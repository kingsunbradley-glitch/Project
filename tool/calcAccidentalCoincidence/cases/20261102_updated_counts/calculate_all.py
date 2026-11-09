#!/usr/bin/env python3
"""新计数分别复算含/不含 SF 的关联链；只输出结果，不改旧案例。"""

from collections import Counter
import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
sys.path.insert(0, str(ROOT))
from common import chain_result, safe_json


def load_method(name):
    spec = importlib.util.spec_from_file_location(f"updated_counts_{name}", ROOT / name / "calculate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.calculate


def calculate(config_path=None):
    source = Path(config_path) if config_path is not None else CASE / "input.json"
    config = json.loads(source.read_text(encoding="utf-8"))
    methods = {name: load_method(name) for name in ("GSI", "Berkeley", "Lazarev", "Dubna")}
    time, bins = config["measurement_time_s"], config["position_bins"]
    counts = config["background_counts"]
    rate = {key: n / time / bins for key, n in counts.items()}
    er, observed = counts["ER"], config["observed_chains"]
    cases = {}
    for with_sf in (False, True):
        sequence = ["ER", "alpha1", "alpha23", "alpha23"] + (["SF"] if with_sf else [])
        groups = list(dict.fromkeys(sequence))
        windows = config["step_windows_s"][:len(sequence) - 1]
        gsi_input = {"mode": "fixed", "measurement_time_s": time, "position_bins": bins,
                     "background_counts": [counts[key] for key in groups],
                     "sequence": [groups.index(key) for key in sequence],
                     "windows_s": windows, "observed_chains": observed}
        gsi = methods["GSI"](gsi_input)
        berkeley_input = {"model": "sequential", "background_time_s": time, "position_bins": bins,
                          "step_background_counts": [counts[key] for key in sequence[1:]],
                          "step_windows_s": windows, "implantations": er, "observed_chains": observed}
        berkeley = methods["Berkeley"](berkeley_input)

        pair_counts = [counts[left] * rate[right] * window
                       for left, right, window in zip(sequence[:-1], sequence[1:], windows)]
        lazarev_input = {"pair_counts": pair_counts,
                         "intermediate_counts": [counts[key] for key in sequence[1:-1]],
                         "observed_chains": observed}
        lazarev = methods["Lazarev"](lazarev_input)
        lazarev["input_status"] = "未提供实测相邻配对数；采用均匀泊松背景率估计期望配对数"

        # 原计数式之外，显式施加全部相邻时间门；同类的两个 alpha23 作为 K=2。
        multiplicities = Counter(sequence[1:])
        span = math.fsum(windows)
        dubna_input = {"counts": {"window_s": span, "categories": [
            {"name": key, "rate_hz": rate[key], "count": n} for key, n in multiplicities.items()]}}
        dubna_raw = methods["Dubna"](dubna_input)
        log_p_raw = dubna_raw["counts"]["log10_P_exact_configuration"] * math.log(10)
        time_acceptance = math.prod(math.factorial(n) for n in multiplicities.values())
        time_acceptance *= math.prod(windows) / span ** len(windows)
        dubna = chain_result(math.log(er) + log_p_raw + math.log(time_acceptance), observed)
        dubna["time_acceptance"] = time_acceptance
        dubna["model_definition"] = (
            "ER 为锚，在固定最大总跨度内恰好出现指定类别个数，再显式要求全部相邻时间门；"
            "允许额外 ER，不允许额外所列衰变类别事件。非原始计数式未经时间筛选的值。"
        )

        # 与独立、直接的解析表达式核对。
        total = math.fsum(rate[key] for key in groups)
        expected_gsi = er * math.prod(rate[key] / total * -math.expm1(-total * window)
                                     for key, window in zip(sequence[1:], windows))
        expected_low = er * math.prod(rate[key] * window for key, window in zip(sequence[1:], windows))
        expected_berkeley = er * math.prod(-math.expm1(-rate[key] * window)
                                          for key, window in zip(sequence[1:], windows))
        expected_dubna = expected_low * math.exp(-span * math.fsum(rate[key] for key in multiplicities))
        for actual, expected in [(gsi["n_b"], expected_gsi), (lazarev["n_b"], expected_low),
                                 (berkeley["n_b"], expected_berkeley), (dubna["n_b"], expected_dubna)]:
            assert math.isclose(actual, expected, rel_tol=1e-12)

        cases["with_SF" if with_sf else "without_SF"] = {
            "sequence": sequence,
            "inputs": {"GSI": gsi_input, "Berkeley": berkeley_input,
                       "Lazarev": lazarev_input, "Dubna_raw_counts": dubna_input},
            "GSI": gsi, "Berkeley": berkeley, "Lazarev_model_pairs": lazarev,
            "Dubna_counts_plus_time_gates": dubna,
            "Dubna_raw_counts_reference_only": dubna_raw,
            "MonteCarlo": {"status": "missing_event_list_or_insertion_statistics"},
        }
    return {"input": config, "rates_per_pixel_hz": rate, "cases": cases,
            "assumptions": ["655.2 h 和 6144 个位置单元沿用前次，背景时间/位置均匀",
                            "ER、alpha1、alpha23、SF 作为互不重叠类别；alpha23 在链中引用两次",
                            "求至少一条链的概率；各 n_b 代入泊松链数分布",
                            "固定背景率点估计，不传播背景率统计不确定度"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="指定新的计数配置；默认读取同目录 input.json")
    args = parser.parse_args()
    print(json.dumps(safe_json(calculate(args.config)), ensure_ascii=False, indent=2, allow_nan=False))
