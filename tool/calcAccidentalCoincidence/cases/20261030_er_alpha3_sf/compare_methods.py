#!/usr/bin/env python3
"""用已确认的带 SF 输入比较各模型；保留四个相邻时间门，不用宽总窗结果代替。"""

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
    spec = importlib.util.spec_from_file_location(f"sf_compare_{name}", ROOT / name / "calculate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.calculate


def calculate_all():
    methods = {name: load_method(name) for name in ("Berkeley", "GSI", "Lazarev", "Dubna")}
    cases = {}
    for sf_count in (188600, 194):
        berkeley_input = json.loads((CASE / f"SF_{sf_count}.json").read_text(encoding="utf-8"))
        time = berkeley_input["background_time_s"]
        bins = berkeley_input["position_bins"]
        er = berkeley_input["implantations"]
        a, b, b_again, sf = berkeley_input["step_background_counts"]
        if b != b_again:
            raise ValueError("本例 alpha2 和 alpha3 应使用同一筛选类别的计数")
        windows = berkeley_input["step_windows_s"]
        observed = berkeley_input["observed_chains"]
        rates = [n / time / bins for n in (er, a, b, sf)]
        _, r_a, r_b, r_sf = rates

        gsi_input = {"mode": "fixed", "background_counts": [er, a, b, sf],
                     "measurement_time_s": time, "position_bins": bins,
                     "sequence": [0, 1, 2, 2, 3], "windows_s": windows,
                     "observed_chains": observed}
        gsi = methods["GSI"](gsi_input)
        # 独立直接计算下一事件的积分乘积，防止重复 alpha23 被计入总率两次。
        total_rate = math.fsum(rates)
        expected_gsi = er * math.prod(r / total_rate * -math.expm1(-total_rate * dt)
                                     for r, dt in zip((r_a, r_b, r_b, r_sf), windows))
        assert math.isclose(gsi["n_b"], expected_gsi, rel_tol=1e-12)

        # 缺实测配对统计，按均匀独立泊松背景推算配对期望数。
        pairs = [er * r_a * windows[0], a * r_b * windows[1],
                 b * r_b * windows[2], b * r_sf * windows[3]]
        lazarev_input = {"pair_counts": pairs, "intermediate_counts": [a, b, b],
                         "observed_chains": observed}
        lazarev = methods["Lazarev"](lazarev_input)
        lazarev["input_status"] = "各相邻配对数为率模型估计，尚未提供实测配对统计"
        assert math.isclose(lazarev["n_b"], gsi["n_b_low_rate"], rel_tol=1e-12)

        # Dubna 固定窗计数并不等于逐级搜索，须显式乘以有序时间条件的接受度。
        span = math.fsum(windows)
        dubna_input = {"counts": {"window_s": span, "categories": [
            {"name": "alpha1", "rate_hz": r_a, "count": 1},
            {"name": "alpha23", "rate_hz": r_b, "count": 2},
            {"name": "SF", "rate_hz": r_sf, "count": 1}]}}
        dubna_counts = methods["Dubna"](dubna_input)
        log_p_counts = dubna_counts["counts"]["log10_P_exact_configuration"] * math.log(10)
        # 一个 A、两个同类且已排序 B、一个 SF 的条件时间密度为 2/span^4。
        acceptance = 2 * math.prod(windows) / span ** 4
        dubna = chain_result(math.log(er) + log_p_counts + math.log(acceptance), observed)
        dubna["time_acceptance"] = acceptance
        dubna["model_definition"] = (
            "ER 为锚；固定最大总跨度内恰好 A1、B2、SF1，并显式要求 ER-A-B-B-SF 顺序和全部相邻时间门；"
            "允许额外 ER，不允许额外 A/B/SF。属于对 Dubna 计数式添加时间选择的扩展。"
        )
        direct_filtered = gsi["n_b_low_rate"] * math.exp(-(r_a + r_b + r_sf) * span)
        assert math.isclose(dubna["n_b"], direct_filtered, rel_tol=1e-12)

        cases[str(sf_count)] = {
            "inputs": {"Berkeley": berkeley_input, "GSI": gsi_input,
                       "Lazarev": lazarev_input, "Dubna_raw_counts": dubna_input},
            "Berkeley_sequential": methods["Berkeley"](berkeley_input),
            "GSI_fixed_no_intervening_events": gsi,
            "Lazarev_model_pair_counts": lazarev,
            "Dubna_counts_with_all_sequential_time_gates": dubna,
            "Dubna_raw_counts_reference_not_matching_search": dubna_counts,
            "MonteCarlo": {"status": "missing_event_list_or_insertion_statistics"},
            "assumptions": [
                "所有计数共用 655.2 h 全探测器曝光，6144 个均匀位置单元",
                "ER、alpha1、alpha23、SF 为互不重叠背景类别，alpha23 在链中使用两次",
                "固定背景率估计，不传播本底计数的统计不确定度",
                "至少一条概率为将 n_b 代入泊松链数分布的结果，不是链为假的后验概率",
            ],
        }
    return cases


if __name__ == "__main__":
    print(json.dumps(safe_json(calculate_all()), ensure_ascii=False, indent=2, allow_nan=False))
