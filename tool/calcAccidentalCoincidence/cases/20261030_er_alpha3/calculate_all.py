#!/usr/bin/env python3
"""本次 ER-alpha1-alpha23-alpha23 参数复算；只向终端输出，不改各方法默认输入。"""

import importlib.util
import json
import math
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BASE))
from common import chain_result, safe_json


def method(name):
    spec = importlib.util.spec_from_file_location(f"case_{name}", BASE / name / "calculate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.calculate


def calculate():
    config = json.loads(Path(__file__).with_name("input.json").read_text())
    time = config["measurement_time_s"]
    bins = config["position_bins"]
    er, a, b = [config["background_counts"][key] for key in ("ER", "alpha1", "alpha23")]
    counts = [er, a, b]
    rates = [n / time for n in counts]
    local = [r / bins for r in rates]
    dt1, dt2, dt3 = config["windows_s"]
    window = dt1 + dt2 + dt3
    observed = config["observed_chains"]

    gsi_input = {"mode": "fixed", "measurement_time_s": time, "position_bins": bins,
                 "background_counts": counts, "sequence": [0, 1, 2, 2],
                 "windows_s": [dt1, dt2, dt3], "observed_chains": observed}
    gsi = method("GSI")(gsi_input)

    # 没有实际配对统计，仅在平稳均匀泊松模型下估计期望配对数。
    pairs = [er * local[1] * dt1, a * local[2] * dt2, b * local[2] * dt3]
    lazarev_input = {"pair_counts": pairs, "intermediate_counts": [a, b], "observed_chains": observed}
    lazarev = method("Lazarev")(lazarev_input)
    lazarev["input_status"] = "配对数为本底率模型估计，非实测；结果应与 GSI 低率式一致"

    # 用户实际使用的 Berkeley：每一步在前一步之后的各自窗口内至少出现一次。
    berkeley_sequential_input = {
        "model": "sequential", "step_background_counts": [a, b, b],
        "background_time_s": time, "position_bins": bins,
        "step_windows_s": [dt1, dt2, dt3], "implantations": er,
        "observed_chains": observed,
    }
    berkeley_sequential = method("Berkeley")(berkeley_sequential_input)

    # 旧的宽总窗计数仅留作不同问题的对照，不是本次逐级筛选结果。
    berkeley_input = {"alpha_rate_hz": rates[1] + rates[2], "position_bins": bins,
                      "window_s": window, "n_alpha": 3, "implantations": er,
                      "count_mode": "exact", "observed_chains": observed}
    berkeley = method("Berkeley")(berkeley_input)
    berkeley_at_least = method("Berkeley")({**berkeley_input, "count_mode": "at_least"})

    # 以每个 ER 为锚，在整个总窗内恰有一个 alpha1 和两个 alpha23。
    # ER 本身已是给定条件，不再乘一个 ER 计数概率。
    dubna_input = {"counts": {"window_s": window, "categories": [
        {"name": "alpha1", "rate_hz": local[1], "count": 1},
        {"name": "alpha23", "rate_hz": local[2], "count": 2}]}}
    dubna = method("Dubna")(dubna_input)
    log_p = dubna["counts"]["log10_P_exact_configuration"] * math.log(10)
    dubna["whole_search_poisson_plugin"] = chain_result(math.log(er) + log_p, observed)

    # 固定总窗内恰有 A1/B2 后再要求所有相邻时间门的几何接受度。
    # 有序 B 对在总窗中的均匀密度为 2/window^2；A 均匀密度为 1/window。
    acceptance = 2 * dt1 * dt2 * dt3 / window ** 3
    filtered_log_mean = math.log(er) + log_p + math.log(acceptance)
    filtered = {"conditional_time_acceptance": acceptance,
                **chain_result(filtered_log_mean, observed),
                "meaning": "在固定总窗内 alpha1 恰一个、alpha23 恰两个，且三步均通过相邻窗；允许其他 ER，不允许额外 alpha"}

    return {
        "input": config, "rates_total_hz": dict(zip(("ER", "alpha1", "alpha23"), rates)),
        "rates_per_pixel_hz": dict(zip(("ER", "alpha1", "alpha23"), local)),
        "berkeley_dubna_total_window_s": window,
        "method_inputs": {"GSI": gsi_input, "Lazarev": lazarev_input,
                          "Berkeley": berkeley_sequential_input,
                          "Berkeley_total_window_reference": berkeley_input, "Dubna": dubna_input},
        "Berkeley_sequential": berkeley_sequential,
        "GSI": gsi, "Lazarev_model_pairs": lazarev, "Berkeley_exact3": berkeley,
        "Berkeley_at_least3": berkeley_at_least,
        "Dubna_fixed_window_counts": dubna,
        "fixed_window_counts_with_all_time_gates": filtered,
        "MonteCarlo": {"status": "not_computable_from_aggregate_counts",
                       "required": "真实列表事件的时间、位置和类别，或实际随机插入试验的插入数与命中数",
                       "reason": "汇总计数无法确定时间/位置非均匀性；未用人工数据代替实验结果"},
        "Dubna_other_tests": {"status": "missing_data",
                              "required": "实际链成员时间、候选核素半衰期、位置差和分辨率/相关矩阵"},
    }


if __name__ == "__main__":
    print(json.dumps(safe_json(calculate()), ensure_ascii=False, indent=2, allow_nan=False))
