#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Schmidt et al. (1984), Eqs. (2), (3), (5)--(8): 随机符合链概率。

运行：python calculate.py --config input.json
所有时间单位为 s，计数率为 Hz；输入率为所分析区域的总率。
依赖 scipy。公式、事件筛选约定和统计适用条件见 README.md。
"""

import argparse
import json
import math
from pathlib import Path

from scipy.stats import poisson


EXAMPLES = {
    "th214-5s": {
        "mode": "fixed", "measurement_time_s": 6080,
        "background_counts": [28, 38], "windows_s": [5],
        "observed_chains": 12, "position_bins": 1,
    },
    "th214-40s": {
        "mode": "fixed", "measurement_time_s": 6080,
        "background_counts": [28, 38], "windows_s": [40],
        "observed_chains": 18, "position_bins": 1,
    },
    "element109": {
        "mode": "partial", "measurement_time_s": 250 * 3600,
        "background_counts": [288, 1, 1, 1],
        "background_time_s": [250 * 3600, 8 * 3600, 75, 250 * 3600],
        "windows_s": [0.005, 0.0273, 12.9],
        "observed_chains": 1, "position_bins": 50,
    },
}


def positive_number(value, name, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} 必须是数值")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f"{name} 必须为有限{'非负' if allow_zero else '正'}数")
    return value


def integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} 必须是 >= {minimum} 的整数")
    return value


def numbers(values, name, allow_zero=False):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} 必须是非空数组")
    return [positive_number(v, name, allow_zero) for v in values]


def poisson_tail(log_mean, observed):
    """返回式 (2) 的概率及其自然对数；下溢时仍保留 log(P)。"""
    if observed == 0:
        return 1.0, 0.0
    if log_mean == -math.inf:
        return 0.0, -math.inf
    mean = math.exp(log_mean)
    probability = float(poisson.sf(observed - 1, mean))
    if probability > 0:
        return probability, math.log(probability)
    # 极小尾概率：从 n=observed 项起递推，避免 1-CDF 的消减误差。
    term = total = 1.0
    index = observed
    while True:
        index += 1
        term *= mean / index
        updated = total + term
        if updated == total:
            break
        total = updated
    log_probability = -mean + observed * log_mean - math.lgamma(observed + 1)
    log_probability += math.log(total)
    return math.exp(log_probability), log_probability


def calculate(config):
    """从配置字典计算；可作为 Python API 导入，不修改输入。"""
    allowed = {
        "mode", "measurement_time_s", "rates_hz", "background_counts",
        "background_time_s", "windows_s", "observed_chains", "position_bins",
        "sequence",
    }
    if not isinstance(config, dict):
        raise ValueError("配置顶层必须是 JSON 对象")
    if set(config) - allowed:
        raise ValueError(f"未知配置项：{sorted(set(config) - allowed)}")
    mode = config.get("mode", "fixed")
    if mode not in ("fixed", "partial"):
        raise ValueError("mode 只能为 fixed 或 partial")
    time = positive_number(config["measurement_time_s"], "measurement_time_s")
    observed = integer(config.get("observed_chains", 1), "observed_chains", 0)
    bins = integer(config.get("position_bins", 1), "position_bins", 1)
    windows = numbers(config["windows_s"], "windows_s", allow_zero=True)
    if ("rates_hz" in config) == ("background_counts" in config):
        raise ValueError("rates_hz 和 background_counts 必须且只能提供一个")
    if "rates_hz" in config:
        if "background_time_s" in config:
            raise ValueError("使用 rates_hz 时不要设置 background_time_s")
        rates = numbers(config["rates_hz"], "rates_hz", allow_zero=True)
    else:
        counts = config["background_counts"]
        if not isinstance(counts, list) or not counts:
            raise ValueError("background_counts 必须为非空整数数组")
        counts = [integer(n, "background_counts", 0) for n in counts]
        exposures = config.get("background_time_s", time)
        if not isinstance(exposures, list):
            exposures = [exposures] * len(counts)
        exposures = numbers(exposures, "background_time_s")
        if len(exposures) != len(counts):
            raise ValueError("background_time_s 数组长度必须与 background_counts 相同")
        rates = [n / t for n, t in zip(counts, exposures)]
    sequence = config.get("sequence", list(range(len(rates))))
    if not isinstance(sequence, list) or len(sequence) < 2:
        raise ValueError("sequence 必须包含至少两个事件类别下标")
    sequence = [integer(i, "sequence", 0) for i in sequence]
    if any(i >= len(rates) for i in sequence):
        raise ValueError("sequence 下标超出计数率数组范围")
    if len(windows) != len(sequence) - 1:
        raise ValueError("时间窗数量必须为链步骤数减一")
    if mode == "partial" and len(set(sequence)) != len(sequence):
        raise ValueError("partial 模式不支持重复事件类别；请使用 fixed 模式")
    if not all(math.isfinite(r) for r in rates):
        raise ValueError("计数率超出浮点范围，请检查输入单位")

    # 输入全区域计数率。M 个等效位置区均匀分布：每区 r_i=R_i/M；
    # 最后对 M 区求和，起始事件总数仍为 R_1*T。
    local_rates = [r / bins for r in rates]
    if mode == "fixed":
        loss_rates = [math.fsum(local_rates)] * len(windows)
    else:
        loss_rates = [local_rates[sequence[0]] + local_rates[i] for i in sequence[1:]]
    occupancy = [rate * dt for rate, dt in zip(loss_rates, windows)]
    if not all(math.isfinite(x) for x in occupancy):
        raise ValueError("计数率×时间窗超出浮点范围，请检查输入")

    zero = any(rates[i] == 0 for i in sequence) or any(dt == 0 for dt in windows)
    if zero:
        log_mean = log_approx = -math.inf
    else:
        # 式 (6)/(8) 的对数；随后乘 (1-exp(-x))/x 得到式 (5)/(7)。
        log_approx = math.log(time) + math.fsum(math.log(rates[i]) for i in sequence)
        log_approx -= (len(sequence) - 1) * math.log(bins)
        log_approx += math.fsum(math.log(dt) for dt in windows)
        corrections = [
            math.log(-math.expm1(-x) / x) if x else 0.0 for x in occupancy
        ]
        log_mean = log_approx + math.fsum(corrections)
    mean, approx_mean = math.exp(log_mean), math.exp(log_approx)
    probability, log_probability = poisson_tail(log_mean, observed)
    approx_probability, _ = poisson_tail(log_approx, observed)
    log_small_p = (observed * log_mean - math.lgamma(observed + 1)
                   if observed else 0.0)
    warnings = []
    if max(occupancy) >= 0.1:
        warnings.append("部分 λΔt >= 0.1：低计数率近似可能有明显偏差，以非展开式结果为主。")
    if mean >= 0.1:
        warnings.append("n_b >= 0.1：不宜直接把 n_b 或式 (3) 当作概率。")
    if any(r == 0 for r in rates):
        warnings.append("含零背景率：代入值给出 n_b=0，但零观测不证明真实背景率为零。")
    span = sum(windows) if mode == "fixed" else max(windows)
    if span / time >= 0.01:
        warnings.append("时间窗相对测量时长不可忽略；本计算未修正起止边界或采集间断。")
    if probability == 0 and math.isfinite(log_probability):
        warnings.append("P_err 小于浮点表示范围；请使用 log10_P_err，不能解读为严格零。")
    return {
        "mode": mode, "mean_equation": 5 if mode == "fixed" else 7,
        "measurement_time_s": time, "rates_hz": rates,
        "windows_s": windows, "observed_chains": observed, "position_bins": bins,
        "sequence": sequence,
        "occupancy": occupancy,
        "n_b": mean, "log10_n_b": log_mean / math.log(10),
        "n_b_low_rate": approx_mean,
        "P_err": probability, "log10_P_err": log_probability / math.log(10),
        "P_err_low_rate": approx_probability,
        "log10_P_eq3": log_small_p / math.log(10), "warnings": warnings,
    }


def json_safe(value):
    """严格 JSON 不支持 Infinity；用字符串表示对数零值 -Infinity。"""
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--config", type=Path, help="JSON 参数文件；默认脚本同目录 input.json")
    source.add_argument("--example", choices=EXAMPLES, help="运行论文例子（非用户实验数据）")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出配置与计算结果")
    args = parser.parse_args()
    try:
        if args.example:
            config = EXAMPLES[args.example]
            label = f"论文示例 {args.example}（非用户实验数据）"
        else:
            path = args.config or Path(__file__).with_name("input.json")
            config = json.loads(path.read_text(encoding="utf-8-sig"))
            label = str(path.resolve())
        result = calculate(config)
    except (OSError, ValueError, KeyError, OverflowError) as exc:
        parser.error(f"输入或数值范围错误：{exc}")
    if args.json:
        print(json.dumps(json_safe({"source": label, **result}), indent=2,
                         ensure_ascii=False, allow_nan=False))
        return
    print(f"输入：{label}")
    print(f"模式：{result['mode']}；n_b 使用原文式 ({result['mean_equation']})")
    print(f"测量时长 = {result['measurement_time_s']:.8g} s；位置分区 = {result['position_bins']}")
    print(f"全区域背景率 (Hz) = {result['rates_hz']}")
    print(f"链步骤对应类别下标 = {result['sequence']}")
    print(f"时间窗 (s) = {result['windows_s']}；观测链数 = {result['observed_chains']}")
    print(f"随机链期望数 n_b       = {result['n_b']:.10e}")
    print(f"P_err = P(N >= n_m)    = {result['P_err']:.10e}")
    print(f"log10(P_err)          = {result['log10_P_err']:.10g}")
    print(f"低计数率近似 n_b       = {result['n_b_low_rate']:.10e}")
    print(f"近似 n_b 代入式 (2)    = {result['P_err_low_rate']:.10e}")
    print(f"式 (3) 的 log10(P)    = {result['log10_P_eq3']:.10g}（仅供 n_b << 1 时参考）")
    print(f"max(λΔt)              = {max(result['occupancy']):.6g}")
    print("概率含义：纯随机背景模型下出现至少该数量链的概率；不是该链为假的后验概率。")
    for warning in result["warnings"]:
        print(f"注意：{warning}")


if __name__ == "__main__":
    main()
