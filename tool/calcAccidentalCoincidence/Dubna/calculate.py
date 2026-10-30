#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dubna/Zlokazov：分类泊松计数、时间、串联衰变及相关位置差检验。"""

import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.stats import gamma

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (array, binomial_interval, check_keys, confidence, integer,
                    number, probability_fields, run_cli)


def counts_test(config):
    check_keys(config, {"window_s", "categories"})
    window = number(config["window_s"], "window_s", zero=True)
    categories = config["categories"]
    if not isinstance(categories, list) or not categories:
        raise ValueError("categories 必须为非空数组")
    rows, log_joint = [], 0.0
    for category in categories:
        check_keys(category, {"name", "rate_hz", "count"})
        rate = number(category["rate_hz"], "rate_hz", zero=True)
        count = integer(category["count"], "count")
        mean = number(rate * window, "rate_hz*window_s", zero=True)
        log_p = count * math.log(mean) - mean - math.lgamma(count + 1) if mean else (0.0 if count == 0 else -math.inf)
        log_joint += log_p
        rows.append({"name": category.get("name", ""), "mean": mean, "count": count,
                     **probability_fields(log_p, "Q_exact")})
    return {"categories": rows, **probability_fields(log_joint, "P_exact_configuration"),
            "meaning": "各类独立背景在固定窗内恰好出现指定个数的联合概率；以终端事件为条件时不再乘终端事件概率"}


def time_test(config):
    check_keys(config, {"background_rate_hz", "event_times_s", "gaps_s"})
    rate = number(config["background_rate_hz"], "background_rate_hz")
    if ("event_times_s" in config) == ("gaps_s" in config):
        raise ValueError("event_times_s 和 gaps_s 必须且只能提供一个")
    if "event_times_s" in config:
        times = array(config["event_times_s"], "event_times_s", zero=True)
        gaps = np.diff(times).tolist()
        if not gaps or any(t <= 0 for t in gaps):
            raise ValueError("event_times_s 至少有两个严格递增时间")
    else:
        gaps = array(config["gaps_s"], "gaps_s")
    n = len(gaps)  # 五个链成员对应四个间隔，不能当成五个样本。
    mean = math.fsum(gaps) / n
    harmonic = math.fsum(1 / j for j in range(1, n + 1))
    return {
        "number_of_gaps": n, "observed_mean_gap_s": mean,
        "background_mean_gap_s": 1 / rate,
        "background_expected_min_gap_s": 1 / (n * rate),
        "background_expected_max_gap_s": harmonic / rate,
        "chain_fitted_expected_min_gap_s": mean / n,
        "chain_fitted_expected_max_gap_s": harmonic * mean,
        "P_mean_gap_le_observed_if_iid_exponential": float(gamma.cdf(n * mean * rate, a=n)),
        "meaning": "最小/最大值是次序统计量的期望，不是置信界；均值尾概率另要求未筛选的固定数量独立指数间隔",
    }


def completion_cdf(half_lives, time):
    """独立串联指数等待时间和的 CDF；矩阵指数避免相等衰变常数的奇点。"""
    rates = math.log(2) / np.asarray(half_lives, dtype=float)
    n = len(rates)
    generator = np.zeros((n + 1, n + 1))
    generator[np.arange(n), np.arange(n)] = -rates
    generator[np.arange(n), np.arange(1, n + 1)] = rates
    # 直接读吸收态概率，避免小 t 时用 1-survival 造成数值消减。
    return float(np.clip(expm(generator * time)[0, n], 0, 1))


def decay_test(config):
    check_keys(config, {"time_s", "half_lives_s", "half_life_bounds_s", "detection_efficiency"})
    time = number(config["time_s"], "time_s", zero=True)
    efficiency = number(config.get("detection_efficiency", 1), "detection_efficiency", zero=True)
    if efficiency > 1:
        raise ValueError("detection_efficiency 不能大于 1")
    if ("half_lives_s" in config) == ("half_life_bounds_s" in config):
        raise ValueError("half_lives_s 和 half_life_bounds_s 必须且只能提供一个")
    if "half_lives_s" in config:
        lives = array(config["half_lives_s"], "half_lives_s")
        result = {"P_completed_by_t": completion_cdf(lives, time)}
    else:
        bounds = config["half_life_bounds_s"]
        if not isinstance(bounds, list) or not bounds:
            raise ValueError("half_life_bounds_s 必须为非空的 [下限,上限] 数组")
        checked = [array(b, "half_life_bounds_s") for b in bounds]
        if any(len(b) != 2 or b[0] > b[1] for b in checked):
            raise ValueError("半衰期界限必须为 [正下限,正上限]")
        # 随着每个半衰期增加，完成时间和随机增大，故 CDF 的极值位于两个对角角点。
        result = {"P_completed_min": completion_cdf([b[1] for b in checked], time),
                  "P_completed_max": completion_cdf([b[0] for b in checked], time)}
    return {**result, "detection_efficiency": efficiency,
            "detected_probabilities": {k: v * efficiency for k, v in result.items()},
            "meaning": "在假设串联衰变模型下，截至 t 完成指定全部步骤的 CDF；效率是整链总效率，不自动包含分支比或漏失路径"}


def position_test(config):
    check_keys(config, {"standardized_differences", "correlation_matrix", "samples", "seed", "confidence"})
    x = np.asarray(config["standardized_differences"], dtype=float)
    if x.ndim != 1 or not x.size or not np.isfinite(x).all():
        raise ValueError("standardized_differences 必须为有限的标准化位置差数组")
    covariance = np.asarray(config["correlation_matrix"], dtype=float)
    n = x.size
    if (covariance.shape != (n, n) or not np.isfinite(covariance).all()
            or not np.allclose(covariance, covariance.T, atol=1e-12, rtol=0)
            or not np.allclose(np.diag(covariance), 1, atol=1e-12, rtol=0)):
        raise ValueError("correlation_matrix 必须为对称、对角为 1 的匹配尺寸矩阵")
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues.min() < -1e-10:
        raise ValueError("相关矩阵必须半正定")
    eigenvalues = np.maximum(eigenvalues, 0)
    trials = integer(config.get("samples", 100000), "samples", 1000)
    seed = integer(config.get("seed", 20261029), "seed")
    level = confidence(config.get("confidence", 0.67))
    statistic = float(x @ x)
    rng = np.random.default_rng(seed)
    # X'X 与 sum(lambda_i * Z_i^2) 同分布；保留差分间的相关性。
    simulated = np.concatenate([
        np.sum(rng.standard_normal((min(10000, trials - start), n)) ** 2 * eigenvalues, axis=1)
        for start in range(0, trials, 10000)
    ])
    tail = (1 - level) / 2
    interval = np.quantile(simulated, [tail, 1 - tail]).tolist()
    hits = int(np.count_nonzero(simulated >= statistic))
    return {
        "S": statistic, "expected_S": float(n),
        "variance_S": float(2 * np.sum(eigenvalues ** 2)),
        "confidence": level, "central_reference_interval": interval,
        "inside_reference_interval": bool(interval[0] <= statistic <= interval[1]),
        "P_S_ge_observed": hits / trials,
        "tail_probability_CI95_MC": binomial_interval(hits, trials),
        "samples": trials, "seed": seed,
        "meaning": "共同位置、已知分辨率及给定相关矩阵下的相容性；不是随机链假符合概率",
    }


def calculate(config, base_dir=None):
    functions = {"counts": counts_test, "time": time_test, "decay": decay_test, "position": position_test}
    check_keys(config, functions)
    if not any(name in config for name in functions):
        raise ValueError("至少需要 counts、time、decay、position 中的一项")
    return {"method": "Dubna", **{name: function(config[name]) for name, function in functions.items() if name in config}}


if __name__ == "__main__":
    run_cli(calculate, __file__, "Dubna / Zlokazov")
