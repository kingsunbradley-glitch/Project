#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Berkeley：固定总窗泊松计数，或逐级窗口至少一次的乘积。"""

import math
from pathlib import Path
import sys

from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import array, check_keys, integer, number, poisson_log_tail, probability_fields, run_cli, chain_result


def sequential(config):
    """每步从前一个选定事件重新计时；同一筛选类别可重复出现在后续步骤。"""
    check_keys(config, {"model", "step_rates_hz", "step_background_counts", "background_time_s",
                        "step_windows_s", "position_bins", "implantations", "observed_chains"})
    bins = integer(config.get("position_bins", 1), "position_bins", 1)
    seeds = integer(config["implantations"], "implantations")
    windows = array(config["step_windows_s"], "step_windows_s", zero=True)
    if ("step_rates_hz" in config) == ("step_background_counts" in config):
        raise ValueError("step_rates_hz 和 step_background_counts 必须且只能提供一个")
    if "step_rates_hz" in config:
        if "background_time_s" in config:
            raise ValueError("使用 step_rates_hz 时不需要 background_time_s")
        rates = array(config["step_rates_hz"], "step_rates_hz", zero=True)
    else:
        counts = config["step_background_counts"]
        if not isinstance(counts, list) or not counts:
            raise ValueError("step_background_counts 必须为非空整数数组")
        counts = [integer(n, "step_background_counts") for n in counts]
        exposure = number(config["background_time_s"], "background_time_s")
        rates = [number(n / exposure, "step_rate_hz", zero=True) for n in counts]
    if len(rates) != len(windows):
        raise ValueError("每个逐级窗口必须对应一个背景率或背景计数")
    means, probabilities, log_probabilities = [], [], []
    for rate, window in zip(rates, windows):
        log_mu = math.log(rate) + math.log(window) - math.log(bins) if rate and window else -math.inf
        mu = math.exp(log_mu)
        log_p = math.log(-math.expm1(-mu)) if mu else log_mu
        means.append(mu)
        probabilities.append(math.exp(log_p))
        log_probabilities.append(log_p)
    log_p = math.fsum(log_probabilities)
    log_mean = math.log(seeds) + log_p if seeds else -math.inf
    return {"method": "Berkeley", "model": "sequential", "step_rates_hz": rates,
            "step_windows_s": windows, "mu_per_step": means, "P_per_step": probabilities,
            **probability_fields(log_p, "P_per_implantation"),
            **chain_result(log_mean, config.get("observed_chains", 1)),
            "meaning": "选定起始事件后，按每步重启的窗口取至少一次概率的乘积；不是在单个总窗中恰好出现若干事件"}


def calculate(config, base_dir=None):
    check_keys(config, {"model", "alpha_rate_hz", "position_bins", "window_s", "n_alpha",
                        "step_rates_hz", "step_background_counts", "background_time_s", "step_windows_s",
                        "implantations", "count_mode", "observed_chains"})
    model = config.get("model", "total_window")
    if model == "sequential":
        return sequential(config)
    if model != "total_window":
        raise ValueError("model 必须为 total_window 或 sequential")
    check_keys(config, {"model", "alpha_rate_hz", "position_bins", "window_s", "n_alpha",
                        "implantations", "count_mode", "observed_chains"})
    rate = number(config["alpha_rate_hz"], "alpha_rate_hz", zero=True)
    bins = integer(config.get("position_bins", 1), "position_bins", 1)
    window = number(config["window_s"], "window_s", zero=True)
    n = integer(config["n_alpha"], "n_alpha", 1)
    seeds = integer(config["implantations"], "implantations")
    mode = config.get("count_mode", "exact")
    if mode not in ("exact", "at_least"):
        raise ValueError("count_mode 必须为 exact 或 at_least")
    log_mu = math.log(rate) + math.log(window) - math.log(bins) if rate and window else -math.inf
    mu = math.exp(log_mu)
    log_p = (n * log_mu - mu - math.lgamma(n + 1) if mode == "exact"
             else poisson_log_tail(log_mu, n))
    log_mean = math.log(seeds) + log_p if seeds else -math.inf
    return {
        "method": "Berkeley", "model": model, "mu_per_pixel": mu, "n_alpha": n, "count_mode": mode,
        **probability_fields(log_p, "P_per_implantation"),
        **chain_result(log_mean, config.get("observed_chains", 1)),
        "meaning": "P_per_implantation 是单次注入后固定总窗的概率；n_b=注入数×该概率",
    }


if __name__ == "__main__":
    run_cli(calculate, __file__, "Berkeley")
