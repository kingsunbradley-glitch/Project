#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MCRP：逐个插入随机标记事件，在原始列表数据中搜索满足条件的完整链。"""

import csv
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (array, binomial_interval, chain_result, check_keys, confidence,
                    integer, number, run_cli)


def load_events(path, bins):
    events = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"time_s", "position_bin", "event_type"} <= set(reader.fieldnames or []):
            raise ValueError("CSV 需要 time_s,position_bin,event_type 列，可选 energy_mev")
        for line, row in enumerate(reader, 2):
            try:
                time = number(float(row["time_s"]), "time_s", zero=True)
                position = integer(int(row["position_bin"]), "position_bin")
                kind = row["event_type"].strip()
                energy = float(row["energy_mev"]) if row.get("energy_mev", "").strip() else None
                if position >= bins or not kind or (energy is not None and not math.isfinite(energy)):
                    raise ValueError("位置越界、空事件类别或无效能量")
                events.append((time, position, kind, energy))
            except (ValueError, TypeError, AttributeError) as exc:
                raise ValueError(f"{path}:{line}: {exc}") from exc
    if not events:
        raise ValueError("CSV 中没有事件")
    return events


def has_chain(times_by_stage, anchor_time, anchor_stage, windows):
    """允许其他事件插入；严格时间顺序；每次只判断是否存在一条完整链。"""
    def search(stage, time, direction):
        next_stage = stage + direction
        if not 0 <= next_stage < len(times_by_stage):
            return True
        lower_gap, upper_gap = windows[min(stage, next_stage)]
        lower = time + lower_gap if direction == 1 else time - upper_gap
        upper = time + upper_gap if direction == 1 else time - lower_gap
        candidates = times_by_stage[next_stage]
        first = np.searchsorted(candidates, lower, side="left")
        last = np.searchsorted(candidates, upper, side="right")
        for candidate in candidates[first:last]:
            if (candidate - time) * direction > 0 and search(next_stage, candidate, direction):
                return True
        return False

    return search(anchor_stage, anchor_time, -1) and search(anchor_stage, anchor_time, 1)


def simulate(config, base_dir):
    bins = integer(config.get("position_bins", 1), "position_bins", 1)
    trials = integer(config["trials"], "trials", 1)
    seed = integer(config.get("seed", 20261029), "seed")
    stages = config["stages"]
    if not isinstance(stages, list) or len(stages) < 2:
        raise ValueError("stages 至少需要两类事件")
    filters = []
    for stage in stages:
        check_keys(stage, {"event_type", "energy_mev"})
        if not isinstance(stage["event_type"], str) or not stage["event_type"].strip():
            raise ValueError("event_type 必须为非空字符串")
        energy = array(stage["energy_mev"], "energy_mev", zero=True) if "energy_mev" in stage else None
        if energy is not None and (len(energy) != 2 or energy[0] > energy[1]):
            raise ValueError("energy_mev 必须为 [能量下限,能量上限]")
        filters.append((stage["event_type"].strip(), energy))
    anchor = integer(config.get("injection_index", 0), "injection_index")
    if anchor >= len(stages):
        raise ValueError("injection_index 超过 stages 范围（从 0 编号）")
    windows = config["windows_s"]
    if not isinstance(windows, list) or len(windows) != len(stages) - 1:
        raise ValueError("windows_s 必须为 K-1 个 [下限,上限]")
    windows = [array(w, "windows_s", zero=True) for w in windows]
    if any(len(w) != 2 or w[0] > w[1] for w in windows):
        raise ValueError("时间窗必须为 [非负下限,上限]")
    intervals = config["live_intervals_s"]
    if not isinstance(intervals, list) or not intervals:
        raise ValueError("live_intervals_s 必须为非空数组")
    intervals = [array(i, "live_intervals_s", zero=True) for i in intervals]
    if any(len(i) != 2 or i[0] >= i[1] for i in intervals):
        raise ValueError("每个采集区段必须为 [开始,结束]，且开始 < 结束")
    if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
        raise ValueError("采集区段必须按时间排序且不重叠")
    time_weights = np.asarray(array(config.get("injection_interval_weights", [b-a for a, b in intervals]), "injection_interval_weights", zero=True))
    position_weights = np.asarray(array(config.get("position_weights", [1] * bins), "position_weights", zero=True))
    if len(time_weights) != len(intervals) or len(position_weights) != bins:
        raise ValueError("时间区段权重/位置权重的长度不匹配")
    if not np.isfinite(time_weights.sum()) or not np.isfinite(position_weights.sum()) or not time_weights.sum() or not position_weights.sum():
        raise ValueError("权重总和必须为有限正数")
    time_weights = time_weights / time_weights.sum()
    position_weights = position_weights / position_weights.sum()
    events = load_events(base_dir / config["events_csv"], bins)
    if any(energy is None for _, _, _, energy in events) and any(f[1] for f in filters):
        raise ValueError("存在能量门时，CSV 所有事件都需要 energy_mev")
    # 每个连续采集区段、位置区、链步骤各自保留真实列表数据的时间结构。
    tables = []
    for start, end in intervals:
        by_position = []
        for position in range(bins):
            by_position.append([
                np.sort(np.array([t for t, p, k, e in events
                                  if start <= t < end and p == position and k == kind
                                  and (energy is None or energy[0] <= e <= energy[1])]))
                for kind, energy in filters
            ])
        tables.append(by_position)
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(trials):
        interval = int(rng.choice(len(intervals), p=time_weights))
        position = int(rng.choice(bins, p=position_weights))
        time = rng.uniform(*intervals[interval])
        hits += int(has_chain(tables[interval][position], time, anchor, windows))
    return hits, trials, {"seed": seed, "injection_index": anchor, "background_events": len(events)}


def calculate(config, base_dir=None):
    check_keys(config, {"mode", "inserted_events", "correlated_events", "confidence", "reference_events",
                        "observed_chains", "events_csv", "trials", "seed", "stages", "windows_s",
                        "position_bins", "position_weights", "live_intervals_s", "injection_index",
                        "injection_interval_weights"})
    mode = config.get("mode", "simulate")
    level = confidence(config.get("confidence", 0.95))
    if mode == "counts":
        check_keys(config, {"mode", "inserted_events", "correlated_events", "confidence",
                            "reference_events", "observed_chains"})
        trials = integer(config["inserted_events"], "inserted_events", 1)
        hits = integer(config["correlated_events"], "correlated_events")
        extra = {}
    elif mode == "simulate":
        if "inserted_events" in config or "correlated_events" in config:
            raise ValueError("simulate 模式不能同时提供已有插入/命中计数")
        hits, trials, extra = simulate(config, Path(base_dir or "."))
    else:
        raise ValueError("mode 必须为 counts 或 simulate")
    interval = binomial_interval(hits, trials, level)
    result = {"method": "MonteCarlo", "mode": mode, "inserted_events": trials,
              "correlated_events": hits, "P_per_insertion": hits / trials,
              "confidence": level, "P_per_insertion_interval": interval, **extra,
              "meaning": "每个随机标记事件至少形成一条满足搜索条件的链的频率；一个标记最多计一次命中"}
    if "reference_events" in config:
        count = integer(config["reference_events"], "reference_events")
        log_mean = math.log(count) + math.log(hits / trials) if count and hits else -math.inf
        estimate = chain_result(log_mean, config.get("observed_chains", 1))
        # 这是将模拟得到的概率代入泊松模型后的估计，不能误认为来自独立全实验模拟。
        estimate["n_b_interval_from_MC"] = [count * p for p in interval]
        result["poisson_plugin_estimate"] = estimate
    elif "observed_chains" in config:
        raise ValueError("observed_chains 需要同时设置 reference_events")
    if hits == 0:
        result["note"] = "零命中只说明本次样本未检出，请使用区间上限；代入估计中的零不是严格概率上界"
    return result


if __name__ == "__main__":
    run_cli(calculate, __file__, "Monte Carlo random insertion")
