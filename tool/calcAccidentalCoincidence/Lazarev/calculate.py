#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""截图式 (4.28)：由相邻类别的偶然配对计数估计随机完整链数。"""

import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import array, chain_result, check_keys, integer, run_cli


def calculate(config, base_dir=None):
    check_keys(config, {"pair_counts", "intermediate_counts", "observed_chains"})
    pairs = array(config["pair_counts"], "pair_counts", zero=True)
    counts = array(config.get("intermediate_counts", []), "intermediate_counts") if len(pairs) > 1 else []
    if len(counts) != len(pairs) - 1 or (len(pairs) == 1 and config.get("intermediate_counts", [])):
        raise ValueError("K 类事件需要 K-1 个 pair_counts 和 K-2 个 intermediate_counts")
    observed = integer(config.get("observed_chains", 1), "observed_chains")
    log_mean = (math.fsum(math.log(c) for c in pairs) - math.fsum(math.log(n) for n in counts)
                if all(pairs) else -math.inf)
    return {"method": "Lazarev", "pair_counts": pairs, "intermediate_counts": counts,
            **chain_result(log_mean, observed),
            "meaning": "由偶然配对统计得到 n_b，再计算整个搜索中至少出现观测链数的泊松上尾概率"}


if __name__ == "__main__":
    run_cli(calculate, __file__, "Lazarev")
