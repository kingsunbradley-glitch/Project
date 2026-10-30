"""新增方法共用的输入检查、概率计算与命令行入口。"""

import argparse
import json
import math
from pathlib import Path

from scipy.stats import beta, poisson


def check_keys(config, allowed):
    if not isinstance(config, dict):
        raise ValueError("配置必须为 JSON 对象")
    unknown = set(config) - set(allowed) - {"description"}
    if unknown:
        raise ValueError(f"未知参数：{sorted(unknown)}")


def number(value, name, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} 必须是数值")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (not zero and value == 0):
        raise ValueError(f"{name} 必须为有限{'非负' if zero else '正'}数")
    return value


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} 必须是 >= {minimum} 的整数")
    return value


def array(values, name, zero=False):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} 必须是非空数组")
    return [number(v, name, zero) for v in values]


def confidence(value):
    value = number(value, "confidence")
    if value >= 1:
        raise ValueError("confidence 必须在 0 与 1 之间")
    return value


def probability_fields(log_p, name="P"):
    return {name: math.exp(log_p), f"log10_{name}": log_p / math.log(10)}


def poisson_log_tail(log_mean, observed):
    integer(observed, "observed_chains")
    if observed == 0:
        return 0.0
    if log_mean == -math.inf:
        return -math.inf
    if log_mean > math.log(float.fromhex("0x1.fffffffffffffp+1023")):
        raise ValueError("期望数超出数值范围")
    mean = math.exp(log_mean)
    probability = float(poisson.sf(observed - 1, mean))
    if probability > 0:
        return math.log(probability)
    if not math.isfinite(probability):
        raise ValueError("泊松参数超出支持范围")
    term = total = 1.0
    k = observed
    while True:
        k += 1
        term *= mean / k
        updated = total + term
        if updated == total:
            break
        total = updated
    return -mean + observed * log_mean - math.lgamma(observed + 1) + math.log(total)


def chain_result(log_mean, observed):
    integer(observed, "observed_chains")
    return {
        "n_b": math.exp(log_mean), "log10_n_b": log_mean / math.log(10),
        "observed_chains": observed,
        **probability_fields(poisson_log_tail(log_mean, observed), "P_at_least_observed_chains"),
    }


def binomial_interval(hits, trials, level=0.95):
    """Clopper--Pearson 双侧精确二项区间，含零命中边界。"""
    integer(trials, "trials", 1)
    integer(hits, "hits")
    if hits > trials:
        raise ValueError("hits 不能超过 trials")
    tail = (1 - confidence(level)) / 2
    lower = 0.0 if hits == 0 else float(beta.ppf(tail, hits, trials - hits + 1))
    upper = 1.0 if hits == trials else float(beta.isf(tail, hits + 1, trials - hits))
    return [lower, upper]


def safe_json(value):
    if isinstance(value, dict):
        return {k: safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def run_cli(calculate, script_file, method):
    parser = argparse.ArgumentParser(description=f"{method}；默认读取脚本同目录 input.json")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--json", action="store_true", help="只输出结果 JSON")
    args = parser.parse_args()
    path = args.config or Path(script_file).with_name("input.json")
    try:
        config = json.loads(path.read_text(encoding="utf-8-sig"))
        result = calculate(config, base_dir=path.resolve().parent)
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as exc:
        parser.error(str(exc))
    if not args.json:
        print(f"方法：{method}\n输入：{path.resolve()}")
        if config.get("description"):
            print(config["description"])
    print(json.dumps(safe_json(result), ensure_ascii=False, indent=2, allow_nan=False))
