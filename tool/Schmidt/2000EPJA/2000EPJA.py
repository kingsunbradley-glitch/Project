#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Schmidt test for decay times from an exponential distribution.
Based on K.H. Schmidt, Eur. Phys. J. A 8, 141–145 (2000).

功能：
1. 读取一列衰变时间；
2. 计算平均寿命 tau 和半衰期 T1/2；
3. 计算 Schmidt 文献中的 sigma_Theta_exp；
4. 与 Schmidt Table 1 的 90% 接受区间比较；
5. 可选 Monte Carlo 给出 p_low 和 p_high；
6. 可选绘制 logarithmic decay-time distribution。

运行例子：
    python 2000EPJA.py "Caldata.txt" --unit ms --mc 1000000 --plot
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np


# Schmidt Table 1:
# n : (E_n, sigma_n, gamma_n, lower_5%, upper_95%)
#
# lower_5% 和 upper_95% 是 90% 接受区间。
# 若 sigma_Theta_exp < lower，说明时间分布过窄；
# 若 sigma_Theta_exp > upper，说明时间分布过宽。
SCHMIDT_TABLE = {
    1:  (0.00, 0.00, 0.00, None, None),
    2:  (0.69, 0.58, 1.42, 0.04, 1.83),
    3:  (0.89, 0.55, 1.24, 0.19, 1.91),
    4:  (0.98, 0.50, 1.13, 0.31, 1.92),
    5:  (1.04, 0.47, 1.12, 0.41, 1.90),
    6:  (1.08, 0.44, 1.10, 0.48, 1.89),
    7:  (1.11, 0.42, 0.99, 0.52, 1.87),
    8:  (1.13, 0.40, 0.96, 0.58, 1.85),
    9:  (1.15, 0.38, 0.95, 0.62, 1.84),
    10: (1.16, 0.37, 0.90, 0.65, 1.82),
    11: (1.17, 0.35, 0.84, 0.67, 1.81),
    12: (1.18, 0.34, 0.84, 0.70, 1.79),
    13: (1.19, 0.33, 0.82, 0.72, 1.77),
    14: (1.19, 0.32, 0.78, 0.73, 1.77),
    15: (1.20, 0.31, 0.78, 0.75, 1.76),
    16: (1.20, 0.30, 0.76, 0.77, 1.75),
    17: (1.21, 0.30, 0.74, 0.78, 1.74),
    18: (1.22, 0.29, 0.72, 0.79, 1.73),
    19: (1.22, 0.28, 0.69, 0.80, 1.72),
    20: (1.22, 0.28, 0.68, 0.81, 1.71),
    30: (1.24, 0.23, 0.57, 0.89, 1.64),
    40: (1.25, 0.20, 0.55, 0.94, 1.60),
    50: (1.25, 0.20, 0.55, 0.98, 1.57),
    60: (1.26, 0.17, 0.44, 1.00, 1.54),
    70: (1.26, 0.15, 0.45, 1.02, 1.53),
    80: (1.27, 0.15, 0.45, 1.04, 1.51),
    90: (1.27, 0.14, 0.40, 1.05, 1.50),
    100: (1.27, 0.13, 0.37, 1.06, 1.49),
}


def read_times(path: str | Path) -> np.ndarray:
    """
    读取衰变时间。

    支持这种格式：
        T1
        6.37
        4.75
        1.20ß

    非数字表头会自动跳过。
    行尾奇怪字符也会自动忽略，比如 1.20ß 会读成 1.20。
    """
    text = Path(path).read_text(encoding="utf-8", errors="ignore")

    vals: list[float] = []

    for line in text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        # 只读取一行开头的数字部分
        m = re.match(
            r"^\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?",
            line,
        )

        if not m:
            continue

        x = float(m.group(0))

        if x > 0:
            vals.append(x)

    if len(vals) == 0:
        raise ValueError(f"No positive time values found in {path}")

    return np.asarray(vals, dtype=float)


def sigma_theta_exp(times: np.ndarray) -> float:
    """
    Schmidt eq. (12)

    Theta_i = ln(t_i)

    sigma_Theta_exp =
        sqrt( sum_i (Theta_i - mean(Theta))^2 / n )

    注意这里除以 n，不是 n-1。
    """
    theta = np.log(times)
    theta_mean = theta.mean()

    sigma = np.sqrt(np.mean((theta - theta_mean) ** 2))

    return float(sigma)


def mean_lifetime_and_half_life(times: np.ndarray) -> tuple[float, float]:
    """
    在完整时间窗口的假设下：

        tau = mean(t_i)
        T1/2 = ln(2) * tau
    """
    tau = float(np.mean(times))
    t12 = math.log(2.0) * tau

    return tau, t12

def schmidt_limits(
    n: int,
    nmc: int = 1_000_000,
    seed: int | None = 12345,
) -> tuple[float | None, float | None, str]:
    """
    Schmidt 90% 接受区间。

    文献 Table 1 中有的 n，直接使用表值。
    文献 Table 1 中没有的 n，不做线性插值，
    而是严格按照文献方法 Monte Carlo 计算 5% 和 95% 分位数。

    对 n=1，sigma_Theta_exp 恒为 0，无法做宽度检验。
    """
    if n in SCHMIDT_TABLE:
        lower = SCHMIDT_TABLE[n][3]
        upper = SCHMIDT_TABLE[n][4]
        return lower, upper, "Schmidt Table 1"

    if n < 2:
        return None, None, "not available for n < 2"

    rng = np.random.default_rng(seed)

    batch = min(200_000, max(10_000, nmc))
    sigmas = np.empty(nmc, dtype=float)

    done = 0

    while done < nmc:
        m = min(batch, nmc - done)

        # 时间尺度不影响 ln(t) 的标准差，所以用 tau=1 的指数分布即可
        t = rng.exponential(scale=1.0, size=(m, n))

        theta = np.log(t)

        sigma = np.sqrt(
            np.mean(
                (theta - theta.mean(axis=1, keepdims=True)) ** 2,
                axis=1,
            )
        )

        sigmas[done:done + m] = sigma
        done += m

    lower = float(np.quantile(sigmas, 0.05))
    upper = float(np.quantile(sigmas, 0.95))

    return lower, upper, f"Monte Carlo according to Schmidt method, {nmc} trials"

def mc_pvalues(
    n: int,
    sigma_obs: float,
    nmc: int,
    seed: int | None = None,
) -> tuple[float, float]:
    """
    Monte Carlo 计算：
        p_low  = P(sigma <= sigma_obs)
        p_high = P(sigma >= sigma_obs)

    因为 Schmidt 检验对时间尺度不敏感，所以直接用 Exp(1) 生成随机数即可。
    """
    rng = np.random.default_rng(seed)

    batch = min(200_000, max(10_000, nmc))

    n_lower_or_equal = 0
    n_higher_or_equal = 0
    done = 0

    while done < nmc:
        m = min(batch, nmc - done)

        t = rng.exponential(scale=1.0, size=(m, n))
        theta = np.log(t)

        sigma = np.sqrt(
            np.mean(
                (theta - theta.mean(axis=1, keepdims=True)) ** 2,
                axis=1,
            )
        )

        n_lower_or_equal += int(np.count_nonzero(sigma <= sigma_obs))
        n_higher_or_equal += int(np.count_nonzero(sigma >= sigma_obs))

        done += m

    p_low = n_lower_or_equal / nmc
    p_high = n_higher_or_equal / nmc

    return p_low, p_high


def plot_log_time(
    times: np.ndarray,
    out: str | Path,
    unit: str = "arb.",
) -> None:
    """
    绘制 logarithmic decay-time distribution。

    横轴是 t 的 log 坐标。

    对单一指数衰变：
        dN / dln(t) = N * lambda * t * exp(-lambda * t)

    其中：
        lambda = 1 / tau
    """
    import matplotlib.pyplot as plt

    n = len(times)

    tau = float(np.mean(times))
    lam = 1.0 / tau

    tmin = min(times) / 3.0
    tmax = max(times) * 3.0

    grid = np.logspace(np.log10(tmin), np.log10(tmax), 600)

    y = n * lam * grid * np.exp(-lam * grid)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))

    ax.plot(
        grid,
        y,
        lw=1.8,
        label=r"$N\lambda t e^{-\lambda t}$",
    )

    ax.scatter(
        times,
        np.zeros_like(times),
        marker="v",
        s=55,
        label="data",
    )

    ax.set_xscale("log")
    ax.set_xlabel(f"t / {unit}")
    ax.set_ylabel(r"counts per $\ln(t)$ bin, arb.")

    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Schmidt test for exponential decay-time data"
    )

    parser.add_argument(
        "input",
        help="text file containing one positive decay time per line",
    )

    parser.add_argument(
        "--unit",
        default="arb.",
        help="time unit label only, e.g. ms, s, us",
    )

    parser.add_argument(
        "--mc",
        type=int,
        default=200_000,
        help="Monte Carlo trials for p-values; use 0 to skip",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="random seed for Monte Carlo",
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="save log-time distribution plot",
    )

    parser.add_argument(
        "--out",
        default="schmidt_logtime.png",
        help="plot filename",
    )

    args = parser.parse_args()

    times = read_times(args.input)

    n = len(times)

    sig = sigma_theta_exp(times)

    tau, t12 = mean_lifetime_and_half_life(times)

    lower, upper, limit_source = schmidt_limits(
    n=n,
    nmc=args.mc if args.mc and args.mc > 0 else 1_000_000,
    seed=args.seed,
)

    print("Schmidt exponential-distribution test")
    print("-------------------------------------")
    print(f"input file       : {args.input}")
    print(f"N events         : {n}")
    print(f"times            : {np.array2string(times, precision=6, separator=', ')} {args.unit}")
    print(f"mean lifetime tau: {tau:.6g} {args.unit}")
    print(f"T1/2 = ln2*tau   : {t12:.6g} {args.unit}")
    print(f"sigma_Theta_exp  : {sig:.6f}")

    if lower is None or upper is None:
        print("Schmidt 90% limits: not defined for this N")
    else:
        print(f"90% limits       : [{lower:.6f}, {upper:.6f}]  ({limit_source})")

        if sig < lower:
            print(
                "judgement        : TOO NARROW; below 5% lower tail. "
                "Possible non-decay/periodic source or missing time range."
            )
        elif sig > upper:
            print(
                "judgement        : TOO BROAD; above 5% upper tail. "
                "Possible mixture of lifetimes/background."
            )
        else:
            print(
                "judgement        : compatible with one exponential species "
                "at Schmidt 90% level."
            )

    if args.mc and n >= 2:
        p_low, p_high = mc_pvalues(
            n=n,
            sigma_obs=sig,
            nmc=args.mc,
            seed=args.seed,
        )

        print(f"MC p_low=P(sigma<=obs)  : {p_low:.6g}  ({args.mc} trials)")
        print(f"MC p_high=P(sigma>=obs) : {p_high:.6g}  ({args.mc} trials)")

    if args.plot:
        plot_log_time(
            times=times,
            out=args.out,
            unit=args.unit,
        )

        print(f"plot saved       : {args.out}")


if __name__ == "__main__":
    main()