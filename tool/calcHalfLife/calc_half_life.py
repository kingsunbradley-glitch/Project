#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐事件衰变时间的半衰期计算（Python 3.8+，仅标准库）。

N < 21：Schmidt et al. (1984) 式 (17) 的精确中心置信区间。
N >= 21：无分箱极大似然估计与渐近似然比置信区间。
两分支均取恰好 0.68 C.L.，半衰期为 ln(2) * mean(t_i)。
默认输入为本脚本同目录的 input.dat，每行一个正的衰变时间差。
适用模型：单一核素、独立指数衰变、无背景、无观测时间截断或删失。
"""

import argparse
import math
from pathlib import Path
from statistics import NormalDist


CONFIDENCE = 0.68
MLE_MIN_EVENTS = 21


def read_times(path):
    """读取单列时间，支持空行、# 注释和 Fortran D 指数。"""
    times = []
    with path.open(encoding="utf-8-sig") as stream:
        for number, line in enumerate(stream, 1):
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            if len(fields) != 1:
                raise ValueError(f"{path}:{number}: 每行应只有一个衰变时间差")
            try:
                value = float(fields[0].replace("D", "e").replace("d", "e"))
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: 无效数值 {fields[0]!r}") from exc
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{path}:{number}: 时间差必须为有限正数")
            times.append(value)
    if not times:
        raise ValueError(f"{path}: 没有有效实验数据")
    return times


def _bisect(function, a, b):
    """对已夹住的连续函数零点二分求解。"""
    fa, fb = function(a), function(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if (fa < 0) == (fb < 0):
        raise ValueError("求根区间未夹住零点")
    for _ in range(160):
        mid = (a + b) / 2
        if mid == a or mid == b:
            return mid
        fm = function(mid)
        if fm == 0:
            return mid
        if (fa < 0) == (fm < 0):
            a, fa = mid, fm
        else:
            b = mid
    return (a + b) / 2


def _gamma_survival(n, x):
    """整数形状 Gamma(n, 1) 的上尾概率，供 N<21 的式 (17) 使用。"""
    term = 1.0
    terms = [term]
    for k in range(1, n):
        term *= x / k
        terms.append(term)
    return math.exp(-x) * math.fsum(terms)


def definition_interval(n, confidence=CONFIDENCE):
    """返回式 (17) 的 tau/tau_hat 精确中心置信区间。

    若 S=sum(t_i)，则 X=S/tau ~ Gamma(n, 1)。令 Q_n 为其上尾概率，
    式 (17) 要求 Q_n(S/tau_low)=alpha/2，
    Q_n(S/tau_high)=1-alpha/2，其中 alpha=1-confidence。
    直接求解有限和，不使用式 (18) 的近似，也不积分归一化似然函数。
    """
    if not isinstance(n, int) or not 1 <= n < MLE_MIN_EVENTS:
        raise ValueError("定义法要求整数事件数 1 <= n < 21")
    if not 0 < confidence < 1:
        raise ValueError("要求 0 < confidence < 1")
    tail = (1 - confidence) / 2

    def solve_x(probability):
        upper = float(n)
        while _gamma_survival(n, upper) > probability:
            upper *= 2
        return _bisect(lambda x: _gamma_survival(n, x) - probability, 0.0, upper)

    return n / solve_x(tail), n / solve_x(1 - tail)


def likelihood_interval(n, confidence=CONFIDENCE):
    """返回 tau/tau_hat 的似然比区间，使用一参数渐近 C.L. 标定。

    无背景、无截断、独立指数衰变：tau_hat = mean(t_i)。
    令 y = log(tau/tau_hat)，则
        log L_max - log L(tau) = n * (y + exp(-y) - 1)。
    阈值为 NormalDist().inv_cdf((1+C.L.)/2)**2 / 2。
    对恰好 0.68 使用其对应阈值，不能将 68.2689% 直接当作 68%。
    这是渐近似然比区间，不宣称有限样本覆盖率恰为 0.68。
    """
    if not isinstance(n, int) or n < 1 or not 0 < confidence < 1:
        raise ValueError("要求 n >= 1 且 0 < confidence < 1")
    threshold = NormalDist().inv_cdf((1 + confidence) / 2) ** 2 / 2

    def equation(y):
        return n * (y + math.expm1(-y)) - threshold

    left, right = -1.0, 1.0
    while equation(left) < 0:
        left *= 2
    while equation(right) < 0:
        right *= 2
    return math.exp(_bisect(equation, left, 0.0)), math.exp(_bisect(equation, 0.0, right))


def calculate(times):
    """返回事件数、平均寿命、半衰期及半衰期区间。"""
    if not times or any(not math.isfinite(t) or t <= 0 for t in times):
        raise ValueError("实验数据必须为非空的有限正时间差序列")
    n = len(times)
    # 先按最大值缩放，避免 sum(t_i) 溢出或逐项除以 n 后下溢。
    scale = max(times)
    tau = scale * (math.fsum(t / scale for t in times) / n)
    half_life = math.log(2) * tau
    interval = definition_interval if n < MLE_MIN_EVENTS else likelihood_interval
    low_factor, high_factor = interval(n)
    low, high = half_life * low_factor, half_life * high_factor
    if not (0 < low <= half_life <= high < math.inf):
        raise ValueError("数值超出浮点表示范围，请变换输入时间单位")
    return n, tau, half_life, low, high


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(__file__).resolve().with_name("input.dat"),
        help="输入文件，默认读取脚本同目录 input.dat",
    )
    parser.add_argument("--unit", default="输入时间单位", help="时间单位标签，如 ms 或 s")
    args = parser.parse_args()
    try:
        n, tau, half_life, low, high = calculate(read_times(args.input))
    except (OSError, ValueError, OverflowError) as exc:
        parser.exit(1, f"错误：{exc}\n")
    print(f"输入文件：{args.input.resolve()}")
    print(f"事件数 N = {n}")
    if n < MLE_MIN_EVENTS:
        print("方法：论文定义法；式 (17) 精确中心置信区间")
    else:
        print("方法：无分箱极大似然估计；渐近似然比置信区间")
    print(f"C.L. = {CONFIDENCE:.2f}")
    print(f"平均寿命 tau = {tau:.10g} {args.unit}")
    print(f"T_1/2 = {half_life:.10g} +{high-half_life:.10g} -{half_life-low:.10g} {args.unit}")
    print(f"区间 = [{low:.10g}, {high:.10g}] {args.unit}")


if __name__ == "__main__":
    main()
