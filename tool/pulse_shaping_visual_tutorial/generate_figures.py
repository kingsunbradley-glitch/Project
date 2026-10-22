#!/usr/bin/env python3
"""生成梯形、三角形与 RC-CR（RCCR）数字脉冲成形教学图。

本脚本只依赖 NumPy 和 Matplotlib。图中的输入是电荷灵敏前放的指数衰减脉冲；
先做理想极零相消，将它还原为阶跃，再分别送入三种成形器。
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

# 某些计算节点的用户目录不可写；把 Matplotlib 字体缓存放到临时目录。
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pulse_shaping_mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


# Okabe-Ito 色盲友好配色
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"
GRAY = "#5F6368"
LIGHT_GRAY = "#E8EAED"
BLACK = "#202124"


def configure_style() -> None:
    """设置适合中文教学图的统一样式。"""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial Unicode MS",
                "Hiragino Sans GB",
                "Heiti SC",
                "DejaVu Sans",
            ],
            "axes.unicode_minus": False,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.linewidth": 0.9,
            "lines.linewidth": 2.2,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def preamplifier_pulse(t: np.ndarray, t0: float, tau_decay: float) -> np.ndarray:
    """理想电荷灵敏前放脉冲：瞬时上升，随后指数衰减。"""
    out = np.zeros_like(t)
    mask = t >= t0
    out[mask] = np.exp(-(t[mask] - t0) / tau_decay)
    return out


def pole_zero_to_step(signal: np.ndarray, dt: float, tau_decay: float) -> np.ndarray:
    """用已知衰减常数做极零相消，并积分得到等幅阶跃。"""
    pole = np.exp(-dt / tau_decay)
    impulse = np.empty_like(signal)
    impulse[0] = signal[0]
    impulse[1:] = signal[1:] - pole * signal[:-1]
    return np.cumsum(impulse)


def moving_average(signal: np.ndarray, length: int) -> np.ndarray:
    """因果长度-L 滑动平均。"""
    kernel = np.ones(length, dtype=float) / length
    return np.convolve(signal, kernel, mode="full")[: signal.size]


def delay(signal: np.ndarray, samples: int) -> np.ndarray:
    """将序列延迟指定采样点，前部补零。"""
    if samples == 0:
        return signal.copy()
    out = np.zeros_like(signal)
    out[samples:] = signal[:-samples]
    return out


def trapezoid_shaper(step: np.ndarray, rise: int, flat: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """两个相隔 rise+flat 的滑动平均之差。flat=0 时退化为三角成形。"""
    first = moving_average(step, rise)
    second = delay(first, rise + flat)
    shaped = first - second
    return shaped, first, second


def cr_high_pass(signal: np.ndarray, dt: float, tau: float) -> np.ndarray:
    """一阶 CR 高通（后向欧拉离散）。"""
    alpha = tau / (tau + dt)
    out = np.zeros_like(signal)
    for i in range(1, signal.size):
        out[i] = alpha * (out[i - 1] + signal[i] - signal[i - 1])
    return out


def rc_low_pass(signal: np.ndarray, dt: float, tau: float) -> np.ndarray:
    """一阶 RC 低通（后向欧拉离散）。"""
    beta = dt / (tau + dt)
    out = np.zeros_like(signal)
    for i in range(1, signal.size):
        out[i] = out[i - 1] + beta * (signal[i] - out[i - 1])
    return out


def rccr_shaper(step: np.ndarray, dt: float, tau_cr: float, tau_rc: float) -> tuple[np.ndarray, np.ndarray]:
    """CR 微分 + RC 积分级联；理想线性系统中两级的次序可交换。"""
    differentiated = cr_high_pass(step, dt, tau_cr)
    shaped = rc_low_pass(differentiated, dt, tau_rc)
    return shaped, differentiated


def normalize(signal: np.ndarray) -> np.ndarray:
    peak = np.max(np.abs(signal))
    return signal / peak if peak else signal.copy()


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, text: str) -> None:
    ax.text(-0.10, 1.04, text, transform=ax.transAxes, weight="bold", fontsize=14, va="bottom")


def save_figure(fig: plt.Figure, output: Path, dpi: int) -> None:
    fig.savefig(output.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def draw_box(ax: plt.Axes, xy: tuple[float, float], size: tuple[float, float], text: str, color: str) -> None:
    x, y = xy
    w, h = size
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.025,rounding_size=0.03",
        facecolor=color,
        edgecolor="none",
        alpha=0.16,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=BLACK, weight="bold")


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, lw=1.5, color=GRAY))


def make_overview(
    t: np.ndarray,
    preamp: np.ndarray,
    step: np.ndarray,
    trap: np.ndarray,
    triangle: np.ndarray,
    rccr: np.ndarray,
    output_dir: Path,
    dpi: int,
) -> None:
    fig = plt.figure(figsize=(12.2, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 1.6])

    flow = fig.add_subplot(gs[0, :])
    flow.set_xlim(0, 1)
    flow.set_ylim(0, 1)
    flow.axis("off")
    draw_box(flow, (0.02, 0.28), (0.17, 0.42), "前放指数脉冲\n（能量 → 幅度）", SKY)
    draw_arrow(flow, (0.20, 0.49), (0.28, 0.49))
    draw_box(flow, (0.29, 0.28), (0.17, 0.42), "极零相消\n还原为阶跃", ORANGE)
    draw_arrow(flow, (0.47, 0.49), (0.55, 0.49))
    draw_box(flow, (0.56, 0.61), (0.17, 0.27), "梯形成形\n平顶读能量", BLUE)
    draw_box(flow, (0.56, 0.365), (0.17, 0.27), "三角成形\n尖顶读峰值", GREEN)
    draw_box(flow, (0.56, 0.12), (0.17, 0.27), "RC–CR 成形\n圆滑单峰", VERMILION)
    draw_arrow(flow, (0.47, 0.49), (0.55, 0.745))
    draw_arrow(flow, (0.47, 0.49), (0.55, 0.50))
    draw_arrow(flow, (0.47, 0.49), (0.55, 0.255))
    draw_arrow(flow, (0.74, 0.745), (0.82, 0.745))
    draw_arrow(flow, (0.74, 0.50), (0.82, 0.50))
    draw_arrow(flow, (0.74, 0.255), (0.82, 0.255))
    draw_box(flow, (0.83, 0.28), (0.15, 0.42), "提取幅度\n标定为能量", PURPLE)
    flow.set_title("从探测器前放波形到能量：三条成形路线", fontsize=16, weight="bold", pad=6)

    ax1 = fig.add_subplot(gs[1, 0])
    ax1.plot(t, preamp, color=SKY, label="前放输出：快升 + 指数衰减")
    ax1.plot(t, step, color=ORANGE, linestyle="--", label="极零相消后的阶跃")
    ax1.axvline(20, color=GRAY, lw=1, ls=":")
    ax1.annotate("粒子到达", xy=(20, 1.0), xytext=(28, 0.75), arrowprops=dict(arrowstyle="->", color=GRAY))
    ax1.set(xlabel="时间 (μs)", ylabel="归一化幅度", xlim=(5, 105), ylim=(-0.05, 1.15))
    ax1.legend(frameon=False, loc="center right")
    ax1.set_title("先把长尾巴校正成阶跃")
    clean_axes(ax1)
    panel_label(ax1, "A")

    ax2 = fig.add_subplot(gs[1, 1])
    ax2.plot(t, normalize(trap), color=BLUE, label="梯形：有平顶")
    ax2.plot(t, normalize(triangle), color=GREEN, ls="--", label="三角：平顶宽度 G = 0")
    ax2.plot(t, normalize(rccr), color=VERMILION, ls="-.", label="RC–CR：指数圆滑单峰")
    ax2.axhline(0, color=BLACK, lw=0.8)
    ax2.set(xlabel="时间 (μs)", ylabel="归一化成形输出", xlim=(12, 90), ylim=(-0.06, 1.12))
    ax2.legend(frameon=False)
    ax2.set_title("同一个阶跃，不同滤波器“雕刻”出不同形状")
    clean_axes(ax2)
    panel_label(ax2, "B")

    save_figure(fig, output_dir / "00_总览_三种成形路线", dpi)


def make_trapezoid_triangle(
    t: np.ndarray,
    first: np.ndarray,
    second: np.ndarray,
    trap: np.ndarray,
    tri_first: np.ndarray,
    tri_second: np.ndarray,
    triangle: np.ndarray,
    rise_samples: int,
    flat_samples: int,
    dt: float,
    output_dir: Path,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10.8, 10.0), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.plot(t, first, color=BLUE, label=r"当前滑动平均 $M_L[n]$")
    ax.plot(t, second, color=ORANGE, ls="--", label=r"延迟副本 $M_L[n-(L+G)]$")
    ax.fill_between(t, second, first, where=first >= second, color=BLUE, alpha=0.16, label="两者之差")
    ax.set(ylabel="滑动平均", ylim=(-0.05, 1.12))
    ax.set_title("梯形成形的本质：两个矩形平均窗相减")
    ax.legend(frameon=False, loc="lower right")
    clean_axes(ax)
    panel_label(ax, "A")

    ax = axes[1]
    ax.plot(t, trap, color=BLUE)
    t0 = 20.0
    t1 = t0 + rise_samples * dt
    t2 = t1 + flat_samples * dt
    t3 = t2 + rise_samples * dt
    ax.axvspan(t0, t1, color=SKY, alpha=0.12)
    ax.axvspan(t1, t2, color=ORANGE, alpha=0.16)
    ax.axvspan(t2, t3, color=SKY, alpha=0.12)
    ax.annotate("上升时间 L", xy=((t0 + t1) / 2, 0.45), ha="center", color=BLUE, weight="bold")
    ax.annotate("平顶 G\n这里读能量更稳", xy=((t1 + t2) / 2, 1.0), xytext=((t1 + t2) / 2, 0.55),
                ha="center", arrowprops=dict(arrowstyle="->", color=ORANGE), color=BLACK)
    ax.annotate("下降时间 L", xy=((t2 + t3) / 2, 0.45), ha="center", color=BLUE, weight="bold")
    ax.set(ylabel="梯形输出", ylim=(-0.06, 1.15))
    ax.set_title("平顶给电荷收集时间留出余量，减小弹道亏损")
    clean_axes(ax)
    panel_label(ax, "B")

    ax = axes[2]
    ax.plot(t, tri_first, color=GREEN, alpha=0.45, label="当前平均")
    ax.plot(t, tri_second, color=ORANGE, alpha=0.45, ls="--", label="恰好延迟 L 的平均")
    ax.plot(t, triangle, color=GREEN, lw=3.0, label="相减结果：三角形")
    peak_index = int(np.argmax(triangle))
    ax.scatter([t[peak_index]], [triangle[peak_index]], s=55, color=VERMILION, zorder=4)
    ax.annotate("G = 0，两段斜坡直接相接\n峰值就是幅度估计",
                xy=(t[peak_index], triangle[peak_index]), xytext=(t[peak_index] + 13, 0.72),
                arrowprops=dict(arrowstyle="->", color=VERMILION), ha="left")
    ax.set(xlabel="时间 (μs)", ylabel="三角输出", xlim=(12, 82), ylim=(-0.06, 1.15))
    ax.set_title("三角算法就是“平顶宽度为零”的梯形算法")
    ax.legend(frameon=False, loc="upper right")
    clean_axes(ax)
    panel_label(ax, "C")

    save_figure(fig, output_dir / "01_梯形与三角_窗口相减", dpi)


def make_rccr(
    t: np.ndarray,
    step: np.ndarray,
    differentiated: np.ndarray,
    rccr: np.ndarray,
    tau: float,
    output_dir: Path,
    dpi: int,
) -> None:
    fig = plt.figure(figsize=(11.2, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.58, 1, 1])

    flow = fig.add_subplot(gs[0])
    flow.set_xlim(0, 1)
    flow.set_ylim(0, 1)
    flow.axis("off")
    draw_box(flow, (0.03, 0.27), (0.20, 0.48), "阶跃输入", ORANGE)
    draw_arrow(flow, (0.24, 0.51), (0.34, 0.51))
    draw_box(flow, (0.35, 0.27), (0.20, 0.48), "CR 高通\n抑制直流，突出变化", VERMILION)
    draw_arrow(flow, (0.56, 0.51), (0.66, 0.51))
    draw_box(flow, (0.67, 0.27), (0.20, 0.48), "RC 低通\n压制高频，变圆滑", BLUE)
    draw_arrow(flow, (0.88, 0.51), (0.97, 0.51))
    flow.text(0.50, 0.04, "RCCR 常写作 RC–CR；理想线性级联中两级次序不影响总传递函数",
              ha="center", va="bottom", color=GRAY)
    flow.set_title("RC–CR：一边去掉慢基线，一边压低快噪声", fontsize=16, weight="bold")

    ax1 = fig.add_subplot(gs[1])
    ax1.plot(t, step, color=ORANGE, label="阶跃输入")
    ax1.plot(t, normalize(differentiated), color=VERMILION, ls="--", label="经过 CR：指数衰减")
    ax1.set(ylabel="归一化幅度", xlim=(12, 85), ylim=(-0.05, 1.12))
    ax1.legend(frameon=False)
    ax1.set_title("CR 像“变化探测器”：恒定部分回到零，只保留边沿后的短脉冲")
    clean_axes(ax1)
    panel_label(ax1, "A")

    ax2 = fig.add_subplot(gs[2])
    ax2.plot(t, normalize(differentiated), color=VERMILION, alpha=0.42, ls="--", label="CR 输出")
    ax2.plot(t, normalize(rccr), color=BLUE, lw=3.0, label="再经过 RC：半高斯单峰")
    peak = int(np.argmax(rccr))
    ax2.scatter([t[peak]], [normalize(rccr)[peak]], s=60, color=BLACK, zorder=4)
    ax2.annotate(
        f"τCR = τRC = {tau:.0f} μs\n峰值约出现在阶跃后 1τ",
        xy=(t[peak], 1.0),
        xytext=(t[peak] + 16, 0.72),
        arrowprops=dict(arrowstyle="->", color=BLACK),
    )
    ax2.text(66, 0.32, r"连续时间形状：$y(t)\propto (t/\tau)e^{-t/\tau}$", color=BLUE, ha="center")
    ax2.set(xlabel="时间 (μs)", ylabel="归一化幅度", xlim=(12, 85), ylim=(-0.05, 1.12))
    ax2.legend(frameon=False)
    ax2.set_title("RC 把尖锐指数脉冲摊开：峰越宽，噪声越小，但堆积概率越高")
    clean_axes(ax2)
    panel_label(ax2, "B")

    save_figure(fig, output_dir / "02_RCCR_RC-CR_两级作用", dpi)


def make_parameter_effects(step: np.ndarray, t: np.ndarray, dt: float, output_dir: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.6), constrained_layout=True)

    ax = axes[0]
    for flat, color, style in [(0, GREEN, "--"), (8, BLUE, "-"), (18, PURPLE, "-.")]:
        shaped, _, _ = trapezoid_shaper(step, rise=12, flat=flat)
        ax.plot(t, normalize(shaped), color=color, ls=style, label=f"G = {flat * dt:g} μs")
    ax.set(xlabel="时间 (μs)", ylabel="归一化输出", xlim=(14, 82), ylim=(-0.05, 1.12))
    ax.set_title("梯形：G 控制平顶")
    ax.legend(frameon=False)
    clean_axes(ax)
    panel_label(ax, "A")

    ax = axes[1]
    for rise, color, style in [(7, ORANGE, ":"), (12, GREEN, "--"), (20, BLUE, "-")]:
        shaped, _, _ = trapezoid_shaper(step, rise=rise, flat=0)
        ax.plot(t, normalize(shaped), color=color, ls=style, label=f"L = {rise * dt:g} μs")
    ax.set(xlabel="时间 (μs)", ylabel="归一化输出", xlim=(14, 82), ylim=(-0.05, 1.12))
    ax.set_title("三角：L 控制宽度")
    ax.legend(frameon=False)
    clean_axes(ax)
    panel_label(ax, "B")

    ax = axes[2]
    for tau, color, style in [(5.0, ORANGE, ":"), (9.0, VERMILION, "--"), (15.0, BLUE, "-")]:
        shaped, _ = rccr_shaper(step, dt, tau, tau)
        ax.plot(t, normalize(shaped), color=color, ls=style, label=f"τ = {tau:.0f} μs")
    ax.set(xlabel="时间 (μs)", ylabel="归一化输出", xlim=(14, 92), ylim=(-0.05, 1.12))
    ax.set_title("RC–CR：τ 控制圆滑程度")
    ax.legend(frameon=False)
    clean_axes(ax)
    panel_label(ax, "C")

    fig.suptitle("参数不是越大越好：降噪能力与脉冲堆积之间要折中", fontsize=16, weight="bold")
    save_figure(fig, output_dir / "03_参数如何改变波形", dpi)


def numerical_checks(
    trap: np.ndarray,
    triangle: np.ndarray,
    rccr: np.ndarray,
    t: np.ndarray,
    rise_samples: int,
    flat_samples: int,
    dt: float,
    tau: float,
) -> None:
    """对教学图依赖的关键性质做轻量数值检查。"""
    start_index = int(round(20.0 / dt))
    plateau = trap[start_index + rise_samples : start_index + rise_samples + flat_samples]
    assert plateau.size == flat_samples
    assert np.ptp(plateau) < 1e-12, "梯形平顶不平"
    assert np.isclose(np.max(triangle), 1.0, atol=1e-12), "三角峰值异常"
    rccr_peak_time = t[int(np.argmax(rccr))] - 20.0
    assert abs(rccr_peak_time - tau) < 2.5 * dt, "RC-CR 峰位偏离约 1τ 的预期"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("figures"), help="图片输出目录")
    parser.add_argument("--dpi", type=int, default=220, help="PNG 分辨率")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    configure_style()
    dt = 0.25  # μs
    t = np.arange(0.0, 120.0, dt)
    t0 = 20.0
    tau_decay = 45.0
    rise_samples = 56  # 14 μs
    flat_samples = 32  # 8 μs
    tau_rccr = 9.0

    preamp = preamplifier_pulse(t, t0, tau_decay)
    step = pole_zero_to_step(preamp, dt, tau_decay)
    trap, first, second = trapezoid_shaper(step, rise_samples, flat_samples)
    triangle, tri_first, tri_second = trapezoid_shaper(step, rise_samples, 0)
    rccr, differentiated = rccr_shaper(step, dt, tau_rccr, tau_rccr)

    numerical_checks(trap, triangle, rccr, t, rise_samples, flat_samples, dt, tau_rccr)
    make_overview(t, preamp, step, trap, triangle, rccr, args.output_dir, args.dpi)
    make_trapezoid_triangle(
        t,
        first,
        second,
        trap,
        tri_first,
        tri_second,
        triangle,
        rise_samples,
        flat_samples,
        dt,
        args.output_dir,
        args.dpi,
    )
    make_rccr(t, step, differentiated, rccr, tau_rccr, args.output_dir, args.dpi)
    make_parameter_effects(step, t, dt, args.output_dir, args.dpi)
    print(f"已生成 4 组 PNG + SVG 教学图：{args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
