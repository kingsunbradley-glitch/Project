#!/usr/bin/env python3
"""生成梯形、三角与 RC-CR 脉冲成形的中文 GIF 动画。"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pulse_shaping_mplconfig"))

import matplotlib

matplotlib.use("Agg")
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np

from generate_figures import (
    BLACK,
    BLUE,
    GRAY,
    GREEN,
    LIGHT_GRAY,
    ORANGE,
    SKY,
    VERMILION,
    clean_axes,
    configure_style,
    normalize,
    pole_zero_to_step,
    preamplifier_pulse,
    rccr_shaper,
    trapezoid_shaper,
)


def make_demo_signals() -> dict[str, object]:
    dt = 0.25  # μs
    t = np.arange(0.0, 120.0, dt)
    preamp = preamplifier_pulse(t, t0=20.0, tau_decay=45.0)
    step = pole_zero_to_step(preamp, dt=dt, tau_decay=45.0)
    rise = 56  # 14 μs
    flat = 32  # 8 μs
    trap, trap_now, trap_delayed = trapezoid_shaper(step, rise=rise, flat=flat)
    triangle, tri_now, tri_delayed = trapezoid_shaper(step, rise=rise, flat=0)
    rccr, cr = rccr_shaper(step, dt=dt, tau_cr=9.0, tau_rc=9.0)
    return {
        "dt": dt,
        "t": t,
        "step": step,
        "rise": rise,
        "flat": flat,
        "trap": trap,
        "trap_now": trap_now,
        "trap_delayed": trap_delayed,
        "triangle": triangle,
        "tri_now": tri_now,
        "tri_delayed": tri_delayed,
        "rccr": rccr,
        "cr": cr,
    }


def set_dot(dot: plt.Line2D, x: float, y: float) -> None:
    dot.set_data([x], [y])


def fast_rise_decay_pulse(
    t: np.ndarray,
    t0: float,
    tau_rise: float,
    tau_decay: float,
) -> np.ndarray:
    """有限快上升、随后慢指数衰减的归一化脉冲。"""
    pulse = np.zeros_like(t)
    mask = t >= t0
    elapsed = t[mask] - t0
    pulse[mask] = (1.0 - np.exp(-elapsed / tau_rise)) * np.exp(-elapsed / tau_decay)
    peak = np.max(pulse)
    return pulse / peak if peak else pulse


def make_three_input_signals() -> dict[str, object]:
    """构造正常脉冲，以及叠加等振幅高、低频正弦噪声的两个版本。"""
    dt = 0.05  # μs，足够解析高频正弦
    t = np.arange(0.0, 90.0 + dt, dt)
    t0 = 18.0
    tau_rise = 0.8
    tau_decay = 22.0
    noise_amplitude = 0.08
    high_frequency = 0.80  # MHz；横轴单位为 μs，所以也是 cycle/μs
    low_frequency = 0.04   # MHz
    normal = fast_rise_decay_pulse(t, t0=t0, tau_rise=tau_rise, tau_decay=tau_decay)
    high_noise = noise_amplitude * np.sin(2.0 * np.pi * high_frequency * t)
    low_noise = noise_amplitude * np.sin(2.0 * np.pi * low_frequency * t)
    return {
        "dt": dt,
        "t": t,
        "t0": t0,
        "tau_rise": tau_rise,
        "tau_decay": tau_decay,
        "noise_amplitude": noise_amplitude,
        "high_frequency": high_frequency,
        "low_frequency": low_frequency,
        "signals": [normal, normal + high_noise, normal + low_noise],
    }


def save_animation(animation: FuncAnimation, output: Path, fps: int, dpi: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    animation.save(output, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(animation._fig)


def make_input_signal_animation(output_dir: Path, fps: int, dpi: int) -> None:
    """演示正常脉冲及分别叠加高、低频正弦噪声后的输入。"""
    signal_data = make_three_input_signals()
    t = np.asarray(signal_data["t"])
    t0 = float(signal_data["t0"])
    tau_rise = float(signal_data["tau_rise"])
    tau_decay = float(signal_data["tau_decay"])
    noise_amplitude = float(signal_data["noise_amplitude"])
    high_frequency = float(signal_data["high_frequency"])
    low_frequency = float(signal_data["low_frequency"])
    signals = [np.asarray(signal) for signal in signal_data["signals"]]
    colors = [BLACK, VERMILION, BLUE]
    titles = [
        "① 正常脉冲：快速上升，随后逐渐衰减",
        "② 正常脉冲 + 高频正弦噪声：短时间内快速抖动",
        "③ 正常脉冲 + 低频正弦噪声：基线缓慢起伏",
    ]
    formulas = [
        r"$p(t)$",
        rf"$p(t)+{noise_amplitude:.2f}\sin(2\pi\times {high_frequency:.2f}\,t)$",
        rf"$p(t)+{noise_amplitude:.2f}\sin(2\pi\times {low_frequency:.2f}\,t)$",
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10.2, 7.6), sharex=True, constrained_layout=True)
    fig.suptitle("三种输入信号：脉冲相同，只改变叠加的正弦噪声", fontsize=16, weight="bold")
    animated: list[tuple[plt.Line2D, plt.Line2D, plt.Line2D, np.ndarray]] = []

    for ax, signal, color, title, formula in zip(axes, signals, colors, titles, formulas):
        ax.plot(t, signal, color=LIGHT_GRAY, lw=2.0, label="完整参考波形")
        reveal, = ax.plot([], [], color=color, lw=2.5, label="当前已采集信号")
        dot, = ax.plot([], [], "o", color=ORANGE, ms=6, zorder=5)
        cursor = ax.axvline(0, color=GRAY, ls=":", lw=1.2)
        ax.axhline(0, color=BLACK, lw=0.8)
        ax.axvline(t0, color=GRAY, lw=1.0, ls="--")
        ax.text(t0 + 0.8, -0.12, "脉冲开始", color=GRAY, va="bottom", fontsize=9)
        ax.text(
            0.985,
            0.82,
            formula,
            transform=ax.transAxes,
            ha="right",
            va="top",
            color=color,
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec=LIGHT_GRAY, alpha=0.92),
        )
        ax.set(xlim=(0, 90), ylim=(-0.16, 1.16), ylabel="归一化幅度")
        ax.set_title(title, loc="left", weight="bold")
        clean_axes(ax)
        animated.append((reveal, dot, cursor, signal))

    axes[0].text(
        0.985,
        0.56,
        rf"上升常数 $\tau_r={tau_rise:g}$ μs；衰减常数 $\tau_d={tau_decay:g}$ μs",
        transform=axes[0].transAxes,
        ha="right",
        color=GRAY,
    )
    axes[-1].set_xlabel("时间 (μs)")
    axes[-1].text(
        0.985,
        0.40,
        "两种噪声幅度相同：A = 0.08\n所以差别来自频率，而不是振幅",
        transform=axes[-1].transAxes,
        ha="right",
        color=GRAY,
    )
    clock = fig.text(0.87, 0.955, "", ha="right", va="top", color=GRAY, weight="bold")
    frame_indices = np.linspace(0, t.size - 1, 120).astype(int)

    def update(frame: int) -> list[object]:
        i = int(frame_indices[frame])
        current_time = t[i]
        clock.set_text(f"采样进行到：t = {current_time:4.1f} μs")
        artists: list[object] = [clock]
        for reveal, dot, cursor, signal in animated:
            reveal.set_data(t[: i + 1], signal[: i + 1])
            set_dot(dot, current_time, signal[i])
            cursor.set_xdata([current_time, current_time])
            artists.extend([reveal, dot, cursor])
        return artists

    animation = FuncAnimation(fig, update, frames=len(frame_indices), interval=1000 / fps, blit=False, repeat=True)
    save_animation(animation, output_dir / "00_三种输入信号演示.gif", fps, dpi)


def make_three_by_three_shaping_animation(output_dir: Path, fps: int, dpi: int) -> None:
    """把三种输入分别送入梯形、三角和 RC-CR，生成动态处理矩阵。"""
    signal_data = make_three_input_signals()
    t = np.asarray(signal_data["t"])
    dt = float(signal_data["dt"])
    t0 = float(signal_data["t0"])
    tau_decay = float(signal_data["tau_decay"])
    inputs = [np.asarray(signal) for signal in signal_data["signals"]]

    rise = int(round(8.0 / dt))
    flat = int(round(4.0 / dt))
    tau_rccr = 5.0
    shaped_by_input: list[list[np.ndarray]] = []
    for raw in inputs:
        corrected = pole_zero_to_step(raw, dt=dt, tau_decay=tau_decay)
        corrected -= np.mean(corrected[t < t0 - 2.0])
        trap, _, _ = trapezoid_shaper(corrected, rise=rise, flat=flat)
        triangle, _, _ = trapezoid_shaper(corrected, rise=rise, flat=0)
        rccr, _ = rccr_shaper(corrected, dt=dt, tau_cr=tau_rccr, tau_rc=tau_rccr)
        shaped_by_input.append([trap, triangle, rccr])

    # 三列必须共用正常信号的标尺，否则单独归一化会掩盖噪声残留与幅度偏差。
    clean_scales = [np.max(np.abs(values)) for values in shaped_by_input[0]]
    normalized_outputs = [
        [values / scale for values, scale in zip(outputs, clean_scales)]
        for outputs in shaped_by_input
    ]

    column_titles = ["正常脉冲", "+ 高频正弦噪声", "+ 低频正弦噪声"]
    row_titles = ["输入", "梯形\nL=8 μs\nG=4 μs", "三角\nL=8 μs", "RCCR\nτ=5 μs"]
    row_colors = [BLACK, BLUE, GREEN, VERMILION]
    matrix = [
        inputs,
        [normalized_outputs[col][0] for col in range(3)],
        [normalized_outputs[col][1] for col in range(3)],
        [normalized_outputs[col][2] for col in range(3)],
    ]

    fig, axes = plt.subplots(4, 3, figsize=(12.8, 9.2), sharex=True, constrained_layout=True)
    fig.suptitle(
        "三种输入 × 三种成形算法：逐样本处理结果\n"
        "成形前统一用 τd=22 μs 做极零相消；三列输出均按正常信号峰值归一化",
        fontsize=14,
        weight="bold",
    )
    animated: list[tuple[plt.Line2D, plt.Line2D, plt.Line2D, np.ndarray]] = []
    baseline_mask = (t >= 8.0) & (t <= 16.0)

    for row in range(4):
        for col in range(3):
            ax = axes[row, col]
            values = np.asarray(matrix[row][col])
            color = row_colors[row]
            ax.plot(t, values, color=LIGHT_GRAY, lw=1.8)
            reveal, = ax.plot([], [], color=color, lw=2.2)
            dot, = ax.plot([], [], "o", color=ORANGE, ms=4.8, zorder=5)
            cursor = ax.axvline(0, color=GRAY, ls=":", lw=1.0)
            ax.axhline(0, color=BLACK, lw=0.7)
            ax.axvline(t0, color=GRAY, lw=0.8, ls="--")
            ax.set_xlim(0, 90)
            ax.set_ylim((-0.18, 1.18) if row == 0 else (-0.22, 1.22))
            clean_axes(ax)
            if row == 0:
                ax.set_title(column_titles[col], weight="bold")
            if col == 0:
                ax.set_ylabel(row_titles[row])
            else:
                ax.tick_params(labelleft=False)
            if row == 3:
                ax.set_xlabel("时间 (μs)")
            if row > 0:
                rms = float(np.sqrt(np.mean(values[baseline_mask] ** 2)))
                ax.text(
                    0.97,
                    0.88,
                    f"成形后基线 RMS={rms:.3f}",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=8,
                    color=GRAY,
                    bbox=dict(boxstyle="round,pad=0.20", fc="white", ec=LIGHT_GRAY, alpha=0.90),
                )
            animated.append((reveal, dot, cursor, values))

    clock = fig.text(0.975, 0.985, "", ha="right", va="top", color=GRAY, weight="bold")
    frame_indices = np.linspace(0, t.size - 1, 120).astype(int)

    def update(frame: int) -> list[object]:
        i = int(frame_indices[frame])
        current_time = t[i]
        clock.set_text(f"t = {current_time:4.1f} μs")
        artists: list[object] = [clock]
        for reveal, dot, cursor, values in animated:
            reveal.set_data(t[: i + 1], values[: i + 1])
            set_dot(dot, current_time, values[i])
            cursor.set_xdata([current_time, current_time])
            artists.extend([reveal, dot, cursor])
        return artists

    animation = FuncAnimation(fig, update, frames=len(frame_indices), interval=1000 / fps, blit=False, repeat=True)
    save_animation(animation, output_dir / "03_三种输入经过三种成形算法.gif", fps, dpi)


def make_window_animation(data: dict[str, object], output_dir: Path, fps: int, dpi: int) -> None:
    """并排展示两个移动窗如何实时形成梯形和三角形。"""
    t = np.asarray(data["t"])
    step = np.asarray(data["step"])
    dt = float(data["dt"])
    rise = int(data["rise"])
    flat = int(data["flat"])
    trap = np.asarray(data["trap"])
    triangle = np.asarray(data["triangle"])
    trap_now = np.asarray(data["trap_now"])
    trap_delayed = np.asarray(data["trap_delayed"])
    tri_now = np.asarray(data["tri_now"])
    tri_delayed = np.asarray(data["tri_delayed"])

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), sharex=True, constrained_layout=True)
    fig.suptitle("移动平均窗边走边算：梯形与三角波形是怎样长出来的", fontsize=16, weight="bold")

    configurations = [
        (axes[0, 0], axes[1, 0], "梯形：延迟 L+G，留下平顶", flat, trap, trap_now, trap_delayed, BLUE),
        (axes[0, 1], axes[1, 1], "三角：G=0，两段斜坡相接", 0, triangle, tri_now, tri_delayed, GREEN),
    ]
    animated: list[dict[str, object]] = []

    for top, bottom, title, gap, output, now, delayed, color in configurations:
        top.plot(t, step, color=BLACK, lw=2.0, label="校正后的阶跃")
        top.set(xlim=(8, 72), ylim=(-0.08, 1.27), ylabel="输入 / 窗口")
        top.set_title(title, weight="bold")
        clean_axes(top)

        current_rect = Rectangle((8, -0.03), rise * dt, 1.12, facecolor=SKY, edgecolor=BLUE, alpha=0.25, lw=1.4)
        delayed_rect = Rectangle((8, -0.03), rise * dt, 1.12, facecolor=ORANGE, edgecolor=ORANGE, alpha=0.20, lw=1.4)
        top.add_patch(current_rect)
        top.add_patch(delayed_rect)
        cursor_top = top.axvline(8, color=GRAY, ls=":", lw=1.2)
        values = top.text(
            0.02,
            0.90,
            "",
            transform=top.transAxes,
            va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=LIGHT_GRAY, alpha=0.92),
        )
        top.text(0.98, 0.13, "蓝窗：当前平均", transform=top.transAxes, ha="right", color=BLUE)
        top.text(0.98, 0.04, "橙窗：延迟平均", transform=top.transAxes, ha="right", color=ORANGE)

        bottom.plot(t, output, color=LIGHT_GRAY, lw=2.0, label="尚未计算")
        reveal, = bottom.plot([], [], color=color, lw=3.0, label="已实时算出的输出")
        dot, = bottom.plot([], [], "o", color=VERMILION, ms=7, zorder=5)
        cursor_bottom = bottom.axvline(8, color=GRAY, ls=":", lw=1.2)
        bottom.axhline(0, color=BLACK, lw=0.8)
        bottom.set(xlim=(8, 72), ylim=(-0.08, 1.16), xlabel="时间 (μs)", ylabel="成形输出")
        bottom.legend(frameon=False, loc="upper right")
        clean_axes(bottom)

        animated.append(
            {
                "gap": gap,
                "output": output,
                "now": now,
                "delayed": delayed,
                "current_rect": current_rect,
                "delayed_rect": delayed_rect,
                "cursor_top": cursor_top,
                "cursor_bottom": cursor_bottom,
                "values": values,
                "reveal": reveal,
                "dot": dot,
            }
        )

    frame_indices = np.linspace(np.searchsorted(t, 10.0), np.searchsorted(t, 68.0), 100).astype(int)

    def update(frame: int) -> list[object]:
        i = int(frame_indices[frame])
        current_time = t[i]
        artists: list[object] = []
        for item in animated:
            gap = int(item["gap"])
            current_start = current_time - (rise - 1) * dt
            delayed_end = current_time - (rise + gap) * dt
            delayed_start = delayed_end - (rise - 1) * dt
            item["current_rect"].set_x(current_start)
            item["delayed_rect"].set_x(delayed_start)
            item["cursor_top"].set_xdata([current_time, current_time])
            item["cursor_bottom"].set_xdata([current_time, current_time])

            now_value = item["now"][i]
            delayed_value = item["delayed"][i]
            output_value = item["output"][i]
            item["values"].set_text(
                f"t = {current_time:4.1f} μs\n"
                f"当前平均 = {now_value:4.2f}\n"
                f"延迟平均 = {delayed_value:4.2f}\n"
                f"两者之差 = {output_value:4.2f}"
            )
            item["reveal"].set_data(t[: i + 1], item["output"][: i + 1])
            set_dot(item["dot"], current_time, output_value)
            artists.extend(
                [
                    item["current_rect"],
                    item["delayed_rect"],
                    item["cursor_top"],
                    item["cursor_bottom"],
                    item["values"],
                    item["reveal"],
                    item["dot"],
                ]
            )
        return artists

    animation = FuncAnimation(fig, update, frames=len(frame_indices), interval=1000 / fps, blit=False, repeat=True)
    save_animation(animation, output_dir / "01_移动窗口形成梯形与三角.gif", fps, dpi)


def make_rccr_animation(data: dict[str, object], output_dir: Path, fps: int, dpi: int) -> None:
    """逐样本显示阶跃经过 CR 和 RC 两级后的传播过程。"""
    t = np.asarray(data["t"])
    step = np.asarray(data["step"])
    cr = normalize(np.asarray(data["cr"]))
    rccr = normalize(np.asarray(data["rccr"]))

    fig, axes = plt.subplots(3, 1, figsize=(9.6, 7.4), sharex=True, constrained_layout=True)
    fig.suptitle("RCCR（RC–CR）信号传播：先抓变化，再把噪声磨圆", fontsize=16, weight="bold")
    series = [
        (step, ORANGE, "① 阶跃输入：高度携带能量信息"),
        (cr, VERMILION, "② CR 高通：恒定部分回零，只保留边沿后的指数脉冲"),
        (rccr, BLUE, "③ RC 低通：尖锐脉冲被摊成圆滑单峰"),
    ]
    animated: list[tuple[plt.Line2D, plt.Line2D, plt.Line2D, np.ndarray]] = []

    for ax, (values, color, title) in zip(axes, series):
        ax.plot(t, values, color=LIGHT_GRAY, lw=2.4)
        reveal, = ax.plot([], [], color=color, lw=3.0)
        dot, = ax.plot([], [], "o", color=BLACK, ms=6, zorder=5)
        cursor = ax.axvline(8, color=GRAY, ls=":", lw=1.2)
        ax.axhline(0, color=BLACK, lw=0.8)
        ax.set(xlim=(8, 82), ylim=(-0.07, 1.12), ylabel="归一化幅度")
        ax.set_title(title, loc="left", weight="bold")
        clean_axes(ax)
        animated.append((reveal, dot, cursor, values))

    axes[-1].set_xlabel("时间 (μs)")
    axes[-1].text(
        0.98,
        0.80,
        r"当 $\tau_{CR}=\tau_{RC}=\tau$：$y(t)\propto(t/\tau)e^{-t/\tau}$",
        transform=axes[-1].transAxes,
        ha="right",
        color=BLUE,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=LIGHT_GRAY, alpha=0.92),
    )
    clock = fig.text(0.86, 0.955, "", ha="right", va="top", color=GRAY, weight="bold")
    frame_indices = np.linspace(np.searchsorted(t, 10.0), np.searchsorted(t, 80.0), 105).astype(int)

    def update(frame: int) -> list[object]:
        i = int(frame_indices[frame])
        current_time = t[i]
        clock.set_text(f"当前样本：t = {current_time:4.1f} μs")
        artists: list[object] = [clock]
        for reveal, dot, cursor, values in animated:
            reveal.set_data(t[: i + 1], values[: i + 1])
            set_dot(dot, current_time, values[i])
            cursor.set_xdata([current_time, current_time])
            artists.extend([reveal, dot, cursor])
        return artists

    animation = FuncAnimation(fig, update, frames=len(frame_indices), interval=1000 / fps, blit=False, repeat=True)
    save_animation(animation, output_dir / "02_RCCR两级动态响应.gif", fps, dpi)


def verify_gif(path: Path, expected_min_frames: int) -> None:
    from PIL import Image

    with Image.open(path) as image:
        frames = getattr(image, "n_frames", 1)
        assert frames >= expected_min_frames, f"{path.name} 只有 {frames} 帧"
        assert image.width >= 700 and image.height >= 450, f"{path.name} 尺寸过小"
        print(f"{path.name}: {frames} 帧, {image.width}×{image.height}, {path.stat().st_size / 1024 / 1024:.1f} MiB")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("animations"), help="GIF 输出目录")
    parser.add_argument("--fps", type=int, default=12, help="每秒帧数")
    parser.add_argument("--dpi", type=int, default=92, help="GIF 分辨率")
    args = parser.parse_args()

    configure_style()
    data = make_demo_signals()
    make_input_signal_animation(args.output_dir, args.fps, args.dpi)
    make_window_animation(data, args.output_dir, args.fps, args.dpi)
    make_rccr_animation(data, args.output_dir, args.fps, args.dpi)
    make_three_by_three_shaping_animation(args.output_dir, args.fps, args.dpi)
    verify_gif(args.output_dir / "00_三种输入信号演示.gif", 110)
    verify_gif(args.output_dir / "01_移动窗口形成梯形与三角.gif", 90)
    verify_gif(args.output_dir / "02_RCCR两级动态响应.gif", 95)
    verify_gif(args.output_dir / "03_三种输入经过三种成形算法.gif", 110)
    print(f"动画已生成：{args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
