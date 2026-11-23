#!/usr/bin/env python3
"""用 Origin 风格绘制 Wang ECC3 的 4n、5n 激发函数与实验数据。"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MultipleLocator, NullFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from plot_excitation_functions import clip_log_polyline

# 用户配置区
CALCULATION_CSV = HERE / "ExFunCalcECC.csv"
ECM_OFFSET_MEV = 134.0  # E* = Ecm - 134 MeV
# 激发能 E* (MeV)、截面 (pb)、能量误差 (MeV)、截面下/上误差 (pb)。
RED_MEASUREMENT = (50, 0.14, 1.07, 0.1, 0.19)
WANG_RED_ENERGY_RANGE = (48.4, 52.3)  # WangECC3 图中红点横向误差棒端点
BLUE_MEASUREMENT = (49.1, 0.18, 1.6, 0.12, 0.44)
# 激发能 E* (MeV)、截面上限 (pb)、能量误差 (MeV)。
UPPER_LIMITS = (
    (45, 1.1, 0.85),
    (42, 0.32, 0.85),
)
BASS_E_STAR_MEV = 44.7
# 用户配置区结束


def draw_excitation_functions(
    energy_star, sigma_4n_pb, sigma_5n_pb, *, label_positions, stroke_scale=1.0,
    show_legend=False, red_measurement=RED_MEASUREMENT,
):
    """绘制以 MeV、pb 为单位的计算曲线和共用实验数据。"""
    font_family = "DejaVu Serif"
    font_paths = [
        Path("/mnt/c/Windows/Fonts") / name
        for name in ("times.ttf", "timesi.ttf", "timesbd.ttf", "timesbi.ttf")
    ]
    if all(path.is_file() for path in font_paths):
        for path in font_paths:
            font_manager.fontManager.addfont(path)
        font_family = "Times New Roman"
    plt.rcParams.update({
        "font.family": "serif", "font.serif": [font_family],
        "mathtext.fontset": "custom", "mathtext.rm": f"{font_family}:bold",
        "mathtext.it": f"{font_family}:bold:italic", "mathtext.bf": f"{font_family}:bold",
        "font.size": 13, "font.weight": "bold",
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.linewidth": 0.8 * stroke_scale,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42, "svg.fonttype": "none",
    })

    red, blue, black = "#D32F2F", "#1565C0", "#000000"
    fig, ax = plt.subplots(figsize=(190 / 25.4, 120 / 25.4))
    fig.subplots_adjust(left=0.12, right=0.975, bottom=0.17, top=0.965)
    ax.set(xlim=(31, 66), ylim=(1e-3, 3), yscale="log")

    for sigma_pb, color, linestyle in ((sigma_4n_pb, black, "--"), (sigma_5n_pb, red, "-")):
        # 裁切函数跳过非正值两侧的线段，缺失值也按断点处理。
        valid_sigma = np.where(np.isfinite(sigma_pb), sigma_pb, 0.0)
        x, y = clip_log_polyline(energy_star, valid_sigma, ax.get_xlim(), ax.get_ylim())
        ax.plot(x, y, color=color, linestyle=linestyle, lw=1.2 * stroke_scale)

    errorbar_style = dict(capsize=3, elinewidth=stroke_scale, capthick=stroke_scale)
    measurement_handles = []
    for point, color, zorder in ((red_measurement, red, 5), (BLUE_MEASUREMENT, blue, 6)):
        center, sigma, energy_error, err_minus, err_plus = point
        handle = ax.errorbar(
            center, sigma, xerr=energy_error, yerr=[[err_minus], [err_plus]],
            fmt="o", color=color, ms=4.5 * stroke_scale, zorder=zorder, **errorbar_style,
        )
        measurement_handles.append(handle.lines[0])
    for center, limit, energy_error in UPPER_LIMITS:
        ax.errorbar(
            center, limit, xerr=energy_error, fmt="none", ecolor=black,
            zorder=4, **errorbar_style,
        )
        # 箭长仅用于表示上限。
        ax.annotate(
            "", xy=(center, limit / 2.1), xytext=(center, limit),
            arrowprops=dict(arrowstyle="-|>", color=black, lw=stroke_scale, mutation_scale=9,
                            shrinkA=0, shrinkB=0),
            zorder=4,
        )

    ax.annotate(
        "", xy=(BASS_E_STAR_MEV, 1.22e-3), xytext=(BASS_E_STAR_MEV, 2.8e-3),
        arrowprops=dict(arrowstyle="-|>", color=black, lw=1.4 * stroke_scale, mutation_scale=12),
        zorder=4,
    )
    ax.text(BASS_E_STAR_MEV + 0.5, 2.65e-3, r"$V_{\mathrm{Bass}}$", fontsize=16, va="center")

    pos_4n, pos_5n = label_positions
    ax.text(*pos_4n, "4n", fontsize=17, color=black)
    ax.text(*pos_5n, "5n", fontsize=17, color=red)
    ax.set_xlabel(r"Excitation energy $E^{*}$ (MeV)")
    ax.set_ylabel(r"Cross section $\sigma$ (pb)")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", length=5, width=0.8 * stroke_scale, pad=6)
    ax.tick_params(which="minor", length=2.5, width=0.6 * stroke_scale)
    if show_legend:
        red_handle, blue_handle = measurement_handles
        limit_handle = Line2D([], [], color=black, marker=r"$\downarrow$",
                              linestyle="none", markersize=13)
        ax.legend(
            [red_handle, limit_handle, blue_handle],
            ["This work (5n)", "This work (4n)", "JINR (5n)"],
            loc="upper right", fontsize=12, frameon=False,
        )
    return fig


def main():
    energy_cm, sigma_4n, sigma_5n = np.loadtxt(
        CALCULATION_CSV, delimiter="\t", skiprows=1, usecols=(0, 5, 6), unpack=True
    )
    mb_to_pb = 1e9
    center, sigma, _, err_minus, err_plus = RED_MEASUREMENT
    low, high = WANG_RED_ENERGY_RANGE
    red_xerr = [[center - low], [high - center]]
    fig = draw_excitation_functions(
        energy_cm - ECM_OFFSET_MEV, sigma_4n * mb_to_pb, sigma_5n * mb_to_pb,
        label_positions=((36.5, 0.018), (51, 0.055)),
        stroke_scale=1.6, show_legend=True,
        red_measurement=(center, sigma, red_xerr, err_minus, err_plus),
    )
    for extension in ("pdf", "svg", "png"):
        fig.savefig(HERE / f"excitation_functions_wang_ecc3.{extension}", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
