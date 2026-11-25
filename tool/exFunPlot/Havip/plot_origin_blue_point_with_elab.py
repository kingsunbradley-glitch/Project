#!/usr/bin/env python3
"""在 Origin 风格激发函数图顶部增加计算表中的 Elab 横轴。"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import LogLocator, MultipleLocator, NullFormatter

from plot_excitation_functions import (
    CALCULATION_CSV,
    ORIGIN_HIDDEN_UPPER_LIMIT,
    clip_log_polyline,
    configured_experimental_rows,
    read_calculation,
)


HERE = Path(__file__).resolve().parent

# 用户配置区。其余实验点和误差见 plot_excitation_functions.py。
BLUE_ENERGY_RANGE_MEV = (47.5, 50.7)
BLUE_SIGMA_PB = 0.18
BLUE_SIGMA_ERRORS_PB = (0.12, 0.44)  # 下误差、上误差
BASS_ENERGY_MEV = 43.8
E_STAR_MAJOR_STEP_MEV = 5
E_LAB_RANGE_MEV = (190, 240)
E_LAB_MAJOR_STEP_MEV = 5


def interp_with_linear_extrapolation(values, source, target):
    """线性插值，并在数据表两端沿相邻两点直线外推。"""
    values = np.asarray(values, dtype=float)
    converted = np.interp(values, source, target)
    low = values < source[0]
    high = values > source[-1]
    converted = np.where(
        low,
        target[0] + (values - source[0]) * (target[1] - target[0]) / (source[1] - source[0]),
        converted,
    )
    converted = np.where(
        high,
        target[-1]
        + (values - source[-1]) * (target[-1] - target[-2]) / (source[-1] - source[-2]),
        converted,
    )
    return converted


def energy_converters(data):
    """返回由 ExFunCalc.csv 构造的 E* 与 Elab 双向换算函数。"""
    e_star = data["E*"]
    e_lab = data["Elab"]

    def e_star_to_e_lab(values):
        return interp_with_linear_extrapolation(values, e_star, e_lab)

    def e_lab_to_e_star(values):
        return interp_with_linear_extrapolation(values, e_lab, e_star)

    return e_star_to_e_lab, e_lab_to_e_star


def add_elab_axis(ax, e_star_to_e_lab, e_lab_to_e_star):
    """按数据表换算关系，在顶部增加 190--240 MeV 的 Elab 横轴。"""
    top_ax = ax.secondary_xaxis(
        "top", functions=(e_star_to_e_lab, e_lab_to_e_star)
    )
    top_ax.xaxis.set_major_locator(MultipleLocator(E_LAB_MAJOR_STEP_MEV))
    top_ax.xaxis.set_minor_locator(MultipleLocator(1))
    top_ax.set_xlabel(r"Laboratory energy $E_{\mathrm{lab}}$ (MeV)", labelpad=8)
    top_ax.tick_params(which="major", direction="in", length=5, width=0.8, pad=5)
    top_ax.tick_params(which="minor", direction="in", length=2.5, width=0.6)
    top_ax.spines["top"].set_linewidth(0.8)
    for side in ("bottom", "left", "right"):
        top_ax.spines[side].set_visible(False)
    return top_ax


# 加载真正的 Times New Roman，保持与原图一致。
for name in ("times.ttf", "timesi.ttf", "timesbd.ttf", "timesbi.ttf"):
    font_manager.fontManager.addfont("/mnt/c/Windows/Fonts/" + name)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman:bold",
        "mathtext.it": "Times New Roman:bold:italic",
        "mathtext.bf": "Times New Roman:bold",
        "font.size": 13,
        "font.weight": "bold",
        "axes.labelsize": 16,
        "axes.labelweight": "bold",
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

data = read_calculation(CALCULATION_CSV)
rows = configured_experimental_rows()
e_star_to_e_lab, e_lab_to_e_star = energy_converters(data)
e_star_range = tuple(e_lab_to_e_star(E_LAB_RANGE_MEV))

fig, ax = plt.subplots(figsize=(190 / 25.4, 120 / 25.4))
# 顶部为 Elab 的刻度标签和轴标题留出空间。
fig.subplots_adjust(left=0.12, right=0.975, bottom=0.17, top=0.80)
ax.set(xlim=e_star_range, ylim=(1.1e-3, 3), yscale="log")

for isotope, color, style, label in (
    ("274Ds", "#000000", "--", r"$4n$: $^{274}$Ds (calc.)"),
    ("273Ds", "#D32F2F", "-", r"$5n$: $^{273}$Ds (calc.)"),
):
    x, y = clip_log_polyline(data["E*"], data[isotope], ax.get_xlim(), ax.get_ylim())
    ax.plot(x, y, color=color, linestyle=style, lw=1.2, label=label)

# 两个道的上限及能量误差相同时，共用黑色符号。
seen = set()
hidden_channel, hidden_energy = ORIGIN_HIDDEN_UPPER_LIMIT
for row in rows:
    x, y = float(row["E_star_MeV"]), float(row["sigma_pb"])
    xe = np.array(
        [[float(row["E_star_err_minus_MeV"])], [float(row["E_star_err_plus_MeV"])]]
    )
    if row["kind"] == "upper_limit":
        if row["channel"] == hidden_channel and np.isclose(x, hidden_energy):
            continue
        key = (x, y, *xe.ravel())
        if key in seen:
            continue
        seen.add(key)
        ax.errorbar(
            x,
            y,
            xerr=xe,
            fmt="none",
            ecolor="#000000",
            elinewidth=1.0,
            capsize=3,
            zorder=4,
        )
        # 箭长仅用于表示上限，不表示误差。
        ax.annotate(
            "",
            xy=(x, y / 2.1),
            xytext=(x, y),
            arrowprops=dict(
                arrowstyle="-|>",
                color="#000000",
                lw=1.0,
                mutation_scale=9,
                shrinkA=0,
                shrinkB=0,
            ),
            zorder=4,
        )
    else:
        ax.errorbar(
            x,
            y,
            xerr=xe,
            yerr=[[float(row["err_minus_pb"])], [float(row["err_plus_pb"])]],
            fmt="o",
            color="#D32F2F",
            ms=4.5,
            capsize=3,
            elinewidth=1.0,
            zorder=5,
        )

# 蓝点的横向误差棒覆盖给定能区。
energy_low, energy_high = BLUE_ENERGY_RANGE_MEV
energy_center = (energy_low + energy_high) / 2
ax.errorbar(
    energy_center,
    BLUE_SIGMA_PB,
    xerr=[[energy_center - energy_low], [energy_high - energy_center]],
    yerr=np.array(BLUE_SIGMA_ERRORS_PB).reshape(2, 1),
    fmt="o",
    color="#1565C0",
    ms=4.5,
    capsize=3,
    elinewidth=1.0,
    zorder=6,
)

# 势垒位置；箭头是能量标记，不代表截面测量。
ax.annotate(
    "",
    xy=(BASS_ENERGY_MEV, 1.22e-3),
    xytext=(BASS_ENERGY_MEV, 2.8e-3),
    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.4, mutation_scale=12),
)
ax.text(
    BASS_ENERGY_MEV + 0.5,
    2.65e-3,
    r"$V_{\mathrm{BASS}}$",
    fontsize=16,
    fontweight="bold",
    va="center",
)
ax.text(36.3, 0.18, "4n", fontsize=17, fontweight="bold", color="black")
ax.text(54.2, 0.16, "5n", fontsize=17, fontweight="bold", color="#D32F2F")

ax.set_xlabel(r"Excitation energy $E^{*}$ (MeV)")
ax.set_ylabel(r"Cross section $\sigma$ (pb)")
ax.xaxis.set_major_locator(MultipleLocator(E_STAR_MAJOR_STEP_MEV))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.tick_params(which="major", length=5, width=0.8, pad=6)
ax.tick_params(which="minor", length=2.5, width=0.6)
# 顶边只显示 Elab 副轴刻度，避免与底部 E* 轴的顶部镜像刻度重叠。
ax.tick_params(axis="x", which="both", top=False)

add_elab_axis(ax, e_star_to_e_lab, e_lab_to_e_star)

out = HERE / "figures"
out.mkdir(exist_ok=True)
for ext in ("pdf", "svg", "png"):
    fig.savefig(out / f"excitation_functions_origin_blue_point_with_elab.{ext}", dpi=600)
plt.close(fig)
print(out / "excitation_functions_origin_blue_point_with_elab.pdf")
