#!/usr/bin/env python3
"""以统一的 Origin 风格绘制各模型的全部蒸发道激发函数。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MultipleLocator, NullFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
EXFUN_ROOT = HERE.parent

# =============================================================================
# 用户配置区
# =============================================================================

# 坐标范围。按要求，纵轴只画到 10^-2 pb，不显示更低的截面。
X_LIMITS_MEV = (30.0, 62.0)
Y_LIMITS_PB = (1.0e-2, 3.0)

# WangECC3 的 Q 值；激发能按 E* = Ecm + Q 计算。
WANG_ECC3_Q_MEV = -132.481

# 红、蓝实验点。
# energy / sigma: 点坐标；xerr_left/right: 左右能量误差；
# yerr_down/up: 向下、向上截面误差。单位依次为 MeV 和 pb。
EXPERIMENTAL_POINTS = (
    {
        "label": "This work (5n)",
        "channel": "5n",
        "color": "#D32F2F",
        "energy": 50.0,
        "sigma": 0.14,
        "xerr_left": 1.0,
        "xerr_right": 1.2,
        "yerr_down": 0.08,
        "yerr_up": 0.18,
    },
    {
        "label": "JINR (5n)",
        "channel": "5n",
        "color": "#1565C0",
        "energy": 49.1,
        "sigma": 0.18,
        "xerr_left": 1.6,
        "xerr_right": 1.6,
        "yerr_down": 0.12,
        "yerr_up": 0.44,
    },
)

# 截面上限。energy / limit 是上限坐标，xerr_left/right 是左右能量误差。
# 向下箭头只表示“上限”，箭头长度不是截面误差。
UPPER_LIMITS = (
    {
        "label": "This work (4n)",
        "channel": "4n",
        "energy": 45.0,
        "limit": 1.10,
        "xerr_left": 1.0,
        "xerr_right": 1.0,
    },
    {
        "label": "This work (4n)",
        "channel": "4n",
        "energy": 42.0,
        "limit": 0.34,
        "xerr_left": 1.0,
        "xerr_right": 1.0,
    },
)

# 如需势垒箭头，把 None 改为激发能数值，例如 44.7。
BARRIER_ENERGY_MEV = None
BARRIER_LABEL = r"$V_{\mathrm{Bass}}$"

# 图幅和输出格式。
FIGURE_SIZE_MM = (190.0, 120.0)
OUTPUT_FORMATS = ("pdf", "svg", "png")
PNG_DPI = 600

# 用户配置区结束
# =============================================================================


CHANNEL_STYLES = {
    "2n": {"color": "#7B3294", "linestyle": ":"},
    "3n": {"color": "#0072B2", "linestyle": "-."},
    "4n": {"color": "#000000", "linestyle": "--"},
    "5n": {"color": "#D32F2F", "linestyle": "-"},
    "6n": {"color": "#008837", "linestyle": (0, (5, 1.5, 1.3, 1.5))},
}

MODEL_SPECS = {
    "havip": {
        "title": "Havip",
        "path": EXFUN_ROOT / "Havip" / "ExFunCalc.csv",
        "output": "havip_all_channels_origin",
    },
    "dns_li": {
        "title": "DNS_Li",
        "path": EXFUN_ROOT / "DNS_Li" / "excitation_functions_2n_5n.csv",
        "output": "dns_li_all_channels_origin",
    },
    "wang_ecc3": {
        "title": "WangECC3",
        "path": EXFUN_ROOT / "WangECC3" / "ExFunCalcECC.csv",
        "output": "wang_ecc3_all_channels_origin",
    },
}


def _float_column(rows, name):
    """把 CSV 列转成浮点数组；空单元格作为 NaN 保留为断点。"""
    values = []
    for row in rows:
        text = row[name].strip()
        values.append(float(text) if text else np.nan)
    return np.asarray(values, dtype=float)


def read_havip(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    energy = _float_column(rows, "E*")
    return energy, {
        "3n": _float_column(rows, "275Ds"),
        "4n": _float_column(rows, "274Ds"),
        "5n": _float_column(rows, "273Ds"),
        "6n": _float_column(rows, "272Ds"),
    }


def read_dns_li(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    energy = _float_column(rows, "E_star_MeV")
    return energy, {
        channel: _float_column(rows, f"sigma_{channel}_pb")
        for channel in ("2n", "3n", "4n", "5n")
    }


def read_wang_ecc3(path):
    # 原表为 Ecm 和 mb；E* = Ecm + Q，1 mb = 10^9 pb。
    table = np.loadtxt(path, delimiter="\t", skiprows=1)
    energy = table[:, 0] + WANG_ECC3_Q_MEV
    return energy, {
        channel: table[:, column] * 1.0e9
        for channel, column in zip(("2n", "3n", "4n", "5n"), (3, 4, 5, 6), strict=True)
    }


MODEL_READERS = {
    "havip": read_havip,
    "dns_li": read_dns_li,
    "wang_ecc3": read_wang_ecc3,
}


def validate_calculation(model_key, energy, curves):
    if energy.ndim != 1 or energy.size < 2 or not np.isfinite(energy).all():
        raise ValueError(f"{model_key}: invalid excitation-energy column")
    if not np.all(np.diff(energy) > 0):
        raise ValueError(f"{model_key}: excitation energies must increase strictly")
    for channel, sigma in curves.items():
        if sigma.shape != energy.shape:
            raise ValueError(f"{model_key}/{channel}: energy and cross-section lengths differ")
        finite = sigma[np.isfinite(sigma)]
        if finite.size == 0 or np.any(finite < 0):
            raise ValueError(f"{model_key}/{channel}: cross sections are empty or negative")


def validate_experiment():
    for point in EXPERIMENTAL_POINTS:
        values = [point[key] for key in (
            "energy", "sigma", "xerr_left", "xerr_right", "yerr_down", "yerr_up"
        )]
        if not np.isfinite(values).all() or any(value < 0 for value in values[2:]):
            raise ValueError(f"Invalid experimental point: {point['label']}")
        if point["sigma"] <= point["yerr_down"]:
            raise ValueError(f"Lower error reaches zero on log axis: {point['label']}")
    for limit in UPPER_LIMITS:
        values = [limit[key] for key in ("energy", "limit", "xerr_left", "xerr_right")]
        if not np.isfinite(values).all() or limit["limit"] <= 0 or any(value < 0 for value in values[2:]):
            raise ValueError(f"Invalid upper limit: {limit['label']}")


def configure_origin_style():
    family = "DejaVu Serif"
    windows_font_dir = Path("/mnt/c/Windows/Fonts")
    times_fonts = [windows_font_dir / name for name in (
        "times.ttf", "timesi.ttf", "timesbd.ttf", "timesbi.ttf"
    )]
    if all(path.is_file() for path in times_fonts):
        for path in times_fonts:
            font_manager.fontManager.addfont(path)
        family = "Times New Roman"

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [family],
        "mathtext.fontset": "stix",
        "font.size": 13,
        "font.weight": "bold",
        "axes.labelsize": 16,
        "axes.labelweight": "bold",
        "axes.linewidth": 1.25,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def plot_experiment(ax):
    handles = []
    error_style = {"capsize": 3.5, "elinewidth": 1.25, "capthick": 1.25}
    for point in EXPERIMENTAL_POINTS:
        handle = ax.errorbar(
            point["energy"], point["sigma"],
            xerr=np.array([[point["xerr_left"]], [point["xerr_right"]]]),
            yerr=np.array([[point["yerr_down"]], [point["yerr_up"]]]),
            fmt="o", color=point["color"], markersize=5.5, zorder=6,
            **error_style,
        )
        handles.append((handle.lines[0], point["label"]))

    seen_labels = set()
    for limit in UPPER_LIMITS:
        ax.errorbar(
            limit["energy"], limit["limit"],
            xerr=np.array([[limit["xerr_left"]], [limit["xerr_right"]]]),
            fmt="none", ecolor="black", zorder=5, **error_style,
        )
        arrow_end = max(limit["limit"] / 2.1, Y_LIMITS_PB[0] * 1.04)
        ax.annotate(
            "", xy=(limit["energy"], arrow_end), xytext=(limit["energy"], limit["limit"]),
            arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.25,
                        "mutation_scale": 10, "shrinkA": 0, "shrinkB": 0},
            zorder=5,
        )
        if limit["label"] not in seen_labels:
            upper_handle = Line2D([], [], color="black", marker=r"$\downarrow$",
                                  linestyle="none", markersize=13)
            handles.append((upper_handle, limit["label"]))
            seen_labels.add(limit["label"])
    return handles


def draw_model(model_key, output_dir):
    spec = MODEL_SPECS[model_key]
    energy, curves = MODEL_READERS[model_key](spec["path"])
    validate_calculation(model_key, energy, curves)

    width, height = (value / 25.4 for value in FIGURE_SIZE_MM)
    fig, ax = plt.subplots(figsize=(width, height))
    fig.subplots_adjust(left=0.12, right=0.975, bottom=0.17, top=0.965)
    ax.set(xlim=X_LIMITS_MEV, ylim=Y_LIMITS_PB, yscale="log")

    curve_handles = []
    for channel, sigma in curves.items():
        style = CHANNEL_STYLES[channel]
        # 零值和空值在对数轴上保留为断点，不用任意小正数替换。
        visible_sigma = np.where(np.isfinite(sigma) & (sigma > 0), sigma, np.nan)
        handle, = ax.plot(
            energy, visible_sigma, color=style["color"], linestyle=style["linestyle"],
            linewidth=1.8, label=channel, zorder=2,
        )
        curve_handles.append(handle)

    experiment_handles = plot_experiment(ax)

    if BARRIER_ENERGY_MEV is not None:
        arrow_bottom = Y_LIMITS_PB[0] * 1.05
        arrow_top = Y_LIMITS_PB[0] * 2.3
        ax.annotate(
            "", xy=(BARRIER_ENERGY_MEV, arrow_bottom),
            xytext=(BARRIER_ENERGY_MEV, arrow_top),
            arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.4,
                        "mutation_scale": 12},
        )
        ax.text(BARRIER_ENERGY_MEV + 0.6, arrow_top, BARRIER_LABEL,
                fontsize=16, va="center")

    ax.text(0.018, 0.968, spec["title"], transform=ax.transAxes,
            fontsize=16, ha="left", va="top")
    ax.set_xlabel(r"Excitation energy $E^{*}$ (MeV)")
    ax.set_ylabel(r"Cross section $\sigma$ (pb)")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", length=6.5, width=1.25, pad=6)
    ax.tick_params(which="minor", length=3.2, width=0.9)

    channel_legend = ax.legend(
        handles=curve_handles, labels=[handle.get_label() for handle in curve_handles],
        loc="upper left", bbox_to_anchor=(0.0, 0.89), ncol=2,
        fontsize=12, frameon=False, handlelength=2.6,
        columnspacing=1.2, labelspacing=0.35,
    )
    ax.add_artist(channel_legend)
    ax.legend(
        handles=[item[0] for item in experiment_handles],
        labels=[item[1] for item in experiment_handles],
        loc="upper right", fontsize=12, frameon=False,
        handlelength=1.8, labelspacing=0.5,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / spec["output"]
    for extension in OUTPUT_FORMATS:
        save_options = {"dpi": PNG_DPI} if extension == "png" else {}
        fig.savefig(stem.with_suffix(f".{extension}"), **save_options)
    plt.close(fig)
    return stem, curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model", choices=("all", *MODEL_SPECS), default="all",
        help="只画指定模型；默认画全部模型",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    args = parser.parse_args()

    validate_experiment()
    configure_origin_style()
    selected = MODEL_SPECS if args.model == "all" else (args.model,)
    for model_key in selected:
        stem, curves = draw_model(model_key, args.output_dir)
        below = [channel for channel, values in curves.items()
                 if np.nanmax(values) < Y_LIMITS_PB[0]]
        note = f"; entirely below y range: {', '.join(below)}" if below else ""
        print(f"Saved {stem}.pdf/.svg/.png{note}")


if __name__ == "__main__":
    main()
