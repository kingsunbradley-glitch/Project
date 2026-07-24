#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plot even-N log10(T1/2 SF / s) systematics for Cf, Fm, and No."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "matplotlib-codex-cache"),
)

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_PATH = SCRIPT_DIR / "Fmsys.dat"
OUT_STEMS = (
    SCRIPT_DIR / "Cf_Fm_No_LogTsf_N",
    SCRIPT_DIR / "Cf_Fm_No_evenN_LogTsf_N",
)
FIGURE_FORMATS = ("pdf", "png", "tiff")

ELEMENTS = {
    98: {
        "symbol": "Cf",
        "label": "Cf (Z=98)",
        "color": "#CE5A0C",
        "marker": "o",
    },
    100: {
        "symbol": "Fm",
        "label": "Fm (Z=100)",
        "color": "#3210D8",
        "marker": "s",
    },
    102: {
        "symbol": "No",
        "label": "No (Z=102)",
        "color": "#009E73",
        "marker": "^",
    },
}


def apply_origin_publication_style() -> None:
    """Origin-like boxed axes and publication-friendly typography."""
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 9.0,
            "axes.labelsize": 10.5,
            "axes.linewidth": 1.15,
            "axes.edgecolor": "black",
            "xtick.labelsize": 8.8,
            "ytick.labelsize": 8.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
            "xtick.minor.size": 2.8,
            "ytick.minor.size": 2.8,
            "xtick.major.width": 1.05,
            "ytick.major.width": 1.05,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "legend.fontsize": 8.5,
            "legend.frameon": False,
            "lines.linewidth": 1.35,
            "lines.markersize": 5.4,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_even_neutron_ground_state_data(path: Path) -> pd.DataFrame:
    """Read Fmsys.dat, keep even-N ground states, and compute log10(Tsf/s)."""
    df = pd.read_csv(path, sep=r"\s+", dtype={"A": "string"})

    required = {"A", "Z", "N", "T1/2_SF", "T1/2_SF_limit"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")

    df = df[df["A"].str.fullmatch(r"\d+", na=False)].copy()
    for col in ["Z", "N", "T1/2_SF"]:
        df[col] = pd.to_numeric(df[col], errors="raise")

    df = df[df["N"] % 2 == 0].copy()
    df["T1/2_SF_limit"] = df["T1/2_SF_limit"].fillna("value").astype(str)
    df["logTsf"] = np.log10(df["T1/2_SF"])

    return df.sort_values(["Z", "N", "A"], kind="stable")


def plot_connected_segments(ax: plt.Axes, data: pd.DataFrame, color: str) -> None:
    """Draw solid connecting lines, leaving the tabulated 242Fm point unconnected."""
    if data.empty:
        return

    if int(data["Z"].iloc[0]) == 100:
        segments = [data[data["N"] < 142], data[data["N"] > 142]]
    else:
        segments = [data]

    for segment in segments:
        if len(segment) < 2:
            continue
        ax.plot(
            segment["N"],
            segment["logTsf"],
            "-",
            color=color,
            linewidth=1.35,
            zorder=2,
        )


def add_242fm_curve_connections(ax: plt.Axes, fm_data: pd.DataFrame) -> None:
    """Connect the two N=142 242Fm values to the main Fm systematics curve."""
    n_value = 142
    fm_color = ELEMENTS[100]["color"]
    anchor_points = fm_data[fm_data["N"] > n_value].sort_values("N", kind="stable")
    if anchor_points.empty:
        return

    anchor = anchor_points.iloc[0]
    anchor_x = anchor["N"]
    anchor_y = anchor["logTsf"]

    star_log_tsf = np.log10(4.0e-6)
    ax.plot(
        [n_value, anchor_x],
        [star_log_tsf, anchor_y],
        "-",
        color="red",
        linewidth=1.35,
        zorder=2.4,
    )

    tabulated_242fm = fm_data[fm_data["N"] == n_value]
    if tabulated_242fm.empty:
        return

    tabulated = tabulated_242fm.iloc[0]
    ax.plot(
        [tabulated["N"], anchor_x],
        [tabulated["logTsf"], anchor_y],
        ":",
        color=fm_color,
        linewidth=1.35,
        zorder=2.3,
    )


def add_limit_arrows(ax: plt.Axes, data: pd.DataFrame, color: str) -> None:
    """Draw down arrows for upper limits and up arrows for lower limits."""
    arrow_span = 0.90
    for _, row in data.iterrows():
        limit = row["T1/2_SF_limit"].strip().lower()
        if limit == "upper":
            xy = (row["N"], row["logTsf"] - arrow_span)
            xytext = (row["N"], row["logTsf"])
        elif limit == "lower":
            xy = (row["N"], row["logTsf"] + arrow_span)
            xytext = (row["N"], row["logTsf"])
        else:
            continue

        ax.annotate(
            "",
            xy=xy,
            xytext=xytext,
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "lw": 0.95,
                "mutation_scale": 8.5,
                "shrinkA": 4.0,
                "shrinkB": 0.0,
            },
            zorder=4,
        )


def plot_element(ax: plt.Axes, data: pd.DataFrame, spec: dict[str, str]) -> Line2D:
    """Plot one isotope chain and return a legend handle."""
    plot_connected_segments(ax, data, spec["color"])
    ax.plot(
        data["N"],
        data["logTsf"],
        linestyle="none",
        marker=spec["marker"],
        markersize=5.8,
        markerfacecolor=spec["color"],
        markeredgecolor=spec["color"],
        markeredgewidth=0.9,
        color=spec["color"],
        zorder=3,
    )
    add_limit_arrows(ax, data, spec["color"])

    return Line2D(
        [0],
        [0],
        color=spec["color"],
        marker=spec["marker"],
        linestyle="-",
        linewidth=1.35,
        markersize=5.8,
        markerfacecolor=spec["color"],
        markeredgecolor=spec["color"],
        label=spec["label"],
    )


def add_242fm_star(ax: plt.Axes) -> Line2D:
    """Add the special 242Fm point at 4 us."""
    n_value = 142
    log_tsf = np.log10(4.0e-6)

    ax.plot(
        [n_value],
        [log_tsf],
        linestyle="none",
        marker="*",
        markersize=12.5,
        markerfacecolor="red",
        markeredgecolor="red",
        markeredgewidth=0.9,
        color="red",
        zorder=6,
    )

    return Line2D(
        [0],
        [0],
        color="red",
        marker="*",
        linestyle="none",
        markersize=12.5,
        markerfacecolor="red",
        markeredgecolor="red",
        label=r"$^{242}\mathrm{Fm}$",
    )


def draw_figure(df: pd.DataFrame) -> tuple[plt.Figure, dict[str, int]]:
    apply_origin_publication_style()

    fig, ax = plt.subplots(figsize=(5.1, 3.65), constrained_layout=True)

    handles: list[Line2D] = []
    counts: dict[str, int] = {}
    fm_data: pd.DataFrame | None = None
    for z, spec in ELEMENTS.items():
        data = df[df["Z"] == z].sort_values("N", kind="stable")
        counts[spec["symbol"]] = len(data)
        if z == 100:
            fm_data = data
        handles.append(plot_element(ax, data, spec))

    if fm_data is not None:
        add_242fm_curve_connections(ax, fm_data)
    handles.append(add_242fm_star(ax))
    counts["242Fm_star"] = 1

    ax.axvline(152, linestyle=(0, (4, 2)), linewidth=0.95, color="0.35", zorder=1)
    ax.text(
        152,
        -7,
        r"$N=152$",
        ha="center",
        va="center",
        fontsize=11,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.7},
    )
    element_labels = [
        (151, -3, "No", ELEMENTS[102]["color"]),
        (148, 1.5, "Fm", ELEMENTS[100]["color"]),
        (146, 6.5, "Cf", ELEMENTS[98]["color"]),
    ]
    for x, y, text, color in element_labels:
        ax.text(
            x,
            y,
            text,
            color=color,
            fontsize=13.0,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5},
            zorder=7,
        )

    ax.set_xlabel(r"Neutron number")
    ax.set_ylabel(r"$\log_{10}(\mathrm{T}_{1/2}^{\mathrm{SF}}/\mathrm{s})$")
    ax.set_xlim(138, 162)
    ax.set_ylim(-10, 14)

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.tick_params(which="both", direction="in", top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.15)
        spine.set_color("black")

    # Explicit legend is disabled; element labels are drawn inside the axes.
    # ax.legend(
    #     handles=handles,
    #     loc="upper left",
    #     handlelength=2.2,
    #     borderaxespad=0.6,
    #     labelspacing=0.45,
    # )

    return fig, counts


def save_figure(fig: plt.Figure, out_stems: tuple[Path, ...]) -> list[Path]:
    paths = []
    for out_stem in out_stems:
        for suffix in FIGURE_FORMATS:
            path = out_stem.with_suffix(f".{suffix}")
            fig.savefig(path)
            paths.append(path)
    return paths


def main() -> None:
    df = load_even_neutron_ground_state_data(DATA_PATH)
    fig, counts = draw_figure(df)
    paths = save_figure(fig, OUT_STEMS)
    plt.close(fig)

    print("Plotted points:")
    for key in ("Cf", "Fm", "No", "242Fm_star"):
        print(f"  {key}: {counts[key]}")
    for path in paths:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
