#!/usr/bin/env python3
"""Plot Qalpha and alpha partial half-life trends from AnaHs.csv."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache")
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator, LogLocator, MultipleLocator, NullFormatter


ELEMENT_Z = {
    "Rf": 104,
    "Sg": 106,
    "Hs": 108,
    "Ds": 110,
}

ELEMENT_ORDER = ("Rf", "Sg", "Hs", "Ds")
ELEMENT_COLORS = {
    "Rf": "#3B6EA8",
    "Sg": "#C95C1A",
    "Hs": "#2F8A5B",
    "Ds": "#8A5BA8",
}
ELEMENT_MARKERS = {
    "Rf": "o",
    "Sg": "s",
    "Hs": "^",
    "Ds": "D",
}


def parse_isotope(value: str) -> tuple[int, str, bool]:
    match = re.match(r"^\s*(\d+)([A-Z][a-z]?)(.*)$", str(value).strip())
    if not match:
        raise ValueError(f"Cannot parse isotope name: {value!r}")

    mass_number = int(match.group(1))
    symbol = match.group(2)
    suffix = match.group(3)
    if symbol not in ELEMENT_Z:
        raise ValueError(f"Missing Z mapping for element: {symbol}")

    return mass_number, symbol, "^m" in suffix


def load_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df.columns = [column.strip() for column in df.columns]
    df["Isotope"] = df["Isotope"].astype(str).str.strip()

    parsed = df["Isotope"].map(parse_isotope)
    df["A"] = parsed.map(lambda item: item[0])
    df["Element"] = parsed.map(lambda item: item[1])
    df["is_isomer"] = parsed.map(lambda item: item[2])
    df["Z"] = df["Element"].map(ELEMENT_Z)
    df["N"] = df["A"] - df["Z"]

    numeric_columns = [
        "Talpha_1/2_s",
        "Talpha_1/2_err_plus_s",
        "Talpha_1/2_err_minus_s",
        "Qalpha_MeV",
        "resolution/keV",
    ]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df["Qalpha_err_MeV"] = df["resolution/keV"] / 1000.0
    return df.sort_values(["Z", "N", "is_isomer", "A"]).reset_index(drop=True)


def apply_origin_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.25,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(ax: plt.Axes, log_y: bool = False) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.25)
        spine.set_color("black")

    ax.tick_params(
        which="major",
        direction="in",
        top=True,
        right=True,
        length=6,
        width=1.15,
        pad=5,
    )
    ax.tick_params(
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=3.5,
        width=0.95,
    )
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    if log_y:
        ax.yaxis.set_major_locator(LogLocator(base=10.0))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.axvline(162, color="0.25", lw=1.1, ls=(0, (5, 3)), zorder=0)
    ax.text(
        162.25,
        0.965,
        r"$N=162$",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=11,
        color="0.2",
    )
    ax.axvline(152, color="0.25", lw=1.1, ls=(0, (5, 3)), zorder=0)
    ax.text(
        152.25,
        0.965,
        r"$N=152$",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=11,
        color="0.2",
    )

def draw_element_series(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_column: str,
    yerr_columns: tuple[str, str],
    label_elements: bool,
) -> None:
    for element in ELEMENT_ORDER:
        group = df[df["Element"] == element].sort_values(["N", "is_isomer"])
        if group.empty:
            continue

        base = group[~group["is_isomer"]]
        isomer = group[group["is_isomer"]]
        color = ELEMENT_COLORS[element]
        marker = ELEMENT_MARKERS[element]
        label = f"{element} (Z={ELEMENT_Z[element]})" if label_elements else "_nolegend_"

        if not base.empty:
            ax.plot(
                base["N"],
                base[y_column],
                color=color,
                lw=1.35,
                alpha=0.85,
                zorder=1,
            )
            ax.errorbar(
                base["N"],
                base[y_column],
                yerr=np.vstack(
                    [
                        base[yerr_columns[0]].to_numpy(),
                        base[yerr_columns[1]].to_numpy(),
                    ]
                ),
                fmt=marker,
                ms=5.8,
                mfc=color,
                mec=color,
                mew=1.0,
                color=color,
                ecolor=color,
                elinewidth=1.0,
                capsize=3.3,
                capthick=1.0,
                linestyle="none",
                label=label,
                zorder=3,
            )

        if not isomer.empty:
            ax.errorbar(
                isomer["N"],
                isomer[y_column],
                yerr=np.vstack(
                    [
                        isomer[yerr_columns[0]].to_numpy(),
                        isomer[yerr_columns[1]].to_numpy(),
                    ]
                ),
                fmt=marker,
                ms=6.2,
                mfc="white",
                mec=color,
                mew=1.35,
                color=color,
                ecolor=color,
                elinewidth=1.0,
                capsize=3.3,
                capthick=1.0,
                linestyle="none",
                label="_nolegend_",
                zorder=4,
            )
            for _, row in isomer.iterrows():
                ax.annotate(
                    "m",
                    xy=(row["N"], row[y_column]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=9,
                    color=color,
                )


def make_combined_figure(df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 8.0), sharex=True)

    draw_element_series(
        axes[0],
        df,
        "Qalpha_MeV",
        ("Qalpha_err_MeV", "Qalpha_err_MeV"),
        label_elements=True,
    )
    style_axis(axes[0], log_y=False)
    axes[0].set_ylabel(r"$Q_{\alpha}$ (MeV)")
    axes[0].set_ylim(8.1, 11.75)
    axes[0].text(
        0.02,
        0.94,
        r"   (a)",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=16,
    )
    axes[0].legend(loc="upper right", frameon=False, handlelength=2.2)

    draw_element_series(
        axes[1],
        df,
        "Talpha_1/2_s",
        ("Talpha_1/2_err_minus_s", "Talpha_1/2_err_plus_s"),
        label_elements=False,
    )
    axes[1].set_yscale("log")
    style_axis(axes[1], log_y=True)
    axes[1].set_xlabel(r"Neutron number $N$")
    axes[1].set_ylabel(r"$T_{1/2}^{\alpha}$ (s)")
    axes[1].set_ylim(7e-5, 2e4)
    axes[1].text(
        0.02,
        0.94,
        r"   (b)",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=16,
    )

    axes[1].set_xlim(149, 172)
    fig.subplots_adjust(hspace=0.06)
    return fig


def make_single_figure(
    df: pd.DataFrame,
    y_column: str,
    yerr_columns: tuple[str, str],
    ylabel: str,
    panel_label: str,
    output_stem: Path,
    log_y: bool = False,
    ylim: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    draw_element_series(ax, df, y_column, yerr_columns, label_elements=True)
    if log_y:
        ax.set_yscale("log")
    style_axis(ax, log_y=log_y)
    ax.set_xlabel(r"Neutron number $N$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(149, 172)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.text(
        0.02,
        0.94,
        panel_label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )
    ax.legend(loc="best", frameon=False, handlelength=2.2)
    fig.savefig(output_stem.with_suffix(".png"))
    fig.savefig(output_stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Draw Origin-style alpha-decay trend plots from AnaHs.csv."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).with_name("AnaHs.csv"),
        help="Input CSV table.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory.",
    )
    args = parser.parse_args()

    apply_origin_style()
    df = load_data(args.csv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    combined = make_combined_figure(df)
    combined_stem = args.outdir / "AnaHs_alpha_trends_origin"
    combined.savefig(combined_stem.with_suffix(".png"))
    combined.savefig(combined_stem.with_suffix(".pdf"))
    plt.close(combined)

    make_single_figure(
        df,
        "Qalpha_MeV",
        ("Qalpha_err_MeV", "Qalpha_err_MeV"),
        r"$Q_{\alpha}$ (MeV)",
        r"",
        args.outdir / "AnaHs_Qalpha_trend_origin",
        log_y=False,
        ylim=(8.1, 11.75),
    )
    make_single_figure(
        df,
        "Talpha_1/2_s",
        ("Talpha_1/2_err_minus_s", "Talpha_1/2_err_plus_s"),
        r"$T_{1/2}^{\alpha}$ (s)",
        r"",
        args.outdir / "AnaHs_Talpha_trend_origin",
        log_y=True,
        ylim=(7e-5, 2e4),
    )

    print(f"Wrote {combined_stem.with_suffix('.png')}")
    print(f"Wrote {combined_stem.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
