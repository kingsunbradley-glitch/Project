#!/usr/bin/env python3
"""Draw a publication-ready waveform from a two-column sample/ADC text file."""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator


def read_waveform(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two numeric columns in {path}")

    sample = data[:, 0]
    adc = data[:, 1]
    return sample, adc


def set_origin_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "ytick.minor.size": 2.0,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.6,
            "ytick.minor.width": 0.6,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def nice_step(span: float, target_ticks: int = 5) -> float:
    if span <= 0:
        return 1.0

    raw_step = span / target_ticks
    exponent = math.floor(math.log10(raw_step))
    base = 10**exponent
    fraction = raw_step / base

    for nice_fraction in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= nice_fraction:
            return nice_fraction * base
    return 10.0 * base


def plot_waveform(
    input_path: Path,
    output_prefix: Path,
    sample_interval_ns: float,
    dpi: int,
) -> None:
    sample, adc = read_waveform(input_path)
    time_us = sample * sample_interval_ns / 1000.0

    set_origin_style()

    fig, ax = plt.subplots(figsize=(3.45, 2.35), constrained_layout=True)
    ax.plot(
        time_us,
        adc,
        color="#1f4e79",
        linewidth=0.9,
        solid_joinstyle="round",
        solid_capstyle="round",
    )

    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel("ADC channel")

    x_step = nice_step(time_us.max() - time_us.min())
    x_min = math.floor(time_us.min() / x_step) * x_step
    x_max = math.ceil(time_us.max() / x_step) * x_step
    ax.set_xlim(x_min, x_max)
    y_margin = 0.05 * (adc.max() - adc.min())
    ax.set_ylim(adc.min() - y_margin, adc.max() + y_margin)

    ax.xaxis.set_major_locator(MultipleLocator(x_step))
    ax.xaxis.set_minor_locator(MultipleLocator(x_step / 2.0))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=True, right=True)

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "eps", "png"):
        fig.savefig(output_prefix.with_suffix(f".{ext}"), dpi=dpi)

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot waveform data with an Origin-like publication style."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path(__file__).with_name("wave_info.txt"),
        help="Input text file with sample and ADC columns.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        default=Path(__file__).with_name("waveform_publication"),
        help="Output path prefix; .pdf, .svg, .eps, and .png are written.",
    )
    parser.add_argument(
        "--sample-interval-ns",
        type=float,
        default=10.0,
        help="Sampling interval in ns per point.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG output DPI.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot_waveform(
        input_path=args.input,
        output_prefix=args.output_prefix,
        sample_interval_ns=args.sample_interval_ns,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
