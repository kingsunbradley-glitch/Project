#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare measured and predicted even-N spontaneous-fission half-lives."""

from __future__ import annotations

import os
import re
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

from plot_fmsys_logtsf import (
    DATA_PATH,
    ELEMENTS,
    FIGURE_FORMATS,
    add_242fm_star,
    apply_origin_publication_style,
    load_even_neutron_ground_state_data,
    plot_element,
)


SCRIPT_DIR = Path(__file__).resolve().parent
THEORY_PATH = SCRIPT_DIR / "TsfFromCPC.txt"
OUT_STEM = SCRIPT_DIR / "Cf_Fm_No_Exp_vs_Theory_LogTsf_N"

# The source table labels its fourth column "EBT", while its reference text and
# the published model name use "ETB" (effective tunneling barrier).
THEORY_STYLES = {
    "EBT": {
        "label": "ETB",
        "linestyle": (0, (5.0, 2.2)),
        "marker": "x",
        "markersize": 4.4,
    },
    "MAA": {
        "label": "MAA",
        "linestyle": (0, (1.2, 2.0)),
        "marker": "D",
        "markersize": 3.7,
    },
}

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
THEORY_ROW = re.compile(
    rf"^\s*(\d+)\s+(\d+)\s+(\d+)\s+({NUMBER})\s+({NUMBER}|--)\s*$"
)


def load_theory_data(path: Path) -> pd.DataFrame:
    """Read the CPC table; EBT and MAA are already log10(TSF/s)."""
    rows: list[dict[str, int | float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = THEORY_ROW.match(line)
        if match is None:
            continue

        mass, proton, neutron, ebt, maa = match.groups()
        rows.append(
            {
                "A": int(mass),
                "Z": int(proton),
                "N": int(neutron),
                "EBT": float(ebt),
                "MAA": np.nan if maa == "--" else float(maa),
            }
        )

    if not rows:
        raise ValueError(f"No theory rows found in {path}")

    theory = pd.DataFrame(rows)
    duplicate = theory.duplicated(["A", "Z", "N"], keep=False)
    if duplicate.any():
        keys = theory.loc[duplicate, ["A", "Z", "N"]].drop_duplicates()
        raise ValueError(f"Duplicate theory keys in {path}:\n{keys.to_string(index=False)}")

    return theory


def match_experiment_to_theory(
    experiment: pd.DataFrame, theory: pd.DataFrame
) -> pd.DataFrame:
    """Keep theory values only at nuclei present in the experimental data."""
    experiment = experiment.copy()
    experiment["A"] = pd.to_numeric(experiment["A"], errors="raise")

    matched = experiment.merge(
        theory,
        on=["A", "Z", "N"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    missing_keys = matched.loc[matched["_merge"] != "both", ["A", "Z", "N"]]
    if not missing_keys.empty:
        raise ValueError(
            "Experimental nuclei missing from the theory table:\n"
            f"{missing_keys.to_string(index=False)}"
        )

    missing_values = matched[["EBT", "MAA"]].isna().any(axis=1)
    if missing_values.any():
        keys = matched.loc[missing_values, ["A", "Z", "N", "EBT", "MAA"]]
        raise ValueError(
            "Matched nuclei with unavailable theory values:\n"
            f"{keys.to_string(index=False)}"
        )

    return matched.drop(columns="_merge").sort_values(
        ["Z", "N", "A"], kind="stable"
    )


def plot_theory_chain(
    ax: plt.Axes, data: pd.DataFrame, color: str, model: str
) -> None:
    """Plot one prediction only at the matched experimental neutron numbers."""
    style = THEORY_STYLES[model]
    marker_kwargs: dict[str, str | float] = {
        "marker": style["marker"],
        "markersize": style["markersize"],
        "markeredgecolor": color,
        "markeredgewidth": 0.9,
    }
    if style["marker"] == "D":
        marker_kwargs["markerfacecolor"] = "white"

    ax.plot(
        data["N"],
        data[model],
        color=color,
        linestyle=style["linestyle"],
        linewidth=1.05,
        zorder=2.1,
        **marker_kwargs,
    )


def add_242fm_experimental_connections(ax: plt.Axes, fm_data: pd.DataFrame) -> None:
    """Connect both measured 242Fm values to the measured Fm systematics."""
    anchor_points = fm_data[fm_data["N"] > 142].sort_values("N", kind="stable")
    if anchor_points.empty:
        return

    anchor = anchor_points.iloc[0]
    ax.plot(
        [142, anchor["N"]],
        [np.log10(4.0e-6), anchor["logTsf"]],
        "-",
        color="red",
        linewidth=1.35,
        zorder=2.4,
    )

    tabulated = fm_data[fm_data["N"] == 142]
    if tabulated.empty:
        return

    row = tabulated.iloc[0]
    ax.plot(
        [row["N"], anchor["N"]],
        [row["logTsf"], anchor["logTsf"]],
        "-",
        color=ELEMENTS[100]["color"],
        linewidth=1.35,
        zorder=2.3,
    )


def emphasize_242fm_maa_prediction(ax: plt.Axes, fm_data: pd.DataFrame) -> None:
    """Redraw the nearby MAA marker above the large 4-us experimental star."""
    prediction = fm_data[fm_data["N"] == 142]
    if prediction.empty:
        return

    row = prediction.iloc[0]
    style = THEORY_STYLES["MAA"]
    ax.plot(
        [row["N"]],
        [row["MAA"]],
        linestyle="none",
        marker=style["marker"],
        markersize=4.2,
        markerfacecolor="white",
        markeredgecolor=ELEMENTS[100]["color"],
        markeredgewidth=1.0,
        zorder=6.5,
    )


def make_source_legend(star_handle: Line2D) -> list[Line2D]:
    """Create neutral handles: color denotes element and style denotes source."""
    handles = [
        Line2D(
            [0],
            [0],
            color="0.15",
            linestyle="-",
            linewidth=1.35,
            label="Experiment",
        )
    ]

    for model, style in THEORY_STYLES.items():
        marker_face = "white" if style["marker"] == "D" else "0.15"
        handles.append(
            Line2D(
                [0],
                [0],
                color="0.15",
                linestyle=style["linestyle"],
                linewidth=1.05,
                marker=style["marker"],
                markersize=style["markersize"],
                markerfacecolor=marker_face,
                markeredgecolor="0.15",
                markeredgewidth=0.9,
                label=style["label"],
            )
        )

    star_handle.set_label(r"$^{242}\mathrm{Fm}$ (4 $\mu$s)")
    handles.append(star_handle)
    return handles


def draw_figure(data: pd.DataFrame) -> tuple[plt.Figure, dict[str, int]]:
    """Draw experiment, ETB, and MAA on the original systematics axes."""
    apply_origin_publication_style()
    fig, ax = plt.subplots(figsize=(5.1, 3.65), constrained_layout=True)

    counts: dict[str, int] = {}
    fm_data: pd.DataFrame | None = None
    for proton, spec in ELEMENTS.items():
        element = data[data["Z"] == proton].sort_values("N", kind="stable")
        counts[spec["symbol"]] = len(element)

        for model in THEORY_STYLES:
            plot_theory_chain(ax, element, spec["color"], model)
        plot_element(ax, element, spec)

        if proton == 100:
            fm_data = element

    if fm_data is not None:
        add_242fm_experimental_connections(ax, fm_data)
    star_handle = add_242fm_star(ax)
    if fm_data is not None:
        emphasize_242fm_maa_prediction(ax, fm_data)
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
    for x, y, label, color in element_labels:
        ax.text(
            x,
            y,
            label,
            color=color,
            fontsize=13.0,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5},
            zorder=7,
        )

    ax.set_xlabel("Neutron number")
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

    ax.legend(
        handles=make_source_legend(star_handle),
        loc="upper left",
        ncols=1,
        handlelength=2.4,
        borderaxespad=0.65,
        labelspacing=0.35,
    )

    return fig, counts


def save_figure(fig: plt.Figure, out_stem: Path) -> list[Path]:
    """Save publication vector and high-resolution raster formats."""
    paths = []
    for suffix in FIGURE_FORMATS:
        path = out_stem.with_suffix(f".{suffix}")
        fig.savefig(path)
        paths.append(path)
    return paths


def main() -> None:
    experiment = load_even_neutron_ground_state_data(DATA_PATH)
    theory = load_theory_data(THEORY_PATH)
    matched = match_experiment_to_theory(experiment, theory)

    fig, counts = draw_figure(matched)
    paths = save_figure(fig, OUT_STEM)
    plt.close(fig)

    print("Matched and plotted nuclei:")
    for key in ("Cf", "Fm", "No"):
        print(f"  {key}: {counts[key]}")
    print(f"  Total: {len(matched)}")
    print(f"  ETB predictions: {matched['EBT'].notna().sum()}")
    print(f"  MAA predictions: {matched['MAA'].notna().sum()}")
    print(f"  242Fm (4 us) star: {counts['242Fm_star']}")
    for path in paths:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
