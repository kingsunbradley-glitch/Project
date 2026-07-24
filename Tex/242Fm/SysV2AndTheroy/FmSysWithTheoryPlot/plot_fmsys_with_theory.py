#!/home/evalie/miniconda3/envs/py_env/bin/python
# -*- coding: utf-8 -*-

"""Plot experimental and theoretical SF half-life systematics.

Run with the py_env interpreter:
    /home/evalie/miniconda3/envs/py_env/bin/python plot_fmsys_with_theory.py
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "matplotlib-fmsys-with-theory"),
)

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_PATH = SCRIPT_DIR / "FmSysWithTheory.dat"
OUTPUT_STEM = SCRIPT_DIR / "FmSys_Exp_BSkG3_EBT_MAA_Origin"
OUTPUT_FORMATS = ("pdf", "svg", "png")

ELEMENTS = {
    98: {"symbol": "Cf", "color": "#D55E00"},
    100: {"symbol": "Fm", "color": "#3657D6"},
    102: {"symbol": "No", "color": "#009E73"},
}

SOURCES = {
    "Experiment": {
        "column": "T1/2_SF",
        "linestyle": "-",
        "linewidth": 1.35,
        "marker": "o",
        "markersize": 5.0,
        "markerfacecolor": "element",
        "zorder": 4.0,
    },
    "BSkG3": {
        "column": "BSkG3_s",
        "linestyle": (0, (5.0, 2.0)),
        "linewidth": 1.05,
        "marker": "s",
        "markersize": 4.0,
        "markerfacecolor": "white",
        "zorder": 2.8,
    },
    "EBT": {
        "column": "EBT_s",
        "linestyle": (0, (5.0, 1.7, 1.2, 1.7)),
        "linewidth": 1.05,
        "marker": "x",
        "markersize": 4.6,
        "markerfacecolor": "element",
        "zorder": 2.6,
    },
    "MAA": {
        "column": "MAA_s",
        "linestyle": (0, (1.2, 1.8)),
        "linewidth": 1.05,
        "marker": "D",
        "markersize": 3.8,
        "markerfacecolor": "white",
        "zorder": 2.4,
    },
}


def apply_origin_style() -> None:
    """Apply an Origin-like boxed-axis publication style."""
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 9.0,
            "axes.labelsize": 11.0,
            "axes.linewidth": 1.15,
            "axes.edgecolor": "black",
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
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
            "legend.fontsize": 8.2,
            "legend.frameon": False,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_data(path: Path) -> tuple[pd.DataFrame, int]:
    """Load the merged table and retain even-even ground-state nuclei."""
    data = pd.read_csv(
        path,
        sep=r"\s+",
        dtype={"A": "string", "T1/2_SF_limit": "string"},
        na_values=["--"],
        keep_default_na=True,
    )
    required = {
        "A",
        "Z",
        "N",
        "T1/2_SF",
        "T1/2_SF_limit",
        "BSkG3_s",
        "EBT_s",
        "MAA_s",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")

    is_ground_state = data["A"].str.fullmatch(r"\d+", na=False)
    proton_number = pd.to_numeric(data["Z"], errors="raise")
    neutron_number = pd.to_numeric(data["N"], errors="raise")
    is_even_even = (
        is_ground_state & (proton_number % 2 == 0) & (neutron_number % 2 == 0)
    )
    omitted_rows = int((~is_even_even).sum())
    data = data.loc[is_even_even].copy()

    numeric_columns = [
        "A",
        "Z",
        "N",
        "T1/2_SF",
        "BSkG3_s",
        "EBT_s",
        "MAA_s",
    ]
    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="raise")

    time_columns = ["T1/2_SF", "BSkG3_s", "EBT_s", "MAA_s"]
    for column in time_columns:
        invalid = data[column].notna() & (data[column] <= 0.0)
        if invalid.any():
            nuclei = data.loc[invalid, ["A", "Z", "N"]]
            raise ValueError(
                f"Non-positive values in {column}:\n{nuclei.to_string(index=False)}"
            )
        data[f"log_{column}"] = np.log10(data[column])

    return data.sort_values(["Z", "N", "A"], kind="stable"), omitted_rows


def add_limit_arrows(ax: plt.Axes, chain: pd.DataFrame, color: str) -> None:
    """Draw downward upper-limit arrows and upward lower-limit arrows."""
    arrow_span = 0.72
    for _, row in chain.iterrows():
        limit = str(row["T1/2_SF_limit"]).strip().lower()
        y_value = float(row["log_T1/2_SF"])
        if limit == "upper":
            endpoint = y_value - arrow_span
        elif limit == "lower":
            endpoint = y_value + arrow_span
        else:
            continue

        ax.annotate(
            "",
            xy=(row["N"], endpoint),
            xytext=(row["N"], y_value),
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "lw": 0.9,
                "mutation_scale": 8.0,
                "shrinkA": 4.0,
                "shrinkB": 0.0,
            },
            zorder=5,
        )


def plot_source(
    ax: plt.Axes,
    chain: pd.DataFrame,
    color: str,
    source: str,
    style: dict[str, object],
) -> None:
    """Plot one experimental or theoretical curve for an isotope chain."""
    log_column = f"log_{style['column']}"
    marker_face = color if style["markerfacecolor"] == "element" else "white"
    ax.plot(
        chain["N"],
        chain[log_column],
        color=color,
        linestyle=style["linestyle"],
        linewidth=style["linewidth"],
        marker=style["marker"],
        markersize=style["markersize"],
        markerfacecolor=marker_face,
        markeredgecolor=color,
        markeredgewidth=0.85,
        zorder=style["zorder"],
    )
    if source == "Experiment":
        add_limit_arrows(ax, chain, color)


def add_242fm_star(ax: plt.Axes, fm_chain: pd.DataFrame) -> Line2D:
    """Add the 4-us 242Fm point and connect it to the even-even Fm chain."""
    neutron_number = 142
    log_tsf = np.log10(4.0e-6)
    anchors = fm_chain.loc[fm_chain["N"] > neutron_number].sort_values(
        "N", kind="stable"
    )
    if not anchors.empty:
        anchor = anchors.iloc[0]
        ax.plot(
            [neutron_number, anchor["N"]],
            [log_tsf, anchor["log_T1/2_SF"]],
            color="red",
            linestyle="-",
            linewidth=1.35,
            zorder=4.6,
        )

    ax.plot(
        [neutron_number],
        [log_tsf],
        color="red",
        linestyle="none",
        marker="*",
        markersize=12.5,
        markerfacecolor="red",
        markeredgecolor="red",
        markeredgewidth=0.9,
        zorder=6.0,
    )
    return Line2D(
        [0],
        [0],
        color="red",
        linestyle="none",
        marker="*",
        markersize=10.5,
        markerfacecolor="red",
        markeredgecolor="red",
        label=r"$^{242}\mathrm{Fm}$ (4 $\mu$s)",
    )


def emphasize_242fm_maa(ax: plt.Axes, fm_chain: pd.DataFrame) -> None:
    """Redraw the nearby MAA marker above the large red star."""
    prediction = fm_chain.loc[fm_chain["N"] == 142]
    if prediction.empty:
        return

    row = prediction.iloc[0]
    style = SOURCES["MAA"]
    ax.plot(
        [row["N"]],
        [row["log_MAA_s"]],
        color=ELEMENTS[100]["color"],
        linestyle="none",
        marker=style["marker"],
        markersize=4.2,
        markerfacecolor="white",
        markeredgecolor=ELEMENTS[100]["color"],
        markeredgewidth=1.0,
        zorder=6.5,
    )


def make_element_legend() -> list[Line2D]:
    """Return color handles for the three element chains."""
    return [
        Line2D(
            [0],
            [0],
            color=spec["color"],
            linewidth=2.0,
            label=rf"{spec['symbol']} ($Z={proton}$)",
        )
        for proton, spec in ELEMENTS.items()
    ]


def make_source_legend(star_handle: Line2D) -> list[Line2D]:
    """Return neutral handles for experiment and prediction models."""
    handles: list[Line2D] = []
    for source, style in SOURCES.items():
        marker_face = "0.15" if style["markerfacecolor"] == "element" else "white"
        handles.append(
            Line2D(
                [0],
                [0],
                color="0.15",
                linestyle=style["linestyle"],
                linewidth=style["linewidth"],
                marker=style["marker"],
                markersize=style["markersize"],
                markerfacecolor=marker_face,
                markeredgecolor="0.15",
                markeredgewidth=0.85,
                label=source,
            )
        )
    handles.append(star_handle)
    return handles


def draw_figure(data: pd.DataFrame) -> plt.Figure:
    """Create the complete Origin-style systematics figure."""
    apply_origin_style()
    fig, ax = plt.subplots(figsize=(6.6, 4.6), constrained_layout=True)

    fm_chain = pd.DataFrame()
    for proton, element in ELEMENTS.items():
        chain = data.loc[data["Z"] == proton].sort_values("N", kind="stable")
        if chain.empty:
            continue
        for source, style in SOURCES.items():
            plot_source(ax, chain, element["color"], source, style)
        if proton == 100:
            fm_chain = chain

    star_handle = add_242fm_star(ax, fm_chain)
    emphasize_242fm_maa(ax, fm_chain)

    ax.axvline(
        152,
        color="0.42",
        linewidth=0.9,
        linestyle=(0, (4.0, 2.2)),
        zorder=1,
    )
    ax.text(
        152.25,
        -7.15,
        r"$N=152$",
        color="0.30",
        fontsize=9.0,
        ha="left",
        va="bottom",
    )

    ax.set_xlabel(r"Neutron number, $N$")
    ax.set_ylabel(r"$\log_{10}\!\left(T_{1/2}^{\mathrm{SF}}/\mathrm{s}\right)$")
    ax.set_xlim(138, 162)
    ax.set_ylim(-8, 14)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", direction="in", top=True, right=True)

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.15)

    element_legend = ax.legend(
        handles=make_element_legend(),
        loc="upper left",
        handlelength=2.0,
        borderaxespad=0.65,
        labelspacing=0.35,
    )
    ax.add_artist(element_legend)
    ax.legend(
        handles=make_source_legend(star_handle),
        loc="lower left",
        ncols=3,
        columnspacing=1.0,
        handlelength=2.5,
        borderaxespad=0.65,
        labelspacing=0.35,
    )

    return fig


def save_figure(fig: plt.Figure) -> list[Path]:
    """Save vector versions and a high-resolution raster version."""
    outputs: list[Path] = []
    for suffix in OUTPUT_FORMATS:
        path = OUTPUT_STEM.with_suffix(f".{suffix}")
        fig.savefig(path)
        outputs.append(path)
    return outputs


def main() -> None:
    data, omitted_rows = load_data(DATA_PATH)
    figure = draw_figure(data)
    outputs = save_figure(figure)
    plt.close(figure)

    print(f"Input: {DATA_PATH}")
    print(f"Even-even ground-state nuclei plotted: {len(data)}")
    print(f"Non-even-even rows omitted: {omitted_rows}")
    print("242Fm (4 us) star: 1")
    for source, style in SOURCES.items():
        count = int(data[style["column"]].notna().sum())
        print(f"{source} points: {count}")
    for path in outputs:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
