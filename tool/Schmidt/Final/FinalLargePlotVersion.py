"""Draw the 5 x 2 energy/decay-time panels for the five decay-chain members.

The event values are transcribed from Table B supplied with this script.  The
energy panels contain only explicitly reported alpha-particle energies;
parenthesized energy deposits and ``SF`` entries are therefore omitted.  Their
positive decay times are still included in the decay-time panels.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


# -----------------------------------------------------------------------------
# Global style
# -----------------------------------------------------------------------------
LEFT_PANEL_LABELS = ["(a)", "(b)", "(c)", "(d)", "(e)"]
RIGHT_PANEL_LABELS = ["(f)", "(g)", "(h)", "(i)", "(j)"]
FONT_WEIGHT = "bold"
FONT_SIZE_NUCLEUS = 44
FONT_SIZE_TICK = 30
FONT_SIZE_AXIS = 32
FONT_SIZE_LEGEND = 26
FONT_SIZE_GLOBAL = 34

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "serif"],
        "mathtext.fontset": "stix",
        "font.weight": FONT_WEIGHT,
        "axes.labelweight": FONT_WEIGHT,
        "axes.titleweight": FONT_WEIGHT,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
    }
)


INSTITUTES = ["IMP", "JINR", "RIKEN", "GSI"]
COLORS = ["#D81525", "blue", "green", "#B8860B"]


def seconds(values):
    """Convert the table's seconds to the common plotting unit, milliseconds."""
    return np.asarray(values, dtype=float) * 1000.0


def dataset(name, energy, time, fit_max_ms=None):
    """Build a dataset, optionally splitting its curve fits at a time cutoff."""
    energy_arrays = [np.asarray(values, dtype=float) for values in energy]
    time_arrays = [np.asarray(values, dtype=float) for values in time]
    all_positive_times = np.concatenate(
        [values[values > 0] for values in time_arrays]
    )
    if fit_max_ms is None:
        fit_data = all_positive_times
        second_fit_data = np.array([])
    else:
        fit_data = all_positive_times[all_positive_times < fit_max_ms]
        second_fit_data = all_positive_times[all_positive_times > fit_max_ms]
    return {
        "name": name,
        "energy": energy_arrays,
        "time": time_arrays,
        "fit": fit_data,
        "second_fit": second_fit_data,
        "fit_cut_ms": fit_max_ms,
    }


# Data order in every list: IMP, JINR, RIKEN, GSI.
# Cn and Ds times are tabulated in ms; Hs, Sg, and Rf times are converted from s.
ALL_DATASETS = [
    dataset(
        r"$^{277}$Cn",
        energy=[[], [], [11.09, 11.32, 11.07], [11.45, 11.17]],
        time=[[], [], [1.10, 1.22, 0.370], [0.280, 1.406]],
    ),
    dataset(
        r"$^{273}$Ds",
        energy=[
            [11.165, 10.858],
            [11.017, 10.929],
            [11.14, 11.15, 11.03],
            [11.08, 11.20],
        ],
        time=[
            [0.108, 8.320],
            [0.184, 41.703],
            [0.520, 0.0399, 0.373],
            [0.110, 0.310],
        ],
        fit_max_ms=1.0,  # Fit the decay curve with only the seven t < 1 ms events.
    ),
    dataset(
        r"$^{269}$Hs",
        energy=[
            [9.103, 8.916],
            [],  # The JINR value (0.342) is a parenthesized energy deposit.
            [9.17, 9.25, 9.15],
            [9.23, 9.18],
        ],
        time=[
            seconds([5.006, 9.363]),
            seconds([4.0978]),
            seconds([14.2, 0.270, 36.0]),
            seconds([19.7, 22.0]),
        ],
    ),
    dataset(
        r"$^{265}$Sg",
        energy=[
            [8.660, 8.690],
            [8.397, 8.509],
            [8.71, 8.70, 8.66],
            [],  # GSI reports only parenthesized deposits (4.60) and (0.2).
        ],
        time=[
            seconds([6.340, 11.663]),
            seconds([9.0766, 58.19]),
            seconds([23.0, 79.0, 13.8]),
            seconds([7.4, 18.8]),
        ],
    ),
    dataset(
        r"$^{261}$Rf",
        energy=[[], [], [], [8.52]],  # All other Rf decays are reported as SF.
        time=[
            seconds([5.298, 7.241]),
            seconds([3.0344, 0.3923]),
            seconds([2.97, 8.30, 3.73]),
            seconds([4.7, 14.5]),
        ],
    ),
]


def theoretical_curve(time_ms, decay_constant):
    """Exponential-decay density expressed per logarithmic time interval."""
    return decay_constant * time_ms * np.exp(-decay_constant * time_ms)


def draw_figure():
    """Create and return the five-row energy/time figure."""
    figure = plt.figure(figsize=(17, 20))
    grid = GridSpec(
        len(ALL_DATASETS),
        2,
        figure=figure,
        width_ratios=[1.5, 1],
        wspace=0,
        hspace=0,
    )

    energy_min, energy_max = 7.8, 11.8
    log_time_min, log_time_max = -3.4, 6.2
    energy_bins = np.linspace(energy_min, energy_max, 201)  # 20-keV bins
    log_time_bins = np.linspace(log_time_min, log_time_max, 90)

    first_energy_axis = None

    for row, data in enumerate(ALL_DATASETS):
        if row == 0:
            energy_axis = figure.add_subplot(grid[row, 0])
            time_axis = figure.add_subplot(grid[row, 1], sharey=energy_axis)
            first_energy_axis = energy_axis
        else:
            energy_axis = figure.add_subplot(
                grid[row, 0], sharex=first_energy_axis
            )
            time_axis = figure.add_subplot(grid[row, 1], sharey=energy_axis)

        # Left column: alpha-energy distribution.
        energy_axis.hist(
            data["energy"],
            bins=energy_bins,
            stacked=True,
            color=COLORS,
            edgecolor="none",
        )
        energy_axis.set_xlim(energy_min, energy_max)
        energy_axis.text(
            0.05,
            0.70,
            data["name"],
            transform=energy_axis.transAxes,
            fontsize=FONT_SIZE_NUCLEUS,
            fontweight=FONT_WEIGHT,
        )
        energy_axis.text(
            0.95,
            0.92,
            LEFT_PANEL_LABELS[row],
            transform=energy_axis.transAxes,
            ha="right",
            va="top",
            fontsize=FONT_SIZE_AXIS,
            fontweight=FONT_WEIGHT,
        )
        energy_axis.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        energy_axis.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        energy_axis.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))
        energy_axis.spines["top"].set_visible(True)
        energy_axis.spines["right"].set_visible(False)
        energy_axis.tick_params(
            direction="in",
            top=True,
            right=False,
            which="both",
            labelsize=FONT_SIZE_TICK,
        )

        # Right column: decay-time distribution on a log10(ms) coordinate.
        log_time_data = [
            np.log10(values[values > 0]) if values.size else np.array([])
            for values in data["time"]
        ]
        histogram_counts, _, _ = time_axis.hist(
            log_time_data,
            bins=log_time_bins,
            stacked=True,
            color=COLORS,
            edgecolor="none",
        )
        time_axis.set_xlim(log_time_min, log_time_max)
        time_axis.text(
            0.05,
            0.92,
            RIGHT_PANEL_LABELS[row],
            transform=time_axis.transAxes,
            ha="left",
            va="top",
            fontsize=FONT_SIZE_AXIS,
            fontweight=FONT_WEIGHT,
        )

        log_time_axis = np.linspace(log_time_min, log_time_max, 500)
        time_ms = 10**log_time_axis
        maximum_count = max(float(np.max(histogram_counts)), 1.0)

        curve_groups = [(data["fit"], "-", None)]
        if data["second_fit"].size:
            cutoff = data["fit_cut_ms"]
            curve_groups = [
                (data["fit"], "-", rf"$t<{cutoff:g}\,\mathrm{{ms}}$ fit"),
                (
                    data["second_fit"],
                    "--",
                    rf"$t>{cutoff:g}\,\mathrm{{ms}}$ fit",
                ),
            ]

        for fit_values, line_style, curve_label in curve_groups:
            mean_time_ms = np.mean(fit_values)
            decay_constant = 1.0 / mean_time_ms
            curve = theoretical_curve(time_ms, decay_constant)
            if np.max(curve) > 0:
                time_axis.plot(
                    log_time_axis,
                    curve * maximum_count / np.max(curve),
                    color="#6A5ACD",
                    linestyle=line_style,
                    linewidth=2,
                    label=curve_label,
                )

        time_axis.xaxis.set_major_locator(ticker.MultipleLocator(2))
        time_axis.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        time_axis.spines["top"].set_visible(True)
        time_axis.spines["left"].set_linestyle("--")
        time_axis.spines["left"].set_linewidth(3)
        time_axis.spines["left"].set_color("gray")
        time_axis.tick_params(
            direction="in",
            top=True,
            right=True,
            which="both",
            labelsize=FONT_SIZE_TICK,
            labelright=True,
        )
        time_axis.tick_params(axis="y", labelleft=False, left=False)

        for tick_label in energy_axis.get_xticklabels() + energy_axis.get_yticklabels():
            tick_label.set_fontweight(FONT_WEIGHT)
        for tick_label in time_axis.get_xticklabels() + time_axis.get_yticklabels():
            tick_label.set_fontweight(FONT_WEIGHT)

        if row < len(ALL_DATASETS) - 1:
            energy_axis.tick_params(labelbottom=False)
            time_axis.tick_params(labelbottom=False)
        else:
            energy_axis.set_xlabel(
                r"$E_{\alpha}$ (MeV)",
                fontsize=FONT_SIZE_AXIS,
                fontweight=FONT_WEIGHT,
            )
            time_axis.set_xlabel(
                r"$\log_{10}[t\ (\mathrm{ms})]$",
                fontsize=FONT_SIZE_AXIS,
                fontweight=FONT_WEIGHT,
            )

        # The institute-color legend belongs specifically to the 277Cn time panel.
        if row == 0:
            legend_handles = [
                Patch(facecolor=color, edgecolor="none", label=label)
                for label, color in zip(INSTITUTES, COLORS)
            ]
            time_axis.legend(
                handles=legend_handles,
                frameon=False,
                fontsize=FONT_SIZE_LEGEND,
                loc="upper right",
                prop={"weight": FONT_WEIGHT, "size": FONT_SIZE_LEGEND},
            )

        # Energy and time panels in a row share their count scale.
        energy_axis.set_ylim(bottom=0)
        upper_limit = max(energy_axis.get_ylim()[1], 1.0)
        energy_axis.set_ylim(0, upper_limit * 1.6)

    figure.text(
        0.08,
        0.5,
        "Counts / 20 keV",
        va="center",
        rotation="vertical",
        fontsize=FONT_SIZE_GLOBAL,
        fontweight=FONT_WEIGHT,
    )
    figure.text(
        0.97,
        0.5,
        "Counts",
        va="center",
        rotation="vertical",
        fontsize=FONT_SIZE_GLOBAL,
        fontweight=FONT_WEIGHT,
    )
    figure.subplots_adjust(
        left=0.12,
        right=0.95,
        top=0.95,
        bottom=0.05,
        hspace=0,
        wspace=0,
    )
    return figure


def main():
    output_directory = Path(__file__).resolve().parent
    figure = draw_figure()
    png_path = output_directory / "decay_properties_5x2.png"
    pdf_path = output_directory / "decay_properties_5x2.pdf"
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"Saved {png_path}")
    print(f"Saved {pdf_path}")


if __name__ == "__main__":
    main()
