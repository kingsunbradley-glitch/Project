#!/Users/evalie/miniconda3/envs/py_env/bin/python
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(
        Path(os.environ.get("TMPDIR", "/private/tmp"))
        / "matplotlib-242fm"
    ),
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import chi2


def strip_comment(line: str) -> str:
    """Remove comments beginning with # and surrounding spaces."""
    return line.split("#", 1)[0].strip()


def split_fields(line: str) -> list[str]:
    """Split a comma-separated or whitespace-separated line."""
    if "," in line:
        return [
            field.strip()
            for field in next(csv.reader([line]))
        ]

    return re.split(r"\s+", line.strip())


def read_delta_t(
    path: Path,
    column: str = "Delta_t1",
) -> np.ndarray:
    """
    Read positive decay-time values from a text file.

    The first non-comment line must contain the column names.
    Delta_t1 is assumed to be in microseconds.
    """

    lines = [
        strip_comment(line)
        for line in path.read_text(
            encoding="utf-8"
        ).splitlines()
    ]

    lines = [
        line
        for line in lines
        if line
    ]

    if not lines:
        raise ValueError(
            f"{path} is empty."
        )

    header = split_fields(lines[0])

    if column not in header:
        raise ValueError(
            f"Column {column!r} was not found in {path}."
        )

    column_index = header.index(column)
    values: list[float] = []

    for line_number, line in enumerate(
        lines[1:],
        start=2,
    ):
        fields = split_fields(line)

        if len(fields) != len(header):
            raise ValueError(
                f"Line {line_number}: expected "
                f"{len(header)} columns, "
                f"got {len(fields)}."
            )

        try:
            value = float(
                fields[column_index]
            )
        except ValueError as exc:
            raise ValueError(
                f"Line {line_number}: "
                f"{column} must be numeric, "
                f"got {fields[column_index]!r}."
            ) from exc

        if value <= 0.0:
            raise ValueError(
                f"Line {line_number}: "
                f"{column} must be positive, "
                f"got {value}."
            )

        values.append(value)

    if not values:
        raise ValueError(
            f"No {column} values found in {path}."
        )

    return np.asarray(
        values,
        dtype=float,
    )


def theoretical_curve_log10(
    t: np.ndarray,
    lambda_est: float,
) -> np.ndarray:
    r"""
    Probability density in x = log10(t).

    The exponential decay-time distribution is

        f(t) = lambda * exp(-lambda * t).

    With

        x = log10(t),

    one has

        dt/dx = ln(10) * t.

    Therefore,

        dP/dx
        = ln(10) * lambda * t * exp(-lambda * t).

    The returned curve is a probability density per unit
    log10(t), and its integral over the full x range is 1.
    """

    return (
        math.log(10.0)
        * lambda_est
        * t
        * np.exp(-lambda_est * t)
    )


def half_life_interval(
    times_us: np.ndarray,
    confidence: float,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
]:
    """
    Calculate the maximum-likelihood half-life and its exact
    central confidence interval for exponential decay times.
    """

    n_events = times_us.size
    total_time = float(
        np.sum(times_us)
    )

    mean_life = (
        total_time / n_events
    )

    half_life = (
        math.log(2.0) * mean_life
    )

    alpha = 1.0 - confidence
    degrees_of_freedom = 2 * n_events

    mean_life_low = (
        2.0
        * total_time
        / chi2.ppf(
            1.0 - alpha / 2.0,
            degrees_of_freedom,
        )
    )

    mean_life_high = (
        2.0
        * total_time
        / chi2.ppf(
            alpha / 2.0,
            degrees_of_freedom,
        )
    )

    half_life_low = (
        math.log(2.0)
        * mean_life_low
    )

    half_life_high = (
        math.log(2.0)
        * mean_life_high
    )

    err_low = (
        half_life
        - half_life_low
    )

    err_high = (
        half_life_high
        - half_life
    )

    return (
        half_life,
        err_low,
        err_high,
        half_life_low,
        half_life_high,
    )


def plot_delta_t(
    times_us: np.ndarray,
    output_path: Path,
    bin_width: float,
    log_min_arg: float | None,
    log_max_arg: float | None,
    half_life: float,
    err_low: float,
    err_high: float,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
]:
    """
    Plot the logarithmic decay-time histogram together with
    the Schmidt exponential distribution.
    """

    # ---------------------------------------------------------
    # Plot style
    # ---------------------------------------------------------

    plt.rcParams["font.family"] = "serif"

    plt.rcParams["font.serif"] = [
        "STIXGeneral",
        "Times New Roman",
        "serif",
    ]

    plt.rcParams["mathtext.fontset"] = "stix"

    plt.rcParams["axes.linewidth"] = 1.0

    plt.rcParams["xtick.major.width"] = 1.0
    plt.rcParams["ytick.major.width"] = 1.0

    plt.rcParams["xtick.minor.width"] = 0.5
    plt.rcParams["ytick.minor.width"] = 0.5

    plt.rcParams["xtick.major.size"] = 3.5
    plt.rcParams["ytick.major.size"] = 3.5

    plt.rcParams["xtick.minor.size"] = 2.0
    plt.rcParams["ytick.minor.size"] = 2.0

    # ---------------------------------------------------------
    # Lifetime parameters
    # ---------------------------------------------------------

    n_events = times_us.size

    mean_life = float(
        np.mean(times_us)
    )

    lambda_est = (
        1.0 / mean_life
    )

    # In logarithmic time coordinates, the Schmidt curve
    # reaches its maximum at t = tau.
    peak_log = math.log10(
        mean_life
    )

    data_log = np.log10(
        times_us
    )

    # ---------------------------------------------------------
    # Automatic plotting range
    # ---------------------------------------------------------

    log_min = (
        math.floor(
            min(
                peak_log - 2.0,
                float(data_log.min()) - 0.5,
            )
            * 10.0
        )
        / 10.0
    )

    log_max = (
        math.ceil(
            max(
                peak_log + 2.0,
                float(data_log.max()) + 0.5,
            )
            * 10.0
        )
        / 10.0
    )

    if log_min_arg is not None:
        log_min = log_min_arg

    if log_max_arg is not None:
        log_max = log_max_arg

    if (
        log_min_arg is None
        and log_max_arg is None
    ):
        log_min -= 0.5
        log_max -= 0.5

    if log_max <= log_min:
        raise ValueError(
            "log-max must be larger than log-min."
        )

    if bin_width <= 0.0:
        raise ValueError(
            "bin-width must be positive."
        )

    # ---------------------------------------------------------
    # Equal-width histogram bins in log10(t)
    # ---------------------------------------------------------

    number_of_bins = int(
        math.ceil(
            (log_max - log_min)
            / bin_width
        )
    )

    bins = (
        log_min
        + np.arange(
            number_of_bins + 1,
            dtype=float,
        )
        * bin_width
    )

    # ---------------------------------------------------------
    # Draw histogram
    # ---------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(6, 5)
    )

    counts, bin_edges, _ = ax.hist(
        data_log,
        bins=bins,
        color="#D81525",
        edgecolor="black",
        linewidth=0.3,
    )

    # ---------------------------------------------------------
    # Schmidt theoretical curve
    #
    # x = log10(t)
    #
    # Probability density:
    #
    # dP/dx
    # = ln(10) * lambda * t * exp(-lambda*t)
    #
    # Expected counts per histogram bin:
    #
    # dN_bin/dx
    # = N * Delta_x * dP/dx
    #
    # The factor Delta_x is required because the histogram
    # vertical axis is Counts per bin rather than density.
    # ---------------------------------------------------------

    x_axis = np.linspace(
        log_min,
        log_max,
        2000,
    )

    t_real = (
        10.0**x_axis
    )

    probability_density = (
        theoretical_curve_log10(
            t_real,
            lambda_est,
        )
    )

    expected_counts_per_bin = (
        n_events
        * bin_width
        * probability_density
    )

    ax.plot(
        x_axis,
        expected_counts_per_bin,
        "-",
        linewidth=2.0,
        color="#6A5ACD",
    )

    # ---------------------------------------------------------
    # Axes
    # ---------------------------------------------------------

    ax.set_xlim(
        log_min,
        log_max,
    )

    ax.set_xlabel(
        r"$\log_{10}[t\;(\mathrm{\mu s})]$",
        fontsize=18,
    )

    ax.set_ylabel(
        "Counts",
        fontsize=18,
    )

    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(1.0)
    )

    ax.xaxis.set_minor_locator(
        ticker.MultipleLocator(0.5)
    )

    ax.yaxis.set_major_locator(
        ticker.MaxNLocator(
            integer=True
        )
    )

    ax.tick_params(
        direction="in",
        top=True,
        right=True,
        which="both",
        labelsize=14,
    )

    for spine in ax.spines.values():
        spine.set_visible(True)

    # ---------------------------------------------------------
    # Labels
    # ---------------------------------------------------------

    ax.text(
        0.05,
        0.90,
        r"$^{242}$Fm",
        transform=ax.transAxes,
        fontsize=22,
        ha="left",
        va="top",
    )

    ax.text(
        0.05,
        0.78,
        (
            rf"$\mathrm{{T}}_{{1/2}}^{{\mathrm{{SF}}}}"
            rf" = {half_life:.2f}"
            rf"^{{+{err_high:.2f}}}"
            rf"_{{-{err_low:.2f}}}"
            rf"\ \mathrm{{\mu s}}$"
        ),
        transform=ax.transAxes,
        fontsize=16,
        ha="left",
        va="top",
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(fig)

    # ---------------------------------------------------------
    # Diagnostic normalization values
    # ---------------------------------------------------------

    histogram_area = float(
        np.sum(
            counts
            * np.diff(bin_edges)
        )
    )

    full_curve_area = float(
        n_events
        * bin_width
    )

    displayed_curve_area = float(
        np.trapezoid(
            expected_counts_per_bin,
            x_axis,
        )
    )

    print(
        f"Histogram geometric area in displayed bins "
        f"= {histogram_area:.6f}"
    )

    print(
        f"Theoretical full-range area N*bin_width "
        f"= {full_curve_area:.6f}"
    )

    print(
        f"Theoretical area inside displayed x-range "
        f"= {displayed_curve_area:.6f}"
    )

    return (
        lambda_est,
        peak_log,
        mean_life,
        log_min,
        log_max,
    )


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description=(
            "Plot Delta_t1 data in logarithmic time "
            "coordinates using the Schmidt exponential "
            "distribution."
        )
    )

    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=here / "input.txt",
        help=(
            "Input text file. "
            "Default: input.txt in the script directory."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=(
            here
            / "delta_t_distribution.pdf"
        ),
        help=(
            "Output figure path. "
            "Default: delta_t_distribution.pdf"
        ),
    )

    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.2,
        help=(
            "Histogram-bin width in log10(t). "
            "Default: 0.2."
        ),
    )

    parser.add_argument(
        "--log-min",
        type=float,
        default=None,
        help=(
            "Optional lower x-axis limit "
            "in log10(us)."
        ),
    )

    parser.add_argument(
        "--log-max",
        type=float,
        default=None,
        help=(
            "Optional upper x-axis limit "
            "in log10(us)."
        ),
    )

    parser.add_argument(
        "--confidence",
        type=float,
        default=0.6827,
        help=(
            "Two-sided confidence level. "
            "Default: 0.6827."
        ),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        input_path = (
            args.input.resolve()
        )

        output_path = (
            args.output.resolve()
        )

        times_us = read_delta_t(
            input_path
        )

        if not (
            0.0
            < args.confidence
            < 1.0
        ):
            raise ValueError(
                "--confidence must be "
                "between 0 and 1."
            )

        if args.bin_width <= 0.0:
            raise ValueError(
                "--bin-width must be positive."
            )

        (
            half_life,
            err_low,
            err_high,
            interval_low,
            interval_high,
        ) = half_life_interval(
            times_us,
            args.confidence,
        )

        (
            lambda_est,
            peak_log,
            mean_life,
            log_min,
            log_max,
        ) = plot_delta_t(
            times_us=times_us,
            output_path=output_path,
            bin_width=args.bin_width,
            log_min_arg=args.log_min,
            log_max_arg=args.log_max,
            half_life=half_life,
            err_low=err_low,
            err_high=err_high,
        )

    except (
        ValueError,
        OSError,
    ) as exc:
        print(
            f"Error: {exc}",
            file=sys.stderr,
        )

        return 1

    print()

    print(
        f"Read {times_us.size} "
        f"Delta_t1 value(s) "
        f"from {input_path}."
    )

    print(
        "Delta_t1 values are treated "
        "as microseconds (us)."
    )

    print(
        "Values (us): "
        + ", ".join(
            f"{value:g}"
            for value in times_us
        )
    )

    print(
        f"Mean lifetime tau "
        f"= {mean_life:.3f} us."
    )

    print(
        f"Decay constant lambda "
        f"= {lambda_est:.6f} us^-1."
    )

    print(
        f"Log-time curve peaks at "
        f"log10(t/us) = {peak_log:.3f}."
    )

    print(
        f"Histogram bin width "
        f"= {args.bin_width:.3f} "
        f"in log10(t)."
    )

    print(
        f"Plot x-range: "
        f"log10(t/us) = "
        f"[{log_min:.3f}, "
        f"{log_max:.3f}]."
    )

    print(
        f"T_1/2 = "
        f"{half_life:.3f} "
        f"+{err_high:.3f} "
        f"-{err_low:.3f} us "
        f"({args.confidence * 100.0:.2f}% C.L.)"
    )

    print(
        f"Confidence interval: "
        f"[{interval_low:.3f}, "
        f"{interval_high:.3f}] us"
    )

    print(
        f"Wrote {output_path}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )