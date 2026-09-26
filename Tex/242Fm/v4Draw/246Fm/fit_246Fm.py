#!/Users/evalie/miniconda3/envs/py_env/bin/python
"""Fit the local 246Fm decay-time and energy data table."""

from __future__ import annotations

import argparse
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(os.environ.get("TMPDIR", "/private/tmp")) / "matplotlib-246fm"),
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2


@dataclass(frozen=True)
class FitResults:
    event_count: int
    mean_lifetime_s: float
    decay_constant_hz: float
    half_life_s: float
    half_life_low_s: float
    half_life_high_s: float
    energy_gaussian_mean_kev: float
    energy_gaussian_mean_error_kev: float
    energy_gaussian_sigma_kev: float
    energy_gaussian_sigma_error_kev: float
    energy_gaussian_amplitude: float
    energy_arithmetic_mean_kev: float


def read_root_scan(
    path: Path, input_time_unit: str = "auto"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read Row, SumE and DeltaT from ROOT Scan or whitespace tables."""
    scan_pattern = re.compile(
        r"^\s*\*\s*(\d+)\s*\*\s*"
        r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*\*\s*"
        r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*\*\s*$"
    )
    table_pattern = re.compile(
        r"^\s*(\d+)\s+"
        r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*$"
    )
    rows: list[int] = []
    energies_kev: list[float] = []
    delta_t_ns: list[float] = []

    file_text = path.read_text(encoding="utf-8")
    for line in file_text.splitlines():
        match = scan_pattern.match(line) or table_pattern.match(line)
        if match is None:
            continue
        rows.append(int(match.group(1)))
        energies_kev.append(float(match.group(2)))
        delta_t_ns.append(float(match.group(3)))

    if not rows:
        raise ValueError(f"No numeric ROOT Scan rows found in {path}.")

    if input_time_unit == "auto":
        unit_matches = re.findall(r"DeltaT_(ns|us|ms|s)\b", file_text, flags=re.IGNORECASE)
        units = {match.lower() for match in unit_matches}
        if len(units) != 1:
            raise ValueError(
                "Could not determine one input time unit from the DeltaT header; "
                "use --input-time-unit."
            )
        input_time_unit = units.pop()

    unit_to_seconds = {"ns": 1.0e-9, "us": 1.0e-6, "ms": 1.0e-3, "s": 1.0}
    energy = np.asarray(energies_kev, dtype=float)
    times_s = np.asarray(delta_t_ns, dtype=float) * unit_to_seconds[input_time_unit]
    if np.any(times_s <= 0.0):
        raise ValueError("All DeltaT values must be positive.")
    return np.asarray(rows, dtype=int), energy, times_s, input_time_unit


def fit_values(energies_kev: np.ndarray, times_s: np.ndarray, confidence: float) -> FitResults:
    """Calculate the exponential MLE and a Gaussian histogram fit."""
    n = times_s.size

    # Exponential MLE: tau_hat = sum(t_i) / N and T_1/2 = ln(2) tau_hat.
    total_time_s = float(np.sum(times_s))
    mean_lifetime_s = total_time_s / n
    half_life_s = math.log(2.0) * mean_lifetime_s
    decay_constant_hz = 1.0 / mean_lifetime_s

    # Exact two-sided interval obtained by inverting the chi-square pivot.
    alpha = 1.0 - confidence
    dof = 2 * n
    mean_lifetime_low_s = 2.0 * total_time_s / chi2.ppf(1.0 - alpha / 2.0, dof)
    mean_lifetime_high_s = 2.0 * total_time_s / chi2.ppf(alpha / 2.0, dof)
    half_life_low_s = math.log(2.0) * mean_lifetime_low_s
    half_life_high_s = math.log(2.0) * mean_lifetime_high_s

    # Gaussian chi-square fit to a Freedman-Diaconis histogram.  Poisson count
    # uncertainties reproduce the usual spectrum-fitting treatment.
    edges = np.histogram_bin_edges(energies_kev, bins="fd")
    counts, _ = np.histogram(energies_kev, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0

    def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
        return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

    initial = [float(np.max(counts)), float(np.mean(energies_kev)), float(np.std(energies_kev, ddof=1))]
    parameters, covariance = curve_fit(
        gaussian,
        centers,
        counts,
        p0=initial,
        sigma=np.sqrt(counts),
        absolute_sigma=True,
        bounds=([0.0, float(np.min(energies_kev)), 0.0], [np.inf, float(np.max(energies_kev)), np.inf]),
        maxfev=100_000,
    )
    parameter_errors = np.sqrt(np.diag(covariance))

    return FitResults(
        event_count=n,
        mean_lifetime_s=mean_lifetime_s,
        decay_constant_hz=decay_constant_hz,
        half_life_s=half_life_s,
        half_life_low_s=half_life_low_s,
        half_life_high_s=half_life_high_s,
        energy_gaussian_mean_kev=float(parameters[1]),
        energy_gaussian_mean_error_kev=float(parameter_errors[1]),
        energy_gaussian_sigma_kev=float(parameters[2]),
        energy_gaussian_sigma_error_kev=float(parameter_errors[2]),
        energy_gaussian_amplitude=float(parameters[0]),
        energy_arithmetic_mean_kev=float(np.mean(energies_kev)),
    )


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["STIXGeneral", "Times New Roman", "serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.0,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def plot_reverse_clock(times_s: np.ndarray, result: FitResults, output_stem: Path) -> None:
    """Plot the logarithmic-time histogram and its reverse-clock curve."""
    data_log = np.log10(times_s)
    peak_log = math.log10(result.mean_lifetime_s)

    # Preserve plot_delta_t.py's automatic range and 0.2-decade bin width.
    log_min = math.floor(min(peak_log - 2.0, float(data_log.min()) - 0.5) * 10.0) / 10.0
    log_max = math.ceil(max(peak_log + 2.0, float(data_log.max()) + 0.5) * 10.0) / 10.0
    log_min -= 0.5
    log_max -= 0.5
    bin_width = 0.3
    bins = np.arange(log_min, log_max + bin_width, bin_width)

    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    counts, _, _ = ax.hist(
        data_log,
        bins=bins,
        color="#D81525",
        edgecolor="black",
        linewidth=0.3,
        label="Data",
    )

    x_axis = np.linspace(log_min, log_max, 800)
    t_real_s = 10.0**x_axis
    curve = (
        math.log(10.0)
        * result.decay_constant_hz
        * t_real_s
        * np.exp(-result.decay_constant_hz * t_real_s)
    )
    if np.max(curve) > 0.0 and np.max(counts) > 0.0:
        curve *= np.max(counts) / np.max(curve)
    ax.plot(x_axis, curve, color="#6A5ACD", linewidth=2.0, label="Fit curve")

    err_low = result.half_life_s - result.half_life_low_s
    err_high = result.half_life_high_s - result.half_life_s
    exponent = math.floor(math.log10(result.half_life_s))
    if exponent == 0:
        half_life_label = (
            rf"$T_{{1/2}}^{{\rm }}="
            rf"({result.half_life_s:.3f}"
            rf"^{{+{err_high:.3f}}}_{{-{err_low:.3f}}})\ \rm s$"
        )
    else:
        scale = 10.0**exponent
        half_life_label = (
            rf"$T_{{1/2}}^{{\rm SF}}="
            rf"({result.half_life_s / scale:.3f}"
            rf"^{{+{err_high / scale:.3f}}}_{{-{err_low / scale:.3f}}})"
            rf"\times10^{{{exponent}}}\ \rm s$"
        )
    ax.text(0.05, 0.91, r"$^{246}$Fm", transform=ax.transAxes, fontsize=22, ha="left")
    ax.text(
        0.05,
        0.80,
        half_life_label,
        transform=ax.transAxes,
        fontsize=15,
        ha="left",
    )
    ax.set_xlim(log_min, log_max)
    ax.set_xlabel(r"$\log_{10}[t\;(\mathrm{s})]$", fontsize=18)
    ax.set_ylabel("Counts", fontsize=18)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.tick_params(top=True, right=True, which="both", labelsize=14)
    ax.legend(loc="upper right", frameon=False, fontsize=11)

    for suffix in (".png", ".pdf"):
        fig.savefig(output_stem.with_suffix(suffix), dpi=300, facecolor="white")
    plt.close(fig)


def plot_energy(energies_kev: np.ndarray, result: FitResults, output_stem: Path) -> None:
    """Plot the energy histogram and the unbinned Gaussian MLE."""
    edges = np.histogram_bin_edges(energies_kev, bins="fd")
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.hist(
        energies_kev,
        bins=edges,
        color="#4C9F70",
        edgecolor="black",
        linewidth=0.5,
        label="Data",
    )

    x = np.linspace(float(edges[0]), float(edges[-1]), 800)
    gaussian = result.energy_gaussian_amplitude * np.exp(
        -0.5 * ((x - result.energy_gaussian_mean_kev) / result.energy_gaussian_sigma_kev) ** 2
    )
    ax.plot(x, gaussian, color="#5B3A9D", linewidth=2.0, label="Gaussian fit")

    ax.text(0.05, 0.91, r"$^{246}$Fm", transform=ax.transAxes, fontsize=22, ha="left")
    ax.text(
        0.05,
        0.80,
        rf"$\mu={result.energy_gaussian_mean_kev:.2f}\pm"
        rf"{result.energy_gaussian_mean_error_kev:.2f}\ \rm keV$"
        "\n"
        rf"$\sigma={result.energy_gaussian_sigma_kev:.2f}\pm"
        rf"{result.energy_gaussian_sigma_error_kev:.2f}\ \rm keV$",
        transform=ax.transAxes,
        fontsize=15,
        ha="left",
        va="top",
    )
    ax.set_xlabel(r"Energy (keV)", fontsize=18)
    ax.set_ylabel("Counts", fontsize=18)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.tick_params(top=True, right=True, which="both", labelsize=14)
    ax.legend(loc="upper right", frameon=False, fontsize=11)

    for suffix in (".png", ".pdf"):
        fig.savefig(output_stem.with_suffix(suffix), dpi=300, facecolor="white")
    plt.close(fig)


def write_results(
    path: Path,
    source: Path,
    input_time_unit: str,
    confidence: float,
    result: FitResults,
) -> None:
    err_low = result.half_life_s - result.half_life_low_s
    err_high = result.half_life_high_s - result.half_life_s
    text = f"""246Fm fit results
=================
Input: {source.resolve()}
Events: {result.event_count}
Input time unit: {input_time_unit}

Decay-time fit (unbinned exponential MLE)
Time unit: s
Confidence level: {confidence * 100:.2f}%
Mean lifetime tau: {result.mean_lifetime_s:.9e} s
Decay constant lambda: {result.decay_constant_hz:.9e} s^-1
Half-life T_1/2: {result.half_life_s:.9e} +{err_high:.9e} -{err_low:.9e} s
Half-life interval: [{result.half_life_low_s:.9e}, {result.half_life_high_s:.9e}] s

Energy fit (Gaussian fit to a Freedman-Diaconis histogram with Poisson errors)
Gaussian mean: {result.energy_gaussian_mean_kev:.6f} +/- {result.energy_gaussian_mean_error_kev:.6f} keV
Gaussian mean: {result.energy_gaussian_mean_kev / 1000.0:.9f} +/- {result.energy_gaussian_mean_error_kev / 1000.0:.9f} MeV
Gaussian sigma: {result.energy_gaussian_sigma_kev:.6f} +/- {result.energy_gaussian_sigma_error_kev:.6f} keV
Arithmetic mean of all energies: {result.energy_arithmetic_mean_kev:.6f} keV
"""
    path.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Fit 246Fm half-life and energy distributions.")
    parser.add_argument("input", nargs="?", type=Path, default=here / "246Fm.dat")
    parser.add_argument(
        "--input-time-unit",
        choices=("auto", "ns", "us", "ms", "s"),
        default="auto",
        help="Input DeltaT unit; default reads it from the column header.",
    )
    parser.add_argument("--confidence", type=float, default=0.6827)
    parser.add_argument("--output-dir", type=Path, default=here)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 0.0 < args.confidence < 1.0:
        raise ValueError("--confidence must be between 0 and 1.")

    _, energies_kev, times_s, input_time_unit = read_root_scan(
        args.input.resolve(), args.input_time_unit
    )
    result = fit_values(energies_kev, times_s, args.confidence)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    configure_style()
    plot_reverse_clock(times_s, result, output_dir / "246Fm_half_life_MLE")
    plot_energy(energies_kev, result, output_dir / "246Fm_energy_gaussian")
    write_results(
        output_dir / "246Fm_fit_results.txt",
        args.input,
        input_time_unit,
        args.confidence,
        result,
    )

    err_low = result.half_life_s - result.half_life_low_s
    err_high = result.half_life_high_s - result.half_life_s
    print(f"Events: {result.event_count}")
    print(f"Input time unit: {input_time_unit}; converted to s")
    print(f"T_1/2 = {result.half_life_s:.9e} +{err_high:.9e} -{err_low:.9e} s")
    print(
        f"Gaussian energy mean = {result.energy_gaussian_mean_kev:.6f} "
        f"+/- {result.energy_gaussian_mean_error_kev:.6f} keV"
    )
    print(
        f"Gaussian sigma = {result.energy_gaussian_sigma_kev:.6f} "
        f"+/- {result.energy_gaussian_sigma_error_kev:.6f} keV"
    )
    print(f"Arithmetic energy mean = {result.energy_arithmetic_mean_kev:.6f} keV")
    print(f"Outputs written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
