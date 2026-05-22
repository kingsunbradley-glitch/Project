#!/Users/evalie/miniconda3/envs/py_env/bin/python
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ.get("TMPDIR", "/private/tmp")) / "matplotlib-242fm"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import chi2


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def split_fields(line: str) -> list[str]:
    if "," in line:
        return [field.strip() for field in next(csv.reader([line]))]
    return re.split(r"\s+", line.strip())


def read_delta_t(path: Path, column: str = "Delta_t1") -> np.ndarray:
    lines = [strip_comment(line) for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        raise ValueError(f"{path} is empty.")

    header = split_fields(lines[0])
    if column not in header:
        raise ValueError(f"Column {column!r} was not found in {path}.")

    column_index = header.index(column)
    values: list[float] = []
    for line_number, line in enumerate(lines[1:], start=2):
        fields = split_fields(line)
        if len(fields) != len(header):
            raise ValueError(f"Line {line_number}: expected {len(header)} columns, got {len(fields)}.")
        try:
            value = float(fields[column_index])
        except ValueError as exc:
            raise ValueError(f"Line {line_number}: {column} must be numeric, got {fields[column_index]!r}.") from exc
        if value <= 0:
            raise ValueError(f"Line {line_number}: {column} must be positive, got {value}.")
        values.append(value)

    if not values:
        raise ValueError(f"No {column} values found in {path}.")
    return np.array(values, dtype=float)


def theoretical_curve_log10(t: np.ndarray, lambda_est: float) -> np.ndarray:
    # For x=log10(t), dP/dx = ln(10) * lambda * t * exp(-lambda * t).
    return math.log(10.0) * lambda_est * t * np.exp(-lambda_est * t)


def half_life_interval(times_us: np.ndarray, confidence: float) -> tuple[float, float, float, float, float]:
    n = times_us.size
    total_time = float(np.sum(times_us))
    mean_life = total_time / n
    half_life = math.log(2.0) * mean_life

    alpha = 1.0 - confidence
    dof = 2 * n
    mean_life_low = 2.0 * total_time / chi2.ppf(1.0 - alpha / 2.0, dof)
    mean_life_high = 2.0 * total_time / chi2.ppf(alpha / 2.0, dof)
    half_life_low = math.log(2.0) * mean_life_low
    half_life_high = math.log(2.0) * mean_life_high

    err_low = half_life - half_life_low
    err_high = half_life_high - half_life
    return half_life, err_low, err_high, half_life_low, half_life_high


def plot_delta_t(
    times_us: np.ndarray,
    output_path: Path,
    bins_count: int,
    log_min_arg: float | None,
    log_max_arg: float | None,
    half_life: float,
    err_low: float,
    err_high: float,
    interval_low: float,
    interval_high: float,
) -> tuple[float, float, float, float, float]:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["STIXGeneral", "Times New Roman", "serif"]
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

    mean_life = float(np.mean(times_us))
    lambda_est = 1.0 / mean_life
    peak_log = math.log10(mean_life)

    data_log = np.log10(times_us)
    log_min = math.floor(min(peak_log - 2.0, float(data_log.min()) - 0.5) * 10.0) / 10.0
    log_max = math.ceil(max(peak_log + 2.0, float(data_log.max()) + 0.5) * 10.0) / 10.0
    if log_min_arg is not None:
        log_min = log_min_arg
    if log_max_arg is not None:
        log_max = log_max_arg
    if log_min_arg is None and log_max_arg is None:
        log_min -= 0.5
        log_max -= 0.5
    if log_max <= log_min:
        raise ValueError("log-max must be larger than log-min.")

    bins = np.linspace(log_min, log_max, bins_count + 1)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(
        data_log,
        bins=bins,
        color="#D81525",
        edgecolor="black",
        linewidth=0.3,
    )

    x_axis = np.linspace(log_min, log_max, 500)
    t_real = 10**x_axis
    y_curve = theoretical_curve_log10(t_real, lambda_est)

    counts, _ = np.histogram(data_log, bins=bins)
    max_count = np.max(counts) if counts.size > 0 else 1
    if max_count <= 0:
        max_count = 1
    if y_curve.max() > 0:
        ax.plot(x_axis, y_curve * (max_count / y_curve.max()), "-", linewidth=2, color="#6A5ACD")

    ax.set_xlim(log_min, log_max)
    ax.set_xlabel(r"$\log_{10}[t\;(\mu\mathrm{s})]$", fontsize=18)
    ax.set_ylabel("Counts", fontsize=18)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=14)

    for spine in ax.spines.values():
        spine.set_visible(True)

    ax.text(0.05, 0.9, r"$^{242}$Fm", transform=ax.transAxes, fontsize=22, ha="left")
    ax.text(
        -1.75,
        2.0,
        "\n".join(
            [
                rf"$T_{{1/2}}^{{SF}} = {half_life:.2f}^{{+{err_high:.2f}}}_{{-{err_low:.2f}}}\ \mu\mathrm{{s}}$",
            ]
        ),
        fontsize=16,
        ha="left",
        va="top",
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return lambda_est, peak_log, mean_life, log_min, log_max


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Plot Delta_t1 data from input.txt in microseconds.")
    parser.add_argument("input", nargs="?", type=Path, default=here / "input.txt")
    parser.add_argument("-o", "--output", type=Path, default=here / "delta_t_distribution.pdf")
    parser.add_argument("--bins", type=int, default=40, help="Number of equal-width bins in log10(t).")
    parser.add_argument("--log-min", type=float, default=None, help="Optional lower x-axis limit in log10(us).")
    parser.add_argument("--log-max", type=float, default=None, help="Optional upper x-axis limit in log10(us).")
    parser.add_argument("--confidence", type=float, default=0.6827, help="Two-sided confidence level.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        times_us = read_delta_t(args.input.resolve())
        if not 0.0 < args.confidence < 1.0:
            raise ValueError("--confidence must be between 0 and 1.")
        if args.bins < 1:
            raise ValueError("--bins must be a positive integer.")

        half_life, err_low, err_high, interval_low, interval_high = half_life_interval(times_us, args.confidence)
        lambda_est, peak_log, mean_life, log_min, log_max = plot_delta_t(
            times_us,
            args.output.resolve(),
            args.bins,
            args.log_min,
            args.log_max,
            half_life,
            err_low,
            err_high,
            interval_low,
            interval_high,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Read {times_us.size} Delta_t1 value(s) from {args.input.resolve()}.")
    print("Delta_t1 values are treated as microseconds (us).")
    print(f"Values (us): {', '.join(f'{value:g}' for value in times_us)}")
    print(f"Mean lifetime tau = {mean_life:.3f} us; lambda = {lambda_est:.6f} us^-1.")
    print(f"Log-time curve peaks at log10(t/us) = {peak_log:.3f}.")
    print(f"Plot x-range: log10(t/us) = [{log_min:.3f}, {log_max:.3f}].")
    print(f"T_1/2 = {half_life:.3f} +{err_high:.3f} -{err_low:.3f} us ({args.confidence * 100:.2f}% C.L.)")
    print(f"Confidence interval: [{interval_low:.3f}, {interval_high:.3f}] us")
    print(f"Wrote {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
