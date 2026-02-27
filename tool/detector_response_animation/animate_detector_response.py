#!/usr/bin/env python3
"""
Animate detector scan:
1) Fixed discrete levels on top.
2) Detector response function scans from low to high energy.
3) Resolution worsens with energy (larger FWHM at high energy).
4) Output split into low-E discrete, high-E discrete, and high-E dense (unresolved) parts.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def gaussian(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    norm = sigma * np.sqrt(2.0 * np.pi)
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / norm


def fwhm_vs_energy(energy: np.ndarray, e_min: float, e_max: float) -> np.ndarray:
    # Increase FWHM with energy: e.g. 70 keV -> 130 keV across the scan range.
    fwhm_low = 0.070
    fwhm_high = 0.130
    frac = (energy - e_min) / (e_max - e_min)
    frac = np.clip(frac, 0.0, 1.0)
    return fwhm_low + (fwhm_high - fwhm_low) * frac


def broaden_group(
    x: np.ndarray,
    level_e: np.ndarray,
    level_w: np.ndarray,
    e_min: float,
    e_max: float,
) -> np.ndarray:
    out = np.zeros_like(x)
    for e0, w in zip(level_e, level_w):
        sigma = fwhm_vs_energy(np.array([e0]), e_min, e_max)[0] / 2.355
        out += w * gaussian(x, e0, sigma)
    return out


def masked_prefix(y: np.ndarray, end_idx: int) -> np.ndarray:
    out = np.full_like(y, np.nan)
    out[: end_idx + 1] = y[: end_idx + 1]
    return out


def main() -> None:
    # Three groups:
    # 1) low-energy discrete lines
    low_e = np.array([8.52, 8.66, 8.82], dtype=float)
    low_w = np.array([1.00, 0.78, 0.62], dtype=float)

    # 2) high-energy discrete lines
    high_disc_e = np.array([9.42, 9.68, 9.96], dtype=float)
    high_disc_w = np.array([0.76, 0.64, 0.58], dtype=float)

    # 3) high-energy dense lines that become unresolved
    high_dense_e = np.array([10.23, 10.26, 10.29, 10.32, 10.35, 10.38, 10.41, 10.44], dtype=float)
    high_dense_w = np.array([0.30, 0.36, 0.41, 0.44, 0.42, 0.39, 0.34, 0.28], dtype=float)

    all_e = np.concatenate([low_e, high_disc_e, high_dense_e])
    all_w = np.concatenate([low_w, high_disc_w, high_dense_w])

    e_min = float(all_e.min() - 0.45)
    e_max = float(all_e.max() + 0.45)
    e_axis = np.linspace(e_min, e_max, 1400)

    low_resp = broaden_group(e_axis, low_e, low_w, e_min, e_max)
    high_disc_resp = broaden_group(e_axis, high_disc_e, high_disc_w, e_min, e_max)
    high_dense_resp = broaden_group(e_axis, high_dense_e, high_dense_w, e_min, e_max)
    total_resp = low_resp + high_disc_resp + high_dense_resp

    fig, (ax_levels, ax_kernel, ax_out) = plt.subplots(
        3, 1, figsize=(11, 9), sharex=True, gridspec_kw={"height_ratios": [1, 1, 2]}
    )

    ax_levels.vlines(low_e, 0.0, low_w, color="tab:blue", lw=2.5, label="Low-E discrete levels")
    ax_levels.vlines(
        high_disc_e, 0.0, high_disc_w, color="tab:orange", lw=2.5, label="High-E discrete levels"
    )
    ax_levels.vlines(
        high_dense_e, 0.0, high_dense_w, color="tab:green", lw=2.0, label="High-E dense levels"
    )
    ax_levels.set_ylabel("Relative intensity")
    ax_levels.set_title("Fixed Input Levels")
    ax_levels.set_ylim(0.0, 1.18 * all_w.max())
    ax_levels.grid(alpha=0.25)
    ax_levels.legend(loc="upper right", frameon=False)

    kernel_line, = ax_kernel.plot(e_axis, np.zeros_like(e_axis), color="tab:red", lw=2.2)
    scan_line_kernel = ax_kernel.axvline(e_axis[0], color="k", ls="--", lw=1.0, alpha=0.75)
    ax_kernel.set_ylabel("R(E | Escan)")
    ax_kernel.set_title("Detector Response Function During Scan")
    ax_kernel.set_ylim(0.0, 1.05)
    ax_kernel.grid(alpha=0.25)

    low_line, = ax_out.plot(
        e_axis, np.full_like(e_axis, np.nan), color="tab:blue", lw=2.0, label="_nolegend_"
    )
    high_disc_line, = ax_out.plot(
        e_axis,
        np.full_like(e_axis, np.nan),
        color="tab:orange",
        lw=2.0,
        label="_nolegend_",
    )
    high_dense_line, = ax_out.plot(
        e_axis,
        np.full_like(e_axis, np.nan),
        color="tab:green",
        lw=2.0,
        label="_nolegend_",
    )
    total_line, = ax_out.plot(
        e_axis, np.full_like(e_axis, np.nan), color="k", lw=2.4, label="Total detected signal"
    )

    scan_line_out = ax_out.axvline(e_axis[0], color="k", ls="--", lw=1.0, alpha=0.75)

    ax_out.set_xlabel("Energy (MeV)")
    ax_out.set_ylabel("Counts (a.u.)")
    ax_out.set_title("Output Built by Low->High Detector Scan")
    ax_out.set_ylim(0.0, 1.10 * total_resp.max())
    ax_out.grid(alpha=0.25)
    ax_out.legend(loc="upper right", frameon=False)

    info_text = ax_out.text(0.02, 0.96, "", transform=ax_out.transAxes, va="top")

    total_frames = len(e_axis)
    kernel_scale = 1.0

    def update(frame: int):
        e_scan = e_axis[frame]
        fwhm = fwhm_vs_energy(np.array([e_scan]), e_min, e_max)[0]
        sigma = fwhm / 2.355

        kernel = gaussian(e_axis, e_scan, sigma)
        kernel /= kernel.max()
        kernel_line.set_ydata(kernel * kernel_scale)
        scan_line_kernel.set_xdata([e_scan, e_scan])

        low_line.set_ydata(masked_prefix(low_resp, frame))
        high_disc_line.set_ydata(masked_prefix(high_disc_resp, frame))
        high_dense_line.set_ydata(masked_prefix(high_dense_resp, frame))
        total_line.set_ydata(masked_prefix(total_resp, frame))
        scan_line_out.set_xdata([e_scan, e_scan])

        info_text.set_text(
            f"Escan = {e_scan:.3f} MeV\n"
            f"FWHM(Escan) = {fwhm * 1000:.1f} keV\n"
            f"sigma(Escan) = {sigma * 1000:.1f} keV"
        )

        return [
            kernel_line,
            scan_line_kernel,
            low_line,
            high_disc_line,
            high_dense_line,
            total_line,
            scan_line_out,
            info_text,
        ]

    anim = FuncAnimation(fig, update, frames=total_frames, interval=1, blit=False, repeat=True)
    _ = anim
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
