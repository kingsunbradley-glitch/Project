#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def sample_decay_times(n, half_life_us, rng):
    tau = half_life_us / np.log(2)
    return rng.exponential(scale=tau, size=n)


def main():
    rng = np.random.default_rng(12345)

    # ============================================================
    # 模拟两组 alpha
    # ============================================================
    n1 = 100000
    n2 = 100000

    e1_center = 8600.0      # keV
    e2_center = 8000.0      # keV
    e_sigma = 40.0          # keV

    t12_1_us = 1.0          # 1 us
    t12_2_us = 1000.0       # 1 ms = 1000 us

    e1 = rng.normal(e1_center, e_sigma, n1)
    e2 = rng.normal(e2_center, e_sigma, n2)

    t1_us = sample_decay_times(n1, t12_1_us, rng)
    t2_us = sample_decay_times(n2, t12_2_us, rng)

    energy = np.concatenate([e1, e2])
    time_us = np.concatenate([t1_us, t2_us])

    mask = time_us > 0
    energy = energy[mask]
    time_us = time_us[mask]

    logt = np.log10(time_us)

    # ============================================================
    # 二维直方图
    # ============================================================
    x_range = [-3, 6]
    y_range = [7800, 8800]
    x_bins = 160
    y_bins = 160

    H, xedges, yedges = np.histogram2d(
        logt,
        energy,
        bins=[x_bins, y_bins],
        range=[x_range, y_range]
    )

    x_centers = 0.5 * (xedges[:-1] + xedges[1:])
    y_centers = 0.5 * (yedges[:-1] + yedges[1:])

    x_projection = H.sum(axis=1)
    y_projection = H.sum(axis=0)

    # ============================================================
    # 同一张图：上方 X projection，中间 2D，右侧 Y projection
    # ============================================================
    fig = plt.figure(figsize=(9, 8))

    gs = GridSpec(
        2, 2,
        width_ratios=[4, 1.2],
        height_ratios=[1.2, 4],
        hspace=0.05,
        wspace=0.05
    )

    ax_xproj = fig.add_subplot(gs[0, 0])
    ax_2d = fig.add_subplot(gs[1, 0], sharex=ax_xproj)
    ax_yproj = fig.add_subplot(gs[1, 1], sharey=ax_2d)

    # ============================================================
    # 主图：二维热力图
    # ============================================================
    mesh = ax_2d.pcolormesh(
        xedges,
        yedges,
        H.T,
        shading="auto"
    )

    ax_2d.set_xlabel(r"$\log_{10}(T/\mu s)$")
    ax_2d.set_ylabel(r"$E_\alpha$ (keV)")

    ax_2d.axvline(np.log10(t12_1_us), linestyle="--", linewidth=1.2)
    ax_2d.axvline(np.log10(t12_2_us), linestyle="--", linewidth=1.2)

    ax_2d.axhline(e1_center, linestyle="--", linewidth=1.2)
    ax_2d.axhline(e2_center, linestyle="--", linewidth=1.2)

    # ============================================================
    # X 投影：logT 分布
    # ============================================================
    ax_xproj.step(
        x_centers,
        x_projection,
        where="mid",
        linewidth=1.3
    )

    ax_xproj.axvline(np.log10(t12_1_us), linestyle="--", linewidth=1.2)
    ax_xproj.axvline(np.log10(t12_2_us), linestyle="--", linewidth=1.2)

    ax_xproj.set_ylabel("Counts")
    ax_xproj.tick_params(labelbottom=False)
    ax_xproj.grid(True, alpha=0.3)

    # ============================================================
    # Y 投影：能量谱
    # ============================================================
    ax_yproj.step(
        y_projection,
        y_centers,
        where="mid",
        linewidth=1.3
    )

    ax_yproj.axhline(e1_center, linestyle="--", linewidth=1.2)
    ax_yproj.axhline(e2_center, linestyle="--", linewidth=1.2)

    ax_yproj.set_xlabel("Counts")
    ax_yproj.tick_params(labelleft=False)
    ax_yproj.grid(True, alpha=0.3)

    # ============================================================
    # colorbar
    # ============================================================
    cbar = fig.colorbar(
        mesh,
        ax=[ax_2d, ax_yproj],
        fraction=0.046,
        pad=0.04
    )
    cbar.set_label("Counts/bin")

    fig.suptitle(
        r"$E_\alpha$ vs $\log_{10}(T)$ with X/Y Projections",
        y=0.98
    )

    plt.savefig("E_vs_logT_with_XY_projection.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()