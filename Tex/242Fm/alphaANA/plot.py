#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def find_column(df, keywords):
    """
    根据关键词模糊寻找列名。
    keywords: list[str]
    """
    for col in df.columns:
        col_str = str(col)
        if all(k in col_str for k in keywords):
            return col

    print("\nAvailable columns:")
    for c in df.columns:
        print(repr(c))

    raise KeyError(f"Cannot find column with keywords: {keywords}")


def main():
    xlsx_path = Path("/Users/evalie/Documents/242FmSF分析.xlsx")
    #sheet_name = "Alpha_shortLife"
    sheet_name = "Alpha2"
    if not xlsx_path.exists():
        raise FileNotFoundError(f"Cannot find file: {xlsx_path}")

    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

    print("Columns in sheet:")
    for c in df.columns:
        print(repr(c))

    # 自动寻找列名
    e_col = find_column(df, ["DSSD_E[1]", "SSD_E"])
    t_col = find_column(df, ["Delta_Ts[1]", "1e6"])

    print(f"\nUsing energy column: {e_col}")
    print(f"Using time column:   {t_col}")

    # 转成数值
    df[e_col] = pd.to_numeric(df[e_col], errors="coerce")
    df[t_col] = pd.to_numeric(df[t_col], errors="coerce")

    data = df[[e_col, t_col]].dropna()

    energy = data[e_col].to_numpy()
    time_ms = data[t_col].to_numpy()

    positive = time_ms > 0

    print(f"\nTotal valid points: {len(data)}")
    print(f"Energy range: {energy.min():.3f} - {energy.max():.3f} keV")
    print(f"Time range:   {time_ms.min():.6g} - {time_ms.max():.6g} ms")

    # ============================================================
    # 1. DSSD_E[1]+SSD_E 一维能量直方图
    # ============================================================
    plt.figure(figsize=(8, 6))

    plt.hist(
        energy,
        bins=100,
        histtype="step",
        linewidth=1.5
    )

    plt.xlabel(r"$DSSD\_E[1]+SSD\_E$ (keV)")
    plt.ylabel("Counts")
    plt.title("Alpha_shortLife: Energy Spectrum")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("Alpha_shortLife_E_hist.png", dpi=300)
    plt.show()

    # ============================================================
    # 2. Delta_Ts[1]/1e6 一维时间直方图，线性时间轴
    # ============================================================
    plt.figure(figsize=(8, 6))

    plt.hist(
        time_ms,
        bins=100,
        histtype="step",
        linewidth=1.5
    )

    plt.xlabel(r"$\Delta T_s[1]/10^6$ (ms)")
    plt.ylabel("Counts")
    plt.title("Alpha_shortLife: Decay Time Spectrum")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("Alpha_shortLife_T_hist_linear.png", dpi=300)
    plt.show()

    # ============================================================
    # 3. Delta_Ts[1]/1e6 一维时间直方图，log10 时间轴
    #    这个更适合看 us-ms-s 跨很多数量级的分布
    # ============================================================
    plt.figure(figsize=(8, 6))

    log_time = np.log10(time_ms[positive])

    plt.hist(
        log_time,
        bins=100,
        histtype="step",
        linewidth=1.5
    )

    plt.xlabel(r"$\log_{10}(\Delta T_s[1]/10^6\ \mathrm{ms})$")
    plt.ylabel("Counts")
    plt.title("Alpha_shortLife: Decay Time Spectrum in log10 Time")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("Alpha_shortLife_logT_hist.png", dpi=300)
    plt.show()

    # ============================================================
    # 4. 能量-时间关联图，线性时间轴
    # ============================================================
    plt.figure(figsize=(8, 6))

    plt.scatter(
        time_ms,
        energy,
        s=18,
        alpha=0.75
    )

    plt.xlabel(r"$\Delta T_s[1]/10^6$ (ms)")
    plt.ylabel(r"$DSSD\_E[1]+SSD\_E$ (keV)")
    plt.title("Alpha_shortLife: Energy vs Decay Time")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("Alpha_shortLife_E_vs_T_scatter.png", dpi=300)
    plt.show()

    # ============================================================
    # 5. 能量-时间关联图，log10 时间轴
    # ============================================================
    plt.figure(figsize=(8, 6))

    plt.scatter(
        log_time,
        energy[positive],
        s=18,
        alpha=0.75
    )

    plt.xlabel(r"$\log_{10}(\Delta T_s[1]/10^6\ \mathrm{ms})$")
    plt.ylabel(r"$DSSD\_E[1]+SSD\_E$ (keV)")
    plt.title("Alpha_shortLife: Energy vs log10 Decay Time")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("Alpha_shortLife_E_vs_logT_scatter.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()