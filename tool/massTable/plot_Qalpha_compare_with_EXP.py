#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Qα (alpha-decay Q value) from three mass models (FRDM, HFB32, WS4+RBF)
and Experimental data (EXP) for Rf/Sg/Hs/Ds isotopic chains.
X-axis is forced to show only integers.

Usage:
  python plot_Qalpha_compare_integer_axis.py --nmin 145 --nmax 180
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def _try_parse_numeric(parts):
    try:
        return [float(x) for x in parts]
    except Exception:
        return None

def read_qalpha_table(path: str, model_name: str) -> pd.DataFrame:
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#") or s.startswith("[") or not s[0].isdigit():
                    continue
                parts = s.split()
                nums = _try_parse_numeric(parts)
                if nums is None:
                    continue
                if len(nums) >= 4:
                    Z = int(round(nums[0]))
                    N = int(round(nums[1]))
                    A = int(round(nums[2]))
                    Q = float(nums[3])
                    rows.append((Z, N, A, Q))
                elif len(nums) == 3:
                    A = int(round(nums[0]))
                    Z = int(round(nums[1]))
                    Q = float(nums[2])
                    N = A - Z
                    rows.append((Z, N, A, Q))
    except FileNotFoundError:
        print(f"[WARN] File not found: {path}")
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha"])
    df["model"] = model_name
    return df

def read_exp_table(path: str) -> pd.DataFrame:
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#") or s.startswith("[") or not s[0].isdigit():
                    continue
                parts = s.split()
                nums = _try_parse_numeric(parts)
                if nums is None or len(nums) < 5:
                    continue
                A = int(round(nums[0]))
                Z = int(round(nums[1]))
                N = int(round(nums[2]))
                Q_kev = float(nums[3])
                Err_kev = float(nums[4])
                Q_mev = Q_kev / 1000.0
                Err_mev = Err_kev / 1000.0
                rows.append((Z, N, A, Q_mev, Err_mev))
    except FileNotFoundError:
        print(f"[WARN] File not found: {path}")
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha", "Error"])
    df["model"] = "EXP"
    return df

def plot_element(ax, df_models: pd.DataFrame, df_exp: pd.DataFrame, Z: int, sym: str, nmin: int, nmax: int):
    # Plot Models
    if not df_models.empty:
        d = df_models[(df_models["Z"] == Z) & (df_models["N"] >= nmin) & (df_models["N"] <= nmax)]
        for model_name in sorted(d["model"].dropna().unique()):
            dm = d[d["model"] == model_name].sort_values("N")
            if dm.empty:
                continue
            ax.plot(dm["N"], dm["Qalpha"], marker="o", markersize=4, linewidth=1.5,
                    label=model_name, alpha=0.8)

    # Plot EXP
    if not df_exp.empty:
        de = df_exp[(df_exp["Z"] == Z) & (df_exp["N"] >= nmin) & (df_exp["N"] <= nmax)].sort_values("N")
        if not de.empty:
            ax.errorbar(de["N"], de["Qalpha"], yerr=de["Error"], 
                        fmt='s', color='black', label='EXP', 
                        capsize=3, markersize=5, zorder=10)

    # Reference lines
    ax.axvline(152, linestyle=":", color="gray", linewidth=1, alpha=0.5)
    ax.axvline(162, linestyle=":", color="gray", linewidth=1, alpha=0.5)

    ax.set_title(f"{sym} (Z={Z})")
    ax.set_xlabel("Neutron number N")
    ax.set_ylabel(r"$Q_\alpha$ (MeV)")
    
    # --- Force Integer Ticks on X-Axis ---
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    ax.grid(True, linestyle=":", linewidth=0.8)
    ax.legend()

def main():
    ap = argparse.ArgumentParser(description="Plot Qalpha comparison with integer X-axis.")
    ap.add_argument("--frdm", default="Q@a_FRDM.dat")
    ap.add_argument("--hfb",  default="Q@a_HFB.dat")
    ap.add_argument("--ws4",  default="Q@a_WS4+RBF.dat")
    ap.add_argument("--exp",  default="Q@a_EXP.dat")
    ap.add_argument("--nmin", type=int, default=145)
    ap.add_argument("--nmax", type=int, default=180)
    args = ap.parse_args()

    df_frdm = read_qalpha_table(args.frdm, "FRDM(2012)")
    df_hfb  = read_qalpha_table(args.hfb,  "HFB-32")
    df_ws4  = read_qalpha_table(args.ws4,  "WS4+RBF")
    df_exp  = read_exp_table(args.exp)

    dfs = [df for df in (df_frdm, df_hfb, df_ws4) if not df.empty]
    df_models = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=["Z", "N", "Qalpha", "model"])

    elements = [(104, "Rf"), (106, "Sg"), (108, "Hs"), (110, "Ds")]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, (Z, sym) in enumerate(elements):
        if i < len(axes):
            plot_element(axes[i], df_models, df_exp, Z, sym, args.nmin, args.nmax)

    plt.tight_layout()
    print("Plot generated.")
    plt.show()

if __name__ == "__main__":
    main()