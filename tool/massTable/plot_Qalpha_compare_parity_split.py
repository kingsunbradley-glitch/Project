#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Qα (alpha-decay Q value) for Rf/Sg/Hs/Ds isotopic chains.
Layout: 4 Columns (Elements) x 2 Rows (Parity).
Top Row: Even-N nuclei.
Bottom Row: Odd-N nuclei.

Usage:
  python plot_Qalpha_compare_parity_split.py --nmin 145 --nmax 180
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

def plot_single_ax(ax, df_models, df_exp, Z, parity_label):
    """
    Helper to plot data on a single axis for a specific Z and parity.
    parity_label: 'Even' or 'Odd'
    """
    # Filter by parity
    is_even = (parity_label == 'Even')
    
    # Filter Models
    if not df_models.empty:
        # Filter by Z
        dm = df_models[df_models["Z"] == Z]
        # Filter by Parity
        dm = dm[dm["N"] % 2 == (0 if is_even else 1)]
        
        for model_name in sorted(dm["model"].dropna().unique()):
            d = dm[dm["model"] == model_name].sort_values("N")
            if not d.empty:
                ax.plot(d["N"], d["Qalpha"], marker="o", markersize=4, linewidth=1.5,
                        label=model_name, alpha=0.8)

    # Filter EXP
    if not df_exp.empty:
        de = df_exp[df_exp["Z"] == Z]
        de = de[de["N"] % 2 == (0 if is_even else 1)].sort_values("N")
        if not de.empty:
            ax.errorbar(de["N"], de["Qalpha"], yerr=de["Error"], 
                        fmt='s', color='black', label='EXP', 
                        capsize=3, markersize=5, zorder=10)

def main():
    ap = argparse.ArgumentParser(description="Plot Qalpha comparison 4x2 grid (Even Top / Odd Bottom).")
    ap.add_argument("--frdm", default="Q@a_FRDM.dat")
    ap.add_argument("--hfb",  default="Q@a_HFB.dat")
    ap.add_argument("--ws4",  default="Q@a_WS4+RBF.dat")
    ap.add_argument("--exp",  default="Q@a_EXP.dat")
    ap.add_argument("--nmin", type=int, default=145)
    ap.add_argument("--nmax", type=int, default=180)
    args = ap.parse_args()

    # Read Data
    df_frdm = read_qalpha_table(args.frdm, "FRDM(2012)")
    df_hfb  = read_qalpha_table(args.hfb,  "HFB-24")
    df_ws4  = read_qalpha_table(args.ws4,  "WS4+RBF")
    df_exp  = read_exp_table(args.exp)
    
    # Filter N range globally for convenience
    dfs = [df for df in (df_frdm, df_hfb, df_ws4) if not df.empty]
    df_models = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=["Z", "N", "Qalpha", "model"])
    
    if not df_models.empty:
        df_models = df_models[(df_models["N"] >= args.nmin) & (df_models["N"] <= args.nmax)]
    if not df_exp.empty:
        df_exp = df_exp[(df_exp["N"] >= args.nmin) & (df_exp["N"] <= args.nmax)]

    elements = [(104, "Rf"), (106, "Sg"), (108, "Hs"), (110, "Ds")]
    
    # 2 Rows (Parity), 4 Columns (Elements)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    for i, (Z, sym) in enumerate(elements):
        # Top Row: Even
        ax_even = axes[0, i]
        plot_single_ax(ax_even, df_models, df_exp, Z, 'Even')
        ax_even.set_title(f"{sym} (Z={Z}) - Even N")
        ax_even.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax_even.grid(True, linestyle=":", linewidth=0.8)
        
        # Bottom Row: Odd
        ax_odd = axes[1, i]
        plot_single_ax(ax_odd, df_models, df_exp, Z, 'Odd')
        ax_odd.set_title(f"{sym} (Z={Z}) - Odd N")
        ax_odd.set_xlabel("Neutron number N")
        ax_odd.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax_odd.grid(True, linestyle=":", linewidth=0.8)

        # Labels and Legends only on the first column to save space
        if i == 0:
            ax_even.set_ylabel(r"$Q_\alpha$ (MeV)")
            ax_odd.set_ylabel(r"$Q_\alpha$ (MeV)")
            ax_even.legend()
            ax_odd.legend()

    plt.tight_layout()
    print("Plot generated.")
    plt.show()

if __name__ == "__main__":
    main()