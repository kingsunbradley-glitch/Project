#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Qα from multiple mass models (FRDM, HFB, WS4 + New Models)
and Experimental data (EXP).
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
                # Handle different column formats if necessary
                if len(nums) >= 4:
                    # Format: Z N A Q
                    Z = int(round(nums[0]))
                    N = int(round(nums[1]))
                    A = int(round(nums[2]))
                    Q = float(nums[3])
                    rows.append((Z, N, A, Q))
                elif len(nums) == 3:
                    # Format: A Z Q (derive N)
                    A = int(round(nums[0]))
                    Z = int(round(nums[1]))
                    Q = float(nums[2])
                    N = A - Z
                    rows.append((Z, N, A, Q))
    except FileNotFoundError:
        print(f"[WARN] File not found: {path} (Skipping model {model_name})")
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
                if nums is None or len(nums) < 4:
                    continue
                # EXP Format usually: A Z N Q(keV) Err(keV) ...
                if len(nums) >= 5:
                    A = int(round(nums[0]))
                    Z = int(round(nums[1]))
                    N = int(round(nums[2]))
                    Q_kev = float(nums[3])
                    Err_kev = float(nums[4])
                    Q_mev = Q_kev / 1000.0
                    Err_mev = Err_kev / 1000.0
                    rows.append((Z, N, A, Q_mev, Err_mev))
    except FileNotFoundError:
        print(f"[WARN] EXP file not found: {path}")
        return pd.DataFrame()
    
    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha", "Error"])
    df["model"] = "EXP"
    return df

def plot_element(ax, df_models: pd.DataFrame, df_exp: pd.DataFrame, Z: int, sym: str, nmin: int, nmax: int):
    # Define styles for each model to distinguish them
    styles = {
        "FRDM(2012)": {"c": "blue", "m": "o"},
        "HFB-32":     {"c": "green", "m": "v"},
        "WS4+RBF":    {"c": "red", "m": "^"},
        "UNEDF1":     {"c": "purple", "m": "s"},
        "SKMS":       {"c": "cyan", "m": "p"},
        "SLY4":       {"c": "magenta", "m": "*"},
        "SV-MIN":     {"c": "orange", "m": "D"},
    }
    
    # Plot Models
    if not df_models.empty:
        # Filter for current Z and N range
        d = df_models[(df_models["Z"] == Z) & (df_models["N"] >= nmin) & (df_models["N"] <= nmax)]
        
        # Plot in a specific order if desired
        model_order = ["FRDM(2012)", "HFB-32", "WS4+RBF", "UNEDF1", "SKMS", "SLY4", "SV-MIN"]
        
        for model_name in model_order:
            if model_name not in d["model"].unique():
                continue
            
            dm = d[d["model"] == model_name].sort_values("N")
            if dm.empty:
                continue
            
            # Use defined style or default
            st = styles.get(model_name, {"c": "gray", "m": "."})
            
            ax.plot(dm["N"], dm["Qalpha"], 
                    marker=st["m"], color=st["c"],
                    markersize=4, linewidth=1.5,
                    label=model_name, alpha=0.8)

    # Plot Experimental Data
    if not df_exp.empty:
        de = df_exp[(df_exp["Z"] == Z) & (df_exp["N"] >= nmin) & (df_exp["N"] <= nmax)].sort_values("N")
        if not de.empty:
            ax.errorbar(de["N"], de["Qalpha"], yerr=de["Error"], 
                        fmt='s', color='black', label='EXP', 
                        capsize=3, markersize=5, zorder=10, 
                        markerfacecolor='none', markeredgewidth=1.5)

    # Formatting
    ax.axvline(152, linestyle=":", color="gray", linewidth=1, alpha=0.5)
    ax.axvline(162, linestyle=":", color="gray", linewidth=1, alpha=0.5)
    ax.set_title(f"{sym} (Z={Z})")
    ax.set_xlabel("Neutron number N")
    ax.set_ylabel(r"$Q_\alpha$ (MeV)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, linestyle=":", linewidth=0.8)
    
    # Legend
    ax.legend(fontsize='large', loc='best')

def main():
    parser = argparse.ArgumentParser(description="Plot Qalpha comparison including new models.")
    parser.add_argument("--nmin", type=int, default=145, help="Min Neutron number")
    parser.add_argument("--nmax", type=int, default=180, help="Max Neutron number")
    args = parser.parse_args()

    # 1. Define all model files
    model_files = {
        #"FRDM(2012)": "Q@a_FRDM.dat",
        #"HFB-32":     "Q@a_HFB.dat",
        #"WS4+RBF":    "Q@a_WS4+RBF.dat",
        "UNEDF1":     "Q@a_UNEDF1.dat",
        #"SKMS":       "Q@a_SKMS.dat", #    all over estimate
        "SLY4":       "Q@a_SLY4.dat",
        "SV-MIN":     "Q@a_SV-MIN.dat"
    }

    # 2. Read all models
    dfs = []
    print("Loading model data...")
    for name, path in model_files.items():
        df = read_qalpha_table(path, name)
        if not df.empty:
            dfs.append(df)
            print(f"  - Loaded {name}: {len(df)} records")
    
    if dfs:
        df_models = pd.concat(dfs, ignore_index=True)
    else:
        df_models = pd.DataFrame(columns=["Z", "N", "Qalpha", "model"])

    # 3. Read EXP data
    df_exp = read_exp_table("Q@a_EXP.dat")
    print(f"  - Loaded EXP: {len(df_exp)} records")

    # 4. Plotting
    elements = [(104, "Rf"), (106, "Sg"), (108, "Hs"), (110, "Ds")]
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    print("Generating plots...")
    for i, (Z, sym) in enumerate(elements):
        if i < len(axes):
            plot_element(axes[i], df_models, df_exp, Z, sym, args.nmin, args.nmax)

    plt.tight_layout()
    output_filename = "Qalpha_comparison_all_models.png"
    plt.savefig(output_filename, dpi=300)
    #print(f"Plot saved to {output_filename}")
    plt.show()
    # plt.show() # Uncomment if running locally with display

if __name__ == "__main__":
    main()