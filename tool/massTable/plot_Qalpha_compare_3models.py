#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  python plot_Qalpha_compare_3models.py --nmin 150 --nmax 170
"""
Compare Qα (alpha-decay Q value) from three mass models (FRDM, HFB-24, WS4+RBF)
for Rf/Sg/Hs/Ds isotopic chains.

Supported input formats (whitespace-separated, header lines allowed):
  Format A (FRDM/HFB style):
      Z  N  A  Qalpha
  Format B (WS4+RBF file from imqmd.com/mass):
      A  Z  Qalpha

Lines starting with '#' are ignored, and non-numeric lines are skipped automatically.

# Usage:
#  python plot_Qalpha_compare_3models.py --nmin 150 --nmax 170

Optionally specify filenames:
  python plot_Qalpha_compare_3models.py --frdm "Q@a_FRDM.dat" --hfb "Q@a_HFB.dat" --ws4 "Q@a_WS4+RBF.dat" --nmin 150 --nmax 170

The script will call plt.show() at the end.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def _try_parse_numeric(parts):
    # Return list[float] if possible; else None.
    try:
        return [float(x) for x in parts]
    except Exception:
        return None


def read_qalpha_table(path: str, model_name: str) -> pd.DataFrame:
    """
    Read Qα table and return DataFrame with columns: Z, N, A, Qalpha, model
    Supports:
      - >=4 columns: Z N A Qalpha  (extra columns ignored)
      - 3 columns:   A Z Qalpha    (WS4+RBF)
    """
    rows = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            nums = _try_parse_numeric(parts)
            if nums is None:
                continue

            # Format A: Z N A Q
            if len(nums) >= 4:
                Z = int(round(nums[0]))
                N = int(round(nums[1]))
                A = int(round(nums[2]))
                Q = float(nums[3])
                if A <= 0 or Z <= 0:
                    continue
                if A != Z + N:
                    N = A - Z
                if N < 0:
                    continue
                rows.append((Z, N, A, Q))
                continue

            # Format B: A Z Q
            if len(nums) == 3:
                A = int(round(nums[0]))
                Z = int(round(nums[1]))
                Q = float(nums[2])
                if A <= 0 or Z <= 0 or A < Z:
                    continue
                N = A - Z
                rows.append((Z, N, A, Q))
                continue

    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha"])
    df["model"] = model_name
    return df


def plot_element(df_all: pd.DataFrame, Z: int, sym: str, nmin: int, nmax: int):
    d = df_all[(df_all["Z"] == Z) & (df_all["N"] >= nmin) & (df_all["N"] <= nmax)].copy()
    if d.empty:
        print(f"[WARN] No data for {sym} (Z={Z}) in N=[{nmin},{nmax}] across all models.")
        return

    fig, ax = plt.subplots(figsize=(8.5, 5.2))

    for model_name in sorted(d["model"].dropna().unique()):
        dm = d[d["model"] == model_name].sort_values("N")
        if dm.empty:
            continue
        ax.plot(dm["N"], dm["Qalpha"], marker="o", markersize=3, linewidth=1,
                label=f"{model_name} (n={len(dm)})")

    # Reference lines
    ax.axvline(152, linestyle=":", linewidth=1)
    ax.axvline(162, linestyle=":", linewidth=1)

    ax.set_title(f"{sym} (Z={Z})  Qα comparison  N=[{nmin},{nmax}]")
    ax.set_xlabel("Neutron number N")
    ax.set_ylabel("Qα (MeV)")
    ax.grid(True, linestyle=":", linewidth=0.8)
    ax.legend()


def main():
    ap = argparse.ArgumentParser(description="Plot Qα comparison for Rf/Sg/Hs/Ds from FRDM, HFB-24, WS4+RBF.")
    ap.add_argument("--frdm", default="Q@a_FRDM.dat", help="FRDM Qα table")
    ap.add_argument("--hfb",  default="Q@a_HFB.dat",  help="HFB-24 Qα table")
    ap.add_argument("--ws4",  default="Q@a_WS4+RBF.dat", help="WS4+RBF Qα table (often A Z Qα format)")
    ap.add_argument("--nmin", type=int, required=True, help="start neutron number (inclusive)")
    ap.add_argument("--nmax", type=int, required=True, help="end neutron number (inclusive)")
    args = ap.parse_args()

    if args.nmin > args.nmax:
        raise SystemExit("ERROR: nmin must be <= nmax")

    df_frdm = read_qalpha_table(args.frdm, "FRDM(2012)")
    df_hfb  = read_qalpha_table(args.hfb,  "HFB-24")
    df_ws4  = read_qalpha_table(args.ws4,  "WS4+RBF")

    print(f"[INFO] Read rows: FRDM={len(df_frdm)}  HFB={len(df_hfb)}  WS4={len(df_ws4)}")
    if len(df_ws4) == 0:
        print("[HINT] WS4 file not parsed. Check format/path. Expected either 'A Z Qalpha' or 'Z N A Qalpha'.")

    dfs = [df for df in (df_frdm, df_hfb, df_ws4) if not df.empty]
    if not dfs:
        raise SystemExit("ERROR: No data read from any input file.")
    df_all = pd.concat(dfs, ignore_index=True)

    for Z, sym in [(104, "Rf"), (106, "Sg"), (108, "Hs"), (110, "Ds")]:
        plot_element(df_all, Z, sym, args.nmin, args.nmax)

    plt.show()


if __name__ == "__main__":
    main()
