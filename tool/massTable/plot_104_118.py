#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Qα vs N trends for Z=104–118.
- 3 models: line + small markers
- EXP: scatter with error bars
- One element (one Z) per figure, and each figure is shown via plt.show() (no saving).

Default input filenames (in current working directory):
  - Q@a_WS4+RBF.dat
  - Q@a_FRDM.dat
  - Q@a_HFB.dat
  - Q@a_EXP.dat

Usage:
  python plot_Qalpha_Z104_118_show_per_element.py
  python plot_Qalpha_Z104_118_show_per_element.py --nmin 145 --nmax 190
"""

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


SYMBOLS = {
    104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs",
    109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh",
    114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
}


def _try_parse_numeric(parts):
    try:
        return [float(x) for x in parts]
    except Exception:
        return None


def read_qalpha_table(path: str, model_name: str) -> pd.DataFrame:
    """
    Robust reader for model files. Accepts either:
      - Z N A Q(MeV) (>=4 numbers)
      - A Z Q(MeV)   (3 numbers)  -> N inferred as A-Z
    Ignores comment/header lines.
    """
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#") or s.startswith("[") or s.startswith("-"):
                    continue
                # Allow lines like: "A   Z   Q_alpha" -> skip (not numeric start)
                if not s[0].isdigit():
                    continue
                parts = re.split(r"\s+", s)
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
    """
    EXP file expected columns (numbers per row):
      A Z N Q Error
    Q and Error are assumed in keV and converted to MeV (as in your reference script).
    """
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#") or s.startswith("[") or s.startswith("-"):
                    continue
                if not s[0].isdigit():
                    continue
                parts = re.split(r"\s+", s)
                nums = _try_parse_numeric(parts)
                if nums is None or len(nums) < 5:
                    continue
                A = int(round(nums[0]))
                Z = int(round(nums[1]))
                N = int(round(nums[2]))
                Q_kev = float(nums[3])
                Err_kev = float(nums[4])
                rows.append((Z, N, A, Q_kev / 1000.0, Err_kev / 1000.0))
    except FileNotFoundError:
        print(f"[WARN] File not found: {path}")
        return pd.DataFrame()

    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha", "Error"])
    df["model"] = "EXP"
    return df


def plot_one_element(Z: int,
                     df_models: pd.DataFrame,
                     df_exp: pd.DataFrame,
                     nmin: int | None,
                     nmax: int | None,
                     marker_model: float,
                     marker_exp: float,
                     lw: float):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    sym = SYMBOLS.get(Z, "")
    title = f"Z={Z}" + (f" ({sym})" if sym else "")
    ax.set_title(title)

    dm = df_models[df_models["Z"] == Z].copy()
    if nmin is not None:
        dm = dm[dm["N"] >= nmin]
    if nmax is not None:
        dm = dm[dm["N"] <= nmax]

    for model_name in sorted(dm["model"].dropna().unique()):
        d = dm[dm["model"] == model_name].sort_values("N")
        if not d.empty:
            ax.plot(
                d["N"], d["Qalpha"],
                marker="o",
                markersize=marker_model,
                linewidth=lw,
                alpha=0.85,
                label=model_name
            )

    de = df_exp[df_exp["Z"] == Z].copy()
    if nmin is not None:
        de = de[de["N"] >= nmin]
    if nmax is not None:
        de = de[de["N"] <= nmax]

    if not de.empty:
        de = de.sort_values("N")
        ax.errorbar(
            de["N"], de["Qalpha"], yerr=de["Error"],
            fmt="s",
            color="black",
            markersize=marker_exp,
            capsize=3,
            linestyle="None",
            label="EXP",
            zorder=10,
        )

    ax.set_xlabel("Neutron number N")
    ax.set_ylabel(r"$Q_\alpha$ (MeV)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.8)
    ax.legend(frameon=False)

    plt.tight_layout()
    plt.show()   # <-- per-element show
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Show Qalpha vs N per element (Z=104–118), no saving.")
    ap.add_argument("--frdm", default="Q@a_FRDM.dat")
    ap.add_argument("--hfb",  default="Q@a_HFB.dat")
    ap.add_argument("--ws4",  default="Q@a_WS4+RBF.dat")
    ap.add_argument("--exp",  default="Q@a_EXP.dat")
    ap.add_argument("--zmin", type=int, default=104)
    ap.add_argument("--zmax", type=int, default=118)
    ap.add_argument("--nmin", type=int, default=None)
    ap.add_argument("--nmax", type=int, default=None)

    # Marker/line style knobs (match your reference style but smaller if you want)
    ap.add_argument("--ms-model", type=float, default=3.0, help="marker size for model curves")
    ap.add_argument("--ms-exp",   type=float, default=4.0, help="marker size for EXP points")
    ap.add_argument("--lw",       type=float, default=1.2, help="line width for model curves")
    args = ap.parse_args()

    df_frdm = read_qalpha_table(args.frdm, "FRDM(2012)")
    df_hfb  = read_qalpha_table(args.hfb,  "HFB-24")
    df_ws4  = read_qalpha_table(args.ws4,  "WS4+RBF")
    df_exp  = read_exp_table(args.exp)

    dfs = [df for df in (df_frdm, df_hfb, df_ws4) if not df.empty]
    df_models = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(
        columns=["Z", "N", "A", "Qalpha", "model"]
    )

    # Loop: one element per figure, each plt.show()
    for Z in range(args.zmin, args.zmax + 1):
        has_any = False
        if not df_models.empty and (df_models["Z"] == Z).any():
            has_any = True
        if not df_exp.empty and (df_exp["Z"] == Z).any():
            has_any = True
        if not has_any:
            continue
        plot_one_element(
            Z,
            df_models=df_models,
            df_exp=df_exp,
            nmin=args.nmin,
            nmax=args.nmax,
            marker_model=args.ms_model,
            marker_exp=args.ms_exp,
            lw=args.lw,
        )

    print("[OK] Done. Closed windows? Script finished.")


if __name__ == "__main__":
    main()
