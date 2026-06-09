#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make Q_alpha and alpha partial half-life systematics for transfermium nuclei.

Default selection:
    Z = 100--104, N = 148--156

Data source:
    IAEA LiveChart ground_states CSV API. The script uses the evaluated columns
    qa, half_life_sec, decay_1/2/3, and decay_1_%/2_%/3_%.

Outputs:
    raw_livechart_ground_states.csv
    alpha_systematics_Z100_104_N148_156.csv
    alpha_systematics_Z100_104_N148_156.pdf
    alpha_systematics_Z100_104_N148_156.png

Recommended workflow:
    1) Run once with LiveChart auto-download.
    2) Inspect the generated CSV.
    3) Add corrections or preferred ENSDF/NuDat/NUBASE values in overrides.csv.
    4) Re-run with --overrides overrides.csv.
"""

from __future__ import annotations

import argparse
import io
import math
import os
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import pandas as pd

LIVECHART_URL = "https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all"
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 LiveChartAlphaSystematics/1.0"

ELEMENT_Z = {
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
}

# Distinct markers; color is left to matplotlib's default cycle.
MARKERS = {
    100: "o",
    101: "s",
    102: "^",
    103: "D",
    104: "v",
}


def read_livechart_csv():
    url = "https://www-nds.iaea.org/relnsd/v0/data?fields=ground_states&nuclides=all"
    tmp = Path("raw_livechart_ground_states.csv")

    subprocess.run(
        [
            "curl", "-L",
            "-A", "Mozilla/5.0",
            url,
            "-o", str(tmp),
        ],
        check=True,
    )

    return pd.read_csv(tmp)


def to_num(x) -> float:
    """Convert a value to float, tolerating blanks and simple textual flags."""
    if x is None:
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "null"}:
        return np.nan
    # Keep the first number in strings like '<5', '~100', '100 AP'.
    m = re.search(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return np.nan
    return float(m.group(0))


def normalize_decay_mode(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    s = str(x).strip().upper()
    s = s.replace("Α", "A")  # Greek capital alpha if present
    s = s.replace("ALPHA", "A")
    s = s.replace("α", "A")
    return s


def branch_flag_from_string(x) -> str:
    """Return exact/approx/lower/upper/blank from a branch-ratio field."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "blank"
    s = str(x).strip().upper()
    if not s or s == "NAN":
        return "blank"
    if any(tok in s for tok in ["<", "LT"]):
        return "upper"   # b_alpha < value -> T_alpha > T/value
    if any(tok in s for tok in [">", "GT"]):
        return "lower"   # b_alpha > value -> T_alpha < T/value
    if any(tok in s for tok in ["~", "AP", "CA", "≈"]):
        return "approx"
    return "exact"


def find_alpha_branch(row: pd.Series, assume_single_alpha_100: bool = True) -> Tuple[float, float, str, str]:
    """
    Find alpha branching ratio among decay_1/2/3.

    Returns
    -------
    alpha_branch : float
        Branching ratio as fraction, e.g. 0.32 for 32%.
    alpha_branch_unc : float
        Absolute uncertainty as fraction if present.
    alpha_branch_flag : str
        exact, approx, lower, upper, assumed, or missing.
    note : str
        Human-readable note.
    """
    decay_modes = []
    for i in (1, 2, 3):
        mode = normalize_decay_mode(row.get(f"decay_{i}", ""))
        br_raw = row.get(f"decay_{i}_%", np.nan)
        unc_raw = row.get(f"unc_{i}", np.nan)
        if mode:
            decay_modes.append(mode)
        if mode == "A":
            br_percent = to_num(br_raw)
            unc_percent = to_num(unc_raw)
            flag = branch_flag_from_string(br_raw)
            if np.isfinite(br_percent):
                return br_percent / 100.0, (unc_percent / 100.0 if np.isfinite(unc_percent) else np.nan), flag, ""
            if assume_single_alpha_100 and len([m for m in decay_modes if m]) == 1:
                return 1.0, np.nan, "assumed", "alpha branch blank; treated as 100% because alpha is the only listed mode"
            return np.nan, np.nan, "missing", "alpha mode listed but alpha branch ratio is blank"
    return np.nan, np.nan, "missing", "no alpha decay mode listed among decay_1/2/3"


def apply_overrides(df: pd.DataFrame, overrides_path: Optional[Path]) -> pd.DataFrame:
    """Apply optional manual overrides from a CSV file."""
    if overrides_path is None:
        return df
    if not overrides_path.exists():
        raise FileNotFoundError(f"Overrides file not found: {overrides_path}")

    ov = pd.read_csv(overrides_path)
    # Allow either lower-case or upper-case column names by normalizing.
    ov.columns = [c.strip() for c in ov.columns]
    key_candidates = ["Z", "N", "A", "symbol"]
    missing_keys = [c for c in ["Z", "N"] if c not in ov.columns]
    if missing_keys:
        raise ValueError(f"Overrides file must contain at least Z,N columns; missing {missing_keys}")

    df = df.copy()
    df["_key"] = df["Z"].astype(int).astype(str) + "_" + df["N"].astype(int).astype(str)
    ov["_key"] = ov["Z"].astype(int).astype(str) + "_" + ov["N"].astype(int).astype(str)
    ov = ov.set_index("_key")

    override_map = {
        "qa_keV_override": "Qalpha_keV",
        "Qalpha_keV_override": "Qalpha_keV",
        "Qalpha_MeV_override": "Qalpha_MeV",
        "T12_total_s_override": "T12_total_s",
        "half_life_sec_override": "T12_total_s",
        "alpha_branch_override": "alpha_branch",
        "alpha_BR_override": "alpha_branch",
        "alpha_branch_unc_override": "alpha_branch_unc",
        "alpha_branch_flag_override": "alpha_branch_flag",
        "state_override": "state",
        "source_override": "source",
        "note_override": "note",
    }

    for idx, row in df.iterrows():
        key = row["_key"]
        if key not in ov.index:
            continue
        orow = ov.loc[key]
        if isinstance(orow, pd.DataFrame):
            orow = orow.iloc[0]
        for ocol, dcol in override_map.items():
            if ocol in ov.columns:
                value = orow.get(ocol, np.nan)
                if pd.isna(value) or str(value).strip() == "":
                    continue
                if dcol in {"Qalpha_keV", "Qalpha_MeV", "T12_total_s", "alpha_branch", "alpha_branch_unc"}:
                    value = to_num(value)
                    if not np.isfinite(value):
                        continue
                df.at[idx, dcol] = value
        if "Qalpha_MeV" in df.columns and np.isfinite(to_num(df.at[idx, "Qalpha_MeV"])):
            df.at[idx, "Qalpha_keV"] = float(df.at[idx, "Qalpha_MeV"]) * 1000.0
    return df.drop(columns=["_key"])


def prepare_dataframe(raw: pd.DataFrame, zmin: int, zmax: int, nmin: int, nmax: int,
                      assume_single_alpha_100: bool = True) -> pd.DataFrame:
    """Select nuclei and compute alpha partial half-life."""
    df = raw.copy()
    # LiveChart uses lower-case z,n,symbol,qa,half_life_sec.
    required = ["z", "n", "symbol", "qa", "half_life_sec"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"LiveChart data missing expected columns: {missing}\nColumns present: {list(df.columns)}")

    df["z"] = pd.to_numeric(df["z"], errors="coerce")
    df["n"] = pd.to_numeric(df["n"], errors="coerce")
    sel = df[(df["z"].between(zmin, zmax)) & (df["n"].between(nmin, nmax))].copy()
    sel = sel.sort_values(["z", "n"])

    rows = []
    for _, row in sel.iterrows():
        z = int(row["z"])
        n = int(row["n"])
        a = z + n
        symbol = str(row.get("symbol", ELEMENT_Z.get(z, ""))).strip()
        alpha_branch, alpha_branch_unc, alpha_branch_flag, br_note = find_alpha_branch(row, assume_single_alpha_100)
        t12_total_s = to_num(row.get("half_life_sec", np.nan))
        qalpha_keV = to_num(row.get("qa", np.nan))
        unc_qalpha_keV = to_num(row.get("unc_qa", np.nan))
        unc_t12_total_s = to_num(row.get("unc_hls", np.nan))
        state_energy_keV = to_num(row.get("energy", 0.0))

        t12_alpha_s = np.nan
        log10_t12_alpha_s = np.nan
        half_life_limit = ""
        if np.isfinite(t12_total_s) and np.isfinite(alpha_branch) and alpha_branch > 0:
            t12_alpha_s = t12_total_s / alpha_branch
            log10_t12_alpha_s = np.log10(t12_alpha_s)
            if alpha_branch_flag == "upper":
                half_life_limit = "lower"  # b < bmax -> T_alpha > T/bmax
            elif alpha_branch_flag == "lower":
                half_life_limit = "upper"  # b > bmin -> T_alpha < T/bmin

        decay_summary = []
        for i in (1, 2, 3):
            mode = normalize_decay_mode(row.get(f"decay_{i}", ""))
            br = row.get(f"decay_{i}_%", "")
            if mode:
                decay_summary.append(f"{mode}:{br}%")

        rows.append({
            "Z": z,
            "N": n,
            "A": a,
            "symbol": symbol,
            "nuclide": f"{a}{symbol}",
            "state": "g.s." if (not np.isfinite(state_energy_keV) or abs(state_energy_keV) < 1e-9) else f"E={state_energy_keV:g} keV",
            "Qalpha_keV": qalpha_keV,
            "unc_Qalpha_keV": unc_qalpha_keV,
            "Qalpha_MeV": qalpha_keV / 1000.0 if np.isfinite(qalpha_keV) else np.nan,
            "T12_total_s": t12_total_s,
            "unc_T12_total_s": unc_t12_total_s,
            "alpha_branch": alpha_branch,
            "alpha_branch_unc": alpha_branch_unc,
            "alpha_branch_percent": alpha_branch * 100.0 if np.isfinite(alpha_branch) else np.nan,
            "alpha_branch_flag": alpha_branch_flag,
            "T12_alpha_s": t12_alpha_s,
            "log10_T12_alpha_s": log10_t12_alpha_s,
            "T12_alpha_limit": half_life_limit,
            "decay_summary": "; ".join(decay_summary),
            "half_life_operator": row.get("operator_hl", ""),
            "source": "IAEA LiveChart ground_states",
            "ENSDFpublicationcut-off": row.get("ENSDFpublicationcut-off", row.get("ENSDF_publication_cut-off", "")),
            "ENSDFauthors": row.get("ENSDFauthors", row.get("ENSDF_authors", "")),
            "Extraction_date": row.get("Extraction_date", ""),
            "note": br_note,
        })
    return pd.DataFrame(rows)


def recompute_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute Qalpha_MeV, alpha_percent, and T_alpha after overrides."""
    df = df.copy()
    df["Qalpha_MeV"] = pd.to_numeric(df["Qalpha_keV"], errors="coerce") / 1000.0
    df["alpha_branch"] = pd.to_numeric(df["alpha_branch"], errors="coerce")
    df["T12_total_s"] = pd.to_numeric(df["T12_total_s"], errors="coerce")
    df["alpha_branch_percent"] = df["alpha_branch"] * 100.0
    df["T12_alpha_s"] = np.where(
        (df["T12_total_s"] > 0) & (df["alpha_branch"] > 0),
        df["T12_total_s"] / df["alpha_branch"],
        np.nan,
    )
    df["log10_T12_alpha_s"] = np.where(df["T12_alpha_s"] > 0, np.log10(df["T12_alpha_s"]), np.nan)
    return df


def setup_matplotlib() -> None:
    """Simple journal-style matplotlib defaults."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "savefig.bbox": "tight",
        "savefig.dpi": 450,
    })


def plot_systematics(df: pd.DataFrame, out_prefix: Path, label_points: bool = False) -> None:
    """Make two-panel Q_alpha and alpha partial half-life systematics plot."""
    setup_matplotlib()
    fig, (ax_q, ax_t) = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True, gridspec_kw={"hspace": 0.08})

    for z, g in df.groupby("Z", sort=True):
        g = g.sort_values("N")
        label = f"Z={z} ({ELEMENT_Z.get(int(z), g['symbol'].iloc[0])})"
        marker = MARKERS.get(int(z), "o")

        qmask = np.isfinite(g["Qalpha_MeV"])
        if qmask.any():
            ax_q.plot(g.loc[qmask, "N"], g.loc[qmask, "Qalpha_MeV"],
                      linestyle="--", linewidth=1.0, alpha=0.65)
            ax_q.scatter(g.loc[qmask, "N"], g.loc[qmask, "Qalpha_MeV"],
                         marker=marker, s=46, edgecolors="black", linewidths=0.7, label=label, zorder=3)
            # Optional Qalpha error bars if uncertainties are available.
            err = pd.to_numeric(g.loc[qmask, "unc_Qalpha_keV"], errors="coerce") / 1000.0
            if np.isfinite(err).any():
                ax_q.errorbar(g.loc[qmask, "N"], g.loc[qmask, "Qalpha_MeV"], yerr=err,
                              fmt="none", linewidth=0.8, capsize=2.5, zorder=2)

        tmask = np.isfinite(g["log10_T12_alpha_s"])
        if tmask.any():
            ax_t.plot(g.loc[tmask, "N"], g.loc[tmask, "log10_T12_alpha_s"],
                      linestyle="--", linewidth=1.0, alpha=0.65)
            exact = g.loc[tmask & (g["T12_alpha_limit"].fillna("") == "")]
            lower = g.loc[tmask & (g["T12_alpha_limit"] == "lower")]
            upper = g.loc[tmask & (g["T12_alpha_limit"] == "upper")]
            if len(exact):
                ax_t.scatter(exact["N"], exact["log10_T12_alpha_s"],
                             marker=marker, s=46, edgecolors="black", linewidths=0.7, zorder=3)
            if len(lower):
                ax_t.scatter(lower["N"], lower["log10_T12_alpha_s"],
                             marker="^", s=54, edgecolors="black", linewidths=0.7, zorder=3)
            if len(upper):
                ax_t.scatter(upper["N"], upper["log10_T12_alpha_s"],
                             marker="v", s=54, edgecolors="black", linewidths=0.7, zorder=3)

        if label_points:
            for _, r in g.iterrows():
                if np.isfinite(r.get("Qalpha_MeV", np.nan)):
                    ax_q.text(r["N"] + 0.03, r["Qalpha_MeV"], r["nuclide"], fontsize=7, va="center")
                if np.isfinite(r.get("log10_T12_alpha_s", np.nan)):
                    ax_t.text(r["N"] + 0.03, r["log10_T12_alpha_s"], r["nuclide"], fontsize=7, va="center")

    for ax in (ax_q, ax_t):
        ax.axvline(152, linestyle=":", linewidth=1.2)
        ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.45)
        ax.minorticks_on()

    ax_q.set_ylabel(r"$Q_{\alpha}$ (MeV)")
    ax_t.set_ylabel(r"$\log_{10}\,T_{1/2,\alpha}$ (s)")
    ax_t.set_xlabel(r"Neutron number $N$")
    ax_q.legend(frameon=False, ncol=3, loc="best")

    ax_q.text(0.02, 0.94, r"(a)", transform=ax_q.transAxes, fontweight="bold")
    ax_t.text(0.02, 0.94, r"(b)", transform=ax_t.transAxes, fontweight="bold")
    ax_q.text(152.1, ax_q.get_ylim()[1] - 0.06 * (ax_q.get_ylim()[1] - ax_q.get_ylim()[0]),
              r"$N=152$", fontsize=10, va="top")

    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"))
    plt.close(fig)


def write_template(path: Path) -> None:
    text = """Z,N,A,symbol,nuclide,qa_keV_override,T12_total_s_override,alpha_branch_override,alpha_branch_unc_override,alpha_branch_flag_override,state_override,source_override,note_override
100,152,252,Fm,252Fm,,,,,,,"ENSDF/NuDat manual check","example row; leave cells blank unless overriding"
102,152,254,No,254No,,,,,,,"ENSDF/NuDat manual check","example row; alpha_branch is fraction, e.g. 0.90 for 90%"
104,152,256,Rf,256Rf,,,,,,,"ENSDF/NuDat manual check","branch flags: exact, approx, lower, upper"
"""
    path.write_text(text, encoding="utf-8")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Plot Q_alpha and alpha partial half-life systematics.")
    parser.add_argument("--zmin", type=int, default=100)
    parser.add_argument("--zmax", type=int, default=104)
    parser.add_argument("--nmin", type=int, default=148)
    parser.add_argument("--nmax", type=int, default=156)
    parser.add_argument("--outdir", default="alpha_systematics_out", help="Output directory")
    parser.add_argument("--raw-csv", default=None, help="Use an existing LiveChart ground_states CSV instead of downloading")
    parser.add_argument("--overrides", default=None, help="Optional manual override CSV")
    parser.add_argument("--no-assume-single-alpha-100", action="store_true",
                        help="Do not assume 100%% alpha branch when alpha is the only listed mode but branch is blank")
    parser.add_argument("--label-points", action="store_true", help="Annotate each data point by nuclide")
    parser.add_argument("--write-template", action="store_true", help="Only write overrides_template.csv and exit")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.write_template:
        write_template(outdir / "overrides_template.csv")
        print(f"Wrote {outdir / 'overrides_template.csv'}")
        return 0

    if args.raw_csv:
        raw_path = Path(args.raw_csv)
        raw = pd.read_csv(raw_path)
    else:
        raw = read_livechart_csv()
        raw_path = outdir / "raw_livechart_ground_states.csv"
        raw.to_csv(raw_path, index=False)
        print(f"Saved raw LiveChart data: {raw_path}")

    df = prepare_dataframe(
        raw,
        zmin=args.zmin,
        zmax=args.zmax,
        nmin=args.nmin,
        nmax=args.nmax,
        assume_single_alpha_100=not args.no_assume_single_alpha_100,
    )

    if args.overrides:
        df = apply_overrides(df, Path(args.overrides))
        df = recompute_derived_columns(df)

    prefix = f"alpha_systematics_Z{args.zmin}_{args.zmax}_N{args.nmin}_{args.nmax}"
    csv_path = outdir / f"{prefix}.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved processed table: {csv_path}")

    out_prefix = outdir / prefix
    plot_systematics(df, out_prefix, label_points=args.label_points)
    print(f"Saved figures: {out_prefix.with_suffix('.pdf')} and {out_prefix.with_suffix('.png')}")

    missing_alpha = df[df["alpha_branch"].isna()][["nuclide", "decay_summary", "note"]]
    if len(missing_alpha):
        print("\nEntries without numeric alpha branch were not included in the alpha partial half-life panel:")
        print(missing_alpha.to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
