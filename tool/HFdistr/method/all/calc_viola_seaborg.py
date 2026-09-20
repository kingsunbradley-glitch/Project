#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_viola_seaborg.py

Compute alpha-decay half-lives with the Viola-Seaborg formula in the form

    log10(T_alpha/s) = (a*Z + b) / sqrt(Q_alpha) + c*Z + d

with parameters from the user's screenshot:
    a = 1.787
    b = -21.40
    c = -0.2549
    d = -28.42

Input: CSV compatible with the Poenaru-style file used in this project, e.g.

    Name,A,Z,Ealpha_MeV,Texp_s,Br_alpha
    273Ds_a,273,110,10.858,0.00832,1.0

The script preserves all original columns and appends calculated columns.
If Q_alpha is not present, it converts E_alpha to Q_alpha using

    Q_alpha = E_alpha * A / (A - 4)

unless --exp-is-qalpha is given, in which case the Ealpha/EXP column is used
as Q_alpha directly.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

# Viola-Seaborg parameters from the screenshot
VS_A = 1.787
VS_B = -21.40
VS_C = -0.2549
VS_D = -28.42


ALIASES = {
    "A": ["A", "Mass", "MassNumber", "A_parent", "Ap"],
    "Z": ["Z", "Z_parent", "Zp", "Proton", "ProtonNumber"],
    "N": ["N", "N_parent", "Np", "Neutron", "NeutronNumber"],
    "Qalpha": ["Qalpha", "Q_alpha", "Qalpha_MeV", "Q_alpha_MeV", "Q", "Q_MeV"],
    "Ealpha": ["Ealpha", "E_alpha", "Ealpha_MeV", "E_alpha_MeV", "EXP", "E", "E_MeV"],
    "Texp": ["Texp_s", "Talpha_exp_s", "T_exp_s", "T0.5", "T12_s", "T1/2_s", "tau_s", "life_s"],
    "Br": ["Br_alpha", "Bralpha", "BR_alpha", "branch", "Branch", "Br"],
}


def find_col(df: pd.DataFrame, key: str) -> Optional[str]:
    """Find a column using exact and case-insensitive aliases."""
    cols = list(df.columns)
    for cand in ALIASES[key]:
        if cand in cols:
            return cand
    lower_map = {str(c).strip().lower(): c for c in cols}
    for cand in ALIASES[key]:
        hit = lower_map.get(cand.lower())
        if hit is not None:
            return hit
    return None


def as_float(value, colname: str, row_index: int) -> float:
    try:
        x = float(value)
    except Exception as exc:
        raise ValueError(f"Row {row_index}: cannot convert {colname}={value!r} to float") from exc
    if not math.isfinite(x):
        raise ValueError(f"Row {row_index}: {colname} is not finite: {value!r}")
    return x


def viola_seaborg_log10T(A: float, Z: float, Qalpha: float) -> float:
    """Return log10(T_alpha/s) using the Viola-Seaborg formula."""
    if A <= 4:
        raise ValueError(f"A must be > 4, got A={A}")
    if Z <= 0:
        raise ValueError(f"Z must be > 0, got Z={Z}")
    if Qalpha <= 0:
        raise ValueError(f"Qalpha must be > 0 MeV, got Qalpha={Qalpha}")
    return (VS_A * Z + VS_B) / math.sqrt(Qalpha) + VS_C * Z + VS_D


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute alpha half-lives with the Viola-Seaborg formula, preserving input.csv format."
    )
    parser.add_argument("input", help="Input CSV file, e.g. input.csv")
    parser.add_argument("-o", "--output", default="output_viola_seaborg.csv", help="Output CSV file")
    parser.add_argument(
        "--exp-is-qalpha",
        action="store_true",
        help="Treat the Ealpha/EXP column as Qalpha directly, not as measured alpha-particle energy.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print detected columns and output path.",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    df = pd.read_csv(in_path)

    col_A = find_col(df, "A")
    col_Z = find_col(df, "Z")
    col_N = find_col(df, "N")
    col_Q = find_col(df, "Qalpha")
    col_E = find_col(df, "Ealpha")
    col_T = find_col(df, "Texp")
    col_Br = find_col(df, "Br")

    if col_A is None or col_Z is None:
        raise ValueError(f"Need A and Z columns. Found columns: {list(df.columns)}")
    if col_Q is None and col_E is None:
        raise ValueError(
            "Need either a Qalpha column or an Ealpha/EXP column. "
            "Supported Ealpha names include Ealpha_MeV, Ealpha, EXP."
        )

    if not args.quiet:
        print(
            "Detected columns: "
            f"A={col_A}, Z={col_Z}, N={col_N}, Qalpha={col_Q}, "
            f"Ealpha={col_E}, Texp={col_T}, Br={col_Br}"
        )

    out_rows = []
    for i, row in df.iterrows():
        A = as_float(row[col_A], col_A, i)
        Z = as_float(row[col_Z], col_Z, i)
        N_calc = as_float(row[col_N], col_N, i) if col_N is not None else A - Z

        if col_Q is not None:
            Qalpha = as_float(row[col_Q], col_Q, i)
            q_source = col_Q
        else:
            Ealpha = as_float(row[col_E], col_E, i)
            if args.exp_is_qalpha:
                Qalpha = Ealpha
                q_source = f"{col_E} treated as Qalpha"
            else:
                Qalpha = Ealpha * A / (A - 4.0)
                q_source = f"{col_E} converted by A/(A-4)"

        log10T = viola_seaborg_log10T(A, Z, Qalpha)
        Tcalc = 10.0 ** log10T

        result = {
            "N_calc": N_calc,
            "Qalpha_VS_MeV": Qalpha,
            "Qalpha_source": q_source,
            "log10T_VS_s": log10T,
            "T_VS_s": Tcalc,
        }

        if col_T is not None:
            Texp = as_float(row[col_T], col_T, i)
            Br = as_float(row[col_Br], col_Br, i) if col_Br is not None else 1.0
            if Br <= 0:
                raise ValueError(f"Row {i}: Br_alpha must be > 0, got {Br}")
            Talpha_exp = Texp / Br
            result["Talpha_exp_s"] = Talpha_exp
            result["log10Talpha_exp_s"] = math.log10(Talpha_exp) if Talpha_exp > 0 else float("nan")
            result["Delta_log10_exp_minus_VS"] = result["log10Talpha_exp_s"] - log10T if Talpha_exp > 0 else float("nan")
            result["HF_VS"] = Talpha_exp / Tcalc if Tcalc > 0 else float("nan")
        out_rows.append(result)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(out_rows)], axis=1)
    out.to_csv(args.output, index=False)

    if not args.quiet:
        print(f"Wrote: {args.output}")
        preview_cols = [c for c in ["Name", col_A, col_Z, col_E, col_Q, "Qalpha_VS_MeV", "log10T_VS_s", "T_VS_s", "HF_VS"] if c and c in out.columns]
        if preview_cols:
            print(out[preview_cols].to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
