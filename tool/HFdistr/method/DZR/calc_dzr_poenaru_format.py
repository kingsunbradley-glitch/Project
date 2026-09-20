#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate alpha-decay half-lives with the original Royer formula
and the Deng-Zhang-Royer (DZR, PRC 101, 034307, 2020) improved formula.

Designed to be compatible with the Poenaru-style/simple input.csv format, e.g.

Name,A,Z,Ealpha_MeV,Texp_s,Br_alpha
273Ds_a,273,110,10.858,0.00832,1.0

Also accepts aliases:
  Ealpha: Ealpha, E_alpha, Ealpha_MeV, E_alpha_MeV, EXP, EXP_MeV, Ea, Ea_MeV
  Qalpha: Qalpha, Q_alpha, Qalpha_MeV, Q_alpha_MeV, Q, Q_MeV
  Texp: Texp_s, T_exp_s, T0.5, T12_s, half_life_s
  l: l, L, DeltaL, Delta_L, lmin
  branch: Br_alpha, Br, Branch, b_alpha

Usage:
  python calc_dzr_poenaru_format_fixed.py input.csv -o output_dzr.csv

By default, if only Ealpha is given, Qalpha = Ealpha * A/(A-4).
Use --exp-is-qalpha if your EXP/Ealpha column is already Qalpha.
"""

import argparse
import math
from pathlib import Path
import pandas as pd

EALPHA_ALIASES = [
    "Ealpha", "E_alpha", "Ealpha_MeV", "E_alpha_MeV", "E_ALPHA_MEV",
    "Ea", "Ea_MeV", "EXP", "EXP_MeV", "Eexp", "E_exp", "Eexp_MeV",
]
QALPHA_ALIASES = [
    "Qalpha", "Q_alpha", "Qalpha_MeV", "Q_alpha_MeV", "Q_ALPHA_MEV",
    "Q", "Q_MeV", "QMeV",
]
TEXP_ALIASES = [
    "Texp_s", "T_exp_s", "Texp", "T_exp", "T0.5", "T1/2", "T12_s",
    "half_life_s", "HalfLife_s", "tau_s", "time_s",
]
L_ALIASES = ["l", "L", "ell", "DeltaL", "Delta_L", "DL", "lmin", "l_min"]
BR_ALIASES = ["Br_alpha", "Br", "BR", "Branch", "b_alpha", "branch_alpha"]

# Original Royer Eq. (2) parameter sets
ROYER_PARAMS = {
    "even-even": (-25.31, -1.1629, 1.5864),
    "evenZ-oddN": (-26.65, -1.0859, 1.5848),
    "oddZ-evenN": (-25.68, -1.1423, 1.5920),
    "odd-odd": (-29.48, -1.1130, 1.6971),
}

# Deng-Zhang-Royer Eq. (7) global parameters and blocking h
DZR_A = -26.8125
DZR_B = -1.1255
DZR_C = 1.6057
DZR_D = 0.0513
DZR_H = {
    "even-even": 0.0,
    "evenZ-oddN": 0.3625,
    "oddZ-evenN": 0.2812,
    "odd-odd": 0.7486,
}


def norm(s: str) -> str:
    """Normalize a column name for robust matching."""
    return str(s).strip().lower().replace(" ", "").replace("-", "_")


def find_col(df: pd.DataFrame, aliases):
    table = {norm(c): c for c in df.columns}
    for a in aliases:
        key = norm(a)
        if key in table:
            return table[key]
    return None


def as_float(x, default=math.nan):
    if pd.isna(x):
        return default
    if isinstance(x, str):
        x = x.strip()
        if x == "":
            return default
    return float(x)


def case_from_ZN(Z: int, N: int) -> str:
    if Z % 2 == 0 and N % 2 == 0:
        return "even-even"
    if Z % 2 == 0 and N % 2 == 1:
        return "evenZ-oddN"
    if Z % 2 == 1 and N % 2 == 0:
        return "oddZ-evenN"
    return "odd-odd"


def logT_royer(A: int, Z: int, Qalpha: float, case: str) -> float:
    a, b, c = ROYER_PARAMS[case]
    return a + b * (A ** (1.0 / 6.0)) * math.sqrt(Z) + c * Z / math.sqrt(Qalpha)


def logT_dzr(A: int, Z: int, Qalpha: float, ell: float, case: str) -> float:
    return (
        DZR_A
        + DZR_B * (A ** (1.0 / 6.0)) * math.sqrt(Z)
        + DZR_C * Z / math.sqrt(Qalpha)
        + DZR_D * ell * (ell + 1.0)
        + DZR_H[case]
    )


def main():
    parser = argparse.ArgumentParser(
        description="Calculate alpha-decay half-lives with Royer and DZR formulas."
    )
    parser.add_argument("input_csv", help="Input CSV file")
    parser.add_argument("-o", "--output", default="output_dzr.csv", help="Output CSV file")
    parser.add_argument(
        "--exp-is-qalpha",
        action="store_true",
        help="Treat the EXP/Ealpha column as Qalpha directly, not as measured Ealpha.",
    )
    args = parser.parse_args()

    path = Path(args.input_csv)
    if not path.exists():
        raise FileNotFoundError(f"Cannot find input file: {path}")

    df = pd.read_csv(path)
    df.columns = [str(c).strip().replace("\ufeff", "") for c in df.columns]

    A_col = find_col(df, ["A", "Mass", "MassNumber"])
    Z_col = find_col(df, ["Z", "Proton", "ProtonNumber"])
    N_col = find_col(df, ["N", "Neutron", "NeutronNumber"])
    E_col = find_col(df, EALPHA_ALIASES)
    Q_col = find_col(df, QALPHA_ALIASES)
    T_col = find_col(df, TEXP_ALIASES)
    l_col = find_col(df, L_ALIASES)
    br_col = find_col(df, BR_ALIASES)

    if A_col is None or Z_col is None:
        raise ValueError(f"Need A and Z columns. Found columns: {list(df.columns)}")
    if Q_col is None and E_col is None:
        raise ValueError(
            "Need either Qalpha column or Ealpha/EXP column. "
            f"Found columns: {list(df.columns)}\n"
            f"Recognized Ealpha names include: {EALPHA_ALIASES}\n"
            f"Recognized Qalpha names include: {QALPHA_ALIASES}"
        )

    out_rows = []
    for _, row in df.iterrows():
        A = int(round(as_float(row[A_col])))
        Z = int(round(as_float(row[Z_col])))
        N = int(round(as_float(row[N_col]))) if N_col else A - Z
        case = case_from_ZN(Z, N)

        if Q_col is not None:
            Qalpha = as_float(row[Q_col])
        else:
            Ealpha = as_float(row[E_col])
            if E_col and "kev" in norm(E_col):
                Ealpha /= 1000.0
            Qalpha = Ealpha if args.exp_is_qalpha else Ealpha * A / (A - 4.0)

        ell = 0.0 if l_col is None or pd.isna(row[l_col]) else as_float(row[l_col], 0.0)

        logR = logT_royer(A, Z, Qalpha, case)
        logD = logT_dzr(A, Z, Qalpha, ell, case)
        TR = 10.0 ** logR
        TD = 10.0 ** logD

        # Experimental alpha partial half-life: Talpha = T_total / Br_alpha.
        Talpha_exp = math.nan
        HFR = math.nan
        HFD = math.nan
        if T_col is not None and not pd.isna(row[T_col]):
            Texp = as_float(row[T_col])
            br = 1.0
            if br_col is not None and not pd.isna(row[br_col]):
                br = as_float(row[br_col], 1.0)
            if br > 0:
                Talpha_exp = Texp / br
                HFR = Talpha_exp / TR
                HFD = Talpha_exp / TD

        extra = {
            "N_calc": N,
            "Qalpha_DZR_MeV": Qalpha,
            "l_used": ell,
            "parity_case": case,
            "log10T_Royer_s": logR,
            "T_Royer_s": TR,
            "log10T_DZR_s": logD,
            "T_DZR_s": TD,
            "Talpha_exp_s": Talpha_exp,
            "HF_Royer": HFR,
            "HF_DZR": HFD,
        }
        out_rows.append(extra)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(out_rows)], axis=1)
    out.to_csv(args.output, index=False)
    print(f"Input columns: {list(df.columns)}")
    print(f"Detected: A={A_col}, Z={Z_col}, N={N_col}, Ealpha={E_col}, Qalpha={Q_col}, Texp={T_col}, Br={br_col}, l={l_col}")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
