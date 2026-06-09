#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calc_poenaru_target_format.py

Calculate alpha-decay partial half-lives using the Poenaru-type formula
as written in the target paper screenshot:

    T_1/2^alpha = 10^[(B1+B2*y+B3*z+B4*y^2+B5*y*z+B6*z^2)*Ks/ln(10) - 20.446]  seconds

with

    y = (N - Ni)/(N_{i+1} - Ni),  Ni < N <= N_{i+1}
    z = (Z - Zi)/(Z_{i+1} - Zi),  Zi < Z <= Z_{i+1}

    Ni = 51, 83, 127, 185, ...
    Zi = 51, 83, 115, 121, ...

    Q = E_alpha * A / A_d,  A_d = A - 4

    Ks = 2.52956 * Z_d * sqrt(A_d/(A*Q)) *
         [acos(sqrt(x)) - sqrt(x*(1-x))]

    x = 0.4253 * Q * (1.5874 + A_d^(1/3)) / Z_d
    Z_d = Z - 2

Target-paper parameters:
    B1 = 0.988662
    B2 = 0.016314
    B3 = 0.020433
    B4 = 0.027896
    B5 = -0.003033
    B6 = -0.003033

Input:
    csv / xlsx / xls

Required columns, with flexible aliases:
    A       : A, Mass, mass
    Z       : Z, Proton, proton
    energy  : Qalpha_MeV / Qalpha / Q / Q_MeV
              or Ealpha_MeV / Ealpha / E / Ea / E_alpha

Optional experimental half-life columns:
    Texp_s / T12_s / T1/2_s / t12_s / half_life_s / T
    Br_alpha / Br / branch / b_alpha

Important:
    If only Ealpha is given, this script converts it to Qalpha by
        Qalpha = Ealpha * A/(A-4)
    as in the target paper.

Example:
    python calc_poenaru_target_format.py input.csv -o output.csv
    python calc_poenaru_target_format.py input.csv --energy-mode Ealpha
    python calc_poenaru_target_format.py input.csv --param-set original1980
"""

import argparse
import math
import os
from typing import Optional

import numpy as np
import pandas as pd


N_MAGIC_PLUS_ONE = [51, 83, 127, 185, 229, 283]
Z_MAGIC_PLUS_ONE = [51, 83, 115, 121, 173, 185]

PARAM_SETS = {
    # Parameters exactly corresponding to the target-paper text shown by the user:
    "target": {
        "B1": 0.988662,
        "B2": 0.016314,
        "B3": 0.020433,
        "B4": 0.027896,
        "B5": -0.003033,
        "B6": -0.003033,
    },
    # Kept only for comparison/debugging if the older script was based on another table.
    "original1980": {
        "B1": 0.988662,
        "B2": 0.016314,
        "B3": 0.020433,
        "B4": 0.027896,
        "B5": -0.003033,
        "B6": -0.16820,
    },
}


ALIASES = {
    "A": ["A", "Mass", "mass", "A_parent", "Ap"],
    "Z": ["Z", "Proton", "proton", "Z_parent", "Zp"],
    "Qalpha": ["Qalpha_MeV", "Qalpha", "Q_alpha", "Q", "Q_MeV", "Qa", "Q_a"],
    "Ealpha": ["Ealpha_MeV", "Ealpha", "E_alpha", "E", "Ea", "E_a", "Ealpha_keV", "E_alpha_keV"],
    "Texp_s": ["Texp_s", "T12_s", "T1/2_s", "t12_s", "half_life_s", "T", "T_s"],
    "Br_alpha": ["Br_alpha", "Br", "BR", "branch", "Branch", "b_alpha", "AlphaBranch"],
}


def read_table(filename: str) -> pd.DataFrame:
    ext = os.path.splitext(filename)[1].lower()
    if ext == ".csv":
        return pd.read_csv(filename)
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(filename)
    raise ValueError("Input file must be .csv, .xlsx, or .xls")


def find_column(df: pd.DataFrame, logical_name: str) -> Optional[str]:
    """Return the first actual column matching a logical name alias."""
    candidates = ALIASES[logical_name]

    # Exact match first
    for c in candidates:
        if c in df.columns:
            return c

    # Case-insensitive fallback, ignoring surrounding spaces
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for c in candidates:
        key = c.strip().lower()
        if key in lower_map:
            return lower_map[key]

    return None


def find_reduced_variable(value: int, magic_list: list[int]) -> tuple[float, int, int]:
    for i in range(len(magic_list) - 1):
        low = magic_list[i]
        high = magic_list[i + 1]
        if low < value <= high:
            return (value - low) / (high - low), low, high

    raise ValueError(
        f"value={value} is outside implemented magic-plus-one intervals: {magic_list}"
    )


def qalpha_from_ealpha(A: int, Ealpha_MeV: float) -> float:
    Ad = A - 4
    return Ealpha_MeV * A / Ad


def normalize_energy_value(colname: str, value: float) -> float:
    """
    Convert energy to MeV if the column name clearly says keV.
    Otherwise assume MeV.
    """
    if "kev" in colname.lower():
        return value / 1000.0
    return value


def calc_poenaru(A: int, Z: int, Qalpha_MeV: float, params: dict) -> dict:
    N = A - Z
    Ad = A - 4
    Zd = Z - 2

    y, Ni, Nip1 = find_reduced_variable(N, N_MAGIC_PLUS_ONE)
    z, Zi, Zip1 = find_reduced_variable(Z, Z_MAGIC_PLUS_ONE)

    x = 0.4253 * Qalpha_MeV * (1.5874 + Ad ** (1.0 / 3.0)) / Zd
    if not (0.0 < x < 1.0):
        raise ValueError(
            f"Invalid x={x:.8g}; require 0<x<1. "
            f"Check A={A}, Z={Z}, Qalpha_MeV={Qalpha_MeV}."
        )

    bracket = math.acos(math.sqrt(x)) - math.sqrt(x * (1.0 - x))
    Ks = 2.52956 * Zd * math.sqrt(Ad / (A * Qalpha_MeV)) * bracket

    B1 = params["B1"]
    B2 = params["B2"]
    B3 = params["B3"]
    B4 = params["B4"]
    B5 = params["B5"]
    B6 = params["B6"]

    F_yz = B1 + B2*y + B3*z + B4*y*y + B5*y*z + B6*z*z
    log10_Tcalc_s = F_yz * Ks / math.log(10.0) - 20.446
    Tcalc_s = 10.0 ** log10_Tcalc_s

    return {
        "N_calc": N,
        "Ad_calc": Ad,
        "Zd_calc": Zd,
        "Qalpha_calc_MeV": Qalpha_MeV,
        "x_calc": x,
        "Ks_calc": Ks,
        "y_calc": y,
        "z_calc": z,
        "N_interval": f"({Ni},{Nip1}]",
        "Z_interval": f"({Zi},{Zip1}]",
        "F_yz_calc": F_yz,
        "log10_T12a_Poenaru_s": log10_Tcalc_s,
        "T12a_Poenaru_s": Tcalc_s,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate Poenaru alpha-decay half-lives in the target-paper format."
    )
    parser.add_argument("input", help="Input table: .csv, .xlsx, or .xls")
    parser.add_argument("-o", "--output", default=None, help="Output CSV filename")
    parser.add_argument(
        "--param-set",
        choices=sorted(PARAM_SETS.keys()),
        default="target",
        help="Parameter set. Default: target, i.e. B5=B6=-0.003033.",
    )
    parser.add_argument(
        "--energy-mode",
        choices=["auto", "Qalpha", "Ealpha"],
        default="auto",
        help=(
            "auto: prefer Qalpha column if present; otherwise use Ealpha. "
            "Qalpha: force Qalpha aliases. Ealpha: force Ealpha aliases and convert to Q."
        ),
    )
    args = parser.parse_args()

    df = read_table(args.input)

    A_col = find_column(df, "A")
    Z_col = find_column(df, "Z")
    if A_col is None:
        raise ValueError(f"Missing A column. Accepted aliases: {ALIASES['A']}")
    if Z_col is None:
        raise ValueError(f"Missing Z column. Accepted aliases: {ALIASES['Z']}")

    Q_col = find_column(df, "Qalpha")
    E_col = find_column(df, "Ealpha")
    T_col = find_column(df, "Texp_s")
    Br_col = find_column(df, "Br_alpha")

    if args.energy_mode == "Qalpha":
        if Q_col is None:
            raise ValueError(f"--energy-mode Qalpha selected, but no Qalpha column found. Accepted aliases: {ALIASES['Qalpha']}")
        energy_source = "Qalpha"
    elif args.energy_mode == "Ealpha":
        if E_col is None:
            raise ValueError(f"--energy-mode Ealpha selected, but no Ealpha column found. Accepted aliases: {ALIASES['Ealpha']}")
        energy_source = "Ealpha"
    else:
        if Q_col is not None:
            energy_source = "Qalpha"
        elif E_col is not None:
            energy_source = "Ealpha"
        else:
            raise ValueError(
                "Need either a Qalpha or Ealpha column.\n"
                f"Qalpha aliases: {ALIASES['Qalpha']}\n"
                f"Ealpha aliases: {ALIASES['Ealpha']}"
            )

    params = PARAM_SETS[args.param_set]
    results = []

    for idx, row in df.iterrows():
        A = int(row[A_col])
        Z = int(row[Z_col])

        if energy_source == "Qalpha":
            raw = float(row[Q_col])
            Qalpha_MeV = normalize_energy_value(Q_col, raw)
            Ealpha_used_MeV = np.nan
        else:
            raw = float(row[E_col])
            Ealpha_used_MeV = normalize_energy_value(E_col, raw)
            Qalpha_MeV = qalpha_from_ealpha(A, Ealpha_used_MeV)

        res = calc_poenaru(A, Z, Qalpha_MeV, params)
        res["energy_source"] = energy_source
        res["Ealpha_used_MeV"] = Ealpha_used_MeV
        res["param_set"] = args.param_set
        for k, v in params.items():
            res[k] = v

        if T_col is not None and not pd.isna(row[T_col]):
            Texp_s = float(row[T_col])
            Br_alpha = 1.0
            if Br_col is not None and not pd.isna(row[Br_col]):
                Br_alpha = float(row[Br_col])

            if not (0.0 < Br_alpha <= 1.0):
                raise ValueError(f"Invalid alpha branch at row {idx}: Br_alpha={Br_alpha}")

            T12a_exp_s = Texp_s / Br_alpha
            HF = T12a_exp_s / res["T12a_Poenaru_s"]

            res["Texp_s_used"] = Texp_s
            res["Br_alpha_used"] = Br_alpha
            res["T12a_exp_s"] = T12a_exp_s
            res["HF_Poenaru"] = HF
            res["log10_HF_Poenaru"] = math.log10(HF) if HF > 0 else np.nan
        else:
            res["Texp_s_used"] = np.nan
            res["Br_alpha_used"] = np.nan
            res["T12a_exp_s"] = np.nan
            res["HF_Poenaru"] = np.nan
            res["log10_HF_Poenaru"] = np.nan

        results.append(res)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(results)], axis=1)

    if args.output is None:
        base, _ = os.path.splitext(args.input)
        args.output = f"{base}_Poenaru_target.csv"

    out.to_csv(args.output, index=False)
    print(f"Done. Output written to: {args.output}")
    print(f"Formula parameter set: {args.param_set}")
    print(f"Energy source: {energy_source}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
