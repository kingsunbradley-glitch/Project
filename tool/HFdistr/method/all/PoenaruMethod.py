#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate alpha-decay partial half-lives using the Poenaru-Ivascu-Mazilu
semiempirical formula, and calculate hindrance factors.

Input file can be .csv or .xlsx.

Required columns:
    A, Z

Energy input:
    Either Ealpha_MeV or Qalpha_MeV

Experimental half-life input:
    Optional:
        Texp_s      : experimental total half-life in seconds
        Br_alpha    : alpha branching ratio, default = 1.0

Output:
    input filename with suffix _Poenaru_HF.csv
"""

import argparse
import math
import os
import pandas as pd
import numpy as np


# Poenaru-Ivascu-Mazilu 1980 parameters
B1 = 0.988662
B2 = 0.016314
B3 = 0.020433
B4 = 0.027896
B5 = -0.003033
B6 = -0.16820


# Magic-plus-one numbers used in the paper
N_MAGIC_PLUS_ONE = [51, 83, 127, 185, 229]
Z_MAGIC_PLUS_ONE = [51, 83, 115, 121, 173]


def find_reduced_variable(value, magic_list):
    """
    Calculate reduced shell variable.

    For N:
        y = (N - Ni) / (N_{i+1} - Ni),  Ni < N <= N_{i+1}

    For Z:
        z = (Z - Zi) / (Z_{i+1} - Zi),  Zi < Z <= Z_{i+1}
    """
    for i in range(len(magic_list) - 1):
        low = magic_list[i]
        high = magic_list[i + 1]
        if low < value <= high:
            return (value - low) / (high - low), low, high

    raise ValueError(
        f"Value {value} is outside the implemented magic-plus-one range: {magic_list}"
    )


def calc_qalpha_from_ealpha(A, Ealpha):
    """
    Recoil correction:
        Qalpha = Ealpha * A / (A - 4)
    """
    return Ealpha * A / (A - 4)


def calc_poenaru_log10_t(A, Z, Qalpha):
    """
    Calculate log10(T_alpha / s) using the Poenaru-Ivascu-Mazilu formula.
    """
    N = A - Z
    Ad = A - 4
    Zd = Z - 2

    y, Ni, Nip1 = find_reduced_variable(N, N_MAGIC_PLUS_ONE)
    z, Zi, Zip1 = find_reduced_variable(Z, Z_MAGIC_PLUS_ONE)

    # x parameter in the Coulomb action term
    x = 0.4253 * Qalpha * (1.5874 + Ad ** (1.0 / 3.0)) / Zd

    if not (0.0 < x < 1.0):
        raise ValueError(
            f"Invalid x = {x:.6g}. Need 0 < x < 1. "
            f"Check A={A}, Z={Z}, Qalpha={Qalpha}."
        )

    bracket = math.acos(math.sqrt(x)) - math.sqrt(x * (1.0 - x))

    Ks = (
        2.52956
        * Zd
        * math.sqrt(Ad / (A * Qalpha))
        * bracket
    )

    F = (
        B1
        + B2 * y
        + B3 * z
        + B4 * y * y
        + B5 * y * z
        + B6 * z * z
    )

    log10_T = F * Ks / math.log(10.0) - 20.446

    return {
        "N": N,
        "Ad": Ad,
        "Zd": Zd,
        "Qalpha_MeV": Qalpha,
        "x": x,
        "Ks": Ks,
        "y": y,
        "z": z,
        "N_interval": f"({Ni},{Nip1}]",
        "Z_interval": f"({Zi},{Zip1}]",
        "F_yz": F,
        "log10_Tcalc_s": log10_T,
        "Tcalc_s": 10.0 ** log10_T,
    }


def read_table(filename):
    ext = os.path.splitext(filename)[1].lower()

    if ext == ".csv":
        return pd.read_csv(filename)
    elif ext in [".xlsx", ".xls"]:
        return pd.read_excel(filename)
    else:
        raise ValueError("Input file must be .csv, .xlsx, or .xls")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate alpha-decay half-lives and hindrance factors."
    )
    parser.add_argument("input", help="Input table: .csv, .xlsx, or .xls")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV filename. Default: input_Poenaru_HF.csv"
    )

    args = parser.parse_args()

    df = read_table(args.input)

    required = ["A", "Z"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    if "Qalpha_MeV" not in df.columns and "Ealpha_MeV" not in df.columns:
        raise ValueError("Need either Qalpha_MeV or Ealpha_MeV column.")

    results = []

    for idx, row in df.iterrows():
        A = int(row["A"])
        Z = int(row["Z"])

        if "Qalpha_MeV" in df.columns and not pd.isna(row.get("Qalpha_MeV")):
            Qalpha = float(row["Qalpha_MeV"])
        else:
            Ealpha = float(row["Ealpha_MeV"])
            Qalpha = calc_qalpha_from_ealpha(A, Ealpha)

        res = calc_poenaru_log10_t(A, Z, Qalpha)

        # Experimental partial alpha half-life
        # If Br_alpha is missing, assume alpha branch = 1.
        if "Texp_s" in df.columns and not pd.isna(row.get("Texp_s")):
            Texp_s = float(row["Texp_s"])

            if "Br_alpha" in df.columns and not pd.isna(row.get("Br_alpha")):
                Br_alpha = float(row["Br_alpha"])
            else:
                Br_alpha = 1.0

            if Br_alpha <= 0.0 or Br_alpha > 1.0:
                raise ValueError(
                    f"Invalid Br_alpha={Br_alpha} at row {idx}. "
                    "It should be in the range 0 < Br_alpha <= 1."
                )

            Talpha_exp_s = Texp_s / Br_alpha
            HF = Talpha_exp_s / res["Tcalc_s"]

            res["Br_alpha"] = Br_alpha
            res["Talpha_exp_s"] = Talpha_exp_s
            res["HF"] = HF
            res["log10_HF"] = math.log10(HF)
        else:
            res["Br_alpha"] = np.nan
            res["Talpha_exp_s"] = np.nan
            res["HF"] = np.nan
            res["log10_HF"] = np.nan

        results.append(res)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(results)], axis=1)

    if args.output is None:
        base, _ = os.path.splitext(args.input)
        output = "output_Poenaru_HF.csv"
    else:
        output = args.output

    out.to_csv(output, index=False)

    print(f"Done. Output written to: {output}")
    print()
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()