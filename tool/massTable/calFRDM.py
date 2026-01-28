#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#python calFRDM.py -i FRDM_12.dat -o Q@a_FRDM.dat --zmin 104 --zmax 118
"""Compute Q_alpha from FRDM2012 table (mass excess).

Default input:  FRDM_12.dat
Default output: Q@a_FRDM.dat

Output columns: Z  N  A  Qalpha(MeV)
"""

import argparse

ALPHA_MASS_EXCESS_MEV = 2.4249158692  # 4He atomic mass excess (MeV)

def parse_frdm_masses(path: str):
    masses = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if "ε2" in s or s.startswith("-") or s.startswith("N") or s.startswith("Z"):
                continue
            parts = s.split()
            if not parts or not parts[0].lstrip("+-").isdigit():
                continue
            if len(parts) < 15:
                continue
            try:
                Z = int(parts[0])
                A = int(parts[2])
                Mth = float(parts[14])   # theoretical mass excess column in this export
            except Exception:
                continue
            masses[(Z, A)] = Mth
    return masses

def compute_qalpha(masses: dict, zmin: int, zmax: int):
    rows = []
    for (Z, A), Mp in masses.items():
        if Z < zmin or Z > zmax:
            continue
        Md = masses.get((Z - 2, A - 4))
        if Md is None:
            continue
        Q = Mp - Md - ALPHA_MASS_EXCESS_MEV
        N = A - Z
        rows.append((Z, N, A, Q))
    rows.sort(key=lambda x: (x[0], x[1]))
    return rows

def write_out(rows, out_path: str):
    with open(out_path, "w") as f:
        f.write("# Z  N  A  Qalpha(MeV)\n")
        for Z, N, A, Q in rows:
            f.write(f"{Z:3d} {N:3d} {A:4d} {Q:10.6f}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="FRDM_12.dat", help="FRDM table file (your export)")
    ap.add_argument("-o", "--output", default="Q@a_FRDM.dat", help="output file")
    ap.add_argument("--zmin", type=int, default=104)
    ap.add_argument("--zmax", type=int, default=118)
    args = ap.parse_args()

    masses = parse_frdm_masses(args.input)
    rows = compute_qalpha(masses, args.zmin, args.zmax)
    write_out(rows, args.output)

    print(f"[FRDM] parsed masses: {len(masses)}")
    print(f"[FRDM] wrote {len(rows)} Qalpha rows -> {args.output}")

if __name__ == "__main__":
    main()
