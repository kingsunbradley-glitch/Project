#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#python calHFB.py  -i HFB_24.dat  -o Q@a_HFB.dat  --zmin 104 --zmax 118
"""Compute Q_alpha from HFB-24 table (mass excess).

Default input:  HFB_24.dat
Default output: Q@a_HFB.dat

Output columns: Z  N  A  Qalpha(MeV)
"""

import argparse

ALPHA_MASS_EXCESS_MEV = 2.4249158692  # 4He atomic mass excess (MeV)

def parse_hfb_masses(path: str):
    masses = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("-") or s.startswith("Z"):
                continue
            parts = s.split()
            if not parts or not parts[0].lstrip("+-").isdigit():
                continue
            if len(parts) < 10:
                continue
            try:
                Z = int(parts[0])
                A = int(parts[1])
                Mcal = float(parts[9])   # Mcal column
            except Exception:
                continue
            if abs(Mcal - 999.99) < 1e-6:
                continue
            masses[(Z, A)] = Mcal
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
    ap.add_argument("-i", "--input", default="HFB_32.dat", help="HFB-24 table file (your export)")
    ap.add_argument("-o", "--output", default="Q@a_HFB.dat", help="output file")
    ap.add_argument("--zmin", type=int, default=104)
    ap.add_argument("--zmax", type=int, default=118)
    args = ap.parse_args()

    masses = parse_hfb_masses(args.input)
    rows = compute_qalpha(masses, args.zmin, args.zmax)
    write_out(rows, args.output)

    print(f"[HFB-24] parsed masses: {len(masses)}")
    print(f"[HFB-24] wrote {len(rows)} Qalpha rows -> {args.output}")

if __name__ == "__main__":
    main()
