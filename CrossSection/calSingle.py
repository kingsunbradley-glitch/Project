#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
python3 calc_trans_eff.py \
  --n_det 20000 \
  --eps 0.55 \
  --I 3 \
  --q 12 \
  --t 10 \
  --d 500 \
  --A 169 \
  --sigma 500
  '''

"""
Calculate separator transmission efficiency for heavy-ion fusion products.

Formula:
    N_beam = I * t / (q * e)

    n_target = d / A * N_A

    N_prod = N_beam * n_target * sigma

    eta_app = N_det / N_prod

    eta_trans = (N_det / eps_det) / N_prod

Units:
    I          : microampere, uA
    t          : minute
    q          : charge state, e.g. 12 for 12+
    d          : target thickness, ug/cm^2
    A          : target mass number, g/mol approximately
    sigma      : microbarn, ub
    N_det      : detected counts
    eps_det    : detection efficiency, e.g. 0.55
"""

import argparse


E_CHARGE = 1.602176634e-19      # C
N_A = 6.02214076e23             # mol^-1


def calculate(
    n_det: float,
    eps_det: float,
    current_uA: float,
    charge_state: float,
    time_min: float,
    target_thickness_ug_cm2: float,
    target_A: float,
    sigma_ub: float,
):
    # Beam current: uA -> A
    current_A = current_uA * 1e-6

    # Time: min -> s
    time_s = time_min * 60.0

    # Target thickness: ug/cm^2 -> g/cm^2
    target_thickness_g_cm2 = target_thickness_ug_cm2 * 1e-6

    # Cross section: microbarn -> cm^2
    # 1 barn = 1e-24 cm^2
    # 1 microbarn = 1e-30 cm^2
    sigma_cm2 = sigma_ub * 1e-30

    # Total beam particles
    n_beam = current_A * time_s / (charge_state * E_CHARGE)

    # Target areal density
    n_target = target_thickness_g_cm2 / target_A * N_A

    # Number of produced nuclei
    n_prod = n_beam * n_target * sigma_cm2

    # Efficiency without detector-efficiency correction
    eta_app = n_det / n_prod

    # Detector-efficiency corrected detected products
    n_dssd = n_det / eps_det

    # Transmission efficiency
    eta_trans = n_dssd / n_prod

    return {
        "time_s": time_s,
        "current_A": current_A,
        "sigma_cm2": sigma_cm2,
        "target_thickness_g_cm2": target_thickness_g_cm2,
        "n_beam": n_beam,
        "n_target": n_target,
        "n_prod": n_prod,
        "n_dssd": n_dssd,
        "eta_app": eta_app,
        "eta_trans": eta_trans,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate separator transmission efficiency."
    )

    parser.add_argument("--n_det", type=float, required=True,
                        help="Detected counts, e.g. 20000")

    parser.add_argument("--eps", type=float, required=True,
                        help="Detection efficiency, e.g. 0.55")

    parser.add_argument("--I", type=float, required=True,
                        help="Beam current in uA, e.g. 3")

    parser.add_argument("--q", type=float, required=True,
                        help="Charge state, e.g. 12 for 12+")

    parser.add_argument("--t", type=float, required=True,
                        help="Measurement time in minutes, e.g. 10")

    parser.add_argument("--d", type=float, required=True,
                        help="Target thickness in ug/cm^2, e.g. 500")

    parser.add_argument("--A", type=float, default=169,
                        help="Target mass number in g/mol, default = 169 for 169Tm")

    parser.add_argument("--sigma", type=float, required=True,
                        help="Reaction cross section in microbarn, e.g. 500")

    args = parser.parse_args()

    result = calculate(
        n_det=args.n_det,
        eps_det=args.eps,
        current_uA=args.I,
        charge_state=args.q,
        time_min=args.t,
        target_thickness_ug_cm2=args.d,
        target_A=args.A,
        sigma_ub=args.sigma,
    )

    print("\n========== Input ==========")
    print(f"N_det              = {args.n_det:.6g}")
    print(f"eps_det            = {args.eps:.6g}")
    print(f"I                  = {args.I:.6g} uA")
    print(f"q                  = {args.q:.6g}+")
    print(f"t                  = {args.t:.6g} min")
    print(f"d                  = {args.d:.6g} ug/cm^2")
    print(f"A_target           = {args.A:.6g} g/mol")
    print(f"sigma              = {args.sigma:.6g} ub")

    print("\n========== Intermediate ==========")
    print(f"time               = {result['time_s']:.6e} s")
    print(f"N_beam             = {result['n_beam']:.6e}")
    print(f"n_target           = {result['n_target']:.6e} atoms/cm^2")
    print(f"sigma              = {result['sigma_cm2']:.6e} cm^2")
    print(f"N_prod             = {result['n_prod']:.6e}")
    print(f"N_DSSD corrected   = {result['n_dssd']:.6e}")

    print("\n========== Result ==========")
    print(f"Apparent efficiency, without detector correction:")
    print(f"eta_app            = {result['eta_app']:.6e}")
    print(f"eta_app            = {result['eta_app'] * 100:.3f} %")

    print(f"\nTransmission efficiency, detector corrected:")
    print(f"eta_trans          = {result['eta_trans']:.6e}")
    print(f"eta_trans          = {result['eta_trans'] * 100:.3f} %")
    print()


if __name__ == "__main__":
    main()