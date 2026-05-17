#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
python singleReduced.py \
  --Am 273 \
  --Zm 110 \
  --Ea 10.857 \
  --T 17.34e-3 \
  --dT-plus 31.55e-3 \
  --dT-minus 6.83e-3 \
  --BR 1.0 \
  --DeltaL 5
  '''

import math
import argparse

# =========================
# Constants
# =========================
HPLK = 4.13566727e-21      # Planck constant, MeV*s
HBARC = 197.32696          # MeV*fm
AMU = 931.494013           # MeV/c^2
E2 = 1.4399644             # e^2, MeV*fm
LN2 = math.log(2.0)


# =========================
# Potentials
# =========================
def v_nuclear(r, A_d):
    """
    Nuclear Igo potential.
    A_d: daughter mass number
    r: fm
    return: MeV
    """
    return -1100.0 * math.exp(-((r - 1.17 * A_d ** (1.0 / 3.0)) / 0.574))


def v_coulomb(r, Z_d):
    """
    Coulomb potential between alpha and daughter.
    Z_alpha = 2.
    return: MeV
    """
    return 2.0 * Z_d * E2 / r


def reduced_mass(A_d):
    """
    alpha-daughter reduced mass.
    return: MeV/c^2
    """
    return 4.0 * A_d / (4.0 + A_d) * AMU


def v_centrifugal(r, A_d, delta_L):
    """
    Centrifugal potential:
        V_L = hbar^2 L(L+1) / (2 mu r^2)
    return: MeV
    """
    mu = reduced_mass(A_d)
    return delta_L * (delta_L + 1.0) * HBARC ** 2 / (2.0 * mu * r ** 2)


def v_total(r, A_d, Z_d, delta_L):
    """
    Total alpha-daughter potential.
    return: MeV
    """
    return (
        v_nuclear(r, A_d)
        + v_coulomb(r, Z_d)
        + v_centrifugal(r, A_d, delta_L)
    )


# =========================
# Q value correction
# =========================
def electron_screening_energy(Z_m):
    """
    Electron screening correction used in the ROOT macro.
    return: MeV
    """
    return (65.3 * Z_m ** (7.0 / 5.0) - 80.0 * Z_m ** (2.0 / 5.0)) * 1.0e-6


def effective_q_value(A_m, Z_m, E_alpha):
    """
    Effective total alpha-decay energy:
        Q_eff = E_alpha * A_m / A_d + E_screening
    E_alpha: measured alpha energy in MeV
    return: MeV
    """
    A_d = A_m - 4.0
    return E_alpha * (A_m / A_d) + electron_screening_energy(Z_m)


# =========================
# Numerical tools
# =========================
def bisection_root(func, a, b, tol=1.0e-10, max_iter=200):
    fa = func(a)
    fb = func(b)

    if abs(fa) < tol:
        return a
    if abs(fb) < tol:
        return b
    if fa * fb > 0:
        raise RuntimeError("Root is not bracketed.")

    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = func(c)

        if abs(fc) < tol or abs(b - a) < tol:
            return c

        if fa * fc <= 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    return 0.5 * (a + b)


def analytic_outer_turning_point(A_d, Z_d, delta_L, Q_eff):
    """
    Approximate outer turning point using Coulomb + centrifugal potential.
    This is used only to define a safe search range.
    """
    mu = reduced_mass(A_d)
    K = HBARC ** 2 * delta_L * (delta_L + 1.0) / (2.0 * mu)

    return (
        Z_d * E2 / Q_eff
        + math.sqrt((Z_d * E2) ** 2 + Q_eff * K) / Q_eff
    )


def find_turning_points(A_d, Z_d, delta_L, Q_eff, n_scan=30000):
    """
    Find inner and outer classical turning points by scanning V(r)-Q.
    return: Rin, Rout in fm
    """
    def f(r):
        return v_total(r, A_d, Z_d, delta_L) - Q_eff

    r_min = 0.05
    r_outer_guess = analytic_outer_turning_point(A_d, Z_d, delta_L, Q_eff)
    r_max = max(r_outer_guess + 5.0, 50.0)

    roots = []

    x0 = r_min
    y0 = f(x0)

    for i in range(1, n_scan + 1):
        x1 = r_min + (r_max - r_min) * i / n_scan
        y1 = f(x1)

        if y0 == 0.0:
            root = x0
            if not roots or abs(root - roots[-1]) > 1.0e-4:
                roots.append(root)

        elif y0 * y1 < 0.0:
            root = bisection_root(f, x0, x1)
            if not roots or abs(root - roots[-1]) > 1.0e-4:
                roots.append(root)

        x0, y0 = x1, y1

    if len(roots) < 2:
        raise RuntimeError(
            f"Could not find two turning points. Found roots = {roots}"
        )

    return roots[0], roots[-1]


def simpson_integral(func, a, b, n=20000):
    """
    Simpson integration.
    n must be even.
    """
    if n % 2 == 1:
        n += 1

    h = (b - a) / n

    s = func(a) + func(b)

    for i in range(1, n):
        x = a + i * h
        s += 4.0 * func(x) if i % 2 == 1 else 2.0 * func(x)

    return s * h / 3.0


# =========================
# Main physics calculation
# =========================
def calculate_reduced_width(
    A_m,
    Z_m,
    E_alpha,
    T_half,
    BR=1.0,
    delta_L=0.0,
    n_int=20000,
):
    """
    Calculate reduced alpha-decay width delta^2.

    A_m: parent mass number
    Z_m: parent proton number
    E_alpha: alpha energy, MeV
    T_half: half-life, s
    BR: alpha branch ratio, e.g. 0.83
    delta_L: orbital angular momentum carried by alpha
    n_int: integration subdivisions

    return: dictionary
    """
    if T_half <= 0:
        raise ValueError("T_half must be positive.")
    if BR < 0:
        raise ValueError("BR must be non-negative.")
    if E_alpha <= 0:
        raise ValueError("E_alpha must be positive.")
    if delta_L < 0:
        raise ValueError("delta_L must be non-negative.")

    A_d = A_m - 4.0
    Z_d = Z_m - 2.0
    Q_eff = effective_q_value(A_m, Z_m, E_alpha)

    mu = reduced_mass(A_d)

    Rin, Rout = find_turning_points(A_d, Z_d, delta_L, Q_eff)

    def integrand(r):
        val = v_total(r, A_d, Z_d, delta_L) - Q_eff
        return math.sqrt(max(val, 0.0))

    integral = simpson_integral(integrand, Rin, Rout, n=n_int)

    logP = -2.0 * math.sqrt(2.0 * mu) * integral / HBARC
    log10P = logP / math.log(10.0)

    if logP > -745.0:
        P = math.exp(logP)
    else:
        P = 0.0

    lambda_alpha = LN2 / T_half * BR

    # delta^2 = h * lambda / P
    # Use log form for numerical stability.
    if BR == 0:
        delta2_keV = 0.0
    else:
        log_delta2 = math.log(HPLK * LN2 * BR / T_half * 1000.0) - logP
        if log_delta2 < 709.0:
            delta2_keV = math.exp(log_delta2)
        else:
            delta2_keV = float("inf")

    return {
        "A_d": A_d,
        "Z_d": Z_d,
        "Q_eff": Q_eff,
        "Rin": Rin,
        "Rout": Rout,
        "integral": integral,
        "logP": logP,
        "log10P": log10P,
        "P": P,
        "lambda_alpha": lambda_alpha,
        "delta2_keV": delta2_keV,
    }


def calculate_with_errors(
    A_m,
    Z_m,
    E_alpha,
    T_half,
    dT_plus=0.0,
    dT_minus=0.0,
    BR=1.0,
    dBR_plus=0.0,
    dBR_minus=0.0,
    dE_plus=0.0,
    dE_minus=0.0,
    delta_L=0.0,
    n_int=20000,
):
    """
    Asymmetric error estimate.

    Convention:
        T_half = T^{+dT_plus}_{-dT_minus}
        BR     = BR^{+dBR_plus}_{-dBR_minus}
        Ealpha = Ealpha^{+dE_plus}_{-dE_minus}

    Because delta^2 roughly increases when:
        T_half decreases,
        BR increases,
        E_alpha decreases,

    lower delta^2 is evaluated with:
        T_half + dT_plus,
        BR - dBR_minus,
        E_alpha + dE_plus

    upper delta^2 is evaluated with:
        T_half - dT_minus,
        BR + dBR_plus,
        E_alpha - dE_minus
    """
    central = calculate_reduced_width(
        A_m, Z_m, E_alpha, T_half, BR, delta_L, n_int
    )

    # Lower bound of delta^2
    T_for_low = T_half + dT_plus
    BR_for_low = max(BR - dBR_minus, 0.0)
    E_for_low = E_alpha + dE_plus

    low = calculate_reduced_width(
        A_m, Z_m, E_for_low, T_for_low, BR_for_low, delta_L, n_int
    )

    # Upper bound of delta^2
    T_for_high = T_half - dT_minus
    if T_for_high <= 0:
        raise ValueError("T_half - dT_minus must be positive.")

    BR_for_high = BR + dBR_plus
    E_for_high = E_alpha - dE_minus
    if E_for_high <= 0:
        raise ValueError("E_alpha - dE_minus must be positive.")

    high = calculate_reduced_width(
        A_m, Z_m, E_for_high, T_for_high, BR_for_high, delta_L, n_int
    )

    delta_c = central["delta2_keV"]
    delta_low = low["delta2_keV"]
    delta_high = high["delta2_keV"]

    err_minus = delta_c - delta_low
    err_plus = delta_high - delta_c

    return central, err_minus, err_plus, low, high


def main():
    parser = argparse.ArgumentParser(
        description="Calculate reduced alpha-decay width delta^2 with Delta L."
    )

    parser.add_argument("--Am", type=float, required=True, help="Parent mass number, e.g. 273")
    parser.add_argument("--Zm", type=float, required=True, help="Parent proton number, e.g. 110")
    parser.add_argument("--Ea", type=float, required=True, help="Alpha energy in MeV")
    parser.add_argument("--T", type=float, required=True, help="Half-life in seconds")
    parser.add_argument("--BR", type=float, default=1.0, help="Alpha branching ratio")

    parser.add_argument(
        "--DeltaL", "-L",
        type=float,
        default=0.0,
        help="Orbital angular momentum carried by alpha"
    )

    parser.add_argument("--dT-plus", type=float, default=0.0, help="Upper error of half-life in seconds")
    parser.add_argument("--dT-minus", type=float, default=0.0, help="Lower error of half-life in seconds")

    parser.add_argument("--dBR-plus", type=float, default=0.0, help="Upper error of branching ratio")
    parser.add_argument("--dBR-minus", type=float, default=0.0, help="Lower error of branching ratio")

    parser.add_argument("--dE-plus", type=float, default=0.0, help="Upper error of alpha energy in MeV")
    parser.add_argument("--dE-minus", type=float, default=0.0, help="Lower error of alpha energy in MeV")

    parser.add_argument("--nint", type=int, default=20000, help="Number of Simpson integration bins")

    args = parser.parse_args()

    central, err_minus, err_plus, low, high = calculate_with_errors(
        A_m=args.Am,
        Z_m=args.Zm,
        E_alpha=args.Ea,
        T_half=args.T,
        dT_plus=args.dT_plus,
        dT_minus=args.dT_minus,
        BR=args.BR,
        dBR_plus=args.dBR_plus,
        dBR_minus=args.dBR_minus,
        dE_plus=args.dE_plus,
        dE_minus=args.dE_minus,
        delta_L=args.DeltaL,
        n_int=args.nint,
    )

    print("\n================ Reduced alpha-decay width ================")
    print(f"Parent nucleus: A = {args.Am:g}, Z = {args.Zm:g}")
    print(f"Daughter nucleus: A = {central['A_d']:.0f}, Z = {central['Z_d']:.0f}")
    print(f"E_alpha = {args.Ea:.6g} MeV")
    print(f"T_1/2 = {args.T:.6g} s")
    print(f"BR = {args.BR:.6g}")
    print(f"Delta L = {args.DeltaL:g}")

    print("\n---------------- Barrier information ----------------")
    print(f"Q_eff = {central['Q_eff']:.9f} MeV")
    print(f"Rin = {central['Rin']:.6f} fm")
    print(f"Rout = {central['Rout']:.6f} fm")
    print(f"WKB integral = {central['integral']:.9f}")
    print(f"log10(P) = {central['log10P']:.9f}")
    print(f"P = {central['P']:.9e}")

    print("\n---------------- Decay width ----------------")
    print(f"lambda_alpha = {central['lambda_alpha']:.9e} s^-1")
    print(
        f"delta^2 = {central['delta2_keV']:.9g} "
        f"- {err_minus:.9g} + {err_plus:.9g} keV"
    )
    print("============================================================\n")


if __name__ == "__main__":
    main()