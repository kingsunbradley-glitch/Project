#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bass fusion-barrier calculator (independent, LISE++-style implementation).

Implements PUBLISHED Bass77 / Bass80 proximity potentials and the Woods-Saxon
example shown in O. Tarasov's LISE++ v9.10 fusion tutorial. This is NOT the
proprietary LISE++ Bass cross-section algorithm.

Dependencies: numpy, scipy, matplotlib (matplotlib needed only for plotting).
Example:
    python calc_bass.py --projectile 40 18 --target 238 92 \
        --energies 217 220 225 --model all --lmax 140 --capture
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import warnings

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

E2_MEV_FM = 1.4399764
HBARC_MEV_FM = 197.3269804
U_MEV_C2 = 931.49410242

PERIODIC = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
    "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def nucleus_label(nucleus):
    z = nucleus.Z
    return f"{nucleus.A}{PERIODIC[z-1]}" if 1 <= z <= len(PERIODIC) else f"A{nucleus.A}Z{z}"


MODEL_LABELS = {
    "bass77": "Bass77",
    "bass80": "Bass80",
    "woods-saxon": "Woods-Saxon",
}


@dataclass(frozen=True)
class Nucleus:
    A: int
    Z: int

    def __post_init__(self):
        if self.A < 1 or self.Z < 1 or self.Z > self.A:
            raise ValueError("Each nucleus must satisfy 1 <= Z <= A")


@dataclass(frozen=True)
class Reaction:
    projectile: Nucleus
    target: Nucleus

    @property
    def mu_mev(self) -> float:
        """Reduced mass * c**2, approximating nuclear masses by A * atomic mass unit."""
        a, b = self.projectile.A, self.target.A
        return (a * b / (a + b)) * U_MEV_C2

    def lab_to_cm(self, energy_lab: float) -> float:
        return energy_lab * self.target.A / (self.projectile.A + self.target.A)

    def cm_to_lab(self, energy_cm: float) -> float:
        return energy_cm * (self.projectile.A + self.target.A) / self.target.A


@dataclass(frozen=True)
class Barrier:
    radius_fm: float
    value_mev: float
    pocket_radius_fm: Optional[float]
    pocket_value_mev: Optional[float]


class Potential:
    def __init__(self, reaction: Reaction, model: str = "bass80",
                 ws_depth: float = 105.0, ws_r0: float = 1.12,
                 ws_a: float = 0.75):
        if model not in MODEL_LABELS:
            raise ValueError(f"Unknown model: {model}")
        if min(ws_depth, ws_r0, ws_a) <= 0:
            raise ValueError("Woods-Saxon parameters must be positive")
        self.rx = reaction
        self.model = model
        self.ws_depth = ws_depth
        self.ws_r0 = ws_r0
        self.ws_a = ws_a
        self.R1 = self.nuclear_radius(reaction.projectile.A)
        self.R2 = self.nuclear_radius(reaction.target.A)
        self.Rtouch = self.R1 + self.R2
        self.Rbar = self.R1 * self.R2 / self.Rtouch
        self.coulomb_scale = E2_MEV_FM * reaction.projectile.Z * reaction.target.Z
        self.lprefactor = HBARC_MEV_FM**2 / (2 * reaction.mu_mev)

    def nuclear_radius(self, A: int) -> float:
        x = A ** (1.0 / 3.0)
        if self.model == "bass77":
            return 1.16 * x - 1.39 / x
        if self.model == "bass80":
            Rs = 1.28 * x - 0.76 + 0.8 / x
            return Rs - 0.98 / Rs
        return self.ws_r0 * x  # used as a WS radius contribution

    def _bass_phi_and_derivative(self, r):
        s = np.asarray(r) - self.Rtouch
        if self.model == "bass77":
            a, b, d1, d2 = 0.030, 0.0061, 3.30, 0.65
        else:
            a, b, d1, d2 = 0.033, 0.007, 3.50, 0.65
        e1, e2 = a * np.exp(np.clip(s / d1, -700, 700)), b * np.exp(np.clip(s / d2, -700, 700))
        den = e1 + e2
        return 1.0 / den, (e1 / d1 + e2 / d2) / (den * den)

    def nuclear(self, r):
        r = np.asarray(r, dtype=float)
        if self.model == "woods-saxon":
            x = (r - self.Rtouch) / self.ws_a
            return -self.ws_depth * expit(-x)
        phi, _ = self._bass_phi_and_derivative(r)
        return -self.Rbar * phi

    def coulomb(self, r):
        return self.coulomb_scale / np.asarray(r, dtype=float)

    def centrifugal(self, r, ell=0):
        return self.lprefactor * ell * (ell + 1) / np.asarray(r, dtype=float)**2

    def total(self, r, ell=0):
        return self.nuclear(r) + self.coulomb(r) + self.centrifugal(r, ell)

    def derivative(self, r, ell=0):
        r = np.asarray(r, dtype=float)
        if self.model == "woods-saxon":
            f = expit((r - self.Rtouch) / self.ws_a)
            dvn = self.ws_depth * f * (1 - f) / self.ws_a
        else:
            _, dphi_neg = self._bass_phi_and_derivative(r)
            dvn = self.Rbar * dphi_neg
        return (dvn - self.coulomb_scale / r**2
                - 2 * self.lprefactor * ell * (ell + 1) / r**3)

    def _stationary_points(self, ell=0):
        """Return all (position, 'min'|'max') via derivative sign crossings."""
        start = 0.5  # include inner extrema even when the Bass tail is extrapolated inward
        stop = self.Rtouch + 30.0
        grid = np.linspace(start, stop, 3501)
        derivatives = self.derivative(grid, ell)
        out = []
        for i in range(len(grid) - 1):
            a, b = float(derivatives[i]), float(derivatives[i + 1])
            if a == 0.0 or a * b < 0.0:
                x = (grid[i] if a == 0.0 else
                     brentq(lambda rr: float(self.derivative(rr, ell)),
                            grid[i], grid[i + 1], xtol=1e-12))
                kind = "min" if a < b else "max"
                if not out or abs(x - out[-1][0]) > 1e-6:
                    out.append((float(x), kind))
        return out

    def barrier(self, ell=0) -> Optional[Barrier]:
        """Outermost physical potential-barrier maximum, if present.

        A barrier requires a preceding pocket minimum. When a pocket disappears
        at high angular momentum, returns None rather than inventing a barrier.
        """
        points = self._stationary_points(ell)
        maximums = [r for r, kind in points if kind == "max"]
        if not maximums:
            return None
        r_b = max(maximums)
        inner_minima = [r for r, kind in points if kind == "min" and r < r_b]
        if not inner_minima:
            return None
        r_p = max(inner_minima)
        return Barrier(r_b, float(self.total(r_b, ell)),
                       r_p, float(self.total(r_p, ell)))


def all_barriers(pot: Potential, lmax: int):
    return {ell: pot.barrier(ell) for ell in range(lmax + 1)}


def capture_cross_section_mb(rx: Reaction, barriers: dict,
                             energy_cm: float, hw: float = 3.0,
                             method: str = "hill-wheeler") -> float:
    """Simplified pocket-limited partial-wave capture; NOT a LISE++ rate.

    sigma = pi/k^2 * sum_l (2l+1) T_l, with 1 fm^2 = 10 mb.
    Missing P_CN, fast fission, evaporation, nuclear deformation etc.
    """
    if energy_cm <= 0 or hw <= 0:
        raise ValueError("E_cm and hbar*omega must both be > 0")
    k2 = 2 * rx.mu_mev * energy_cm / HBARC_MEV_FM**2
    summation = 0.0
    for ell, barrier in barriers.items():
        if barrier is None:
            continue
        if method == "hill-wheeler":
            t = float(expit((2 * np.pi / hw) * (energy_cm - barrier.value_mev)))
        elif method == "classical":
            t = 1.0 if energy_cm >= barrier.value_mev else 0.0
        else:
            raise ValueError("method must be hill-wheeler or classical")
        summation += (2 * ell + 1) * t
    return 10.0 * np.pi / k2 * summation


def write_csv(path: Path, headers: list[str], rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def make_plots(output: Path, rx: Reaction, potentials: dict, all_l: dict,
               selected_l: list[int]):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # One figure per plot, not a multipanel canvas.
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    for model, pot in potentials.items():
        ref = pot.barrier(0)
        rmin = max(0.5, pot.Rtouch - 2)
        rmax = pot.Rtouch + 12
        grid = np.linspace(rmin, rmax, 1400)
        ax.plot(grid, pot.total(grid), label=f"{MODEL_LABELS[model]} (l=0)")
        if ref:
            ax.plot(ref.radius_fm, ref.value_mev, "o", ms=4)
    ax.set_xlabel("Center-to-center distance r (fm)")
    ax.set_ylabel("V(r) (MeV)")
    ax.set_title(f"{nucleus_label(rx.projectile)} + {nucleus_label(rx.target)}: entrance potentials")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output / "potentials_l0.png", dpi=200)
    plt.close(fig)

    # Individual potential contributions for the Bass80 curve if present.
    comp_model = "bass80" if "bass80" in potentials else next(iter(potentials))
    comp = potentials[comp_model]
    grid = np.linspace(max(0.5, comp.Rtouch - 2), comp.Rtouch + 12, 1400)
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    ax.plot(grid, comp.coulomb(grid), label="V_C (point Coulomb)")
    ax.plot(grid, comp.nuclear(grid), label="V_N (nuclear)")
    ax.plot(grid, comp.total(grid), label="V_total, l=0")
    b = comp.barrier(0)
    if b:
        ax.plot(b.radius_fm, b.value_mev, "o", ms=6,
                label=f"Barrier {b.value_mev:.3f} MeV at {b.radius_fm:.3f} fm")
    ax.axhline(0, lw=0.7)
    ax.set_xlabel("Center-to-center distance r (fm)")
    ax.set_ylabel("Potential (MeV)")
    ax.set_title(f"{MODEL_LABELS[comp_model]}: Coulomb + nuclear components")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(output / "potential_components.png", dpi=200)
    plt.close(fig)

    chosen_model = "bass80" if "bass80" in potentials else next(iter(potentials))
    pot = potentials[chosen_model]
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    grid = np.linspace(max(0.5, pot.Rtouch - 2), pot.Rtouch + 12, 1400)
    for ell in selected_l:
        ax.plot(grid, pot.total(grid, ell), label=f"l={ell}")
        barrier = pot.barrier(ell)
        if barrier:
            ax.plot(barrier.radius_fm, barrier.value_mev, "o", ms=4)
    ax.set_xlabel("Center-to-center distance r (fm)")
    ax.set_ylabel("V(r,l) (MeV)")
    ax.set_title(f"{MODEL_LABELS[chosen_model]}: angular-momentum dependence")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output / "potentials_vs_l.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    for model, rows in all_l.items():
        lvalues = [ell for ell, b in rows.items() if b is not None]
        values = [rows[ell].value_mev for ell in lvalues]
        ax.plot(lvalues, values, label=MODEL_LABELS[model])
    ax.set_xlabel("Angular momentum l (hbar)")
    ax.set_ylabel("Barrier B_l (MeV)")
    ax.set_title("Outer barrier vs angular momentum (pocket required)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output / "barrier_vs_l.png", dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--projectile", type=int, nargs=2, metavar=("A", "Z"),
                        default=[40, 18], help="Projectile nucleus, default 40 18")
    parser.add_argument("--target", type=int, nargs=2, metavar=("A", "Z"),
                        default=[238, 92], help="Target nucleus, default 238 92")
    parser.add_argument("--energies", type=float, nargs="*", default=[217, 220, 225],
                        help="Laboratory beam energies in MeV (total, not MeV/u)")
    parser.add_argument("--model", choices=["all", *MODEL_LABELS.keys()], default="all")
    parser.add_argument("--lmax", type=int, default=140)
    parser.add_argument("--show-l", type=int, nargs="*", default=[0, 20, 40, 60])
    parser.add_argument("--ws-depth", type=float, default=105.0)
    parser.add_argument("--ws-r0", type=float, default=1.12)
    parser.add_argument("--ws-a", type=float, default=0.75)
    parser.add_argument("--capture", action="store_true",
                        help="Calculate an illustrative pocket-limited capture estimate")
    parser.add_argument("--hw", type=float, default=3.0,
                        help="hbar omega (MeV) for approximate Hill-Wheeler transmission")
    parser.add_argument("--out-dir", type=Path, default=Path("bass_results"))
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if args.lmax < 0 or any(x < 0 for x in args.show_l):
        parser.error("lmax and all listed l must be nonnegative")
    if any(e <= 0 for e in args.energies):
        parser.error("All beam energies must be positive")
    if args.hw <= 0:
        parser.error("--hw must be positive")

    rx = Reaction(Nucleus(*args.projectile), Nucleus(*args.target))
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    names = list(MODEL_LABELS) if args.model == "all" else [args.model]
    potentials = {name: Potential(rx, name, args.ws_depth, args.ws_r0, args.ws_a)
                  for name in names}
    all_l = {name: all_barriers(pot, args.lmax) for name, pot in potentials.items()}

    print(f"Reaction: {nucleus_label(rx.projectile)} + {nucleus_label(rx.target)}")
    print("Assumptions: spherical nuclei, point Coulomb, nonrelativistic Ecm, "
          "mu≈(A1*A2/(A1+A2))*u; models are independently implemented.")
    print("WARNING: Bass potential pocket and l_crit can rely on deep overlap extrapolation; "
          "these are model-dependent, NOT validated LISE++ outputs.")
    print("\nModel          R_touch(fm)   R_B(fm)    V_B_cm(MeV)  V_B_lab(MeV)  l_crit(discrete)")
    summary_rows = []
    for name, pot in potentials.items():
        b = all_l[name][0]
        good = [ell for ell, val in all_l[name].items() if val is not None]
        lcrit = max(good) + 1 if good and max(good) < args.lmax else None
        summary_rows.append([name, pot.R1, pot.R2, pot.Rtouch,
                             b.radius_fm if b else "", b.value_mev if b else "",
                             rx.cm_to_lab(b.value_mev) if b else "",
                             lcrit if lcrit is not None else "not reached"])
        if b:
            lc = str(lcrit) if lcrit is not None else f">{args.lmax}"
            print(f"{MODEL_LABELS[name]:<14s} {pot.Rtouch:>9.3f}  "
                  f"{b.radius_fm:>9.3f}    {b.value_mev:>11.3f}  "
                  f"{rx.cm_to_lab(b.value_mev):>12.3f}  {lc:>6s}")
        else:
            print(f"{MODEL_LABELS[name]}: no outer barrier with pocket found")
    write_csv(out / "summary.csv",
              ["model", "R1_fm", "R2_fm", "Rtouch_fm", "RB_fm", "VB_cm_MeV", "VB_lab_MeV", "lcrit_discrete"],
              summary_rows)
    if len(all_l) > 0:
        write_csv(out / "barrier_vs_l.csv",
                  ["model", "l", "RB_fm", "VB_cm_MeV", "pocket_R_fm", "pocket_V_MeV"],
                  ([name, ell,
                    b.radius_fm if b else "", b.value_mev if b else "",
                    b.pocket_radius_fm if b else "", b.pocket_value_mev if b else ""]
                   for name, bs in all_l.items() for ell, b in bs.items()))
    rmin = max(0.4, min(p.Rtouch for p in potentials.values()) - 2)
    rmax = max(p.Rtouch for p in potentials.values()) + 12
    grid = np.linspace(rmin, rmax, 1201)
    write_csv(out / "potential_vs_r.csv",
              ["model", "r_fm", "V_C_MeV", "V_N_MeV", "V_l0_MeV"],
              ([name, float(r), float(p.coulomb(r)), float(p.nuclear(r)), float(p.total(r))]
               for name, p in potentials.items() for r in grid))
    print("\nBeam energies (total laboratory kinetic energies):")
    energy_rows = []
    for elab in args.energies:
        ecm = rx.lab_to_cm(elab)
        d = {name: (ecm - all_l[name][0].value_mev if all_l[name][0] else np.nan)
             for name in potentials}
        energy_rows.extend([[name, elab, ecm, d[name]] for name in potentials])
        print(f"  E_lab={elab:>8.3f} MeV ; E_cm={ecm:>8.3f} MeV ; "
              + ", ".join(f"E-VB[{MODEL_LABELS[name]}]={d[name]:+.3f}" for name in potentials))
    write_csv(out / "beam_energies.csv", ["model", "E_lab_MeV", "E_cm_MeV", "E_cm_minus_VB_MeV"], energy_rows)

    if args.capture:
        if args.lmax < 50:
            warnings.warn("--lmax is small; capture cross section may be truncated")
        print("\nIllustrative pocket-limited CAPTURE estimates, NOT LISE++ compound/ER cross sections:")
        capture_rows = []
        for name, bs in all_l.items():
            for elab in args.energies:
                ecm = rx.lab_to_cm(elab)
                classic = capture_cross_section_mb(rx, bs, ecm, args.hw, "classical")
                hw = capture_cross_section_mb(rx, bs, ecm, args.hw, "hill-wheeler")
                capture_rows.append([name, elab, ecm, args.hw, classic, hw])
                print(f"  {MODEL_LABELS[name]:<13s}, E_lab={elab:6.1f}: "
                      f"classical={classic:10.4g} mb ; Hill-Wheeler={hw:10.4g} mb")
        write_csv(out / "capture_estimate.csv",
                  ["model", "E_lab_MeV", "E_cm_MeV", "hbar_omega_MeV", "capture_classical_mb", "capture_HillWheeler_mb"],
                  capture_rows)
    if not args.no_plots:
        make_plots(out, rx, potentials, all_l, args.show_l)
    print(f"\nResults written to {out.resolve()}")


if __name__ == "__main__":
    main()
