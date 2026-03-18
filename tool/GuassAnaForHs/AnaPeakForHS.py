#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
high_group_m1_m2_fit.py

Fit the high-energy alpha group only (default: energy >= 9000 keV) and compare:

    M1: single peak
    M2: double peak

This version is intended for literature values with reported uncertainties.
It does NOT introduce any extra inter-lab offset or extra broadening.

Accepted input formats (whitespace-delimited text with header), e.g.:

    energy_keV sigma_keV source label
    9103 20 Decay paperA
    9230 30 Decay paperB
    9030 50 EVR   paperC

Also supported:
    energy_MeV sigma_MeV ...
    energy sigma ...

Extra columns are ignored.
"""

import argparse
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

SQRT2PI = math.sqrt(2.0 * math.pi)


def norm_pdf(x, mu, sigma):
    sigma = np.asarray(sigma, dtype=float)
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (SQRT2PI * sigma)


def negloglike_m1(theta, x, s):
    mu = theta[0]
    pdf = norm_pdf(x, mu, s)
    pdf = np.clip(pdf, 1e-300, None)
    return -np.sum(np.log(pdf))


def unpack_m2(theta):
    mu1 = float(theta[0])
    gap = math.exp(float(theta[1]))
    mu2 = mu1 + gap
    a = float(theta[2])
    w1 = 1.0 / (1.0 + math.exp(-a))
    w2 = 1.0 - w1
    return mu1, mu2, w1, w2


def negloglike_m2(theta, x, s):
    mu1, mu2, w1, w2 = unpack_m2(theta)
    pdf = w1 * norm_pdf(x, mu1, s) + w2 * norm_pdf(x, mu2, s)
    pdf = np.clip(pdf, 1e-300, None)
    return -np.sum(np.log(pdf))


def fit_m1(x, s):
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    w = 1.0 / (s ** 2)
    mu0 = float(np.sum(w * x) / np.sum(w))
    res = minimize(
        negloglike_m1,
        x0=np.array([mu0]),
        args=(x, s),
        method="L-BFGS-B",
    )
    return {
        "mu": float(res.x[0]),
        "logL": -float(res.fun),
        "success": bool(res.success),
    }


def fit_m2(x, s, n_starts=40, seed=1234):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)

    xmin, xmax = float(np.min(x)), float(np.max(x))
    xr = max(xmax - xmin, 1e-6)
    q25, q50, q75 = np.quantile(x, [0.25, 0.5, 0.75])

    starts = []
    for m1 in [q25, q50 - 0.1 * xr, q25 - 0.05 * xr]:
        for gap in [0.02 * xr, 0.05 * xr, 0.10 * xr, 0.15 * xr]:
            gap = max(gap, 1e-3)
            for w1 in [0.25, 0.4, 0.5, 0.6, 0.75]:
                a = math.log(w1 / (1.0 - w1))
                starts.append([m1, math.log(gap), a])

    while len(starts) < n_starts:
        m1 = rng.uniform(xmin - 0.05 * xr, q75)
        gap = rng.uniform(max(1e-3, 0.01 * xr), max(2e-3, 0.25 * xr))
        w1 = rng.uniform(0.1, 0.9)
        starts.append([m1, math.log(gap), math.log(w1 / (1.0 - w1))])

    best_res = None
    best_fun = None
    for st in starts:
        try:
            res = minimize(
                negloglike_m2,
                x0=np.array(st, dtype=float),
                args=(x, s),
                method="L-BFGS-B",
            )
            if best_fun is None or res.fun < best_fun:
                best_fun = float(res.fun)
                best_res = res
        except Exception:
            continue

    if best_res is None:
        raise RuntimeError("M2 fit failed for all starting points.")

    mu1, mu2, w1, w2 = unpack_m2(best_res.x)
    return {
        "mu1": mu1,
        "mu2": mu2,
        "w1": w1,
        "w2": w2,
        "logL": -float(best_res.fun),
        "success": bool(best_res.success),
    }


def aic(logL, k):
    return 2 * k - 2 * logL


def bic(logL, k, n):
    return k * math.log(n) - 2 * logL


def infer_columns(df):
    cols = {c.lower(): c for c in df.columns}

    energy_col = None
    sigma_col = None
    unit = None

    if "energy_kev" in cols:
        energy_col = cols["energy_kev"]
        unit = "keV"
    elif "energy_mev" in cols:
        energy_col = cols["energy_mev"]
        unit = "MeV"
    elif "energy" in cols:
        energy_col = cols["energy"]
        unit = "auto"

    if "sigma_kev" in cols:
        sigma_col = cols["sigma_kev"]
    elif "sigma_mev" in cols:
        sigma_col = cols["sigma_mev"]
    elif "sigma" in cols:
        sigma_col = cols["sigma"]
    else:
        raise ValueError("Could not find sigma column.")

    if energy_col is None:
        raise ValueError("Could not find energy column.")

    return energy_col, sigma_col, unit


def load_data(path, emin_kev=9000.0):
    df = pd.read_csv(path, sep=r"\s+", comment="#")
    energy_col, sigma_col, unit = infer_columns(df)

    e = pd.to_numeric(df[energy_col], errors="coerce").to_numpy(dtype=float)
    s = pd.to_numeric(df[sigma_col], errors="coerce").to_numpy(dtype=float)

    mask = np.isfinite(e) & np.isfinite(s) & (s > 0)
    df = df.loc[mask].copy()
    e = e[mask]
    s = s[mask]

    if unit == "keV":
        e_kev = e
        s_kev = s
    elif unit == "MeV":
        e_kev = e * 1000.0
        s_kev = s * 1000.0
    else:
        if np.nanmedian(e) < 100.0:
            e_kev = e * 1000.0
            s_kev = s * 1000.0 if np.nanmedian(s) < 10.0 else s
        else:
            e_kev = e
            s_kev = s

    df["energy_kev_internal"] = e_kev
    df["sigma_kev_internal"] = s_kev

    df = df[df["energy_kev_internal"] >= float(emin_kev)].copy()
    if len(df) < 3:
        raise ValueError("Too few points remain after the high-energy selection.")

    return df.reset_index(drop=True)


def simulate_under_m1(mu, s, rng):
    return rng.normal(loc=mu, scale=s, size=len(s))


def bootstrap_pvalue(x, s, fit1, lr_obs, n_boot=200, seed=1234):
    rng = np.random.default_rng(seed)
    lrs = []
    n_fail = 0

    for _ in range(int(n_boot)):
        xb = simulate_under_m1(fit1["mu"], s, rng)
        try:
            f1 = fit_m1(xb, s)
            f2 = fit_m2(xb, s, seed=int(rng.integers(1, 10_000_000)))
            lr = 2.0 * (f2["logL"] - f1["logL"])
            lrs.append(lr)
        except Exception:
            n_fail += 1

    lrs = np.asarray(lrs, dtype=float)
    if len(lrs) == 0:
        return {"pvalue": np.nan, "lrs": lrs, "n_fail": n_fail}

    pvalue = float(np.mean(lrs >= lr_obs))
    return {"pvalue": pvalue, "lrs": lrs, "n_fail": n_fail}


def jackknife_lr(x, s):
    rows = []
    n = len(x)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        xi = x[mask]
        si = s[mask]
        try:
            f1 = fit_m1(xi, si)
            f2 = fit_m2(xi, si, seed=1000 + i)
            lr = 2.0 * (f2["logL"] - f1["logL"])
        except Exception:
            lr = np.nan
        rows.append({"removed_idx": i, "removed_energy_kev": x[i], "lr": lr})
    return pd.DataFrame(rows)


def plot_results(df, fit1, fit2, lr_obs, boot, jack, savefig=None, show=False):
    x = df["energy_kev_internal"].to_numpy()
    s = df["sigma_kev_internal"].to_numpy()

    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.0])

    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])

    labels_done = set()
    for _, row in df.iterrows():
        src = str(row["source"]) if "source" in df.columns else "data"
        lab = src if src not in labels_done else None
        labels_done.add(src)
        ax0.errorbar(
            row["energy_kev_internal"],
            0.0,
            xerr=row["sigma_kev_internal"],
            fmt="o",
            ms=6,
            capsize=3,
            alpha=0.9,
            label=lab,
        )

    xx = np.linspace(float(np.min(x) - 4 * np.max(s)), float(np.max(x) + 4 * np.max(s)), 1200)

    yy1 = norm_pdf(xx, fit1["mu"], np.median(s))
    yy1 = yy1 / np.max(yy1)

    y21 = fit2["w1"] * norm_pdf(xx, fit2["mu1"], np.median(s))
    y22 = fit2["w2"] * norm_pdf(xx, fit2["mu2"], np.median(s))
    yy2 = y21 + y22
    yy2 = yy2 / np.max(yy2)
    y21n = y21 / np.max(yy2)
    y22n = y22 / np.max(yy2)

    ax0.plot(xx, yy1, lw=2, label=f"M1 single peak ({fit1['mu']:.1f} keV)")
    ax0.plot(xx, yy2, "--", lw=2, label=f"M2 double peak ({fit2['mu1']:.1f}, {fit2['mu2']:.1f} keV)")
    ax0.plot(xx, y21n, ":", lw=1.5, label="M2 component 1")
    ax0.plot(xx, y22n, ":", lw=1.5, label="M2 component 2")

    ax0.set_xlabel("Energy (keV)")
    ax0.set_ylabel("Normalized shape / data markers")
    ax0.set_title("High-energy group: M1 (single peak) vs M2 (double peak)")
    ax0.legend(frameon=False, fontsize=9)

    txt = (
        f"logL(M1) = {fit1['logL']:.4f}\n"
        f"logL(M2) = {fit2['logL']:.4f}\n"
        f"LR = 2*(logL2-logL1) = {lr_obs:.4f}\n"
        f"M1: mu = {fit1['mu']:.3f} keV\n"
        f"M2: mu1 = {fit2['mu1']:.3f} keV, mu2 = {fit2['mu2']:.3f} keV\n"
        f"weights = [{fit2['w1']:.3f}, {fit2['w2']:.3f}]"
    )
    ax0.text(
        0.02, 0.98, txt,
        transform=ax0.transAxes, va="top", ha="left",
        fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

    lrs = boot["lrs"]
    if len(lrs) > 0:
        ax1.hist(lrs, bins=min(20, max(8, len(lrs) // 5)), alpha=0.8)
        ax1.axvline(lr_obs, linestyle="--", linewidth=2, label=f"obs LR={lr_obs:.3f}")
        ax1.legend(frameon=False, fontsize=9)
    ax1.set_xlabel("Bootstrap LR under M1 null")
    ax1.set_ylabel("Count")
    ax1.set_title(f"Parametric bootstrap (p ≈ {boot['pvalue']:.3f})")

    if len(jack) > 0:
        jj = jack.sort_values("lr").reset_index(drop=True)
        ax2.plot(np.arange(len(jj)), jj["lr"].to_numpy(), "o-", ms=4)
        ax2.axhline(lr_obs, linestyle="--", linewidth=2, label=f"full-data LR={lr_obs:.3f}")
        ax2.legend(frameon=False, fontsize=9)
    ax2.set_xlabel("Removed-point rank")
    ax2.set_ylabel("LR after removing one point")
    ax2.set_title("Jackknife influence check")

    fig.tight_layout()

    if savefig:
        fig.savefig(savefig, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Fit high-energy literature points: M1 single peak vs M2 double peak.")
    parser.add_argument("datafile", help="Input text file with energy/sigma columns.")
    parser.add_argument("--emin-kev", type=float, default=9000.0, help="Lower cut for high-energy group in keV (default: 9000).")
    parser.add_argument("--bootstrap", type=int, default=200, help="Number of bootstrap replicas (default: 200).")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed.")
    parser.add_argument("--savefig", type=str, default=None, help="Optional output figure path.")
    parser.add_argument("--show", action="store_true", help="Show the plot window.")
    args = parser.parse_args()

    df = load_data(args.datafile, emin_kev=args.emin_kev)
    x = df["energy_kev_internal"].to_numpy(dtype=float)
    s = df["sigma_kev_internal"].to_numpy(dtype=float)

    fit1 = fit_m1(x, s)
    fit2 = fit_m2(x, s, seed=args.seed)

    k1, k2 = 1, 3
    n = len(x)
    lr = 2.0 * (fit2["logL"] - fit1["logL"])

    print("=== High-energy group: M1 (single peak) vs M2 (double peak) ===")
    print(f"N = {n}")
    print(f"logL(M1) = {fit1['logL']:.6f}")
    print(f"logL(M2) = {fit2['logL']:.6f}")
    print(f"LR = 2*(logL2-logL1) = {lr:.6f}")
    print(f"AIC(M1) = {aic(fit1['logL'], k1):.6f}")
    print(f"AIC(M2) = {aic(fit2['logL'], k2):.6f}")
    print(f"BIC(M1) = {bic(fit1['logL'], k1, n):.6f}")
    print(f"BIC(M2) = {bic(fit2['logL'], k2, n):.6f}")
    print(f"M1: mu = {fit1['mu']:.6f} keV")
    print(f"M2: mu1 = {fit2['mu1']:.6f} keV, mu2 = {fit2['mu2']:.6f} keV")
    print(f"weights M2 = [{fit2['w1']:.4f}, {fit2['w2']:.4f}]")

    jack = jackknife_lr(x, s)
    if len(jack) > 0:
        print("\n=== Jackknife (drop 1 point each time) ===")
        print(f"LR range after removing one point: {np.nanmin(jack['lr']):.6f} to {np.nanmax(jack['lr']):.6f}")
        print(jack.sort_values('lr').head(8).to_string(index=False))

    boot = bootstrap_pvalue(x, s, fit1, lr, n_boot=args.bootstrap, seed=args.seed)
    print("\n=== Parametric bootstrap under M1 null ===")
    print(f"bootstrap successful = {len(boot['lrs'])}, failed = {boot['n_fail']}")
    print(f"Bootstrap p-value ~= P(LR_boot >= LR_obs | M1 true) = {boot['pvalue']:.6f}")
    if len(boot["lrs"]) > 0:
        print(f"LR_boot median = {np.median(boot['lrs']):.6f}, 95% quantile = {np.quantile(boot['lrs'], 0.95):.6f}, max = {np.max(boot['lrs']):.6f}")

    plot_results(df, fit1, fit2, lr, boot, jack, savefig=args.savefig, show=args.show)


if __name__ == "__main__":
    main()
