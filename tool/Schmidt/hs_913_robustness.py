#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
稳健性检验：在已知 8.95 + 9.13 两群存在的前提下，检验 9.13 是否需要再裂成两个子峰。

做三件事：
1) 主拟合：M2(low+high) vs M3(low+high1+high2)
2) Jackknife：每次删掉一个点，重拟合，检查 9.13 裂峰结论是否被少数点驱动
3) Parametric bootstrap under M2：如果真实只有 8.95+9.13 两峰，那么像当前这么大的改进会不会经常偶然出现

说明：
- 不同实验室分辨率直接进似然，不是手工调权重。
- EVR 为参考，Decay 允许相对零点偏移 delta_decay。
- 允许一个额外展宽 tau 吸收未建模展宽。
- 默认参数相对保守：允许更大的 delta/tau/split，以避免“边界卡死”造成的假阳性。
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm


def fwhm_kev_to_sigma_mev(fwhm_kev: float) -> float:
    return (fwhm_kev / 1000.0) / 2.355


def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)


def softmax(v: np.ndarray) -> np.ndarray:
    m = np.max(v)
    e = np.exp(v - m)
    return e / np.sum(e)


def logsumexp_rows(a: np.ndarray) -> np.ndarray:
    m = np.max(a, axis=0)
    return m + np.log(np.sum(np.exp(a - m), axis=0))


def read_two_column_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    with open(path, 'r', encoding='utf-8') as f:
        header = [h.strip() for h in f.readline().rstrip('\n').split('\t')]
        ncol = len(header)
        rows = []
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < ncol:
                parts += [''] * (ncol - len(parts))
            elif len(parts) > ncol:
                parts = parts[: ncol - 1] + [''.join(parts[ncol - 1 :])]
            parts = [p.strip() for p in parts]
            for col, token in zip(header, parts):
                if token == '':
                    continue
                try:
                    rows.append({'energy': float(token), 'lab': col})
                except ValueError:
                    pass
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError('没有读到有效数据。')
    return out


def attach_resolution(df: pd.DataFrame, sigma_map: dict[str, float]) -> pd.DataFrame:
    out = df.copy()
    out['sigma_inst'] = out['lab'].map(sigma_map)
    if out['sigma_inst'].isna().any():
        missing = out.loc[out['sigma_inst'].isna(), 'lab'].unique()
        raise ValueError(f'这些实验室没有分辨率: {missing}')
    out['is_decay'] = (out['lab'].str.lower() == 'decay').astype(float)
    return out


def unpack_m2(theta: np.ndarray) -> dict:
    mu_low = theta[0]
    gap = np.exp(theta[1])
    mu_high = mu_low + gap
    w = softmax(np.array([0.0, theta[2]]))
    delta_decay = theta[3]
    tau = softplus(theta[4])
    return {
        'mu_low': mu_low,
        'mu_high': mu_high,
        'w': w,
        'delta_decay': delta_decay,
        'tau': tau,
    }


def unpack_m3(theta: np.ndarray) -> dict:
    mu_low = theta[0]
    gap = np.exp(theta[1])
    high_center = mu_low + gap
    high_split = np.exp(theta[2])
    mu_h1 = high_center - high_split / 2.0
    mu_h2 = high_center + high_split / 2.0
    w = softmax(np.array([0.0, theta[3], theta[4]]))
    delta_decay = theta[5]
    tau = softplus(theta[6])
    return {
        'mu_low': mu_low,
        'high_center': high_center,
        'high_split': high_split,
        'mu_h1': mu_h1,
        'mu_h2': mu_h2,
        'w': w,
        'delta_decay': delta_decay,
        'tau': tau,
    }


def nll_m2(theta, x, sigma_inst, is_decay) -> float:
    p = unpack_m2(theta)
    delta = p['delta_decay'] * is_decay
    sigma = np.sqrt(sigma_inst**2 + p['tau']**2)
    ll0 = np.log(p['w'][0] + 1e-300) + norm.logpdf(x, loc=p['mu_low'] + delta, scale=sigma)
    ll1 = np.log(p['w'][1] + 1e-300) + norm.logpdf(x, loc=p['mu_high'] + delta, scale=sigma)
    return -np.sum(logsumexp_rows(np.vstack([ll0, ll1])))


def nll_m3(theta, x, sigma_inst, is_decay) -> float:
    p = unpack_m3(theta)
    delta = p['delta_decay'] * is_decay
    sigma = np.sqrt(sigma_inst**2 + p['tau']**2)
    ll0 = np.log(p['w'][0] + 1e-300) + norm.logpdf(x, loc=p['mu_low'] + delta, scale=sigma)
    ll1 = np.log(p['w'][1] + 1e-300) + norm.logpdf(x, loc=p['mu_h1'] + delta, scale=sigma)
    ll2 = np.log(p['w'][2] + 1e-300) + norm.logpdf(x, loc=p['mu_h2'] + delta, scale=sigma)
    return -np.sum(logsumexp_rows(np.vstack([ll0, ll1, ll2])))


def fit_best(fun, starts, bounds, args, maxiter=200):
    best = None
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    for s in starts:
        s = np.clip(np.asarray(s, dtype=float), lo, hi)
        res = minimize(fun, s, args=args, method='L-BFGS-B', bounds=bounds,
                       options={'maxiter': maxiter})
        if best is None or res.fun < best.fun:
            best = res
    return best


def build_starts(max_delta, max_tau, max_split, rng, base2=None, base3=None, nrand=4):
    starts_m2 = []
    if base2 is not None:
        starts_m2.append(base2.copy())
    starts_m2 += [
        np.array([8.90, np.log(0.23), 1.0, 0.0, np.log(0.010)]),
        np.array([8.92, np.log(0.21), 1.5, 0.0, np.log(0.005)]),
        np.array([8.88, np.log(0.25), 0.5, 0.0, np.log(0.015)]),
    ]
    for _ in range(nrand):
        if base2 is not None:
            starts_m2.append(base2 + rng.normal(0, [0.003, 0.03, 0.25, 0.0015, 0.12]))
        starts_m2.append(np.array([
            rng.uniform(8.84, 8.95),
            rng.uniform(np.log(0.12), np.log(0.35)),
            rng.uniform(-2.0, 2.0),
            rng.uniform(-max_delta, max_delta),
            rng.uniform(np.log(0.001), np.log(max_tau)),
        ]))

    starts_m3 = []
    if base3 is not None:
        starts_m3.append(base3.copy())
    if base2 is not None and base3 is None:
        starts_m3.append(np.array([base2[0], base2[1], np.log(0.01), 0.9, 0.6, base2[3], base2[4]]))
    starts_m3 += [
        np.array([8.90, np.log(0.23), np.log(0.015), 0.8, 0.8, 0.0, np.log(0.010)]),
        np.array([8.92, np.log(0.21), np.log(0.020), 1.2, 0.5, 0.0, np.log(0.005)]),
        np.array([8.88, np.log(0.25), np.log(0.030), 0.5, 1.0, 0.0, np.log(0.015)]),
    ]
    for _ in range(nrand):
        if base3 is not None:
            starts_m3.append(base3 + rng.normal(0, [0.003, 0.03, 0.08, 0.3, 0.3, 0.0015, 0.12]))
        starts_m3.append(np.array([
            rng.uniform(8.84, 8.95),
            rng.uniform(np.log(0.12), np.log(0.35)),
            rng.uniform(np.log(0.002), np.log(max_split)),
            rng.uniform(-2.0, 2.0),
            rng.uniform(-2.0, 2.0),
            rng.uniform(-max_delta, max_delta),
            rng.uniform(np.log(0.001), np.log(max_tau)),
        ]))
    return starts_m2, starts_m3


def fit_models(df, max_delta_kev=30.0, max_extra_kev=40.0, max_split_kev=100.0,
               seed=0, nrand=4, base2=None, base3=None, maxiter=200):
    rng = np.random.default_rng(seed)
    x = df['energy'].to_numpy()
    sigma_inst = df['sigma_inst'].to_numpy()
    is_decay = df['is_decay'].to_numpy()

    max_delta = max_delta_kev / 1000.0
    max_tau = max_extra_kev / 1000.0
    max_split = max_split_kev / 1000.0

    bounds_m2 = [
        (8.75, 9.02),
        (np.log(0.08), np.log(0.40)),
        (-5.0, 5.0),
        (-max_delta, max_delta),
        (-20.0, np.log(max_tau)),
    ]
    bounds_m3 = [
        (8.75, 9.02),
        (np.log(0.08), np.log(0.40)),
        (np.log(0.002), np.log(max_split)),
        (-5.0, 5.0),
        (-5.0, 5.0),
        (-max_delta, max_delta),
        (-20.0, np.log(max_tau)),
    ]

    starts_m2, starts_m3 = build_starts(max_delta, max_tau, max_split, rng, base2=base2, base3=base3, nrand=nrand)
    res2 = fit_best(nll_m2, starts_m2, bounds_m2, args=(x, sigma_inst, is_decay), maxiter=maxiter)
    res3 = fit_best(nll_m3, starts_m3, bounds_m3, args=(x, sigma_inst, is_decay), maxiter=maxiter)
    return res2, res3


def simulate_from_m2(df_template: pd.DataFrame, p2: dict, rng) -> pd.DataFrame:
    sigma_inst = df_template['sigma_inst'].to_numpy()
    is_decay = df_template['is_decay'].to_numpy()
    u = rng.uniform(size=len(df_template))
    is_high = u > p2['w'][0]
    mu = np.where(is_high, p2['mu_high'], p2['mu_low']) + p2['delta_decay'] * is_decay
    sigma = np.sqrt(sigma_inst**2 + p2['tau']**2)
    x = rng.normal(mu, sigma)
    out = df_template.copy()
    out['energy'] = x
    return out


def component_pdf_m2(xx: np.ndarray, sigma_ref: float, p: dict):
    s = np.sqrt(sigma_ref**2 + p['tau']**2)
    c0 = p['w'][0] * norm.pdf(xx, loc=p['mu_low'], scale=s)
    c1 = p['w'][1] * norm.pdf(xx, loc=p['mu_high'], scale=s)
    return c0 + c1, [c0, c1]


def component_pdf_m3(xx: np.ndarray, sigma_ref: float, p: dict):
    s = np.sqrt(sigma_ref**2 + p['tau']**2)
    c0 = p['w'][0] * norm.pdf(xx, loc=p['mu_low'], scale=s)
    c1 = p['w'][1] * norm.pdf(xx, loc=p['mu_h1'], scale=s)
    c2 = p['w'][2] * norm.pdf(xx, loc=p['mu_h2'], scale=s)
    return c0 + c1 + c2, [c0, c1, c2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('datafile')
    ap.add_argument('--fwhm-evr-kev', type=float, default=50.0)
    ap.add_argument('--fwhm-decay-kev', type=float, default=20.0)
    ap.add_argument('--max-delta-kev', type=float, default=30.0)
    ap.add_argument('--max-extra-kev', type=float, default=40.0)
    ap.add_argument('--max-high-split-kev', type=float, default=100.0)
    ap.add_argument('--bootstrap', type=int, default=60)
    ap.add_argument('--seed', type=int, default=1234)
    ap.add_argument('--savefig', type=str, default=None)
    ap.add_argument('--show', action='store_true')
    args = ap.parse_args()

    sigma_map = {
        'EVR': fwhm_kev_to_sigma_mev(args.fwhm_evr_kev),
        'Decay': fwhm_kev_to_sigma_mev(args.fwhm_decay_kev),
    }
    df = attach_resolution(read_two_column_table(args.datafile), sigma_map)

    # Main fit
    res2, res3 = fit_models(df, args.max_delta_kev, args.max_extra_kev, args.max_high_split_kev,
                            seed=args.seed, nrand=6, maxiter=220)
    p2 = unpack_m2(res2.x)
    p3 = unpack_m3(res3.x)
    ll2 = -float(res2.fun)
    ll3 = -float(res3.fun)
    lr_obs = 2.0 * (ll3 - ll2)

    print('=== Main fit: M2 (8.95 + 9.13) vs M3 (8.95 + 9.13a + 9.13b) ===')
    print(f'logL(M2) = {ll2:.4f}')
    print(f'logL(M3) = {ll3:.4f}')
    print(f'LR = 2*(logL3-logL2) = {lr_obs:.4f}')
    print(f"M2: mu_low={p2['mu_low']:.6f} MeV, mu_high={p2['mu_high']:.6f} MeV, delta_decay={p2['delta_decay']*1000:.2f} keV, tau={p2['tau']*1000:.2f} keV")
    print(f"M3: mu_low={p3['mu_low']:.6f} MeV, mu_h1={p3['mu_h1']:.6f} MeV, mu_h2={p3['mu_h2']:.6f} MeV, split={p3['high_split']*1000:.2f} keV, delta_decay={p3['delta_decay']*1000:.2f} keV, tau={p3['tau']*1000:.2f} keV")
    print(f"weights M2={np.round(p2['w'],4)}")
    print(f"weights M3={np.round(p3['w'],4)}")
    print()

    # Jackknife
    jack = []
    for i in range(len(df)):
        sub = df.drop(index=df.index[i]).reset_index(drop=True)
        r2_i, r3_i = fit_models(sub, args.max_delta_kev, args.max_extra_kev, args.max_high_split_kev,
                                seed=args.seed + i + 1, nrand=2, base2=res2.x, base3=res3.x, maxiter=120)
        p3_i = unpack_m3(r3_i.x)
        jack.append({
            'removed_idx': int(i),
            'removed_lab': df.iloc[i]['lab'],
            'removed_energy': float(df.iloc[i]['energy']),
            'lr': float(2.0 * (-r3_i.fun + r2_i.fun)),
            'split_kev': float(p3_i['high_split'] * 1000.0),
            'delta_kev': float(p3_i['delta_decay'] * 1000.0),
            'tau_kev': float(p3_i['tau'] * 1000.0),
        })
    jack_df = pd.DataFrame(jack).sort_values('lr')
    print('=== Jackknife (drop 1 point each time) ===')
    print(f"LR range after removing one point: {jack_df['lr'].min():.4f} to {jack_df['lr'].max():.4f}")
    print('Most influential removals (smallest LR first):')
    print(jack_df.head(8).to_string(index=False, float_format=lambda v: f'{v:.3f}'))
    print()

    # Parametric bootstrap under M2 null
    rng = np.random.default_rng(args.seed)
    lrs_boot = []
    fail = 0
    for b in range(args.bootstrap):
        df_b = simulate_from_m2(df, p2, rng)
        try:
            rb2, rb3 = fit_models(df_b, args.max_delta_keV if hasattr(args,'max_delta_keV') else args.max_delta_kev,
                                  args.max_extra_keV if hasattr(args,'max_extra_keV') else args.max_extra_kev,
                                  args.max_high_split_keV if hasattr(args,'max_high_split_keV') else args.max_high_split_kev,
                                  seed=args.seed + 1000 + b, nrand=1, base2=res2.x, base3=res3.x, maxiter=100)
            lrs_boot.append(float(2.0 * (-rb3.fun + rb2.fun)))
        except Exception:
            fail += 1
    lrs_boot = np.asarray(lrs_boot, dtype=float)
    p_boot = np.mean(lrs_boot >= lr_obs) if len(lrs_boot) else np.nan
    print('=== Parametric bootstrap under M2 null ===')
    print(f'bootstrap successful = {len(lrs_boot)}, failed = {fail}')
    if len(lrs_boot):
        print(f'Bootstrap p-value ~= P(LR_boot >= LR_obs | M2 true) = {p_boot:.4f}')
        print(f'LR_boot median = {np.median(lrs_boot):.4f}, 95% quantile = {np.quantile(lrs_boot, 0.95):.4f}, max = {np.max(lrs_boot):.4f}')
    print()

    # Plot
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0])

    ax1 = fig.add_subplot(gs[0, :])
    xx = np.linspace(df['energy'].min() - 0.05, df['energy'].max() + 0.05, 1000)
    pdf2, comps2 = component_pdf_m2(xx, sigma_map['EVR'], p2)
    pdf3, comps3 = component_pdf_m3(xx, sigma_map['EVR'], p3)
    labs = sorted(df['lab'].unique())
    y_positions = {lab: 0.03 + i * 0.02 for i, lab in enumerate(labs)}
    markers = {'EVR': 'o', 'Decay': 's'}
    for lab in labs:
        sub = df[df['lab'] == lab]
        y = np.full(len(sub), y_positions[lab])
        xerr = np.full(len(sub), sigma_map[lab] * 2.355 / 2.0)
        ax1.errorbar(sub['energy'], y, xerr=xerr, fmt=markers.get(lab, 'o'), capsize=2,
                     linestyle='none', label=f'{lab} data')
    ax1.plot(xx, pdf2, linewidth=2.2, label='M2 total')
    ax1.plot(xx, pdf3, linewidth=2.2, linestyle='--', label='M3 total')
    for i, c in enumerate(comps3, start=1):
        ax1.plot(xx, c, linestyle=':', linewidth=1.5, label=f'M3 comp{i}')
    ax1.set_title('Main fit: is the 9.13-MeV group itself split?')
    ax1.set_xlabel('Energy (MeV)')
    ax1.set_ylabel('Relative density / rug level')
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=9, ncol=3)

    ax2 = fig.add_subplot(gs[1, 0])
    if len(lrs_boot):
        ax2.hist(lrs_boot, bins=15)
        ax2.axvline(lr_obs, linestyle='--', linewidth=2, label=f'Observed LR={lr_obs:.2f}')
        ax2.set_title(f'Bootstrap under M2 null (p~{p_boot:.3f})')
        ax2.legend(fontsize=9)
    ax2.set_xlabel('LR statistic')
    ax2.set_ylabel('Count')
    ax2.grid(alpha=0.25)

    ax3 = fig.add_subplot(gs[1, 1])
    order = np.arange(len(jack_df))
    colors = np.where(jack_df['removed_lab'].eq('Decay'), 'tab:orange', 'tab:blue')
    ax3.scatter(order, jack_df['lr'], c=colors)
    ax3.axhline(lr_obs, linestyle='--', linewidth=2, label=f'Full-data LR={lr_obs:.2f}')
    ax3.set_title('Jackknife influence on split evidence')
    ax3.set_xlabel('Removed point rank (sorted by LR)')
    ax3.set_ylabel('LR after removing one point')
    ax3.grid(alpha=0.25)
    ax3.legend(fontsize=9)

    fig.tight_layout()
    if args.savefig:
        fig.savefig(args.savefig, dpi=200)
        print(f'Saved figure to: {args.savefig}')
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
