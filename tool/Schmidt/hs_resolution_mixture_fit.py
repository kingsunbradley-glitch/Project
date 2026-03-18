#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
联合分析来自不同实验室的数据，并显式考虑各自的实验分辨率。

功能：
1. 读取两列数据（如 Decay / EVR），自动忽略空白。
2. 将实验室分辨率（FWHM）写入似然函数。
3. 比较：
   - H1: 单峰 + 实验室零点偏移 + 额外展宽
   - H2: 双峰 + 实验室零点偏移 + 额外展宽
4. 输出 MLE、AIC、BIC。
5. 绘图展示原始点、单峰/双峰拟合曲线，以及 H2 下每个点属于两个峰的后验概率。

默认假设：
- 输入能量单位为 MeV。
- 分辨率输入单位为 keV（例如 EVR=50 keV, Decay=20 keV）。
- 以 EVR 为参考实验室，其零点偏移固定为 0；Decay 允许一个相对偏移 delta_decay。

用法：
python hs_resolution_mixture_fit.py /mnt/data/HsData.txt \
    --fwhm-evr-kev 50 --fwhm-decay-kev 20 --show

如果你想保存图：
python hs_resolution_mixture_fit.py /mnt/data/HsData.txt \
    --savefig hs_fit.png
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm


# =========================
# 工具函数
# =========================

def fwhm_kev_to_sigma_mev(fwhm_kev: float) -> float:
    """将 FWHM(keV) 转为 sigma(MeV)。"""
    return (fwhm_kev / 1000.0) / 2.355


def softplus(x: float | np.ndarray) -> float | np.ndarray:
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)


def safe_logsumexp(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """对 log(a)+log(b) 的两项进行稳定求和。这里输入是两个 log-array。"""
    m = np.maximum(a, b)
    return m + np.log(np.exp(a - m) + np.exp(b - m))


@dataclass
class FitResult:
    name: str
    params: dict
    nll: float
    loglike: float
    k: int
    aic: float
    bic: float
    success: bool
    message: str


# =========================
# 数据读取
# =========================

def read_two_column_table(path: str | Path) -> pd.DataFrame:
    """
    读取带表头的 tab 分隔数据。

    这里不用普通的 whitespace 读取，是因为像 "\t\t8.92" 这种行会被错误折叠，
    从而把列错位。我们显式按 tab 处理，并把多余的 tab 折叠进最后一列。

    输出长表：energy, lab
    """
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        header = [h.strip() for h in f.readline().rstrip("\n").split("\t")]
        ncol = len(header)
        long_rows = []

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < ncol:
                parts += [""] * (ncol - len(parts))
            elif len(parts) > ncol:
                # 保留前 ncol-1 列，其余都并入最后一列，适配类似 "\t\t8.92" 的情况
                parts = parts[: ncol - 1] + ["".join(parts[ncol - 1 :])]

            parts = [p.strip() for p in parts]
            for col, token in zip(header, parts):
                if token == "":
                    continue
                try:
                    val = float(token)
                except ValueError:
                    continue
                long_rows.append({"energy": val, "lab": col})

    out = pd.DataFrame(long_rows)
    if out.empty:
        raise ValueError("没有读到有效数据，请检查文件格式。")
    return out


# =========================
# 模型定义
# =========================

def attach_resolution(df: pd.DataFrame, sigma_map: dict[str, float]) -> pd.DataFrame:
    df = df.copy()
    df["sigma_inst"] = df["lab"].map(sigma_map)
    if df["sigma_inst"].isna().any():
        missing = df.loc[df["sigma_inst"].isna(), "lab"].unique()
        raise ValueError(f"这些实验室没有提供分辨率: {missing}")
    # 以 EVR 为参考实验室，偏移固定为 0；其它实验室拟合相对偏移
    df["is_decay"] = (df["lab"].str.lower() == "decay").astype(float)
    return df


def unpack_h1(theta: np.ndarray) -> dict:
    mu = theta[0]
    delta_decay = theta[1]
    tau = softplus(theta[2])
    return {
        "mu": mu,
        "delta_decay": delta_decay,
        "tau": tau,
    }


def unpack_h2(theta: np.ndarray) -> dict:
    center = theta[0]
    delta_decay = theta[1]
    log_sep = theta[2]
    logit_pi = theta[3]
    tau = softplus(theta[4])

    sep = np.exp(log_sep)
    pi = 1.0 / (1.0 + np.exp(-logit_pi))
    mu1 = center - sep / 2.0
    mu2 = center + sep / 2.0
    return {
        "center": center,
        "delta_decay": delta_decay,
        "sep": sep,
        "pi": pi,
        "mu1": mu1,
        "mu2": mu2,
        "tau": tau,
    }


def nll_h1(theta: np.ndarray, x: np.ndarray, sigma_inst: np.ndarray, is_decay: np.ndarray) -> float:
    p = unpack_h1(theta)
    mu_eff = p["mu"] + p["delta_decay"] * is_decay
    sigma = np.sqrt(sigma_inst**2 + p["tau"]**2)
    ll = norm.logpdf(x, loc=mu_eff, scale=sigma)
    return -np.sum(ll)


def nll_h2(theta: np.ndarray, x: np.ndarray, sigma_inst: np.ndarray, is_decay: np.ndarray) -> float:
    p = unpack_h2(theta)
    delta = p["delta_decay"] * is_decay
    sigma = np.sqrt(sigma_inst**2 + p["tau"]**2)

    ll1 = np.log(p["pi"] + 1e-300) + norm.logpdf(x, loc=p["mu1"] + delta, scale=sigma)
    ll2 = np.log(1.0 - p["pi"] + 1e-300) + norm.logpdf(x, loc=p["mu2"] + delta, scale=sigma)
    ll = safe_logsumexp(ll1, ll2)
    return -np.sum(ll)


# =========================
# 拟合器
# =========================

def fit_model(name: str, fun, starts: list[np.ndarray], bounds, args, unpacker, k: int) -> FitResult:
    best = None
    for theta0 in starts:
        res = minimize(fun, theta0, args=args, method="L-BFGS-B", bounds=bounds)
        if best is None or res.fun < best.fun:
            best = res

    params = unpacker(best.x)
    n = len(args[0])
    nll = float(best.fun)
    loglike = -nll
    aic = 2 * k - 2 * loglike
    bic = k * np.log(n) - 2 * loglike
    return FitResult(
        name=name,
        params=params,
        nll=nll,
        loglike=loglike,
        k=k,
        aic=aic,
        bic=bic,
        success=bool(best.success),
        message=str(best.message),
    )


def build_starts_h1(x: np.ndarray) -> list[np.ndarray]:
    xm = np.mean(x)
    xs = np.std(x, ddof=1) if len(x) > 1 else 0.03
    starts = [
        np.array([xm, 0.0, math.log(max(xs, 1e-4))]),
        np.array([np.median(x), 0.0, math.log(0.01)]),
        np.array([xm, 0.005, math.log(0.02)]),
        np.array([xm, -0.005, math.log(0.02)]),
    ]
    return starts


def build_starts_h2(x: np.ndarray) -> list[np.ndarray]:
    q25, q50, q75 = np.quantile(x, [0.25, 0.5, 0.75])
    starts = [
        np.array([q50, 0.0, math.log(max(q75 - q25, 0.005)), 0.0, math.log(0.01)]),
        np.array([q50, 0.0, math.log(0.015), 0.0, math.log(0.01)]),
        np.array([q50, 0.005, math.log(0.020), math.log(0.7 / 0.3), math.log(0.015)]),
        np.array([q50, -0.005, math.log(0.010), math.log(0.3 / 0.7), math.log(0.015)]),
        np.array([q50, 0.0, math.log(0.030), 0.0, math.log(0.020)]),
    ]
    return starts


# =========================
# H2 后验归属概率
# =========================

def posterior_membership_h2(df: pd.DataFrame, fit: FitResult) -> pd.DataFrame:
    p = fit.params
    out = df.copy()
    delta = p["delta_decay"] * out["is_decay"].values
    sigma = np.sqrt(out["sigma_inst"].values**2 + p["tau"]**2)
    x = out["energy"].values

    comp1 = p["pi"] * norm.pdf(x, loc=p["mu1"] + delta, scale=sigma)
    comp2 = (1.0 - p["pi"]) * norm.pdf(x, loc=p["mu2"] + delta, scale=sigma)
    denom = comp1 + comp2 + 1e-300
    out["p_peak1"] = comp1 / denom
    out["p_peak2"] = comp2 / denom
    return out


# =========================
# 绘图
# =========================

def predicted_density(xgrid: np.ndarray, lab: str, fit: FitResult, sigma_map: dict[str, float]) -> np.ndarray:
    sigma_inst = sigma_map[lab]
    is_decay = 1.0 if lab.lower() == "decay" else 0.0

    if fit.name == "H1_single":
        p = fit.params
        mu = p["mu"] + p["delta_decay"] * is_decay
        sigma = np.sqrt(sigma_inst**2 + p["tau"]**2)
        return norm.pdf(xgrid, loc=mu, scale=sigma)

    elif fit.name == "H2_double":
        p = fit.params
        delta = p["delta_decay"] * is_decay
        sigma = np.sqrt(sigma_inst**2 + p["tau"]**2)
        return (
            p["pi"] * norm.pdf(xgrid, loc=p["mu1"] + delta, scale=sigma)
            + (1.0 - p["pi"]) * norm.pdf(xgrid, loc=p["mu2"] + delta, scale=sigma)
        )
    else:
        raise ValueError("未知模型名称")


def plot_results(df: pd.DataFrame,
                 fit_h1: FitResult,
                 fit_h2: FitResult,
                 sigma_map: dict[str, float],
                 savefig: str | None = None,
                 show: bool = False) -> None:
    labs = list(df["lab"].unique())
    xmin = df["energy"].min() - 0.08
    xmax = df["energy"].max() + 0.08
    xgrid = np.linspace(xmin, xmax, 1200)

    fig = plt.figure(figsize=(10.5, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0])

    # -------- 上图：原始点 + 拟合密度 --------
    ax1 = fig.add_subplot(gs[0, 0])
    y_pos = {lab: i for i, lab in enumerate(labs[::-1], start=1)}

    for lab in labs:
        dd = df[df["lab"] == lab].copy()
        y = np.full(len(dd), y_pos[lab], dtype=float)
        jitter = np.linspace(-0.06, 0.06, len(dd)) if len(dd) > 1 else np.array([0.0])
        ax1.scatter(dd["energy"], y + jitter, s=42, alpha=0.9, label=f"{lab} data")

        dens1 = predicted_density(xgrid, lab, fit_h1, sigma_map)
        dens2 = predicted_density(xgrid, lab, fit_h2, sigma_map)

        # 把密度缩放到实验室对应的水平线上方，便于展示
        ax1.plot(xgrid, y_pos[lab] + 0.35 * dens1 / dens1.max(), lw=2.0, ls="--",
                 label=f"{lab} H1" if lab == labs[0] else None)
        ax1.plot(xgrid, y_pos[lab] + 0.35 * dens2 / dens2.max(), lw=2.2,
                 label=f"{lab} H2" if lab == labs[0] else None)

    ax1.set_yticks([y_pos[lab] for lab in labs])
    ax1.set_yticklabels(labs)
    ax1.set_xlabel("Energy (MeV)")
    ax1.set_ylabel("Dataset")
    ax1.set_title("Different resolutions are written into the likelihood, not added by hand on the plot")
    ax1.grid(alpha=0.25)
    ax1.legend(frameon=False, ncol=3, fontsize=10)

    txt = (
        f"H1: logL={fit_h1.loglike:.2f}, AIC={fit_h1.aic:.2f}, BIC={fit_h1.bic:.2f}\n"
        f"H2: logL={fit_h2.loglike:.2f}, AIC={fit_h2.aic:.2f}, BIC={fit_h2.bic:.2f}\n"
        f"ΔBIC = BIC(H2)-BIC(H1) = {fit_h2.bic - fit_h1.bic:.2f}  (negative favors H2)"
    )
    ax1.text(0.02, 0.98, txt, transform=ax1.transAxes, va="top", ha="left",
             fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    # -------- 下图：H2 下每个点属于两个峰的后验概率 --------
    ax2 = fig.add_subplot(gs[1, 0])
    post = posterior_membership_h2(df, fit_h2)
    post = post.sort_values("energy").reset_index(drop=True)
    idx = np.arange(len(post))

    ax2.scatter(post["energy"], post["p_peak1"], s=48, alpha=0.9, label="P(peak 1 | data)")
    ax2.scatter(post["energy"], post["p_peak2"], s=48, alpha=0.9, label="P(peak 2 | data)")
    ambiguous = post[(post["p_peak2"] > 0.1) & (post["p_peak2"] < 0.9)]
    for _, row in ambiguous.iterrows():
        ax2.text(row["energy"] + 0.002, row["p_peak2"] + 0.02, row["lab"], fontsize=9, alpha=0.85)

    ax2.axvline(fit_h2.params["mu1"], ls="--", lw=1.5)
    ax2.axvline(fit_h2.params["mu2"], ls="--", lw=1.5)
    ax2.set_ylim(-0.03, 1.03)
    ax2.set_xlabel("Energy (MeV)")
    ax2.set_ylabel("Posterior membership")
    ax2.set_title("If H2 is adopted, which points actually support each peak?")
    ax2.grid(alpha=0.25)
    ax2.legend(frameon=False)

    if savefig:
        fig.savefig(savefig, dpi=220, bbox_inches="tight")
        print(f"[Saved figure] {savefig}")
    if show:
        plt.show()
    else:
        plt.close(fig)


# =========================
# 主程序
# =========================

def main():
    parser = argparse.ArgumentParser(description="不同实验室分辨率下的单峰/双峰联合拟合")
    parser.add_argument("input", type=str, help="输入数据文件，例如 HsData.txt")
    parser.add_argument("--fwhm-evr-kev", type=float, default=50.0, help="EVR 分辨率 FWHM (keV)")
    parser.add_argument("--fwhm-decay-kev", type=float, default=20.0, help="Decay 分辨率 FWHM (keV)")
    parser.add_argument("--savefig", type=str, default=None, help="保存图像路径")
    parser.add_argument("--show", action="store_true", help="直接显示图")
    args = parser.parse_args()

    df = read_two_column_table(args.input)

    sigma_map = {
        "EVR": fwhm_kev_to_sigma_mev(args.fwhm_evr_kev),
        "Decay": fwhm_kev_to_sigma_mev(args.fwhm_decay_kev),
    }
    df = attach_resolution(df, sigma_map)

    x = df["energy"].values
    sigma_inst = df["sigma_inst"].values
    is_decay = df["is_decay"].values

    starts_h1 = build_starts_h1(x)
    bounds_h1 = [
        (x.min() - 0.2, x.max() + 0.2),  # mu
        (-0.08, 0.08),                   # delta_decay
        (-12.0, math.log(0.15)),         # tau via softplus(param)
    ]
    fit_h1 = fit_model(
        name="H1_single",
        fun=nll_h1,
        starts=starts_h1,
        bounds=bounds_h1,
        args=(x, sigma_inst, is_decay),
        unpacker=unpack_h1,
        k=3,
    )

    starts_h2 = build_starts_h2(x)
    bounds_h2 = [
        (x.min() - 0.2, x.max() + 0.2),  # center
        (-0.08, 0.08),                   # delta_decay
        (math.log(1e-4), math.log(0.25)),# sep
        (-6.0, 6.0),                     # pi logit
        (-12.0, math.log(0.15)),         # tau via softplus(param)
    ]
    fit_h2 = fit_model(
        name="H2_double",
        fun=nll_h2,
        starts=starts_h2,
        bounds=bounds_h2,
        args=(x, sigma_inst, is_decay),
        unpacker=unpack_h2,
        k=5,
    )

    # -------- 打印结果 --------
    print("=" * 72)
    print("Input summary")
    print(df.groupby("lab").agg(N=("energy", "size"), mean=("energy", "mean"), std=("energy", "std")))
    print("\nInstrument sigma (MeV):")
    for k, v in sigma_map.items():
        print(f"  {k:>5s}: {v:.6f} MeV  ({v*1000:.2f} keV sigma)")

    print("\n" + "=" * 72)
    print("H1: single peak + lab offset + extra width")
    print(f"success = {fit_h1.success}, message = {fit_h1.message}")
    for k, v in fit_h1.params.items():
        if isinstance(v, float):
            print(f"  {k:>12s} = {v:.6f}")
    print(f"  logL = {fit_h1.loglike:.4f}")
    print(f"  AIC  = {fit_h1.aic:.4f}")
    print(f"  BIC  = {fit_h1.bic:.4f}")

    print("\n" + "=" * 72)
    print("H2: double peak + lab offset + extra width")
    print(f"success = {fit_h2.success}, message = {fit_h2.message}")
    for k, v in fit_h2.params.items():
        if isinstance(v, float):
            print(f"  {k:>12s} = {v:.6f}")
    print(f"  logL = {fit_h2.loglike:.4f}")
    print(f"  AIC  = {fit_h2.aic:.4f}")
    print(f"  BIC  = {fit_h2.bic:.4f}")

    print("\n" + "=" * 72)
    d_aic = fit_h2.aic - fit_h1.aic
    d_bic = fit_h2.bic - fit_h1.bic
    print(f"Delta AIC = AIC(H2)-AIC(H1) = {d_aic:.4f}  (negative favors H2)")
    print(f"Delta BIC = BIC(H2)-BIC(H1) = {d_bic:.4f}  (negative favors H2)")

    if d_bic < -6:
        print("BIC interpretation: 对 H2（双峰）有较强支持。")
    elif d_bic < -2:
        print("BIC interpretation: 对 H2（双峰）有一定支持，但不算非常强。")
    elif d_bic < 2:
        print("BIC interpretation: H1 与 H2 差别不大，当前数据不足以强判。")
    else:
        print("BIC interpretation: 当前数据更支持 H1（单峰或至少不需要双峰）。")

    post = posterior_membership_h2(df, fit_h2)
    print("\nPosterior membership under H2:")
    print(post[["lab", "energy", "p_peak1", "p_peak2"]].sort_values("energy").to_string(index=False))

    plot_results(df, fit_h1, fit_h2, sigma_map, savefig=args.savefig, show=args.show)


if __name__ == "__main__":
    main()
