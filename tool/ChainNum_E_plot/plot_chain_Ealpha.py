#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

FILE = Path("./273DsDataSum.txt")

# ---------- 参数 ----------
SCALE_1000 = True          # 能量 *1000 -> keV
X_NUM_MIN = 8000
X_LEFT_MAX = 9500          # 左段上限
X_RIGHT_MIN = 10500        # 右段下限
X_NUM_MAX = 11500
XTICK_STEP = 500

X_SF = 7900                # SF 放最左侧（独立坐标）
# -------------------------

NAMES = ["ChainNum", "DsE", "HsE", "SgE", "RfE"]

def detect_delim(first_line: str) -> str:
    if "\t" in first_line:
        return "\t"
    if "," in first_line:
        return ","
    return None

def read_table(fname: Path) -> pd.DataFrame:
    if not fname.exists():
        raise FileNotFoundError(f"找不到文件: {fname.resolve()}")

    with fname.open("r", encoding="utf-8", errors="ignore") as f:
        first = ""
        for line in f:
            s = line.rstrip("\n")
            if s.strip() and not s.lstrip().startswith("#"):
                first = s
                break
    if not first:
        raise ValueError("文件为空或只有注释。")

    delim = detect_delim(first)
    if delim is None:
        raise ValueError(
            "检测不到明确分隔符（tab 或逗号）。\n"
            "如果是用空格“对齐”的文本，空字段会丢失导致列左移；请从 Excel 重新导出为 TSV/CSV。"
        )

    has_header = any(ch.isalpha() for ch in first)

    df = pd.read_csv(
        fname,
        sep=delim,
        engine="python",
        header=0 if has_header else None,
        names=None if has_header else NAMES,
        comment="#",
        skip_blank_lines=True,
        dtype=str
    )

    df = df.iloc[:, :5].copy()
    df.columns = NAMES

    for c in NAMES:
        df[c] = df[c].astype("string").str.strip()
        df.loc[df[c] == "", c] = pd.NA

    return df

def to_x(series: pd.Series, scale_1000: bool, x_sf: float) -> np.ndarray:
    s = series.astype("string").str.strip()
    is_sf = s.str.upper().eq("SF").fillna(False)

    x = np.full(len(s), np.nan, dtype=float)
    x[is_sf.to_numpy()] = x_sf

    num = pd.to_numeric(s.where(~is_sf), errors="coerce").to_numpy(dtype=float)
    if scale_1000:
        num *= 1000.0
    num[(~np.isfinite(num)) | (num <= 0)] = np.nan

    x[~is_sf.to_numpy()] = num[~is_sf.to_numpy()]
    return x

def add_break_marks(axL, axR, d=0.012):
    # 断轴“//”符号（zorder 高一点，盖在网格线上）
    kwargsL = dict(transform=axL.transAxes, color="k", clip_on=False, linewidth=1.2, zorder=10)
    axL.plot((1-d, 1+d), (-d, +d), **kwargsL)
    axL.plot((1-d, 1+d), (1-d, 1+d), **kwargsL)

    kwargsR = dict(transform=axR.transAxes, color="k", clip_on=False, linewidth=1.2, zorder=10)
    axR.plot((-d, +d), (-d, +d), **kwargsR)
    axR.plot((-d, +d), (1-d, 1+d), **kwargsR)

def draw_continuous_ygrid(fig, axL, axR, y_positions, linestyle="--", linewidth=0.8, alpha=0.7):
    """
    在 figure 坐标系里画横向网格线，使其跨过断轴中间的空隙不断开。
    """
    # 横线从左子图左边界到右子图右边界
    x0 = axL.get_position().x0
    x1 = axR.get_position().x1
    inv = fig.transFigure.inverted()

    for yy in y_positions:
        # 把 (任意x, yy) 从 data 坐标转到 display，再转到 figure 坐标
        y_disp = axL.transData.transform((0, yy))[1]
        y_fig = inv.transform((0, y_disp))[1]

        fig.add_artist(Line2D(
            [x0, x1], [y_fig, y_fig],
            transform=fig.transFigure,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            zorder=0,          # 放到底层当网格
            clip_on=False
        ))

def main():
    df = read_table(FILE)

    y = np.arange(len(df))
    y_labels = df["ChainNum"].astype("string").fillna("").to_list()

    fig, (axL, axR) = plt.subplots(
        1, 2, sharey=True, figsize=(13, 6),
        gridspec_kw={"width_ratios": [3.2, 2.2], "wspace": 0.03}
    )

    # Rf 放最左（图例顺序也这样）
    cols = [
        ("RfE", r"$^{261}\mathrm{Rf}$", "D"),
        ("SgE", r"$^{265}\mathrm{Sg}$", "^"),
        ("HsE", r"$^{269}\mathrm{Hs}$", "s"),
        ("DsE", r"$^{273}\mathrm{Ds}$", "o"),
    ]

    handles, labels = [], []

    for col, lab, mk in cols:
        x = to_x(df[col], SCALE_1000, X_SF)

        mL = np.isfinite(x) & ((x == X_SF) | ((x >= X_NUM_MIN) & (x <= X_LEFT_MAX)))
        mR = np.isfinite(x) & ((x >= X_RIGHT_MIN) & (x <= X_NUM_MAX))

        h = axL.scatter(x[mL], y[mL], s=45, marker=mk, label=lab, zorder=3)
        axR.scatter(x[mR], y[mR], s=45, marker=mk, label=lab, zorder=3)

        handles.append(h)
        labels.append(lab)

    # y 轴刻度：ChainNum 文本
    axL.set_yticks(y)
    axL.set_yticklabels(y_labels)
    axL.set_ylabel("ChainNum")
    axL.invert_yaxis()

    # 左图 x：SF + 8000-9500，并在 SF 与数值区之间画 dashed 分隔线
    axL.set_xlim(X_SF - 80, X_LEFT_MAX)
    xticks_L = [X_SF] + list(np.arange(X_NUM_MIN, X_LEFT_MAX + 1, XTICK_STEP))
    xlabels_L = ["SF"] + [str(int(t)) for t in xticks_L[1:]]
    axL.set_xticks(xticks_L)
    axL.set_xticklabels(xlabels_L)

    x_sep_sf = (X_SF + X_NUM_MIN) / 2.0
    axL.axvline(x_sep_sf, linestyle="--", linewidth=1.0, alpha=0.9, zorder=2)  # dashed

    # 右图 x：10500-11500
    axR.set_xlim(X_RIGHT_MIN, X_NUM_MAX)
    xticks_R = list(np.arange(X_RIGHT_MIN, X_NUM_MAX + 1, XTICK_STEP))
    axR.set_xticks(xticks_R)
    axR.set_xticklabels([str(int(t)) for t in xticks_R])

    # 断轴拼接外观
    axL.spines["right"].set_visible(False)
    axR.spines["left"].set_visible(False)
    axR.tick_params(labelleft=False, left=False)
    add_break_marks(axL, axR)

    # 不用 ax.grid，改用“跨断轴不断开”的 figure-level 横向虚线网格
    # 先做 layout，再画网格线，避免 tight_layout 改位置后不对齐
    fig.supxlabel(r"$E_{\alpha}\ (\mathrm{keV})$")

    fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=4,
        frameon=False
    )

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.canvas.draw()  # 让坐标变换/布局先稳定

    draw_continuous_ygrid(fig, axL, axR, y_positions=y, linestyle="--", linewidth=0.8, alpha=0.7)

    out = "chain_Ealpha_2D_brokenX_SFleft.png"
    fig.savefig(out, dpi=200)
    plt.show()
    print(f"Saved: {out}")

if __name__ == "__main__":
    main()
