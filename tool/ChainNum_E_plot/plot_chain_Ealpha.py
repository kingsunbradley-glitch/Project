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
X_LEFT_MAX = 9600          # 左段上限
X_RIGHT_MIN = 10400        # 右段下限
X_NUM_MAX = 11500          # 右段上限
XTICK_STEP = 500

X_SF = 7900                # SF 放最左侧

# 要添加的辅助线 (MeV)
AUX_LINES_MEV = [9.13, 8.95, 8.84, 8.69, 8.51,8.28]
# -------------------------

NAMES = ["ChainNum", "DsE", "HsE", "SgE", "RfE"]

def detect_delim(first_line: str) -> str:
    if "\t" in first_line: return "\t"
    if "," in first_line: return ","
    return None

def read_table(fname: Path) -> pd.DataFrame:
    if not fname.exists():
        print(f"Warning: {fname} not found. Using dummy data.")
        # Dummy data for testing
        return pd.DataFrame({
            "ChainNum": ["1", "2", "3", "4"],
            "DsE": ["11.1", "10.8", "SF", "11.4"],
            "HsE": ["9.2", "9.3", "9.1", "SF"],
            "SgE": ["8.5", "SF", "8.6", "8.4"],
            "RfE": ["SF", "7.1", "7.2", "7.0"]
        })

    with fname.open("r", encoding="utf-8", errors="ignore") as f:
        first = ""
        for line in f:
            s = line.rstrip("\n")
            if s.strip() and not s.lstrip().startswith("#"):
                first = s
                break
    if not first: raise ValueError("文件为空或只有注释。")

    delim = detect_delim(first)
    if delim is None: raise ValueError("检测不到明确分隔符。")

    has_header = any(ch.isalpha() for ch in first)
    df = pd.read_csv(fname, sep=delim, engine="python", header=0 if has_header else None,
                     names=None if has_header else NAMES, comment="#", skip_blank_lines=True, dtype=str)
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
    if scale_1000: num *= 1000.0
    num[(~np.isfinite(num)) | (num <= 0)] = np.nan
    x[~is_sf.to_numpy()] = num[~is_sf.to_numpy()]
    return x

# ✅ 修正后的断轴符号函数（自动修正斜率）
def add_break_marks(axL, axR, d=0.012):
    posL = axL.get_position()
    posR = axR.get_position()
    width_ratio = posL.width / posR.width
    
    dL = d
    dR = d * width_ratio # 修正右侧宽度

    kwargs = dict(color="k", clip_on=False, linewidth=1.2, zorder=10)
    # 左轴
    axL.plot((1-dL, 1+dL), (-dL, +dL), transform=axL.transAxes, **kwargs)
    axL.plot((1-dL, 1+dL), (1-dL, 1+dL), transform=axL.transAxes, **kwargs)
    # 右轴
    axR.plot((-dR, +dR), (-dL, +dL), transform=axR.transAxes, **kwargs)
    axR.plot((-dR, +dR), (1-dL, 1+dL), transform=axR.transAxes, **kwargs)

def draw_continuous_ygrid(fig, axL, axR, y_positions, linestyle="--", linewidth=0.8, alpha=0.7, zorder=1.1):
    x0 = axL.get_position().x0
    x1 = axR.get_position().x1
    inv = fig.transFigure.inverted()
    for yy in y_positions:
        y_disp = axL.transData.transform((0, yy))[1]
        y_fig = inv.transform((0, y_disp))[1]
        fig.add_artist(Line2D([x0, x1], [y_fig, y_fig], transform=fig.transFigure,
                              linestyle=linestyle, linewidth=linewidth, alpha=alpha,
                              color="0.6", zorder=zorder, clip_on=False))

def main():
    try:
        df = read_table(FILE)
    except Exception as e:
        print(e)
        return

    y = np.arange(len(df))
    y_labels = df["ChainNum"].astype("string").fillna("").to_list()

    fig, (axL, axR) = plt.subplots(
        1, 2, sharey=True, figsize=(13, 6),
        gridspec_kw={"width_ratios": [3.2, 2.2], "wspace": 0.03}
    )

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

    # -----------------------------------------------------------
    # 【新增】绘制辅助线
    # -----------------------------------------------------------
    for mev_val in AUX_LINES_MEV:
        kev_val = mev_val * 1000.0
        # 你的这些值(8510-9130)都在左轴范围内，所以只画在 axL
        if X_NUM_MIN <= kev_val <= X_LEFT_MAX:
            # 使用 axvline 画垂直线
            axL.axvline(x=kev_val, color='tab:blue', linestyle='-.', linewidth=1.0, alpha=0.8, zorder=2.5)
            
            # (可选) 如果你想在顶部标出数值，可以取消下面这行的注释
            # axL.text(kev_val, -1, f"{mev_val}", rotation=90, va='top', ha='center', color='tab:blue', fontsize=8, transform=axL.get_xaxis_transform())

    # Y 轴设置
    axL.set_yticks(y)
    axL.set_yticklabels(y_labels)
    axL.set_ylabel("ChainNum")
    axL.invert_yaxis()

    # 左轴 X 设置
    axL.set_xlim(X_SF - 80, X_LEFT_MAX)
    xticks_L = [X_SF] + list(np.arange(X_NUM_MIN, X_LEFT_MAX + 1, XTICK_STEP))
    xlabels_L = ["SF"] + [str(int(t)) for t in xticks_L[1:]]
    axL.set_xticks(xticks_L)
    axL.set_xticklabels(xlabels_L)
    
    # SF 分隔线
    x_sep_sf = (X_SF + X_NUM_MIN) / 2.0
    axL.axvline(x_sep_sf, linestyle="--", linewidth=1.0, alpha=0.9, zorder=2, color='gray')

    # 右轴 X 设置 (保留之前要求的：视窗10400起，刻度10500起)
    axR.set_xlim(X_RIGHT_MIN, X_NUM_MAX)
    xticks_R = list(np.arange(10500, X_NUM_MAX + 1, XTICK_STEP))
    axR.set_xticks(xticks_R)
    axR.set_xticklabels([str(int(t)) for t in xticks_R])

    # 样式
    axL.spines["right"].set_visible(False)
    axR.spines["left"].set_visible(False)
    axR.tick_params(labelleft=False, left=False)
    
    # 添加修正后的断轴符号
    add_break_marks(axL, axR)

    fig.supxlabel(r"$E_{\alpha}\ (\mathrm{keV})$")
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=4, frameon=False)
    
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.canvas.draw()
    
    draw_continuous_ygrid(fig, axL, axR, y_positions=y)

    out = "chain_Ealpha_2D_withlines.png"
    fig.savefig(out, dpi=200)
    plt.show()
    print(f"Saved: {out}")

if __name__ == "__main__":
    main()