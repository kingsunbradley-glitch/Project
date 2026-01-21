import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from matplotlib.colors import Normalize, LogNorm, PowerNorm

# =========================
# 0) 全局可调参数（你主要改这里）
# =========================

CMAP_MODE = "viridis"
NORM_MODE = "power"
POWER_GAMMA = 0.7
LOG_EPS = 1e-3
ROBUST_PERCENTILES = (2, 98)

ANNOT_ON = True
ANNOT_FONTSIZE = 10
AUTO_ANNOT_COLOR = True
ANNOT_FMT = ".0f"
ANNOT_LUMA_THRESHOLD = 0.55

GRID_LINEWIDTH = 0.5
GRID_LINECOLOR = "0.95"
NAN_COLOR = "white"

BASE_W, PER_COL_W = 6.0, 1.0
BASE_H, PER_ROW_H = 3.0, 0.4

OUTPUT_DIR = "heatmaps_wide"
INPUT_FILE = "processed_result.txt"

# 打印并列最小/最大时最多列多少条（防止刷屏）
MAX_TIE_PRINT = 200  # <-- 需要就调大/调小


# =========================
# 1) 期刊风格参数
# =========================
def set_style():
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    mpl.rcParams["font.size"] = 11
    mpl.rcParams["axes.labelsize"] = 12
    mpl.rcParams["axes.titlesize"] = 14
    mpl.rcParams["xtick.labelsize"] = 11
    mpl.rcParams["ytick.labelsize"] = 11

    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"

    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["savefig.bbox"] = "tight"


# =========================
# 2) colormap + norm 选择（核心美化）
# =========================
def get_cmap_and_norm(data_2d: np.ndarray):
    p_lo, p_hi = ROBUST_PERCENTILES
    vmin = float(np.nanpercentile(data_2d, p_lo))
    vmax = float(np.nanpercentile(data_2d, p_hi))

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = 0.0, 1.0

    if CMAP_MODE in ["mako", "rocket", "vlag", "icefire", "crest", "flare"]:
        cmap = sns.color_palette(CMAP_MODE, as_cmap=True)
    else:
        try:
            cmap = mpl.colormaps[CMAP_MODE]
        except Exception:
            cmap = plt.get_cmap(CMAP_MODE)

    try:
        cmap = cmap.copy()
    except Exception:
        pass
    cmap.set_bad(color=NAN_COLOR)

    if NORM_MODE == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    elif NORM_MODE == "power":
        norm = PowerNorm(gamma=POWER_GAMMA, vmin=max(0.0, vmin), vmax=vmax)
    elif NORM_MODE == "log":
        positive = data_2d[np.isfinite(data_2d) & (data_2d > 0)]
        min_pos = float(np.min(positive)) if positive.size > 0 else LOG_EPS
        vmin_log = max(min_pos, LOG_EPS)
        norm = LogNorm(vmin=vmin_log, vmax=max(vmax, vmin_log * 1.01))
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    return cmap, norm, (vmin, vmax)


# =========================
# 3) 自动黑/白字标注（修复错位版）
# =========================
def apply_auto_annot_color(ax, data, cmap, norm, threshold=0.55):
    n_rows, n_cols = data.shape
    for t in ax.texts:
        x, y = t.get_position()
        col = int(np.floor(x))
        row = int(np.floor(y))

        if row < 0 or row >= n_rows or col < 0 or col >= n_cols:
            continue

        val = data.iat[row, col]
        if pd.isna(val):
            continue

        r, g, b, a = cmap(norm(val))
        luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
        t.set_color("black" if luminance > threshold else "white")


# =========================
# 4) 主函数：绘图
# =========================
def plot_wide_heatmap(input_file=INPUT_FILE):
    set_style()

    try:
        df = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError:
        print(f"找不到文件 {input_file}，请先运行此前的 sort.py 生成数据。")
        return

    # 确保 N 是数值
    df["N"] = pd.to_numeric(df["N"], errors="coerce")

    # 先算全局最小/最大（画完图再统一打印）
    df_valid = df.dropna(subset=["N"]).copy()
    min_val = max_val = None
    min_rows = max_rows = None

    if not df_valid.empty:
        min_val = df_valid["N"].min()
        max_val = df_valid["N"].max()

        # 并列最小/最大行（只保留你要的列）
        min_rows = df_valid.loc[df_valid["N"] == min_val, ["D1", "Q1", "Q2", "N"]]
        max_rows = df_valid.loc[df_valid["N"] == max_val, ["D1", "Q1", "Q2", "N"]]

        # 排序一下输出更可读（可改）
        min_rows = min_rows.sort_values(["D1", "Q1", "Q2"]).reset_index(drop=True)
        max_rows = max_rows.sort_values(["D1", "Q1", "Q2"]).reset_index(drop=True)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    d1_values = sorted(df["D1"].unique())

    for d1 in d1_values:
        subset = df[df["D1"] == d1]
        heatmap_data = subset.pivot(index="Q2", columns="Q1", values="N")
        heatmap_data.sort_index(ascending=False, inplace=True)
        heatmap_data.sort_index(axis=1, ascending=True, inplace=True)

        n_rows, n_cols = heatmap_data.shape
        fig_width = BASE_W + n_cols * PER_COL_W
        fig_height = BASE_H + n_rows * PER_ROW_H

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        vals = heatmap_data.to_numpy(dtype=float)
        cmap, norm, _ = get_cmap_and_norm(vals)

        hm = sns.heatmap(
            heatmap_data,
            ax=ax,
            annot=ANNOT_ON,
            fmt=ANNOT_FMT,
            annot_kws={"size": ANNOT_FONTSIZE, "weight": "normal"},
            cmap=cmap,
            norm=norm,
            square=False,
            linewidths=GRID_LINEWIDTH,
            linecolor=GRID_LINECOLOR,
            cbar_kws={"label": "\n" + r"$\bf{Average\ Counts\ per\ Hour}$", "shrink": 0.75, "pad": 0.02},
            mask=heatmap_data.isnull(),
        )

        ax.set_title(f"YX   :   D1 = {d1:.3f} T·m", pad=14, fontweight="bold")
        ax.set_xlabel("Q1 / (T/m)", labelpad=10)
        ax.set_ylabel("Q2 / (T/m)", labelpad=10)

        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, length=4)
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)

        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(direction="in", length=3)

        if ANNOT_ON and AUTO_ANNOT_COLOR:
            apply_auto_annot_color(
                ax,
                heatmap_data,
                cmap,
                norm,
                threshold=ANNOT_LUMA_THRESHOLD,
            )

        filename = os.path.join(OUTPUT_DIR, f"heatmap_wide_D1_{d1:.3f}.png")
        plt.savefig(filename)
        plt.close(fig)
        print(f"已保存: {filename}")

    # =========================
    # 画完图后：打印全局最小/最大 + 并列全部列出
    # =========================
    print("\n========== Global N summary ==========")
    if min_rows is None or max_rows is None:
        print("N 全是 NaN/无有效数值，无法统计最小/最大。")
        print("======================================\n")
        return

    # --- MIN ---
    print(f"MIN N = {min_val}  (count={len(min_rows)})")
    if len(min_rows) <= MAX_TIE_PRINT:
        for i, r in min_rows.iterrows():
            print(f"  - D1={r['D1']}  Q1={r['Q1']}  Q2={r['Q2']}")
    else:
        print(f"  (并列最小有 {len(min_rows)} 条，超过 MAX_TIE_PRINT={MAX_TIE_PRINT}，只打印前 {MAX_TIE_PRINT} 条)")
        for i, r in min_rows.head(MAX_TIE_PRINT).iterrows():
            print(f"  - D1={r['D1']}  Q1={r['Q1']}  Q2={r['Q2']}")

    # --- MAX ---
    print(f"\nMAX N = {max_val}  (count={len(max_rows)})")
    if len(max_rows) <= MAX_TIE_PRINT:
        for i, r in max_rows.iterrows():
            print(f"  - D1={r['D1']}  Q1={r['Q1']}  Q2={r['Q2']}")
    else:
        print(f"  (并列最大有 {len(max_rows)} 条，超过 MAX_TIE_PRINT={MAX_TIE_PRINT}，只打印前 {MAX_TIE_PRINT} 条)")
        for i, r in max_rows.head(MAX_TIE_PRINT).iterrows():
            print(f"  - D1={r['D1']}  Q1={r['Q1']}  Q2={r['Q2']}")

    print("======================================\n")


if __name__ == "__main__":
    plot_wide_heatmap(INPUT_FILE)
