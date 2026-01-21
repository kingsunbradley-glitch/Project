import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, PowerNorm

# =========================
# Tunables
# =========================
CMAP_MODE = "viridis"
NORM_MODE = "power"          # linear / power / log
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

INPUT_FILE = "sortedData.txt"
OUTPUT_DIR = "heatmaps_trimmed"

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

def get_cmap_and_norm(data_2d: np.ndarray):
    p_lo, p_hi = ROBUST_PERCENTILES
    finite = data_2d[np.isfinite(data_2d)]

    if finite.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin = float(np.nanpercentile(finite, p_lo))
        vmax = float(np.nanpercentile(finite, p_hi))
        if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or vmin == vmax:
            vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
            if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or vmin == vmax:
                vmin, vmax = 0.0, 1.0

    cmap = mpl.colormaps.get(CMAP_MODE, mpl.colormaps["viridis"])
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
        positive = finite[finite > 0]
        min_pos = float(np.min(positive)) if positive.size > 0 else LOG_EPS
        vmin_log = max(min_pos, LOG_EPS)
        norm = LogNorm(vmin=vmin_log, vmax=max(vmax, vmin_log * 1.01))
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    return cmap, norm

def apply_auto_annot_color(text_objs, cmap, norm, threshold=0.55):
    for t in text_objs:
        val = getattr(t, "_heatmap_val", None)
        if val is None or (not np.isfinite(val)):
            continue
        r, g, b, a = cmap(norm(val))
        luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
        t.set_color("black" if luminance > threshold else "white")

def main():
    set_style()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_FILE, sep="\t")
    for c in ["D1 (T.m)", "Q1 (T/m)", "Q2 (T/m)", "A80eq_fixN_2"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # invalid values (<=0) -> NaN
    df.loc[df["A80eq_fixN_2"] <= 0, "A80eq_fixN_2"] = np.nan

    for d1 in sorted(df["D1 (T.m)"].dropna().unique()):
        subset = df[df["D1 (T.m)"] == d1]
        heatmap_data = subset.pivot(index="Q2 (T/m)", columns="Q1 (T/m)", values="A80eq_fixN_2")

        # KEY: trim columns/rows that are entirely NaN (removes large blank margins)
        heatmap_data = heatmap_data.dropna(axis=1, how="all").dropna(axis=0, how="all")
        if heatmap_data.empty:
            continue

        heatmap_data = heatmap_data.sort_index(ascending=False)
        heatmap_data = heatmap_data.sort_index(axis=1, ascending=True)

        vals = heatmap_data.to_numpy(dtype=float)
        masked = np.ma.masked_invalid(vals)
        cmap, norm = get_cmap_and_norm(vals)

        n_rows, n_cols = masked.shape
        fig_width = BASE_W + n_cols * PER_COL_W
        fig_height = BASE_H + n_rows * PER_ROW_H
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")

        ax.set_xticks(np.arange(n_cols))
        ax.set_yticks(np.arange(n_rows))
        ax.set_xticklabels([f"{x:g}" for x in heatmap_data.columns.to_list()])
        ax.set_yticklabels([f"{y:g}" for y in heatmap_data.index.to_list()])
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, length=4)

        ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
        ax.grid(which="minor", color=GRID_LINECOLOR, linewidth=GRID_LINEWIDTH)
        ax.tick_params(which="minor", bottom=False, left=False)

        ax.set_title(f"YX   :   D1 = {d1:.3f} T·m", pad=14, fontweight="bold")
        ax.set_xlabel("Q1 / (T/m)", labelpad=10)
        ax.set_ylabel("Q2 / (T/m)", labelpad=10)

        cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label("A80eq_fixN_2 (avg)")
        cbar.ax.tick_params(direction="in", length=3)

        if ANNOT_ON:
            text_objs = []
            fmt = "{:" + ANNOT_FMT + "}"
            for i in range(n_rows):
                for j in range(n_cols):
                    v = vals[i, j]
                    if not np.isfinite(v):
                        continue
                    t = ax.text(j, i, fmt.format(v), ha="center", va="center",
                                fontsize=ANNOT_FONTSIZE, weight="normal")
                    t._heatmap_val = float(v)
                    text_objs.append(t)

            if AUTO_ANNOT_COLOR:
                apply_auto_annot_color(text_objs, cmap, norm, threshold=ANNOT_LUMA_THRESHOLD)

        out_png = os.path.join(OUTPUT_DIR, f"heatmap_trimmed_D1_{d1:.3f}.png")
        fig.savefig(out_png)
        plt.close(fig)
        print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
