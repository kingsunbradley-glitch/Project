#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate No and Fm SF partial half-life tables and plot log10(Tsf/s) vs neutron number N.

Input:
  1. Fm branch table:
     /mnt/data/Fm_SF_branch_T12.dat
     columns: A  SF_branch_percent  T12_value  T12_unit

  2. No data are transcribed in this script from the screenshot.

Output:
  No_SF_branch_T12.dat
  No_SF_Tsf_s.dat
  Fm_No_SF_Tsf_s.dat
  Fmsys.dat
  Fm_No_all_LogTsf_s.{pdf,png,tiff}
  Fm_No_evenN_LogTsf_s.{pdf,png,tiff}
"""

from pathlib import Path
import re
import math
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# Paths
# =========================
script_dir = Path(__file__).resolve().parent
out_dir = script_dir
fm_branch_path = script_dir / "Fm_SF_branch_T12.dat"

no_branch_path = out_dir / "No_SF_branch_T12.dat"
no_tsf_path = out_dir / "No_SF_Tsf_s.dat"
combined_path = out_dir / "Fm_No_SF_Tsf_s.dat"
fmsys_path = out_dir / "Fmsys.dat"
combined_all_stem = out_dir / "Fm_No_all_LogTsf_s"
combined_even_stem = out_dir / "Fm_No_evenN_LogTsf_s"
figure_formats = ("pdf", "png", "tiff")


# =========================
# Basic parsers
# =========================
unit_to_s = {
    "ns": 1e-9,
    "us": 1e-6,
    "μs": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
    "min": 60.0,
    "h": 3600.0,
    "d": 86400.0,
}


def apply_origin_publication_style():
    """
    Origin-like publication style: boxed axes, inward ticks, clean white
    background, print-safe fonts, and vector-friendly PDF text.
    """
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.5,
        "axes.labelsize": 9.5,
        "axes.linewidth": 1.15,
        "axes.edgecolor": "black",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
        "xtick.minor.size": 2.8,
        "ytick.minor.size": 2.8,
        "xtick.major.width": 1.05,
        "ytick.major.width": 1.05,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.35,
        "lines.markersize": 5.0,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def parse_branch_number(x):
    """
    直接取分支比中的数字部分：
      <1.5   -> 1.5
      <=100  -> 100
      >97    -> 97
      100    -> 100

    含 ? 的值跳过，返回 NaN。
    """
    s = str(x).strip()
    if "?" in s or s == "":
        return float("nan")

    s = s.replace("<=", "")
    s = s.replace(">=", "")
    s = s.replace("<", "")
    s = s.replace(">", "")
    return float(s)


def branch_limit_type(x):
    """
    Return the limit carried by the SF branch value itself.
    """
    s = str(x).strip()
    if s.startswith("<"):
        return "upper"
    if s.startswith(">"):
        return "lower"
    return ""


def tsf_limit_from_branch(sf_limit):
    """
    Tsf = T1/2 / brSF, so a branch upper limit gives a Tsf lower limit,
    and a branch lower limit gives a Tsf upper limit.
    """
    if sf_limit == "upper":
        return "lower"
    if sf_limit == "lower":
        return "upper"
    return ""


def parse_value_unc(x):
    """
    解析半衰期及误差，返回:
      value, err_low, err_high

    支持：
      4.6(2)        -> 4.6 ± 0.2
      2.467(16)     -> 2.467 ± 0.016
      43.0+22-15    -> 43.0 +22 -15
      1.23+12-11    -> 1.23 +12 -11
      ≈5            -> 5, no error
    """
    s = str(x).strip().replace(" ", "")
    s = s.replace("~", "≈")

    if s.startswith("≈"):
        return float(s[1:]), float("nan"), float("nan")

    # 非对称误差，例如 43.0+22-15
    m = re.fullmatch(
        r"([0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?\d+)?)"
        r"\+([0-9]+(?:\.[0-9]+)?)"
        r"-([0-9]+(?:\.[0-9]+)?)",
        s
    )
    if m:
        value = float(m.group(1))
        err_high = float(m.group(2))
        err_low = float(m.group(3))
        return value, err_low, err_high

    # 括号误差，例如 2.467(16)
    m = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?\d+)?)\((\d+)\)", s)
    if m:
        val_str = m.group(1)
        unc_digits = m.group(2)

        value = float(val_str)

        if "." in val_str:
            decimals = len(val_str.split(".")[1].split("E")[0].split("e")[0])
        else:
            decimals = 0

        err = int(unc_digits) * 10 ** (-decimals)
        return value, err, err

    return float(s), float("nan"), float("nan")


def make_tsf_table(rows, Z, element):
    """
    rows: list of (A_raw, SF_branch_percent_raw, T12_value_raw, T12_unit)
    """
    output = []

    for A_raw, sf_raw, t_raw, unit in rows:
        A = int(re.match(r"\d+", str(A_raw)).group(0))
        N = A - Z

        sf = parse_branch_number(sf_raw)
        sf_limit = branch_limit_type(sf_raw)
        tsf_limit = tsf_limit_from_branch(sf_limit)

        t_value, t_err_low, t_err_high = parse_value_unc(t_raw)

        unit = str(unit).strip()
        if unit not in unit_to_s:
            raise ValueError(f"Unknown unit: {unit}")

        factor = unit_to_s[unit]

        t_s = t_value * factor
        t_err_low_s = t_err_low * factor if not math.isnan(t_err_low) else float("nan")
        t_err_high_s = t_err_high * factor if not math.isnan(t_err_high) else float("nan")

        if not math.isnan(sf) and sf > 0:
            tsf_s = t_s / (sf / 100.0)

            tsf_err_low_s = (
                t_err_low_s / (sf / 100.0)
                if not math.isnan(t_err_low_s)
                else float("nan")
            )
            tsf_err_high_s = (
                t_err_high_s / (sf / 100.0)
                if not math.isnan(t_err_high_s)
                else float("nan")
            )

            log_tsf = math.log10(tsf_s)

            if not math.isnan(tsf_err_low_s) and tsf_s - tsf_err_low_s > 0:
                log_err_low = log_tsf - math.log10(tsf_s - tsf_err_low_s)
            else:
                log_err_low = float("nan")

            if not math.isnan(tsf_err_high_s):
                log_err_high = math.log10(tsf_s + tsf_err_high_s) - log_tsf
            else:
                log_err_high = float("nan")

        else:
            tsf_s = float("nan")
            tsf_err_low_s = float("nan")
            tsf_err_high_s = float("nan")
            log_tsf = float("nan")
            log_err_low = float("nan")
            log_err_high = float("nan")

        output.append({
            "Element": element,
            "Z": Z,
            "A_raw": A_raw,
            "A": A,
            "N": N,
            "SF_branch_percent_raw": sf_raw,
            "SF_branch_percent_numeric": sf,
            "SF_branch_limit": sf_limit,
            "T12_raw": t_raw,
            "T12_unit": unit,
            "T12_s": t_s,
            "T12_err_low_s": t_err_low_s,
            "T12_err_high_s": t_err_high_s,
            "Tsf_s": tsf_s,
            "Tsf_err_low_s": tsf_err_low_s,
            "Tsf_err_high_s": tsf_err_high_s,
            "Tsf_limit": tsf_limit,
            "Log10_Tsf_s": log_tsf,
            "Log10_Tsf_err_low": log_err_low,
            "Log10_Tsf_err_high": log_err_high,
        })

    return pd.DataFrame(output)


def read_fm_branch_table(path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            A_raw, sf_raw, t_raw, unit = parts[:4]

            # 含 ? 的分支比不计算
            if "?" in sf_raw:
                continue

            try:
                sf = parse_branch_number(sf_raw)
            except Exception:
                continue

            if math.isnan(sf) or sf <= 0:
                continue

            rows.append((A_raw, sf_raw, t_raw, unit))

    return rows


def yerr_arrays(df):
    low = []
    high = []

    for lo, hi in zip(df["Log10_Tsf_err_low"], df["Log10_Tsf_err_high"]):
        low.append(0 if pd.isna(lo) else lo)
        high.append(0 if pd.isna(hi) else hi)

    return [low, high]


def write_fmsys_table(plot_df, path):
    """
    Write the compact plotting interface table.
    Units: T1/2 and T1/2_SF are in seconds; brSF is in percent.
    T1/2_SF_limit marks whether T1/2_SF is a value, upper limit, or lower limit.
    """
    fmsys = plot_df[
        [
            "A_raw",
            "Z",
            "N",
            "T12_s",
            "SF_branch_percent_numeric",
            "Tsf_s",
            "Tsf_limit",
        ]
    ].rename(columns={
        "A_raw": "A",
        "T12_s": "T1/2",
        "SF_branch_percent_numeric": "brSF",
        "Tsf_s": "T1/2_SF",
        "Tsf_limit": "T1/2_SF_limit",
    })
    fmsys["T1/2_SF_limit"] = fmsys["T1/2_SF_limit"].replace("", "value")

    fmsys = fmsys.sort_values(["Z", "N", "A"], kind="stable")
    fmsys.to_csv(path, sep="\t", index=False, float_format="%.8g")


def add_242fm_star(ax):
    """
    242Fm special point:
      Tsf = 3.65 +1.08 -0.68 us
          = 3.65e-6 +1.08e-6 -0.68e-6 s
    """
    T_s = 3.65e-6
    T_low_s = 0.68e-6
    T_high_s = 1.08e-6

    y = math.log10(T_s)
    yerr_low = y - math.log10(T_s - T_low_s)
    yerr_high = math.log10(T_s + T_high_s) - y

    ax.errorbar(
        [142],
        [y],
        yerr=[[yerr_low], [yerr_high]],
        fmt="*",
        markersize=10,
        capsize=3.5,
        capthick=0.9,
        elinewidth=0.9,
        color="red",
        ecolor="red",
        label=r"$^{242}\mathrm{Fm}$",
        zorder=5,
    )


def plot_line_segments(ax, df, color, break_mask):
    segment = []

    for idx, row in df.iterrows():
        if bool(break_mask.loc[idx]):
            if len(segment) >= 2:
                seg = pd.DataFrame(segment)
                ax.plot(
                    seg["N"],
                    seg["Log10_Tsf_s"],
                    "-",
                    color=color,
                    linewidth=1.35,
                    zorder=2,
                )
            segment = []
            continue

        segment.append(row)

    if len(segment) >= 2:
        seg = pd.DataFrame(segment)
        ax.plot(
            seg["N"],
            seg["Log10_Tsf_s"],
            "-",
            color=color,
            linewidth=1.35,
            zorder=2,
        )


def add_limit_markers(ax, df, color):
    for _, row in df.iterrows():
        limit = row["Tsf_limit"]
        if limit not in {"lower", "upper"}:
            continue

        x = row["N"]
        y = row["Log10_Tsf_s"]
        dy = 0.55 if limit == "lower" else -0.55

        ax.annotate(
            "",
            xy=(x, y + dy),
            xytext=(x, y),
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "lw": 0.85,
                "mutation_scale": 7.5,
                "shrinkA": 4,
                "shrinkB": 0,
            },
            zorder=4,
        )


def plot_series(ax, df, marker, color, label, break_mask=None):
    if len(df) == 0:
        return

    if break_mask is None:
        break_mask = pd.Series(False, index=df.index)

    plot_line_segments(ax, df, color, break_mask)

    ax.errorbar(
        df["N"],
        df["Log10_Tsf_s"],
        yerr=yerr_arrays(df),
        fmt=marker,
        linestyle="none",
        color=color,
        markerfacecolor=color,
        markeredgecolor=color,
        markeredgewidth=1.0,
        markersize=5.2,
        elinewidth=0.9,
        capsize=3.0,
        capthick=0.9,
        label=label,
        zorder=3,
    )
    add_limit_markers(ax, df, color)


def save_publication_figure(fig, out_stem):
    saved_paths = []
    for suffix in figure_formats:
        path = out_stem.with_suffix(f".{suffix}")
        fig.savefig(path)
        saved_paths.append(path)
    return saved_paths


def draw_combined(plot_df, out_stem):
    apply_origin_publication_style()
    fig, ax = plt.subplots(figsize=(5.0, 3.55), constrained_layout=True)

    fm = plot_df[plot_df["Element"] == "Fm"].sort_values("N")
    no = plot_df[plot_df["Element"] == "No"].sort_values("N")

    fm242_break = fm["A_raw"].astype(str).eq("242")
    plot_series(ax, fm, "o", "#300EDA", "Fm", break_mask=fm242_break)
    plot_series(ax, no, "s", "#E6A90F", "No")

    add_242fm_star(ax)

    ax.axvline(152, linestyle=(0, (4, 2)), linewidth=0.9, color="0.35", zorder=1)
    ax.text(
        152,
        0.3,
        r"$N=152$",
        transform=ax.get_xaxis_transform(),
        va="center",
        ha="center",
        fontsize=12,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
    )

    ax.set_xlabel(r"Neutron number")
    #ax.set_ylabel(r"$\mathrm{\log_{10}(T^{SF}_{1/2}/s)}$")
    ax.set_ylabel(r"$\log_{10}\!\left(T_{1/2}^{\mathrm{SF}}/\mathrm{s}\right)$")
    ax.set_xlim(138, 162)
    ax.set_ylim(-9, 11)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.tick_params(which="both", direction="in", top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.15)
        spine.set_color("black")

    ax.legend(loc="upper left", handlelength=2.3, borderaxespad=0.6)

    saved_paths = save_publication_figure(fig, out_stem)
    plt.close(fig)
    return saved_paths


def main():
    # =========================
    # No data transcribed from screenshot
    # Columns:
    #   A, SF_branch_percent, T12_value, T12_unit
    # =========================
    no_rows = [
        ("250",   "100",       "4.6(2)",       "us"),
        ("250m",  "100",       "43.0+22-15",   "us"),
        ("251",   "<1.4E-3",   "0.8(1)",       "s"),
        ("252",   "31",        "2.467(16)",    "s"),
        ("254",   "0.17",      "51.2(4)",      "s"),
        ("254m1", "0.02",      "264.9(14)",    "ms"),
        ("254m3", "<=0.012",   "184(3)",       "us"),
        ("256",   "0.55",      "2.91(5)",      "s"),
        ("257",   "<1.5",      "24.5(5)",      "s"),
        ("258",   "<=100",     "1.23+12-11",   "ms"),
        ("259",   "<10",       "58(5)",        "min"),
        ("260",   "100",       "106(8)",       "ms"),
        ("262",   "100",       "≈5",           "ms"),
    ]

    # 先生成 No 的 A、SF分支比、半衰期
    no_branch = pd.DataFrame(
        no_rows,
        columns=["A", "SF_branch_percent", "T12_value", "T12_unit"]
    )
    no_branch.to_csv(no_branch_path, sep="\t", index=False)

    # No: 计算 Tsf
    no_tsf = make_tsf_table(no_rows, Z=102, element="No")
    no_tsf.to_csv(no_tsf_path, sep="\t", index=False, float_format="%.8g")

    # Fm: 从已有 dat 读入并计算 Tsf
    if not fm_branch_path.exists():
        raise FileNotFoundError(
            f"Cannot find {fm_branch_path}. "
            "Please put Fm_SF_branch_T12.dat in the same directory as this script."
        )

    fm_rows = read_fm_branch_table(fm_branch_path)
    fm_tsf = make_tsf_table(fm_rows, Z=100, element="Fm")

    # 合并表
    combined = pd.concat([fm_tsf, no_tsf], ignore_index=True)
    combined.to_csv(combined_path, sep="\t", index=False, float_format="%.8g")

    # 所有有数据点，包括 isomer
    all_plot = combined[combined["Log10_Tsf_s"].notna()].copy()
    write_fmsys_table(all_plot, fmsys_path)

    all_figures = draw_combined(
        all_plot,
        combined_all_stem,
    )

    # 偶中子核，只取基态 A_raw 为纯数字
    even_plot = combined[
        combined["A_raw"].astype(str).str.fullmatch(r"\d+")
        & (combined["N"] % 2 == 0)
        & combined["Log10_Tsf_s"].notna()
    ].copy()

    even_figures = draw_combined(
        even_plot,
        combined_even_stem,
    )

    print(f"Saved: {no_branch_path}")
    print(f"Saved: {no_tsf_path}")
    print(f"Saved: {combined_path}")
    print(f"Saved: {fmsys_path}")
    for path in [*all_figures, *even_figures]:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
