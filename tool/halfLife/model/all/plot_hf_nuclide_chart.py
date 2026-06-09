#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_hf_nuclide_chart.py

根据 input CSV 生成交互式核素图：
- x 轴：中子数 N
- y 轴：质子数 Z / 元素符号
- 颜色：log10(HF)，hover 中同时显示原始 HF
- 下拉菜单：选择不同 half-life calculation method

输入 CSV 至少需要这些列：
    name, method, HF
推荐同时包含：
    Ealpha, Texp, Tcal

运行示例：
    python plot_hf_nuclide_chart.py main.csv -o hf_nuclide_chart.html

如果缺少 plotly：
    pip install pandas plotly
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

try:
    import plotly.graph_objects as go
except ImportError as exc:
    raise SystemExit(
        "缺少 plotly。请先安装：pip install plotly pandas"
    ) from exc


# 元素表：索引即 Z。覆盖到超重核 Og。
ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYMBOL_TO_Z: Dict[str, int] = {sym: z for z, sym in enumerate(ELEMENTS) if sym}


def parse_nuclide_name(name: str) -> Tuple[int, str, int, int]:
    """
    将 '273Ds' / 'Ds273' / '273-Ds' 解析成 A, symbol, Z, N。
    当前主格式是 '246Fm'。
    """
    s = str(name).strip()

    m = re.match(r"^\s*(\d+)\s*[-_ ]?\s*([A-Za-z]{1,3})\s*$", s)
    if not m:
        m = re.match(r"^\s*([A-Za-z]{1,3})\s*[-_ ]?\s*(\d+)\s*$", s)
        if not m:
            raise ValueError(f"无法解析核素名：{name!r}。请使用类似 273Ds 或 Ds273 的格式。")
        sym = m.group(1).capitalize()
        A = int(m.group(2))
    else:
        A = int(m.group(1))
        sym = m.group(2).capitalize()

    if sym not in SYMBOL_TO_Z:
        raise ValueError(f"未知元素符号：{sym!r}，来自核素名 {name!r}")

    Z = SYMBOL_TO_Z[sym]
    N = A - Z
    return A, sym, Z, N


def fmt_float(x, digits: int = 4) -> str:
    """hover 用的数字格式。"""
    try:
        if pd.isna(x):
            return ""
        x = float(x)
    except Exception:
        return str(x)

    if x == 0:
        return "0"
    ax = abs(x)
    if ax >= 1e4 or ax < 1e-3:
        return f"{x:.{digits}e}"
    return f"{x:.{digits}g}"


def build_chart(csv_path: Path, out_html: Path, title: str = "Nuclide chart colored by HF") -> None:
    df = pd.read_csv(csv_path)

    required = {"name", "method", "HF"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV 缺少必要列：{sorted(missing)}。当前列为：{list(df.columns)}")

    # 解析核素名到 A, Z, N。
    parsed = df["name"].apply(parse_nuclide_name)
    df[["A", "symbol", "Z", "N"]] = pd.DataFrame(parsed.tolist(), index=df.index)

    # 转换数值列。
    for col in ["HF", "Ealpha", "Texp", "Tcal"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # HF 必须为正值才能取 log10。
    bad = df["HF"].isna() | (df["HF"] <= 0)
    if bad.any():
        bad_rows = df.loc[bad, ["name", "method", "HF"]].head(10).to_string(index=False)
        raise ValueError("发现 HF <= 0 或无法转换为数字的行，无法绘制 log10(HF)：\n" + bad_rows)

    df["log10HF"] = df["HF"].apply(math.log10)
    df["label"] = df["A"].astype(str) + df["symbol"]

    # 如果同一 method + 核素重复，保留第一条并给出提示。
    duplicated = df.duplicated(subset=["method", "name"], keep=False)
    if duplicated.any():
        dup_info = (
            df.loc[duplicated, ["method", "name"]]
            .value_counts()
            .reset_index(name="count")
            .head(20)
            .to_string(index=False)
        )
        print("[Warning] 发现重复的 method + name 组合；绘图时保留第一条：")
        print(dup_info)
        df = df.drop_duplicates(subset=["method", "name"], keep="first")

    methods = list(dict.fromkeys(df["method"].astype(str).tolist()))
    methods = sorted(methods, key=lambda x: x.lower())

    z_min, z_max = int(df["Z"].min()), int(df["Z"].max())
    n_min, n_max = int(df["N"].min()), int(df["N"].max())

    tickvals = list(range(z_min, z_max + 1))
    ticktext = [
        f"{z} {ELEMENTS[z]}" if 0 <= z < len(ELEMENTS) and ELEMENTS[z] else str(z)
        for z in tickvals
    ]

    global_cmin = float(df["log10HF"].min())
    global_cmax = float(df["log10HF"].max())

    fig = go.Figure()

    for i, method in enumerate(methods):
        d = df[df["method"].astype(str) == method].copy()
        d = d.sort_values(["Z", "N"])

        hover_lines = []
        for _, r in d.iterrows():
            lines = [
                f"<b>{r['label']}</b>",
                f"Method: {method}",
                f"Z = {int(r['Z'])}, N = {int(r['N'])}, A = {int(r['A'])}",
                f"HF = {fmt_float(r['HF'], 5)}",
                f"log10(HF) = {fmt_float(r['log10HF'], 4)}",
            ]
            if "Ealpha" in d.columns:
                lines.append(f"Eα = {fmt_float(r.get('Ealpha'), 5)} MeV")
            if "Texp" in d.columns:
                lines.append(f"Texp = {fmt_float(r.get('Texp'), 5)} s")
            if "Tcal" in d.columns:
                lines.append(f"Tcal = {fmt_float(r.get('Tcal'), 5)} s")
            hover_lines.append("<br>".join(lines))

        fig.add_trace(
            go.Scatter(
                x=d["N"],
                y=d["Z"],
                mode="markers+text",
                text=d["label"],
                textposition="middle center",
                textfont=dict(size=10, color="black"),
                hovertext=hover_lines,
                hoverinfo="text",
                visible=(i == 0),
                marker=dict(
                    symbol="square",
                    size=42,
                    color=d["log10HF"],
                    cmin=global_cmin,
                    cmax=global_cmax,
                    colorscale="Viridis",
                    colorbar=dict(
                        title="log10(HF)",
                        ticks="outside",
                    ),
                    line=dict(color="rgba(0,0,0,0.55)", width=1.0),
                ),
                name=method,
            )
        )

    buttons = []
    for i, method in enumerate(methods):
        visible = [False] * len(methods)
        visible[i] = True
        buttons.append(
            dict(
                label=method,
                method="update",
                args=[
                    {"visible": visible},
                    {
                        "title": (
                            f"{title}<br>"
                            f"<sup>Method: {method}; color = log10(HF), hover shows raw HF</sup>"
                        )
                    },
                ],
            )
        )

    fig.update_layout(
        title=(
            f"{title}<br>"
            f"<sup>Method: {methods[0]}; color = log10(HF), hover shows raw HF</sup>"
            if methods else title
        ),
        width=1200,
        height=850,
        template="plotly_white",
        font=dict(size=14),
        margin=dict(l=80, r=120, t=110, b=80),
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                x=0.01,
                xanchor="left",
                y=1.12,
                yanchor="top",
                showactive=True,
                bgcolor="white",
                bordercolor="gray",
                borderwidth=1,
            )
        ],
        annotations=[
            dict(
                text="Select method:",
                x=0.01,
                y=1.17,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                font=dict(size=13),
            )
        ],
    )

    fig.update_xaxes(
        title="Neutron number N",
        range=[n_min - 1, n_max + 1],
        dtick=2,
        showgrid=True,
        gridwidth=1,
        zeroline=False,
        mirror=True,
        ticks="outside",
        showline=True,
    )
    fig.update_yaxes(
        title="Proton number Z",
        range=[z_min - 1, z_max + 1],
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
        showgrid=True,
        gridwidth=1,
        zeroline=False,
        mirror=True,
        ticks="outside",
        showline=True,
        scaleanchor="x",   # 让核素格子尽量接近正方形
        scaleratio=1,
    )

    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
    print(f"Saved: {out_html}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an interactive nuclide chart colored by HF with method dropdown."
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default="main.csv",
        help="Input CSV path. Default: main.csv",
    )
    parser.add_argument(
        "-o", "--out",
        default="hf_nuclide_chart.html",
        help="Output HTML file. Default: hf_nuclide_chart.html",
    )
    parser.add_argument(
        "--title",
        default="Nuclide chart colored by HF",
        help="Figure title.",
    )
    args = parser.parse_args()

    build_chart(Path(args.csv), Path(args.out), title=args.title)


if __name__ == "__main__":
    main()
