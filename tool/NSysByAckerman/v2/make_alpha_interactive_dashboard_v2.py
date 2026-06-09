#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a transferable ROOT-like interactive HTML dashboard for alpha-systematics CSV files.

This v2 script is compatible with the CSV produced by make_alpha_systematics.py v2,
for example:
    out/alpha_systematics_Z100_104_N148_156.csv

It also keeps backward compatibility with the older Ackermann-style CSV when the
corresponding columns are present.

Main updates vs the earlier dashboard
-------------------------------------
* Draw Q_alpha error bars from unc_Qalpha_keV / Qalpha_err_MeV.
* Draw alpha partial half-life error bars from unc_T12_alpha_s / T_alpha_half_err_s.
* Drop entries carrying estimated/systematics '#' markers by default when the CSV
  still contains has_estimated_marker or raw string fields with '#'.
* Highlight N=152 by default for the Z=100--104, N=148--156 systematics.
* Hide variables whose y-column is absent, so TSF is shown only when present.
* Use the new v2 column names: T12_total_s, alpha_branch, T12_alpha_s,
  alpha_branch_flag, T12_alpha_limit, decay_summary, etc.

Default usage in tool/152NSys:
    python make_alpha_interactive_dashboard_v2.py

Explicit usage:
    python make_alpha_interactive_dashboard_v2.py \
        --csv out/alpha_systematics_Z100_104_N148_156.csv \
        --html out/alpha_interactive_dashboard.html

Open the resulting HTML with a browser, or forward that single HTML file.
The HTML uses Plotly from CDN, so the browser needs internet access to load Plotly.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype


Z_BY_SYMBOL = {
    "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104,
    "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
    "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116,
    "Ts": 117, "Og": 118,
}

ELEMENT_COLORS = {
    "Es": "#6B7280", "Fm": "#111827", "Md": "#8C564B", "No": "#E377C2",
    "Lr": "#7F7F7F", "Rf": "#BCBD22", "Db": "#17BECF", "Sg": "#8C564B",
    "Bh": "#FF7F0E", "Hs": "#9467BD", "Mt": "#2CA02C", "Ds": "#1F77B4",
    "Rg": "#D62728", "Cn": "#4B5563", "Nh": "#0F766E", "Fl": "#9333EA",
    "Mc": "#F97316", "Lv": "#DC2626", "Ts": "#2563EB", "Og": "#16A34A",
}

ELEMENT_SYMBOLS = {
    "Es": "circle", "Fm": "square", "Md": "diamond", "No": "triangle-up",
    "Lr": "triangle-down", "Rf": "cross", "Db": "x", "Sg": "pentagon",
    "Bh": "hexagon", "Hs": "diamond-wide", "Mt": "triangle-left", "Ds": "square",
    "Rg": "circle", "Cn": "triangle-right", "Nh": "star", "Fl": "hourglass",
    "Mc": "bowtie", "Lv": "asterisk", "Ts": "hash", "Og": "circle-open-dot",
}


def first_existing_column(df: pd.DataFrame, names: Iterable[str]) -> Optional[str]:
    for name in names:
        if name in df.columns:
            return name
    return None


def to_numeric_series(df: pd.DataFrame, name: Optional[str]) -> pd.Series:
    if name is None:
        return pd.Series([math.nan] * len(df), index=df.index, dtype="float64")
    return pd.to_numeric(df[name], errors="coerce")


def clean_value(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return float(value)
    return str(value)


def has_hash_marker(value: Any) -> bool:
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except TypeError:
        pass
    return "#" in str(value)


def row_has_estimated_marker(row: pd.Series) -> bool:
    if "has_estimated_marker" in row.index:
        v = row.get("has_estimated_marker")
        if isinstance(v, bool):
            return v
        if str(v).strip().lower() in {"true", "1", "yes", "y"}:
            return True
    critical = [
        "qa", "unc_qa", "Qalpha_keV", "unc_Qalpha_keV", "Qalpha_MeV", "Qalpha_err_MeV",
        "half_life", "half_life_sec", "unc_hls", "T12_total_s", "unc_T12_total_s",
        "decay_1_%", "decay_2_%", "decay_3_%", "unc_1", "unc_2", "unc_3",
        "alpha_branch", "alpha_branch_unc", "alpha_branch_percent",
    ]
    return any(col in row.index and has_hash_marker(row.get(col)) for col in critical)


def format_with_unc(value: Any, unc: Any, unit: str = "", qualifier: Any = "") -> str:
    if value is None:
        return ""
    try:
        vf = float(value)
    except Exception:  # noqa: BLE001
        return ""
    if not math.isfinite(vf):
        return ""
    prefix = ""
    q = "" if qualifier is None else str(qualifier).strip().lower()
    if q in {"<", "upper", "upper limit"}:
        prefix = "< "
    elif q in {">", "lower", "lower limit"}:
        prefix = "> "
    elif q in {"approx", "~", "≈"}:
        prefix = "≈ "
    text = f"{prefix}{vf:.6g}"
    try:
        uf = float(unc)
        if math.isfinite(uf) and uf > 0:
            text += f" ± {uf:.2g}"
    except Exception:  # noqa: BLE001
        pass
    return f"{text} {unit}".strip()


def guess_symbol(z: Any) -> str:
    try:
        zi = int(z)
    except Exception:  # noqa: BLE001
        return ""
    for sym, zz in Z_BY_SYMBOL.items():
        if zz == zi:
            return sym
    return ""


def normalize_dataframe(df: pd.DataFrame, keep_estimated: bool = False) -> Tuple[pd.DataFrame, int]:
    """Normalize v2 alpha-systematics CSV or older Ackermann-style CSV to canonical columns."""
    df = df.copy()
    if not keep_estimated:
        mask_est = df.apply(row_has_estimated_marker, axis=1)
        n_excluded = int(mask_est.sum())
        df = df.loc[~mask_est].copy()
    else:
        n_excluded = 0

    col_z = first_existing_column(df, ["Z", "z"])
    col_n = first_existing_column(df, ["N", "n"])
    col_a = first_existing_column(df, ["A", "a", "mass"])
    col_symbol = first_existing_column(df, ["symbol", "element", "Element"])
    col_nuclide = first_existing_column(df, ["nuclide", "isotope", "isotope_label", "Nuclide"])

    if col_z is None or col_n is None:
        raise ValueError(f"CSV must contain Z/N or z/n columns. Columns present: {list(df.columns)}")

    out = pd.DataFrame(index=df.index)
    out["Z"] = to_numeric_series(df, col_z).astype("Int64")
    out["N"] = to_numeric_series(df, col_n).astype("Int64")
    if col_a is not None:
        out["A"] = to_numeric_series(df, col_a).astype("Int64")
    else:
        out["A"] = (out["Z"].astype("float64") + out["N"].astype("float64")).round().astype("Int64")

    if col_symbol is not None:
        out["symbol"] = df[col_symbol].astype(str).str.strip()
    else:
        out["symbol"] = [guess_symbol(z) for z in out["Z"]]

    if col_nuclide is not None:
        out["nuclide"] = df[col_nuclide].astype(str).str.strip()
    else:
        out["nuclide"] = [f"{'' if pd.isna(a) else int(a)}{sym}" for a, sym in zip(out["A"], out["symbol"])]

    out["state"] = df[first_existing_column(df, ["state", "State"])] if first_existing_column(df, ["state", "State"]) else ""

    # Q_alpha and uncertainty.
    q_mev_col = first_existing_column(df, ["Qalpha_MeV", "Q_alpha_MeV", "Qalpha", "Q_alpha"])
    q_kev_col = first_existing_column(df, ["Qalpha_keV", "Q_alpha_keV", "qa", "qalpha_keV"])
    if q_mev_col:
        out["Qalpha_MeV"] = to_numeric_series(df, q_mev_col)
    elif q_kev_col:
        out["Qalpha_MeV"] = to_numeric_series(df, q_kev_col) / 1000.0
    else:
        out["Qalpha_MeV"] = math.nan

    qerr_mev_col = first_existing_column(df, ["Qalpha_err_MeV", "Q_alpha_err_MeV", "unc_Qalpha_MeV"])
    qerr_kev_col = first_existing_column(df, ["unc_Qalpha_keV", "Qalpha_err_keV", "unc_qa", "Q_alpha_err_keV"])
    if qerr_mev_col:
        out["Qalpha_err_MeV"] = to_numeric_series(df, qerr_mev_col)
    elif qerr_kev_col:
        out["Qalpha_err_MeV"] = to_numeric_series(df, qerr_kev_col) / 1000.0
    else:
        out["Qalpha_err_MeV"] = math.nan

    # Total half-life and alpha branch.
    out["T_half_s"] = to_numeric_series(df, first_existing_column(df, ["T12_total_s", "T_half_s", "half_life_sec", "T12_s"]))
    out["T_half_err_s"] = to_numeric_series(df, first_existing_column(df, ["unc_T12_total_s", "T_half_err_s", "unc_hls", "T12_total_err_s"]))
    out["alpha_branch"] = to_numeric_series(df, first_existing_column(df, ["alpha_branch", "alpha_BR", "b_alpha", "BR_alpha"]))
    # Older dashboard used percent branch columns.
    b_alpha_percent = to_numeric_series(df, first_existing_column(df, ["alpha_branch_percent", "b_alpha_percent", "Ialpha_percent"]))
    out.loc[out["alpha_branch"].isna() & b_alpha_percent.notna(), "alpha_branch"] = b_alpha_percent[out["alpha_branch"].isna() & b_alpha_percent.notna()] / 100.0
    out["alpha_branch_percent"] = out["alpha_branch"] * 100.0
    out["alpha_branch_unc"] = to_numeric_series(df, first_existing_column(df, ["alpha_branch_unc", "alpha_BR_unc", "b_alpha_unc"]))
    b_alpha_unc_percent = to_numeric_series(df, first_existing_column(df, ["alpha_branch_unc_percent", "b_alpha_percent_unc", "unc_b_alpha_percent"]))
    out.loc[out["alpha_branch_unc"].isna() & b_alpha_unc_percent.notna(), "alpha_branch_unc"] = b_alpha_unc_percent[out["alpha_branch_unc"].isna() & b_alpha_unc_percent.notna()] / 100.0
    out["alpha_branch_unc_percent"] = out["alpha_branch_unc"] * 100.0

    # Alpha partial half-life and uncertainty.
    talpha_col = first_existing_column(df, ["T12_alpha_s", "T_alpha_half_s", "Talpha_s", "T_alpha_s"])
    out["T_alpha_half_s"] = to_numeric_series(df, talpha_col)
    need_talpha = out["T_alpha_half_s"].isna() & out["T_half_s"].notna() & out["alpha_branch"].notna() & (out["alpha_branch"] > 0)
    out.loc[need_talpha, "T_alpha_half_s"] = out.loc[need_talpha, "T_half_s"] / out.loc[need_talpha, "alpha_branch"]

    talpha_err_col = first_existing_column(df, ["unc_T12_alpha_s", "T_alpha_half_err_s", "Talpha_err_s", "unc_T_alpha_half_s"])
    out["T_alpha_half_err_s"] = to_numeric_series(df, talpha_err_col)
    need_err = out["T_alpha_half_err_s"].isna() & out["T_alpha_half_s"].notna()
    rel2 = pd.Series(0.0, index=out.index)
    used = pd.Series(False, index=out.index)
    ok_t = out["T_half_s"].gt(0) & out["T_half_err_s"].ge(0)
    rel2.loc[ok_t] += (out.loc[ok_t, "T_half_err_s"] / out.loc[ok_t, "T_half_s"]) ** 2
    used.loc[ok_t] = True
    ok_b = out["alpha_branch"].gt(0) & out["alpha_branch_unc"].ge(0)
    rel2.loc[ok_b] += (out.loc[ok_b, "alpha_branch_unc"] / out.loc[ok_b, "alpha_branch"]) ** 2
    used.loc[ok_b] = True
    calc_err = out["T_alpha_half_s"] * rel2.pow(0.5)
    out.loc[need_err & used, "T_alpha_half_err_s"] = calc_err[need_err & used]

    out["T_alpha_limit"] = df[first_existing_column(df, ["T12_alpha_limit", "T_alpha_half_qualifier", "alpha_half_life_limit"])] if first_existing_column(df, ["T12_alpha_limit", "T_alpha_half_qualifier", "alpha_half_life_limit"]) else ""
    out["alpha_branch_flag"] = df[first_existing_column(df, ["alpha_branch_flag", "alpha_BR_flag"])] if first_existing_column(df, ["alpha_branch_flag", "alpha_BR_flag"]) else ""

    # Optional SF columns, retained only if present.
    out["b_SF_percent"] = to_numeric_series(df, first_existing_column(df, ["b_SF_percent", "SF_branch_percent", "sf_branch_percent"]))
    out["T_SF_half_s"] = to_numeric_series(df, first_existing_column(df, ["T_SF_half_s", "T12_SF_s", "TSF_s"]))
    out["T_SF_half_err_s"] = to_numeric_series(df, first_existing_column(df, ["T_SF_half_err_s", "unc_T12_SF_s", "TSF_err_s"]))

    # Information fields for hover/click panel.
    for new_col, candidates in {
        "decay_summary": ["decay_summary", "decay_modes_raw", "Modes"],
        "source": ["source", "data_source"],
        "note": ["note", "comments", "Note"],
        "Jpi_raw": ["Jpi_raw", "jp", "Jpi", "spin_parity"],
        "Extraction_date": ["Extraction_date", "extraction_date"],
        "ENSDFauthors": ["ENSDFauthors", "ENSDF_authors"],
    }.items():
        col = first_existing_column(df, candidates)
        out[new_col] = df[col] if col else ""

    return out.reset_index(drop=True), n_excluded


def dataframe_to_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        record = {col: clean_value(row.get(col)) for col in df.columns}
        isotope = str(record.get("nuclide") or f"{record.get('A', '')}{record.get('symbol', '')}")
        symbol = str(record.get("symbol") or "")
        record["isotope_label"] = isotope
        record["element"] = symbol
        record["point_label"] = f"{isotope}  (Z={record.get('Z')}, N={record.get('N')})"
        record["Qalpha_display"] = format_with_unc(record.get("Qalpha_MeV"), record.get("Qalpha_err_MeV"), "MeV")
        record["Talpha_display"] = format_with_unc(record.get("T_alpha_half_s"), record.get("T_alpha_half_err_s"), "s", record.get("T_alpha_limit"))
        record["Ttotal_display"] = format_with_unc(record.get("T_half_s"), record.get("T_half_err_s"), "s")
        record["branch_display"] = format_with_unc(record.get("alpha_branch_percent"), record.get("alpha_branch_unc_percent"), "%", record.get("alpha_branch_flag"))
        record["TSF_display"] = format_with_unc(record.get("T_SF_half_s"), record.get("T_SF_half_err_s"), "s")
        records.append(record)
    return records


def available_variables(df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    variables: Dict[str, Dict[str, Any]] = {}
    if "Qalpha_MeV" in df and pd.to_numeric(df["Qalpha_MeV"], errors="coerce").notna().any():
        variables["Qalpha"] = {
            "label": "Qα",
            "title": "Ground-state Qα systematics",
            "ycol": "Qalpha_MeV",
            "err_plus": "Qalpha_err_MeV",
            "err_minus": "Qalpha_err_MeV",
            "unit": "MeV",
            "axis_title": "Q<sub>α</sub> (MeV)",
            "logy": False,
        }
    if "T_alpha_half_s" in df and pd.to_numeric(df["T_alpha_half_s"], errors="coerce").notna().any():
        variables["Talpha"] = {
            "label": "Tα1/2",
            "title": "Ground-state partial α half-lives",
            "ycol": "T_alpha_half_s",
            "err_plus": "T_alpha_half_err_s",
            "err_minus": "T_alpha_half_err_s",
            "limit_col": "T_alpha_limit",
            "unit": "s",
            "axis_title": "T<sup>α</sup><sub>1/2</sub> (s)",
            "logy": True,
        }
    if "T_SF_half_s" in df and pd.to_numeric(df["T_SF_half_s"], errors="coerce").notna().any():
        variables["TSF"] = {
            "label": "TSF1/2",
            "title": "Ground-state partial SF half-lives",
            "ycol": "T_SF_half_s",
            "err_plus": "T_SF_half_err_s",
            "err_minus": "T_SF_half_err_s",
            "unit": "s",
            "axis_title": "T<sup>SF</sup><sub>1/2</sub> (s)",
            "logy": True,
        }
    if not variables:
        raise ValueError("No plottable columns found. Need Qalpha_MeV and/or T12_alpha_s/T_alpha_half_s.")
    return variables


def make_html(records: List[Dict[str, Any]], variables: Dict[str, Dict[str, Any]], shell_n: int, title: str) -> str:
    present_elements = sorted({r["element"] for r in records if r.get("element")}, key=lambda s: Z_BY_SYMBOL.get(s, 999))
    isotope_order = sorted(
        [r.get("isotope_label") for r in records if r.get("isotope_label")],
        key=lambda x: (int("".join(ch for ch in str(x) if ch.isdigit()) or 0), str(x)),
    )
    default_var = "Qalpha" if "Qalpha" in variables else next(iter(variables))
    payload = {
        "records": records,
        "elements": present_elements,
        "isotopes": isotope_order,
        "colors": {el: ELEMENT_COLORS.get(el, "#000000") for el in present_elements},
        "symbols": {el: ELEMENT_SYMBOLS.get(el, "circle") for el in present_elements},
        "variables": variables,
        "defaultVar": default_var,
        "shellN": shell_n,
        "title": title,
    }
    return HTML_TEMPLATE.replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False, allow_nan=False))


def read_input_csv(path: Path) -> pd.DataFrame:
    """Read v2 CSV files and padded Ackermann v1 CSV output consistently."""
    df = pd.read_csv(path, skipinitialspace=True)
    stripped_columns = [str(col).strip() for col in df.columns]
    seen = set()
    duplicates = []
    for col in stripped_columns:
        if col in seen and col not in duplicates:
            duplicates.append(col)
        seen.add(col)
    if duplicates:
        raise ValueError(f"Duplicate CSV columns after trimming whitespace: {duplicates}")
    df.columns = stripped_columns

    for col in df.columns:
        if is_string_dtype(df[col]) or is_object_dtype(df[col]):
            df[col] = df[col].map(lambda value: value.strip() if isinstance(value, str) else value)
    return df


HTML_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Alpha systematics interactive dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  :root {
    --bg: #f3f4f6; --card: #ffffff; --text: #111827;
    --muted: #6b7280; --border: #d1d5db;
  }
  html, body { margin: 0; background: var(--bg); color: var(--text); font-family: "Times New Roman", "STIXGeneral", "DejaVu Serif", serif; }
  .layout { display: grid; grid-template-columns: minmax(760px, 1fr) 370px; gap: 14px; min-height: 100vh; padding: 14px; box-sizing: border-box; }
  .card { background: var(--card); border: 1px solid var(--border); border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.06); }
  .plot-card { min-width: 0; overflow: hidden; }
  #plot { width: 100%; height: calc(100vh - 30px); min-height: 680px; }
  .controller { max-height: calc(100vh - 30px); overflow: auto; padding: 12px; box-sizing: border-box; }
  .title { font-size: 20px; font-weight: 700; margin: 0 0 4px; }
  .hint { margin: 0 0 12px; color: var(--muted); font-size: 13px; line-height: 1.35; }
  .section { border-top: 1px solid var(--border); padding-top: 12px; margin-top: 12px; }
  .section-title { font-weight: 700; margin-bottom: 8px; font-size: 15px; }
  .row { display: flex; gap: 8px; align-items: center; margin: 6px 0; }
  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
  .grid3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; }
  button, input, select { font: inherit; border: 1px solid #9ca3af; border-radius: 6px; background: #ffffff; padding: 6px 8px; box-sizing: border-box; }
  input, select { width: 100%; }
  button { cursor: pointer; }
  .button-main.active, .action.active { background: #111827; color: #ffffff; border-color: #111827; font-weight: 700; }
  .action { width: 100%; }
  .element-list { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }
  .element-button { display: flex; align-items: center; justify-content: space-between; border-left-width: 6px; color: var(--text); }
  .element-button.inactive { color: #9ca3af; background: #f9fafb; font-weight: 400; opacity: 0.70; }
  .small { color: var(--muted); font-size: 12px; }
  .details { font-size: 13px; line-height: 1.35; background: #f9fafb; border: 1px solid var(--border); border-radius: 8px; padding: 10px; overflow-x: auto; }
  .details b { font-size: 15px; }
  .kv { display: grid; grid-template-columns: 135px 1fr; gap: 3px 8px; margin-top: 8px; }
  .kv span:nth-child(odd) { color: var(--muted); }
  @media (max-width: 1100px) { .layout { grid-template-columns: 1fr; } .controller { max-height: unset; } #plot { height: 700px; } }
</style>
</head>
<body>
<main class="layout">
  <section class="card plot-card"><div id="plot"></div></section>
  <aside class="card controller">
    <p class="title" id="dash-title">Alpha systematics viewer</p>
    <p class="hint">ROOT-like 操作：工具栏选 zoom/pan；鼠标滚轮缩放；框选放大；拖拽平移；双击复位。右侧可手动设置横轴/纵轴范围。</p>

    <div class="section">
      <div class="section-title">Physical quantity</div>
      <div class="grid3" id="var-buttons"></div>
    </div>

    <div class="section">
      <div class="section-title">Element visibility</div>
      <div class="grid2">
        <button class="action" id="show-all" type="button">Show all</button>
        <button class="action" id="hide-all" type="button">Hide all</button>
      </div>
      <div class="element-list" id="element-list"></div>
    </div>

    <div class="section">
      <div class="section-title">Nuclide filter</div>
      <select id="nuclide-select"></select>
      <div class="row"><input id="search" type="text" placeholder="Search, e.g. 254No or No" /></div>
      <div class="grid2">
        <button class="action" id="clear-nuclide" type="button">Clear nuclide</button>
        <button class="action" id="focus-selected" type="button">Focus selected</button>
      </div>
      <div class="small">选择某个核素后只显示该点；Search 支持部分匹配。</div>
    </div>

    <div class="section">
      <div class="section-title">Axis range</div>
      <div class="grid2">
        <input id="xmin" type="number" step="1" placeholder="N min" />
        <input id="xmax" type="number" step="1" placeholder="N max" />
        <input id="ymin" type="number" step="any" placeholder="Y min" />
        <input id="ymax" type="number" step="any" placeholder="Y max" />
      </div>
      <div class="grid2" style="margin-top:8px;">
        <button class="action" id="apply-range" type="button">Apply range</button>
        <button class="action" id="reset-range" type="button">Reset range</button>
      </div>
      <div class="small">半衰期图为 log 纵轴，Y min/Y max 仍填写真实秒数，例如 1e-3 到 1e3。</div>
    </div>

    <div class="section">
      <div class="section-title">Clicked point</div>
      <div class="details" id="details">Click a point to inspect one nuclide.</div>
    </div>
  </aside>
</main>

<script>
const payload = __PAYLOAD__;
const records = payload.records;
const variables = payload.variables;
const elements = payload.elements;
const colors = payload.colors;
const symbols = payload.symbols;
const shellN = payload.shellN;
let currentVar = payload.defaultVar;
let activeElements = new Set(elements);
let selectedNuclide = "";
let searchText = "";
let manualRange = null;

document.getElementById("dash-title").textContent = payload.title;

const config = {
  responsive: true,
  displaylogo: false,
  scrollZoom: true,
  doubleClick: "reset",
  modeBarButtonsToAdd: ["drawline", "eraseshape"],
  toImageButtonOptions: {format: "png", filename: "alpha_systematics_dashboard", height: 900, width: 1200, scale: 3}
};

function finite(v) { return v !== null && v !== undefined && Number.isFinite(Number(v)); }

function pointPassesFilter(r, ycol) {
  if (!activeElements.has(r.element)) return false;
  if (!finite(r.N) || !finite(r[ycol])) return false;
  if (variables[currentVar].logy && Number(r[ycol]) <= 0) return false;
  if (selectedNuclide && r.isotope_label !== selectedNuclide) return false;
  if (searchText) {
    const hay = `${r.isotope_label} ${r.element} ${r.Z} ${r.N} ${r.A}`.toLowerCase();
    if (!hay.includes(searchText.toLowerCase())) return false;
  }
  return true;
}

function makeHover(r, varKey) {
  return [
    `<b>${r.isotope_label}</b>` + ` (Z=${r.Z}, N=${r.N})`,
    `A=${r.A}`,
    `Qα: ${r.Qalpha_display || ""}`,
    `T total: ${r.Ttotal_display || ""}`,
    `bα: ${r.branch_display || ""}`,
    `Tα1/2: ${r.Talpha_display || ""}`,
    `TSF1/2: ${r.TSF_display || ""}`,
    `State: ${r.state || ""}`,
    `Jπ: ${r.Jpi_raw || ""}`,
    `Decay: ${r.decay_summary || ""}`,
    `Source: ${r.source || ""}`,
    `Note: ${r.note || ""}`
  ].join("<br>");
}

function groupRecords(varKey) {
  const ycol = variables[varKey].ycol;
  const grouped = new Map();
  elements.forEach(el => grouped.set(el, []));
  records.forEach(r => { if (pointPassesFilter(r, ycol)) grouped.get(r.element)?.push(r); });
  return grouped;
}

function limitSymbol(r, baseSymbol, v) {
  const lim = String(r[v.limit_col || ""] || "").toLowerCase();
  if (lim === "lower") return "triangle-up";
  if (lim === "upper") return "triangle-down";
  return baseSymbol;
}

function makeTraces(varKey) {
  const v = variables[varKey];
  const grouped = groupRecords(varKey);
  const traces = [];
  elements.forEach(el => {
    const rows = (grouped.get(el) || []).sort((a, b) => Number(a.N) - Number(b.N));
    if (!rows.length) return;
    const x = rows.map(r => Number(r.N));
    const y = rows.map(r => Number(r[v.ycol]));
    const trace = {
      type: "scatter",
      mode: "markers+lines",
      name: el,
      x: x,
      y: y,
      customdata: rows,
      text: rows.map(r => makeHover(r, varKey)),
      hovertemplate: "%{text}<extra></extra>",
      marker: {
        color: colors[el] || "#000000",
        symbol: rows.map(r => limitSymbol(r, symbols[el] || "circle", v)),
        size: 9,
        line: {color: colors[el] || "#000000", width: 1.2}
      },
      line: {color: colors[el] || "#000000", width: 1.2},
      connectgaps: false
    };
    if (v.err_plus && rows.some(r => finite(r[v.err_plus]) && Number(r[v.err_plus]) > 0)) {
      const plus = [];
      const minus = [];
      rows.forEach(r => {
        const yy = Number(r[v.ycol]);
        const ep = finite(r[v.err_plus]) ? Math.max(0, Number(r[v.err_plus])) : 0;
        let em = finite(r[v.err_minus]) ? Math.max(0, Number(r[v.err_minus])) : ep;
        if (v.logy && yy - em <= 0) em = Math.max(0, yy * 0.999);
        plus.push(ep);
        minus.push(em);
      });
      trace.error_y = {type: "data", symmetric: false, array: plus, arrayminus: minus, visible: true, color: colors[el] || "#000000", thickness: 1.2, width: 3};
    }
    traces.push(trace);
  });
  return traces;
}

function autoRange(varKey) {
  const v = variables[varKey];
  const xs = [];
  const ys = [];
  records.forEach(r => {
    if (pointPassesFilter(r, v.ycol)) { xs.push(Number(r.N)); ys.push(Number(r[v.ycol])); }
  });
  const xmin = xs.length ? Math.min(...xs) - 1 : 148;
  const xmax = xs.length ? Math.max(...xs) + 1 : 156;
  let ymin, ymax;
  if (!ys.length) { ymin = v.logy ? 1e-9 : 0; ymax = v.logy ? 1e9 : 1; }
  else if (v.logy) {
    const logs = ys.filter(y => y > 0).map(y => Math.log10(y));
    ymin = Math.pow(10, Math.floor(Math.min(...logs) - 0.5));
    ymax = Math.pow(10, Math.ceil(Math.max(...logs) + 0.5));
  } else {
    const lo = Math.min(...ys); const hi = Math.max(...ys); const pad = Math.max((hi - lo) * 0.08, 0.05);
    ymin = lo - pad; ymax = hi + pad;
  }
  return {x: [xmin, xmax], y: [ymin, ymax]};
}

function axisRangeForPlotly(varKey, range) {
  const v = variables[varKey];
  const yr = range.y.slice();
  if (v.logy) { yr[0] = Math.log10(Math.max(yr[0], 1e-300)); yr[1] = Math.log10(Math.max(yr[1], 1e-300)); }
  return {x: range.x, y: yr};
}

function makeLayout(varKey) {
  const v = variables[varKey];
  const range = manualRange || autoRange(varKey);
  const pr = axisRangeForPlotly(varKey, range);
  return {
    template: "plotly_white",
    title: {text: v.title, x: 0.02, xanchor: "left", font: {size: 22}},
    margin: {l: 90, r: 34, t: 62, b: 76},
    font: {family: "Times New Roman, STIXGeneral, DejaVu Serif, serif", size: 16},
    dragmode: "zoom",
    hovermode: "closest",
    xaxis: {title: {text: "Neutron number <i>N</i>", font: {size: 23}}, range: pr.x, dtick: 1, ticks: "inside", mirror: "ticks", showline: true, linecolor: "#000000", linewidth: 1.5, showgrid: true, gridcolor: "rgba(0,0,0,0.10)", zeroline: false, rangeslider: {visible: true, thickness: 0.08}},
    yaxis: {title: {text: v.axis_title, font: {size: 23}}, type: v.logy ? "log" : "linear", range: pr.y, ticks: "inside", mirror: "ticks", showline: true, linecolor: "#000000", linewidth: 1.5, showgrid: true, gridcolor: "rgba(0,0,0,0.10)", zeroline: false},
    legend: {x: 1.01, y: 1, xanchor: "left", yanchor: "top", bgcolor: "rgba(255,255,255,0)", borderwidth: 0, font: {size: 14}},
    shapes: [{type: "line", xref: "x", yref: "paper", x0: shellN, x1: shellN, y0: 0, y1: 1, line: {color: "#7f7f7f", width: 1.3, dash: "dot"}}],
    annotations: [{x: shellN + 0.12, y: 0.95, xref: "x", yref: "paper", text: `<i>N</i>=${shellN}`, showarrow: false, font: {size: 16, color: "#111827"}, xanchor: "left"}]
  };
}

function updatePlot(keepRange=false) {
  if (!keepRange) manualRange = null;
  Plotly.react("plot", makeTraces(currentVar), makeLayout(currentVar), config).then(syncInputs);
  updateVariableButtons();
  updateElementButtons();
}

function syncInputs() {
  const range = manualRange || autoRange(currentVar);
  document.getElementById("xmin").value = Math.round(range.x[0]);
  document.getElementById("xmax").value = Math.round(range.x[1]);
  document.getElementById("ymin").value = Number(range.y[0]).toExponential(4);
  document.getElementById("ymax").value = Number(range.y[1]).toExponential(4);
}

function updateVariableButtons() {
  Object.keys(variables).forEach(k => {
    const btn = document.getElementById(`btn-${k}`);
    if (btn) btn.classList.toggle("active", currentVar === k);
  });
}

function updateElementButtons() {
  elements.forEach(el => {
    const btn = document.getElementById(`el-${el}`);
    if (btn) btn.classList.toggle("inactive", !activeElements.has(el));
  });
}

function buildController() {
  const varHost = document.getElementById("var-buttons");
  Object.keys(variables).forEach(k => {
    const btn = document.createElement("button");
    btn.className = "button-main";
    btn.id = `btn-${k}`;
    btn.type = "button";
    btn.textContent = variables[k].label;
    btn.addEventListener("click", () => { currentVar = k; manualRange = null; updatePlot(true); });
    varHost.appendChild(btn);
  });

  const elHost = document.getElementById("element-list");
  elements.forEach(el => {
    const count = records.filter(r => r.element === el).length;
    const btn = document.createElement("button");
    btn.type = "button"; btn.id = `el-${el}`; btn.className = "element-button";
    btn.style.borderLeftColor = colors[el] || "#000";
    btn.innerHTML = `<span>${el}</span><span class="small">${count}</span>`;
    btn.addEventListener("click", () => { if (activeElements.has(el)) activeElements.delete(el); else activeElements.add(el); updatePlot(false); });
    elHost.appendChild(btn);
  });

  const sel = document.getElementById("nuclide-select");
  const all = document.createElement("option"); all.value = ""; all.textContent = "All nuclides"; sel.appendChild(all);
  payload.isotopes.forEach(iso => { const opt = document.createElement("option"); opt.value = iso; opt.textContent = iso; sel.appendChild(opt); });
}

function showDetails(r) {
  const html = `
    <b>${r.isotope_label}</b> <span class="small">Z=${r.Z}, N=${r.N}, A=${r.A}</span>
    <div class="kv">
      <span>Qα</span><span>${r.Qalpha_display || ""}</span>
      <span>T total</span><span>${r.Ttotal_display || ""}</span>
      <span>bα</span><span>${r.branch_display || ""}</span>
      <span>Tα1/2</span><span>${r.Talpha_display || ""}</span>
      <span>TSF1/2</span><span>${r.TSF_display || ""}</span>
      <span>State</span><span>${r.state || ""}</span>
      <span>Jπ</span><span>${r.Jpi_raw || ""}</span>
      <span>Decay</span><span>${r.decay_summary || ""}</span>
      <span>Source</span><span>${r.source || ""}</span>
      <span>ENSDF</span><span>${r.ENSDFauthors || ""}</span>
      <span>Note</span><span>${r.note || ""}</span>
    </div>`;
  document.getElementById("details").innerHTML = html;
}

function attachEvents() {
  document.getElementById("show-all").addEventListener("click", () => { activeElements = new Set(elements); updatePlot(false); });
  document.getElementById("hide-all").addEventListener("click", () => { activeElements = new Set(); updatePlot(false); });
  document.getElementById("nuclide-select").addEventListener("change", e => { selectedNuclide = e.target.value; updatePlot(false); });
  document.getElementById("search").addEventListener("input", e => { searchText = e.target.value.trim(); selectedNuclide = ""; document.getElementById("nuclide-select").value = ""; updatePlot(false); });
  document.getElementById("clear-nuclide").addEventListener("click", () => { selectedNuclide = ""; searchText = ""; document.getElementById("nuclide-select").value = ""; document.getElementById("search").value = ""; updatePlot(false); });
  document.getElementById("focus-selected").addEventListener("click", () => {
    const target = selectedNuclide || document.getElementById("search").value.trim();
    if (!target) return;
    const hits = records.filter(r => r.isotope_label.toLowerCase().includes(target.toLowerCase()));
    if (!hits.length) return;
    const r = hits[0]; selectedNuclide = r.isotope_label; document.getElementById("nuclide-select").value = selectedNuclide; searchText = ""; document.getElementById("search").value = ""; activeElements = new Set([r.element]);
    const v = variables[currentVar]; const y = finite(r[v.ycol]) ? Number(r[v.ycol]) : null;
    if (y !== null) { manualRange = {x: [Number(r.N) - 2, Number(r.N) + 2], y: [v.logy ? y / 5 : y - Math.abs(y) * 0.12 - 0.1, v.logy ? y * 5 : y + Math.abs(y) * 0.12 + 0.1]}; }
    updatePlot(true); showDetails(r);
  });
  document.getElementById("apply-range").addEventListener("click", () => {
    const xmin = Number(document.getElementById("xmin").value); const xmax = Number(document.getElementById("xmax").value);
    const ymin = Number(document.getElementById("ymin").value); const ymax = Number(document.getElementById("ymax").value);
    if (![xmin, xmax, ymin, ymax].every(Number.isFinite) || xmax <= xmin || ymax <= ymin) return;
    manualRange = {x: [xmin, xmax], y: [ymin, ymax]}; updatePlot(true);
  });
  document.getElementById("reset-range").addEventListener("click", () => { manualRange = null; updatePlot(false); });
  document.getElementById("plot").on("plotly_click", ev => { const p = ev.points && ev.points[0]; if (p && p.customdata) showDetails(p.customdata); });
}

buildController();
attachEvents();
updatePlot(true);
</script>
</body>
</html>
'''


def find_default_csv(here: Path) -> Path:
    alpha_name = "alpha_systematics_Z100_104_N148_156.csv"
    candidates = [
        here / "alpha_systematics_out" / alpha_name,
        here / "out" / alpha_name,
        here.parent / "alpha_systematics_out" / alpha_name,
        here.parent / "out" / alpha_name,
        here / "out" / "ackermann_ground_states.csv",
        here.parent / "out" / "ackermann_ground_states.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate ROOT-like interactive HTML from alpha-systematics CSV.")
    here = Path(__file__).resolve().parent
    parser.add_argument("--csv", type=Path, default=None, help="Input CSV from make_alpha_systematics.py v2")
    parser.add_argument("--html", type=Path, default=None, help="Output HTML file")
    parser.add_argument("--keep-estimated", action="store_true", help="Keep rows carrying '#' estimated/systematics markers")
    parser.add_argument("--shell-n", type=int, default=152, help="Vertical reference line neutron number")
    parser.add_argument("--title", default="Alpha systematics viewer", help="Dashboard title")
    args = parser.parse_args()

    csv_path = (args.csv if args.csv is not None else find_default_csv(here)).resolve()
    out_html = (args.html if args.html is not None else csv_path.parent / "alpha_interactive_dashboard.html").resolve()
    if not csv_path.exists():
        raise FileNotFoundError(
            f"CSV not found: {csv_path}\n"
            "Run make_alpha_systematics.py first, or pass --csv path/to/alpha_systematics_Z100_104_N148_156.csv"
        )

    raw = read_input_csv(csv_path)
    df, n_excluded = normalize_dataframe(raw, keep_estimated=args.keep_estimated)
    variables = available_variables(df)
    records = dataframe_to_records(df)
    html = make_html(records, variables, shell_n=args.shell_n, title=args.title)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")

    print(f"Input CSV: {csv_path}")
    print(f"Records:   {len(records)}")
    if n_excluded:
        print(f"Excluded estimated '#' rows: {n_excluded}  (use --keep-estimated to include them)")
    print(f"Variables: {', '.join(variables.keys())}")
    print(f"HTML:      {out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
