#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a transferable ROOT-like interactive HTML dashboard from
Ackermann ground-state CSV output.

Default usage from tool/NSysByAckerman:
    python make_ackermann_interactive_dashboard.py

Input default:
    out/ackermann_ground_states.csv

Output default:
    out/ackermann_interactive_dashboard.html

Interactive features
--------------------
* Switch among Q_alpha, total T_1/2, T_1/2^alpha, and T_1/2^SF.
* Toggle whether Q_alpha / half-life values marked by # are shown.
* Filter by element buttons, exact nuclide selector, or text search.
* ROOT-like axis interaction: wheel zoom, box zoom, pan, double-click reset.
* Manual axis ranges: N_min/N_max and Y_min/Y_max.
* Click a point to show its full table information.

The HTML uses Plotly from CDN. It is a single file except for that JS library.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd


Z_BY_SYMBOL = {
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
}

# Same general palette idea as the uploaded PlotHalfLifeByElement.py, but extended.
ELEMENT_COLORS = {
    "Es": "#6B7280",
    "Fm": "#111827",
    "Md": "#8C564B",
    "No": "#E377C2",
    "Lr": "#7F7F7F",
    "Rf": "#BCBD22",
    "Db": "#17BECF",
    "Sg": "#8C564B",
    "Bh": "#FF7F0E",
    "Hs": "#9467BD",
    "Mt": "#2CA02C",
    "Ds": "#1F77B4",
    "Rg": "#D62728",
    "Cn": "#4B5563",
    "Nh": "#0F766E",
    "Fl": "#9333EA",
    "Mc": "#F97316",
    "Lv": "#DC2626",
    "Ts": "#2563EB",
    "Og": "#16A34A",
}

ELEMENT_SYMBOLS = {
    "Es": "circle",
    "Fm": "square",
    "Md": "diamond",
    "No": "triangle-up",
    "Lr": "triangle-down",
    "Rf": "cross",
    "Db": "x",
    "Sg": "pentagon",
    "Bh": "hexagon",
    "Hs": "diamond-wide",
    "Mt": "triangle-left",
    "Ds": "square",
    "Rg": "circle",
    "Cn": "triangle-right",
    "Nh": "star",
    "Fl": "hourglass",
    "Mc": "bowtie",
    "Lv": "asterisk",
    "Ts": "hash",
    "Og": "circle-open-dot",
}

PARITY_CLASSES = [
    {"key": "even-N", "label": "even-N", "color": "#111827"},
    {"key": "even-Z", "label": "even-Z", "color": "#16A34A"},
    {"key": "odd-Z", "label": "odd-Z", "color": "#DC2626"},
    {"key": "odd-N", "label": "odd-N", "color": "#2563EB"},
]

VARIABLES = {
    "Qalpha": {
        "label": "Qα",
        "title": "Ground-state Qα systematics",
        "ycol": "Qalpha_MeV",
        "err_plus": "Qalpha_err_MeV",
        "err_minus": "Qalpha_err_MeV",
        "estimated_col": "Qalpha_estimated",
        "unit": "MeV",
        "axis_title": "Q<sub>α</sub> (MeV)",
        "logy": False,
        "format": ".4g",
    },
    "Ttotal": {
        "label": "T1/2",
        "title": "Ground-state total half-lives",
        "ycol": "T_half_s",
        "err_plus": "T_half_err_plus_s",
        "err_minus": "T_half_err_minus_s",
        "estimated_col": "T_half_estimated",
        "unit": "s",
        "axis_title": "T<sub>1/2</sub> (s)",
        "logy": True,
        "format": ".4g",
    },
    "Talpha": {
        "label": "Tα1/2",
        "title": "Ground-state partial α half-lives",
        "ycol": "T_alpha_half_s",
        "err_plus": "T_alpha_half_err_plus_s",
        "err_minus": "T_alpha_half_err_minus_s",
        "estimated_col": "T_alpha_half_estimated",
        "unit": "s",
        "axis_title": "T<sup>α</sup><sub>1/2</sub> (s)",
        "logy": True,
        "format": ".4g",
    },
    "TSF": {
        "label": "TSF1/2",
        "title": "Ground-state partial SF half-lives",
        "ycol": "T_SF_half_s",
        "err_plus": "T_SF_half_err_plus_s",
        "err_minus": "T_SF_half_err_minus_s",
        "estimated_col": "T_SF_half_estimated",
        "unit": "s",
        "axis_title": "T<sup>SF</sup><sub>1/2</sub> (s)",
        "logy": True,
        "format": ".4g",
    },
}


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
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return float(value)
    if isinstance(value, int):
        return int(value)
    return str(value)


def truthy(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def finite_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    try:
        out = float(value)
    except Exception:  # noqa: BLE001
        return None
    return out if math.isfinite(out) else None


def parity_class(z: Any, n: Any) -> str:
    zf = finite_float(z)
    nf = finite_float(n)
    if zf is None or nf is None:
        return ""
    zi = int(zf)
    ni = int(nf)
    if zi % 2 == 0 and ni % 2 == 0:
        return "even-Z even-N"
    if zi % 2 == 0 and ni % 2 == 1:
        return "even-Z odd-N"
    if zi % 2 == 1 and ni % 2 == 0:
        return "odd-Z even-N"
    return "odd-Z odd-N"


def parity_filters(z: Any, n: Any) -> List[str]:
    zf = finite_float(z)
    nf = finite_float(n)
    if zf is None or nf is None:
        return []
    zi = int(zf)
    ni = int(nf)
    filters: List[str] = []
    if ni % 2 == 0:
        filters.append("even-N")
    if zi % 2 == 0:
        filters.append("even-Z")
    else:
        filters.append("odd-Z")
    if ni % 2 == 1:
        filters.append("odd-N")
    return filters


def format_with_unc(
    qual: Any,
    value: Any,
    err_plus: Any,
    err_minus: Any,
    unit: str,
    estimated: Any = False,
) -> str:
    vf = finite_float(value)
    if vf is None:
        return ""
    q = "" if qual is None or pd.isna(qual) else str(qual)
    if q in ("<", "<=", ">", ">=", "approx"):
        prefix = {"<": "< ", "<=": "≤ ", ">": "> ", ">=": "≥ ", "approx": "≈ "}[q]
    elif "uncertain" in q:
        prefix = "~ "
    else:
        prefix = ""
    text = f"{prefix}{vf:.6g}"
    ep = finite_float(err_plus)
    em = finite_float(err_minus)
    if ep is not None and ep > 0 and em is not None and em > 0:
        if math.isclose(ep, em, rel_tol=1e-9, abs_tol=0.0):
            text += f" ± {ep:.2g}"
        else:
            text += f" +{ep:.2g}/-{em:.2g}"
    elif ep is not None and ep > 0:
        text += f" +{ep:.2g}"
    elif em is not None and em > 0:
        text += f" -{em:.2g}"
    if unit:
        text = f"{text} {unit}"
    if truthy(estimated):
        text += " (#)"
    return text.strip()


def dataframe_to_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    needed_numeric = [
        "Z", "N", "A",
        "T_half_s", "T_half_err_plus_s", "T_half_err_minus_s", "T_half_err_s",
        "b_alpha_percent", "b_alpha_err_plus_percent", "b_alpha_err_minus_percent", "b_alpha_err_percent",
        "T_alpha_half_s", "T_alpha_half_err_plus_s", "T_alpha_half_err_minus_s", "T_alpha_half_err_s",
        "b_SF_percent", "b_SF_err_plus_percent", "b_SF_err_minus_percent", "b_SF_err_percent",
        "T_SF_half_s", "T_SF_half_err_plus_s", "T_SF_half_err_minus_s", "T_SF_half_err_s",
        "Qalpha_MeV", "Qalpha_err_MeV",
        "discovery_year",
    ]
    for col in needed_numeric:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    records: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        record = {col: clean_value(row.get(col)) for col in df.columns}
        isotope = str(record.get("isotope") or f"{record.get('A', '')}{record.get('symbol', '')}")
        symbol = str(record.get("symbol") or "")
        record["isotope_label"] = isotope
        record["element"] = symbol
        record["parity_class"] = parity_class(record.get("Z"), record.get("N"))
        record["parity_filters"] = parity_filters(record.get("Z"), record.get("N"))
        record["point_label"] = f"{isotope}  (Z={record.get('Z')}, N={record.get('N')})"
        record["Qalpha_display"] = format_with_unc(
            record.get("Qalpha_qualifier"),
            record.get("Qalpha_MeV"),
            record.get("Qalpha_err_MeV"),
            record.get("Qalpha_err_MeV"),
            "MeV",
            record.get("Qalpha_estimated"),
        )
        record["Ttotal_display"] = format_with_unc(
            record.get("T_half_qualifier"),
            record.get("T_half_s"),
            record.get("T_half_err_plus_s"),
            record.get("T_half_err_minus_s"),
            "s",
            record.get("T_half_estimated"),
        )
        record["branch_display"] = format_with_unc(
            record.get("b_alpha_qualifier"),
            record.get("b_alpha_percent"),
            record.get("b_alpha_err_plus_percent"),
            record.get("b_alpha_err_minus_percent"),
            "%",
            record.get("b_alpha_estimated"),
        )
        record["Talpha_display"] = format_with_unc(
            record.get("T_alpha_half_qualifier"),
            record.get("T_alpha_half_s"),
            record.get("T_alpha_half_err_plus_s"),
            record.get("T_alpha_half_err_minus_s"),
            "s",
            record.get("T_alpha_half_estimated"),
        )
        record["sf_branch_display"] = format_with_unc(
            record.get("b_SF_qualifier"),
            record.get("b_SF_percent"),
            record.get("b_SF_err_plus_percent"),
            record.get("b_SF_err_minus_percent"),
            "%",
            record.get("b_SF_estimated"),
        )
        record["TSF_display"] = format_with_unc(
            record.get("T_SF_half_qualifier"),
            record.get("T_SF_half_s"),
            record.get("T_SF_half_err_plus_s"),
            record.get("T_SF_half_err_minus_s"),
            "s",
            record.get("T_SF_half_estimated"),
        )
        records.append(record)
    return records


def make_html(records: List[Dict[str, Any]]) -> str:
    present_elements = sorted({r["element"] for r in records if r.get("element")}, key=lambda s: Z_BY_SYMBOL.get(s, 999))
    isotope_order = sorted(
        [r.get("isotope_label") for r in records if r.get("isotope_label")],
        key=lambda x: (int("".join(ch for ch in str(x) if ch.isdigit()) or 0), str(x)),
    )
    payload = {
        "records": records,
        "elements": present_elements,
        "isotopes": isotope_order,
        "colors": {el: ELEMENT_COLORS.get(el, "#000000") for el in present_elements},
        "symbols": {el: ELEMENT_SYMBOLS.get(el, "circle") for el in present_elements},
        "parityClasses": PARITY_CLASSES,
        "variables": VARIABLES,
    }

    return HTML_TEMPLATE.replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False, allow_nan=False))


HTML_TEMPLATE = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Ackermann ground-state decay dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  :root {
    --bg: #f3f4f6;
    --card: #ffffff;
    --text: #111827;
    --muted: #6b7280;
    --border: #d1d5db;
    --accent: #111827;
  }
  html, body {
    margin: 0;
    background: var(--bg);
    color: var(--text);
    font-family: "Times New Roman", "STIXGeneral", "DejaVu Serif", serif;
  }
  .layout {
    display: grid;
    grid-template-columns: minmax(760px, 1fr) 360px;
    gap: 14px;
    min-height: 100vh;
    padding: 14px;
    box-sizing: border-box;
  }
  .card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.06);
  }
  .plot-card { min-width: 0; overflow: hidden; }
  #plot {
    width: 100%;
    height: calc(100vh - 30px);
    min-height: 680px;
  }
  .controller {
    max-height: calc(100vh - 30px);
    overflow: auto;
    padding: 12px;
    box-sizing: border-box;
  }
  .title {
    font-size: 20px;
    font-weight: 700;
    margin: 0 0 4px;
  }
  .hint {
    margin: 0 0 12px;
    color: var(--muted);
    font-size: 13px;
    line-height: 1.35;
  }
  .section {
    border-top: 1px solid var(--border);
    padding-top: 12px;
    margin-top: 12px;
  }
  .section-title {
    font-weight: 700;
    margin-bottom: 8px;
    font-size: 15px;
  }
  .row { display: flex; gap: 8px; align-items: center; margin: 6px 0; }
  .check-row { display: flex; gap: 8px; align-items: center; margin: 6px 0; }
  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
  .grid3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; }
  button, input, select {
    font: inherit;
    border: 1px solid #9ca3af;
    border-radius: 6px;
    background: #ffffff;
    padding: 6px 8px;
    box-sizing: border-box;
  }
  input, select { width: 100%; }
  input[type="checkbox"] { width: auto; }
  button { cursor: pointer; }
  .button-main.active, .action.active {
    background: #111827;
    color: #ffffff;
    border-color: #111827;
    font-weight: 700;
  }
  .action { width: 100%; }
  .element-list {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 6px;
  }
  .element-button {
    display: flex;
    align-items: center;
    justify-content: space-between;
    border-left-width: 6px;
    color: var(--text);
  }
  .element-button.inactive {
    color: #9ca3af;
    background: #f9fafb;
    font-weight: 400;
    opacity: 0.70;
  }
  .small { color: var(--muted); font-size: 12px; }
  .details {
    font-size: 13px;
    line-height: 1.35;
    background: #f9fafb;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 10px;
    overflow-x: auto;
  }
  .details b { font-size: 15px; }
  .kv { display: grid; grid-template-columns: 135px 1fr; gap: 3px 8px; margin-top: 8px; }
  .kv span:nth-child(odd) { color: var(--muted); }
  @media (max-width: 1100px) {
    .layout { grid-template-columns: 1fr; }
    .controller { max-height: unset; }
    #plot { height: 700px; }
  }
</style>
</head>
<body>
<main class="layout">
  <section class="card plot-card">
    <div id="plot"></div>
  </section>

  <aside class="card controller">
    <p class="title">Ackermann ground-state viewer</p>
    <p class="hint">ROOT-like 操作：工具栏选 zoom/pan；鼠标滚轮缩放；框选放大；拖拽平移；双击复位。右侧也可以手动设置横轴/纵轴范围。</p>

    <div class="section">
      <div class="section-title">Physical quantity</div>
      <div class="grid2" id="var-buttons"></div>
    </div>

    <div class="section">
      <div class="section-title">Data quality</div>
      <label class="check-row">
        <input id="show-estimated" type="checkbox" />
        <span>Show estimated # values</span>
      </label>
    </div>

    <div class="section">
      <div class="section-title">N/Z parity</div>
      <div class="grid2">
        <button class="action" id="show-all-parity" type="button">Show all</button>
        <button class="action" id="hide-all-parity" type="button">Hide all</button>
      </div>
      <div class="element-list" id="parity-list"></div>
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
      <div class="row">
        <input id="search" type="text" placeholder="Search, e.g. 273Ds or Ds" />
      </div>
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
      <div class="small">半衰期图为 log 纵轴，Y min/Y max 仍然填写真实秒数，例如 1e-3 到 1e3。</div>
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
const parityClasses = payload.parityClasses;
let currentVar = "Qalpha";
let activeElements = new Set(elements);
let activeParityClasses = new Set(parityClasses.map(c => c.key));
let selectedNuclide = "";
let searchText = "";
let manualRange = null;
let showEstimated = false;

const config = {
  responsive: true,
  displaylogo: false,
  scrollZoom: true,
  doubleClick: "reset",
  modeBarButtonsToAdd: ["drawline", "eraseshape"],
  toImageButtonOptions: {
    format: "png",
    filename: "ackermann_ground_state_systematics",
    height: 900,
    width: 1200,
    scale: 3
  }
};

function finite(v) {
  return v !== null && v !== undefined && Number.isFinite(Number(v));
}

function truthy(v) {
  return v === true || v === 1 || ["true", "1", "yes", "y"].includes(String(v).toLowerCase());
}

function isEstimatedForVar(r, varKey) {
  const col = variables[varKey].estimated_col;
  return col ? truthy(r[col]) : false;
}

function rowParityFilters(r) {
  if (Array.isArray(r.parity_filters)) return r.parity_filters;
  return r.parity_class ? [r.parity_class] : [];
}

function pointPassesParity(r) {
  const filters = new Set(rowParityFilters(r));
  const zSelected = ["even-Z", "odd-Z"].filter(k => activeParityClasses.has(k));
  const nSelected = ["even-N", "odd-N"].filter(k => activeParityClasses.has(k));
  if (!zSelected.length && !nSelected.length) return false;
  const zPass = !zSelected.length || zSelected.some(k => filters.has(k));
  const nPass = !nSelected.length || nSelected.some(k => filters.has(k));
  return zPass && nPass;
}

function pointPassesFilter(r, varKey) {
  const v = variables[varKey];
  if (!activeElements.has(r.element)) return false;
  if (!pointPassesParity(r)) return false;
  if (!finite(r.N) || !finite(r[v.ycol])) return false;
  if (v.logy && Number(r[v.ycol]) <= 0) return false;
  if (!showEstimated && isEstimatedForVar(r, varKey)) return false;
  if (selectedNuclide && r.isotope_label !== selectedNuclide) return false;
  if (searchText) {
    const hay = `${r.isotope_label} ${r.element} ${r.Z} ${r.N} ${r.A}`.toLowerCase();
    if (!hay.includes(searchText.toLowerCase())) return false;
  }
  return true;
}

function formatForHover(v, unit) {
  if (!finite(v)) return "";
  return `${Number(v).toExponential(4)} ${unit}`;
}

function makeHover(r, varKey) {
  const v = variables[varKey];
  const y = r[v.ycol];
  const estimated = isEstimatedForVar(r, varKey) ? "yes" : "no";
  return [
    `<b>${r.isotope_label}</b>` + ` (Z=${r.Z}, N=${r.N})`,
    `A=${r.A}`,
    `Qα: ${r.Qalpha_display || ""}`,
    `T total: ${r.Ttotal_display || ""}`,
    `bα: ${r.branch_display || ""}`,
    `Tα1/2: ${r.Talpha_display || ""}`,
    `bSF: ${r.sf_branch_display || ""}`,
    `TSF1/2: ${r.TSF_display || ""}`,
    `Class: ${r.parity_class || ""}`,
    `Estimated marker for this quantity: ${estimated}`,
    `Jπ: ${r.Jpi_raw || ""}`,
    `Modes: ${r.decay_modes_raw || ""}`
  ].join("<br>");
}

function groupRecords(varKey) {
  const grouped = new Map();
  elements.forEach(el => grouped.set(el, []));
  records.forEach(r => {
    if (pointPassesFilter(r, varKey)) grouped.get(r.element)?.push(r);
  });
  return grouped;
}

function makeTraces(varKey) {
  const v = variables[varKey];
  const grouped = groupRecords(varKey);
  const traces = [];
  elements.forEach(el => {
    const rows = (grouped.get(el) || []).sort((a, b) => Number(a.N) - Number(b.N));
    if (!rows.length) return;
    const trace = {
      type: "scatter",
      mode: "markers+lines",
      name: el,
      x: rows.map(r => Number(r.N)),
      y: rows.map(r => Number(r[v.ycol])),
      customdata: rows,
      text: rows.map(r => makeHover(r, varKey)),
      hovertemplate: "%{text}<extra></extra>",
      marker: {
        color: colors[el] || "#000000",
        symbol: symbols[el] || "circle",
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
      trace.error_y = {
        type: "data",
        symmetric: false,
        array: plus,
        arrayminus: minus,
        visible: true,
        color: colors[el] || "#000000",
        thickness: 1.2,
        width: 3
      };
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
    if (pointPassesFilter(r, varKey)) {
      xs.push(Number(r.N));
      ys.push(Number(r[v.ycol]));
    }
  });
  const xmin = xs.length ? Math.min(...xs) - 1 : 130;
  const xmax = xs.length ? Math.max(...xs) + 1 : 180;
  let ymin, ymax;
  if (!ys.length) {
    ymin = v.logy ? 1e-9 : 0;
    ymax = v.logy ? 1e9 : 1;
  } else if (v.logy) {
    const logs = ys.filter(y => y > 0).map(y => Math.log10(y));
    const lo = Math.floor(Math.min(...logs) - 0.5);
    const hi = Math.ceil(Math.max(...logs) + 0.5);
    ymin = Math.pow(10, lo);
    ymax = Math.pow(10, hi);
  } else {
    const lo = Math.min(...ys);
    const hi = Math.max(...ys);
    const pad = Math.max((hi - lo) * 0.08, 0.05);
    ymin = lo - pad;
    ymax = hi + pad;
  }
  return {x: [xmin, xmax], y: [ymin, ymax]};
}

function axisRangeForPlotly(varKey, range) {
  const v = variables[varKey];
  const yr = range.y.slice();
  if (v.logy) {
    yr[0] = Math.log10(Math.max(yr[0], 1e-300));
    yr[1] = Math.log10(Math.max(yr[1], 1e-300));
  }
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
    xaxis: {
      title: {text: "Neutron number <i>N</i>", font: {size: 23}},
      range: pr.x,
      dtick: 2,
      ticks: "inside",
      mirror: "ticks",
      showline: true,
      linecolor: "#000000",
      linewidth: 1.5,
      showgrid: true,
      gridcolor: "rgba(0,0,0,0.10)",
      zeroline: false,
      rangeslider: {visible: true, thickness: 0.08}
    },
    yaxis: {
      title: {text: v.axis_title, font: {size: 23}},
      type: v.logy ? "log" : "linear",
      range: pr.y,
      ticks: "inside",
      mirror: "ticks",
      showline: true,
      linecolor: "#000000",
      linewidth: 1.5,
      showgrid: true,
      gridcolor: "rgba(0,0,0,0.10)",
      zeroline: false
    },
    legend: {
      x: 1.01,
      y: 1,
      xanchor: "left",
      yanchor: "top",
      bgcolor: "rgba(255,255,255,0)",
      borderwidth: 0,
      font: {size: 14}
    },
    shapes: [
      {type: "line", xref: "x", yref: "paper", x0: 152, x1: 152, y0: 0, y1: 1, line: {color: "#7f7f7f", width: 1.3, dash: "dot"}},
      {type: "line", xref: "x", yref: "paper", x0: 162, x1: 162, y0: 0, y1: 1, line: {color: "#7f7f7f", width: 1.3, dash: "dot"}}
    ],
    annotations: [
      {x: 152.25, y: 0.88, xref: "x", yref: "paper", text: "<i>N</i>=152", showarrow: false, font: {size: 16, color: "#111827"}, xanchor: "left"},
      {x: 162.25, y: 0.95, xref: "x", yref: "paper", text: "<i>N</i>=162", showarrow: false, font: {size: 16, color: "#111827"}, xanchor: "left"}
    ]
  };
}

function updatePlot(keepRange=false) {
  if (!keepRange) manualRange = null;
  const traces = makeTraces(currentVar);
  const layout = makeLayout(currentVar);
  Plotly.react("plot", traces, layout, config).then(() => {
    syncInputs();
  });
  updateVariableButtons();
  updateElementButtons();
  updateParityButtons();
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

function updateParityButtons() {
  parityClasses.forEach(cls => {
    const btn = document.getElementById(`parity-${cls.key}`);
    if (btn) btn.classList.toggle("inactive", !activeParityClasses.has(cls.key));
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
    btn.addEventListener("click", () => {
      currentVar = k;
      manualRange = null;
      updatePlot(true);
    });
    varHost.appendChild(btn);
  });

  const parityHost = document.getElementById("parity-list");
  parityClasses.forEach(cls => {
    const count = records.filter(r => rowParityFilters(r).includes(cls.key)).length;
    const btn = document.createElement("button");
    btn.type = "button";
    btn.id = `parity-${cls.key}`;
    btn.className = "element-button";
    btn.style.borderLeftColor = cls.color || "#000";
    btn.innerHTML = `<span>${cls.label}</span><span class="small">${count}</span>`;
    btn.addEventListener("click", () => {
      if (activeParityClasses.has(cls.key)) activeParityClasses.delete(cls.key);
      else activeParityClasses.add(cls.key);
      updatePlot(false);
    });
    parityHost.appendChild(btn);
  });

  const elHost = document.getElementById("element-list");
  elements.forEach(el => {
    const count = records.filter(r => r.element === el).length;
    const btn = document.createElement("button");
    btn.type = "button";
    btn.id = `el-${el}`;
    btn.className = "element-button";
    btn.style.borderLeftColor = colors[el] || "#000";
    btn.innerHTML = `<span>${el}</span><span class="small">${count}</span>`;
    btn.addEventListener("click", () => {
      if (activeElements.has(el)) activeElements.delete(el);
      else activeElements.add(el);
      updatePlot(false);
    });
    elHost.appendChild(btn);
  });

  const sel = document.getElementById("nuclide-select");
  const all = document.createElement("option");
  all.value = "";
  all.textContent = "All nuclides";
  sel.appendChild(all);
  payload.isotopes.forEach(iso => {
    const opt = document.createElement("option");
    opt.value = iso;
    opt.textContent = iso;
    sel.appendChild(opt);
  });
}

function showDetails(r) {
  const html = `
    <b>${r.isotope_label}</b> <span class="small">Z=${r.Z}, N=${r.N}, A=${r.A}</span>
    <div class="kv">
      <span>Qα</span><span>${r.Qalpha_display || ""}</span>
      <span>T total</span><span>${r.Ttotal_display || ""}</span>
      <span>bα</span><span>${r.branch_display || ""}</span>
      <span>Tα1/2</span><span>${r.Talpha_display || ""}</span>
      <span>bSF</span><span>${r.sf_branch_display || ""}</span>
      <span>TSF1/2</span><span>${r.TSF_display || ""}</span>
      <span>Class</span><span>${r.parity_class || ""}</span>
      <span># fields</span><span>${r.nonexperimental_fields || ""}</span>
      <span>Jπ</span><span>${r.Jpi_raw || ""}</span>
      <span>Modes</span><span>${r.decay_modes_raw || ""}</span>
      <span>Branches</span><span>${r.branching_raw || ""}</span>
      <span>Discovery</span><span>${r.discovery_year || ""}</span>
    </div>`;
  document.getElementById("details").innerHTML = html;
}

function attachEvents() {
  document.getElementById("show-estimated").addEventListener("change", e => {
    showEstimated = e.target.checked;
    updatePlot(false);
  });
  document.getElementById("show-all-parity").addEventListener("click", () => {
    activeParityClasses = new Set(parityClasses.map(c => c.key));
    updatePlot(false);
  });
  document.getElementById("hide-all-parity").addEventListener("click", () => {
    activeParityClasses = new Set();
    updatePlot(false);
  });
  document.getElementById("show-all").addEventListener("click", () => {
    activeElements = new Set(elements);
    updatePlot(false);
  });
  document.getElementById("hide-all").addEventListener("click", () => {
    activeElements = new Set();
    updatePlot(false);
  });
  document.getElementById("nuclide-select").addEventListener("change", e => {
    selectedNuclide = e.target.value;
    updatePlot(false);
  });
  document.getElementById("search").addEventListener("input", e => {
    searchText = e.target.value.trim();
    selectedNuclide = "";
    document.getElementById("nuclide-select").value = "";
    updatePlot(false);
  });
  document.getElementById("clear-nuclide").addEventListener("click", () => {
    selectedNuclide = "";
    searchText = "";
    document.getElementById("nuclide-select").value = "";
    document.getElementById("search").value = "";
    updatePlot(false);
  });
  document.getElementById("focus-selected").addEventListener("click", () => {
    const target = selectedNuclide || document.getElementById("search").value.trim();
    if (!target) return;
    const hits = records.filter(r => r.isotope_label.toLowerCase().includes(target.toLowerCase()));
    if (!hits.length) return;
    const r = hits[0];
    if (isEstimatedForVar(r, currentVar)) {
      showEstimated = true;
      document.getElementById("show-estimated").checked = true;
    }
    selectedNuclide = r.isotope_label;
    document.getElementById("nuclide-select").value = selectedNuclide;
    searchText = "";
    document.getElementById("search").value = "";
    activeElements = new Set([r.element]);
    const filters = rowParityFilters(r);
    if (filters.length) activeParityClasses = new Set(filters);
    const v = variables[currentVar];
    const y = finite(r[v.ycol]) ? Number(r[v.ycol]) : null;
    if (y !== null) {
      const ypadLow = v.logy ? y / 5 : y - Math.abs(y) * 0.12 - 0.1;
      const ypadHigh = v.logy ? y * 5 : y + Math.abs(y) * 0.12 + 0.1;
      manualRange = {x: [Number(r.N) - 2, Number(r.N) + 2], y: [ypadLow, ypadHigh]};
    }
    updatePlot(true);
    showDetails(r);
  });
  document.getElementById("apply-range").addEventListener("click", () => {
    const xmin = Number(document.getElementById("xmin").value);
    const xmax = Number(document.getElementById("xmax").value);
    const ymin = Number(document.getElementById("ymin").value);
    const ymax = Number(document.getElementById("ymax").value);
    if (![xmin, xmax, ymin, ymax].every(Number.isFinite) || xmax <= xmin || ymax <= ymin) return;
    manualRange = {x: [xmin, xmax], y: [ymin, ymax]};
    updatePlot(true);
  });
  document.getElementById("reset-range").addEventListener("click", () => {
    manualRange = null;
    updatePlot(false);
  });

  document.getElementById("plot").on("plotly_click", ev => {
    const p = ev.points && ev.points[0];
    if (p && p.customdata) showDetails(p.customdata);
  });
}

buildController();
attachEvents();
updatePlot(true);
</script>
</body>
</html>
'''


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate ROOT-like interactive HTML from Ackermann ground-state CSV.")
    here = Path(__file__).resolve().parent
    base = here.parent
    parser.add_argument("--csv", type=Path, default=base / "out" / "ackermann_ground_states.csv", help="Input CSV from make_ackermann_ground_state_systematics.py")
    parser.add_argument("--html", type=Path, default=base / "out" / "ackermann_interactive_dashboard.html", help="Output HTML file")
    args = parser.parse_args()

    csv_path = args.csv.resolve()
    out_html = args.html.resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}\nRun make_ackermann_ground_state_systematics.py first.")

    df = pd.read_csv(csv_path)
    records = dataframe_to_records(df)
    html = make_html(records)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    print(f"Input CSV: {csv_path}")
    print(f"Records:   {len(records)}")
    print(f"HTML:      {out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
