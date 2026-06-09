#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse Ackermann DSAS review LaTeX table (Table isotope_list), keep ground-state
entries only, export a CSV, and plot Q_alpha, total half-lives, partial alpha
half-lives, and partial SF half-lives versus neutron number N.

Default usage from tool/NSysByAckerman:
    python make_ackermann_ground_state_systematics.py

Optional:
    python make_ackermann_ground_state_systematics.py \
        --tex source/ppnp-LaTeX-D_ackermann_DSAS_Rev_revised_final.tex \
        --outdir out

Notes
-----
* The table lists total half-lives T_1/2 and branching ratios in percent.
  Partial half-lives are computed as
      T_1/2(mode) = T_1/2(total) / (branch_percent/100).
* If a branching ratio is a limit, the partial half-life is marked as the
  opposite limit in the CSV. Example: b_SF < 1% => T_1/2^SF > T_1/2/0.01.
* Uncertainties in total half-lives and branching ratios are propagated to the
  partial half-lives when available.
* Isomeric/letter-state rows such as ^{250m}Es, ^{254m1}No, or ^{265a}Sg are excluded
  unless explicitly added in MANUAL_EXTRA_ROWS. Parenthesized uncertain ground-state
  assignments are kept but flagged in `is_parenthesized`.
* Q_alpha values marked by # are kept in the CSV but flagged in `Qalpha_estimated`;
  static Q_alpha plots exclude them by default. Use --include-estimated-qalpha to
  include those systematics/estimated values in the static plot.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd


# Element symbols used in the review table, Z=99--118.
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

SECONDS_BY_UNIT = {
    "ns": 1e-9,
    "us": 1e-6,
    "µs": 1e-6,
    "μs": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
    "sec": 1.0,
    "min": 60.0,
    "h": 3600.0,
    "d": 86400.0,
    "y": 365.25 * 86400.0,
    "yr": 365.25 * 86400.0,
}

MANUAL_EXTRA_ROWS = [
    {
        "isotope": "273Ds^m",
        "symbol": "Ds",
        "Z": 110,
        "N": 163,
        "A": 273,
        "is_parenthesized": False,
        "mass_number_uncertain": False,
        "T_half_raw": "0.16^+0.08_-0.05 ms",
        "T_half_s": 0.16e-3,
        "T_half_err_plus_s": 0.08e-3,
        "T_half_err_minus_s": 0.05e-3,
        "T_half_err_s": 0.065e-3,
        "T_half_unit": "ms",
        "T_half_qualifier": "",
        "T_half_estimated": False,
        "Jpi_raw": "",
        "Jpi_estimated": False,
        "decay_modes_raw": "alpha",
        "branching_raw": "1",
        "b_alpha_percent": 100.0,
        "b_alpha_err_plus_percent": None,
        "b_alpha_err_minus_percent": None,
        "b_alpha_err_percent": None,
        "b_alpha_qualifier": "",
        "b_alpha_estimated": False,
        "T_alpha_half_s": 0.16e-3,
        "T_alpha_half_err_plus_s": 0.08e-3,
        "T_alpha_half_err_minus_s": 0.05e-3,
        "T_alpha_half_err_s": 0.065e-3,
        "T_alpha_half_qualifier": "",
        "T_alpha_half_estimated": False,
        "b_SF_percent": None,
        "b_SF_err_plus_percent": None,
        "b_SF_err_minus_percent": None,
        "b_SF_err_percent": None,
        "b_SF_qualifier": "",
        "b_SF_estimated": False,
        "T_SF_half_s": None,
        "T_SF_half_err_plus_s": None,
        "T_SF_half_err_minus_s": None,
        "T_SF_half_err_s": None,
        "T_SF_half_qualifier": "",
        "T_SF_half_estimated": False,
        "Qalpha_raw": "11.33(2) MeV",
        "Qalpha_keV": 11330.0,
        "Qalpha_MeV": 11.33,
        "Qalpha_err_keV": 20.0,
        "Qalpha_err_MeV": 0.02,
        "Qalpha_estimated": False,
        "Qalpha_is_experimental": True,
        "Qalpha_qualifier": "",
        "has_estimated_marker": False,
        "nonexperimental_fields": "",
        "discovery_year": None,
    },
]


@dataclass
class RawRecord:
    isotope_raw: str
    half_life_raw: str
    spin_raw: str
    qalpha_raw: str
    discovery_raw: str
    mode_branch_rows: List[Tuple[str, str]] = field(default_factory=list)


def strip_unescaped_percent(line: str) -> str:
    r"""Remove LaTeX comments while preserving escaped \%."""
    out = []
    escaped = False
    for ch in line:
        if ch == "%" and not escaped:
            break
        out.append(ch)
        escaped = (ch == "\\" and not escaped)
        if ch != "\\":
            escaped = False
    return "".join(out)


def latex_to_text(s: str) -> str:
    """Convert the small subset of LaTeX used in the table to parseable text."""
    if s is None:
        return ""
    s = s.strip()
    replacements = {
        "~": " ",
        "\\,": " ",
        "\\alpha": "alpha",
        "\\beta": "beta",
        "\\gamma": "gamma",
        "\\mu": "u",
        "\\#": "#",
        "\\leq": "<=",
        "\\geq": ">=",
        "\\approx": "approx ",
        "\\sim": "~",
        "$<$": "<",
        "$>$": ">",
        "$E$": "E",
    }
    for a, b in replacements.items():
        s = s.replace(a, b)

    # Drop citations and formatting commands, but keep command arguments where useful.
    s = re.sub(r"\\cite\{[^{}]*\}", "", s)
    s = re.sub(r"\\href\{[^{}]*\}\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\tiny", "", s)
    s = re.sub(r"\\mathrm\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\text\{([^{}]*)\}", r"\1", s)

    # Remove remaining dollar signs and simple LaTeX structural braces.
    s = s.replace("$", "")
    s = s.replace("{", "").replace("}", "")
    s = s.replace("\\", "")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def split_latex_row(line: str) -> Optional[List[str]]:
    line = strip_unescaped_percent(line).strip()
    if not line or "&" not in line:
        return None
    # The table rows are one line each. Remove final LaTeX row terminator.
    line = re.sub(r"\\\\\s*$", "", line).strip()
    cells = [c.strip() for c in line.split("&")]
    # Pad to the expected longtable width.
    if len(cells) < 8:
        cells += [""] * (8 - len(cells))
    return cells[:8]


def extract_isotope_table(tex: str) -> str:
    label_pos = tex.find(r"\label{tab:isotope_list}")
    if label_pos < 0:
        raise RuntimeError("Could not find \\label{tab:isotope_list} in the LaTeX file.")
    start = tex.rfind(r"\begin{longtable}", 0, label_pos)
    end = tex.find(r"\end{longtable}", label_pos)
    if start < 0 or end < 0:
        raise RuntimeError("Could not isolate the isotope_list longtable.")
    return tex[start:end]


def build_raw_records(table_text: str) -> List[RawRecord]:
    records: List[RawRecord] = []
    current: Optional[RawRecord] = None

    for line in table_text.splitlines():
        cells = split_latex_row(line)
        if cells is None:
            continue
        first_cell = cells[0].strip()
        first_clean = latex_to_text(first_cell)
        is_data_start = bool(re.search(r"\^\{?\d+", first_cell)) and any(
            sym in first_cell for sym in Z_BY_SYMBOL
        )

        if is_data_start:
            current = RawRecord(
                isotope_raw=cells[0],
                half_life_raw=cells[1],
                spin_raw=cells[2],
                qalpha_raw=cells[5],
                discovery_raw=cells[6],
            )
            current.mode_branch_rows.append((cells[3], cells[4]))
            records.append(current)
        elif current is not None:
            # Continuation line. Append extra decay/branching information when present.
            if latex_to_text(cells[3]) or latex_to_text(cells[4]):
                current.mode_branch_rows.append((cells[3], cells[4]))

    return records


def parse_isotope(isotope_raw: str) -> Dict[str, object]:
    # Examples: $^{240}$Es, $^{250m}$Es, $^{254m1}$No, $^{265a}$Sg, ($^{286}$Cn)$^v
    m = re.search(r"\^\{\s*(\d+)(?:\(\d+\))?([A-Za-z]\d*)?\s*\}\$?\s*([A-Z][a-z]?)", isotope_raw)
    if not m:
        raise ValueError(f"Could not parse isotope cell: {isotope_raw}")
    A = int(m.group(1))
    state = m.group(2) or ""
    symbol = m.group(3)
    if symbol not in Z_BY_SYMBOL:
        raise ValueError(f"Unknown element symbol {symbol} in {isotope_raw}")
    Z = Z_BY_SYMBOL[symbol]
    return {
        "A": A,
        "symbol": symbol,
        "Z": Z,
        "N": A - Z,
        "state": state,
        "is_ground_state": state == "",
        "is_parenthesized": isotope_raw.strip().startswith("("),
        "mass_number_uncertain": bool(re.search(r"\^\{\s*\d+\(\d+\)", isotope_raw)),
        "isotope": f"{A}{symbol}",
    }


def qualifier_from_text(raw: str) -> str:
    txt = latex_to_text(raw)
    if "<=" in txt or "≤" in txt:
        return "<="
    if ">=" in txt or "≥" in txt:
        return ">="
    if "<" in txt:
        return "<"
    if ">" in txt:
        return ">"
    if "approx" in txt or "~" in txt:
        return "approx"
    return ""


def first_float(text: str) -> Optional[float]:
    m = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?", text)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def bracket_uncertainty(raw: str) -> Optional[float]:
    """Return parenthetic uncertainty in the same units as the leading number."""
    txt = latex_to_text(raw)
    # Find leading number and following (digits) uncertainty.
    m = re.search(r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))\s*\((\d+(?:\.\d+)?)\)", txt)
    if not m:
        return None
    number_s, unc_s = m.group(1), m.group(2)
    if "." in unc_s:
        return float(unc_s)
    decimals = len(number_s.split(".", 1)[1]) if "." in number_s else 0
    return float(unc_s) * (10.0 ** (-decimals))


def asymmetric_uncertainty(raw: str) -> Tuple[Optional[float], Optional[float]]:
    """Return (+unc, -unc) in the same units as the leading number."""
    txt = re.sub(r"\s+", "", latex_to_text(raw))
    num = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    m = re.search(rf"\^\+({num})_-({num})", txt)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def uncertainty_pair(raw: str) -> Tuple[Optional[float], Optional[float]]:
    """Return (+unc, -unc), accepting asymmetric and parenthetic forms."""
    plus, minus = asymmetric_uncertainty(raw)
    if plus is not None or minus is not None:
        return plus, minus
    sym = bracket_uncertainty(raw)
    if sym is None:
        return None, None
    return sym, sym


def symmetric_error(plus: Optional[float], minus: Optional[float]) -> Optional[float]:
    vals = [v for v in (plus, minus) if v is not None and math.isfinite(v)]
    if not vals:
        return None
    return sum(vals) / len(vals)


def parse_half_life(raw: str) -> Tuple[Optional[float], Optional[float], Optional[float], str, str, bool]:
    txt = latex_to_text(raw)
    qual = qualifier_from_text(raw)
    estimated = "#" in txt
    # Remove common uncertainty forms before unit matching.
    txt_no_unc = re.sub(r"\^\+[^_\s]+_-[^\s]+", " ", txt)
    txt_no_unc = re.sub(r"\([^)]*\)", " ", txt_no_unc)
    value = first_float(txt)
    if value is None:
        return None, None, None, "", qual, estimated
    unit_match = re.search(r"\b(ns|us|µs|μs|ms|min|sec|s|h|d|yr|y)\b", txt_no_unc, re.I)
    # If parentheses removal swallowed the unit, fall back to the original clean text.
    if not unit_match:
        unit_match = re.search(r"\b(ns|us|µs|μs|ms|min|sec|s|h|d|yr|y)\b", txt, re.I)
    if not unit_match:
        return None, None, None, "", qual, estimated
    unit = unit_match.group(1).lower().replace("μ", "µ")
    scale = SECONDS_BY_UNIT[unit]
    unc_plus, unc_minus = uncertainty_pair(raw)
    seconds = value * scale
    err_plus_s = unc_plus * scale if unc_plus is not None else None
    err_minus_s = unc_minus * scale if unc_minus is not None else None
    return seconds, err_plus_s, err_minus_s, unit, qual, estimated


def parse_qalpha(raw: str) -> Tuple[Optional[float], Optional[float], bool, str]:
    """Parse Q_alpha and return (keV, err_keV, estimated_flag, qualifier)."""
    txt = latex_to_text(raw)
    if not txt:
        return None, None, False, ""
    qual = qualifier_from_text(raw)
    estimated = "#" in txt
    value = first_float(txt)
    if value is None:
        return None, None, estimated, qual
    err = bracket_uncertainty(raw)

    # Header says keV, but the newest Lv entries are written as 11.240(15), i.e. MeV.
    if value < 100.0:
        value_keV = value * 1000.0
        err_keV = err * 1000.0 if err is not None else None
    else:
        value_keV = value
        err_keV = err
    return value_keV, err_keV, estimated, qual


def normalize_modes(mode_raw: str) -> List[str]:
    txt = latex_to_text(mode_raw)
    txt = txt.replace("beta^-", "beta-").replace("beta^+", "beta+")
    txt = txt.replace("beta-", "beta-").replace("beta+", "beta+")
    modes = []
    pattern = re.compile(r"alpha|EC|SF|IT|beta\s*[-+]|\bbeta\b", re.I)
    for m in pattern.finditer(txt):
        token = m.group(0).replace(" ", "")
        low = token.lower()
        if low == "alpha":
            norm = "alpha"
        elif low == "ec":
            norm = "EC"
        elif low == "sf":
            norm = "SF"
        elif low == "it":
            norm = "IT"
        elif low.startswith("beta-"):
            norm = "beta-"
        elif low.startswith("beta+"):
            norm = "beta+"
        else:
            norm = "beta"
        modes.append(norm)
    return modes


def parse_branch_token(token_raw: str) -> Tuple[Optional[float], Optional[float], Optional[float], str, str, bool]:
    txt = latex_to_text(token_raw)
    txt = txt.strip().strip("()")
    if not txt:
        return None, None, None, "unknown", txt, False
    qual = qualifier_from_text(token_raw)
    estimated = "#" in txt
    value = first_float(txt)
    if value is None:
        return None, None, None, "unknown", txt, estimated
    err_plus, err_minus = uncertainty_pair(token_raw)
    # Entries such as "100 (?)" are numerically useful but should be flagged.
    if "?" in txt and not qual:
        qual = "uncertain"
    elif "?" in txt and qual:
        qual = f"{qual};uncertain"
    return value, err_plus, err_minus, qual, txt, estimated


def parse_mode_branches(rows: Iterable[Tuple[str, str]]) -> Dict[str, Dict[str, object]]:
    result: Dict[str, Dict[str, object]] = {}
    all_modes_raw = []
    all_branches_raw = []

    for mode_raw, branch_raw in rows:
        modes = normalize_modes(mode_raw)
        branch_txt = latex_to_text(branch_raw)
        branch_tokens = [tok.strip() for tok in branch_txt.split("/")] if branch_txt else []
        all_modes_raw.append(latex_to_text(mode_raw))
        all_branches_raw.append(branch_txt)

        for i, mode in enumerate(modes):
            token = branch_tokens[i] if i < len(branch_tokens) else ""
            value, err_plus, err_minus, qual, token_clean, estimated = parse_branch_token(token)
            # Do not overwrite an already parsed numerical value with an unknown continuation.
            if mode in result and result[mode].get("value") is not None and value is None:
                continue
            result[mode] = {
                "value": value,
                "err_plus": err_plus,
                "err_minus": err_minus,
                "qualifier": qual,
                "raw": token_clean,
                "estimated": estimated,
            }

    result["__raw__"] = {
        "decay_modes_raw": "; ".join(x for x in all_modes_raw if x),
        "branching_raw": "; ".join(x for x in all_branches_raw if x),
    }
    return result


def apply_single_alpha_default(branches: Dict[str, Dict[str, object]]) -> None:
    """Treat a sole alpha decay mode with no listed branch as a 100% alpha branch."""
    mode_keys = {key for key in branches if key != "__raw__"}
    alpha = branches.get("alpha")
    if mode_keys != {"alpha"} or alpha is None:
        return
    if alpha.get("value") is not None:
        return
    alpha["value"] = 100.0
    alpha["err_plus"] = None
    alpha["err_minus"] = None
    alpha["qualifier"] = "inferred_single_alpha"
    alpha["estimated"] = False


def partial_half_life(
    total_s: Optional[float],
    branch_percent: Optional[float],
    branch_qual: str,
    total_qual: str = "",
) -> Tuple[Optional[float], str]:
    if total_s is None or branch_percent is None or branch_percent <= 0:
        return None, ""
    t = total_s / (branch_percent / 100.0)

    # Branching-ratio limits invert into partial half-life limits.
    branch_partial_qual = ""
    if branch_qual in ("<", "<="):
        branch_partial_qual = ">" if branch_qual == "<" else ">="
    elif branch_qual in (">", ">="):
        branch_partial_qual = "<" if branch_qual == ">" else "<="
    elif branch_qual in ("approx", "uncertain") or "uncertain" in str(branch_qual):
        branch_partial_qual = branch_qual

    # A limit on total T_1/2 carries through directly when the branch is exact.
    if total_qual and branch_partial_qual:
        return t, f"{total_qual};{branch_partial_qual}"
    if branch_partial_qual:
        return t, branch_partial_qual
    if total_qual:
        return t, total_qual
    return t, ""


def partial_half_life_errors(
    partial_s: Optional[float],
    total_s: Optional[float],
    total_err_plus_s: Optional[float],
    total_err_minus_s: Optional[float],
    branch_percent: Optional[float],
    branch_err_plus_percent: Optional[float],
    branch_err_minus_percent: Optional[float],
) -> Tuple[Optional[float], Optional[float]]:
    """Propagate T_partial = T_total / branch. Returns (+err, -err) in seconds."""
    if partial_s is None or total_s is None or total_s <= 0:
        return None, None
    if branch_percent is None or branch_percent <= 0:
        return None, None

    rel_plus2 = 0.0
    rel_minus2 = 0.0
    used_plus = False
    used_minus = False

    if total_err_plus_s is not None and math.isfinite(total_err_plus_s) and total_err_plus_s >= 0:
        rel_plus2 += (total_err_plus_s / total_s) ** 2
        used_plus = True
    if total_err_minus_s is not None and math.isfinite(total_err_minus_s) and total_err_minus_s >= 0:
        rel_minus2 += (total_err_minus_s / total_s) ** 2
        used_minus = True

    # The branch term is inverted: a smaller branch gives a larger partial half-life.
    if branch_err_minus_percent is not None and math.isfinite(branch_err_minus_percent) and branch_err_minus_percent >= 0:
        rel_plus2 += (branch_err_minus_percent / branch_percent) ** 2
        used_plus = True
    if branch_err_plus_percent is not None and math.isfinite(branch_err_plus_percent) and branch_err_plus_percent >= 0:
        rel_minus2 += (branch_err_plus_percent / branch_percent) ** 2
        used_minus = True

    err_plus = partial_s * math.sqrt(rel_plus2) if used_plus else None
    err_minus = partial_s * math.sqrt(rel_minus2) if used_minus else None
    return err_plus, err_minus


def records_to_dataframe(records: List[RawRecord]) -> pd.DataFrame:
    rows = []
    for rec in records:
        iso = parse_isotope(rec.isotope_raw)
        if not iso["is_ground_state"]:
            continue

        total_s, total_err_plus_s, total_err_minus_s, unit, tqual, t_est = parse_half_life(rec.half_life_raw)
        q_keV, qerr_keV, q_est, qqual = parse_qalpha(rec.qalpha_raw)
        branches = parse_mode_branches(rec.mode_branch_rows)
        apply_single_alpha_default(branches)
        alpha = branches.get("alpha", {"value": None, "err_plus": None, "err_minus": None, "qualifier": "", "raw": "", "estimated": False})
        sf = branches.get("SF", {"value": None, "err_plus": None, "err_minus": None, "qualifier": "", "raw": "", "estimated": False})
        talpha_s, talpha_qual = partial_half_life(total_s, alpha.get("value"), alpha.get("qualifier", ""), tqual)
        tsf_s, tsf_qual = partial_half_life(total_s, sf.get("value"), sf.get("qualifier", ""), tqual)
        talpha_err_plus_s, talpha_err_minus_s = partial_half_life_errors(
            talpha_s,
            total_s,
            total_err_plus_s,
            total_err_minus_s,
            alpha.get("value"),
            alpha.get("err_plus"),
            alpha.get("err_minus"),
        )
        tsf_err_plus_s, tsf_err_minus_s = partial_half_life_errors(
            tsf_s,
            total_s,
            total_err_plus_s,
            total_err_minus_s,
            sf.get("value"),
            sf.get("err_plus"),
            sf.get("err_minus"),
        )
        jpi_est = "#" in latex_to_text(rec.spin_raw)
        estimated_fields = []
        if q_est:
            estimated_fields.append("Qalpha")
        if t_est:
            estimated_fields.append("T_half")
        if alpha.get("estimated"):
            estimated_fields.append("b_alpha")
        if sf.get("estimated"):
            estimated_fields.append("b_SF")
        if jpi_est:
            estimated_fields.append("Jpi")

        year = first_float(latex_to_text(rec.discovery_raw))
        rows.append(
            {
                "isotope": iso["isotope"],
                "symbol": iso["symbol"],
                "Z": iso["Z"],
                "N": iso["N"],
                "A": iso["A"],
                "is_parenthesized": iso["is_parenthesized"],
                "mass_number_uncertain": iso["mass_number_uncertain"],
                "T_half_raw": latex_to_text(rec.half_life_raw),
                "T_half_s": total_s,
                "T_half_err_plus_s": total_err_plus_s,
                "T_half_err_minus_s": total_err_minus_s,
                "T_half_err_s": symmetric_error(total_err_plus_s, total_err_minus_s),
                "T_half_unit": unit,
                "T_half_qualifier": tqual,
                "T_half_estimated": t_est,
                "Jpi_raw": latex_to_text(rec.spin_raw),
                "Jpi_estimated": jpi_est,
                "decay_modes_raw": branches["__raw__"]["decay_modes_raw"],
                "branching_raw": branches["__raw__"]["branching_raw"],
                "b_alpha_percent": alpha.get("value"),
                "b_alpha_err_plus_percent": alpha.get("err_plus"),
                "b_alpha_err_minus_percent": alpha.get("err_minus"),
                "b_alpha_err_percent": symmetric_error(alpha.get("err_plus"), alpha.get("err_minus")),
                "b_alpha_qualifier": alpha.get("qualifier", ""),
                "b_alpha_estimated": alpha.get("estimated", False),
                "T_alpha_half_s": talpha_s,
                "T_alpha_half_err_plus_s": talpha_err_plus_s,
                "T_alpha_half_err_minus_s": talpha_err_minus_s,
                "T_alpha_half_err_s": symmetric_error(talpha_err_plus_s, talpha_err_minus_s),
                "T_alpha_half_qualifier": talpha_qual,
                "T_alpha_half_estimated": t_est or alpha.get("estimated", False),
                "b_SF_percent": sf.get("value"),
                "b_SF_err_plus_percent": sf.get("err_plus"),
                "b_SF_err_minus_percent": sf.get("err_minus"),
                "b_SF_err_percent": symmetric_error(sf.get("err_plus"), sf.get("err_minus")),
                "b_SF_qualifier": sf.get("qualifier", ""),
                "b_SF_estimated": sf.get("estimated", False),
                "T_SF_half_s": tsf_s,
                "T_SF_half_err_plus_s": tsf_err_plus_s,
                "T_SF_half_err_minus_s": tsf_err_minus_s,
                "T_SF_half_err_s": symmetric_error(tsf_err_plus_s, tsf_err_minus_s),
                "T_SF_half_qualifier": tsf_qual,
                "T_SF_half_estimated": t_est or sf.get("estimated", False),
                "Qalpha_raw": latex_to_text(rec.qalpha_raw),
                "Qalpha_keV": q_keV,
                "Qalpha_MeV": q_keV / 1000.0 if q_keV is not None else None,
                "Qalpha_err_keV": qerr_keV,
                "Qalpha_err_MeV": qerr_keV / 1000.0 if qerr_keV is not None else None,
                "Qalpha_estimated": q_est,
                "Qalpha_is_experimental": not q_est,
                "Qalpha_qualifier": qqual,
                "has_estimated_marker": bool(estimated_fields),
                "nonexperimental_fields": ";".join(estimated_fields),
                "discovery_year": int(year) if year is not None else None,
            }
        )

    rows.extend(dict(row) for row in MANUAL_EXTRA_ROWS)
    df = pd.DataFrame(rows)
    df = df.sort_values(["Z", "N", "A"]).reset_index(drop=True)
    return df


def save_csv(df: pd.DataFrame, outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outpath, index=False, quoting=csv.QUOTE_MINIMAL)


def plot_by_element(
    df: pd.DataFrame,
    ycol: str,
    ylabel: str,
    outstem: Path,
    yerr_col: Optional[str] = None,
    yerr_plus_col: Optional[str] = None,
    yerr_minus_col: Optional[str] = None,
    logy: bool = False,
    title: Optional[str] = None,
    exclude_estimated_col: Optional[str] = None,
) -> None:
    plot_df = df.dropna(subset=["N", ycol]).copy()
    if exclude_estimated_col and exclude_estimated_col in plot_df.columns:
        est = plot_df[exclude_estimated_col].fillna(False)
        if est.dtype != bool:
            est = est.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})
        plot_df = plot_df.loc[~est].copy()
    if plot_df.empty:
        print(f"[WARN] No finite data for {ycol}; skip plot.")
        return

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    symbols = sorted(plot_df["symbol"].unique(), key=lambda s: Z_BY_SYMBOL[s])

    for sym in symbols:
        sub = plot_df[plot_df["symbol"] == sym].sort_values("N")
        x = sub["N"].to_numpy()
        y = sub[ycol].to_numpy(dtype=float)
        plus_col = yerr_plus_col or yerr_col
        minus_col = yerr_minus_col or yerr_col
        if plus_col and minus_col and plus_col in sub.columns and minus_col in sub.columns:
            yerr_plus = pd.to_numeric(sub[plus_col], errors="coerce").fillna(0).clip(lower=0).to_numpy(dtype=float)
            yerr_minus = pd.to_numeric(sub[minus_col], errors="coerce").fillna(0).clip(lower=0).to_numpy(dtype=float)
            if logy:
                yerr_minus = [min(em, max(0.0, yy * 0.999)) for yy, em in zip(y, yerr_minus)]
            if not any(ep > 0 for ep in yerr_plus) and not any(em > 0 for em in yerr_minus):
                ax.plot(x, y, marker="o", linewidth=1.2, markersize=4.5, label=sym)
            else:
                ax.errorbar(
                    x,
                    y,
                    yerr=[yerr_minus, yerr_plus],
                    marker="o",
                    linewidth=1.2,
                    markersize=4.5,
                    capsize=2.5,
                    label=sym,
                )
        else:
            ax.plot(x, y, marker="o", linewidth=1.2, markersize=4.5, label=sym)

    ax.set_xlabel("Neutron number N")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(title="Element", loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, fontsize=8, ncol=1)
    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    fig.savefig(outstem.with_suffix(".png"), dpi=300)
    fig.savefig(outstem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract Ackermann review ground-state decay table and plot systematics.")
    here = Path(__file__).resolve().parent
    base = here.parent
    parser.add_argument(
        "--tex",
        type=Path,
        default=base / "source" / "ppnp-LaTeX-D_ackermann_DSAS_Rev_revised_final.tex",
        help="Path to the Ackermann review LaTeX file.",
    )
    parser.add_argument("--outdir", type=Path, default=base / "out", help="Output directory.")
    parser.add_argument("--csv-name", default="ackermann_ground_states.csv", help="Output CSV filename.")
    parser.add_argument(
        "--include-estimated-qalpha",
        action="store_true",
        help="Include Q_alpha values marked by # in the static Q_alpha plot.",
    )
    args = parser.parse_args()

    tex_path = args.tex.resolve()
    outdir = args.outdir.resolve()
    if not tex_path.exists():
        raise FileNotFoundError(f"LaTeX source not found: {tex_path}")

    tex = tex_path.read_text(encoding="utf-8", errors="replace")
    table = extract_isotope_table(tex)
    records = build_raw_records(table)
    df = records_to_dataframe(records)

    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / args.csv_name
    save_csv(df, csv_path)

    plot_by_element(
        df,
        ycol="Qalpha_MeV",
        yerr_col="Qalpha_err_MeV",
        ylabel=r"$Q_\alpha$ (MeV)",
        outstem=outdir / "ackermann_gs_Qalpha_vs_N",
        logy=False,
        title=r"Ground-state $Q_\alpha$ from Ackermann DSAS review table",
        exclude_estimated_col=None if args.include_estimated_qalpha else "Qalpha_estimated",
    )
    plot_by_element(
        df,
        ycol="T_half_s",
        yerr_plus_col="T_half_err_plus_s",
        yerr_minus_col="T_half_err_minus_s",
        ylabel=r"$T_{1/2}$ (s)",
        outstem=outdir / "ackermann_gs_Ttotal_half_vs_N",
        logy=True,
        title=r"Ground-state total half-lives",
        exclude_estimated_col="T_half_estimated",
    )
    plot_by_element(
        df,
        ycol="T_alpha_half_s",
        yerr_plus_col="T_alpha_half_err_plus_s",
        yerr_minus_col="T_alpha_half_err_minus_s",
        ylabel=r"$T^{\alpha}_{1/2}$ (s)",
        outstem=outdir / "ackermann_gs_Talpha_half_vs_N",
        logy=True,
        title=r"Ground-state partial alpha half-lives",
        exclude_estimated_col="T_alpha_half_estimated",
    )
    plot_by_element(
        df,
        ycol="T_SF_half_s",
        yerr_plus_col="T_SF_half_err_plus_s",
        yerr_minus_col="T_SF_half_err_minus_s",
        ylabel=r"$T^{SF}_{1/2}$ (s)",
        outstem=outdir / "ackermann_gs_TSF_half_vs_N",
        logy=True,
        title=r"Ground-state partial spontaneous-fission half-lives",
        exclude_estimated_col="T_SF_half_estimated",
    )

    print(f"Parsed raw isotope/isomer rows: {len(records)}")
    print(f"Ground-state rows written:     {len(df)}")
    print(f"Qalpha rows marked '#':        {int(df['Qalpha_estimated'].sum()) if 'Qalpha_estimated' in df else 0}")
    print(f"CSV: {csv_path}")
    print(f"Figures: {outdir / 'ackermann_gs_Qalpha_vs_N.png'}")
    print(f"         {outdir / 'ackermann_gs_Ttotal_half_vs_N.png'}")
    print(f"         {outdir / 'ackermann_gs_Talpha_half_vs_N.png'}")
    print(f"         {outdir / 'ackermann_gs_TSF_half_vs_N.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
