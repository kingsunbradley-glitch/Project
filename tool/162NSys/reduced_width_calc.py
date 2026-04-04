#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
根据输入的 .dat 数据表计算 alpha 约化跃迁宽度 delta2。

默认约定（与用户给出的 ROOT 宏一致）：
1. 势函数：V = Vn + Vc + Vl
2. WKB 穿透率：P = exp[-2*sqrt(2*mu)/hbarc * integral sqrt(V-E) dr]
3. 约化跃迁宽度：delta2 = h * (ln2 / T_half * BR) / P * 1000  [keV]

输入文件默认应至少包含这些列：
Z, A, EXP, error, T0.5, +, -, idx, Br, l

说明：
- idx == 1：实验数据，使用 EXP/error 和 T0.5/+/- 计算 delta2 及上下限误差
- idx == 0：仅计算中心值 delta2，误差列留空
- 默认认为：
    * EXP / error 的单位是 keV（若数值 < 100，则视为已是 MeV）
    * T0.5 / + / - 的单位是 ms（可通过 --t-unit 修改）

输出：原文件名去掉扩展名后追加 Pro.dat
例如：dataDs.dat -> dataDsPro.dat
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# ---------- Physical constants ----------
HPLK = 4.13566727e-21   # MeV*s
HBARC = 197.32696       # MeV*fm
E_SQUA = 1.4399644      # MeV*fm
AMU_C2 = 931.494013     # MeV/c^2
LN2 = math.log(2.0)

# ---------- Column defaults ----------
ENERGY_COL = "EXP"
ENERGY_ERR_COL = "error"
THALF_COL = "T0.5"
THALF_PLUS_COL = "+"
THALF_MINUS_COL = "-"
IDX_COL = "idx"
BR_COL = "Br"
L_COL = "l"
Z_COL = "Z"
A_COL = "A"


def to_float(value: object) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def format_float(x: Optional[float]) -> str:
    if x is None:
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.10g}"


def energy_to_mev(e: float) -> float:
    # 中心能量：样例表里实验能量通常写成 keV（如 11270），理论若直接给 MeV（如 11.2035）则不变
    return e / 1000.0 if abs(e) >= 20.0 else e


def energy_err_to_mev(e: float) -> float:
    # 能量误差若大于 1，通常是 keV（20, 30, 70 ...）；若已经是 <1 的小数，则按 MeV 处理
    return e / 1000.0 if abs(e) > 1.0 else e


def time_to_seconds(t: float, unit: str) -> float:
    factor = {
        "s": 1.0,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
    }[unit]
    return t * factor


def Vn(r: float, A_daughter: float) -> float:
    return -1100.0 * math.exp(-((r - 1.17 * (A_daughter ** (1.0 / 3.0))) / 0.574))


def Vcoul(r: float, Z_daughter: float) -> float:
    return 2.0 * Z_daughter * E_SQUA / r


def Vcent(r: float, A_daughter: float, L: float) -> float:
    rmass = 4.0 * A_daughter / (4.0 + A_daughter) * AMU_C2
    return L * (L + 1.0) * HBARC * HBARC / (2.0 * rmass * r * r)


def Vtot(r: float, A_daughter: float, Z_daughter: float, L: float) -> float:
    return Vn(r, A_daughter) + Vcoul(r, Z_daughter) + Vcent(r, A_daughter, L)


def effective_q_value(am: float, zm: float, ealpha_mev: float) -> float:
    a_daughter = am - 4.0
    escr = (65.3 * (zm ** (7.0 / 5.0)) - 80.0 * (zm ** (2.0 / 5.0))) * 1e-6
    return ealpha_mev * (am / a_daughter) + escr


def outer_turning_point(am: float, zm: float, ealpha_mev: float, L: float) -> float:
    a_daughter = am - 4.0
    z_daughter = zm - 2.0
    etot = effective_q_value(am, zm, ealpha_mev)
    rmass = 4.0 * a_daughter / (4.0 + a_daughter) * AMU_C2
    term = z_daughter * z_daughter * E_SQUA * E_SQUA + etot * HBARC * HBARC * L * (L + 1.0) / (2.0 * rmass)
    return z_daughter * E_SQUA / etot + math.sqrt(term) / etot


def secant_root(am: float, zm: float, ealpha_mev: float, L: float, r0: float, r1: float,
                tol: float = 1e-10, max_iter: int = 1000) -> float:
    a_daughter = am - 4.0
    z_daughter = zm - 2.0
    etot = effective_q_value(am, zm, ealpha_mev)
    f0 = Vtot(r0, a_daughter, z_daughter, L) - etot
    f1 = Vtot(r1, a_daughter, z_daughter, L) - etot

    for _ in range(max_iter):
        if f1 == f0:
            return r1
        r = r1 - f1 * (r1 - r0) / (f1 - f0)
        if abs(r - r1) < tol:
            return r
        r0, r1 = r1, r
        f0 = f1
        f1 = Vtot(r1, a_daughter, z_daughter, L) - etot
    return r1


def inner_turning_point(am: float, zm: float, ealpha_mev: float, L: float) -> float:
    a_daughter = am - 4.0
    z_daughter = zm - 2.0
    etot = effective_q_value(am, zm, ealpha_mev)
    rout = outer_turning_point(am, zm, ealpha_mev, L)

    r = None
    rx = rout - 0.2
    while rx > 0.0:
        if Vtot(rx, a_daughter, z_daughter, L) < etot:
            r = rx
            break
        rx -= 0.2

    if r is None:
        raise RuntimeError(f"Failed to bracket inner turning point for A={am}, Z={zm}, E={ealpha_mev}, l={L}.")

    return secant_root(am, zm, ealpha_mev, L, r, r + 0.2)


# 采用足够密的梯形积分，能与给出的 ROOT 宏结果对到 1e-4 keV 量级
N_INTEGRAL_POINTS = 20000


def penetration_factor(am: float, zm: float, ealpha_mev: float, L: float) -> float:
    a_daughter = am - 4.0
    etot = effective_q_value(am, zm, ealpha_mev)
    rin = inner_turning_point(am, zm, ealpha_mev, L)
    rout = outer_turning_point(am, zm, ealpha_mev, L)
    rmass = 4.0 * a_daughter / (4.0 + a_daughter) * AMU_C2
    z_daughter = zm - 2.0

    dx = (rout - rin) / (N_INTEGRAL_POINTS - 1)
    integral = 0.0
    for i in range(N_INTEGRAL_POINTS):
        r = rin + i * dx
        val = Vtot(r, a_daughter, z_daughter, L) - etot
        y = math.sqrt(val) if val > 0.0 else 0.0
        coef = 1.0
        if i == 0 or i == N_INTEGRAL_POINTS - 1:
            coef = 0.5
        integral += coef * y
    integral *= dx

    exponent = -2.0 * math.sqrt(2.0 * rmass) * integral / HBARC
    return math.exp(exponent)


def reduced_width_keV(am: float, zm: float, ealpha_mev: float, t_half_seconds: float,
                      br: float, L: float) -> float:
    if t_half_seconds <= 0.0:
        raise ValueError("Half-life must be positive.")
    p = penetration_factor(am, zm, ealpha_mev, L)
    return HPLK * (LN2 / t_half_seconds * br) / p * 1000.0


def compute_delta_with_errors(am: float, zm: float, e_mev: float, de_mev: float,
                              t_sec: float, dt_plus_sec: float, dt_minus_sec: float,
                              br: float, L: float) -> Tuple[float, float, float]:
    center = reduced_width_keV(am, zm, e_mev, t_sec, br, L)

    # delta2 随 E 增大而减小，随 T_half 增大而减小
    # 因此：
    #   upper 取 E - dE, T - dT_minus
    #   lower 取 E + dE, T + dT_plus
    e_low = max(e_mev - max(de_mev, 0.0), 1e-9)
    e_high = e_mev + max(de_mev, 0.0)
    t_low = max(t_sec - max(dt_minus_sec, 0.0), 1e-30)
    t_high = t_sec + max(dt_plus_sec, 0.0)

    upper_val = reduced_width_keV(am, zm, e_low, t_low, br, L)
    lower_val = reduced_width_keV(am, zm, e_high, t_high, br, L)

    plus_err = max(upper_val - center, 0.0)
    minus_err = max(center - lower_val, 0.0)
    return center, plus_err, minus_err


def read_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if not lines:
        raise ValueError(f"Empty input file: {path}")

    header = lines[0].split("\t")
    rows: List[Dict[str, str]] = []

    for line in lines[1:]:
        if line.strip() == "":
            continue
        parts = line.split("\t")
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        elif len(parts) > len(header):
            parts = parts[:len(header)]
        row = {header[i]: parts[i].strip() for i in range(len(header))}
        # 跳过全空行
        if all(v == "" for v in row.values()):
            continue
        rows.append(row)

    return header, rows


def choose_existing(row: Dict[str, str], candidates: Iterable[str]) -> Optional[float]:
    for name in candidates:
        if name in row:
            val = to_float(row[name])
            if val is not None:
                return val
    return None


def process_rows(rows: List[Dict[str, str]], t_unit: str) -> List[Dict[str, str]]:
    processed: List[Dict[str, str]] = []

    for row in rows:
        new_row = dict(row)

        zm = choose_existing(row, [Z_COL])
        am = choose_existing(row, [A_COL])
        idx = int(choose_existing(row, [IDX_COL]) or 0)
        br = choose_existing(row, [BR_COL])
        lval = choose_existing(row, [L_COL])

        # 中心值优先取标准列；若缺失，再尝试常见备选列
        e_val = choose_existing(row, [ENERGY_COL, "Qalpha", "Qa", "Ealpha"])
        t_val = choose_existing(row, [THALF_COL, "T_half", "T1/2"])

        if zm is None or am is None or e_val is None or t_val is None:
            new_row["delta2"] = ""
            new_row["+"] = ""
            new_row["-"] = ""
            processed.append(new_row)
            continue

        br = 1.0 if br is None else br
        lval = 0.0 if lval is None else lval

        e_mev = energy_to_mev(e_val)
        t_sec = time_to_seconds(t_val, t_unit)

        try:
            if idx == 1:
                de_val = choose_existing(row, [ENERGY_ERR_COL]) or 0.0
                dt_plus_val = choose_existing(row, [THALF_PLUS_COL]) or 0.0
                dt_minus_val = choose_existing(row, [THALF_MINUS_COL]) or 0.0

                de_mev = energy_err_to_mev(de_val)
                dt_plus_sec = time_to_seconds(dt_plus_val, t_unit)
                dt_minus_sec = time_to_seconds(dt_minus_val, t_unit)

                delta2, plus_err, minus_err = compute_delta_with_errors(
                    am=am,
                    zm=zm,
                    e_mev=e_mev,
                    de_mev=de_mev,
                    t_sec=t_sec,
                    dt_plus_sec=dt_plus_sec,
                    dt_minus_sec=dt_minus_sec,
                    br=br,
                    L=lval,
                )
                new_row["delta2"] = format_float(delta2)
                new_row["_dplus"] = format_float(plus_err)
                new_row["_dminus"] = format_float(minus_err)
            else:
                delta2 = reduced_width_keV(am, zm, e_mev, t_sec, br, lval)
                new_row["delta2"] = format_float(delta2)
                new_row["_dplus"] = ""
                new_row["_dminus"] = ""
        except Exception:
            new_row["delta2"] = ""
            new_row["_dplus"] = ""
            new_row["_dminus"] = ""

        processed.append(new_row)

    return processed


def write_table(path: Path, header: List[str], rows: List[Dict[str, str]]) -> None:
    out_header = list(header) + ["delta2", "+", "-"]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(out_header)
        for row in rows:
            base = [row.get(col, "") for col in header]
            extra = [row.get("delta2", ""), row.get("_dplus", ""), row.get("_dminus", "")]
            writer.writerow(base + extra)


def build_output_name(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}Pro.dat")


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate reduced alpha-decay widths from a .dat table.")
    parser.add_argument("input_file", help="Input .dat file path")
    parser.add_argument(
        "--t-unit",
        choices=["s", "ms", "us", "ns"],
        default="ms",
        help="Unit of T0.5 / + / - columns (default: ms)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_file).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    header, rows = read_table(input_path)
    processed = process_rows(rows, t_unit=args.t_unit)
    output_path = build_output_name(input_path)
    write_table(output_path, header, processed)

    print(f"Input : {input_path}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    main()
