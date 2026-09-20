#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run all half-life method scripts in this folder with input.csv and merge the
results into main.csv.

Merged columns:
    name, method, Ealpha, Texp, Tcal, HF
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "input.csv"
DEFAULT_OUTPUT_DIR = BASE_DIR / "output"
DEFAULT_MAIN = BASE_DIR / "main.csv"
SUMMARY_FIELDS = ["name", "method", "Ealpha", "Texp", "Tcal", "HF"]

METHOD_NAMES = {
    "PoenaruMethod.py": "Poenaru",
    "calc_dzr_poenaru_format.py": "DZR",
    "calc_ismail2022_formulaE.py": "Ismail2022_FormulaE",
    "calc_qi2009_udl.py": "Qi2009_UDL",
    "calc_viola_seaborg.py": "Viola-Seaborg",
    "calc_xu2022_unified.py": "Xu2022_Unified",
    "calc_poenaru_E. Rurarz_format.py": "E. Rurarz",
}

METHOD_COLUMNS = {
    "DZR": {
        "Tcal": ["T_DZR_s"],
        "HF": ["HF_DZR"],
    },
    "Poenaru": {
        "Tcal": ["Tcalc_s"],
        "HF": ["HF"],
    },
    "Ismail2022_FormulaE": {
        "Tcal": ["T_Ismail2022_FE_s"],
        "HF": ["HF_Ismail2022_FE"],
    },
    "Qi2009_UDL": {
        "Tcal": ["T_Qi2009_UDL_s"],
        "HF": ["HF_Qi2009_UDL"],
    },
    "Viola-Seaborg": {
        "Tcal": ["T_VS_s"],
        "HF": ["HF_VS"],
    },
    "Xu2022_Unified": {
        "Tcal": ["T_Xu2022_s"],
        "HF": ["HF_Xu2022"],
    },
    "E. Rurarz": {
        "Tcal": ["T12a_Poenaru_s"],
        "HF": ["HF_Poenaru"],
    },
}

COMMON_COLUMNS = {
    "name": ["Name", "name"],
    "Ealpha": ["Ealpha_MeV", "E_alpha_MeV", "Ealpha", "E_alpha", "EXP", "EXP_MeV"],
    "Texp": ["Texp_s", "Talpha_exp_s", "T_exp_s", "T0.5", "T12_s", "T1/2_s"],
    "Tcal": [
        "Tcalc_s",
        "T_calc_s",
        "T_DZR_s",
        "T_Ismail2022_FE_s",
        "T_Qi2009_UDL_s",
        "T_VS_s",
        "T_Xu2022_s",
        "T_Royer_s",
        "T12a_Poenaru_s",
    ],
    "HF": [
        "HF",
        "HF_DZR",
        "HF_Ismail2022_FE",
        "HF_Qi2009_UDL",
        "HF_VS",
        "HF_Xu2022",
        "HF_Royer",
        "HF_Poenaru",
    ],
}


def method_name(script: Path) -> str:
    if script.name in METHOD_NAMES:
        return METHOD_NAMES[script.name]
    name = script.stem
    if name.startswith("calc_"):
        name = name[len("calc_") :]
    return name


def find_scripts(base_dir: Path, self_path: Path) -> list[Path]:
    scripts = []
    for script in sorted(base_dir.glob("*.py")):
        if script.resolve() == self_path.resolve():
            continue
        if script.name.startswith("_"):
            continue
        scripts.append(script)
    return scripts


def output_name(method: str) -> str:
    safe = method.replace("/", "_").replace(" ", "_")
    return f"{safe}.csv"


def run_script(script: Path, method: str, input_csv: Path, output_dir: Path) -> Path:
    output_csv = output_dir / output_name(method)
    cmd = [sys.executable, str(script), str(input_csv), "-o", str(output_csv)]
    result = subprocess.run(
        cmd,
        cwd=BASE_DIR,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{script.name} failed with exit code {result.returncode}\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    if not output_csv.exists():
        raise RuntimeError(f"{script.name} finished but did not create {output_csv}")
    return output_csv


def unique_headers(headers: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for header in headers:
        count = seen.get(header, 0)
        out.append(header if count == 0 else f"{header}.{count}")
        seen[header] = count + 1
    return out


def read_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        try:
            headers = unique_headers(next(reader))
        except StopIteration as exc:
            raise RuntimeError(f"Empty CSV: {csv_path}") from exc
        return [dict(zip(headers, row)) for row in reader]


def first_value(row: dict[str, str], columns: list[str]) -> str:
    for col in columns:
        value = row.get(col, "")
        if value != "":
            return value
    return ""


def summary_value(row: dict[str, str], method: str, field: str) -> str:
    method_cols = METHOD_COLUMNS.get(method, {}).get(field, [])
    common_cols = COMMON_COLUMNS[field]
    return first_value(row, method_cols + common_cols)


def merge_outputs(results: list[tuple[str, Path]], main_csv: Path) -> None:
    rows = []
    for method, output_csv in results:
        for row in read_rows(output_csv):
            rows.append(
                {
                    "name": summary_value(row, method, "name"),
                    "method": method,
                    "Ealpha": summary_value(row, method, "Ealpha"),
                    "Texp": summary_value(row, method, "Texp"),
                    "Tcal": summary_value(row, method, "Tcal"),
                    "HF": summary_value(row, method, "HF"),
                }
            )

    if not rows:
        raise RuntimeError("No rows were generated.")

    with main_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run every method script in this folder and merge to main.csv."
    )
    parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_INPUT,
        type=Path,
        help=f"Input CSV. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        type=Path,
        help=f"Directory for each method output CSV. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "-m",
        "--main",
        default=DEFAULT_MAIN,
        type=Path,
        help=f"Merged CSV. Default: {DEFAULT_MAIN}",
    )
    args = parser.parse_args()

    input_csv = args.input.resolve()
    output_dir = args.output_dir.resolve()
    main_csv = args.main.resolve()

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    output_dir.mkdir(parents=True, exist_ok=True)
    scripts = find_scripts(BASE_DIR, Path(__file__))
    if not scripts:
        raise RuntimeError(f"No method scripts found in {BASE_DIR}")

    results = []
    for script in scripts:
        method = method_name(script)
        print(f"Running {script.name} -> {method}")
        output_csv = run_script(script.resolve(), method, input_csv, output_dir)
        results.append((method, output_csv))
        print(f"  wrote {output_csv}")

    merge_outputs(results, main_csv)
    print(f"Merged CSV written to: {main_csv}")


if __name__ == "__main__":
    main()
