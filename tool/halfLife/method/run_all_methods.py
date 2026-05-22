#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run every half-life method script in the subdirectories with the shared
input.csv, then collect all results under output/.

The merged CSV keeps only:
    name, method, Texp, Tcal, HF
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "input.csv"
DEFAULT_OUTPUT_DIR = BASE_DIR / "output"
MERGED_FILENAME = "all_methods_combined.csv"
SUMMARY_FIELDS = ["name", "method", "Texp", "Tcal", "HF"]


def find_method_scripts(base_dir: Path) -> list[tuple[str, Path]]:
    scripts: list[tuple[str, Path]] = []
    for method_dir in sorted(p for p in base_dir.iterdir() if p.is_dir()):
        if method_dir.name == "output" or method_dir.name.startswith("."):
            continue
        for script in sorted(method_dir.glob("*.py")):
            scripts.append((method_dir.name, script))
    return scripts


def run_script(method: str, script: Path, input_csv: Path, output_dir: Path) -> Path:
    output_path = output_dir / f"{method}_{script.stem}.csv"

    with tempfile.TemporaryDirectory(prefix=f"{method}_", dir=output_dir) as tmp:
        work_dir = Path(tmp)
        work_input = work_dir / "input.csv"
        shutil.copy2(input_csv, work_input)

        cmd = [
            sys.executable,
            str(script),
            str(work_input),
            "-o",
            str(output_path),
        ]
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"{method}/{script.name} failed with exit code {result.returncode}\n"
                f"Command: {' '.join(cmd)}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )

        # Some legacy scripts ignore command-line arguments and always write a
        # fixed output filename in cwd. Move that file into output/.
        if not output_path.exists():
            generated = sorted(
                p for p in work_dir.glob("*.csv") if p.name != work_input.name
            )
            if not generated:
                raise RuntimeError(
                    f"{method}/{script.name} finished but did not create a CSV output.\n"
                    f"stdout:\n{result.stdout}\n"
                    f"stderr:\n{result.stderr}"
                )
            shutil.copy2(generated[-1], output_path)

    return output_path


def unique_headers(headers: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for header in headers:
        count = seen.get(header, 0)
        out.append(header if count == 0 else f"{header}.{count}")
        seen[header] = count + 1
    return out


def read_csv_rows(csv_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        try:
            raw_headers = next(reader)
        except StopIteration as exc:
            raise RuntimeError(f"Empty CSV output: {csv_path}") from exc

        headers = unique_headers(raw_headers)
        rows = [dict(zip(headers, row)) for row in reader]
    return headers, rows


def first_value(row: dict[str, str], columns: list[str]) -> str:
    for col in columns:
        if col in row and row[col] != "":
            return row[col]
    return ""


def method_value(row: dict[str, str], method: str, field: str) -> str:
    method_map = {
        "DZR": {
            "Tcal": ["T_DZR_s"],
            "HF": ["HF_DZR"],
        },
        "P2": {
            "Tcal": ["T_calc_s"],
            "HF": ["HF"],
        },
        "Poenaru": {
            "Tcal": ["Tcalc_s"],
            "HF": ["HF"],
        },
        "Viola-Seaborg": {
            "Tcal": ["T_VS_s"],
            "HF": ["HF_VS"],
        },
    }
    fallback_map = {
        "Tcal": ["Tcalc_s", "T_calc_s", "T_DZR_s", "T_VS_s", "T_Royer_s"],
        "HF": ["HF", "HF_DZR", "HF_VS", "HF_Royer"],
    }
    columns = method_map.get(method, {}).get(field, []) + fallback_map[field]
    return first_value(row, columns)


def merge_outputs(results: list[tuple[str, str, Path]], merged_path: Path) -> None:
    merged_rows: list[dict[str, str]] = []

    for method, script_name, csv_path in results:
        _, rows = read_csv_rows(csv_path)
        for row in rows:
            merged_rows.append(
                {
                    "name": first_value(row, ["Name", "name"]),
                    "method": method,
                    "Texp": first_value(row, ["Talpha_exp_s", "Texp_s", "T_exp_s"]),
                    "Tcal": method_value(row, method, "Tcal"),
                    "HF": method_value(row, method, "HF"),
                }
            )

    if not merged_rows:
        raise RuntimeError("No method outputs were generated.")

    with merged_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run all method/*.py scripts with method/input.csv and merge outputs."
    )
    parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_INPUT,
        type=Path,
        help=f"Shared input CSV. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        type=Path,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--merged-name",
        default=MERGED_FILENAME,
        help=f"Merged output filename. Default: {MERGED_FILENAME}",
    )
    args = parser.parse_args()

    input_csv = args.input.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    scripts = find_method_scripts(BASE_DIR)
    if not scripts:
        raise RuntimeError(f"No method scripts found under: {BASE_DIR}")

    results: list[tuple[str, str, Path]] = []
    for method, script in scripts:
        print(f"Running {method}/{script.name} ...")
        output_path = run_script(method, script.resolve(), input_csv, output_dir)
        results.append((method, script.name, output_path))
        print(f"  wrote {output_path}")

    merged_path = output_dir / args.merged_name
    merge_outputs(results, merged_path)
    print(f"Merged output written to: {merged_path}")


if __name__ == "__main__":
    main()
