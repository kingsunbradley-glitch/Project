#!/usr/bin/env python3
"""Propagate the selected uncertainties into the 205Fr transport efficiency."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "205Fr.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "205Fr_with_uncertainty.csv"
DEFAULT_ACCT_VARIANCE_INPUT = SCRIPT_DIR / "205Fr_with_ACCT_variance.csv"

ADDED_COLUMNS = [
    "Counts_stat_rel_unc_pct",
    "ACCT_avg_rel_unc_pct",
    "crossSection_rel_unc_pct",
    "TargetThickness_rel_unc_pct",
    "BranchingRatio_rel_unc_pct",
    "TranEff_total_rel_unc_pct",
    "TranEff_abs_unc_pp",
    "TranEff_lower_pct",
    "TranEff_upper_pct",
]

REQUIRED_COLUMNS = [
    "runnum",
    "Counts",
    "ACCT_avg",
    "crossSection/ub",
    "TargetThickness(mg/cm2)",
    "BranchingRatio",
    "TranEff",
]


def relative_uncertainty_pct(value: str) -> float:
    """Parse a finite relative uncertainty in the inclusive range 0--100%."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"不是有效数字: {value!r}") from exc

    if not math.isfinite(parsed) or not 0.0 <= parsed <= 100.0:
        raise argparse.ArgumentTypeError(
            f"相对误差必须是 0 到 100 之间的有限百分数: {value!r}"
        )
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "传播 Counts 泊松统计误差以及 ACCT_avg、crossSection、"
            "TargetThickness 的相对误差。"
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"输入 CSV（默认: {DEFAULT_INPUT}）",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"输出 CSV（默认: {DEFAULT_OUTPUT}）",
    )
    parser.add_argument(
        "--acct-rel-unc-pct",
        type=relative_uncertainty_pct,
        default=10.0,
        help=(
            "ACCT_avg 固定相对误差，单位 %%（默认: 10；"
            "指定 --acct-variance-input 时不使用）"
        ),
    )
    parser.add_argument(
        "--acct-variance-input",
        type=Path,
        default=None,
        help=(
            "逐 run ACCT 方差 CSV；指定后按 sqrt(ACCT_sample_variance) / "
            "abs(ACCT_raw_mean) 计算 ACCT 相对误差"
            f"（当前数据文件: {DEFAULT_ACCT_VARIANCE_INPUT}）"
        ),
    )
    parser.add_argument(
        "--cross-section-rel-unc-pct",
        type=relative_uncertainty_pct,
        default=1.2,
        help="crossSection 相对误差，单位 %%（默认: 1.2）",
    )
    parser.add_argument(
        "--target-thickness-rel-unc-pct",
        type=relative_uncertainty_pct,
        default=10.0,
        help="TargetThickness 相对误差，单位 %%（默认: 10）",
    )
    parser.add_argument(
        "--branching-ratio-abs-unc-pp",
        type=relative_uncertainty_pct,
        default=0.4,
        help=(
            "BranchingRatio 的绝对误差，单位为分支比百分数的百分点"
            "（205Fr 默认: 0.4，即 98.5(4)%%）"
        ),
    )
    return parser


def normalized_header(header: list[str]) -> tuple[list[str], dict[str, int]]:
    names = [cell.strip().lstrip("\ufeff") for cell in header]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"规范化后存在重复列名: {', '.join(duplicates)}")

    positions = {name: index for index, name in enumerate(names)}
    missing = [name for name in REQUIRED_COLUMNS if name not in positions]
    if missing:
        raise ValueError(f"输入 CSV 缺少必需列: {', '.join(missing)}")

    conflicts = [name for name in ADDED_COLUMNS if name in positions]
    if conflicts:
        raise ValueError(f"输入 CSV 已含输出列: {', '.join(conflicts)}")
    return names, positions


def parse_positive_count(value: str, row_number: int) -> float:
    try:
        count = float(value.strip())
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 Counts 不是有效数字: {value!r}") from exc

    if not math.isfinite(count) or count <= 0.0:
        raise ValueError(f"第 {row_number} 行 Counts 必须为正的有限数: {value!r}")
    return count


def parse_efficiency_pct(value: str, row_number: int) -> float:
    cleaned = value.strip()
    if cleaned.endswith("%"):
        cleaned = cleaned[:-1].strip()
    try:
        efficiency = float(cleaned)
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 TranEff 不是有效百分数: {value!r}") from exc

    if not math.isfinite(efficiency) or efficiency < 0.0:
        raise ValueError(f"第 {row_number} 行 TranEff 必须是非负有限百分数: {value!r}")
    return efficiency


def parse_branching_ratio(value: str, row_number: int) -> float:
    try:
        branching_ratio = float(value.strip())
    except ValueError as exc:
        raise ValueError(
            f"第 {row_number} 行 BranchingRatio 不是有效数字: {value!r}"
        ) from exc
    if not math.isfinite(branching_ratio) or not 0.0 < branching_ratio <= 1.0:
        raise ValueError(
            f"第 {row_number} 行 BranchingRatio 必须在 0 到 1 之间: {value!r}"
        )
    return branching_ratio


def parse_run_number(value: str, row_number: int) -> int:
    try:
        run_number = float(value.strip())
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 runnum 不是有效数字: {value!r}") from exc
    if not math.isfinite(run_number) or not run_number.is_integer():
        raise ValueError(f"第 {row_number} 行 runnum 必须是有限整数: {value!r}")
    return int(run_number)


def load_acct_relative_uncertainties(path: Path) -> dict[int, float]:
    """Return per-run ACCT relative standard deviations in percent."""
    required = ["runnum", "ACCT_raw_mean", "ACCT_sample_variance"]
    with path.open("r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.reader(input_file)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("ACCT 方差 CSV 为空") from exc

        names = [cell.strip().lstrip("\ufeff") for cell in header]
        positions = {name: index for index, name in enumerate(names)}
        missing = [name for name in required if name not in positions]
        if missing:
            raise ValueError(f"ACCT 方差 CSV 缺少必需列: {', '.join(missing)}")

        uncertainties: dict[int, float] = {}
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"ACCT 方差 CSV 第 {row_number} 行有 {len(row)} 列，"
                    f"表头有 {len(header)} 列"
                )
            run_number = parse_run_number(row[positions["runnum"]], row_number)
            if run_number in uncertainties:
                raise ValueError(f"ACCT 方差 CSV 中 run {run_number} 重复")
            try:
                mean = float(row[positions["ACCT_raw_mean"]].strip())
                variance = float(row[positions["ACCT_sample_variance"]].strip())
            except ValueError as exc:
                raise ValueError(
                    f"ACCT 方差 CSV 第 {row_number} 行均值或方差不是有效数字"
                ) from exc
            if not math.isfinite(mean) or mean == 0.0:
                raise ValueError(
                    f"ACCT 方差 CSV 第 {row_number} 行 ACCT_raw_mean 必须为非零有限数"
                )
            if not math.isfinite(variance) or variance < 0.0:
                raise ValueError(
                    f"ACCT 方差 CSV 第 {row_number} 行 ACCT_sample_variance "
                    "必须为非负有限数"
                )
            uncertainties[run_number] = 100.0 * math.sqrt(variance) / abs(mean)

    if not uncertainties:
        raise ValueError("ACCT 方差 CSV 没有数据行")
    return uncertainties


def calculate_rows(
    rows: list[list[str]],
    positions: dict[str, int],
    acct_rel_unc_pct: float,
    cross_section_rel_unc_pct: float,
    target_thickness_rel_unc_pct: float,
    branching_ratio_abs_unc_pp: float,
    acct_rel_uncertainties_by_run: dict[int, float] | None = None,
) -> tuple[list[list[str]], list[float]]:
    calculated: list[list[str]] = []
    total_relative_uncertainties: list[float] = []

    for row_number, row in enumerate(rows, start=2):
        if len(row) != len(positions):
            raise ValueError(
                f"第 {row_number} 行有 {len(row)} 列，表头有 {len(positions)} 列"
            )

        count = parse_positive_count(row[positions["Counts"]], row_number)
        efficiency_pct = parse_efficiency_pct(row[positions["TranEff"]], row_number)
        branching_ratio = parse_branching_ratio(
            row[positions["BranchingRatio"]], row_number
        )
        run_number = parse_run_number(row[positions["runnum"]], row_number)

        if acct_rel_uncertainties_by_run is None:
            row_acct_rel_unc_pct = acct_rel_unc_pct
        else:
            try:
                row_acct_rel_unc_pct = acct_rel_uncertainties_by_run[run_number]
            except KeyError as exc:
                raise ValueError(
                    f"ACCT 方差 CSV 中缺少 run {run_number} 的数据"
                ) from exc

        count_stat_rel_unc_pct = 100.0 / math.sqrt(count)
        branching_ratio_rel_unc_pct = (
            branching_ratio_abs_unc_pp / branching_ratio
        )
        total_rel_unc_pct = math.sqrt(
            count_stat_rel_unc_pct**2
            + row_acct_rel_unc_pct**2
            + cross_section_rel_unc_pct**2
            + target_thickness_rel_unc_pct**2
            + branching_ratio_rel_unc_pct**2
        )
        absolute_uncertainty_pp = efficiency_pct * total_rel_unc_pct / 100.0

        calculated.append(
            row
            + [
                f"{count_stat_rel_unc_pct:.6f}",
                f"{row_acct_rel_unc_pct:.6f}",
                f"{cross_section_rel_unc_pct:.6f}",
                f"{target_thickness_rel_unc_pct:.6f}",
                f"{branching_ratio_rel_unc_pct:.6f}",
                f"{total_rel_unc_pct:.6f}",
                f"{absolute_uncertainty_pp:.4f}",
                f"{efficiency_pct - absolute_uncertainty_pp:.4f}",
                f"{efficiency_pct + absolute_uncertainty_pp:.4f}",
            ]
        )
        total_relative_uncertainties.append(total_rel_unc_pct)

    if not calculated:
        raise ValueError("输入 CSV 没有数据行")
    return calculated, total_relative_uncertainties


def run(args: argparse.Namespace) -> tuple[int, float, float]:
    input_path = args.input.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    if input_path == output_path:
        raise ValueError("输入和输出路径不能相同，以免覆盖原始数据")

    acct_rel_uncertainties_by_run = None
    if args.acct_variance_input is not None:
        acct_variance_path = args.acct_variance_input.expanduser().resolve()
        acct_rel_uncertainties_by_run = load_acct_relative_uncertainties(
            acct_variance_path
        )

    with input_path.open("r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.reader(input_file)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("输入 CSV 为空") from exc
        _, positions = normalized_header(header)
        source_rows = list(reader)

    calculated_rows, totals = calculate_rows(
        source_rows,
        positions,
        args.acct_rel_unc_pct,
        args.cross_section_rel_unc_pct,
        args.target_thickness_rel_unc_pct,
        args.branching_ratio_abs_unc_pp,
        acct_rel_uncertainties_by_run,
    )

    with output_path.open("w", encoding="utf-8-sig", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(header + ADDED_COLUMNS)
        writer.writerows(calculated_rows)

    return len(calculated_rows), min(totals), max(totals)


def main() -> int:
    args = build_parser().parse_args()
    try:
        row_count, minimum, maximum = run(args)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    print(
        f"已写入 {args.output}: {row_count} 行；"
        f"总相对误差范围 {minimum:.6f}%--{maximum:.6f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
