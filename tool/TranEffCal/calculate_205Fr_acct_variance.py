#!/usr/bin/env python3
"""Calculate per-run ACCT statistics from the timestamped raw CSV files."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "205Fr.csv"
DEFAULT_ACCT_DIR = SCRIPT_DIR / "ACCT"
DEFAULT_OUTPUT = SCRIPT_DIR / "205Fr_with_ACCT_variance.csv"

REQUIRED_COLUMNS = ["runnum", "date", "Tsatrt", "Tsustain/s", "ACCT_avg"]
ADDED_COLUMNS = [
    "ACCT_window_start",
    "ACCT_window_end_exclusive",
    "ACCT_sample_count",
    "ACCT_raw_mean",
    "ACCT_sample_variance",
    "ACCT_sample_std",
    "ACCT_rel_std_pct",
    "ACCT_mean_sem",
    "ACCT_rel_sem_pct",
    "ACCT_avg_minus_raw_mean",
]


@dataclass
class RunningStatistics:
    count: int = 0
    mean: float = 0.0
    sum_squared_deviations: float = 0.0
    minimum: float = math.inf
    maximum: float = -math.inf

    def add(self, value: float) -> None:
        """Update mean and variance accumulators using Welford's algorithm."""
        self.count += 1
        delta = value - self.mean
        self.mean += delta / self.count
        self.sum_squared_deviations += delta * (value - self.mean)
        self.minimum = min(self.minimum, value)
        self.maximum = max(self.maximum, value)

    @property
    def sample_variance(self) -> float:
        if self.count < 2:
            raise ValueError("计算样本方差至少需要 2 个 ACCT 数据点")
        return self.sum_squared_deviations / (self.count - 1)

    @property
    def sample_standard_deviation(self) -> float:
        return math.sqrt(self.sample_variance)

    @property
    def standard_error_of_mean(self) -> float:
        return self.sample_standard_deviation / math.sqrt(self.count)


@dataclass
class RunWindow:
    run_number: int
    start: datetime
    end: datetime
    listed_acct_average: float
    source_row: list[str]
    statistics: RunningStatistics = field(default_factory=RunningStatistics)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="按 205Fr.csv 的日期、起始时刻和持续时间计算逐 run ACCT 方差。"
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT, help=f"输入 CSV（默认: {DEFAULT_INPUT}）"
    )
    parser.add_argument(
        "--acct-dir",
        type=Path,
        default=DEFAULT_ACCT_DIR,
        help=f"原始 ACCT CSV 目录（默认: {DEFAULT_ACCT_DIR}）",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"输出 CSV（默认: {DEFAULT_OUTPUT}）",
    )
    parser.add_argument(
        "--year",
        type=int,
        default=2026,
        help="当 date 列只有月和日时采用的年份（默认: 2026）",
    )
    return parser


def normalized_positions(header: list[str]) -> dict[str, int]:
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
    return positions


def parse_finite_float(value: str, label: str, row_number: int) -> float:
    try:
        parsed = float(value.strip())
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 {label} 不是有效数字: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"第 {row_number} 行 {label} 必须是有限数: {value!r}")
    return parsed


def parse_date(value: str, default_year: int, row_number: int) -> tuple[int, int, int]:
    parts = value.strip().split(".")
    try:
        if len(parts) == 2:
            month, day = map(int, parts)
            return default_year, month, day
        if len(parts) == 3:
            year, month, day = map(int, parts)
            return year, month, day
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 date 格式无效: {value!r}") from exc
    raise ValueError(f"第 {row_number} 行 date 应为 M.D 或 YYYY.M.D: {value!r}")


def parse_start_time(
    date_value: str, time_value: str, default_year: int, row_number: int
) -> datetime:
    year, month, day = parse_date(date_value, default_year, row_number)
    time_parts = time_value.strip().split(":")
    try:
        if len(time_parts) == 2:
            hour, minute = map(int, time_parts)
            second = 0
        elif len(time_parts) == 3:
            hour, minute, second = map(int, time_parts)
        else:
            raise ValueError
        return datetime(year, month, day, hour, minute, second)
    except ValueError as exc:
        raise ValueError(f"第 {row_number} 行 Tsatrt 格式无效: {time_value!r}") from exc


def read_run_windows(
    input_path: Path, default_year: int
) -> tuple[list[str], list[RunWindow]]:
    with input_path.open("r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.reader(input_file)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("输入 CSV 为空") from exc
        positions = normalized_positions(header)

        runs: list[RunWindow] = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"第 {row_number} 行有 {len(row)} 列，表头有 {len(header)} 列"
                )
            run_value = parse_finite_float(row[positions["runnum"]], "runnum", row_number)
            if not run_value.is_integer():
                raise ValueError(f"第 {row_number} 行 runnum 必须为整数: {run_value}")
            duration = parse_finite_float(
                row[positions["Tsustain/s"]], "Tsustain/s", row_number
            )
            if duration <= 0.0:
                raise ValueError(f"第 {row_number} 行 Tsustain/s 必须为正数")
            listed_average = parse_finite_float(
                row[positions["ACCT_avg"]], "ACCT_avg", row_number
            )
            start = parse_start_time(
                row[positions["date"]],
                row[positions["Tsatrt"]],
                default_year,
                row_number,
            )
            runs.append(
                RunWindow(
                    run_number=int(run_value),
                    start=start,
                    end=start + timedelta(seconds=duration),
                    listed_acct_average=listed_average,
                    source_row=row,
                )
            )

    if not runs:
        raise ValueError("输入 CSV 没有数据行")

    runs.sort(key=lambda run: run.start)
    for previous, current in zip(runs, runs[1:]):
        if current.start < previous.end:
            raise ValueError(
                f"run {previous.run_number} 与 run {current.run_number} 的时间窗重叠"
            )
    return header, runs


def collect_acct_statistics(acct_dir: Path, runs: list[RunWindow]) -> None:
    acct_files = sorted(acct_dir.glob("*.csv"))
    if not acct_files:
        raise ValueError(f"ACCT 目录中没有 CSV 文件: {acct_dir}")

    starts = [run.start for run in runs]
    earliest = runs[0].start
    latest = runs[-1].end
    relevant_dates = {
        (earliest + timedelta(days=offset)).strftime("%Y/%m/%d")
        for offset in range((latest.date() - earliest.date()).days + 1)
    }

    for acct_path in acct_files:
        with acct_path.open("r", encoding="utf-8", newline="") as acct_file:
            for line_number, row in enumerate(csv.reader(acct_file), start=1):
                if len(row) < 3:
                    raise ValueError(
                        f"{acct_path.name} 第 {line_number} 行不足 3 列"
                    )

                timestamp_text = row[1].strip()
                if timestamp_text[:10] not in relevant_dates:
                    continue
                try:
                    timestamp = datetime.strptime(
                        timestamp_text, "%Y/%m/%d-%H:%M:%S.%f"
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"{acct_path.name} 第 {line_number} 行时间戳无效: {timestamp_text!r}"
                    ) from exc
                if timestamp < earliest or timestamp >= latest:
                    continue

                run_index = bisect.bisect_right(starts, timestamp) - 1
                if run_index < 0 or timestamp >= runs[run_index].end:
                    continue

                try:
                    current = float(row[2].strip())
                except ValueError as exc:
                    raise ValueError(
                        f"{acct_path.name} 第 {line_number} 行 ACCT 值无效: {row[2]!r}"
                    ) from exc
                if not math.isfinite(current):
                    raise ValueError(
                        f"{acct_path.name} 第 {line_number} 行 ACCT 值不是有限数"
                    )
                runs[run_index].statistics.add(current)

    missing = [run.run_number for run in runs if run.statistics.count < 2]
    if missing:
        raise ValueError(f"以下 run 缺少足够的 ACCT 数据点: {missing}")


def write_results(output_path: Path, header: list[str], runs: list[RunWindow]) -> None:
    with output_path.open("w", encoding="utf-8-sig", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(header + ADDED_COLUMNS)
        for run in runs:
            stats = run.statistics
            standard_deviation = stats.sample_standard_deviation
            standard_error = stats.standard_error_of_mean
            if stats.mean == 0.0:
                relative_standard_deviation = math.nan
                relative_standard_error = math.nan
            else:
                relative_standard_deviation = 100.0 * standard_deviation / abs(stats.mean)
                relative_standard_error = 100.0 * standard_error / abs(stats.mean)

            writer.writerow(
                run.source_row
                + [
                    run.start.strftime("%Y-%m-%d %H:%M:%S"),
                    run.end.strftime("%Y-%m-%d %H:%M:%S"),
                    str(stats.count),
                    f"{stats.mean:.10f}",
                    f"{stats.sample_variance:.10f}",
                    f"{standard_deviation:.10f}",
                    f"{relative_standard_deviation:.6f}",
                    f"{standard_error:.10f}",
                    f"{relative_standard_error:.6f}",
                    f"{run.listed_acct_average - stats.mean:.10f}",
                ]
            )


def run(args: argparse.Namespace) -> tuple[list[RunWindow], list[int]]:
    input_path = args.input.expanduser().resolve()
    acct_dir = args.acct_dir.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    if input_path == output_path:
        raise ValueError("输入和输出路径不能相同，以免覆盖原始数据")

    header, runs = read_run_windows(input_path, args.year)
    collect_acct_statistics(acct_dir, runs)
    write_results(output_path, header, runs)

    mismatches = [
        item.run_number
        for item in runs
        if not math.isclose(
            item.listed_acct_average,
            item.statistics.mean,
            rel_tol=0.0,
            abs_tol=5e-7,
        )
    ]
    return runs, mismatches


def main() -> int:
    args = build_parser().parse_args()
    try:
        runs, mismatches = run(args)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    variances = [item.statistics.sample_variance for item in runs]
    print(
        f"已写入 {args.output}: {len(runs)} 个 run；"
        f"ACCT 样本方差范围 {min(variances):.10f}--{max(variances):.10f}"
    )
    if mismatches:
        print(
            "注意：以下 run 的原表 ACCT_avg 与完整时间窗原始均值不一致："
            + ", ".join(map(str, mismatches))
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
