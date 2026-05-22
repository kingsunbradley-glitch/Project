#!/usr/bin/env python3
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


INPUT_COLUMNS = [
    "E_ER",
    "X_Pos",
    "Y_Pos",
    "Delta_t1",
    "E_1_DSSD",
    "E_1_SSD",
]

SF_COLUMNS = [
    "chain_no",
    "Elab",
    "file_no",
    "x_strip",
    "y_strip",
    "ER_E",
    "SF_E",
    "DSSD_E1",
    "SSD_E1",
    "DeltaT_ms",
]


def read_input_header(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                return [field.strip() for field in line.split(",")]
    raise ValueError(f"{path} is empty.")


def read_sf_rows(path: Path) -> list[dict[str, Decimal]]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            values = line.split()
            if len(values) != len(SF_COLUMNS):
                raise ValueError(
                    f"{path}:{line_no}: expected {len(SF_COLUMNS)} columns, got {len(values)}"
                )
            rows.append(dict(zip(SF_COLUMNS, (Decimal(value) for value in values))))
    return rows


def round_int(value: Decimal) -> str:
    return str(int(value.quantize(Decimal("1"), rounding=ROUND_HALF_UP)))


def mev_to_keV(value: Decimal) -> str:
    return round_int(value * Decimal("1000"))


def ms_to_us(value: Decimal) -> str:
    text = format(value * Decimal("1000"), "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def sf_to_input_row(row: dict[str, Decimal]) -> dict[str, str]:
    return {
        "E_ER": mev_to_keV(row["ER_E"]),
        "X_Pos": round_int(row["x_strip"]),
        "Y_Pos": round_int(row["y_strip"]),
        "Delta_t1": ms_to_us(row["DeltaT_ms"]),
        "E_1_DSSD": mev_to_keV(row["DSSD_E1"]),
        "E_1_SSD": mev_to_keV(row["SSD_E1"]),
    }


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Convert table/SF.dat into the input.txt format used by generate_fission_tikz.py."
    )
    parser.add_argument("--sf", type=Path, default=here / "table" / "SF.dat", help="Source SF.dat file.")
    parser.add_argument("--input", type=Path, default=here / "input.txt", help="Existing input.txt with header.")
    parser.add_argument("--output", type=Path, default=None, help="Output path. Defaults to overwriting input.txt.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = (args.output or args.input).resolve()
    sf_path = args.sf.resolve()

    header = read_input_header(input_path)
    missing = [column for column in INPUT_COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{input_path} is missing required columns: {', '.join(missing)}")

    rows = [sf_to_input_row(row) for row in read_sf_rows(sf_path)]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(row[column] for column in header) + "\n")

    print(f"Wrote {output_path} from {sf_path} with {len(rows)} row(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
