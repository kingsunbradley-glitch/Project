#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path


REQUIRED_COLUMNS = [
    "E_ER",
    "X_Pos",
    "Y_Pos",
    "Delta_t1",
    "E_1_DSSD",
    "E_1_SSD",
]
EVENTS_PER_ROW = 3


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def split_fields(line: str) -> list[str]:
    if "," in line:
        return [field.strip() for field in next(csv.reader([line]))]
    return re.split(r"\s+", line.strip())


def read_input(path: Path) -> list[dict[str, str]]:
    lines = [strip_comment(line) for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        raise ValueError(f"{path} is empty.")

    header = split_fields(lines[0])
    missing = [col for col in REQUIRED_COLUMNS if col not in header]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    rows = []
    for line in lines[1:]:
        values = split_fields(line)
        if len(values) != len(header):
            raise ValueError(f"Column count mismatch in line: {line}")
        rows.append(dict(zip(header, values)))
    return rows


def decimal_places(text: str) -> int:
    if "e" in text.lower():
        value = Decimal(text)
        return max(0, -value.as_tuple().exponent)
    if "." not in text:
        return 0
    return len(text.split(".", 1)[1].rstrip())


def decimal_value(text: str, field: str, index: int) -> Decimal:
    try:
        return Decimal(text)
    except InvalidOperation as exc:
        raise ValueError(f"Row {index}: {field} must be numeric, got {text!r}.") from exc


def format_decimal(value: Decimal, places: int) -> str:
    if places == 0:
        return str(value.quantize(Decimal("1")))
    quantum = Decimal("1").scaleb(-places)
    return format(value.quantize(quantum), "f")


def energy_mev(text: str, field: str, index: int) -> tuple[Decimal, int]:
    value = decimal_value(text, field, index)
    places = decimal_places(text)
    if abs(value) >= Decimal("1000"):
        return value / Decimal("1000"), max(3, places + 3)
    return value, places


def format_energy_mev(text: str, field: str, index: int) -> str:
    value, places = energy_mev(text, field, index)
    return format_decimal(value, places)


def fission_energies(row: dict[str, str], index: int) -> tuple[str, str, str]:
    dssd, dssd_places = energy_mev(row["E_1_DSSD"], "E_1_DSSD", index)
    ssd, ssd_places = energy_mev(row["E_1_SSD"], "E_1_SSD", index)
    places = max(dssd_places, ssd_places)
    total = dssd + ssd
    return (
        format_decimal(total, places),
        format_decimal(dssd, dssd_places),
        format_decimal(ssd, ssd_places),
    )


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def nuclide_label(text: str) -> str:
    text = text.strip()
    if any(token in text for token in ("\\", "^", "$")):
        math_text = text.strip("$")
        return rf"$\bm{{{math_text}}}$"

    match = re.fullmatch(r"([A-Za-z]+)\s*(\d+)", text)
    if match:
        element, mass = match.groups()
        return rf"$\bm{{{{}}^{{{mass}}}\mathrm{{{tex_escape(element)}}}}}$"

    match = re.fullmatch(r"(\d+)\s*([A-Za-z]+)", text)
    if match:
        mass, element = match.groups()
        return rf"$\bm{{{{}}^{{{mass}}}\mathrm{{{tex_escape(element)}}}}}$"

    return tex_escape(text)


def event_block(row: dict[str, str], index: int, x_shift: float, y_shift: float, nuclide: str) -> str:
    prefix = f"evt{index}"
    total, dssd, ssd = (tex_escape(value) for value in fission_energies(row, index))
    e_er = tex_escape(format_energy_mev(row["E_ER"], "E_ER", index))
    x_pos = tex_escape(row["X_Pos"])
    y_pos = tex_escape(row["Y_Pos"])
    delta_t = tex_escape(row["Delta_t1"])
    label = nuclide_label(nuclide)

    return rf"""
% Event {index}
\begin{{scope}}[shift={{({x_shift:.2f},{y_shift:.2f})}}]
  \node[nucW] ({prefix}implant) at (5.80,4.00) {{ER}};
  \node[nucSF] ({prefix}sf) at (5.80,0.95) {{{label}}};

  \draw[arr] ({prefix}implant.south) -- ({prefix}sf.north);
  \node[datatxt,anchor=west] at ($({prefix}implant.south)!0.50!({prefix}sf.north)+(0.25,0)$)
    {{$\Delta t={delta_t}\,\mu\mathrm{{s}}$}};
  \draw[arr] ({prefix}sf.south) -- ++(-0.72,-0.72);
  \draw[arr] ({prefix}sf.south) -- ++(0.72,-0.72);
  \node[sflabel,anchor=north] at ($({prefix}sf.south)+(0,-0.95)$) {{SF}};

  \node[datatxt,anchor=west] at ($({prefix}implant.east)+(0.35,0.18)$)
    {{{e_er} MeV\\{x_pos} / {y_pos} mm}};
  \node[datatxt,anchor=west] at ($({prefix}sf.east)+(0.45,0.22)$)
    {{{total} MeV\\\textit{{({dssd}+{ssd})}}}};
\end{{scope}}
""".strip()


def build_document(rows: list[dict[str, str]], nuclide: str) -> str:
    blocks = []
    col_sep = 9.0
    row_sep = 8.2
    for zero_based, row in enumerate(rows):
        col = zero_based % EVENTS_PER_ROW
        page_row = zero_based // EVENTS_PER_ROW
        blocks.append(event_block(row, zero_based + 1, col * col_sep, -page_row * row_sep, nuclide))

    body = "\n\n".join(blocks)
    return rf"""\documentclass[tikz,border=6pt]{{standalone}}
\usepackage{{fontspec}}
\usepackage{{newtxmath}}
\usepackage{{bm}}
\usepackage{{tikz}}
\usetikzlibrary{{calc,arrows.meta}}

\definecolor{{sfgreen}}{{RGB}}{{0,166,86}}
\setmainfont{{Times New Roman}}

\begin{{document}}
\begin{{tikzpicture}}[x=1cm,y=1cm]
\tikzset{{
  nucW/.style={{
    draw=black,
    line width=1.25pt,
    fill=white,
    minimum width=2.05cm,
    minimum height=2.05cm,
    align=center,
    inner sep=2pt,
    font=\bfseries\fontsize{{24}}{{26}}\selectfont
  }},
  nucSF/.style={{
    draw=black,
    line width=1.25pt,
    fill=sfgreen,
    minimum width=2.05cm,
    minimum height=2.05cm,
    align=center,
    inner sep=2pt,
    font=\bfseries\fontsize{{24}}{{26}}\selectfont
  }},
  arr/.style={{
    -{{Latex[length=3.5mm,width=2.5mm]}},
    line width=1.25pt,
    draw=black
  }},
  datatxt/.style={{
    font=\bfseries\fontsize{{18}}{{20}}\selectfont,
    text=black,
    align=left
  }},
  sflabel/.style={{
    font=\bfseries\fontsize{{18}}{{20}}\selectfont,
    text=black,
    align=center
  }}
}}

{body}

\end{{tikzpicture}}
\end{{document}}
"""


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=f"Generate a {EVENTS_PER_ROW}-per-row TikZ figure for 242Fm implantation-fission events."
    )
    parser.add_argument("input", nargs="?", type=Path, default=here / "input.txt")
    parser.add_argument("output", nargs="?", type=Path, default=here / "fissionChain.tex")
    parser.add_argument("--nuclide", default="Fm242", help="Label written inside the nuclide boxes.")
    parser.add_argument("--pdf", action="store_true", help="Run xelatex after writing the .tex file.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()

    try:
        rows = read_input(input_path)
        if not rows:
            raise ValueError("No event rows found after the header.")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(build_document(rows, args.nuclide), encoding="utf-8")
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {output_path} with {len(rows)} event(s), arranged {EVENTS_PER_ROW} per row.")

    if args.pdf:
        try:
            subprocess.run(
                ["xelatex", "-interaction=nonstopmode", output_path.name],
                cwd=output_path.parent,
                check=True,
            )
        except FileNotFoundError:
            print("Error: xelatex was not found. The .tex file was generated successfully.", file=sys.stderr)
            return 1
        except subprocess.CalledProcessError as exc:
            print(f"Error: xelatex failed with exit code {exc.returncode}.", file=sys.stderr)
            return exc.returncode
        print(f"Wrote {output_path.with_suffix('.pdf')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
