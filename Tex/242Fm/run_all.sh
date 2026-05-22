#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON_BIN="${PYTHON:-python3}"

echo "[1/3] Generating fissionChain.tex"
"$PYTHON_BIN" generate_fission_tikz.py

echo "[2/3] Compiling fissionChain.pdf"
xelatex -interaction=nonstopmode -halt-on-error fissionChain.tex

echo "[3/3] Plotting delta_t_distribution.pdf"
if [[ -n "${PYTHON:-}" ]]; then
  "$PYTHON_BIN" plot_delta_t.py
else
  ./plot_delta_t.py
fi

echo "Done."
