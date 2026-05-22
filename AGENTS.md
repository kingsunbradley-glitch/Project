# Repository Guidelines

## Workspace Shape

This repository is a collection of nuclear-physics analysis utilities, ROOT/C++
calibration programs, Python plotting/calculation scripts, and LaTeX manuscript
material. There is no single root-level build or test command. Work inside the
specific subdirectory that owns the task.

Important areas:

- `A80_fixedN_final_package/`: ROOT/C++ A80 run processing package.
- `DSSD_recal_all/`: DSSD pipeline: preselection, normalization, and
  recalibration.
- `DSSD_outcal_all_byCh/`: independent per-channel 3-peak calibration.
- `SSD_calibration/`: SSD internal alpha-source calibration workflow.
- `CrossSection/`: Python + ROOT beam-dose and cross-section calculations.
- `Tex/`: LaTeX manuscripts and TikZ generation utilities.
- `tool/`: assorted analysis, plotting, half-life, mass-table, and fitting
  scripts.
- Directories containing `backup`, `_backup`, or older numbered variants should
  be treated as historical unless the user explicitly targets them.

## Environment

- ROOT is required for most C++ and PyROOT workflows. Use `root-config` from the
  active shell rather than hard-coding include or library paths.
- Current local tools observed during initialization:
  - ROOT `6.34.10`
  - Python `3.13.5`
  - Apple clang `17.0.0`
- C++ standards vary by subproject (`C++11`, `C++17`, and `C++20` all appear).
  Follow the local Makefile or existing compile flags.
- Python scripts are standalone; inspect imports before assuming dependencies.
  Common likely needs include PyROOT, matplotlib, pandas, numpy, scipy, and
  openpyxl depending on the script.

## Build And Run Commands

Run commands from the owning subdirectory.

- `A80_fixedN_final_package/`
  - `make`
  - Builds `process_runs_A80`.
- `DSSD_recal_all/`
  - `make`
  - `make Pre`
  - `make Norm`
  - `make Cali`
  - `make run`
  - Configuration is centralized in `Config.h`.
- `DSSD_outcal_all_byCh/`
  - `make`
  - `make run IN=/path/to/preselected.root TREE=tr_map`
- `SSD_calibration/`
  - `make`
  - `make run`
  - Configuration is centralized in `config.h`.
- `Tex/242Fm/`
  - `python3 generate_fission_tikz.py`
  - `python3 generate_fission_tikz.py --pdf`
  - `./plot_delta_t.py`

There is no global automated test suite. Verification usually means compiling
the targeted program, running the relevant script on the expected input, and
checking generated ROOT files, plots, tables, or PDFs.

## Editing Conventions

- Keep changes tightly scoped to the requested directory and workflow.
- Preserve the existing language and style of the local files. Many docs and
  comments are in Chinese; do not translate them unless asked.
- Prefer the local configuration header (`Config.h` or `config.h`) over
  hard-coded paths or constants in C++ source files.
- Do not mass-format or reorganize ROOT macros, legacy analysis scripts, or
  generated LaTeX unless the user requested cleanup.
- Treat compiled binaries, generated ROOT files, plots, PDFs, and intermediate
  data products as outputs. The root `.gitignore` ignores many of these, but
  some historical outputs are already tracked.
- Never overwrite user data files or calibration result tables without checking
  the targeted workflow and expected outputs.

## Git Hygiene

- The worktree may contain user edits. Do not revert unrelated changes.
- At initialization the worktree had:
  - Modified: `tempForLearn/main.cpp`
  - Untracked: `Tex/242Fm/`
- Avoid broad `git add .` operations in this repository because generated and
  historical files are mixed with source files.
