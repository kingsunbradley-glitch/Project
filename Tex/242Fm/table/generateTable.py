input_file = "SF.dat"
output_file = "table_ER_SF.tex"

COLUMNS = [
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

def read_rows(path):
    rows = []
    with open(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            values = line.split()
            if len(values) != len(COLUMNS):
                raise ValueError(
                    f"{path}:{line_no}: expected {len(COLUMNS)} columns, got {len(values)}"
                )

            rows.append(dict(zip(COLUMNS, map(float, values))))
    return rows

def fmt_x(x):
    if abs(x - round(x)) < 1e-6:
        return f"{int(round(x))}"
    return f"{x:.0f}"

def fmt_y(y):
    if abs(y - round(y)) < 1e-6:
        return f"{int(round(y))}"
    return f"{y:.0f}"

rows = read_rows(input_file)

with open(output_file, "w") as f:
    f.write(r"""\documentclass{article}
\usepackage[margin=1in]{geometry}
\begin{document}

\begin{table*}[htbp]
\centering
\begin{tabular}{ccccccc}
\hline
Chain no. & File no. & X strip & Y strip & Type & Energy (MeV) & $\Delta T$ ($\mu$s) \\
\hline
""")

    for row in rows:
        chain_no = int(row["chain_no"])
        file_no = int(row["file_no"])
        xstrip = fmt_x(row["x_strip"])
        ystrip = fmt_y(row["y_strip"])

        er_e = row["ER_E"]
        sf_e = row["SF_E"]
        dssd_e = row["DSSD_E1"]
        ssd_e = row["SSD_E1"]
        dt_us = row["DeltaT_ms"] * 1000.0

        # ER row
        f.write(
            f"{chain_no} & {file_no} & {xstrip} & {ystrip} "
            f"& ER & {er_e:.2f} & 0 \\\\\n"
        )

        # SF row
        if ssd_e > 0.01:
            sf_energy = f"{sf_e:.1f}({dssd_e:.1f}+{ssd_e:.2f})"
        else:
            sf_energy = f"{sf_e:.1f}"

        f.write(
            f" & & & & "
            f"SF($^{{242}}$Fm) & {sf_energy} & {dt_us:.3f} \\\\\n"
        )

        f.write(r"\hline" + "\n")

    f.write(r"""\end{tabular}
\end{table*}
\end{document}
""")

print(f"Saved LaTeX table to {output_file}")
