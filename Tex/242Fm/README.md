# 242Fm implantation-fission TikZ generator

Edit `input.txt`, keeping the header:

```text
E_ER,X_Pos,Y_Pos,Delta_t1,E_1_DSSD,E_1_SSD
13.13,23,109,1.6350,168.3,20.0
12.84,41,96,0.8746,151.8,22.6
```

The header and data rows may be comma-, space-, or tab-separated. Energy values
larger than or equal to `1000` are treated as keV and displayed as MeV.

Then run:

```sh
python3 generate_fission_tikz.py
```

This writes `fissionChain.tex`. Every three input rows are placed on one TikZ row.

To regenerate all outputs after editing `input.txt`:

```sh
./run_all.sh
```

This runs `generate_fission_tikz.py`, compiles `fissionChain.tex` with `xelatex`,
and writes the `Delta_t1` distribution plot.
By default the plotting step uses `plot_delta_t.py` directly so it keeps the
script's configured Python environment. To force another Python executable, run
for example `PYTHON=python3 ./run_all.sh`.

To generate the PDF directly:

```sh
python3 generate_fission_tikz.py --pdf
```

You can also change the box label:

```sh
python3 generate_fission_tikz.py --nuclide "Fm242"
```

Labels like `Fm242`, `Fm 242`, and `242Fm` are rendered as
`${}^{242}\mathrm{Fm}$` in the TikZ output.

To plot the `Delta_t1` values from `input.txt` as microseconds and print the
half-life estimate:

```sh
./plot_delta_t.py
```

This writes `delta_t_distribution.pdf` and prints the asymmetric confidence
interval for `T_1/2` in `us`.
