# ROOT Figure Redrawing Instructions

## 1. Goal

Use the supplied ROOT file

```text
plot_data.root
```

to redraw two figures with publication-quality ROOT/C++ plotting code.

The original experimental ROOT files are **not required**.

Do not reconstruct the event selection from the original data files. All data needed for plotting have already been extracted into `plot_data.root`.

The task is mainly to experiment with different plotting styles, layouts, fonts, line styles, marker styles, axis formatting, legends, and publication-quality appearance.

The underlying data must not be modified.

---

# 2. Contents of `plot_data.root`

The ROOT file contains the following main objects:

```text
h_SumE_run4_11
h_SumE_run12_40

g_DeltaT_run4_11
g_DeltaT_run12_40

tree_run4_11
tree_run12_40
```

Additional metadata are also stored:

```text
selection_info
Emin_keV
Emax_keV
spectrum_bin_width_keV
```

---

# 3. Meaning of the datasets

There are two experimental run groups:

```text
Run 4-11
Run 12-40
```

They should always be shown separately so that their distributions can be compared.

The colors/styles may be changed freely when testing plotting styles.

However, the same run group must use a consistent style within one figure.

For example:

```text
Run 4-11  -> green / blue / dark green / dashed
Run 12-40 -> red / black / dark red / solid
```

Different versions may use different combinations.

---

# 4. Figure 1: Sum-energy spectrum

## Physics content

The first figure shows the energy spectrum

```text
SumE[1]
```

for events satisfying the additional time condition

```text
DeltaT < 10 s
```

The corresponding histograms are already stored as

```cpp
TH1D *h1 =
    (TH1D*)f->Get("h_SumE_run4_11");

TH1D *h2 =
    (TH1D*)f->Get("h_SumE_run12_40");
```

These histograms use

```text
Energy range: 0 - 250000 keV
Bin width:    20 keV
```

Do not alter the histogram contents.

---

## Required presentation

Use a logarithmic Y axis:

```cpp
gPad->SetLogy();
```

Suggested axis labels:

```text
X: Energy (keV)
Y: Counts
```

or

```text
X: Sum energy (keV)
Y: Counts / 20 keV
```

The second form is preferred if the bin width is kept at 20 keV.

The minimum of the logarithmic axis should be positive, for example:

```cpp
hist->SetMinimum(0.5);
```

Do not use a linear Y axis for the final version unless explicitly testing an alternative style.

---

# 5. Figure 2: correlation between energy and decay time

The second figure is a scatter plot showing

```text
X = SumE
Y = log10(DeltaT [ns])
```

The corresponding `TGraph` objects are

```cpp
TGraph *g1 =
    (TGraph*)f->Get("g_DeltaT_run4_11");

TGraph *g2 =
    (TGraph*)f->Get("g_DeltaT_run12_40");
```

Approximate display ranges:

```text
X: 0 - 250000 keV
Y: 2 - 11
```

This must be a **scatter plot**.

Do NOT use

```text
COLZ
```

and do not convert the data into a density heat map unless specifically generating an additional comparison version.

Use markers such as

```cpp
SetMarkerStyle(20);
SetMarkerSize(...);
```

with different colors or marker styles for the two run groups.

---

# 6. Alternative approach using the stored TTrees

The ROOT file also contains

```text
tree_run4_11
tree_run12_40
```

The branches are approximately:

```text
SumE
DeltaT_ns
DeltaT_s
log10_DeltaT_ns
pass_10s
```

These trees can be used if a new histogram binning or another time-axis representation is needed.

For example:

```cpp
tree->Draw(
    "SumE>>h(12500,0,250000)",
    "pass_10s",
    "goff"
);
```

However, for simple restyling, use the already stored `TH1D` and `TGraph` objects.

The TTrees should only be used when the plot itself needs to be reconstructed.

---

# 7. Important rule: separate data from style

The supplied ROOT file should be treated as a **data container**.

Do not overwrite or modify it.

Create a separate ROOT macro, for example:

```text
draw_plot_data.C
```

or several style variants:

```text
draw_style_A.C
draw_style_B.C
draw_style_C.C
```

The plotting script should contain all stylistic definitions:

```text
canvas size
margins
fonts
font sizes
axis divisions
tick directions
line width
line style
marker size
marker style
colors
legend position
labels
annotations
output format
```

The experimental data must remain untouched.

---

# 8. ROOT file loading

Start the plotting code approximately as follows:

```cpp
TFile *f = TFile::Open("plot_data.root");

if (!f || f->IsZombie()) {
    std::cerr << "Cannot open plot_data.root" << std::endl;
    return;
}
```

Then retrieve the objects:

```cpp
TH1D *h_run4_11 =
    (TH1D*)f->Get("h_SumE_run4_11");

TH1D *h_run12_40 =
    (TH1D*)f->Get("h_SumE_run12_40");

TGraph *g_run4_11 =
    (TGraph*)f->Get("g_DeltaT_run4_11");

TGraph *g_run12_40 =
    (TGraph*)f->Get("g_DeltaT_run12_40");
```

Check that every object exists before drawing.

---

# 9. General figure style

The figures are intended for a nuclear-physics paper or presentation.

The style should resemble modern figures in journals such as:

```text
Physical Review C
Physical Review Letters
Physics Letters B
```

Avoid the default ROOT appearance.

Use a clean scientific style with:

```text
white background
no statistics box
no histogram title
thick enough axes
clear typography
proper margins
compact legends
publication-quality line widths
well-sized tick marks
```

Recommended:

```cpp
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
```

---

# 10. Fonts

Prefer a serif scientific font.

A preferred ROOT font code is:

```cpp
132
```

Use it consistently for:

```text
axis titles
axis labels
legend
annotations
panel labels
```

Example:

```cpp
axis->SetTitleFont(132);
axis->SetLabelFont(132);
```

Avoid mixing multiple font families.

---

# 11. Axis formatting

The axes should be more polished than ROOT defaults.

Suggested settings:

```cpp
SetTitleSize(0.055);
SetLabelSize(0.050);
SetTitleOffset(...);
SetNdivisions(505);
```

Adjust offsets according to the canvas proportions.

Make sure:

```text
axis titles do not overlap labels
labels are not clipped
tick marks are clearly visible
```

Use outward or visually balanced ticks where appropriate.

---

# 12. Line and marker appearance

For Figure 1, use histogram lines rather than filled histograms.

Example:

```cpp
h1->SetLineWidth(2);
h2->SetLineWidth(2);
```

Potential style combinations to test:

### Style A

```text
Run 4-11  : blue dashed
Run 12-40 : black solid
```

### Style B

```text
Run 4-11  : dark green solid
Run 12-40 : red solid
```

### Style C

```text
Run 4-11  : gray dashed
Run 12-40 : black solid
```

For Figure 2, use relatively small markers because many events may overlap.

For example:

```cpp
SetMarkerStyle(20);
SetMarkerSize(0.4);
```

or smaller if necessary.

---

# 13. Legend

The legend should be compact and unobtrusive.

Labels should be:

```text
Run 4–11
Run 12–40
```

Use an en dash if practical.

Recommended ROOT settings:

```cpp
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextFont(132);
```

Do not use a large boxed default ROOT legend.

---

# 14. Figure layout

Initially generate the two figures independently:

```text
Figure 1:
Sum-energy spectrum

Figure 2:
Energy-time scatter plot
```

After the individual versions look good, also test a combined publication-style figure:

```text
(a) Sum-energy spectrum
(b) Energy-time correlation
```

with two vertically stacked panels.

For a combined figure, use panel labels

```text
(a)
(b)
```

near the upper-left corner of each panel.

The panels should share a consistent visual style.

---

# 15. Output files

For every preferred style, save both

```text
PDF
PNG
```

For example:

```cpp
canvas->SaveAs("SumE_styleA.pdf");
canvas->SaveAs("SumE_styleA.png");
```

and

```cpp
canvas->SaveAs("DeltaT_styleA.pdf");
canvas->SaveAs("DeltaT_styleA.png");
```

For combined versions:

```text
combined_styleA.pdf
combined_styleA.png
```

PDF is the primary publication-quality output.

---

# 16. What should NOT be done

Do not:

```text
modify event values
modify histogram contents
re-run the original event selection
depend on the original runXXXXX_Rechain.root files
use COLZ for the scatter plot
use ROOT default titles/statistics boxes
use excessively bright colors
use thick filled histograms
use very large markers
hide one dataset behind the other
```

Do not treat the current colors as part of the physics.

Colors are only graphical identifiers.

---

# 17. Requested workflow

Please proceed in the following way.

First inspect the file:

```cpp
TFile f("plot_data.root");
f.ls();
```

Verify the stored objects.

Then write one clean ROOT macro that reproduces the two figures.

After that, produce several clearly different visual styles without altering the physics data.

A useful set would be:

```text
Style A:
PRC-like traditional serif style

Style B:
clean modern black/blue publication style

Style C:
presentation-friendly high-contrast style
```

Finally, choose one style as the recommended publication version, but keep all style implementations easy to modify.

---

# 18. Code quality

The ROOT macro should be easy to maintain.

Avoid repeating the same formatting code many times.

Create helper functions if useful, for example:

```cpp
void SetAxisStyle(TAxis *axis);

void StyleHistogram(
    TH1D *h,
    int color,
    int lineStyle
);

void StyleGraph(
    TGraph *g,
    int color,
    int markerStyle
);
```

The plotting script should clearly separate:

```text
1. file loading
2. data retrieval
3. style configuration
4. canvas creation
5. drawing
6. legends / annotations
7. output
```

The final code should compile or run directly in ROOT without manual editing.

---

# 19. Main principle

The most important principle is:

> `plot_data.root` contains the physics data; the new ROOT macros should only determine how those data are presented.

The goal is therefore to generate several visually different, publication-quality representations of exactly the same datasets without touching the underlying event information.
