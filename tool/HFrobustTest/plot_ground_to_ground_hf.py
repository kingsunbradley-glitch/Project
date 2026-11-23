#!/usr/bin/env python3
"""Export and plot the even-even ground-to-ground HF for VS and Rurarz."""

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, MultipleLocator, NullFormatter

OUTPUTS = Path(__file__).resolve().parent / "outputs"
STEM = "ground_to_ground_hf"
FORMULAS = {
    "Viola_Seaborg": ("vs", "Viola–Seaborg (VS)"),
    "E_Rurarz": ("rurarz", "E. Rurarz"),
}
UNCERTAIN_FLAGS = {
    "expected_alpha", "uncertain_alpha", "unresolved_alpha_multiplet",
    "energy_balance_inconsistent",
}
LIMIT_DIRECTION = {"<": -1, "≤": -1, ">": 1, "≥": 1}
REGION_COLORS = {"actinide": "#275AB4", "superheavy": "#D32F2F"}


def is_ground(transition, state):
    return transition[f"E{state}_keV"] == 0 and "+" not in str(transition[f"E{state}_raw"])


def export_rows(database):
    results = defaultdict(list)
    for result in database["results"]:
        if result["method"] in FORMULAS:
            results[result["transition_id"]].append(result)

    # Match the atlas definition using all records before applying plot selections.
    ground_branches = Counter()
    maximum_energy = {}
    for transition in database["transitions"]:
        if transition["parity_class"] != "ee" or not is_ground(transition, "i"):
            continue
        parent, energy = transition["parent"], transition["Ealpha_MeV"]
        ground_branches[parent] += 1
        if energy is not None:
            maximum_energy[parent] = max(energy, maximum_energy.get(parent, -math.inf))

    rows = []
    for transition in database["transitions"]:
        if transition["parity_class"] != "ee" or not is_ground(transition, "i"):
            continue
        parent, energy = transition["parent"], transition["Ealpha_MeV"]
        highest_even_even = energy is not None and energy == maximum_energy.get(parent)
        if not (is_ground(transition, "f") or highest_even_even):
            continue

        uncertain_flags = UNCERTAIN_FLAGS.intersection(transition["flags"])
        branch = transition["b_alpha"]
        half_life = transition["T_total_s"]
        single_line = (
            ground_branches[parent] == 1
            and transition["preferred"] and not uncertain_flags and energy is not None
            and transition["I_alpha"] is None and not (transition["I_input_record"] or {}).get("raw")
            and branch is not None and 0 < branch <= 1 and half_life is not None and half_life > 0
        )
        common = {key: transition[key] for key in (
            "A", "Z", "N", "parent", "daughter", "parity_class", "region", "Ei_keV", "Ef_keV",
            "Ealpha_MeV", "Qalpha_MeV", "T_total_s", "b_alpha", "I_alpha", "preferred",
        )}
        for result in results[transition["id"]]:
            hf, qualifier = result["HF"], result["HF_qualifier"]
            source = "evaluated_inputs" if hf is not None else "unavailable"
            if hf is None and single_line and result["Tcalc_s"] and result["Tcalc_s"] > 0:
                hf = half_life / branch / result["Tcalc_s"]
                qualifier = qualifier or "≈"
                source = "single_line_reference"

            if hf is None or not math.isfinite(hf) or hf <= 0:
                reason = "no_positive_HF"
            elif not transition["preferred"]:
                reason = "nonpreferred_record"
            else:
                reason = ""
            rows.append({
                **common,
                "transition_id": transition["id"],
                "ground_to_ground_basis": "zero_level_energies" if is_ground(transition, "f") else "highest_even_even_Ealpha",
                "b_alpha_qualifier": (transition["b_alpha_record"] or {}).get("qualifier", ""),
                "reference_I_alpha": 1 if source == "single_line_reference" else None,
                "method": result["method"], "version": result["version"],
                "mode": result["mode"], "l": result["l"], "candidate_index": result["candidate_index"],
                "formal_l0": result["formal_l0"], "Tcalc_s": result["Tcalc_s"],
                "HF": hf, "HF_qualifier": qualifier, "HF_source": source,
                "log10_HF_err_minus": result["log10_HF_err_minus"],
                "log10_HF_err_plus": result["log10_HF_err_plus"],
                "plotted": not reason, "plot_exclusion": reason,
                "flags": "|".join(transition["flags"]), "source_url": transition["source_url"],
            })

    method_order = {method: index for index, method in enumerate(FORMULAS)}
    return sorted(rows, key=lambda row: (
        row["A"], row["Z"], row["transition_id"], method_order[row["method"]],
        row["mode"], row["candidate_index"] if row["candidate_index"] is not None else -1,
    ))


def origin_style():
    font_dir = Path("/mnt/c/Windows/Fonts")
    if (font_dir / "times.ttf").exists():
        for name in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
            path = font_dir / name
            if path.exists():
                font_manager.fontManager.addfont(path)
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "font.size": 10, "font.weight": "bold", "axes.linewidth": .8, "mathtext.fontset": "stix",
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def plot_formula(ax, rows):
    for row in rows:
        mass, hf, qualifier = row["A"], row["HF"], row["HF_qualifier"]
        color = REGION_COLORS[row["region"]]
        reference = row["HF_source"] == "single_line_reference"
        mixed = qualifier == "mixed_constraints"
        face = "none" if reference or mixed else color

        ax.plot(mass, hf, marker="o", color=color, markerfacecolor=face,
                ms=3.2, mew=.8, linestyle="none", zorder=3)
        if mixed:
            ax.plot(mass, hf, marker="x", color=color, ms=2.5, mew=.7,
                    linestyle="none", zorder=4)

        direction = LIMIT_DIRECTION.get(qualifier)
        if direction:
            ax.annotate("", xy=(mass, hf * 1.45**direction), xytext=(mass, hf),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=.75,
                                        mutation_scale=7, shrinkA=2.5, shrinkB=0), zorder=4)


def draw_figure(rows, title, ylim, output):
    fig, ax = plt.subplots(figsize=(190 / 25.4, 112 / 25.4))
    fig.subplots_adjust(left=.11, right=.97, bottom=.25, top=.87)
    ax.set(yscale="log", xlim=(200, 300), ylim=ylim,
           xlabel=r"Parent mass number $A$", ylabel="Hindrance factor HF")
    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=9)
    ax.axhline(1, color=".55", linestyle=(0, (3, 3)), linewidth=.65, zorder=1)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=12))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(3, 4, 6, 7, 8, 9)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", length=5, width=.8, pad=5)
    ax.tick_params(which="minor", length=2.5, width=.6)
    plot_formula(ax, rows)

    handles = [
        Line2D([], [], marker="o", color=REGION_COLORS["actinide"], ms=4,
               linestyle="none", label="Actinides"),
        Line2D([], [], marker="o", color=REGION_COLORS["superheavy"], ms=4,
               linestyle="none", label="Superheavy nuclei"),
    ]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1, 1.015), ncol=2,
              frameon=False, borderaxespad=0)
    fig.text(
        .11, .05,
        r"Open circle: single-line reference ($I_\alpha=1$)" + "   ·   "
        "Cross inside circle: mixed constraints   ·   Arrow: HF limit",
        fontsize=9,
    )
    for extension in ("pdf", "svg", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), dpi=600)
    plt.close(fig)


def main():
    database = json.loads((OUTPUTS / "database.json").read_text(encoding="utf-8"))
    rows = export_rows(database)
    csv_path = OUTPUTS / f"{STEM}_vs_rurarz.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    origin_style()
    plotted = [row for row in rows if row["plotted"]]
    ylim = (min(row["HF"] for row in plotted) / 10**.24,
            max(row["HF"] for row in plotted) * 10**.24)
    for method, (slug, title) in FORMULAS.items():
        draw_figure([row for row in plotted if row["method"] == method], title, ylim,
                    OUTPUTS / "figures" / f"{STEM}_{slug}")
    print(f"{len({row['transition_id'] for row in rows})} transitions, {len(rows)} result rows, "
          f"{sum(row['plotted'] for row in rows)} plotted results")


if __name__ == "__main__":
    main()
