#!/usr/bin/env python3
"""按 nature-figure 流程绘制 Ds 激发函数及给定实验约束。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
CHANNELS = {
    "275Ds": {"channel": "3n", "color": "#666666", "style": "-."},
    "274Ds": {"channel": "4n", "color": "#0072B2", "style": "--"},
    "273Ds": {"channel": "5n", "color": "#D55E00", "style": "-"},
    "272Ds": {"channel": "6n", "color": "#009E73", "style": ":"},
}


def read_calculation(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        expected = {"Elab", "E*", *CHANNELS}
        if set(reader.fieldnames or []) != expected:
            raise ValueError(f"Unexpected calculation columns: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError("The calculation table is empty")
    data = {key: np.array([float(row[key]) for row in rows]) for key in expected}
    if not all(np.isfinite(values).all() for values in data.values()):
        raise ValueError("Calculation table contains non-finite values")
    if not np.all(np.diff(data["E*"]) > 0):
        raise ValueError("Excitation energies must be strictly increasing")
    for isotope in CHANNELS:
        if np.any(data[isotope] < 0):
            raise ValueError(f"Negative cross section in {isotope}")
    return data


def read_experiment(path):
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    seen = set()
    for row in rows:
        for key in ("E_lab_MeV", "E_star_MeV", "sigma_pb"):
            row[key] = float(row[key])
        for key in ("err_minus_pb", "err_plus_pb"):
            row[key] = float(row[key]) if row[key] else None
        isotope = row["residue"]
        if isotope not in CHANNELS or CHANNELS[isotope]["channel"] != row["channel"]:
            raise ValueError("Inconsistent residue/channel assignment")
        if row["kind"] not in {"measurement", "upper_limit"}:
            raise ValueError("Unknown experimental record kind")
        values = [v for v in row.values() if isinstance(v, float)]
        if not np.isfinite(values).all() or row["sigma_pb"] <= 0:
            raise ValueError("Invalid experimental value for a logarithmic plot")
        if row["kind"] == "measurement":
            errors = (row["err_minus_pb"], row["err_plus_pb"])
            if any(error is None or error < 0 for error in errors):
                raise ValueError("Measured cross sections require nonnegative errors")
            if row["sigma_pb"] - errors[0] <= 0:
                raise ValueError("Measured interval crosses zero on a log axis")
        elif row["err_minus_pb"] is not None or row["err_plus_pb"] is not None:
            raise ValueError("An upper limit must not carry a fabricated error")
        key = (row["channel"], row["E_star_MeV"])
        if key in seen:
            raise ValueError(f"Duplicate experimental record: {key}")
        seen.add(key)
    if not rows:
        raise ValueError("The experimental table is empty")
    return rows


def find_skill(explicit):
    candidates = []
    if explicit:
        candidates.append(explicit)
    else:
        if os.environ.get("CODEX_HOME"):
            candidates.append(Path(os.environ["CODEX_HOME"]) / "skills/nature-figure")
        candidates.append(Path.home() / ".codex/skills/nature-figure")
    for candidate in candidates:
        if (candidate / "scripts/audit_panel_alignment.py").is_file():
            return candidate.resolve()
    raise FileNotFoundError("Install nature-figure or specify --skill-dir")


def clip_log_polyline(x, y, xlim, ylim):
    """在 E*-log(sigma) 平面裁切原折线，防止 PDF 保留图框外的不可见线段。

    仅添加原线段与图框的交点，不插值分析数据，不改变图内的线形。
    零值两侧的线段不连接；完整输入始终保留在总览图中。
    """
    visible_x, visible_y = [], []
    bounds = [(xlim[0], xlim[1]), (np.log10(ylim[0]), np.log10(ylim[1]))]
    for i in range(len(x) - 1):
        if y[i] <= 0 or y[i + 1] <= 0:
            continue
        start = np.array([x[i], np.log10(y[i])])
        end = np.array([x[i + 1], np.log10(y[i + 1])])
        delta = end - start
        enter, leave = 0.0, 1.0
        for dim, (low, high) in enumerate(bounds):
            if delta[dim] == 0:
                if not low <= start[dim] <= high:
                    enter, leave = 1.0, 0.0
                    break
            else:
                intersections = sorted(((low - start[dim]) / delta[dim],
                                        (high - start[dim]) / delta[dim]))
                enter = max(enter, intersections[0])
                leave = min(leave, intersections[1])
        if enter > leave:
            continue
        p, q = start + enter * delta, start + leave * delta
        if not visible_x or not (np.isclose(visible_x[-1], p[0])
                                and np.isclose(visible_y[-1], 10 ** p[1])):
            visible_x.extend([np.nan, p[0]])
            visible_y.extend([np.nan, 10 ** p[1]])
        visible_x.append(q[0])
        visible_y.append(10 ** q[1])
    return visible_x, visible_y


def plot_calculation(ax, data, isotope, *, clip_to_view=False):
    spec = CHANNELS[isotope]
    # 零值不能在对数坐标显示。掩码保留断点，不跨零值连接、不平滑。
    values = np.ma.masked_less_equal(data[isotope], 0)
    x = data["E*"]
    if clip_to_view:
        x, values = clip_log_polyline(x, data[isotope], ax.get_xlim(), ax.get_ylim())
    return ax.plot(
        x, values, color=spec["color"], linestyle=spec["style"],
        linewidth=1.3, label=rf"{spec['channel']}: $^{{{isotope[:3]}}}$Ds",
    )[0]


def plot_experiment(ax, rows, isotope):
    color = CHANNELS[isotope]["color"]
    for row in rows:
        if row["residue"] != isotope:
            continue
        x, y = row["E_star_MeV"], row["sigma_pb"]
        if row["kind"] == "upper_limit":
            # 横帽位于真实上限；向下箭头仅是上限符号，长度不是误差。
            ax.plot(x, y, marker="_", markersize=7, markeredgewidth=1.05,
                    color=color, linestyle="none", zorder=4)
            ax.annotate("", xy=(x, y / 2), xytext=(x, y),
                        arrowprops={"arrowstyle": "-|>", "color": color,
                                    "lw": 1.05, "mutation_scale": 8,
                                    "shrinkA": 0, "shrinkB": 0}, zorder=4)
        else:
            ax.errorbar(x, y, yerr=[[row["err_minus_pb"]], [row["err_plus_pb"]]],
                        fmt="o", color=color, markersize=4.5, capsize=3,
                        elinewidth=1.05, capthick=1.05, zorder=5)


def data_report(data, experiment, calc_path):
    report = {
        "calculation_sha256": hashlib.sha256(calc_path.read_bytes()).hexdigest(),
        "calculation_rows": len(data["E*"]),
        "experimental_records": len(experiment),
        "experimental_measurements": sum(r["kind"] == "measurement" for r in experiment),
        "experimental_upper_limits": sum(r["kind"] == "upper_limit" for r in experiment),
        "energy_basis": "Use each input's E* directly; no shared Elab conversion",
        "zero_policy": "Mask zero cross sections on logarithmic axes; preserve gaps",
        "error_confidence_level": "Not specified in the supplied table",
        "channel_statistics": {},
        "comparison_at_experimental_energies": [],
    }
    for isotope in CHANNELS:
        values = data[isotope]
        peak = int(np.argmax(values))
        report["channel_statistics"][isotope] = {
            "positive_count": int(np.count_nonzero(values > 0)),
            "zero_count": int(np.count_nonzero(values == 0)),
            "sampled_peak_E_star_MeV": float(data["E*"][peak]),
            "sampled_peak_sigma_pb": float(values[peak]),
            "positive_points_outside_overview": int(np.count_nonzero(
                (values > 0) & ((values < 1e-8) | (values > 1)))),
            "positive_points_outside_detail_axes": int(np.count_nonzero(
                (values > 0) & ((values < 1e-3) | (values > 2.5)
                                | (data["E*"] < 30) | (data["E*"] > 65)))),
        }
    for row in experiment:
        matching = np.flatnonzero(np.isclose(data["E*"], row["E_star_MeV"], atol=1e-9, rtol=0))
        if matching.size != 1:
            raise ValueError("Experimental energy is absent from calculation grid; no implicit interpolation")
        i = matching[0]
        report["comparison_at_experimental_energies"].append({
            "channel": row["channel"], "E_star_MeV": row["E_star_MeV"],
            "calculated_sigma_pb": float(data[row["residue"]][i]),
            "experimental_sigma_pb": row["sigma_pb"], "kind": row["kind"],
            "calculation_Elab_MeV": float(data["Elab"][i]),
            "experimental_Elab_MeV": row["E_lab_MeV"],
        })
    return report


def run_audit(script, args, destination):
    result = subprocess.run([sys.executable, str(script), *map(str, args)],
                            capture_output=True, text=True, check=False)
    destination.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Figure audit failed: {script.name}; see {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calculation", type=Path, default=HERE / "ExFunCalc.csv")
    parser.add_argument("--experiment", type=Path, default=HERE / "experimental_data.csv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    parser.add_argument("--skill-dir", type=Path)
    args = parser.parse_args()
    data = read_calculation(args.calculation)
    experiment = read_experiment(args.experiment)
    # skill = find_skill(args.skill_dir)
    # sys.path.insert(0, str(skill / "scripts"))
    # from audit_panel_alignment import require_matplotlib_panel_alignment

    output = args.output_dir.resolve()
    qa = output / "qa"
    qa.mkdir(parents=True, exist_ok=True)
    # run_audit(skill / "scripts/validate_figure.py", [Path(__file__), "--json"], qa / "source.json")
    report = data_report(data, experiment, args.calculation)
    (qa / "data.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "mathtext.fontset": "dejavusans", "font.size": 8,
        "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "xtick.direction": "out", "ytick.direction": "out",
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 92 / 25.4))
    fig.subplots_adjust(left=0.083, right=0.985, bottom=0.23, top=0.83, wspace=0.38)
    titles = ["Calculated channels", r"$4n$: $^{274}$Ds", r"$5n$: $^{273}$Ds"]
    for letter, ax, title in zip("abc", axes, titles):
        ax.set_yscale("log")
        ax.set_title(title, pad=10)
        ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-13, 10),
                    textcoords="offset points", fontweight="bold", fontsize=9,
                    ha="left", va="baseline")
        ax.tick_params(width=0.7, length=3)
        ax.tick_params(which="minor", width=0.5, length=1.6)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=[2, 5]))
        ax.yaxis.set_minor_formatter(NullFormatter())

    for isotope in CHANNELS:
        plot_calculation(axes[0], data, isotope)
    axes[0].set(xlim=(30, 80), ylim=(1e-8, 1))
    axes[0].set_xticks([30, 40, 50, 60, 70, 80])
    axes[0].set_yticks([1e-8, 1e-6, 1e-4, 1e-2, 1])
    axes[0].legend(loc="lower right", borderaxespad=0.45, handlelength=2.2,
                   labelspacing=0.4, handletextpad=0.6)
    for ax, isotope in zip(axes[1:], ("274Ds", "273Ds")):
        ax.set(xlim=(30, 65), ylim=(1e-3, 2.5))
        plot_calculation(ax, data, isotope, clip_to_view=True)
        plot_experiment(ax, experiment, isotope)
        ax.set_xticks([30, 40, 50, 60])
        ax.set_yticks([1e-3, 1e-2, 1e-1, 1])

    fig.supylabel(r"Cross section, $\sigma$ (pb)", x=0.012, y=0.53, fontsize=9)
    fig.supxlabel(r"Excitation energy, $E^{*}$ (MeV)", y=0.105, fontsize=9)
    fig.text(0.5, 0.025, "Lines: calculation    Arrows: upper limits    Point: measurement",
             ha="center", va="bottom", fontsize=8)
    fig.canvas.draw()
    # require_matplotlib_panel_alignment(
    #     fig, json_out=qa / "alignment.json", overlay_svg=qa / "alignment.svg",
    #     tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True,
    # )
    stem = output / "excitation_functions"
    # 不使用 tight 裁切，保证最终物理尺寸为 183 × 92 mm。
    fig.savefig(stem.with_suffix(".pdf"), format="pdf")
    fig.savefig(stem.with_suffix(".svg"), format="svg")
    fig.savefig(stem.with_suffix(".png"), format="png", dpi=600)
    plt.close(fig)
    pdf = stem.with_suffix(".pdf")
    # run_audit(skill / "scripts/audit_pdf_text.py", [pdf, "--min-pt", "5", "--json"], qa / "pdf-text.json")
    # run_audit(skill / "scripts/audit_figure_collisions.py",
            #   [pdf, "--json-out", qa / "collisions.json", "--overlay-pdf", qa / "collisions.pdf"],
            #   qa / "collisions.txt")
    print(f"Saved {stem}.pdf/.svg/.png")
    print(json.dumps(report["comparison_at_experimental_energies"], indent=2))


if __name__ == "__main__":
    main()
