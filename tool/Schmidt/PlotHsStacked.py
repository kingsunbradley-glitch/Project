"""绘制 269Hs 的 EVR 与衰变产生事件的堆叠能谱。"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


def main():
    data_path = Path(__file__).with_name("HsData.txt")
    decay_energies = []
    evr_energies = []
    with data_path.open(encoding="utf-8") as data_file:
        next(data_file)
        for line in data_file:
            values = line.split()
            if not values:
                continue
            # 行首留空表示缺少衰变能量，只有 EVR 数据。
            if line[0].isspace():
                evr_energies.append(float(values[-1]))
            else:
                decay_energies.append(float(values[0]))
                if len(values) > 1:
                    evr_energies.append(float(values[1]))

    # 从 8.810 MeV 起按 20 keV 分箱，并保留高能端留白。
    first_edge_kev = 8790
    bin_width_kev = 20
    energies = evr_energies + decay_energies
    energy_min, energy_max = min(energies), max(energies)
    padding = 0.2 * (energy_max - energy_min)
    last_offset_kev = (energy_max + padding) * 1000 - first_edge_kev
    bin_count = int(np.ceil(last_offset_kev / bin_width_kev))
    bins = (first_edge_kev + np.arange(bin_count + 1) * bin_width_kev) / 1000

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "serif"],
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "mathtext.fontset": "stix",
    })
    fig, ax = plt.subplots(figsize=(8, 3.8))
    stacked_counts, _, (_, decay_bars) = ax.hist(
        [evr_energies, decay_energies],
        bins=bins,
        stacked=True,
        color=["black", "#D81525"],
        label=["Direct (EVR)", r"From $^{273}$Ds $\alpha$ decay"],
        alpha=0.9,
        edgecolor="white",
        linewidth=0.5,
    )

    # 将原图 (c) 的斜线加在 8.915 MeV 事件对应的红色柱段上。
    hatched_energy = 8.915
    if np.any(np.isclose(decay_energies, hatched_energy, rtol=0, atol=1e-6)):
        bin_index = np.searchsorted(bins, hatched_energy, side="right") - 1
        decay_bars[bin_index].set_hatch("///")

    ax.legend(loc="upper left", frameon=False, fontsize=16)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, max(1, stacked_counts[-1].max()) / 0.75)
    ax.set_xlabel(r"$E_{\alpha}$ / MeV", fontsize=18)
    ax.set_ylabel(f"Counts / {bin_width_kev} keV", fontsize=18)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=14)
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight("bold")

    fig.savefig("HsData_Stacked.pdf", bbox_inches="tight")
    fig.savefig("HsData_Stacked.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
