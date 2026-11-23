#!/usr/bin/env python3
"""沿用 WangECC3 版式绘制 DNS_Li 的 4n、5n 激发函数。"""

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from WangECC3.plot_wang_ecc3 import draw_excitation_functions
import matplotlib.pyplot as plt


def main():
    data = np.genfromtxt(HERE / "excitation_functions_2n_5n.csv", delimiter=",", names=True)
    fig = draw_excitation_functions(
        data["E_star_MeV"], data["sigma_4n_pb"], data["sigma_5n_pb"],
        label_positions=((38, 0.045), (55.2, 1.05)),
    )
    for extension in ("pdf", "svg", "png"):
        fig.savefig(HERE / f"excitation_functions_dns_li.{extension}", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
