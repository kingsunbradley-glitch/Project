# multPlot.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import norm

# ======================================================
#  数据整合自 doublePlotDs, doublePlotHs, doublePlotSg
# ======================================================

datasets = {
    "Ds": {
        "time_fit": np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373]),
        "time_test": np.array([8.319, 41.703]),
        "energy_fit": np.array([11.14, 11.15, 11.03, 11.2, 11.08, 11.165, 11.017]),
        "energy_test": np.array([10.857, 10.929]),
        "time_unit": "ms",
    },
    "Hs": {
        "time_fit": np.array([5.005, 4.978, 14.2, 0.27, 36, 9.18, 9.23]) * 1e3,  # s→ms
        "time_test": np.array([9.363]) * 1e3,
        "energy_fit": np.array([9.103, 9.17, 9.25, 9.15, 9.18, 9.23]),
        "energy_test": np.array([8.915]),
        "time_unit": "ms",
    },
    "Sg": {
        "time_fit": np.array([6.339, 11.663, 23, 79.9, 13.8, 18.8, 7.4]) * 1e3,  # s→ms
        "time_test": np.array([9.076, 58.19]) * 1e3,
        "energy_fit": np.array([8.66, 8.671, 8.71, 8.7, 8.66]),
        "energy_test": np.array([8.397, 8.509]),
        "time_unit": "ms",
    }
}

# ======================================================
#  绘图参数
# ======================================================

# 提取整体横轴范围
all_energies = np.concatenate([v["energy_fit"] for v in datasets.values()] + [v["energy_test"] for v in datasets.values()])
energy_min, energy_max = all_energies.min() - 0.1, all_energies.max() + 0.1

all_times = np.concatenate([v["time_fit"] for v in datasets.values()] + [v["time_test"] for v in datasets.values()])
time_min, time_max = 0.01, 1e5  # 用户要求固定范围

# 创建画布：3行2列
fig, axes = plt.subplots(3, 2, figsize=(12, 15), sharex='col', sharey='col')
fig.subplots_adjust(hspace=0, wspace=0)

# ======================================================
#  绘制每一行
# ======================================================

for i, (nuclide, data) in enumerate(datasets.items()):
    ax_energy = axes[i, 0]
    ax_time = axes[i, 1]

    # -------- 能量分布 (左) --------
    energy_data = np.concatenate([data["energy_fit"], data["energy_test"]])
    bins_energy = np.linspace(energy_min, energy_max, 15)
    mu, sigma = np.mean(data["energy_fit"]), np.std(data["energy_fit"])
    x_fit = np.linspace(energy_min, energy_max, 200)
    y_pdf = norm.pdf(x_fit, mu, sigma)
    y_pdf_scaled = y_pdf * len(data["energy_fit"]) * (bins_energy[1] - bins_energy[0])

    # 各实验分组 (按照原 hatch 配色)
    if nuclide == "Ds":
        groups = {
            "IMP": np.array([11.165, 10.857]),
            "Dubna": np.array([11.017, 10.929]),
            "RIKEN": np.array([11.14, 11.15, 11.03]),
            "GSI": np.array([11.2, 11.08])
        }
    elif nuclide == "Hs":
        groups = {
            "IMP": np.array([9.103, 8.915]),
            "Dubna": np.array([]),
            "RIKEN": np.array([9.17, 9.25, 9.15]),
            "GSI": np.array([9.18, 9.23])
        }
    else:  # Sg
        groups = {
            "IMP": np.array([8.66, 8.671]),
            "Dubna": np.array([8.397, 8.509]),
            "RIKEN": np.array([8.71, 8.7, 8.66]),
            "GSI": np.array([])
        }

    colors = {'IMP': 'none', 'Dubna': 'black', 'RIKEN': 'gray', 'GSI': 'white'}
    hatches = {'IMP': '//', 'Dubna': '', 'RIKEN': '', 'GSI': ''}

    valid_groups = [g for g in groups.values() if g.size > 0]
    labels = [k for k, v in groups.items() if v.size > 0]
    cols = [colors[k] for k, v in groups.items() if v.size > 0]
    hts = [hatches[k] for k, v in groups.items() if v.size > 0]

    if valid_groups:
        ax_energy.hist(valid_groups, bins=bins_energy, stacked=True, label=labels,
                       color=cols, hatch=hts, edgecolor='black')
    ax_energy.plot(x_fit, y_pdf_scaled, 'r-', linewidth=1.5)
    ax_energy.legend(frameon=False, fontsize=9)
    ax_energy.set_xlim(energy_min, energy_max)
    ax_energy.set_ylim(bottom=0)
    ax_energy.tick_params(direction='in', top=True, right=True)
    ax_energy.minorticks_on()
    ax_energy.set_ylabel(f"{nuclide}", fontsize=13, rotation=0, labelpad=25, loc='center')

    # -------- 时间分布 (右) --------
    time_data = np.concatenate([data["time_fit"], data["time_test"]])
    t_bar = np.mean(data["time_fit"])
    lam = 1.0 / t_bar
    t_axis = np.logspace(np.log10(time_min), np.log10(time_max), 400)
    y_curve = lam * t_axis * np.exp(-lam * t_axis)

    # 直方图
    bins_time = np.logspace(np.log10(time_min), np.log10(time_max), 16)
    ax_time.hist(time_data, bins=bins_time, color='gray', edgecolor='black')
    # 缩放拟合曲线
    if len(time_data) > 0:
        hist_counts, _ = np.histogram(time_data, bins=bins_time)
        if hist_counts.max() > 0:
            y_curve_scaled = y_curve * (hist_counts.max() / y_curve.max())
            ax_time.plot(t_axis, y_curve_scaled, 'b-', linewidth=2)
    ax_time.set_xscale('log')
    ax_time.set_xlim(time_min, time_max)
    ax_time.set_ylim(bottom=0)
    ax_time.tick_params(direction='in', top=True, right=True)
    ax_time.minorticks_on()

# ======================================================
#  添加底部总标题
# ======================================================
axes[-1, 0].set_xlabel("Energy / MeV", fontsize=13)
axes[-1, 1].set_xlabel("Time t / ms", fontsize=13)
fig.text(0.25, 0.04, "Energy Distribution (MeV)", fontsize=15, ha='center')
fig.text(0.75, 0.04, "Decay-Time Distribution", fontsize=15, ha='center')

# ======================================================
#  保存图像
# ======================================================
plt.savefig("multPlot_combined.png", dpi=300, bbox_inches='tight')
plt.show()
print("✅ multPlot_combined.png 已生成。")
