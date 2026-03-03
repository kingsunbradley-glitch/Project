import numpy as np
import matplotlib.pyplot as plt

# ================= 负号和字体修复 =================
plt.rcParams['axes.unicode_minus'] = False 
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix' 
plt.rcParams['font.size'] = 14  # 整体基础字体调大
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1.5

# ================= 数据提取和转换 =================
# vh构型
x1 = np.array([0.00865, 0.00288, 0, -0.00288, -0.00576, -0.01153]) * 100
y1 = np.array([4, 1, 0, -0.5, -3, -3.5])

# hv构型
x2 = np.array([-0.00583, -0.00292, 0, 0.00583, 0.01166]) * 100
y2 = np.array([-21, -12, 0, 15, 36])

# ================= 拟合计算 =================
# 计算拟合1 (vh)
coeffs1 = np.polyfit(x1, y1, 1)
slope1, intercept1 = coeffs1[0], coeffs1[1]

# 计算拟合2 (hv)
coeffs2 = np.polyfit(x2, y2, 1)
slope2, intercept2 = coeffs2[0], coeffs2[1]

# 生成用于绘制拟合直线的通用 x 轴点
all_x = np.concatenate((x1, x2))
x_fit = np.linspace(min(all_x) - 0.2, max(all_x) + 0.2, 100)
y_fit1 = slope1 * x_fit + intercept1
y_fit2 = slope2 * x_fit + intercept2

# ================= 创建单一画布绘图 =================
fig, ax = plt.subplots(figsize=(9, 7))

# 绘制vh数据及拟合线 (改：黑色数据点，黑色点划线)
ax.scatter(x1, y1, color='black', marker='o', s=80, zorder=3, label=r'$\it{vh}$ Data')
ax.plot(x_fit, y_fit1, color='black', linestyle='-.', linewidth=2, zorder=2, label=r'$\it{vh}$  Linear Fit')

# 绘制hv数据及拟合线 (改：蓝色数据点，蓝色实线)
ax.scatter(x2, y2, color='blue', marker='^', s=80, zorder=3, label=r'$\it{hv}$ Data')
ax.plot(x_fit, y_fit2, color='blue', linestyle='-', linewidth=2, zorder=2, label=r'$\it{hv}$  Linear Fit')

# 坐标轴格式化 (全边框刻度向内)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(width=1.5, length=6,labelsize=24)

# 设置标签
ax.set_xlabel(r'$\Delta D_1$ (%)', fontweight='bold', fontsize=16)
ax.set_ylabel(r'$\Delta X$ (mm)', fontweight='bold', fontsize=16)

# 图例放置在左上角 (改：去除了边框，并大幅度放大了图例字体和图标尺寸)
ax.legend(loc='upper left', frameon=False, fontsize=15, markerscale=1.2)

# (注：已删除生成公式文本框的代码)

plt.tight_layout()
plt.savefig('dispersion_plot.png', dpi=600)
plt.show()
