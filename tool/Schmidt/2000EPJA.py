import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. 导入数据与基础计算 (修复文本导致无法取对数的问题)
# ==========================================
data = pd.read_csv('Caldata.txt', header=0, names=['T1'])
# 强制将数据转换为浮点数，忽略无法转换的非法字符(转为NaN)并剔除
t_exp = pd.to_numeric(data['T1'], errors='coerce').dropna().values
n = len(t_exp)
theta = np.log(t_exp)
theta_mean = np.mean(theta)
sigma_theta_exp = np.sqrt(np.sum((theta - theta_mean)**2) / n)

# ==========================================
# 2. Origin 风格图表全局设置 (彻底修复 Linux 负号与字体问题)
# ==========================================
plt.rcParams['font.family'] = 'serif'
# 首选 Times New Roman，如果 Linux 没有该字体，安全回退到几乎长一样的 DejaVu Serif
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif']
# 使用 'cm' (Computer Modern) 替代 'stix'，它是标准的 LaTeX 数学字体，绝不会丢负号
plt.rcParams['mathtext.fontset'] = 'cm' 
plt.rcParams['axes.unicode_minus'] = False # 强制非数学文本使用 ASCII 连字符作为负号

plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6

fig, ax = plt.subplots(figsize=(8, 6))

# ==========================================
# 3. 在对数时间尺度上绘制实验数据的直方图
# ==========================================
num_bins = 15  # 使用更细的网格

# 根据数据的实际跨度设定边界
theta_min = np.min(theta) - 0.2
theta_max = np.max(theta) + 0.2
bins_theta = np.linspace(theta_min, theta_max, num_bins + 1)
bin_width = bins_theta[1] - bins_theta[0]

counts, _ = np.histogram(theta, bins=bins_theta)
dn_dtheta_exp = counts / bin_width
bins_t = np.exp(bins_theta)

# 绘制柱状图 (bar)
ax.bar(bins_t[:-1], dn_dtheta_exp, width=np.diff(bins_t), align='edge', 
       color='lightgray', edgecolor='black', linewidth=1.5, label='Experimental data', zorder=2)

# ==========================================
# 4. 计算并绘制单源理论曲线
# ==========================================
ln_lambda = -theta_mean
t_smooth = np.logspace(np.log10(bins_t[0]*0.8), np.log10(bins_t[-1]*1.2), 500)
theta_smooth = np.log(t_smooth)
dn_dtheta_theory = n * np.exp(theta_smooth + ln_lambda) * np.exp(-np.exp(theta_smooth + ln_lambda))

ax.plot(t_smooth, dn_dtheta_theory, color='red', linewidth=2.5, label='Expected theoretical distribution', zorder=3)

# ==========================================
# 5. 坐标轴及图例美化
# ==========================================
ax.set_xscale('log')
ax.set_xlabel('t (ms)', fontweight='bold', fontsize=16) 
ax.set_ylabel(r'$dn / d\Theta$', fontweight='bold', fontsize=16)
ax.legend(loc='upper right', frameon=False, fontsize=12)

text_str = f'$n = {n}$\n$\\sigma_{{\\Theta_{{exp}}}} = {sigma_theta_exp:.3f}$'
ax.text(0.04, 0.96, text_str, transform=ax.transAxes, fontsize=16,
        verticalalignment='top', 
        bbox=dict(boxstyle='square,pad=0.5', facecolor='white', edgecolor='black', lw=1.5), zorder=4)

plt.tight_layout()
plt.savefig('Histogram_Plot_Fine.png', dpi=300, bbox_inches='tight')
print("图片已成功生成并保存！")
plt.show()