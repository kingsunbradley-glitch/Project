import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
from scipy.stats import norm

# -----------------------------------------------------------------
# --- 设置全局字体为 STIX (类 Times New Roman) ---
# -----------------------------------------------------------------
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# --- 1. 绘图辅助函数 (来自 plot_4color_time.py) ---
# -----------------------------------------------------------------

def theoretical_curve(t, lambda_est):
    """计算理论衰变曲线"""
    # t 单位为 ms, lambda_est 单位为 ms^-1
    return lambda_est * t * np.exp(-lambda_est * t)

# -----------------------------------------------------------------
# --- 2. 集中定义 273Ds 的数据 ---
# -----------------------------------------------------------------

# --- 能量数据 ---
ds_en_present = np.array([11.165, 10.857])
ds_en_dubna = np.array([11.017, 10.929])
ds_en_riken = np.array([11.14, 11.15, 11.03])
ds_en_gsi = np.array([11.2, 11.08])
#hs_en_present = np.array([9.103, 8.915])
#hs_en_dubna = np.array([])
#hs_en_riken = np.array([9.17, 9.25, 9.15])
#hs_en_gsi = np.array([9.18, 9.23])
#hs_en_D = np.array([8.91,8.92,8.93,8.93,9.03,9.06,9.11,9.13,9.17,9.18])
all_hist_data_energy = [ds_en_present, ds_en_dubna, ds_en_riken, ds_en_gsi]

# --- 时间数据 (ms) ---
#hs_hist_present = np.array([5.005, 9.363]) * 1000.0
#hs_hist_dubna = np.array([4.978]) * 1000.0
#hs_hist_riken = np.array([14.2, 0.27, 36]) * 1000.0
#hs_hist_gsi = np.array([22, 19.7]) * 1000.0
#hs_hist_D = np.array([])
ds_hist_present = np.array([0.107, 8.319])
ds_hist_dubna = np.array([0.184, 41.703])
ds_hist_riken = np.array([0.52, 0.0399, 0.373])
ds_hist_gsi = np.array([0.11, 0.31])
all_hist_data_time = [ds_hist_present, ds_hist_dubna, ds_hist_riken, ds_hist_gsi]

# --- 拟合数据 (ms) ---
ds_fit_data = np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373])

#ds_fit_data = np.array([5.005, 4.978, 14.2, 0.27, 36, 22, 19.7]) * 1000.0

# --- 通用样式 ---
labels = ['Present', 'Dubna', 'RIKEN', 'GSI']
colors = ['black', 'none', 'gray', 'white']
hatches = ['', '//', '', '']  # 为能量数据添加一个空的图例样式

# -----------------------------------------------------------------
# --- 3. 创建单个图像，使用 GridSpec 实现合并效果 ---
# -----------------------------------------------------------------

from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(14, 7))
gs = GridSpec(1, 2, figure=fig, wspace=0, width_ratios=[1, 1])

ax1 = fig.add_subplot(gs[0, 0])  # 左侧：能量
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)  # 右侧：时间

fig.suptitle(r"Energy and Decay Time Distribution for $^{273}$Ds", fontsize=18, y=0.98)

# -----------------------------------------------------------------
# --- 面板 1: 能量 (ax1) ---
# -----------------------------------------------------------------

bins_energy = np.linspace(10.8, 11.25, 20) 

# --- 绘制能量直方图 ---
ax1.hist(all_hist_data_energy,
         bins=bins_energy,
         stacked=True,
         color=colors,
         hatch=hatches,
         edgecolor='black')

# --- 格式化 ax1 (能量) ---
ax1.set_xlim(10.801, 11.281) 
ax1.set_ylim(bottom=0, top=3.5) 
ax1.set_xlabel("Energy / MeV", fontsize=14)
ax1.set_ylabel("Event Counts", fontsize=14)

# --- 刻度 ---
ax1.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune='lower'))  
ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# --- 文本 ---
ax1.text(0.1, 0.8, r"$^{273}$Ds", transform=ax1.transAxes, fontsize=24, weight='normal')

# 保留能量值标注
ax1.text(10.857, 1.2, '10.857', fontsize=16, ha='center', va='bottom', color='red', rotation=90)
ax1.text(10.929, 1.2, '10.929', fontsize=16, ha='center', va='bottom', color='blue', rotation=90)

# -----------------------------------------------------------------
# --- 面板 2: 时间 (ax2) ---
# -----------------------------------------------------------------

# --- 坐标轴设置 ---
ax2.set_xscale('log')
ax2.set_xlim(10**-2.1, 10**2) 
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1)) 
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))

# 重新定义分箱以匹配新的 X 轴范围
bins_time = np.logspace(np.log10(ax2.get_xlim()[0]), np.log10(ax2.get_xlim()[1]), num=40)

# --- 绘制时间直方图 (带图例) ---
n_hist, _, _ = ax2.hist(all_hist_data_time,
                         bins=bins_time,
                         stacked=True,
                         label=labels,
                         color=colors,
                         hatch=hatches,
                         edgecolor='black')

# --- 计算并绘制拟合曲线 ---
if ds_fit_data.size > 0:
    t_bar_ms = np.mean(ds_fit_data)
    lambda_est = 1.0 / t_bar_ms
    
    t_axis = np.logspace(np.log10(ax2.get_xlim()[0]), np.log10(ax2.get_xlim()[1]), 500) 
    y_curve_unscaled = theoretical_curve(t_axis, lambda_est) 
    
    max_hist_count = np.max(n_hist) if n_hist.size > 0 else 1
    if max_hist_count == 0: max_hist_count = 1
        
    if y_curve_unscaled.max() > 0:
        y_scaled_curve = y_curve_unscaled * (max_hist_count / y_curve_unscaled.max())
        ax2.plot(t_axis, y_scaled_curve, '-', linewidth=2 ,color='#6A5ACD')

# --- 格式化 ax2 (时间) ---
ax2.set_xlabel("Time (ms)", fontsize=14)

# 保留时间值标注
ax2.text(8.319, 1.2, '8.319 ', fontsize=16, ha='center', va='bottom', color='red', rotation=90, fontweight='bold')
ax2.text(41.703, 1.2, '41.703 ', fontsize=16, ha='center', va='bottom', color='blue', rotation=90)

# --- 图例 ---
ax2.legend(frameon=False, fontsize=12, loc='upper right', bbox_to_anchor=(0.98, 0.98))

# -----------------------------------------------------------------
# --- 4. 应用 'Origin' 风格并创建 '拼接' 效果 ---
# -----------------------------------------------------------------

# --- 格式化 ax1 (左侧面板) ---
ax1.spines['top'].set_visible(True)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.spines['right'].set_visible(False)  # 隐藏右侧轴
ax1.tick_params(axis='both', which='major', direction='in', top=True, right=False, labelsize=12)
ax1.tick_params(axis='both', which='minor', direction='in', top=True, right=False)
ax1.minorticks_on()
    
# --- 格式化 ax2 (右侧面板) ---
ax2.spines['top'].set_visible(True)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(True)  # 显示左侧轴
ax2.spines['left'].set_linestyle('--')  # 左侧轴设为虚线（分隔线）
ax2.spines['left'].set_linewidth(1.5)  # 加粗虚线
ax2.spines['left'].set_color('gray')  # 设为灰色
ax2.spines['right'].set_visible(True)
ax2.tick_params(axis='both', which='major', direction='in', top=True, right=True, labelsize=12)
ax2.tick_params(axis='both', which='minor', direction='in', top=True, right=True)

# 隐藏 ax2 的 Y 轴刻度标签（因为共享Y轴）
ax2.tick_params(axis='y', labelleft=False, left=False, which='both')

# -----------------------------------------------------------------
# --- 5. 调整布局并保存 ---
# -----------------------------------------------------------------
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("Ds273_Energy_Time_Merged.png", dpi=300, bbox_inches='tight')
print("Plot 'Ds273_Energy_Time_Merged.png' has been saved.")
plt.show()