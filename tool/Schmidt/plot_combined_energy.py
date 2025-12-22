import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
from scipy.stats import norm

# -----------------------------------------------------------------
# --- [新增] 设置全局字体为 STIX (类 Times New Roman) ---
# -----------------------------------------------------------------
# 设置字体家族为 'serif' (衬线字体)
plt.rcParams['font.family'] = 'serif'
# 设置 'serif' 家族的首选字体列表 (STIX 是一种专业学术字体)
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
# [关键] 设置数学字体 ( $...$ ) 使用 'stix' 字体集, 以便与常规文本匹配
plt.rcParams['mathtext.fontset'] = 'stix'
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# --- 1. 数据定义 (来自 3 个文件, 0A-3 部分) ---
# (这部分与前一个脚本完全相同)
# -----------------------------------------------------------------

# --- 273Ds (来自 doublePlotDs1.py) ---
ds_en_present = np.array([11.165, 10.857])
ds_en_previous = np.concatenate([
    np.array([11.017, 10.929]), # Dubna
    np.array([11.14, 11.15, 11.03]), # RIKEN
    np.array([11.2, 11.08])  # GSI
])

# --- 269Hs (来自 doublePlotHs.py) ---
hs_en_present = np.array([9.103, 8.915])
hs_en_previous = np.concatenate([
    np.array([]), # Dubna
    np.array([9.17, 9.25, 9.15]), # RIKEN
    np.array([9.18, 9.23]) # GSI
])

# --- 265Sg (来自 doublePlotSg.py) ---
sg_en_present = np.array([8.66, 8.671])
sg_en_previous = np.concatenate([
    np.array([8.397, 8.509]), # Dubna
    np.array([8.71, 8.7, 8.66]), # RIKEN
    np.array([]) # GSI
])

# --- 合并所有数据 ---
all_energy_present = np.concatenate([ds_en_present, hs_en_present, sg_en_present])
all_energy_previous = np.concatenate([ds_en_previous, hs_en_previous, sg_en_previous])


# -----------------------------------------------------------------
# --- 2. 创建带折线的图像 (两个子图) ---
# -----------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
fig.subplots_adjust(wspace=0.05)

# --- 定义两个区域的分箱 (Binning) ---
bins1 = np.linspace(8.2, 9.5, 25)  # 左侧图 (Sg, Hs)
bins2 = np.linspace(10.8, 11.25, 20) # 右侧图 (Ds) 

# --- 准备直方图数据 ---
all_hist_data = [all_energy_present, all_energy_previous]
labels = ['Present', 'Previous']
colors = ['none', 'gray'] 
hatches = ['//', '']       

# --- 绘制直方图 ---
ax1.hist(all_hist_data,
         bins=bins1,
         stacked=True,
         color=colors,
         hatch=hatches,
         edgecolor='black')

ax2.hist(all_hist_data,
         bins=bins2,
         stacked=True,
         label=labels,
         color=colors,
         hatch=hatches,
         edgecolor='black')

# --- 设置两个子图的 X 轴范围 ---
ax1.set_xlim(8.2, 9.5)  # 左侧图范围
ax2.set_xlim(10.801, 11.3) # 将起始点向右移动极小的量 (0.001) 来 "裁切" 掉第一个 bin 的边缘线

# 调整 ax2 的 X 轴刻度, 移除 '10.8'
ax2.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune='lower'))


# -----------------------------------------------------------------
# --- 3. 格式化和标注 ---
# -----------------------------------------------------------------

# --- 格式化 ax1 (左侧图) ---
ax1.set_ylabel("Event Counts", fontsize=14)
ax1.set_ylim(bottom=0, top=3.5) 
ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax2.legend(frameon=False, fontsize=12, loc='upper right', bbox_to_anchor=(0.98, 0.98))

# --- 添加共享的 X 轴标签 ---
fig.text(0.5, 0.02, f"Energy / MeV", ha='center', va='center', fontsize=14)
fig.suptitle(" Energy Distribution ", fontsize=16)

# --- 添加两条竖直虚线 ---
ax1.axvline(x=8.8, color='k', linestyle='--', linewidth=1)
ax1.axvline(x=9.4, color='k', linestyle='--', linewidth=1)

# --- 添加 Tex 文本 ---
ax1.text(8.4, 2.5, r"$^{265}$Sg", fontsize=24, ha='center', va='bottom')
ax1.text(9.0, 2.5, r"$^{269}$Hs", fontsize=24, ha='center', va='bottom')
ax2.text(10.9, 2.5, r"$^{273}$Ds", fontsize=24, ha='center', va='bottom')

ax1.text(8.940, 1.2, '8.915', fontsize=18, ha='center', va='bottom' , color='blue' , rotation=90)
ax2.text(10.857, 1.2, '10.857', fontsize=18, ha='center', va='bottom' ,color='blue' , rotation=90)
ax2.text(10.929, 1.2, '10.929', fontsize=18, ha='center', va='bottom' , color='blue' , rotation=90)


# -----------------------------------------------------------------
# --- 4. 'Origin' 风格并处理折线 ---
# -----------------------------------------------------------------

# 绘制折线标记 ( // )
d = 0.015  # 标记的大小
kwargs = dict(transform=ax1.transAxes, color='black', clip_on=False, linewidth=1)
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs) # ax1 右上
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)      # ax1 右下

kwargs.update(transform=ax2.transAxes)  # 切换到 ax2 的坐标
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs) # ax2 左上
ax2.plot((-d, +d), (-d, +d), **kwargs)      # ax2 左下

# --- 对两个子图应用 'Origin' 风格 ---
for ax in [ax1, ax2]:
    ax.spines['top'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['bottom'].set_color('black')
    
    # 设置刻度线: 方向朝内, 顶部显示
    ax.tick_params(axis='both', which='major', direction='in', top=True, labelsize=10)
    ax.tick_params(axis='both', which='minor', direction='in', top=True)
    ax.minorticks_on()

# --- 特殊处理折叠轴的轴脊和刻度 ---
ax1.spines['right'].set_visible(False) 
ax1.tick_params(axis='y', right=False) 
ax1.tick_params(axis='x', top=True)    
ax1.tick_params(axis='y', left=True)   

ax2.spines['left'].set_visible(False)  
ax2.spines['right'].set_visible(True)  
# [修改] 明确隐藏 ax2 的左侧 *所有* 刻度线 (主刻度 + 次刻度)
ax2.tick_params(axis='y', which='both', left=False) 
ax2.tick_params(axis='y', right=True)  
ax2.tick_params(axis='x', top=True)    


# --- 5. 保存 ---
plt.savefig("combined_energy_plot_broken_axis_v6.png", dpi=300)
print("Plot 'combined_energy_plot_broken_axis_v6.png' has been saved.")
plt.show()