import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# hatched histogram
# Hs decay time plot - final version
# -----------------------------------------------------------------
# --- 全局字体与线条设置 ---
# -----------------------------------------------------------------
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'

plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.5
plt.rcParams['xtick.major.size'] = 3.5
plt.rcParams['ytick.major.size'] = 3.5
plt.rcParams['xtick.minor.size'] = 2.0
plt.rcParams['ytick.minor.size'] = 2.0

# -----------------------------------------------------------------
# --- 辅助函数：理论曲线 ---
# -----------------------------------------------------------------
def theoretical_curve(t, lambda_est):
    # 对数坐标下的理论分布形状: t * lambda * exp(-lambda * t)
    return lambda_est * t * np.exp(-lambda_est * t)

# -----------------------------------------------------------------
# --- 269Hs 数据准备 ---
# -----------------------------------------------------------------
# 原始数据 (ms)
# IMP: 5.005 (普通), 9.363 (特殊，需打阴影)
imp_val_normal = np.array([5.005]) * 1000.0   
imp_val_special = np.array([9.363]) * 1000.0  

dubna_data = np.array([4.978]) * 1000.0
riken_data = np.array([14.2, 0.27, 36]) * 1000.0
gsi_data = np.array([22, 19.7]) * 1000.0

# 拟合用的完整数据集 (所有事件)
all_fit_data = np.concatenate([imp_val_normal, imp_val_special, dubna_data, riken_data, gsi_data])

# 堆叠直方图数据列表
# 顺序: IMP(普通), IMP(特殊), Dubna, RIKEN, GSI
data_ms_list = [imp_val_normal, imp_val_special, dubna_data, riken_data, gsi_data]

# 转换为 Log10 值
data_log_list = [np.log10(d) for d in data_ms_list]

# 样式定义
colors = ["#D81525", "#D81525", 'blue', 'green', "#B8860b"] # IMP 都是红色
hatches = ['', '\\\\', '', '', ''] # 第二个元素 (IMP特殊) 加斜线
labels = ['IMP', '', 'Dubna', 'RIKEN', 'GSI'] # 空标签避免图例重复

# -----------------------------------------------------------------
# --- 绘图逻辑 ---
# -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5))

# X轴范围 (Log10)
log_min, log_max = 2.0, 6.0
bins = np.linspace(log_min, log_max, 40) # 保持之前的分bin密度

# 1. 绘制直方图
patches = ax.hist(data_log_list, bins=bins, stacked=True, 
                  color=colors, edgecolor='black', linewidth=0.3, label=labels)

# 2. 应用阴影 (Hatching)
bar_containers = patches[2]
for k, container in enumerate(bar_containers):
    hatch_style = hatches[k]
    for rect in container:
        rect.set_hatch(hatch_style)
        rect.set_edgecolor('white')

# 3. 计算并绘制 Fit 曲线 (不带图例)
if all_fit_data.size > 0:
    t_bar = np.mean(all_fit_data)
    lambda_est = 1.0 / t_bar
    
    # 生成曲线的 x 坐标 (log scale)
    x_axis = np.linspace(log_min, log_max, 500)
    t_real = 10**x_axis 
    
    # 计算理论值
    y_curve_unscaled = theoretical_curve(t_real, lambda_est)
    
    # 缩放曲线高度以匹配直方图最大值
    # 先把所有数据展平算一下原本直方图的最高点
    all_log_data_flat = np.concatenate(data_log_list)
    counts, _ = np.histogram(all_log_data_flat, bins=bins)
    max_count = np.max(counts) if len(counts) > 0 else 1
    
    if y_curve_unscaled.max() > 0:
        scale_factor = max_count / y_curve_unscaled.max()
        y_scaled = y_curve_unscaled * scale_factor
        
        # 绘制曲线 (没有 label，这样就不会出现在图例里)
        ax.plot(x_axis, y_scaled, '-', linewidth=2, color='#6A5ACD', zorder=10)

# -----------------------------------------------------------------
# --- 坐标轴与装饰 ---
# -----------------------------------------------------------------
ax.set_xlim(log_min, log_max)
ax.set_xlabel(r"$\log_{10}[t (\mathrm{ms})]$", fontsize=18)
ax.set_ylabel("Event Counts", fontsize=18)

# 刻度设置
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.tick_params(direction='in', top=True, right=True, which='both', labelsize=14)

# 边框设置 (全实线)
ax.spines['left'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(True)
ax.spines['bottom'].set_visible(True)

# 图例处理 (过滤掉空标签，保留唯一项)
handles, lbls = ax.get_legend_handles_labels()
by_label = dict(zip(lbls, handles))
clean_handles = []
clean_labels = []
# 指定顺序
target_order = ['IMP', 'Dubna', 'RIKEN', 'GSI']
for l in target_order:
    if l in by_label:
        clean_handles.append(by_label[l])
        clean_labels.append(l)

# 图例放右上
ax.legend(clean_handles, clean_labels, frameon=False, fontsize=14, loc='upper right')

# 269Hs 标签放左上
ax.text(0.05, 0.9, r"$^{269}$Hs", transform=ax.transAxes, fontsize=22, ha='left')

plt.tight_layout()
plt.savefig('Hs_final_plot_v2.pdf', dpi=300)
plt.show()