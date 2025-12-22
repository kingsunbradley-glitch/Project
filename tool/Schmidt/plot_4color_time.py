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
# --- 1. 绘图辅助函数 ---
# -----------------------------------------------------------------

def theoretical_curve(t, lambda_est):
    """计算理论衰变曲线"""
    # t 单位为 ms, lambda_est 单位为 ms^-1
    return lambda_est * t * np.exp(-lambda_est * t)

def plot_panel(ax, hist_present, hist_dubna, hist_riken, hist_gsi, fit_data, nuclide_label, tau_label_unused, show_legend=False):
    """
    绘制单个面板的函数
    (已修改为接受 4 组数据)
    (已修改为使用 ms 单位)
    """
    
    # --- 1. 定义坐标轴和分箱 ---
    ax.set_xscale('log')
    # [修改] 坐标轴范围: 10^-5 s -> 10^-2 ms; 10^3 s -> 10^6 ms
    ax.set_xlim(10**-2, 10**6) 
    
    ax.set_ylim(0, 3.5)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1)) 
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    
    # [修改] 增加分bin数量, 让直方图更细 (num=60)
    bins = np.logspace(np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1]), num=60)

    # --- 2. 绘制堆叠直方图 (4 种样式) ---
    # [修改] all_hist_data 已经传入 ms 单位
    all_hist_data = [hist_present, hist_dubna, hist_riken, hist_gsi]
    labels = ['Present', 'Dubna', 'RIKEN', 'GSI']
    #colors = ['black','none',  'gray', 'white']
    #colors = ['black', '#2A9BC4', '#A8D5E2', '#E8F1F2']
    colors = ['#485098', '#5A96CC', '#A5C582', '#F9EF7C']
    # [修正] 将 '//' 改为 '/'
    hatches = ['', '', '', ''] 
    
    n_hist, _, _ = ax.hist(all_hist_data,
                             bins=bins,
                             stacked=True,
                             label=labels,
                             color=colors,
                             hatch=hatches,
                             edgecolor='black')
    
    # --- 3. 计算并绘制拟合曲线 ---
    # [修改] fit_data 已经传入 ms 单位
    if fit_data.size > 0:
        t_bar_ms = np.mean(fit_data) # t_bar 单位是 ms
        lambda_est = 1.0 / t_bar_ms  # lambda_est 单位是 ms^-1
        
        t_axis = np.logspace(np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1]), 500) # t_axis 单位是 ms
        y_curve_unscaled = theoretical_curve(t_axis, lambda_est) # ms 和 ms^-1 匹配
        
        max_hist_count = np.max(n_hist) if n_hist.size > 0 else 1
        if max_hist_count == 0: max_hist_count = 1
            
        if y_curve_unscaled.max() > 0:
            y_scaled_curve = y_curve_unscaled * (max_hist_count / y_curve_unscaled.max())
            ax.plot(t_axis, y_scaled_curve, 'k-', linewidth=1.5)
        
    # --- 4. 添加文本标签 ---
    # [修改] 调大核素标签的字号 (fontsize=18)
    ax.text(0.1, 0.8, nuclide_label, transform=ax.transAxes, fontsize=18, weight='bold')
    
    # --- [修改] 移除 tau 标签 (来自 v14 的修改) ---
    # ax.text(0.1, 0.65, tau_label, transform=ax.transAxes, fontsize=12)
    # --- [修改结束] ---

    # --- 5. 应用 'Origin' 风格 ---
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, labelsize=10)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
    
    # --- 6. 图例 ---
    if show_legend:
        # (图例微调)
        ax.legend(frameon=False, fontsize=10, loc='upper right', bbox_to_anchor=(0.98, 0.98))

# -----------------------------------------------------------------
# --- 2. 数据定义 (来自 3 个文件) ---
# [修改] 所有数据点转换为毫秒 (ms)
# -----------------------------------------------------------------

# --- 数据: 273Ds (来自 doublePlotDs1.py) ---
# [修改] 原始数据单位为 ms, 移除 / 1000.0
ds_hist_present = np.array([0.107, 8.319])
ds_hist_dubna = np.array([0.184, 41.703])
ds_hist_riken = np.array([0.52, 0.0399, 0.373])
ds_hist_gsi = np.array([0.11, 0.31])

ds_fit_data = np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373])
ds_t_bar_ms = np.mean(ds_fit_data) # 单位已是 ms
# [修改] 移除 * 1000.0, 因为 ds_t_bar_ms 已经是 ms
ds_tau_label = f"τ = {ds_t_bar_ms:.2f} ms" 

# --- 数据: 269Hs (来自 doublePlotHs.py) ---
# [修改] 原始数据单位为 s, 转换为 ms (乘以 1000.0)
hs_hist_present = np.array([5.005, 9.363]) * 1000.0
hs_hist_dubna = np.array([4.978]) * 1000.0
hs_hist_riken = np.array([14.2, 0.27, 36]) * 1000.0
hs_hist_gsi = np.array([22, 19.7]) * 1000.0

hs_fit_data = np.array([5.005, 4.978, 14.2, 0.27, 36, 22, 19.7]) * 1000.0
hs_t_bar_ms = np.mean(hs_fit_data) # 单位已是 ms
# [修改] 变量已是 ms, 仅修改标签字符串
hs_tau_label = f"τ = {hs_t_bar_ms:.1f} ms"

# --- 数据: 265Sg (来自 doublePlotSg.py) ---
# [修改] 原始数据单位为 s, 转换为 ms (乘以 1000.0)
sg_hist_present = np.array([6.339, 11.663]) * 1000.0
sg_hist_dubna = np.array([9.076, 58.19]) * 1000.0
sg_hist_riken = np.array([23, 79.9, 13.8]) * 1000.0
sg_hist_gsi = np.array([18.8, 7.4]) * 1000.0

sg_fit_data = np.array([6.339, 11.663, 23, 79.9, 13.8, 18.8, 7.4]) * 1000.0
sg_t_bar_ms = np.mean(sg_fit_data) # 单位已是 ms
# [修改] 变量已是 ms, 仅修改标签字符串
sg_tau_label = f"τ = {sg_t_bar_ms:.1f} ms"


# -----------------------------------------------------------------
# --- 3. 创建图像 ---
# -----------------------------------------------------------------

# 创建 3 行 1 列的子图, 共享 X 轴和 Y 轴
fig, (ax_ds, ax_hs, ax_sg) = plt.subplots(3, 1, figsize=(8, 10), sharex=True, sharey=True)

# --- 最终格式化 ---

# --- [新增] 添加 suptitle 并使用 y=0.93 下移 ---
fig.suptitle(" Decay Time Distribution ", fontsize=16, y=0.93)
# --- [修改结束] ---

# --- 绘制面板 ---
# 面板 1: Ds
plot_panel(ax_ds, 
           ds_hist_present, ds_hist_dubna, ds_hist_riken, ds_hist_gsi, 
           ds_fit_data, 
           r"$^{273}$Ds", ds_tau_label, 
           show_legend=True)

# --- [新增] 标注 273Ds 的两个时间值 ---
# [修改] X 坐标更新为 ms 单位
ax_ds.text(8.319, 1.2, '8.319 ', fontsize=12, ha='center', va='bottom' , color='blue' , rotation=90 ) #fontweight='bold'
ax_ds.text(41.703, 1.2, '41.703 ', fontsize=12, ha='center', va='bottom' , color='blue' , rotation=90 )
# --- [修改结束] ---

# 面板 2: Hs
plot_panel(ax_hs, 
           hs_hist_present, hs_hist_dubna, hs_hist_riken, hs_hist_gsi, 
           hs_fit_data, 
           r"$^{269}$Hs", hs_tau_label)

# 面板 3: Sg
plot_panel(ax_sg, 
           sg_hist_present, sg_hist_dubna, sg_hist_riken, sg_hist_gsi, 
           sg_fit_data, 
           r"$^{265}$Sg", sg_tau_label)



# 移除子图之间的垂直间距
fig.subplots_adjust(hspace=0)

# 隐藏顶部和中间图的X轴刻度标签
ax_ds.tick_params(axis='x', labelbottom=False)
ax_hs.tick_params(axis='x', labelbottom=False)

# 添加共享的X轴和Y轴标签
# [修改] X 轴标签单位改为 ms
ax_sg.set_xlabel("Time (ms)", fontsize=14)
fig.text(0.04, 0.5, 'Number of events', va='center', rotation='vertical', fontsize=14)

# --- 保存 ---
# [修改] 更改保存文件名以反映单位变化
plt.savefig("multi_panel_time_plot_ms_units.pdf", dpi=300)
print("Plot 'multi_panel_time_plot_ms_units.pdf' has been saved.")

# In a script, plt.show() would block execution. 
# In this environment, saving is sufficient.
# plt.show()