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
    return lambda_est * t * np.exp(-lambda_est * t)

def plot_panel(ax, hist_present, hist_dubna, hist_riken, hist_gsi, fit_data, nuclide_label, tau_label, show_legend=False):
    """
    绘制单个面板的函数
    (已修改为接受 4 组数据)
    """
    
    # --- 1. 定义坐标轴和分箱 ---
    ax.set_xscale('log')
    ax.set_xlim(10**-5, 10**3) 
    
    ax.set_ylim(0, 3.5)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1)) 
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    
    # [修改] 增加分bin数量, 让直方图更细 (num=40)
    bins = np.logspace(np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1]), num=40)

    # --- 2. 绘制堆叠直方图 (4 种样式) ---
    all_hist_data = [hist_present, hist_dubna, hist_riken, hist_gsi]
    labels = ['Present', 'Dubna', 'RIKEN', 'GSI']
    colors = ['none', 'black', 'gray', 'white']
    hatches = ['/', '', '', ''] 
    
    n_hist, _, _ = ax.hist(all_hist_data,
                             bins=bins,
                             stacked=True,
                             label=labels,
                             color=colors,
                             hatch=hatches,
                             edgecolor='black')
    
    # --- 3. 计算并绘制拟合曲线 ---
    if fit_data.size > 0:
        t_bar = np.mean(fit_data)
        lambda_est = 1.0 / t_bar
        
        t_axis = np.logspace(np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1]), 500)
        y_curve_unscaled = theoretical_curve(t_axis, lambda_est)
        
        max_hist_count = np.max(n_hist) if n_hist.size > 0 else 1
        if max_hist_count == 0: max_hist_count = 1
            
        if y_curve_unscaled.max() > 0:
            y_scaled_curve = y_curve_unscaled * (max_hist_count / max_hist_count) # [修正] 简单缩放
            
            # [修正] 确保曲线不会超过直方图最大值太多
            # 找到直方图在峰值附近的bin的高度
            peak_x = t_axis[np.argmax(y_curve_unscaled)]
            bin_index = np.searchsorted(bins, peak_x) - 1
            if bin_index < 0: bin_index = 0
            
            # n_hist[0] 是 Present, n_hist[1] Dubna...
            # 我们需要总高度
            total_hist_heights = np.array([h[bin_index] for h in n_hist]).sum()
            
            # 使用一个更可靠的缩放：让曲线的峰值匹配直方图的最高点
            y_scaled_curve = y_curve_unscaled * (max_hist_count / y_curve_unscaled.max())
            
            ax.plot(t_axis, y_scaled_curve, 'k-', linewidth=1.5)
        
    # --- 4. 添加文本标签 ---
    # [修改] 调大核素标签的字号 (fontsize=18)
    ax.text(0.1, 0.8, nuclide_label, transform=ax.transAxes, fontsize=18, weight='bold')
    
    # --- [修改] 按要求移除 tau 标签 ---
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
# -----------------------------------------------------------------

# --- 数据: 273Ds (来自 doublePlotDs1.py) ---
ds_hist_present = np.array([0.107, 8.319]) / 1000.0
ds_hist_dubna = np.array([0.184, 41.703]) / 1000.0
ds_hist_riken = np.array([0.52, 0.0399, 0.373]) / 1000.0
ds_hist_gsi = np.array([0.11, 0.31]) / 1000.0
ds_fit_data = np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373]) / 1000.0
ds_t_bar_s = np.mean(ds_fit_data)
ds_tau_label = f"τ = {ds_t_bar_s * 1000.0:.2f} ms" 

# --- 数据: 269Hs (来自 doublePlotHs.py) ---
hs_hist_present = np.array([5.005, 9.363])
hs_hist_dubna = np.array([4.978])
hs_hist_riken = np.array([14.2, 0.27, 36])
hs_hist_gsi = np.array([22, 19.7])
hs_fit_data = np.array([5.005, 4.978, 14.2, 0.27, 36, 22, 19.7]) 
hs_t_bar_s = np.mean(hs_fit_data)
hs_tau_label = f"τ = {hs_t_bar_s:.1f} s"

# --- 数据: 265Sg (来自 doublePlotSg.py) ---
sg_hist_present = np.array([6.339, 11.663])
sg_hist_dubna = np.array([9.076, 58.19])
sg_hist_riken = np.array([23, 79.9, 13.8])
sg_hist_gsi = np.array([18.8, 7.4])
sg_fit_data = np.array([6.339, 11.663, 23, 79.9, 13.8, 18.8, 7.4])
sg_t_bar_s = np.mean(sg_fit_data)
sg_tau_label = f"τ = {sg_t_bar_s:.1f} s"


# -----------------------------------------------------------------
# --- 3. 创建图像 ---
# -----------------------------------------------------------------

# 创建 3 行 1 列的子图, 共享 X 轴和 Y 轴
fig, (ax_ds, ax_hs, ax_sg) = plt.subplots(3, 1, figsize=(8, 10), sharex=True, sharey=True)

# --- 绘制面板 ---
# [修改] 尽管 tau_label 仍在传递, 但 plot_panel 函数内部已将其注释掉
plot_panel(ax_ds, 
           ds_hist_present, ds_hist_dubna, ds_hist_riken, ds_hist_gsi, 
           ds_fit_data, 
           r"$^{273}$Ds", ds_tau_label, 
           show_legend=True)

plot_panel(ax_hs, 
           hs_hist_present, hs_hist_dubna, hs_hist_riken, hs_hist_gsi, 
           hs_fit_data, 
           r"$^{269}$Hs", hs_tau_label)

plot_panel(ax_sg, 
           sg_hist_present, sg_hist_dubna, sg_hist_riken, sg_hist_gsi, 
           sg_fit_data, 
           r"$^{265}$Sg", sg_tau_label)

# --- 最终格式化 ---
fig.subplots_adjust(hspace=0)
ax_ds.tick_params(axis='x', labelbottom=False)
ax_hs.tick_params(axis='x', labelbottom=False)
ax_sg.set_xlabel("Time (s)", fontsize=14)
fig.text(0.04, 0.5, 'Number of events', va='center', rotation='vertical', fontsize=14)

# --- 保存 ---
plt.savefig("multi_panel_time_plot_v14_no_tau.pdf", dpi=300)
print("Plot 'multi_panel_time_plot_v14_no_tau.pdf' has been saved.")
plt.show()