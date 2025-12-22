import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# -----------------------------------------------------------------
# --- Global Font Settings ---
# -----------------------------------------------------------------
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'

# -----------------------------------------------------------------
# --- Helper Functions ---
# -----------------------------------------------------------------
def theoretical_curve(t, lambda_est):
    return lambda_est * t * np.exp(-lambda_est * t)

# -----------------------------------------------------------------
# --- Data Definitions ---
# -----------------------------------------------------------------

# --- 1. 273Ds Data ---
ds_en_present = np.array([11.165, 10.857])
ds_en_dubna = np.array([11.017, 10.929])
ds_en_riken = np.array([11.14, 11.15, 11.03])
ds_en_gsi = np.array([11.2, 11.08])

ds_hist_present = np.array([0.107, 8.319])
ds_hist_dubna = np.array([0.184, 41.703])
ds_hist_riken = np.array([0.52, 0.0399, 0.373])
ds_hist_gsi = np.array([0.11, 0.31])

ds_fit_data = np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373])

data_ds = {
    "name": r"$^{273}$Ds",
    "energy": [ds_en_present, ds_en_dubna, ds_en_riken, ds_en_gsi],
    "time": [ds_hist_present, ds_hist_dubna, ds_hist_riken, ds_hist_gsi],
    "fit": ds_fit_data,
    "fmt": "{:.2f}"
}

# --- 2. 269Hs Data ---
hs_en_present = np.array([9.103, 8.915])
hs_en_dubna = np.array([])
hs_en_riken = np.array([9.17, 9.25, 9.15])
hs_en_gsi = np.array([9.18, 9.23])

# Converted to ms (* 1000)
hs_hist_present = np.array([5.005, 9.363]) * 1000.0
hs_hist_dubna = np.array([4.978]) * 1000.0
hs_hist_riken = np.array([14.2, 0.27, 36]) * 1000.0
hs_hist_gsi = np.array([22, 19.7]) * 1000.0

hs_fit_data = np.array([5.005, 4.978, 14.2, 0.27, 36, 22, 19.7]) * 1000.0

data_hs = {
    "name": r"$^{269}$Hs",
    "energy": [hs_en_present, hs_en_dubna, hs_en_riken, hs_en_gsi ],
    "time": [hs_hist_present, hs_hist_dubna, hs_hist_riken, hs_hist_gsi],
    "fit": hs_fit_data,
    "fmt": "{:.1f}"
}

# --- 3. 265Sg Data ---
sg_en_present = np.array([8.66, 8.671])
sg_en_dubna = np.array([8.397, 8.509])
sg_en_riken = np.array([8.71, 8.7, 8.66])
sg_en_gsi = np.array([])

# Converted to ms (* 1000)
sg_hist_present = np.array([6.339, 11.663]) * 1000.0
sg_hist_dubna = np.array([9.076, 58.19]) * 1000.0
sg_hist_riken = np.array([23, 79.9, 13.8]) * 1000.0
sg_hist_gsi = np.array([18.8, 7.4]) * 1000.0

sg_fit_data = np.array([6.339, 11.663, 23, 79.9, 13.8, 18.8, 7.4]) * 1000.0

data_sg = {
    "name": r"$^{265}$Sg",
    "energy": [sg_en_present, sg_en_dubna, sg_en_riken, sg_en_gsi ],
    "time": [sg_hist_present, sg_hist_dubna, sg_hist_riken, sg_hist_gsi],
    "fit": sg_fit_data,
    "fmt": "{:.1f}"
}

# --- 4. 261Rf Data ---
rf_en_present = np.array([])
rf_en_dubna = np.array([])
rf_en_riken = np.array([])
rf_en_gsi = np.array([8.52])

# Converted to ms (* 1000)
rf_hist_present = np.array([5.297, 7.24]) * 1000.0
rf_hist_dubna = np.array([3.0324, 0.3923]) * 1000.0
rf_hist_riken = np.array([2.97, 8.3, 3.73]) * 1000.0
rf_hist_gsi = np.array([4.7, 14.5]) * 1000.0

# Collect all time data for fitting
rf_fit_data = np.concatenate([
    x for x in [rf_hist_present, rf_hist_dubna, rf_hist_riken, rf_hist_gsi] 
    if len(x) > 0
])

data_rf = {
    "name": r"$^{261}$Rf",
    "energy": [rf_en_present, rf_en_dubna, rf_en_riken, rf_en_gsi],
    "time": [rf_hist_present, rf_hist_dubna, rf_hist_riken, rf_hist_gsi],
    "fit": rf_fit_data,
    "fmt": "{:.1f}"
}

# --- Combine all datasets ---
all_datasets = [data_ds, data_hs, data_sg, data_rf]

# Styles
labels = ['IMP', 'Dubna', 'RIKEN', 'GSI']
colors = ["#D81525", 'blue', 'green', "#B8860b"]
hatches = ['', '', '', '', '']

# -----------------------------------------------------------------
# --- Plotting Main Logic ---
# -----------------------------------------------------------------
fig = plt.figure(figsize=(12, 14)) 
num_rows = len(all_datasets)

# [修改 1] width_ratios=[2, 1] 让能量占 2/3，时间占 1/3
gs = GridSpec(num_rows, 2, figure=fig, wspace=0, hspace=0, width_ratios=[1.5, 1])

# Global Axes Limits
energy_min, energy_max = 7.8, 11.7

# [修改 2 & 3] 时间轴范围改为 log10 值 (-2.0 到 6.0)
log_time_min, log_time_max = -2.2, 6.0 

bins_energy = np.linspace(energy_min, energy_max, 400) 
# 使用线性bins来统计 log 后的数据
bins_time_log = np.linspace(log_time_min, log_time_max, 80)

first_ax_left = None

for i, data in enumerate(all_datasets):
    # Create Subplots
    if i == 0:
        ax_left = fig.add_subplot(gs[i, 0])
        ax_right = fig.add_subplot(gs[i, 1], sharey=ax_left)
        first_ax_left = ax_left 
    else:
        ax_left = fig.add_subplot(gs[i, 0], sharex=first_ax_left) 
        ax_right = fig.add_subplot(gs[i, 1], sharey=ax_left)

    # --- Left Plot: Energy ---
    has_energy_data = any(len(arr) > 0 for arr in data['energy'])
    
    if has_energy_data:
        ax_left.hist(data['energy'],
                     bins=bins_energy,
                     stacked=True,
                     color=colors,
                     hatch=hatches,
                     edgecolor='none', 
                     label=labels)
    
    ax_left.set_xlim(energy_min, energy_max)
    ax_left.text(0.05, 0.82, data['name'], transform=ax_left.transAxes, fontsize=20)
    ax_left.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax_left.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax_left.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))
    
    ax_left.spines['top'].set_visible(True)
    ax_left.spines['right'].set_visible(False)
    ax_left.tick_params(direction='in', top=True, right=False, which='both', labelsize=12)
    
    if i < num_rows - 1:
        ax_left.tick_params(labelbottom=False)
    else:
        # 正确写法
        ax_left.set_xlabel(r"$E_{\alpha}$ (MeV) ", fontsize=14)
    
    # --- Right Plot: Time (Log10 Axis) ---
    ax_right.set_xlim(log_time_min, log_time_max)
    
    # 将时间数据转换为 log10(T)
    log_time_data_list = []
    for arr in data['time']:
        if len(arr) > 0:
            valid_t = arr[arr > 0] # 确保大于0
            if len(valid_t) > 0:
                log_time_data_list.append(np.log10(valid_t))
            else:
                log_time_data_list.append(np.array([]))
        else:
            log_time_data_list.append(np.array([]))

    has_time_data = any(len(arr) > 0 for arr in log_time_data_list)

    n_hist = None
    if has_time_data:
        n_hist, _, _ = ax_right.hist(log_time_data_list,
                                     bins=bins_time_log, # 使用线性bins
                                     stacked=True,
                                     color=colors,
                                     hatch=hatches,
                                     edgecolor='none',
                                     label=labels)
    
    # Fit Curve (Modified for Log10 Axis)
    fit_data = data['fit']
    if fit_data.size > 0:
        t_bar_ms = np.mean(fit_data)
        lambda_est = 1.0 / t_bar_ms
        
        # 横坐标轴现在是 x = log10(t)
        x_axis = np.linspace(log_time_min, log_time_max, 500)
        t_real = 10**x_axis 
        
        y_curve_unscaled = theoretical_curve(t_real, lambda_est)
        
        # Scale
        if n_hist is not None:
             all_log_time_flat = np.concatenate([d for d in log_time_data_list if len(d)>0])
             counts, _ = np.histogram(all_log_time_flat, bins=bins_time_log)
             max_count = np.max(counts)
        else:
            max_count = 1
            
        if max_count == 0: max_count = 1
        
        if y_curve_unscaled.max() > 0:
            scale_factor = max_count / y_curve_unscaled.max()
            y_scaled_curve = y_curve_unscaled * scale_factor
            ax_right.plot(x_axis, y_scaled_curve, '-', linewidth=2, color='#6A5ACD')

    # Frame Styles
    ax_right.spines['top'].set_visible(True)
    ax_right.spines['left'].set_linestyle('--')
    ax_right.spines['left'].set_linewidth(1.5)
    ax_right.spines['left'].set_color('gray')
    ax_right.tick_params(direction='in', top=True, right=True, which='both', labelsize=12)
    ax_right.tick_params(axis='y', labelleft=False, left=False) 
    
    if i < num_rows - 1:
        ax_right.tick_params(labelbottom=False)
    else:
        # [修改 5] 更改标签
        ax_right.set_xlabel(r"$\log_{10}[t (\mathrm{ms})]$", fontsize=14)
        
    if i == 0:
        ax_right.legend(frameon=False, fontsize=10, loc='upper right')

    # -------------------------------------------------------------
    # Set Y-Limits (with 1.6 factor)
    # -------------------------------------------------------------
    ax_left.set_ylim(bottom=0)
    current_ylim = ax_left.get_ylim()
    ax_left.set_ylim(0, max(current_ylim[1], 1) * 1.6)

fig.text(0.06, 0.5, 'Counts', va='center', rotation='vertical', fontsize=16)

plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0, wspace=0)

plt.savefig('combined_plot_log10_axis_rf_fixed.png', dpi=300, bbox_inches='tight')
plt.show()