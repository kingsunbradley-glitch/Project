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

# --- Combine all datasets (Removed Rf) ---
all_datasets = [data_ds, data_hs, data_sg]

# Styles
labels = ['IMP', 'Dubna', 'RIKEN', 'GSI']
colors = ["#D81525", '#126782', '#2A9BC4', '#A8D5E2']
hatches = ['', '', '', '', '']

# -----------------------------------------------------------------
# --- Plotting Main Logic ---
# -----------------------------------------------------------------
# Adjusted figsize for 3 rows (slightly shorter)
fig = plt.figure(figsize=(12, 10.5)) 
# Adjusted GridSpec to 3 rows
gs = GridSpec(3, 2, figure=fig, wspace=0, hspace=0, width_ratios=[1, 1])

# Global Axes Limits
energy_min, energy_max = 7.8, 12.2
time_min, time_max = 10**-2.7, 1e6

# Pre-calculate bins to ensure alignment
bins_energy = np.linspace(energy_min, energy_max, 40) 
bins_time = np.logspace(np.log10(time_min), np.log10(time_max), 40)

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
                     edgecolor='black',
                     label=labels)
    
    ax_left.set_xlim(energy_min, energy_max)
    ax_left.set_ylim(bottom=0)
    ylim = ax_left.get_ylim()
    # Ensure enough space for label
    ax_left.set_ylim(0, max(ylim[1], 1) * 1.3) 

    # Label Isotope Name
    ax_left.text(0.05, 0.82, data['name'], transform=ax_left.transAxes, fontsize=20)
    
    # Grid/Ticks
    ax_left.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax_left.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax_left.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    
    # Frame Styles
    ax_left.spines['top'].set_visible(True)
    ax_left.spines['right'].set_visible(False)
    ax_left.tick_params(direction='in', top=True, right=False, which='both', labelsize=12)
    
    # Hide X labels for top 2 rows (since there are 3 rows total now)
    if i < 2:
        ax_left.tick_params(labelbottom=False)
    else:
        ax_left.set_xlabel("Energy / MeV", fontsize=14)
    
    # --- Right Plot: Time ---
    ax_right.set_xscale('log')
    ax_right.set_xlim(time_min, time_max)
    
    has_time_data = any(len(arr) > 0 for arr in data['time'])
    
    n_hist = None
    if has_time_data:
        n_hist, _, _ = ax_right.hist(data['time'],
                                     bins=bins_time,
                                     stacked=True,
                                     color=colors,
                                     hatch=hatches,
                                     edgecolor='black',
                                     label=labels)
    
    # Fit Curve
    fit_data = data['fit']
    if fit_data.size > 0:
        t_bar_ms = np.mean(fit_data)
        lambda_est = 1.0 / t_bar_ms
        
        t_axis = np.logspace(np.log10(time_min), np.log10(time_max), 500)
        y_curve_unscaled = theoretical_curve(t_axis, lambda_est)
        
        # Scale
        if n_hist is not None:
             all_time_flat = np.concatenate([x for x in data['time'] if len(x)>0])
             counts, _ = np.histogram(all_time_flat, bins=bins_time)
             max_count = np.max(counts)
        else:
            max_count = 1
            
        if max_count == 0: max_count = 1
        
        if y_curve_unscaled.max() > 0:
            scale_factor = max_count / y_curve_unscaled.max()
            y_scaled_curve = y_curve_unscaled * scale_factor
            ax_right.plot(t_axis, y_scaled_curve, '-', linewidth=2, color='#6A5ACD')

            # Add Mean Life Text
            

    # Frame Styles
    ax_right.spines['top'].set_visible(True)
    ax_right.spines['left'].set_linestyle('--')
    ax_right.spines['left'].set_linewidth(1.5)
    ax_right.spines['left'].set_color('gray')
    ax_right.tick_params(direction='in', top=True, right=True, which='both', labelsize=12)
    ax_right.tick_params(axis='y', labelleft=False, left=False) # Hide Y text
    
    # Hide X labels for top 2 rows
    if i < 2:
        ax_right.tick_params(labelbottom=False)
    else:
        ax_right.set_xlabel("Time (ms)", fontsize=14)
        
    # Legend (Only on first plot)
    if i == 0:
        ax_right.legend(frameon=False, fontsize=10, loc='upper right')

# Global Y Label
fig.text(0.06, 0.5, 'Event Counts', va='center', rotation='vertical', fontsize=16)

plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0, wspace=0)

plt.savefig('combined_plot_3rows.png', dpi=300, bbox_inches='tight')
plt.show()