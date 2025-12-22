import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker # <-- 导入模块

# --- 1. USER CONFIGURATION AREA ---
# User's new data:

# 'data_to_fit': The dataset used to *generate* the theoretical curve.
# This will be plotted as a BLACK histogram.
data_to_fit = np.array([0.107, 0.184, 0.11, 0.31, 0.52, 0.0399, 0.373])

# 'data_to_test': Another dataset you want to visually assess against the curve.
# This will be plotted as a LIGHT BLUE histogram.
data_to_test = np.array([8.319, 41.703])

# Time unit for plot labeling
time_unit = "ms" 

# --- 2. CALCULATION FUNCTIONS (based on the paper) ---
# (Functions are identical to the previous version)

def calculate_sigma_theta_exp(data):
    """
    Calculates the standard deviation of the logarithmic decay times (sigma_Theta_exp).
    [cite_start]This is the main test statistic in the paper. (See Equations 12 & 13) [cite: 132, 137, 168-170].
    """
    n = len(data)
    if n < 2:
        return 0, 0, n
    log_t = np.log(data)
    theta_bar_exp = np.mean(log_t)
    sum_sq_dev = np.sum((log_t - theta_bar_exp)**2)
    sigma_theta_exp = np.sqrt(sum_sq_dev / n)
    return sigma_theta_exp, theta_bar_exp, n

def get_limits_from_table1(n):
    """
    Hard-coded lookup for the confidence limits from Table 1 for a given n.
    """
    table_1 = {
        2: (0.04, 1.83), 3: (0.19, 1.91), 4: (0.31, 1.92), 5: (0.41, 1.90),
        6: (0.48, 1.89), 7: (0.52, 1.87), 8: (0.58, 1.85), 9: (0.62, 1.84),
        10: (0.65, 1.82), 20: (0.81, 1.71), 50: (0.98, 1.57), 100: (1.06, 1.49)
    }
    if n in table_1:
        return table_1[n]
    else:
        return None, None

# --- 3. STATISTICAL CALCULATIONS ---
# (Statistical output remains the same)

print("Analyzing datasets...\n")

# Analyze Dataset 1 ('data_to_fit')
sigma1, mean_log1, n1 = calculate_sigma_theta_exp(data_to_fit)
lower1, upper1 = get_limits_from_table1(n1)
print(f"--- Dataset 1 ('data_to_fit') ---")
print(f"Number of events (n): {n1}")
print(f"Sigma_Theta_exp (Eq. 12): {sigma1:.3f}")
if lower1:
    print(f"90% confidence limits from Table 1 for n={n1}: ({lower1}, {upper1})")
    if lower1 <= sigma1 <= upper1:
        print(" -> RESULT: Value is INSIDE the range. Data is consistent with a single decay.")
    else:
        print(" -> RESULT: Value is OUTSIDE the range. Data is not consistent.")

print("\n" + "-"*30 + "\n")

# Analyze Dataset 2 ('data_to_test')
sigma2, mean_log2, n2 = calculate_sigma_theta_exp(data_to_test)
lower2, upper2 = get_limits_from_table1(n2)
print(f"--- Dataset 2 ('data_to_test') ---")
print(f"Number of events (n): {n2}")
print(f"Sigma_Theta_exp (Eq. 12): {sigma2:.3f}")
if lower2:
    print(f"90% confidence limits from Table 1 for n={n2}: ({lower2}, {upper2})")
    if lower2 <= sigma2 <= upper2:
        print(" -> RESULT: Value is INSIDE the range. Data is consistent with a single decay.")
    else:
        print(" -> RESULT: Value is OUTSIDE the range. Data is not consistent.")

print("\n" + "-"*30 + "\n")

# --- Combined Dataset (Test) ---
data_combined = np.concatenate([data_to_fit, data_to_test])
sigma_comb, mean_log_comb, n_comb = calculate_sigma_theta_exp(data_combined)
lower_comb, upper_comb = get_limits_from_table1(n_comb)
print(f"--- Combined Dataset (Test) ---")
print(f"Number of events (n): {n_comb}")
print(f"Sigma_Theta_exp (Eq. 12): {sigma_comb:.3f}")
if lower_comb:
    print(f"90% confidence limits from Table 1 for n={n_comb}: ({lower_comb}, {upper_comb})")
    if lower_comb <= sigma_comb <= upper_comb:
        print(" -> RESULT: Value is INSIDE the range. The combined set IS consistent with a single decay.")
    else:
        print(" -> RESULT: Value is OUTSIDE the range. The combined set is NOT consistent.")

# --- 4. PREPARING THEORETICAL CURVE PLOTTING ---
# (This part is unchanged)
t_bar_1 = np.mean(data_to_fit)
lambda_est_1 = 1.0 / t_bar_1

def theoretical_curve(t, lambda_est):
    return lambda_est * t * np.exp(-lambda_est * t)

all_data = np.concatenate([data_to_fit, data_to_test])
if all_data.size == 0:
    all_data = np.array([0.1, 10])
t_max_data = all_data.max() 
user_requested_start_log = -3.0
t_axis = np.logspace(
    user_requested_start_log,
    np.log10(t_max_data) + 0.5,
    500
)
y_curve = theoretical_curve(t_axis, lambda_est_1)

# --- 5. CREATE THE PLOT (Modified for Histograms) ---
fig, ax = plt.subplots(figsize=(12, 7))

# --- Plot Histograms (Left Y-axis) ---
if all_data.size > 0:
    t_min_data = all_data.min()
    bins = np.logspace(
        np.log10(t_min_data) - 0.1, 
        np.log10(t_max_data) + 0.1, 
        num=16 # 16 bin edges = 15 bins
    )
else:
    bins = np.logspace(np.log10(0.1), np.log10(10), num=16)

if n1 > 0:
    ax.hist(data_to_fit, bins=bins, color='black', alpha=1.0, 
            label=f"'data_to_fit' (n={n1})", edgecolor='white')
if n2 > 0:
    ax.hist(data_to_test, bins=bins, color='lightblue', alpha=0.7, 
            label=f"'data_to_test' (n={n2})", edgecolor='gray')

# Format the histogram axis (ax)
ax.set_xscale('log')
ax.set_xlabel(f"Time t / {time_unit}")
ax.set_ylabel("Event Counts")
ax.set_ylim(bottom=0)
ax.set_xlim(t_axis.min(), t_axis.max()) 

# --- *** MODIFICATION HERE *** ---
# Force Y-axis to use integer ticks
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
# --- *** END OF MODIFICATION *** ---

ax.grid(True, which='both', axis='x', linestyle='--', linewidth=0.5)

# --- Plot Theoretical Curve (Right Y-axis) ---
ax2 = ax.twinx()
line, = ax2.plot(t_axis, y_curve, 'r--', 
                 label=f"Theoretical Curve (Eq 7, based on 'data_to_fit')")
ax2.set_ylabel("Theoretical Density (Arbitrary Units)", color='r')
ax2.tick_params(axis='y', labelcolor='r')
ax2.set_ylim(bottom=0)

# --- Combined Legend ---
handles, labels = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
ax.legend(handles + handles2, labels + labels2, loc='upper right')

ax.set_title("Logarithmic Decay-Time Distribution (Histogram vs. Theory)")
plt.tight_layout()

plt.savefig("decay_histogram_plot_integer_yaxis.png")
print(f"\nPlot 'decay_histogram_plot_integer_yaxis.png' has been saved.")

# Display the figure inline
plt.show()