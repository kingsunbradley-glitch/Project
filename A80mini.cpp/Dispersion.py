import numpy as np
import matplotlib.pyplot as plt

# ================= Font and Format Settings =================
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix' 
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.linewidth'] = 1.5

# ================= Data Extraction =================
# vh configuration (6 points)
x1 = np.array([0.00865, 0.00288, 0, -0.00288, -0.00576, -0.01153]) * 100
y1 = np.array([4, 1, 0, -0.5, -3, -3.5])

# hv configuration (5 points)
x2 = np.array([-0.00583, -0.00292, 0, 0.00583, 0.01166]) * 100
y2 = np.array([-21, -12, 0, 15, 36])

# ================= Linear Fitting =================
# Fit 1 (vh)
coeffs1 = np.polyfit(x1, y1, 1)
slope1, intercept1 = coeffs1[0], coeffs1[1]

# Fit 2 (hv)
coeffs2 = np.polyfit(x2, y2, 1)
slope2, intercept2 = coeffs2[0], coeffs2[1]

# Generate points for fitting lines
all_x = np.concatenate((x1, x2))
x_fit = np.linspace(min(all_x) - 0.2, max(all_x) + 0.2, 100)
y_fit1 = slope1 * x_fit + intercept1
y_fit2 = slope2 * x_fit + intercept2

# ================= Plotting =================
# Set figure size to 16x8 inches. With dpi=100, output image will be exactly 1600x800 pixels.
fig, ax = plt.subplots(figsize=(16, 8))

# Plot vh data and fit (Black)
ax.scatter(x1, y1, color='black', marker='o', s=150, zorder=3, label=r'$\it{vh}$ Data')
ax.plot(x_fit, y_fit1, color='black', linestyle='-.', linewidth=2.5, zorder=2, label=r'$\it{vh}$ Linear Fit')

# Plot hv data and fit (Blue)
ax.scatter(x2, y2, color='blue', marker='^', s=150, zorder=3, label=r'$\it{hv}$ Data')
ax.plot(x_fit, y_fit2, color='blue', linestyle='-', linewidth=2.5, zorder=2, label=r'$\it{hv}$ Linear Fit')

# Axis formatting (ticks pointing inwards on all sides)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(width=1.5, length=8, labelsize=20)

# Set labels with bold fonts
ax.set_xlabel(r'$\Delta D_1$ (%)', fontweight='bold', fontsize=22)
ax.set_ylabel(r'$\Delta X$ (mm)', fontweight='bold', fontsize=22)

# Legend at upper left, no frame
ax.legend(loc='upper left', frameon=False, fontsize=20, markerscale=1.2)

# Adjust layout and save
plt.tight_layout()
plt.savefig('dispersion_plot.png', dpi=300)
# plt.show() # Uncomment if you want to preview in your local environment