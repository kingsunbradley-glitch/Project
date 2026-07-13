import numpy as np
import matplotlib.pyplot as plt

# =========================
# Input data
# =========================
I = np.array([4, 5, 8, 10, 15], dtype=float)

PD1 = np.array([
    249000,
    323000,
    464000,
    570000,
    695000
], dtype=float)

# PD0 only uses the available data points.
# The missing values at I = 8 and I = 10 are ignored.
I_PD0 = np.array([4, 5, 15], dtype=float)

PD0 = np.array([
    11350,
    15000,
    49000
], dtype=float)


# =========================
# Linear fitting
# y = k*x + b
# =========================
k1, b1 = np.polyfit(I, PD1, 1)
k0, b0 = np.polyfit(I_PD0, PD0, 1)


def calculate_r2(x, y, k, b):
    """Calculate coefficient of determination R^2."""
    y_fit = k * x + b

    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    return 1.0 - ss_res / ss_tot


r2_PD1 = calculate_r2(I, PD1, k1, b1)
r2_PD0 = calculate_r2(I_PD0, PD0, k0, b0)


# =========================
# Print fitting results
# =========================
print("PD1 fitting result:")
print(f"PD1 = {k1:.6f} * I + {b1:.6f}")
print(f"R^2 = {r2_PD1:.6f}")

print()

print("PD0 fitting result:")
print(f"PD0 = {k0:.6f} * I + {b0:.6f}")
print(f"R^2 = {r2_PD0:.6f}")


# =========================
# Plot
# =========================
x_fit = np.linspace(3.5, 15.5, 300)

fig, ax1 = plt.subplots(figsize=(8, 6))

# PD1 data and fitting line
ax1.scatter(
    I,
    PD1,
    marker="o",
    s=60,
    label="PD1 data"
)

ax1.plot(
    x_fit,
    k1 * x_fit + b1,
    linewidth=2,
    label=(
        rf"PD1 fit: $y={k1:.2f}I{b1:+.2f}$"
        "\n"
        rf"$R^2={r2_PD1:.4f}$"
    )
)

ax1.set_xlabel("I", fontsize=14)
ax1.set_ylabel("PD1", fontsize=14)

ax1.tick_params(
    axis="both",
    which="both",
    direction="in",
    top=True,
    labelsize=12
)


# Create the second y axis for PD0
ax2 = ax1.twinx()

ax2.scatter(
    I_PD0,
    PD0,
    marker="s",
    s=60,
    label="PD0 data"
)

ax2.plot(
    x_fit,
    k0 * x_fit + b0,
    linestyle="--",
    linewidth=2,
    label=(
        rf"PD0 fit: $y={k0:.2f}I{b0:+.2f}$"
        "\n"
        rf"$R^2={r2_PD0:.5f}$"
    )
)

ax2.set_ylabel("PD0", fontsize=14)

ax2.tick_params(
    axis="y",
    which="both",
    direction="in",
    right=True,
    labelsize=12
)


# Combine the legends of the two axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax1.legend(
    lines1 + lines2,
    labels1 + labels2,
    loc="upper left",
    frameon=False,
    fontsize=10
)

fig.tight_layout()

plt.savefig(
    "PD_linear_fit.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()