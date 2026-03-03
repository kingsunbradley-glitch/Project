import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ==========================================
# Beam optics style version:
# entrance drift + Qv + drift + Qh
# Positive particle along +z
# Qv: y focus / x defocus
# Qh: x focus / y defocus
# ==========================================

# ---------- 1) Beam parameters ----------
q = 1.602e-19
m = 1.672e-27
E_k_MeV = 10.0
E_k_J = E_k_MeV * 1e6 * q
v0 = np.sqrt(2.0 * E_k_J / m)

# Beamline layout (m)
drift0_end = 0.35
qv_end = 0.95
drift1_end = 1.35
qh_end = 1.95
z_end = 2.45

# Quadrupole parameters
G_qv = 5.0
G_qh = 5.0
aperture = 0.08

# Initial conditions: slight divergence in both planes
N_particles = 4
theta_x_mrad = np.array([12, -10, 16, -14])
theta_y_mrad = np.array([24, 18, -20, -15])

# Integrator
max_z = 2.5
t_max = z_end / v0 * 1.18
dt = 2.0e-11
steps = int(round(t_max / dt))

pos = np.zeros((N_particles, 3))
vel = np.zeros((N_particles, 3))

tx = theta_x_mrad * 1e-3
ty = theta_y_mrad * 1e-3
vz = v0 / np.sqrt(1.0 + tx**2 + ty**2)
vel[:, 0] = vz * tx
vel[:, 1] = vz * ty
vel[:, 2] = vz

sample_rate = 25
num_samples = (steps - 1) // sample_rate + 1
trajectory = np.zeros((num_samples, N_particles, 3))

# ---------- 2) Piecewise lattice fields ----------
def get_B_field_beamline(pos_vec):
    """
    Piecewise hard-edge beamline:
      Drift (0 -> drift0_end)
      Qv    (drift0_end -> qv_end)
      Drift (qv_end -> drift1_end)
      Qh    (drift1_end -> qh_end)
      Drift (qh_end -> z_end)
    """
    B = np.zeros_like(pos_vec)
    z = pos_vec[:, 2]
    x = pos_vec[:, 0]
    y = pos_vec[:, 1]

    in_qv = (z >= drift0_end) & (z <= qv_end)
    in_qh = (z >= drift1_end) & (z <= qh_end)

    # Qv: Bx=-G*y, By=-G*x -> y focus, x defocus
    B[in_qv, 0] = -G_qv * y[in_qv]
    B[in_qv, 1] = -G_qv * x[in_qv]

    # Qh: sign flipped -> x focus, y defocus
    B[in_qh, 0] = +G_qh * y[in_qh]
    B[in_qh, 1] = +G_qh * x[in_qh]

    return B


def get_accel(v_vec, p_vec):
    return (q / m) * np.cross(v_vec, get_B_field_beamline(p_vec))

# ---------- 3) RK4 ----------
print(f"开始计算 beamline (drift + Qv + drift + Qh)，共 {steps} 步...")
for i in range(steps):
    if i % sample_rate == 0:
        trajectory[i // sample_rate, :, :] = pos

    k1_v = get_accel(vel, pos)
    k1_p = vel

    k2_v = get_accel(vel + 0.5 * k1_v * dt, pos + 0.5 * k1_p * dt)
    k2_p = vel + 0.5 * k1_v * dt

    k3_v = get_accel(vel + 0.5 * k2_v * dt, pos + 0.5 * k2_p * dt)
    k3_p = vel + 0.5 * k2_v * dt

    k4_v = get_accel(vel + k3_v * dt, pos + k3_p * dt)
    k4_p = vel + k3_v * dt

    vel += (dt / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
    pos += (dt / 6.0) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p)

X = trajectory[:, :, 0]
Y = trajectory[:, :, 1]
Z = trajectory[:, :, 2]

# ---------- 4) Visualization ----------
fig = plt.figure(figsize=(22, 6.5))
fig.subplots_adjust(wspace=0.35, left=0.05, right=0.95, bottom=0.15)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
max_xy = 0.12

# Left: 3D
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('3D Trajectory ', fontsize=14, pad=20)
ax1.set_xlabel('Z axis (m)', labelpad=15)
ax1.set_ylabel('X axis (m)', labelpad=15)
ax1.set_zlabel('Y axis (m)', labelpad=15)
ax1.set_xlim(0, max_z)
ax1.set_ylim(-max_xy, max_xy)
ax1.set_zlim(-max_xy, max_xy)
ax1.scatter([0], [0], [0], color='k', s=50, label='Target Point')

# Draw two transparent square tubes
for zs, ze, label in [(drift0_end, qv_end, 'Qv'), (drift1_end, qh_end, 'Qh')]:
    z_tube = np.linspace(zs, ze, 40)
    y_tube = np.linspace(-aperture, aperture, 2)
    Zt, Yt = np.meshgrid(z_tube, y_tube)
    Xt_plus = np.full_like(Zt, aperture)
    Xt_minus = np.full_like(Zt, -aperture)
    ax1.plot_surface(Zt, Xt_plus, Yt, color='gray', alpha=0.12, linewidth=0)
    ax1.plot_surface(Zt, Xt_minus, Yt, color='gray', alpha=0.12, linewidth=0)
    x_tube = np.linspace(-aperture, aperture, 2)
    Zt2, Xt = np.meshgrid(z_tube, x_tube)
    Yt_plus = np.full_like(Zt2, aperture)
    Yt_minus = np.full_like(Zt2, -aperture)
    ax1.plot_surface(Zt2, Xt, Yt_plus, color='gray', alpha=0.12, linewidth=0)
    ax1.plot_surface(Zt2, Xt, Yt_minus, color='gray', alpha=0.12, linewidth=0)
    ax1.text((zs + ze) * 0.5, 0.09, 0.09, label, color='gray')

lines3d = [ax1.plot([], [], [], color=c, linewidth=2.0)[0] for c in colors]
points3d = [ax1.plot([], [], [], marker='o', color=c, markersize=5)[0] for c in colors]

# Middle: x-z
ax2 = fig.add_subplot(132)
ax2.set_title('Horizontal Plane: x-z', fontsize=14, pad=20)
ax2.set_xlabel('Z axis (m)')
ax2.set_ylabel('X axis (m)')
ax2.set_xlim(0, max_z)
ax2.set_ylim(-max_xy, max_xy)
ax2.grid(True, linestyle=':')
ax2.axvspan(drift0_end, qv_end, color='gray', alpha=0.15, label='Qv')
ax2.axvspan(drift1_end, qh_end, color='gray', alpha=0.28, label='Qh')
ax2.legend(loc='upper left')
lines_x = [ax2.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points_x = [ax2.plot([], [], marker='o', color=c, markersize=5)[0] for c in colors]

# Right: y-z
ax3 = fig.add_subplot(133)
ax3.set_title('Vertical Plane: y-z', fontsize=14, pad=20)
ax3.set_xlabel('Z axis (m)')
ax3.set_ylabel('Y axis (m)')
ax3.set_xlim(0, max_z)
ax3.set_ylim(-max_xy, max_xy)
ax3.grid(True, linestyle=':')
ax3.axvspan(drift0_end, qv_end, color='gray', alpha=0.15, label='Qv')
ax3.axvspan(drift1_end, qh_end, color='gray', alpha=0.28, label='Qh')
ax3.legend(loc='upper left')
lines_y = [ax3.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points_y = [ax3.plot([], [], marker='o', color=c, markersize=5)[0] for c in colors]

ax2.text(1.02, 0.095, 'Qv: x defocus', fontsize=10)
ax2.text(1.02, 0.082, 'Qh: x focus', fontsize=10)
ax3.text(1.02, 0.095, 'Qv: y focus', fontsize=10)
ax3.text(1.02, 0.082, 'Qh: y defocus', fontsize=10)

frames_total = num_samples


def update(frame):
    if frame % 50 == 0:
        print(f"渲染中... {frame}/{frames_total}")
    idx = max(1, frame)
    for n in range(N_particles):
        lines3d[n].set_data(Z[:idx, n], X[:idx, n])
        lines3d[n].set_3d_properties(Y[:idx, n])
        points3d[n].set_data(Z[idx:idx+1, n], X[idx:idx+1, n])
        points3d[n].set_3d_properties(Y[idx:idx+1, n])

        lines_x[n].set_data(Z[:idx, n], X[:idx, n])
        points_x[n].set_data(Z[idx:idx+1, n], X[idx:idx+1, n])

        lines_y[n].set_data(Z[:idx, n], Y[:idx, n])
        points_y[n].set_data(Z[idx:idx+1, n], Y[idx:idx+1, n])

    return lines3d + points3d + lines_x + points_x + lines_y + points_y


# preview
preview_idx = -1
for n in range(N_particles):
    lines3d[n].set_data(Z[:preview_idx, n], X[:preview_idx, n])
    lines3d[n].set_3d_properties(Y[:preview_idx, n])
    points3d[n].set_data(Z[preview_idx:preview_idx+1, n], X[preview_idx:preview_idx+1, n])
    points3d[n].set_3d_properties(Y[preview_idx:preview_idx+1, n])
    lines_x[n].set_data(Z[:preview_idx, n], X[:preview_idx, n])
    points_x[n].set_data(Z[preview_idx:preview_idx+1, n], X[preview_idx:preview_idx+1, n])
    lines_y[n].set_data(Z[:preview_idx, n], Y[:preview_idx, n])
    points_y[n].set_data(Z[preview_idx:preview_idx+1, n], Y[preview_idx:preview_idx+1, n])

fig.savefig('beamline_Qv_drift_Qh_style_preview.png', dpi=180)
plt.close(fig)

# fresh figure for GIF
fig = plt.figure(figsize=(22, 6.5))
fig.subplots_adjust(wspace=0.35, left=0.05, right=0.95, bottom=0.15)
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('3D Trajectory ', fontsize=14, pad=20)
ax1.set_xlabel('Z axis (m)', labelpad=15)
ax1.set_ylabel('X axis (m)', labelpad=15)
ax1.set_zlabel('Y axis (m)', labelpad=15)
ax1.set_xlim(0, max_z)
ax1.set_ylim(-max_xy, max_xy)
ax1.set_zlim(-max_xy, max_xy)
ax1.scatter([0], [0], [0], color='k', s=50, label='Target Point')
for zs, ze, label in [(drift0_end, qv_end, 'Qv'), (drift1_end, qh_end, 'Qh')]:
    z_tube = np.linspace(zs, ze, 40)
    y_tube = np.linspace(-aperture, aperture, 2)
    Zt, Yt = np.meshgrid(z_tube, y_tube)
    Xt_plus = np.full_like(Zt, aperture)
    Xt_minus = np.full_like(Zt, -aperture)
    ax1.plot_surface(Zt, Xt_plus, Yt, color='gray', alpha=0.12, linewidth=0)
    ax1.plot_surface(Zt, Xt_minus, Yt, color='gray', alpha=0.12, linewidth=0)
    x_tube = np.linspace(-aperture, aperture, 2)
    Zt2, Xt = np.meshgrid(z_tube, x_tube)
    Yt_plus = np.full_like(Zt2, aperture)
    Yt_minus = np.full_like(Zt2, -aperture)
    ax1.plot_surface(Zt2, Xt, Yt_plus, color='gray', alpha=0.12, linewidth=0)
    ax1.plot_surface(Zt2, Xt, Yt_minus, color='gray', alpha=0.12, linewidth=0)
    ax1.text((zs + ze) * 0.5, 0.09, 0.09, label, color='gray')
lines3d = [ax1.plot([], [], [], color=c, linewidth=2.0)[0] for c in colors]
points3d = [ax1.plot([], [], [], marker='o', color=c, markersize=5)[0] for c in colors]

ax2 = fig.add_subplot(132)
ax2.set_title('Horizontal Plane: x-z', fontsize=14, pad=20)
ax2.set_xlabel('Z axis (m)')
ax2.set_ylabel('X axis (m)')
ax2.set_xlim(0, max_z)
ax2.set_ylim(-max_xy, max_xy)
ax2.grid(True, linestyle=':')
ax2.axvspan(drift0_end, qv_end, color='gray', alpha=0.15, label='Qv')
ax2.axvspan(drift1_end, qh_end, color='gray', alpha=0.28, label='Qh')
ax2.legend(loc='upper left')
ax2.text(1.02, 0.095, 'Qv: x defocus', fontsize=10)
ax2.text(1.02, 0.082, 'Qh: x focus', fontsize=10)
lines_x = [ax2.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points_x = [ax2.plot([], [], marker='o', color=c, markersize=5)[0] for c in colors]

ax3 = fig.add_subplot(133)
ax3.set_title('Vertical Plane: y-z', fontsize=14, pad=20)
ax3.set_xlabel('Z axis (m)')
ax3.set_ylabel('Y axis (m)')
ax3.set_xlim(0, max_z)
ax3.set_ylim(-max_xy, max_xy)
ax3.grid(True, linestyle=':')
ax3.axvspan(drift0_end, qv_end, color='gray', alpha=0.15, label='Qv')
ax3.axvspan(drift1_end, qh_end, color='gray', alpha=0.28, label='Qh')
ax3.legend(loc='upper left')
ax3.text(1.02, 0.095, 'Qv: y focus', fontsize=10)
ax3.text(1.02, 0.082, 'Qh: y defocus', fontsize=10)
lines_y = [ax3.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points_y = [ax3.plot([], [], marker='o', color=c, markersize=5)[0] for c in colors]

ani = animation.FuncAnimation(fig, update, frames=frames_total, blit=False, interval=30)
print('正在生成 beamline GIF，请稍候...')
ani.save('beamline_Qv_drift_Qh_style.gif', writer='pillow', fps=18)
print('GIF 渲染完成：beamline_Qv_drift_Qh_style.gif')
