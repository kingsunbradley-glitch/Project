import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ============================================================
# 0. Shared visual settings
# ============================================================
COLORS = ['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e']
LABELS = ['P1: +X (1°)', 'P2: +Y (6°)', 'P3: -X (11°)', 'P4: -Y (16°)']
THETA_EMIT_DEG = np.array([1, 6, 11, 16])
PHI_EMIT_DEG = np.array([0, 90, 180, 270])
N_PARTICLES = 4

# Use the quadrupole figure range as the common XY range for every panel.
MAX_Z = 1.35
MAX_R_XY = 0.3

# Shared particle constants
q = 1.602e-19
m = 1.672e-27
E_k_MeV = 10.0
E_k_J = E_k_MeV * 1e6 * 1.602e-19
v0 = np.sqrt(2 * E_k_J / m)


def initial_state():
    pos = np.zeros((N_PARTICLES, 3))
    vel = np.zeros((N_PARTICLES, 3))
    theta = np.deg2rad(THETA_EMIT_DEG)
    phi = np.deg2rad(PHI_EMIT_DEG)
    vel[:, 0] = v0 * np.sin(theta) * np.cos(phi)
    vel[:, 1] = v0 * np.sin(theta) * np.sin(phi)
    vel[:, 2] = v0 * np.cos(theta)
    return pos, vel


# ============================================================
# 1. Solenoid simulation
# ============================================================
B_0 = 2.5
SOL_A = 0.20
SOL_Z1 = 0.3
SOL_Z2 = 1.3
K_drag_visual = 3e29
SOL_DT = 1e-11
SOL_TMAX = 5e-8
SOL_STEPS = int(round(SOL_TMAX / SOL_DT))
SOL_SAMPLE_RATE = 10


def get_B_field_solenoid(pos_vec):
    B_out = np.zeros_like(pos_vec)
    x = pos_vec[:, 0]
    y = pos_vec[:, 1]
    z = pos_vec[:, 2]
    r = np.sqrt(x**2 + y**2)
    r_safe = np.where(r == 0, 1e-10, r)

    term1 = (z - SOL_Z1) / np.sqrt((z - SOL_Z1) ** 2 + SOL_A**2)
    term2 = (z - SOL_Z2) / np.sqrt((z - SOL_Z2) ** 2 + SOL_A**2)
    Bz = (B_0 / 2.0) * (term1 - term2)

    dterm1 = (SOL_A**2) / ((z - SOL_Z1) ** 2 + SOL_A**2) ** 1.5
    dterm2 = (SOL_A**2) / ((z - SOL_Z2) ** 2 + SOL_A**2) ** 1.5
    dBz_dz = (B_0 / 2.0) * (dterm1 - dterm2)
    Br = -(r / 2.0) * dBz_dz

    B_out[:, 0] = Br * (x / r_safe)
    B_out[:, 1] = Br * (y / r_safe)
    B_out[:, 2] = Bz
    return B_out


def get_accel_solenoid(v_vec, p_vec):
    B_vec = get_B_field_solenoid(p_vec)
    a_L = (q / m) * np.cross(v_vec, B_vec)

    v_mag = np.linalg.norm(v_vec, axis=1, keepdims=True)
    v_mag[v_mag == 0] = 1e-10
    a_drag = -(K_drag_visual / (v_mag**2)) * (v_vec / v_mag)
    return a_L + a_drag


def simulate_solenoid():
    pos, vel = initial_state()
    traj = [pos.copy()]

    print(f"[Solenoid] integrating until z reaches {SOL_Z2} m...")
    for i in range(SOL_STEPS):
        k1_v = get_accel_solenoid(vel, pos)
        k1_p = vel
        k2_v = get_accel_solenoid(vel + 0.5 * k1_v * SOL_DT, pos + 0.5 * k1_p * SOL_DT)
        k2_p = vel + 0.5 * k1_v * SOL_DT
        k3_v = get_accel_solenoid(vel + 0.5 * k2_v * SOL_DT, pos + 0.5 * k2_p * SOL_DT)
        k3_p = vel + 0.5 * k2_v * SOL_DT
        k4_v = get_accel_solenoid(vel + k3_v * SOL_DT, pos + k3_p * SOL_DT)
        k4_p = vel + k3_v * SOL_DT

        vel += (SOL_DT / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        pos += (SOL_DT / 6.0) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p)

        if i % SOL_SAMPLE_RATE == 0:
            traj.append(pos.copy())

        if np.all(pos[:, 2] >= SOL_Z2):
            if not np.allclose(traj[-1], pos):
                traj.append(pos.copy())
            break

    return np.array(traj)


# ============================================================
# 2. Quadrupole simulation (Qv)
# ============================================================
G_0 = 1.5
QUAD_APER = 0.4
QUAD_Z1 = 0.3
QUAD_Z2 = 1.3
QUAD_DT = 1e-11
QUAD_TMAX = 6e-8
QUAD_STEPS = int(round(QUAD_TMAX / QUAD_DT))
QUAD_SAMPLE_RATE = 10


def get_B_field_quad(pos_vec):
    B_out = np.zeros_like(pos_vec)
    x = pos_vec[:, 0]
    y = pos_vec[:, 1]
    z = pos_vec[:, 2]
    d = 0.03
    mask = 0.5 * (np.tanh((z - QUAD_Z1) / d) - np.tanh((z - QUAD_Z2) / d))
    B_out[:, 0] = -G_0 * y * mask
    B_out[:, 1] = -G_0 * x * mask
    B_out[:, 2] = 0.0
    return B_out


def get_accel_quad(v_vec, p_vec):
    B_vec = get_B_field_quad(p_vec)
    return (q / m) * np.cross(v_vec, B_vec)


def simulate_quadrupole():
    pos, vel = initial_state()
    traj = [pos.copy()]

    print(f"[Quadrupole] integrating until z reaches {QUAD_Z2} m...")
    for i in range(QUAD_STEPS):
        k1_v = get_accel_quad(vel, pos)
        k1_p = vel
        k2_v = get_accel_quad(vel + 0.5 * k1_v * QUAD_DT, pos + 0.5 * k1_p * QUAD_DT)
        k2_p = vel + 0.5 * k1_v * QUAD_DT
        k3_v = get_accel_quad(vel + 0.5 * k2_v * QUAD_DT, pos + 0.5 * k2_p * QUAD_DT)
        k3_p = vel + 0.5 * k2_v * QUAD_DT
        k4_v = get_accel_quad(vel + k3_v * QUAD_DT, pos + k3_p * QUAD_DT)
        k4_p = vel + k3_v * QUAD_DT

        vel += (QUAD_DT / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        pos += (QUAD_DT / 6.0) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p)

        if i % QUAD_SAMPLE_RATE == 0:
            traj.append(pos.copy())

        if np.all(pos[:, 2] >= QUAD_Z2):
            if not np.allclose(traj[-1], pos):
                traj.append(pos.copy())
            break

    return np.array(traj)


# ============================================================
# 3. Run both simulations
# ============================================================
traj_sol = simulate_solenoid()
traj_quad = simulate_quadrupole()

X_sol, Y_sol, Z_sol = traj_sol[:, :, 0], traj_sol[:, :, 1], traj_sol[:, :, 2]
X_quad, Y_quad, Z_quad = traj_quad[:, :, 0], traj_quad[:, :, 1], traj_quad[:, :, 2]

frames_total = max(len(traj_sol), len(traj_quad))


# ============================================================
# 4. Build the 2x2 animated figure
# ============================================================
fig = plt.figure(figsize=(14, 12))
fig.subplots_adjust(top=0.88, left=0.06, right=0.97, bottom=0.06, wspace=0.12, hspace=0.18)

# --- Top-left: Solenoid 3D ---
ax_s3d = fig.add_subplot(221, projection='3d')
ax_s3d.set_title('Solenoid: 3D Trajectory', fontsize=13, pad=12)
ax_s3d.set_xlabel('Z axis (m)', labelpad=10)
ax_s3d.set_ylabel('X axis (m)', labelpad=10)
ax_s3d.set_zlabel('Y axis (m)', labelpad=10)
ax_s3d.set_xlim(0, MAX_Z)
ax_s3d.set_ylim(-MAX_R_XY, MAX_R_XY)
ax_s3d.set_zlim(-MAX_R_XY, MAX_R_XY)
ax_s3d.view_init(elev=22, azim=-60)

z_cyl = np.linspace(SOL_Z1, SOL_Z2, 60)
theta_cyl = np.linspace(0, 2 * np.pi, 60)
Zc, THc = np.meshgrid(z_cyl, theta_cyl)
ax_s3d.plot_surface(Zc, SOL_A * np.cos(THc), SOL_A * np.sin(THc), color='gray', alpha=0.12, linewidth=0)
ax_s3d.scatter([0], [0], [0], color='k', s=30)

# Beam-direction annotation in plotting coordinates (plot-x is physical +Z)
beam_arrow_start = (0.06, -0.24, 0.24)   # (plot_x=Z, plot_y=X, plot_z=Y)
beam_arrow_dxyz = (0.22, 0.0, 0.0)
ax_s3d.quiver(
    beam_arrow_start[0], beam_arrow_start[1], beam_arrow_start[2],
    beam_arrow_dxyz[0], beam_arrow_dxyz[1], beam_arrow_dxyz[2],
    arrow_length_ratio=0.18,
    color='black', linewidth=1.8
)
ax_s3d.text(
    beam_arrow_start[0] + beam_arrow_dxyz[0] + 0.02,
    beam_arrow_start[1],
    beam_arrow_start[2],
    '+Z = EVR',
    fontsize=10
)

sol_lines_3d = [ax_s3d.plot([], [], [], color=COLORS[i], linewidth=2.0)[0] for i in range(N_PARTICLES)]
sol_pts_3d = [ax_s3d.plot([], [], [], marker='o', color=COLORS[i], markersize=4.5)[0] for i in range(N_PARTICLES)]

# --- Top-right: Solenoid XY ---
ax_sxy = fig.add_subplot(222)
ax_sxy.set_title('Solenoid: X-Y Projection', fontsize=13, pad=12)
ax_sxy.set_xlabel('X axis (m)')
ax_sxy.set_ylabel('Y axis (m)')
ax_sxy.set_xlim(-MAX_R_XY, MAX_R_XY)
ax_sxy.set_ylim(-MAX_R_XY, MAX_R_XY)
ax_sxy.set_aspect('equal')
ax_sxy.grid(True, linestyle=':')
ax_sxy.scatter([0], [0], color='k', s=30)

# Solenoid bore / aperture range in XY projection
sol_aperture = plt.Circle((0, 0), SOL_A, fill=False, linestyle='--', linewidth=1.8, color='gray', alpha=0.95)
ax_sxy.add_artist(sol_aperture)
ax_sxy.text(0.0, SOL_A + 0.015, f'Solenoid aperture (R = {SOL_A:.2f} m)',
            ha='center', va='bottom', fontsize=10, color='gray')

sol_lines_xy = [ax_sxy.plot([], [], color=COLORS[i], linewidth=2.0)[0] for i in range(N_PARTICLES)]
sol_pts_xy = [ax_sxy.plot([], [], marker='o', color=COLORS[i], markersize=4.5)[0] for i in range(N_PARTICLES)]

# --- Bottom-left: Quadrupole 3D ---
ax_q3d = fig.add_subplot(223, projection='3d')
ax_q3d.set_title('Quadrupole (Qv): 3D Trajectory', fontsize=13, pad=12)
ax_q3d.set_xlabel('Z axis (m)', labelpad=10)
ax_q3d.set_ylabel('X axis (m)', labelpad=10)
ax_q3d.set_zlabel('Y axis (m)', labelpad=10)
ax_q3d.set_xlim(0, MAX_Z)
ax_q3d.set_ylim(-MAX_R_XY, MAX_R_XY)
ax_q3d.set_zlim(-MAX_R_XY, MAX_R_XY)
ax_q3d.view_init(elev=22, azim=-60)

pole_angles = [np.pi / 4, 3 * np.pi / 4, 5 * np.pi / 4, 7 * np.pi / 4]
z_pole = np.linspace(QUAD_Z1, QUAD_Z2, 12)
for angle in pole_angles:
    px = QUAD_APER * np.cos(angle) * np.ones_like(z_pole)
    py = QUAD_APER * np.sin(angle) * np.ones_like(z_pole)
    ax_q3d.plot(z_pole, px, py, color='gray', linewidth=7, alpha=0.45)
ax_q3d.scatter([0], [0], [0], color='k', s=30)

# Beam-direction annotation in plotting coordinates (plot-x is physical +Z)
quad_beam_arrow_start = (0.04, -0.24, 0.24)   # (plot_x=Z, plot_y=X, plot_z=Y)
quad_beam_arrow_dxyz = (0.22, 0.0, 0.0)
ax_q3d.quiver(
    quad_beam_arrow_start[0], quad_beam_arrow_start[1], quad_beam_arrow_start[2],
    quad_beam_arrow_dxyz[0], quad_beam_arrow_dxyz[1], quad_beam_arrow_dxyz[2],
    arrow_length_ratio=0.18,
    color='black', linewidth=1.8
)
ax_q3d.text(
    quad_beam_arrow_start[0] + quad_beam_arrow_dxyz[0] + 0.02,
    quad_beam_arrow_start[1],
    quad_beam_arrow_start[2],
    '+Z = EVR',
    fontsize=10
)

quad_lines_3d = [ax_q3d.plot([], [], [], color=COLORS[i], linewidth=2.0)[0] for i in range(N_PARTICLES)]
quad_pts_3d = [ax_q3d.plot([], [], [], marker='o', color=COLORS[i], markersize=4.5)[0] for i in range(N_PARTICLES)]

# --- Bottom-right: Quadrupole XY ---
ax_qxy = fig.add_subplot(224)
ax_qxy.set_title('Quadrupole (Qv): X-Y Projection', fontsize=13, pad=12)
ax_qxy.set_xlabel('X axis (m)')
ax_qxy.set_ylabel('Y axis (m)')
ax_qxy.set_xlim(-MAX_R_XY, MAX_R_XY)
ax_qxy.set_ylim(-MAX_R_XY, MAX_R_XY)
ax_qxy.set_aspect('equal')
ax_qxy.grid(True, linestyle=':')
ax_qxy.scatter([0], [0], color='k', s=30)

for idx, angle in enumerate(pole_angles):
    px, py = QUAD_APER * np.cos(angle), QUAD_APER * np.sin(angle)
    pole_color = 'red' if idx % 2 == 0 else 'blue'
    pole_text = 'N' if idx % 2 == 0 else 'S'
    circle = plt.Circle((px, py), 0.05, color=pole_color, alpha=0.25)
    ax_qxy.add_artist(circle)
    ax_qxy.text(px, py, pole_text, ha='center', va='center', fontsize=10, fontweight='bold')

quad_lines_xy = [ax_qxy.plot([], [], color=COLORS[i], linewidth=2.0)[0] for i in range(N_PARTICLES)]
quad_pts_xy = [ax_qxy.plot([], [], marker='o', color=COLORS[i], markersize=4.5)[0] for i in range(N_PARTICLES)]

# Shared legend above the whole 2x2 figure
legend_handles = [
    Line2D([0], [0], color=COLORS[i], lw=2.5, marker='o', markersize=5, label=LABELS[i])
    for i in range(N_PARTICLES)
]
fig.legend(
    handles=legend_handles,
    labels=LABELS,
    loc='upper center',
    bbox_to_anchor=(0.5, 0.985),
    ncol=4,
    frameon=False,
    fontsize=14,
)


def _set_3d_line_and_point(line, point, Z, X, Y, idx, n):
    line.set_data(Z[:idx, n], X[:idx, n])
    line.set_3d_properties(Y[:idx, n])
    point.set_data(Z[idx:idx+1, n], X[idx:idx+1, n])
    point.set_3d_properties(Y[idx:idx+1, n])


def _set_2d_line_and_point(line, point, X, Y, idx, n):
    line.set_data(X[:idx, n], Y[:idx, n])
    point.set_data(X[idx:idx+1, n], Y[idx:idx+1, n])


def update(frame):
    if frame % 50 == 0:
        print(f"Rendering frame {frame}/{frames_total}")

    idx_sol = min(max(1, frame), len(traj_sol) - 1)
    idx_quad = min(max(1, frame), len(traj_quad) - 1)

    artists = []
    for n in range(N_PARTICLES):
        _set_3d_line_and_point(sol_lines_3d[n], sol_pts_3d[n], Z_sol, X_sol, Y_sol, idx_sol, n)
        _set_2d_line_and_point(sol_lines_xy[n], sol_pts_xy[n], X_sol, Y_sol, idx_sol, n)
        _set_3d_line_and_point(quad_lines_3d[n], quad_pts_3d[n], Z_quad, X_quad, Y_quad, idx_quad, n)
        _set_2d_line_and_point(quad_lines_xy[n], quad_pts_xy[n], X_quad, Y_quad, idx_quad, n)

        artists.extend([
            sol_lines_3d[n], sol_pts_3d[n],
            sol_lines_xy[n], sol_pts_xy[n],
            quad_lines_3d[n], quad_pts_3d[n],
            quad_lines_xy[n], quad_pts_xy[n],
        ])
    return artists


ani = animation.FuncAnimation(fig, update, frames=frames_total, blit=False, interval=30)

print('Saving GIF: combined_solenoid_quadrupole.gif')
ani.save('combined_solenoid_quadrupole.gif', writer='pillow', fps=30)
print('Done! Output: combined_solenoid_quadrupole.gif')
