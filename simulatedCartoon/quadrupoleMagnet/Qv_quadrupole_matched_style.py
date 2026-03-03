import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 1. 物理参数与四极铁磁场设定
# ==========================================
q = 1.602e-19       # 质子电荷 (C)
m = 1.672e-27       # 质子质量 (kg)
E_k_MeV = 10.0      # 动能 10 MeV
E_k_J = E_k_MeV * 1e6 * 1.602e-19
v0 = np.sqrt(2 * E_k_J / m)

# 理想四极铁参数 (Qv / QD: 垂直聚焦，水平散焦)
# 依据文献中的磁场展开：Bx = -G*y, By = -G*x (当粒子带正电且沿+Z飞行时)
# 【调整】将 G_0 降至 1.5 T/m，以避免过聚焦（多次跨越轴线），呈现平滑的聚焦包络
G_0 = 1.5           # 四极铁磁场梯度 (T/m) 
a_aper = 0.4        # 四极铁孔径 (m) (仅用于绘图示意)
z1 = 0.3            # 磁铁入口位置 (m)
z2 = 1.3            # 磁铁出口位置 (m)

# 4个点源发散粒子的初始发射角 (保留用户的设定)
N_particles = 4
theta_emit = np.array([1, 6, 11, 16])      # 极角 (度)
phi_emit = np.array([0, 90, 180, 270])    # 方位角 (度) 对应 +X, +Y, -X, -Y

dt = 1e-11          
t_max = 6e-8        
steps = int(round(t_max / dt))

# ==========================================
# 2. 初始化状态与 RK4 积分
# ==========================================
pos = np.zeros((N_particles, 3))
vel = np.zeros((N_particles, 3))

theta = theta_emit * np.pi / 180.0
phi = phi_emit * np.pi / 180.0
vel[:, 0] = v0 * np.sin(theta) * np.cos(phi)
vel[:, 1] = v0 * np.sin(theta) * np.sin(phi)
vel[:, 2] = v0 * np.cos(theta)

sample_rate = 25 
num_samples = (steps - 1) // sample_rate + 1
trajectory = np.zeros((num_samples, N_particles, 3))

def get_B_field_quad(pos_vec):
    """
    基于文献 Hill's Equation 理论的四极铁磁场
    利用 tanh 函数作为 Fringe Field (边缘场) 平滑过渡，确保积分稳定性
    """
    B_out = np.zeros_like(pos_vec)
    x = pos_vec[:, 0]
    y = pos_vec[:, 1]
    z = pos_vec[:, 2]
    
    # 边缘场平滑度参数
    d = 0.03
    # 纵向分布函数 (平滑化的有效场长)
    mask = 0.5 * (np.tanh((z - z1)/d) - np.tanh((z - z2)/d))
    
    # Qv (垂直方向聚焦) 磁场分量
    B_out[:, 0] = -G_0 * y * mask
    B_out[:, 1] = -G_0 * x * mask
    B_out[:, 2] = 0.0
    return B_out

def get_accel_pure(v_vec, p_vec):
    # 纯哈密顿束流光学保守力场
    B_vec = get_B_field_quad(p_vec)
    a_L = (q / m) * np.cross(v_vec, B_vec)
    return a_L

print(f"开始四极铁 (Qv) 束流光学物理积分计算，共 {steps} 步...")
for i in range(steps):
    if i % sample_rate == 0:
        trajectory[i // sample_rate, :, :] = pos
    
    k1_v = get_accel_pure(vel, pos)
    k1_p = vel
    k2_v = get_accel_pure(vel + 0.5 * k1_v * dt, pos + 0.5 * k1_p * dt)
    k2_p = vel + 0.5 * k1_v * dt
    k3_v = get_accel_pure(vel + 0.5 * k2_v * dt, pos + 0.5 * k2_p * dt)
    k3_p = vel + 0.5 * k2_v * dt
    k4_v = get_accel_pure(vel + k3_v * dt, pos + k3_p * dt)
    k4_p = vel + k3_v * dt
    
    vel += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    pos += (dt / 6.0) * (k1_p + 2*k2_p + 2*k3_p + k4_p)

X = trajectory[:, :, 0]
Y = trajectory[:, :, 1]
Z = trajectory[:, :, 2]

# ==========================================
# 3. 宽屏三图联动可视化设置 (四极铁定制版)
# ==========================================
fig = plt.figure(figsize=(22, 6.5))
fig.subplots_adjust(wspace=0.35, left=0.05, right=0.95, bottom=0.15)

colors = ['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e'] 
labels = [f'P{i+1}: +X (Defocus, {theta_emit[i]}°)' if i % 2 == 0 else f'P{i+1}: +Y (Focus, {theta_emit[i]}°)' for i in range(N_particles)]

max_z = 2.5
# 【调整】因为 22度发散角非常大，需要更大的横向视场
max_r_xy = 1.0 

# --- 左图：3D 轨迹图 ---
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('3D Trajectory in Quadrupole (Qv)', fontsize=14, pad=20)
ax1.set_xlabel('Z axis (m)', labelpad=15)
ax1.set_ylabel('X axis (m)', labelpad=15)
ax1.set_zlabel('Y axis (m)', labelpad=15)
ax1.set_xlim(0, max_z)
ax1.set_ylim(-max_r_xy, max_r_xy)
ax1.set_zlim(-max_r_xy, max_r_xy)

# 在3D空间画出四极铁的四个极头
pole_angles = [np.pi/4, 3*np.pi/4, 5*np.pi/4, 7*np.pi/4]
z_pole = np.linspace(z1, z2, 10)
for angle in pole_angles:
    px = a_aper * np.cos(angle) * np.ones_like(z_pole)
    py = a_aper * np.sin(angle) * np.ones_like(z_pole)
    ax1.plot(z_pole, px, py, color='gray', linewidth=8, alpha=0.5)

ax1.scatter([0],[0],[0], color='k', s=50, label='Source')

lines3d = [ax1.plot([], [], [], color=colors[i], linewidth=2.0)[0] for i in range(N_particles)]
points3d = [ax1.plot([], [], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]

# --- 中图：X-Y 投影图 ---
ax2 = fig.add_subplot(132)
ax2.set_title('X-Y Projection', fontsize=14, pad=20)
ax2.set_xlabel('X axis (m)')
ax2.set_ylabel('Y axis (m)')
ax2.set_xlim(-max_r_xy, max_r_xy)
ax2.set_ylim(-max_r_xy, max_r_xy)
ax2.grid(True, linestyle=':')
ax2.set_aspect('equal')

# 画截面极头 (N极和S极)
for idx, angle in enumerate(pole_angles):
    px, py = a_aper * np.cos(angle), a_aper * np.sin(angle)
    color = 'red' if idx % 2 == 0 else 'blue'
    text = 'N' if idx % 2 == 0 else 'S'
    circle = plt.Circle((px, py), 0.05, color=color, alpha=0.3)
    ax2.add_artist(circle)
    ax2.text(px, py, text, color='k', ha='center', va='center', fontweight='bold')

ax2.scatter([0],[0], color='k', s=50)
lines2d = [ax2.plot([], [], color=colors[i], linewidth=2.0)[0] for i in range(N_particles)]
points2d = [ax2.plot([], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]

# --- 右图：横向位置演化图 (展示 Betatron Oscillation) ---
ax3 = fig.add_subplot(133)
ax3.set_title('Transverse Coordinate vs Z', fontsize=14, pad=20)
ax3.set_xlabel('Z axis (m)')
ax3.set_ylabel('Active Transverse Coordinate (m)')
ax3.set_xlim(0, max_z)
ax3.set_ylim(-max_r_xy, max_r_xy)
ax3.grid(True, linestyle=':')

ax3.axvspan(z1, z2, color='gray', alpha=0.15, label='Quadrupole Field Region')

lines_1d = [ax3.plot([], [], color=colors[i], linewidth=2.0, label=labels[i])[0] for i in range(N_particles)]
points_1d = [ax3.plot([], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]
ax3.legend(loc='upper right', fontsize=9)

frames_total = num_samples 

def update(frame):
    if frame % 50 == 0: print(f"渲染中... {frame}/{frames_total}")
    idx = max(1, frame)
    for n in range(N_particles):
        # 3D update
        lines3d[n].set_data(Z[:idx, n], X[:idx, n])
        lines3d[n].set_3d_properties(Y[:idx, n])
        points3d[n].set_data(Z[idx:idx+1, n], X[idx:idx+1, n])
        points3d[n].set_3d_properties(Y[idx:idx+1, n])
        
        # 2D XY update
        lines2d[n].set_data(X[:idx, n], Y[:idx, n])
        points2d[n].set_data(X[idx:idx+1, n], Y[idx:idx+1, n])
        
        # 1D Transverse vs Z update
        # 粒子 0和2 在 X 轴上运动，粒子 1和3 在 Y 轴上运动
        if n == 0 or n == 2:
            transverse_coord = X
        else:
            transverse_coord = Y
            
        lines_1d[n].set_data(Z[:idx, n], transverse_coord[:idx, n])
        points_1d[n].set_data(Z[idx:idx+1, n], transverse_coord[idx:idx+1, n])
    
    return lines3d + points3d + lines2d + points2d + lines_1d + points_1d

ani = animation.FuncAnimation(fig, update, frames=frames_total, blit=False, interval=30)

# ==========================================
# 4. 生成保存
# ==========================================
print(f"正在生成宽屏动画 (四极铁聚焦 Qv)，请稍候...")
ani.save('Qv_quadrupole.gif', writer='pillow', fps=25)
print("GIF 渲染完成！请查看 Qv_quadrupole.gif")