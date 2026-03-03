import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 1. 物理参数、真实磁场与气体阻尼设定
# ==========================================
q = 1.602e-19       
m = 1.672e-27       
E_k_MeV = 10.0      
E_k_J = E_k_MeV * 1e6 * 1.602e-19
v0 = np.sqrt(2 * E_k_J / m)

# 真实螺线管参数 (处于靶后一段距离)
B_0 = 2.5           # 中心磁场 (T)
a = 0.20            # 螺线管半径 (m)
z1 = 0.3            # 入口位置 (m)
z2 = 1.3            # 出口位置 (m)

# 100Pa He气能损模拟阻尼系数
K_drag_visual = 3e29 

# 4个点源发散粒子的初始发射角
N_particles = 4
theta_emit = np.array([1, 6, 11, 16])     # 极角 (度)
phi_emit = np.array([0, 90, 180, 270])    # 方位角 (度)

dt = 1e-11          
t_max = 5e-8        
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

sample_rate = 20 
num_samples = (steps - 1) // sample_rate + 1
trajectory = np.zeros((num_samples, N_particles, 3))

def get_B_field_real(pos_vec):
    B_out = np.zeros_like(pos_vec)
    x = pos_vec[:, 0]
    y = pos_vec[:, 1]
    z = pos_vec[:, 2]
    r = np.sqrt(x**2 + y**2)
    r_safe = np.where(r == 0, 1e-10, r)
    
    term1 = (z - z1) / np.sqrt((z - z1)**2 + a**2)
    term2 = (z - z2) / np.sqrt((z - z2)**2 + a**2)
    Bz = (B_0 / 2.0) * (term1 - term2)
    
    dterm1 = (a**2) / ((z - z1)**2 + a**2)**1.5
    dterm2 = (a**2) / ((z - z2)**2 + a**2)**1.5
    dBz_dz = (B_0 / 2.0) * (dterm1 - dterm2)
    Br = - (r / 2.0) * dBz_dz
    
    B_out[:, 0] = Br * (x / r_safe)
    B_out[:, 1] = Br * (y / r_safe)
    B_out[:, 2] = Bz
    return B_out

def get_accel_combined(v_vec, p_vec):
    B_vec = get_B_field_real(p_vec)
    a_L = (q / m) * np.cross(v_vec, B_vec)
    
    v_mag = np.linalg.norm(v_vec, axis=1, keepdims=True)
    v_mag[v_mag == 0] = 1e-10 
    a_drag = - (K_drag_visual / (v_mag**2)) * (v_vec / v_mag)
    
    return a_L + a_drag

print(f"开始终极物理引擎计算 (真实磁场 + He气冷却), 共 {steps} 步...")
for i in range(steps):
    if i % sample_rate == 0:
        trajectory[i // sample_rate, :, :] = pos
    
    k1_v = get_accel_combined(vel, pos)
    k1_p = vel
    k2_v = get_accel_combined(vel + 0.5 * k1_v * dt, pos + 0.5 * k1_p * dt)
    k2_p = vel + 0.5 * k1_v * dt
    k3_v = get_accel_combined(vel + 0.5 * k2_v * dt, pos + 0.5 * k2_p * dt)
    k3_p = vel + 0.5 * k2_v * dt
    k4_v = get_accel_combined(vel + k3_v * dt, pos + k3_p * dt)
    k4_p = vel + k3_v * dt
    
    vel += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    pos += (dt / 6.0) * (k1_p + 2*k2_p + 2*k3_p + k4_p)

X = trajectory[:, :, 0]
Y = trajectory[:, :, 1]
Z = trajectory[:, :, 2]
R = np.sqrt(X**2 + Y**2)

# ==========================================
# 4. 宽屏三图联动可视化设置 (图例统一版)
# ==========================================
fig = plt.figure(figsize=(22, 6.5))
fig.subplots_adjust(wspace=0.35, left=0.05, right=0.95, bottom=0.15)

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'] 
# 【新增】定义图例标签，对应初始发射角 1°, 6°, 11°, 16°
labels = ['P1: +X (1°)', 'P2: +Y (6°)', 'P3: -X (11°)', 'P4: -Y (16°)']

max_z = 2.5
max_r_xy = 0.25 

# --- 左图：3D 轨迹图 ---
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('3D Trajectory in Superconducting Solenoid', fontsize=14, pad=20)
ax1.set_xlabel('Z axis (m)', labelpad=15)
ax1.set_ylabel('X axis (m)', labelpad=15)
ax1.set_zlabel('Y axis (m)', labelpad=15)
ax1.set_xlim(0, max_z)
ax1.set_ylim(-max_r_xy, max_r_xy)
ax1.set_zlim(-max_r_xy, max_r_xy)

z_cyl = np.linspace(z1, z2, 50)
theta_cyl = np.linspace(0, 2*np.pi, 50)
Z_c, THETA_c = np.meshgrid(z_cyl, theta_cyl)
ax1.plot_surface(Z_c, a*np.cos(THETA_c), a*np.sin(THETA_c), color='gray', alpha=0.15)
ax1.scatter([0],[0],[0], color='k', s=50, label='Target Point')

# 将原本的 for c in colors 改为 for i in range(N_particles)，方便统一索引
lines3d = [ax1.plot([], [], [], color=colors[i], linewidth=2.0)[0] for i in range(N_particles)]
points3d = [ax1.plot([], [], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]

# --- 中图：X-Y 投影图 ---
ax2 = fig.add_subplot(132)
ax2.set_title('X-Y Projection ', fontsize=14, pad=20)
ax2.set_xlabel('X axis (m)')
ax2.set_ylabel('Y axis (m)')
ax2.set_xlim(-max_r_xy, max_r_xy)
ax2.set_ylim(-max_r_xy, max_r_xy)
ax2.grid(True, linestyle=':')
ax2.set_aspect('equal')
ax2.scatter([0],[0], color='k', s=50)

lines2d = [ax2.plot([], [], color=colors[i], linewidth=2.0)[0] for i in range(N_particles)]
points2d = [ax2.plot([], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]

# --- 右图：极径 R-Z 演化图 ---
ax3 = fig.add_subplot(133)
ax3.set_title(r'Polar Radius ($\rho$ vs z)', fontsize=14, pad=20)
ax3.set_xlabel('Z axis (m)')
ax3.set_ylabel(r'Polar Radius $\rho$ (m)')
ax3.set_xlim(0, max_z)
ax3.set_ylim(0, max_r_xy)
ax3.grid(True, linestyle=':')

# 背景区域
ax3.axvspan(z1, z2, color='gray', alpha=0.15, label='Solenoid B-field Region')

# 【修改】给每一条线绑定 label
lines_r = [ax3.plot([], [], color=colors[i], linewidth=2.0, label=labels[i])[0] for i in range(N_particles)]
points_r = [ax3.plot([], [], marker='o', color=colors[i], markersize=5)[0] for i in range(N_particles)]

# 【修改】最后再调用 legend 生成图例，并且设置字体大小与四极铁保持一致
ax3.legend(loc='upper right', fontsize=9)

frames_total = num_samples 

def update(frame):
    if frame % 50 == 0: print(f"渲染中... {frame}/{frames_total}")
    idx = max(1, frame)
    for n in range(N_particles):
        lines3d[n].set_data(Z[:idx, n], X[:idx, n])
        lines3d[n].set_3d_properties(Y[:idx, n])
        points3d[n].set_data(Z[idx:idx+1, n], X[idx:idx+1, n])
        points3d[n].set_3d_properties(Y[idx:idx+1, n])
        
        lines2d[n].set_data(X[:idx, n], Y[:idx, n])
        points2d[n].set_data(X[idx:idx+1, n], Y[idx:idx+1, n])
        
        lines_r[n].set_data(Z[:idx, n], R[:idx, n])
        points_r[n].set_data(Z[idx:idx+1, n], R[idx:idx+1, n])
    
    return lines3d + points3d + lines2d + points2d + lines_r + points_r

ani = animation.FuncAnimation(fig, update, frames=frames_total, blit=False, interval=30)

# ==========================================
# 5. 生成保存
# ==========================================
print(f"正在生成宽屏动画 (修复版)，请稍候...")
ani.save('ultimate_solenoid_simulation_fixed.gif', writer='pillow', fps=25)
print("GIF 渲染完成！请查看 ultimate_solenoid_simulation_fixed.gif")