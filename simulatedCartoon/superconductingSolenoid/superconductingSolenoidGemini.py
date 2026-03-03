import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 1. 物理常数与模拟参数设置
# ==========================================
q = 1.602e-19       
m = 1.672e-27       
B_z = 2.0           

E_k_MeV = 10.0      
E_k_J = E_k_MeV * 1e6 * 1.602e-19
v0 = np.sqrt(2 * E_k_J / m)

N_particles = 4
theta_emit = np.array([1, 6, 16, 21])     
phi_emit = np.array([0, 90, 180, 270])    

K_drag_visual = 2e29 

dt = 1e-10          
t_max = 5e-8        
steps = int(round(t_max / dt))

# ==========================================
# 2. 初始化状态与 RK4 积分引擎
# ==========================================
pos = np.zeros((N_particles, 3))
vel = np.zeros((N_particles, 3))

theta = theta_emit * np.pi / 180.0
phi = phi_emit * np.pi / 180.0

vel[:, 0] = v0 * np.sin(theta) * np.cos(phi)
vel[:, 1] = v0 * np.sin(theta) * np.sin(phi)
vel[:, 2] = v0 * np.cos(theta)

# 【平滑优化】：降低采样率，增加总帧数，让动画变慢变细腻
sample_rate = 2 
num_samples = (steps - 1) // sample_rate + 1
trajectory = np.zeros((num_samples, N_particles, 3))
B_vec = np.array([0, 0, B_z])

def get_accel(v_vec):
    v_cross_B = np.cross(v_vec, B_vec)
    a_L = (q / m) * v_cross_B
    v_mag = np.linalg.norm(v_vec, axis=1, keepdims=True)
    v_mag[v_mag == 0] = 1e-10 
    a_drag = - (K_drag_visual / (v_mag**2)) * (v_vec / v_mag)
    return a_L + a_drag

for i in range(steps):
    if i % sample_rate == 0:
        trajectory[i // sample_rate, :, :] = pos
    
    k1_v = get_accel(vel)
    k1_p = vel
    
    k2_v = get_accel(vel + 0.5 * k1_v * dt)
    k2_p = vel + 0.5 * k1_v * dt
    
    k3_v = get_accel(vel + 0.5 * k2_v * dt)
    k3_p = vel + 0.5 * k2_v * dt
    
    k4_v = get_accel(vel + k3_v * dt)
    k4_p = vel + k3_v * dt
    
    vel += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    pos += (dt / 6.0) * (k1_p + 2*k2_p + 2*k3_p + k4_p)

# 提取各轴数据，并计算极径 R
X = trajectory[:, :, 0]
Y = trajectory[:, :, 1]
Z = trajectory[:, :, 2]
R = np.sqrt(X**2 + Y**2)  # 新增：计算每个粒子的极径

# ==========================================
# 3. 设置 1x3 大画布与动画
# ==========================================
# 加宽画布以容纳三个图
fig = plt.figure(figsize=(22, 6.5))
fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95)

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'] 
max_z = np.max(Z)
max_r_xy = np.max(R) * 1.1 

# --- 左图：3D 轨迹 ---
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('3D Trajectory in Solenoid', fontsize=14, pad=15)
ax1.set_xlabel('Z axis (m)', fontsize=11)
ax1.set_ylabel('X axis (m)', fontsize=11)
ax1.set_zlabel('Y axis (m)', fontsize=11)
ax1.set_xlim(0, max_z)
ax1.set_ylim(-max_r_xy, max_r_xy)
ax1.set_zlim(-max_r_xy, max_r_xy)
ax1.plot([0, max_z], [0, 0], [0, 0], 'k--', linewidth=1.5, alpha=0.5)
lines3d = [ax1.plot([], [], [], color=c, linewidth=2.0)[0] for c in colors]
points3d = [ax1.plot([], [], [], marker='o', color=c, markersize=6)[0] for c in colors]

# --- 中图：X-Y 投影 ---
ax2 = fig.add_subplot(132)
ax2.set_title('X-Y Projection', fontsize=14, pad=15)
ax2.set_xlabel('X axis (m)', fontsize=11)
ax2.set_ylabel('Y axis (m)', fontsize=11)
ax2.set_xlim(-max_r_xy, max_r_xy)
ax2.set_ylim(-max_r_xy, max_r_xy)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.set_aspect('equal')
lines2d = [ax2.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points2d = [ax2.plot([], [], marker='o', color=c, markersize=6)[0] for c in colors]

# --- 右图：极径 R 随 Z 轴的变化 ---
ax3 = fig.add_subplot(133)
ax3.set_title('Polar Radius ($\\rho$) vs Z-axis', fontsize=14, pad=15)
ax3.set_xlabel('Z axis (m)', fontsize=11)
ax3.set_ylabel('Polar Radius $\\rho$ (m)', fontsize=11)
ax3.set_xlim(0, max_z)
ax3.set_ylim(0, max_r_xy)
ax3.grid(True, linestyle='--', alpha=0.7)
lines_r = [ax3.plot([], [], color=c, linewidth=2.0)[0] for c in colors]
points_r = [ax3.plot([], [], marker='o', color=c, markersize=6)[0] for c in colors]

frames_total = num_samples 

def init():
    for i in range(N_particles):
        lines3d[i].set_data([], [])
        lines3d[i].set_3d_properties([])
        points3d[i].set_data([], [])
        points3d[i].set_3d_properties([])
        
        lines2d[i].set_data([], [])
        points2d[i].set_data([], [])
        
        lines_r[i].set_data([], [])
        points_r[i].set_data([], [])
        
    return lines3d + points3d + lines2d + points2d + lines_r + points_r

def update(frame):
    if frame % 50 == 0:
        print(f"后台渲染中 (三图联动)... 已完成 {frame}/{frames_total} 帧")
        
    idx = frame 
    if idx == 0: idx = 1 
        
    for n in range(N_particles):
        lines3d[n].set_data(Z[:idx, n], X[:idx, n])
        lines3d[n].set_3d_properties(Y[:idx, n])
        points3d[n].set_data(Z[idx:idx+1, n], X[idx:idx+1, n])
        points3d[n].set_3d_properties(Y[idx:idx+1, n])
        
        lines2d[n].set_data(X[:idx, n], Y[:idx, n])
        points2d[n].set_data(X[idx:idx+1, n], Y[idx:idx+1, n])
        
        # 更新第三个图 (极径随 Z 轴变化)
        lines_r[n].set_data(Z[:idx, n], R[:idx, n])
        points_r[n].set_data(Z[idx:idx+1, n], R[idx:idx+1, n])
    
    return lines3d + points3d + lines2d + points2d + lines_r + points_r

ani = animation.FuncAnimation(fig, update, frames=frames_total, 
                              init_func=init, blit=False, interval=20)

# ==========================================
# 4. 生成并保存
# ==========================================
print(f"开始生成 {frames_total} 帧的三图联动宽屏动画...")
# 降低 fps 以实现平滑慢放
ani.save('solenoid_trajectory_3panels.gif', writer='pillow', fps=20)
print("GIF 生成完毕！快去查看 solenoid_trajectory_3panels.gif！")