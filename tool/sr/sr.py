import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad

def calculate_and_visualize_solid_angle(ang_px, ang_nx, ang_py, ang_ny):
    """
    计算并绘制由四个方向角度定义的立体角。
    输入单位为度 (degrees)。
    """
    # 1. 数据准备
    # 将角度转换为弧度
    t_px = np.deg2rad(ang_px)
    t_nx = np.deg2rad(ang_nx)
    t_py = np.deg2rad(ang_py)
    t_ny = np.deg2rad(ang_ny)

    # 计算平均半轴角度 (假设可能为椭圆锥)
    # 如果四个角度相等，theta_x = theta_y，即为正圆
    theta_x = (t_px + t_nx) / 2.0
    theta_y = (t_py + t_ny) / 2.0

    print("-" * 30)
    print(f"输入角度: +x={ang_px}°, -x={ang_nx}°, +y={ang_py}°, -y={ang_ny}°")

    # 2. 立体角计算 (数值积分法，适用于圆或椭圆)
    # 定义椭圆锥边界函数 theta(phi)
    # 基于公式: tan(theta)^2 = 1 / ( (cos(phi)/tan(theta_x))^2 + (sin(phi)/tan(theta_y))^2 )
    def boundary_theta(phi):
        term_x = (np.cos(phi) / np.tan(theta_x)) ** 2
        term_y = (np.sin(phi) / np.tan(theta_y)) ** 2
        tan_theta = 1.0 / np.sqrt(term_x + term_y)
        return np.arctan(tan_theta)

    # 立体角公式: Omega = integral(0 to 2pi) of (1 - cos(theta(phi))) d_phi
    def integrand(phi):
        return 1.0 - np.cos(boundary_theta(phi))

    # 执行积分
    result, error = quad(integrand, 0, 2 * np.pi)
    
    print(f"计算立体角 (Omega): {result:.5f} sr")
    if abs(theta_x - theta_y) < 1e-5:
        # 验证正圆公式
        check = 2 * np.pi * (1 - np.cos(theta_x))
        print(f"正圆公式验证:     {check:.5f} sr")
    print("-" * 30)

    # 3. 3D 可视化
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # --- A. 绘制单位球 (Unit Sphere) ---
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50) # 完整球体
    x_sph = np.outer(np.cos(u), np.sin(v))
    y_sph = np.outer(np.sin(u), np.sin(v))
    z_sph = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(x_sph, y_sph, z_sph, color='gray', alpha=0.1, linewidth=0.5)

    # --- B. 计算锥体数据 ---
    phi_vals = np.linspace(0, 2*np.pi, 100)
    theta_vals = np.array([boundary_theta(p) for p in phi_vals])
    
    # 物理平面距离 (为了视觉效果放在 Z=2)
    dist_phys = 2.0 
    
    # --- C. 绘制投影区域 (在单位球面上, R=1) ---
    # 构建投影面的网格
    # 从中心轴(theta=0) 到 边界(theta=theta_vals) 插值
    r_grid = np.linspace(0, 1, 10) # 径向插值因子
    X_proj_surf = []
    Y_proj_surf = []
    Z_proj_surf = []
    
    for r_factor in r_grid:
        # 当前环的 theta = r_factor * 边界theta (线性近似用于填充颜色)
        curr_thetas = theta_vals * r_factor
        X_proj_surf.append(np.sin(curr_thetas) * np.cos(phi_vals))
        Y_proj_surf.append(np.sin(curr_thetas) * np.sin(phi_vals))
        Z_proj_surf.append(np.cos(curr_thetas))
    
    ax.plot_surface(np.array(X_proj_surf), np.array(Y_proj_surf), np.array(Z_proj_surf), 
                    color='cyan', alpha=0.6, label='Solid Angle Projection')
    
    # 绘制投影边缘 (蓝色)
    x_proj_rim = 1.0 * np.sin(theta_vals) * np.cos(phi_vals)
    y_proj_rim = 1.0 * np.sin(theta_vals) * np.sin(phi_vals)
    z_proj_rim = 1.0 * np.cos(theta_vals)
    ax.plot(x_proj_rim, y_proj_rim, z_proj_rim, color='blue', lw=2, label='Projected Area (r=1)')

    # --- D. 绘制物理物体 (在 Z=dist_phys 处) ---
    # 物理半径 = dist_phys * tan(theta)
    r_phys = dist_phys * np.tan(theta_vals)
    x_phys = r_phys * np.cos(phi_vals)
    y_phys = r_phys * np.sin(phi_vals)
    z_phys = np.full_like(phi_vals, dist_phys)
    
    ax.plot(x_phys, y_phys, z_phys, color='red', lw=3, label=f'Physical Ring (z={dist_phys})')

    # --- E. 绘制辅助线 (从原点到物理环) ---
    # 每隔 45度画一条线
    for i in range(0, len(phi_vals), 12):
        ax.plot([0, x_phys[i]], [0, y_phys[i]], [0, z_phys[i]], 'r--', alpha=0.3, lw=1)

    # 绘制特定四个方向的轴线 (+x, -x, +y, -y)
    # 为了准确，我们需要找到 phi 接近 0, 90, 180, 270 的索引
    axis_phis = [0, np.pi/2, np.pi, 3*np.pi/2]
    axis_labels = ["+x", "+y", "-x", "-y"]
    
    for i, a_phi in enumerate(axis_phis):
        # 计算该方向的边界 theta
        b_theta = boundary_theta(a_phi)
        px = dist_phys * np.tan(b_theta) * np.cos(a_phi)
        py = dist_phys * np.tan(b_theta) * np.sin(a_phi)
        pz = dist_phys
        ax.plot([0, px], [0, py], [0, pz], 'k-', lw=1.5)
        ax.text(px, py, pz, f'{axis_labels[i]} ({np.degrees(b_theta):.1f}°)', fontsize=10)

    # --- F. 设置图形参数 ---
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([0, 2.5])
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    
    # 标注原点
    ax.scatter([0], [0], [0], color='black', s=50, label='Origin (Observer)')
    
    title_str = f"Solid Angle Visualization\nOmega = {result:.4f} sr"
    ax.set_title(title_str, fontsize=14)
    ax.view_init(elev=25, azim=45)
    ax.legend(loc='lower right')
    
    plt.tight_layout()
    plt.show()

# ==========================================
# 在这里输入您的四个方向角度
# ==========================================
# 题目设定: 全部为 6.3 度
angle_plus_x  = 6.5
angle_minus_x = 6.5
angle_plus_y  = 4
angle_minus_y = 4

# 运行函数
calculate_and_visualize_solid_angle(angle_plus_x, angle_minus_x, angle_plus_y, angle_minus_y)