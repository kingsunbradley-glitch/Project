import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 1. 计算核心 (保留您提供的优秀计算逻辑)
# ==========================================
def solid_angle_elliptic_cone(theta_xp_deg, theta_xm_deg,
                              theta_yp_deg, theta_ym_deg,
                              n_phi=2000):
    """
    用四个方向(+x, -x, +y, -y)的开角(单位:度)定义一个椭圆锥，
    数值计算其对应的立体角。
    返回: 立体角值, 方位角数组, 对应的最大天顶角数组
    对称分布的话正负角度相等。
    """
    # 角度 -> 弧度
    txp, txm, typ, tym = np.radians(
        [theta_xp_deg, theta_xm_deg, theta_yp_deg, theta_ym_deg]
    )

    # 在切平面上的半轴长度: 取 tan(theta) 并在 ± 方向上做平均
    # 这定义了物理空间中 z=1 平面上的椭圆形状参数
    ax = 0.5 * (np.tan(txp) + np.tan(txm))  # x 方向半轴
    ay = 0.5 * (np.tan(typ) + np.tan(tym))  # y 方向半轴

    # 方位角采样
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)

    # 椭圆在切平面上的极径 r_tan(phi)
    # 公式: r = (ax * ay) / sqrt((ay * cosφ)^2 + (ax * sinφ)^2)
    r_tan = ax * ay / np.sqrt((ay * np.cos(phi))**2 +
                              (ax * np.sin(phi))**2)

    # r_tan = tan(theta_max)  =>  theta_max(phi)
    # 这是单位球面上投影边界的天顶角
    theta_max = np.arctan(r_tan)

    # 立体角积分: Ω = ∫0^{2π} [1 - cos θ_max(φ)] dφ
    integrand = 1.0 - np.cos(theta_max)
    # 使用梯形法则进行数值积分
    Omega = np.trapz(integrand, phi)

    return Omega, phi, theta_max


# ==========================================
# 2. 新的绘图函数 (采用改进的视觉风格)
# ==========================================
def plot_improved_geometry(phi, theta_max, omega_val, z_ring=2.0):
    """
    使用改进的视觉风格绘制几何体：
    - 单位球 (灰色线框)
    - 立体角投影面 (青色实体表面，位于单位球上)
    - 物理空间的椭圆环 (红色实线，位于 z=z_ring)
    - 连接原点和物理环的视线 (红色虚线)
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # ---- A. 画单位球面 (作为背景参考) ----
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_sph = np.outer(np.cos(u), np.sin(v))
    y_sph = np.outer(np.sin(u), np.sin(v))
    z_sph = np.outer(np.ones_like(u), np.cos(v))
    # 使用淡灰色线框，避免遮挡主体
    ax.plot_wireframe(x_sph, y_sph, z_sph, color='gray', alpha=0.15, linewidth=0.5)

    # ---- B. 画立体角在单位球面上的投影区域 (核心展示) ----
    # 创建用于填充颜色的网格数据
    # 从球心(theta=0)向边界(theta=theta_max)进行插值
    r_interp = np.linspace(0, 1, 30) # 径向插值点
    theta_grid = np.outer(r_interp, theta_max) # (30, n_phi)
    phi_grid = np.outer(np.ones_like(r_interp), phi) # (30, n_phi)

    # 转换为笛卡尔坐标用于绘图 (半径 R=1 的单位球上)
    x_proj_surf = np.sin(theta_grid) * np.cos(phi_grid)
    y_proj_surf = np.sin(theta_grid) * np.sin(phi_grid)
    z_proj_surf = np.cos(theta_grid)

    # 绘制青色实体表面，这个表面的面积数值就是立体角 Omega
    surf = ax.plot_surface(x_proj_surf, y_proj_surf, z_proj_surf,
                           color='cyan', alpha=0.8, linewidth=0,
                           label='Solid Angle Projection (Area = Ω)')

    # 绘制投影面的边界线 (深蓝色)
    x_proj_rim = np.sin(theta_max) * np.cos(phi)
    y_proj_rim = np.sin(theta_max) * np.sin(phi)
    z_proj_rim = np.cos(theta_max)
    ax.plot(x_proj_rim, y_proj_rim, z_proj_rim, color='darkblue', lw=2, label='Projection Boundary (on Unit Sphere)')

    # ---- C. 画物理空间中的实际物体 (椭圆环) ----
    # 物理环的形状必须严格对应计算出的 theta_max
    # 物理半径 R_phys(phi) = z_ring * tan(theta_max(phi))
    r_phys = z_ring * np.tan(theta_max)
    x_phys = r_phys * np.cos(phi)
    y_phys = r_phys * np.sin(phi)
    z_phys = np.full_like(phi, z_ring)

    # 绘制红色物理环
    ax.plot(x_phys, y_phys, z_phys, color='red', lw=3, label=f'Physical Ring Object (at z={z_ring})')

    # ---- D. 画辅助视线 (光锥) ----
    # 连接原点到物理环边缘，展示“张角”的概念
    # 为了不显得杂乱，每隔一定角度画一条线
    skip = len(phi) // 16 # 大约画 16 条线
    for i in range(0, len(phi), skip):
        ax.plot([0, x_phys[i]], [0, y_phys[i]], [0, z_phys[i]],
                color='red', linestyle='--', alpha=0.3, linewidth=1)
    # 手动添加一个图例项用于表示这些虚线
    ax.plot([], [], [], color='red', linestyle='--', alpha=0.3, linewidth=1, label='Sightlines (Cone of Vision)')

    # ---- E. 设置场景和标签 ----
    # 标出原点 O
    ax.scatter([0], [0], [0], color="black", s=60, label="Origin O (Observer)")

    # 自动调整坐标轴范围以清晰显示所有元素
    max_range = np.max([np.max(np.abs(x_phys)), np.max(np.abs(y_phys)), z_ring]) * 1.1
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([0, max_range])

    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")

    # 设置一个便于观察的视角
    ax.view_init(elev=25, azim=45)

    # 标题和图例
    plt.title(f"Solid Angle Visualization\nCalculated Omega = {omega_val:.8f} sr", fontsize=14, fontweight='bold')
    # 由于 surface 图例支持不好，我们需要手动创建一个代理图例 (proxy artist)
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='cyan', alpha=0.8, label=f'Solid Angle Projection (Area = {omega_val:.4f} sr)'),
        plt.Line2D([0], [0], color='darkblue', lw=2, label='Projection Boundary'),
        plt.Line2D([0], [0], color='red', lw=3, label=f'Physical Ring Object (z={z_ring})'),
        plt.Line2D([0], [0], color='red', linestyle='--', alpha=0.3, label='Sightlines'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=8, label='Origin O'),
    ]
    ax.legend(handles=legend_elements, loc="best")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # ==========================================
    # 设置输入参数 (使用您提供的不对称案例)
    # ==========================================
    theta_xp = 6.3   # +x 方向角度 (度)
    theta_xm = 6.3   # -x 方向角度 (度)
    theta_yp = 6.3   # +y 方向角度 (度)
    theta_ym = 6.3   # -y 方向角度 (度)

    print("-" * 30)
    print(f"输入角度: X方向 = {theta_xp}°/{theta_xm}°, Y方向 = {theta_yp}°/{theta_ym}°")

    # 1. 调用您的计算函数
    Omega_num, phi, theta_max = solid_angle_elliptic_cone(
        theta_xp, theta_xm, theta_yp, theta_ym
    )
    print(f"精确计算的立体角 = {Omega_num:.8f} sr")
    print("-" * 30)

    # 2. 调用新的绘图函数
    # 将物理环放置在 z=2.0 的位置，以便与单位球 (R=1) 区分开
    plot_improved_geometry(phi, theta_max, Omega_num, z_ring=2.0)