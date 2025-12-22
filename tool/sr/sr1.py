import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def solid_angle_elliptic_cone(theta_xp_deg, theta_xm_deg,
                              theta_yp_deg, theta_ym_deg,
                              n_phi=2000):
    """
    用四个方向(+x, -x, +y, -y)的开角(单位:度)定义一个椭圆锥，
    数值计算其对应的立体角。

    返回:
        Omega: 立体角 (sr)
        phi:   采样的方位角数组 (rad)
        theta_max: 每个 phi 方向上的最大极角 (rad)
    """
    # 角度 -> 弧度
    txp, txm, typ, tym = np.radians(
        [theta_xp_deg, theta_xm_deg, theta_yp_deg, theta_ym_deg]
    )

    # 在切平面上的半轴长度: 取 tan(theta) 并在 ± 方向上做平均
    ax = 0.5 * (np.tan(txp) + np.tan(txm))  # x 方向半轴
    ay = 0.5 * (np.tan(typ) + np.tan(tym))  # y 方向半轴

    # 方位角采样
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)

    # 椭圆在切平面上的极径 r_tan(phi)
    #   r = (a b) / sqrt((b cosφ)^2 + (a sinφ)^2)
    r_tan = ax * ay / np.sqrt((ay * np.cos(phi))**2 +
                              (ax * np.sin(phi))**2)

    # r_tan = tan(theta_max)  =>  theta_max(phi)
    theta_max = np.arctan(r_tan)

    # 立体角: Ω = ∫0^{2π} [1 - cos θ_max(φ)] dφ
    integrand = 1.0 - np.cos(theta_max)
    Omega = np.trapz(integrand, phi)

    return Omega, phi, theta_max


def solid_angle_circular(theta_deg):
    """
    轴对称圆锥 (四个方向角都相同) 的解析立体角。
    """
    theta = np.radians(theta_deg)
    return 2.0 * np.pi * (1.0 - np.cos(theta))


def plot_geometry_with_ring(phi, theta_max,
                            theta_ring_deg,
                            z_ring=1.5,
                            sphere_radius=1.0):
    """
    画:
      - 单位球面
      - 由 phi, theta_max 定义的立体角在球面上的区域
      - z = z_ring 平面上的圆环 (中心在 (0,0,z_ring))
      - 圆环投影到球面上的曲线
      - O(0,0,0) 和 O1(0,0,z_ring)
    """

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # ==== 1. 画单位球面 ====
    u = np.linspace(0, 2 * np.pi, 80)
    v = np.linspace(0, np.pi, 80)
    xs = sphere_radius * np.outer(np.cos(u), np.sin(v))
    ys = sphere_radius * np.outer(np.sin(u), np.sin(v))
    zs = sphere_radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, alpha=0.15, linewidth=0)

    # ==== 2. 画立体角对应球面区域 ====
    n_theta = 60
    theta = np.linspace(0.0, 1.0, n_theta)[:, None] * theta_max[None, :]

    x_patch = np.sin(theta) * np.cos(phi)[None, :] * sphere_radius
    y_patch = np.sin(theta) * np.sin(phi)[None, :] * sphere_radius
    z_patch = np.cos(theta) * sphere_radius

    ax.plot_surface(x_patch, y_patch, z_patch,
                    alpha=0.6, linewidth=0, antialiased=True)

    # 边界曲线（同时也是“圆环在立体角上的投影”的理论边界）
    xb = np.sin(theta_max) * np.cos(phi) * sphere_radius
    yb = np.sin(theta_max) * np.sin(phi) * sphere_radius
    zb = np.cos(theta_max) * sphere_radius
    ax.plot(xb, yb, zb, "k", linewidth=1.5, label="Solid angle boundary")

    # ==== 3. 画实际空间中的圆环 ====
    theta_ring = np.radians(theta_ring_deg)
    # 圆环所在平面 z = z_ring，上面到原点的方向与 z 轴夹角为 theta_ring
    # tan(theta) = r / z  =>  r = z * tan(theta)
    R_ring = z_ring * np.tan(theta_ring)
    t = np.linspace(0, 2 * np.pi, 400)

    xr = R_ring * np.cos(t)
    yr = R_ring * np.sin(t)
    zr = np.full_like(t, z_ring)

    ax.plot(xr, yr, zr, "r", linewidth=2, label="Ring")

    # 圆心 O1
    ax.scatter([0], [0], [z_ring], color="r", s=40)
    ax.text(0, 0, z_ring, " O1", color="r")

    # ==== 4. 把圆环投影到球面上 ====
    v_ring = np.vstack((xr, yr, zr))           # 3 x N
    norms = np.linalg.norm(v_ring, axis=0)     # 长度
    v_unit = v_ring / norms                    # 单位向量
    x_proj = v_unit[0] * sphere_radius
    y_proj = v_unit[1] * sphere_radius
    z_proj = v_unit[2] * sphere_radius

    ax.plot(x_proj, y_proj, z_proj, "g--", linewidth=2,
            label="Projection of ring on sphere")

    # ==== 5. 标出原点 O ====
    ax.scatter([0], [0], [0], color="k", s=40)
    ax.text(0, 0, 0, " O", color="k")

    # 坐标轴指示
    L = 2.0
    ax.plot([0, L], [0, 0], [0, 0], "k")
    ax.text(L, 0, 0, "x", color="k")
    ax.plot([0, 0], [0, L], [0, 0], "k")
    ax.text(0, L, 0, "y", color="k")
    ax.plot([0, 0], [0, 0], [0, L], "k")
    ax.text(0, 0, L, "z", color="k")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect([1, 1, 1])

    ax.legend(loc="upper left")
    plt.title("Ring, its projection, and corresponding solid angle on the sphere")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # ===== 在这里设置四个方向角（单位: 度）=====
    theta_xp = 6.5   # +x 方向
    theta_xm = 6.5   # -x 方向
    theta_yp = 4   # +y 方向
    theta_ym = 4   # -y 方向

    # 数值计算椭圆锥立体角
    Omega_num, phi, theta_max = solid_angle_elliptic_cone(
        theta_xp, theta_xm, theta_yp, theta_ym
    )
    print(f"数值积分立体角  = {Omega_num:.8f} sr")

    # 如果四个角相等，可以用解析公式对照一下
    if abs(theta_xp - theta_xm) < 1e-6 and abs(theta_yp - theta_ym) < 1e-6 \
       and abs(theta_xp - theta_yp) < 1e-6:
        Omega_ana = solid_angle_circular(theta_xp)
        print(f"解析公式立体角  = {Omega_ana:.8f} sr")

    # 圆环的夹角用四个角平均值（如果不完全相等，就取平均）
    theta_ring_deg = (theta_xp + theta_xm + theta_yp + theta_ym) / 4.0

    # 画出几何图像：球 + 立体角区域 + 圆环 + 投影
    plot_geometry_with_ring(phi, theta_max,
                            theta_ring_deg=theta_ring_deg,
                            z_ring=1.5,
                            sphere_radius=1.0)
