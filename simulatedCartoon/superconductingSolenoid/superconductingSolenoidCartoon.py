# -*- coding: utf-8 -*-
"""
Trajectory visualization of reaction products in a superconducting solenoid
with 100 Pa He gas (continuous dE/dx + small-angle scattering, simplified).

Focus: visualize particle paths (trajectory), not just final spot comparison.
"""

import numpy as np
import matplotlib.pyplot as plt

# 如果你想看3D轨迹，把下面注释打开
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# =========================================================
# 1) 常量
# =========================================================
E_CHARGE = 1.602176634e-19      # C
AMU = 1.66053906660e-27         # kg
C = 299792458.0                 # m/s
MEV_TO_J = 1.602176634e-13      # J/MeV
MASS_U_MEV = 931.49410242       # MeV/c^2/u


# =========================================================
# 2) He气体密度（理想气体）
# =========================================================
def helium_density_g_cm3(P_pa=100.0, T_K=300.0):
    M = 4.002602e-3  # kg/mol
    R = 8.314462618
    rho_kg_m3 = P_pa * M / (R * T_K)
    return rho_kg_m3 * 1e-3  # kg/m^3 -> g/cm^3


# =========================================================
# 3) 初始粒子采样（打靶产物）
# =========================================================
def sample_products(
    N=300,
    A=270,
    E0_MeV=180.0,          # 总动能（例子）
    E_sigma_MeV=8.0,
    q_mean=14.0,
    q_sigma=2.0,
    sigma_xy_mm=1.0,       # 靶点尺寸
    sigma_th_mrad=30.0,    # 发散角
    seed=2026
):
    rng = np.random.default_rng(seed)

    # 初始位置（z=0）
    x0 = rng.normal(0.0, sigma_xy_mm * 1e-3, N)
    y0 = rng.normal(0.0, sigma_xy_mm * 1e-3, N)
    z0 = np.zeros(N)

    # 小角近似：tx=dx/dz, ty=dy/dz
    tx = rng.normal(0.0, sigma_th_mrad * 1e-3, N)
    ty = rng.normal(0.0, sigma_th_mrad * 1e-3, N)

    E0 = np.clip(rng.normal(E0_MeV, E_sigma_MeV, N), 1.0, None)
    q = np.clip(np.rint(rng.normal(q_mean, q_sigma, N)).astype(int), 1, None)

    m_kg = A * AMU
    v = np.sqrt(2.0 * E0 * MEV_TO_J / m_kg)  # 非相对论近似（重离子常用快速可视化可接受）

    denom = np.sqrt(1.0 + tx**2 + ty**2)
    ux, uy, uz = tx / denom, ty / denom, 1.0 / denom

    vx0 = v * ux
    vy0 = v * uy
    vz0 = v * uz

    return {
        "x": x0, "y": y0, "z": z0,
        "vx": vx0, "vy": vy0, "vz": vz0,
        "E": E0, "q": q,
        "A": A
    }


# =========================================================
# 4) 轨迹推进（记录轨迹）
# =========================================================
def propagate_with_trajectory(
    state,
    Bz_T=4.0,                 # 螺线管均匀Bz（简化）
    L_sol_m=1.2,              # 螺线管长度
    L_drift_m=0.5,            # 后漂移段
    dz_m=0.002,               # z步长
    P_he_pa=100.0,
    T_he_K=300.0,
    S_mass_MeV_cm2_g=5000.0,  # 有效质量阻止本领（简化参数，建议后续用SRIM替换）
    include_scattering=True,
    scatter_scale=1.0,        # 可调散射强度缩放，便于视觉展示
    n_show_tracks=40,         # 要显示多少条轨迹
    seed=7
):
    """
    返回:
      final_state: 最终状态
      traj: 被选中粒子的轨迹数组 (nstep+1, n_show)
    """
    rng = np.random.default_rng(seed)

    # 拷贝状态，避免修改原始输入
    x = state["x"].copy()
    y = state["y"].copy()
    z = state["z"].copy()
    vx = state["vx"].copy()
    vy = state["vy"].copy()
    vz = state["vz"].copy()
    E = state["E"].copy()
    q = state["q"].copy()

    A = state["A"]
    N = len(x)
    m_kg = A * AMU
    Mc2_MeV = A * MASS_U_MEV

    rho_g_cm3 = helium_density_g_cm3(P_he_pa, T_he_K)
    X0_He_g_cm2 = 94.32  # He辐射长度，近似用于散射量级估计

    L_total = L_sol_m + L_drift_m
    n_steps = int(np.ceil(L_total / dz_m))
    dz_m = L_total / n_steps

    # 随机挑一些粒子显示轨迹
    n_show_tracks = min(n_show_tracks, N)
    show_idx = rng.choice(N, size=n_show_tracks, replace=False)

    # 轨迹记录
    traj = {
        "idx": show_idx,
        "x": np.zeros((n_steps + 1, n_show_tracks)),
        "y": np.zeros((n_steps + 1, n_show_tracks)),
        "z": np.zeros((n_steps + 1, n_show_tracks)),
        "E": np.zeros((n_steps + 1, n_show_tracks)),
    }
    traj["x"][0, :] = x[show_idx]
    traj["y"][0, :] = y[show_idx]
    traj["z"][0, :] = z[show_idx]
    traj["E"][0, :] = E[show_idx]

    alive = np.ones(N, dtype=bool)

    for i in range(1, n_steps + 1):
        idx = np.where(alive)[0]
        if len(idx) == 0:
            # 所有粒子都“失活”了（比如倒飞）
            traj["x"][i:, :] = traj["x"][i - 1, :]
            traj["y"][i:, :] = traj["y"][i - 1, :]
            traj["z"][i:, :] = traj["z"][i - 1, :]
            traj["E"][i:, :] = traj["E"][i - 1, :]
            break

        # 防止vz<=0
        bad = vz[idx] <= 1e-8
        if np.any(bad):
            alive[idx[bad]] = False
            idx = np.where(alive)[0]
            if len(idx) == 0:
                traj["x"][i:, :] = traj["x"][i - 1, :]
                traj["y"][i:, :] = traj["y"][i - 1, :]
                traj["z"][i:, :] = traj["z"][i - 1, :]
                traj["E"][i:, :] = traj["E"][i - 1, :]
                break

        dt = dz_m / vz[idx]  # 每步沿z前进约dz_m
        in_sol = z[idx] < L_sol_m

        # -----------------------------
        # (a) 磁场作用：均匀Bz下横向速度旋转
        # -----------------------------
        if Bz_T != 0.0:
            phi = (q[idx] * E_CHARGE * Bz_T / m_kg) * dt  # 每步回旋角
            cph = np.cos(phi)
            sph = np.sin(phi)

            vx_old = vx[idx].copy()
            vy_old = vy[idx].copy()

            # 仅在螺线管内部旋转
            vx[idx] = np.where(in_sol, cph * vx_old - sph * vy_old, vx_old)
            vy[idx] = np.where(in_sol, sph * vx_old + cph * vy_old, vy_old)

        # 位置推进
        x[idx] += vx[idx] * dt
        y[idx] += vy[idx] * dt
        z[idx] += vz[idx] * dt

        # -----------------------------
        # (b) 100 Pa He: 能损 + 小角散射（仅螺线管区）
        # -----------------------------
        if np.any(in_sol):
            id_sol = idx[in_sol]

            speed = np.sqrt(vx[id_sol]**2 + vy[id_sol]**2 + vz[id_sol]**2)
            # 实际路径 ds = v*dt ; 转为cm
            ds_cm = speed * (dz_m / np.maximum(vz[id_sol], 1e-12)) * 100.0

            # 连续能损（有效近似）
            dE = S_mass_MeV_cm2_g * rho_g_cm3 * ds_cm
            E[id_sol] = np.maximum(E[id_sol] - dE, 0.2)

            # 用新能量更新速率
            speed_new = np.sqrt(2.0 * E[id_sol] * MEV_TO_J / m_kg)

            # 散射（Highland-like，给趋势图足够）
            if include_scattering:
                x_over_X0 = np.clip(rho_g_cm3 * ds_cm / X0_He_g_cm2, 1e-20, None)
                pc_MeV = np.sqrt(np.clip(2.0 * Mc2_MeV * E[id_sol], 1e-12, None))
                beta = np.clip(speed_new / C, 1e-6, 0.99)

                theta0 = (13.6 / (beta * pc_MeV)) * np.abs(q[id_sol]) * np.sqrt(x_over_X0)
                theta0 *= (1.0 + 0.038 * np.log(np.maximum(x_over_X0, 1e-12)))
                theta0 *= scatter_scale
                theta0 = np.clip(theta0, 0.0, 0.08)

                tx = vx[id_sol] / np.maximum(vz[id_sol], 1e-12)
                ty = vy[id_sol] / np.maximum(vz[id_sol], 1e-12)

                tx += rng.normal(0.0, theta0 / np.sqrt(2), size=len(id_sol))
                ty += rng.normal(0.0, theta0 / np.sqrt(2), size=len(id_sol))

                denom = np.sqrt(1.0 + tx**2 + ty**2)
                ux, uy, uz = tx / denom, ty / denom, 1.0 / denom
                vx[id_sol] = speed_new * ux
                vy[id_sol] = speed_new * uy
                vz[id_sol] = speed_new * uz
            else:
                # 不加散射，仅缩放速率，方向不变
                speed_old = np.sqrt(vx[id_sol]**2 + vy[id_sol]**2 + vz[id_sol]**2)
                scale = speed_new / np.maximum(speed_old, 1e-12)
                vx[id_sol] *= scale
                vy[id_sol] *= scale
                vz[id_sol] *= scale

        # 删除倒飞粒子（少量）
        alive &= (vz > 1e-8)

        # 记录显示轨迹
        traj["x"][i, :] = x[show_idx]
        traj["y"][i, :] = y[show_idx]
        traj["z"][i, :] = z[show_idx]
        traj["E"][i, :] = E[show_idx]

    final_state = {
        "x": x, "y": y, "z": z,
        "vx": vx, "vy": vy, "vz": vz,
        "E": E, "q": q,
        "rho_g_cm3": rho_g_cm3
    }
    return final_state, traj


# =========================================================
# 5) 绘图：轨迹可视化（核心）
# =========================================================
def plot_trajectories(traj, final_state, Bz_T, L_sol_m, L_drift_m, show_3d=True):
    x_mm = traj["x"] * 1e3
    y_mm = traj["y"] * 1e3
    z_m = traj["z"]

    n_show = x_mm.shape[1]
    z_end = L_sol_m + L_drift_m

    # ---- 主图：x-z 和 y-z 轨迹 ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # x-z
    ax = axes[0, 0]
    for j in range(n_show):
        ax.plot(z_m[:, j], x_mm[:, j], linewidth=1.0, alpha=0.8)
    ax.axvspan(0, L_sol_m, alpha=0.12)
    ax.axvline(L_sol_m, linestyle='--', linewidth=1)
    ax.set_title(f"x-z trajectories (B = {Bz_T:.2f} T)")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [mm]")
    ax.grid(True, alpha=0.3)

    # y-z
    ax = axes[0, 1]
    for j in range(n_show):
        ax.plot(z_m[:, j], y_mm[:, j], linewidth=1.0, alpha=0.8)
    ax.axvspan(0, L_sol_m, alpha=0.12)
    ax.axvline(L_sol_m, linestyle='--', linewidth=1)
    ax.set_title("y-z trajectories")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("y [mm]")
    ax.grid(True, alpha=0.3)

    # x-y 投影（终点）
    ax = axes[1, 0]
    ax.scatter(final_state["x"] * 1e3, final_state["y"] * 1e3, s=4, alpha=0.35)
    ax.set_title("Final x-y distribution (detector plane)")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)

    # 能量随z（只显示被跟踪轨迹）
    ax = axes[1, 1]
    for j in range(n_show):
        ax.plot(z_m[:, j], traj["E"][:, j], linewidth=1.0, alpha=0.8)
    ax.axvspan(0, L_sol_m, alpha=0.12, label="Solenoid + He region")
    ax.axvline(L_sol_m, linestyle='--', linewidth=1)
    ax.set_title("Energy along trajectory (tracked particles)")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("E [MeV]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    rho = final_state["rho_g_cm3"]
    fig.suptitle(
        f"Trajectory visualization in superconducting solenoid (He = 100 Pa, rho ≈ {rho:.2e} g/cm^3)\n"
        f"Shaded region: solenoid (0 ~ {L_sol_m:.2f} m), total z = {z_end:.2f} m",
        y=0.98
    )
    fig.tight_layout()
    plt.show()

    # ---- 可选：3D轨迹图 ----
    if show_3d:
        fig3d = plt.figure(figsize=(10, 8))
        ax3d = fig3d.add_subplot(111, projection='3d')

        for j in range(min(n_show, 30)):  # 3D图别画太多，避免太乱
            ax3d.plot(z_m[:, j], x_mm[:, j], y_mm[:, j], linewidth=1.0, alpha=0.9)

        # 画螺线管区域边界（示意）
        # 用一个简单圆柱轮廓示意，不代表真实孔径
        r_cyl_mm = max(10, np.nanpercentile(np.sqrt(x_mm**2 + y_mm**2), 90) * 1.2)
        theta = np.linspace(0, 2*np.pi, 80)
        z1 = np.zeros_like(theta)
        z2 = np.full_like(theta, L_sol_m)
        xc = r_cyl_mm * np.cos(theta)
        yc = r_cyl_mm * np.sin(theta)
        ax3d.plot(z1, xc, yc, linewidth=1.0)
        ax3d.plot(z2, xc, yc, linewidth=1.0)

        for k in [0, 20, 40, 60]:
            ax3d.plot([0, L_sol_m], [xc[k], xc[k]], [yc[k], yc[k]], linewidth=0.8, alpha=0.7)

        ax3d.set_title("3D trajectories (z, x, y)")
        ax3d.set_xlabel("z [m]")
        ax3d.set_ylabel("x [mm]")
        ax3d.set_zlabel("y [mm]")
        plt.tight_layout()
        plt.show()


# =========================================================
# 6) 主程序
# =========================================================
if __name__ == "__main__":
    # ---------- 你最常改的参数 ----------
    N_particles = 600          # 总粒子数（用于统计+背景）
    N_show_tracks = 50         # 真正画轨迹的条数（太多会乱）
    A_product = 270            # 重产物质量数（示意）

    Bz_T = 4.0                 # 超导螺线管磁场
    L_sol_m = 1.2
    L_drift_m = 0.6
    dz_m = 0.002

    # He = 100 Pa（按你的要求）
    P_he_pa = 100.0
    T_he_K = 300.0

    # 有效能损（简化参数）——你后续可用SRIM替换
    S_mass_MeV_cm2_g = 5000.0

    # ---------- 初态采样 ----------
    state0 = sample_products(
        N=N_particles,
        A=A_product,
        E0_MeV=180.0,
        E_sigma_MeV=10.0,
        q_mean=14.0,
        q_sigma=2.0,
        sigma_xy_mm=1.2,
        sigma_th_mrad=35.0,
        seed=2026
    )

    # ---------- 传播并记录轨迹 ----------
    final_state, traj = propagate_with_trajectory(
        state0,
        Bz_T=Bz_T,
        L_sol_m=L_sol_m,
        L_drift_m=L_drift_m,
        dz_m=dz_m,
        P_he_pa=P_he_pa,
        T_he_K=T_he_K,
        S_mass_MeV_cm2_g=S_mass_MeV_cm2_g,
        include_scattering=True,
        scatter_scale=1.0,       # 如果你想更明显看到He散射影响，可临时调到2~3（展示用）
        n_show_tracks=N_show_tracks,
        seed=7
    )

    # ---------- 打印一点信息 ----------
    mean_dE = np.mean(state0["E"] - final_state["E"])
    print("=" * 70)
    print(f"He @ {P_he_pa:.1f} Pa, {T_he_K:.1f} K: rho = {final_state['rho_g_cm3']:.4e} g/cm^3")
    print(f"Mean energy loss in transport = {mean_dE:.4f} MeV")
    print(f"Displayed trajectories = {len(traj['idx'])}, Total particles = {N_particles}")
    print("=" * 70)

    # ---------- 轨迹可视化（核心） ----------
    plot_trajectories(
        traj=traj,
        final_state=final_state,
        Bz_T=Bz_T,
        L_sol_m=L_sol_m,
        L_drift_m=L_drift_m,
        show_3d=True
    )