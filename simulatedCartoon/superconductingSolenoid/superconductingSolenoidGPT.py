# -*- coding: utf-8 -*-
"""
Monte Carlo simulation of heavy reaction products in a superconducting solenoid
with 100 Pa He gas (small energy loss + multiple scattering).

Author: ChatGPT
"""

import numpy as np
import matplotlib.pyplot as plt

# =========================
# 0) Physical constants
# =========================
E_CHARGE = 1.602176634e-19      # C
AMU = 1.66053906660e-27         # kg
C = 299792458.0                 # m/s
MEV_TO_J = 1.602176634e-13      # J/MeV
MASS_U_MEV = 931.49410242       # MeV/c^2 per u

# =========================
# 1) Utility functions
# =========================
def helium_density_g_cm3(P_pa=100.0, T_K=300.0):
    """
    Ideal-gas helium density at pressure P (Pa) and temperature T (K).
    Return: g/cm^3
    """
    M = 4.002602e-3      # kg/mol (He)
    R = 8.314462618      # J/(mol·K)
    rho_kg_m3 = P_pa * M / (R * T_K)
    return rho_kg_m3 * 1e-3   # kg/m^3 -> g/cm^3


def sample_initial_products(
    N,
    A=270,                    # mass number of product
    E0_MeV=180.0,             # total kinetic energy (NOT MeV/u), example value
    dE_sigma_MeV=10.0,        # energy spread
    q_mean=14.0,              # mean charge state
    q_sigma=2.5,              # charge-state spread
    sigma_xy_mm=1.5,          # target spot sigma in x/y
    sigma_th_mrad=35.0,       # angular divergence sigma
    seed=2026,
):
    """
    Generate initial positions/velocities/energies/charge states of recoils.
    """
    rng = np.random.default_rng(seed)

    x0 = rng.normal(0, sigma_xy_mm * 1e-3, N)
    y0 = rng.normal(0, sigma_xy_mm * 1e-3, N)
    z0 = np.zeros(N)

    # small-angle slopes tx = dx/dz, ty = dy/dz
    tx = rng.normal(0, sigma_th_mrad * 1e-3, N)
    ty = rng.normal(0, sigma_th_mrad * 1e-3, N)

    E0 = np.clip(rng.normal(E0_MeV, dE_sigma_MeV, N), 1.0, None)
    q_state = np.clip(np.rint(rng.normal(q_mean, q_sigma, N)).astype(int), 1, None)

    m_kg = A * AMU
    v = np.sqrt(2.0 * E0 * MEV_TO_J / m_kg)  # non-rel approximation

    denom = np.sqrt(1.0 + tx**2 + ty**2)
    ux, uy, uz = tx / denom, ty / denom, 1.0 / denom

    vx0 = v * ux
    vy0 = v * uy
    vz0 = v * uz

    return dict(x0=x0, y0=y0, z0=z0, vx0=vx0, vy0=vy0, vz0=vz0, E0=E0, q=q_state)


def simulate_transport(
    init,
    A=270,
    Bz_T=4.0,                 # solenoid B field (uniform)
    L_sol_m=1.2,              # solenoid length
    L_drift_m=0.8,            # drift after solenoid (to show focusing image plane)
    dz_m=0.002,               # integration z-step
    P_he_pa=100.0,
    T_he_K=300.0,
    S_mass_MeV_cm2_g=5000.0,  # effective mass stopping power (tunable)
    include_scattering=True,
    seed=1,
):
    """
    Vectorized step-by-step transport:
    - exact transverse rotation in uniform Bz (inside solenoid)
    - continuous dE/dx in He (inside solenoid)
    - small-angle multiple scattering in He (inside solenoid)
    """
    rng = np.random.default_rng(seed)

    # Copy initial state
    x = init["x0"].copy()
    y = init["y0"].copy()
    z = init["z0"].copy()
    vx = init["vx0"].copy()
    vy = init["vy0"].copy()
    vz = init["vz0"].copy()
    E = init["E0"].copy()
    q_state = init["q"].copy()

    N = len(x)
    m_kg = A * AMU
    Mc2_MeV = A * MASS_U_MEV

    L_total = L_sol_m + L_drift_m
    n_steps = int(np.ceil(L_total / dz_m))
    dz_m = L_total / n_steps  # make endpoint exact

    rho_g_cm3 = helium_density_g_cm3(P_he_pa, T_he_K)
    X0_He_g_cm2 = 94.32  # approx radiation length of He (g/cm^2)

    z_axis = np.linspace(dz_m, L_total, n_steps)
    sigma_r = np.empty(n_steps)  # beam envelope indicator
    meanE = np.empty(n_steps)

    alive = np.ones(N, dtype=bool)

    for i in range(n_steps):
        if not np.any(alive):
            sigma_r[i:] = np.nan
            meanE[i:] = np.nan
            break

        idx = np.where(alive)[0]

        # force forward propagation in z
        bad = vz[idx] <= 1e-6
        if np.any(bad):
            alive[idx[bad]] = False
            idx = np.where(alive)[0]
            if len(idx) == 0:
                sigma_r[i:] = np.nan
                meanE[i:] = np.nan
                break

        # choose dt so each particle advances ~dz in z
        dt = dz_m / vz[idx]

        # Solenoid field region mask
        in_sol = z[idx] < L_sol_m

        # -------- Lorentz motion in uniform Bz (exact rotation in transverse plane) --------
        if Bz_T != 0.0:
            phi = (q_state[idx] * E_CHARGE * Bz_T / m_kg) * dt  # cyclotron angle per step
            cph = np.cos(phi)
            sph = np.sin(phi)

            vx_old = vx[idx].copy()
            vy_old = vy[idx].copy()

            # rotate only inside solenoid
            vx_new = np.where(in_sol, cph * vx_old - sph * vy_old, vx_old)
            vy_new = np.where(in_sol, sph * vx_old + cph * vy_old, vy_old)

            vx[idx] = vx_new
            vy[idx] = vy_new

        # position update
        x[idx] += vx[idx] * dt
        y[idx] += vy[idx] * dt
        z[idx] += vz[idx] * dt  # approximately = dz_m

        # -------- He gas effects (only inside solenoid) --------
        if np.any(in_sol):
            id_sol = idx[in_sol]

            speed = np.sqrt(vx[id_sol]**2 + vy[id_sol]**2 + vz[id_sol]**2)
            ds_cm = speed * (dz_m / np.maximum(vz[id_sol], 1e-12)) * 100.0

            # (1) Continuous energy loss: dE = (S_mass) * rho * ds
            dE = S_mass_MeV_cm2_g * rho_g_cm3 * ds_cm
            E[id_sol] = np.maximum(E[id_sol] - dE, 0.2)

            # match speed to updated E
            speed_new = np.sqrt(2.0 * E[id_sol] * MEV_TO_J / m_kg)

            if include_scattering:
                # (2) Small-angle multiple scattering (Highland-like, approximate)
                x_over_X0 = np.clip(rho_g_cm3 * ds_cm / X0_He_g_cm2, 1e-20, None)

                # non-rel p*c ~ sqrt(2 M c^2 E)
                pc_MeV = np.sqrt(np.clip(2.0 * Mc2_MeV * E[id_sol], 1e-12, None))
                beta = np.clip(speed_new / C, 1e-6, 0.99)

                theta0 = (13.6 / (beta * pc_MeV)) * np.abs(q_state[id_sol]) * np.sqrt(x_over_X0)
                theta0 *= (1.0 + 0.038 * np.log(np.maximum(x_over_X0, 1e-12)))
                theta0 = np.clip(theta0, 0.0, 0.05)  # avoid pathological large kicks

                tx = vx[id_sol] / np.maximum(vz[id_sol], 1e-12)
                ty = vy[id_sol] / np.maximum(vz[id_sol], 1e-12)

                # isotropic split to x/y
                tx += rng.normal(0.0, theta0 / np.sqrt(2), size=len(id_sol))
                ty += rng.normal(0.0, theta0 / np.sqrt(2), size=len(id_sol))

                denom = np.sqrt(1.0 + tx**2 + ty**2)
                ux, uy, uz = tx / denom, ty / denom, 1.0 / denom

                vx[id_sol] = speed_new * ux
                vy[id_sol] = speed_new * uy
                vz[id_sol] = speed_new * uz
            else:
                # only speed rescale, keep direction
                speed_old = np.sqrt(vx[id_sol]**2 + vy[id_sol]**2 + vz[id_sol]**2)
                scale = speed_new / np.maximum(speed_old, 1e-12)
                vx[id_sol] *= scale
                vy[id_sol] *= scale
                vz[id_sol] *= scale

        # remove backward particles if any
        alive &= (vz > 1e-5)

        sigma_r[i] = np.std(np.hypot(x, y))
        meanE[i] = np.mean(E)

    return {
        "x": x, "y": y, "z": z,
        "vx": vx, "vy": vy, "vz": vz,
        "E": E, "q": q_state,
        "z_axis": z_axis,
        "sigma_r": sigma_r,
        "meanE": meanE,
        "rho_g_cm3": rho_g_cm3,
        "dz_m": dz_m,
        "n_steps": n_steps,
    }


# =========================
# 2) Main demo
# =========================
if __name__ == "__main__":
    # -------- User-tunable parameters --------
    # 你可以按自己的实验条件改这里
    beam_params = dict(
        N=4000,              # Monte Carlo particles
        A=270,               # e.g. ER around superheavy region
        E0_MeV=180.0,        # total kinetic energy of products
        dE_sigma_MeV=10.0,
        q_mean=14.0,
        q_sigma=2.5,
        sigma_xy_mm=1.5,     # source spot size
        sigma_th_mrad=35.0,  # angular spread
        seed=2026,
    )

    transport_common = dict(
        A=beam_params["A"],
        L_sol_m=1.2,                 # solenoid length
        L_drift_m=0.8,               # detector plane after drift
        dz_m=0.002,
        P_he_pa=100.0,               # He pressure = 100 Pa
        T_he_K=300.0,
        S_mass_MeV_cm2_g=5000.0,     # effective stopping power (IMPORTANT: tune with SRIM/ATIMA)
        include_scattering=True,
    )

    B_focus_T = 4.0   # try 2~6 T to see focusing strength change

    # -------- Initial products --------
    init = sample_initial_products(**beam_params)

    # Compare no-field vs solenoid focusing
    res_noB = simulate_transport(init, Bz_T=0.0, seed=1, **transport_common)
    res_B   = simulate_transport(init, Bz_T=B_focus_T, seed=1, **transport_common)

    # -------- Derived quantities --------
    r0_mm     = np.hypot(init["x0"], init["y0"]) * 1e3
    r_noB_mm  = np.hypot(res_noB["x"], res_noB["y"]) * 1e3
    r_B_mm    = np.hypot(res_B["x"], res_B["y"]) * 1e3

    dE_noB = init["E0"] - res_noB["E"]
    dE_B   = init["E0"] - res_B["E"]

    print("=" * 70)
    print("He gas density @100 Pa, 300 K = %.4e g/cm^3" % res_B["rho_g_cm3"])
    print("Mean energy loss (B=0 T)   = %.4f MeV" % np.mean(dE_noB))
    print("Mean energy loss (B=%.1f T) = %.4f MeV" % (B_focus_T, np.mean(dE_B)))
    print("Detector-plane RMS radius:")
    print("  B=0 T      : %.3f mm" % np.std(r_noB_mm))
    print("  B=%.1f T    : %.3f mm" % (B_focus_T, np.std(r_B_mm)))
    print("Mean radius:")
    print("  B=0 T      : %.3f mm" % np.mean(r_noB_mm))
    print("  B=%.1f T    : %.3f mm" % (B_focus_T, np.mean(r_B_mm)))
    print("=" * 70)

    # =========================
    # 3) Visualization
    # =========================
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # (a) Initial distribution
    ax = axes[0, 0]
    ax.scatter(init["x0"] * 1e3, init["y0"] * 1e3, s=2, alpha=0.35)
    ax.set_title("Initial product distribution (z=0)")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)

    # (b) Detector plane x-y (overlay)
    ax = axes[0, 1]
    ax.scatter(res_noB["x"] * 1e3, res_noB["y"] * 1e3, s=2, alpha=0.25, label="B=0 T")
    ax.scatter(res_B["x"] * 1e3,   res_B["y"] * 1e3,   s=2, alpha=0.25, label=f"B={B_focus_T:.1f} T")
    ax.set_title("Distribution at detector plane (after solenoid + drift)")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend(markerscale=3)

    # (c) Beam envelope (RMS radius vs z)
    ax = axes[1, 0]
    ax.plot(res_noB["z_axis"], res_noB["sigma_r"] * 1e3, label="B=0 T")
    ax.plot(res_B["z_axis"],   res_B["sigma_r"] * 1e3,   label=f"B={B_focus_T:.1f} T")
    ax.axvspan(0.0, transport_common["L_sol_m"], alpha=0.12, label="Solenoid region")
    ax.set_title("Beam envelope evolution (RMS radius)")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("RMS(r) [mm]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # (d) Radial distribution at detector plane
    ax = axes[1, 1]
    ax.hist(r_noB_mm, bins=60, histtype="step", linewidth=1.8, label="B=0 T")
    ax.hist(r_B_mm,   bins=60, histtype="step", linewidth=1.8, label=f"B={B_focus_T:.1f} T")
    ax.set_title("Radial distribution at detector plane")
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("Counts")
    ax.grid(True, alpha=0.3)
    ax.legend()

    txt = (
        f"He = 100 Pa\n"
        f"Lsol = {transport_common['L_sol_m']} m\n"
        f"<ΔE> ≈ {np.mean(dE_B):.3f} MeV\n"
        f"RMS(r): {np.std(r_noB_mm):.2f} -> {np.std(r_B_mm):.2f} mm"
    )
    ax.text(
        0.98, 0.95, txt,
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(boxstyle="round", alpha=0.15)
    )

    fig.tight_layout()
    plt.show()