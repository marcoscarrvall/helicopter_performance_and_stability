import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sys
sys.path.insert(0, '/mnt/user-data/uploads')
from data import ApacheAH64 as AH64

# Unpack parameters
rho          = AH64.flight_condition["rho"]
W            = AH64.helicopter["W"]
R            = AH64.main_rotor["R"]
Omega        = AH64.main_rotor["omega"]
sigma        = AH64.main_rotor["sigma"]
k            = AH64.main_rotor["k"]
R_tr         = AH64.tail_rotor["R_tr"]
Omega_tr     = AH64.tail_rotor["omega_tr"]
sigma_tr     = AH64.tail_rotor["sigma_tr"]
k_tr         = AH64.tail_rotor["k_tr"]
fin_blockage = AH64.tail_rotor["fin_blockage"]
l_tr         = AH64.tail_rotor["l_tr"]
V_max        = AH64.specs["v_max"]
V_cruise      = AH64.specs["v_cruise"]
FM           = AH64.main_rotor["FM"]
CDS_f        = AH64.helicopter["flat_plate_area"]
Cd_p         = AH64.main_rotor["cd_p"]
Cd_p_tr      = AH64.tail_rotor["cd_p_tr"]

A            = np.pi * R**2
A_tr         = np.pi * R_tr**2
tip_speed    = Omega    * R
tip_speed_tr = Omega_tr * R_tr

# 1. Hover induced velocity
def compute_hover_induced_velocity(W, rho, A):
    T    = W
    vi_h = np.sqrt(T / (2 * rho * A))
    print(f"Hover induced velocity  vi_h = {vi_h:.4f} m/s")
    return vi_h


# 2. Hover power calculations
def compute_hover_powers(W, vi_h, FM, k, sigma, Cd_p, rho, tip_speed, R):
    T         = W
    P_ideal   = W * vi_h
    P_act     = P_ideal / FM
    P_induced = k * T * vi_h
    P_profile = (sigma * Cd_p / 8) * rho * tip_speed**3 * np.pi * R**2
    P_bem     = P_induced + P_profile

    print(f"\nHover Powers:")
    print(f"  Ideal (ACT)  P_id  = {P_ideal/1000:.3f} kW")
    print(f"  ACT (FM={FM}) P_hov = {P_act/1000:.3f} kW")
    print(f"  BEM (main rotor only) P_bem = {P_bem/1000:.3f} kW")

    return P_ideal, P_act, P_bem

# 3. Find mean lift coefficients
def compute_mean_lift_coefficients(W, P_mr, rho, A_mr, A_tr, Omega_mr, Omega_tr, 
                                    R_mr, R_tr, sigma_mr, sigma_tr, l_tr):

    C_T_hov_mr = W / (rho * A_mr * (Omega_mr * R_mr)**2)
    C_L_mean_hov_mr = 6.6*C_T_hov_mr / sigma_mr

    Q_mr = P_mr / Omega_mr
    T_tr = Q_mr / l_tr
    C_T_hov_tr = T_tr / (rho * A_tr * (Omega_tr * R_tr)**2)
    C_L_mean_hov_tr = 6.6*C_T_hov_tr / sigma_tr

    print(f"\nMean Lift Coefficients at Hover:")
    print(f"  Main Rotor: C_T = {C_T_hov_mr:.4f}, C_L = {C_L_mean_hov_mr:.4f}")
    print(f"  Tail Rotor: C_T = {C_T_hov_tr:.4f}, C_L = {C_L_mean_hov_tr:.4f}")

    return C_T_hov_mr, C_L_mean_hov_mr, C_T_hov_tr, C_L_mean_hov_tr

# A. Glauert induced velocity solver
def solve_vi_forward(V, alpha_d, T_thrust, rho, A, vi_h):
    """
    Solve Glauert's 4th-order implicit equation for forward-flight induced
    velocity . Non-dimensionalised residual:
        vi^2 * [(V cos a)^2 + (V sin a + vi)^2] - vi_h^4 = 0

    Returns:
        vi [m/s]: induced velocity at flight speed V
    """
    def residual(vi):
        VR2 = (V * np.cos(alpha_d))**2 + (V * np.sin(alpha_d) + vi)**2
        return vi**2 * VR2 - vi_h**4

    vi0    = vi_h if V < 1.0 else T_thrust / (2 * rho * A * V)
    vi_sol = fsolve(residual, vi0, full_output=False)[0]
    return abs(vi_sol)

# 4. Forward flight sweep
def compute_forward_flight(V_max, W, rho, A, A_tr, vi_h,
                            CDS_f, k, sigma, Cd_p, tip_speed, Omega, R,
                            sigma_tr, Cd_p_tr, tip_speed_tr, Omega_tr, R_tr,
                            k_tr, fin_blockage, l_tr, N=200):
    """
    Sweep from V=0 to V_max and compute all BEM power components (Eq. 12):
        P = Ppar + Pi + Pp + Pd + Ptr

    Returns:
        V_arr        [m/s]  airspeed array
        P_par_arr    [W]    fuselage parasite power
        P_i_arr      [W]    main rotor induced power
        P_prof_arr   [W]    main rotor profile + H-force power
        P_tr_arr     [W]    tail rotor total power
        P_tot_arr    [W]    total power required
    """
    V_arr = np.append(np.linspace(0.1, V_max, N), V_cruise)
    V_arr = np.sort(V_arr)
    P_par_arr  = np.zeros(len(V_arr))
    P_i_arr    = np.zeros(len(V_arr))
    P_prof_arr = np.zeros(len(V_arr))
    P_tr_arr   = np.zeros(len(V_arr))
    P_tot_arr  = np.zeros(len(V_arr))
    v_i_cruise = 0

    for idx, V in enumerate(V_arr):
        # 4a. Disc angle of attack: sin(alpha_d) = D_fus / W  (Eqs. 5-7)
        D_fus   = 0.5 * rho * V**2 * CDS_f
        alpha_d = np.arcsin(np.clip(D_fus / W, -1, 1))

        # 4b. Main rotor thrust (level flight)
        T_mr = W / np.cos(alpha_d)

        # 4c. Induced velocity (Glauert, Eq. 4)
        vi   = solve_vi_forward(V, alpha_d, T_mr, rho, A, vi_h)

        # 4d. Advance ratio
        mu   = (V * np.cos(alpha_d)) / tip_speed

        # 4e. Main rotor power components
        P_par  = 0.5 * rho * V**3 * CDS_f
        P_i    = k * T_mr * vi
        P_prof = (1/8) * rho * sigma * np.pi * R**2 * tip_speed**3 * Cd_p * (1 + 4.65 * mu**2)

        # 4f. Tail rotor power
        Q_mr    = (P_i + P_prof) / Omega
        T_tr    = Q_mr / l_tr
        vi_h_tr = np.sqrt(abs(T_tr) / (2 * rho * A_tr))
        mu_tr   = V / tip_speed_tr

        def residual_tr(vi_tr):
            VR2 = V**2 + vi_tr**2
            return vi_tr**2 * VR2 - vi_h_tr**4

        vi_tr0 = vi_h_tr if V < 1.0 else T_tr / (2 * rho * A_tr * V)
        vi_tr  = abs(fsolve(residual_tr, max(vi_tr0, 0.01))[0])

        P_i_tr = fin_blockage * k_tr * T_tr * vi_tr
        P_p_tr = (1/8) * rho * sigma_tr * np.pi * R_tr**2 * tip_speed_tr**3 * Cd_p_tr * (1 + 4.65 * mu_tr**2)
        P_tr   = P_i_tr + P_p_tr

        # 4g. Total power
        P_par_arr[idx]  = P_par
        P_i_arr[idx]    = P_i
        P_prof_arr[idx] = P_prof
        P_tr_arr[idx]   = P_tr
        P_tot_arr[idx]  = P_par + P_i + P_prof + P_tr
        if V == V_cruise:
            v_i_cruise = vi
    print(f"\nPower at ccruise speed (V={V_cruise:.2f} m/s):")
    print(f"  P_par = {P_par_arr[idx]/1000:.2f} kW")
    print(f"  P_i = {P_i_arr[idx]/1000:.2f} kW")
    print(f"  P_prof = {P_prof_arr[idx]/1000:.2f} kW")
    print(f"  P_tr = {P_tr_arr[idx]/1000:.2f} kW")
    print(f"  P_total = {P_tot_arr[idx]/1000:.2f} kW")
    print(f"  v_i at cruise speed (V={V_cruise:.2f} m/s) = {v_i_cruise:.4f} m/s")

    return V_arr, P_par_arr, P_i_arr, P_prof_arr, P_tr_arr, P_tot_arr

# 5. Best endurance & best range speeds
def compute_best_speeds(V_arr, P_tot_arr):
    idx_E = np.argmin(P_tot_arr)
    idx_R = np.argmin(P_tot_arr / V_arr)

    results = {
        "V_endurance": V_arr[idx_E], "P_endurance": P_tot_arr[idx_E],
        "V_range":     V_arr[idx_R], "P_range":     P_tot_arr[idx_R],
    }

    V_cruise = AH64.specs["v_cruise"]
    print(f"\nForward Flight Results:")
    print(f"  V_best_endurance = {results['V_endurance']:.2f} m/s  ({results['V_endurance']*3.6:.1f} km/h)  P = {results['P_endurance']/1000:.2f} kW")
    print(f"  V_best_range     = {results['V_range']:.2f} m/s  ({results['V_range']*3.6:.1f} km/h)  P = {results['P_range']/1000:.2f} kW")
    print(f"\n  (Manufacturer cruise speed = {V_cruise:.2f} m/s = {V_cruise*3.6:.1f} km/h)")
    print(f"  (Manufacturer max speed    = {V_max:.2f} m/s = {V_max*3.6:.1f} km/h)")
    return results

# 6. Plot
def plot_results(V_arr, P_par_arr, P_i_arr, P_prof_arr, P_tr_arr, P_tot_arr,
                 best_speeds, V_max, save_path="figures/forward_flight_performance.png"):
    V_kmh  = V_arr * 3.6
    V_E    = best_speeds["V_endurance"]
    P_E    = best_speeds["P_endurance"]
    V_R    = best_speeds["V_range"]
    P_R    = best_speeds["P_range"]

    slope  = P_R / V_R
    V_tang = np.array([0, V_max])

    plt.plot(V_tang * 3.6, slope * V_tang / 1000, "--", label="Best-range tangent", color="gray", linewidth = 0.5)
    plt.plot(V_kmh, P_par_arr  / 1000, "-",  label=r"Parasite ($P_{par}$)",  color="k")
    plt.plot(V_kmh, P_i_arr    / 1000, ":",  label=r"Induced ($P_i$)",        color="k")
    plt.plot(V_kmh, P_prof_arr / 1000, "-.", label=r"Profile ($P_p+P_d$)",    color="k")
    plt.plot(V_kmh, P_tr_arr   / 1000, "--", label=r"Tail rotor ($P_{tr}$)",  color="k")
    plt.plot(V_kmh, P_tot_arr  / 1000,       label=r"Total ($P_{total}$)",    color="k", linewidth=2)

    plt.plot(V_E * 3.6, P_E / 1000, "k^", markersize=8,
             label=rf"$V_{{E}}$ = {V_E*3.6:.0f} km/h, P = {P_E/1000:.0f} kW")
    plt.plot(V_R * 3.6, P_R / 1000, "ks", markersize=8,
             label=rf"$V_{{R}}$ = {V_R*3.6:.0f} km/h, P = {P_R/1000:.0f} kW")

    plt.xlabel("Airspeed V [km/h]")
    plt.ylabel("Power [kW]")
    plt.legend(fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, V_max * 3.6)
    plt.ylim(0)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved → {save_path}")
    plt.show()


# Main
if __name__ == "__main__":
    vi_h = compute_hover_induced_velocity(W, rho, A)

    P_hov_ideal, P_hov_act, P_hov_bem = compute_hover_powers(W, vi_h, FM, k, sigma, Cd_p, rho, tip_speed, R)

    C_T_hov_mr, C_L_mean_hov_mr, C_T_hov_tr, C_L_mean_hov_tr = compute_mean_lift_coefficients(W, P_hov_bem, rho, A, A_tr, Omega, Omega_tr, 
                                    R, R_tr, sigma, sigma_tr, l_tr)

    V_arr, P_par_arr, P_i_arr, P_prof_arr, P_tr_arr, P_tot_arr = compute_forward_flight(
        V_max, W, rho, A, A_tr, vi_h,
        CDS_f, k, sigma, Cd_p, tip_speed, Omega, R,
        sigma_tr, Cd_p_tr, tip_speed_tr, Omega_tr, R_tr,
        k_tr, fin_blockage, l_tr
    )

    best_speeds = compute_best_speeds(V_arr, P_tot_arr)

    plot_results(V_arr, P_par_arr, P_i_arr, P_prof_arr, P_tr_arr, P_tot_arr,
                 best_speeds, V_max)
    