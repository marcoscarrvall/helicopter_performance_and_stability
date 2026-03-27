import math
from matplotlib import ticker
from data import ApacheAH64 as ap
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from matplotlib.patches import Circle, Wedge, Patch



def find_tilt_coefficients(
    V_ff,
    vi,
    alpha_c,
    theta0,
    q_deg,
    p_deg,
    print_results=True,
    theta1c_deg = 2.0,
    theta1s_deg = 1.0,
):
    """
    Compute rotor disc tilt coefficients (a0, a1, b1) for a helicopter
    in forward flight with body pitch rate q and roll rate p.

    Flapping angle (Eq. 13):
        β(ψ) = a0 - a1·cos(ψ) - b1·sin(ψ)

    Parameters
    ----------
    V_ff      : float  Forward flight speed [m/s]
    vi        : float  Induced velocity [m/s]  (from Glauert / Fig. 3)
    alpha_c   : float  Control plane angle of attack [deg]
    theta0    : float  Collective pitch angle [deg]
    q_deg     : float  Body pitch rate [deg/s]
    p_deg     : float  Body roll rate  [deg/s]
    print_results : bool  Print formatted summary (default True)

    Returns
    -------
    dict with keys:
        'a0'      : float  Cone angle [rad]
        'a1'      : float  Longitudinal tilt [rad]
        'b1'      : float  Lateral tilt [rad]
        'gamma'   : float  Lock number [-]
        'mu'      : float  Advance ratio [-]
        'lambda_i': float  Non-dimensional induced inflow [-]
        'lambda_c': float  Non-dimensional control-plane inflow [-]
    """

    # --- Load rotor/atmosphere parameters from data.py ---
    rho   = ap.flight_condition["rho"]   # [kg/m³]
    R     = ap.main_rotor["R"]           # [m]
    Omega = ap.main_rotor["omega"]       # [rad/s]
    c     = ap.main_rotor["c"]           # [m]
    sigma = ap.main_rotor["sigma"]       # [-]
    I     = ap.main_rotor["blade_I"]     # [kg·m²]
    W     = ap.helicopter["W"]           # [N]

    cl_alpha = 2 * math.pi               # [1/rad]  thin-aerofoil lift slope

    # --- Unit conversions ---
    alpha_c_rad = math.radians(alpha_c)
    theta0_rad  = math.radians(theta0)
    theta1c_rad  = math.radians(theta1c_deg)
    theta1s_rad  = math.radians(theta1s_deg)
    q           = math.radians(q_deg)    # [rad/s]
    p           = math.radians(p_deg)    # [rad/s]

    # --- Non-dimensional parameters ---
    gamma    = (rho * cl_alpha * c * R**4) / I                 # Lock number
    mu       = (V_ff * math.cos(alpha_c_rad)) / (Omega * R)    # advance ratio
    lambda_i = vi / (Omega * R)                                 # induced inflow
    lambda_c = (V_ff * math.sin(alpha_c_rad)) / (Omega * R)    # control-plane inflow
    q_norm   = q / Omega
    p_norm   = p / Omega

    # --- Disc tilt coefficients (Eqs. 22, 25, 26) ---
    # Cone angle — unchanged by roll rate (Eq. 22)
    a0 = gamma * (
        (1.0 / 8.0) * theta0_rad * (1.0 + mu**2)
        + (1.0 / 6.0) * mu * theta1c_rad
        + (1.0 / 6.0) * (lambda_i - lambda_c)
        + (1.0 / 12.0) * mu * p_norm
    )

    # Longitudinal tilt — cross-coupled by both q and p (Eq. 25)
    a1 = (
        - (8.0 / 3.0) * mu * theta0_rad
        - theta1c_rad * (1.0 + 1.5 * mu**2)
        + 2.0 * mu * (lambda_c - lambda_i)
        - p_norm
        - (16.0 / gamma) * q_norm
    ) / (1.0 - 0.5 * mu**2)

    # Lateral tilt — cross-coupled by both q and p (Eq. 26)
    b1 = (
        (4.0 / 3.0) * mu * a0
        - q_norm
        - (16.0 / gamma) * p_norm
    ) / (1.0 + 0.5 * mu**2) - theta1s_rad

    # --- Optional printed summary ---
    if print_results:
        SEP  = "=" * 62
        SEP2 = "-" * 62
        print(SEP)
        print("  AE4314-21  |  Disc Tilt Coefficients  |  Part 2 / App. B")
        print(SEP)
        print("\n[ PARAMETERS FROM data.py ]")
        print(f"  Rotor radius           R      = {R:.4f} m")
        print(f"  Rotor speed            Ω      = {Omega:.4f} rad/s")
        print(f"  Blade chord            c      = {c:.4f} m")
        print(f"  Rotor solidity         σ      = {sigma:.4f}")
        print(f"  Blade inertia          I      = {I:.2f} kg·m²")
        print(f"  Air density            ρ      = {rho:.4f} kg/m³")
        print(f"  Helicopter weight      W      = {W:.1f} N")
        print("\n[ FLIGHT CONDITION ]")
        print(f"  Forward speed          V      = {V_ff:.2f} m/s")
        print(f"  Pitch rate             q      = {q_deg:.1f} deg/s  ({q:.5f} rad/s)")
        print(f"  Roll rate              p      = {p_deg:.1f} deg/s  ({p:.5f} rad/s)")
        print(f"  Collective pitch       θ₀     = {theta0:.2f} deg  ({theta0_rad:.5f} rad)")
        print(f"  Control plane AoA      αc     = {alpha_c:.2f} deg  ({alpha_c_rad:.5f} rad)")
        print(f"  Induced velocity       vi     = {vi:.3f} m/s")
        print("\n[ DERIVED NON-DIMENSIONAL PARAMETERS ]")
        print(f"  Lock number            γ      = {gamma:.4f}")
        print(f"  Advance ratio          μ      = {mu:.4f}")
        print(f"  Induced inflow         λi     = {lambda_i:.4f}")
        print(f"  Control-plane inflow   λc     = {lambda_c:.4f}")
        print(f"  Norm. pitch rate       q/Ω    = {q_norm:.5f}")
        print(f"  Norm. roll rate        p/Ω    = {p_norm:.5f}")
        print("\n[ DISC TILT COEFFICIENTS (Eqs. 22, 25, 26) ]")
        print(SEP2)
        print(f"  Cone angle             a0     = {a0:.6f} rad  =  {math.degrees(a0):.4f} deg")
        print(f"  Longitudinal tilt      a1     = {a1:.6f} rad  =  {math.degrees(a1):.4f} deg")
        print(f"  Lateral tilt           b1     = {b1:.6f} rad  =  {math.degrees(b1):.4f} deg")
        print(SEP2)
        print(f"\n{SEP}")
        print("  Done.")
        print(SEP)

    return {"a0": a0, "a1": a1, "b1": b1,
            "gamma": gamma, "mu": mu,
            "lambda_i": lambda_i, "lambda_c": lambda_c}

def plot_flapping_angle(
    V_ff,
    vi,
    alpha_c,
    q_deg,
    p_deg,
    theta0_deg   = 6.0,   # [deg] Collective pitch
    theta1c_deg  = 2.0,   # [deg] Longitudinal cyclic
    theta1s_deg  = 1.0,   # [deg] Lateral cyclic
    psi_points   = 360,
):
    """
    Compute disc tilt coefficients and plot the blade flapping angle β(ψ)
    for one complete rotor revolution (0° → 360°).
 
    Produces two separate figures:
        1. Cartesian plot  — β vs ψ
        2. Polar plot      — β(ψ) with ψ=0 at bottom, ψ=180 at top
 
    Parameters
    ----------
    V_ff         : float  Forward flight speed [m/s]
    vi           : float  Induced velocity [m/s]
    alpha_c      : float  Control plane angle of attack [deg]
    q_deg        : float  Body pitch rate [deg/s]
    p_deg        : float  Body roll rate  [deg/s]
    theta0_deg   : float  Collective pitch [deg]          (default 6°)
    theta1c_deg  : float  Longitudinal cyclic pitch [deg] (default 2°)
    theta1s_deg  : float  Lateral cyclic pitch [deg]      (default 1°)
    psi_points   : int    Number of azimuth points        (default 360)
 
    Returns
    -------
    fig1 : matplotlib Figure  (Cartesian)
    fig2 : matplotlib Figure  (Polar)
    """
 
    # ------------------------------------------------------------------
    # 1. Get disc tilt coefficients
    # ------------------------------------------------------------------
    coeffs = find_tilt_coefficients(
        V_ff      = V_ff,
        vi        = vi,
        alpha_c   = alpha_c,
        theta0    = theta0_deg,
        q_deg     = q_deg,
        p_deg     = p_deg,
        print_results=True,
    )
 
    a0 = coeffs["a0"]   # [rad]
    a1 = coeffs["a1"]   # [rad]
    b1 = coeffs["b1"]   # [rad]
    a0_deg = math.degrees(a0)
 
    # ------------------------------------------------------------------
    # 2. Build azimuth array and compute β(ψ)
    # ------------------------------------------------------------------
    psi_rad     = np.linspace(0, 2 * math.pi, psi_points, endpoint=False)
    psi_deg_arr = np.degrees(psi_rad)
 
    # β(ψ) = a0 - a1·cos(ψ) - b1·sin(ψ)
    beta_rad    = a0 - a1 * np.cos(psi_rad) - b1 * np.sin(psi_rad)
    beta_deg_arr = np.degrees(beta_rad)
 
    # Key reference points
    idx_min  = int(np.argmin(beta_deg_arr))
    idx_max  = int(np.argmax(beta_deg_arr))
    beta_min = beta_deg_arr[idx_min]
    beta_max = beta_deg_arr[idx_max]
    psi_min  = psi_deg_arr[idx_min]
    psi_max  = psi_deg_arr[idx_max]
 
    # Shared title suffix
    condition_str = (
        f"V = {V_ff} m/s,  q = {q_deg}°/s,  p = {p_deg}°/s  |  "
        f"θ₀ = {theta0_deg}°,  θ₁c = {theta1c_deg}°,  θ₁s = {theta1s_deg}°"
    )
 
    # ------------------------------------------------------------------
    # 3a. Cartesian plot
    # ------------------------------------------------------------------
    fig1, ax = plt.subplots(figsize=(10, 5))
    fig1.patch.set_facecolor("white")
    ax.set_facecolor("white")
 
    # Shaded advancing / retreating halves
    ax.axvspan(0,   180, alpha=0.07, color="#1a6faf")
    ax.axvspan(180, 360, alpha=0.07, color="#d94f3d")
 
    # Region labels (placed just above the curve maximum)
    ax.text(90,  beta_max + 0.08, "Advancing side", ha="center",
            fontsize=9, color="#1a6faf")
    ax.text(270, beta_max + 0.08, "Retreating side", ha="center",
            fontsize=9, color="#d94f3d")
 
    # Main flapping curve
    ax.plot(psi_deg_arr, beta_deg_arr,
            color="#1a6faf", linewidth=2.0, zorder=3, label="β(ψ)")
 
    # Min / max markers
    ax.plot(psi_min, beta_min, "v", color="#d94f3d", markersize=8, zorder=4,
            label=f"min β = {beta_min:.3f}° at ψ = {psi_min:.0f}°")
    ax.plot(psi_max, beta_max, "^", color="#2ca02c", markersize=8, zorder=4,
            label=f"max β = {beta_max:.3f}° at ψ = {psi_max:.0f}°")
 
    # Cardinal azimuth lines
    for psi_card in [0, 90, 180, 270, 360]:
        ax.axvline(psi_card, color="grey", linewidth=0.6,
                   linestyle=":", alpha=0.6, zorder=1)
 
    ax.set_xlim(0, 360)
    ax.set_xlabel("Azimuth angle ψ [deg]", fontsize=11)
    ax.set_ylabel("Flapping angle β [deg]", fontsize=11)

    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xticklabels(
        ["0°\n(Tail)", "90°\n(Adv.)", "180°\n(Nose)", "270°\n(Ret.)", "360°"],
        fontsize=9,
    )
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f°"))
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_edgecolor("lightgrey")
    ax.grid(True, color="lightgrey", linewidth=0.6, linestyle="-")
    ax.legend(loc="lower right", fontsize=9)
 
    fig1.tight_layout()
    fig1.savefig("figures/flapping_angle_cartesian.png", dpi=150, bbox_inches="tight")
    print("  Plot saved → flapping_angle_cartesian.png")
 
    # ------------------------------------------------------------------
    # 3b. Polar plot
    #     ψ = 0   at bottom  (Tail)
    #     ψ = 90  at left    (Retreating)
    #     ψ = 180 at top     (Nose)
    #     ψ = 270 at right   (Advancing)
    # ------------------------------------------------------------------
    fig2, ax_p = plt.subplots(figsize=(7, 7), subplot_kw={"projection": "polar"})
    fig2.patch.set_facecolor("white")
    ax_p.set_facecolor("white")
 
    ax_p.set_theta_zero_location("S")   # ψ = 0° at bottom
    ax_p.set_theta_direction(1)         # counter-clockwise: 0→90→180→270 = bottom→right→top→left
 
    # Fill & curve
    ax_p.fill(psi_rad, beta_deg_arr, alpha=0.15, color="#1a6faf")
    ax_p.plot(psi_rad, beta_deg_arr, color="#1a6faf", linewidth=2.0)
 
    # Min / max markers
    ax_p.plot(math.radians(psi_min), beta_min,
              "v", color="#d94f3d", markersize=9, zorder=5,
              label=f"min β = {beta_min:.3f}°")
    ax_p.plot(math.radians(psi_max), beta_max,
              "^", color="#2ca02c", markersize=9, zorder=5,
              label=f"max β = {beta_max:.3f}°")
 
    # 0→bottom, 90→right (advancing), 180→top, 270→left (retreating)
    ax_p.set_thetagrids(
        [0, 90, 180, 270],
        labels=["0°  (Tail)", "90°  (Adv.)", "180°  (Nose)", "270°  (Ret.)"],
        fontsize=9,
    )
    ax_p.tick_params(labelsize=8)
    ax_p.grid(True, color="lightgrey", linewidth=0.6)
    ax_p.spines["polar"].set_edgecolor("lightgrey")
    ax_p.legend(loc="upper right", bbox_to_anchor=(1.3, 1.12), fontsize=9)
 
    fig2.tight_layout()
    fig2.savefig("figures/flapping_angle_polar.png", dpi=150, bbox_inches="tight")
    print("  Plot saved → flapping_angle_polar.png")
 
    plt.show()
    return fig1, fig2


def plot_aoa_disc(
    V_ff,
    vi,
    alpha_c,
    q_deg,
    p_deg,
    theta0_deg  = 6.0,
    theta1c_deg = 2.0,
    theta1s_deg = 1.0,
    n_r         = 200,
    n_psi       = 720,
):
    """
    Plot curves of constant angle of attack across the rotor disc
    (Leishman-style contour map).
 
    α(r̄, ψ) = θ(ψ) − UP/UT   [small-angle, but with reversed-flow masking]
 
    The reversed-flow region (UT < 0) is hatched separately.
    Disc orientation (top-down view):
        Right  → Advancing side (ψ = 90°)
        Top    → Nose           (ψ = 180°)
        Left   → Retreating     (ψ = 270°)
        Bottom → Tail           (ψ = 0°)
    """
 
    # ------------------------------------------------------------------
    # 1. Parameters
    # ------------------------------------------------------------------
    R     = ap.main_rotor["R"]
    Omega = ap.main_rotor["omega"]
 
    alpha_c_rad = math.radians(alpha_c)
    mu          = (V_ff * math.cos(alpha_c_rad)) / (Omega * R)
    lambda_i    = vi / (Omega * R)
    lambda_c    = (V_ff * math.sin(alpha_c_rad)) / (Omega * R)
    lam         = lambda_i + lambda_c
 
    theta0  = math.radians(theta0_deg)
    theta1c = math.radians(theta1c_deg)
    theta1s = math.radians(theta1s_deg)
 
    coeffs = find_tilt_coefficients(
        V_ff=V_ff, vi=vi, alpha_c=alpha_c,
        theta0=theta0_deg, q_deg=q_deg, p_deg=p_deg,
        print_results=False,
    )
    a0 = coeffs["a0"]
    a1 = coeffs["a1"]
    b1 = coeffs["b1"]
 
    # ------------------------------------------------------------------
    # 2. Polar grid
    # ------------------------------------------------------------------
    r_bar   = np.linspace(0.04, 1.0,  n_r)
    psi_arr = np.linspace(0, 2 * np.pi, n_psi, endpoint=False)
    PSI, R_BAR = np.meshgrid(psi_arr, r_bar)   # (n_r, n_psi)
 
    THETA    = theta0  + theta1c * np.cos(PSI) + theta1s * np.sin(PSI)
    BETA     = a0 - a1 * np.cos(PSI) - b1 * np.sin(PSI)
    BETA_DOT = a1 * np.sin(PSI) - b1 * np.cos(PSI)   # β̇/Ω
 
    UT = R_BAR + mu * np.sin(PSI)          # tangential velocity (non-dim)
    UP = lam + BETA_DOT * R_BAR + BETA * mu * np.cos(PSI)   # perpendicular
 
    # Reversed-flow mask: where UT <= 0 the standard α formula is undefined
    REV_FLOW = UT <= 0.0
 
    # Compute α only in forward-flow region; set NaN in reversed flow
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(~REV_FLOW, np.arctan2(UP, UT), np.nan)
 
    ALPHA_DEG = np.where(~REV_FLOW, np.degrees(THETA - phi), np.nan)
 
    # ------------------------------------------------------------------
    # 3. Cartesian coordinates (disc top-down view)
    #    x = r̄·sin(ψ)  → right = advancing (ψ=90°)
    #    y = r̄·cos(ψ)  → up    = nose      (ψ=180°)
    # ------------------------------------------------------------------
    X = R_BAR * np.sin(PSI)
    Y = R_BAR * np.cos(PSI)
 
    # ------------------------------------------------------------------
    # 4. Interpolate onto regular Cartesian grid
    # ------------------------------------------------------------------
    valid = ~np.isnan(ALPHA_DEG.ravel())
    xy_pts     = np.column_stack([X.ravel()[valid], Y.ravel()[valid]])
    alpha_pts  = ALPHA_DEG.ravel()[valid]
 
    grid_res = 400
    xi = np.linspace(-1, 1, grid_res)
    yi = np.linspace(-1, 1, grid_res)
    XI, YI = np.meshgrid(xi, yi)
 
    ALPHA_GRID = griddata(xy_pts, alpha_pts, (XI, YI), method="linear")
 
    # Mask outside rotor disc
    outside = XI**2 + YI**2 > 1.0
    ALPHA_GRID[outside] = np.nan
 
    # Also mark the reversed-flow zone on the Cartesian grid
    # Reversed flow: r̄·sin(ψ) < −μ  →  x < −μ  (left of x = −μ circle)
    # More precisely: UT = r̄ + μ·sin(ψ) < 0  → the region is a circle of
    # radius μ centred at (−μ, 0) in the (x,y) disc plane.
    REV_GRID = (XI + mu)**2 + YI**2 < mu**2
 
    # ------------------------------------------------------------------
    # 5. Contour levels — use actual data range, step of 2°
    # ------------------------------------------------------------------
    vmin = np.nanpercentile(ALPHA_GRID, 2)    # robust min (ignore tiny corners)
    vmax = np.nanpercentile(ALPHA_GRID, 98)
    step = 2
    levels = np.arange(
        math.floor(vmin / step) * step,
        math.ceil(vmax  / step) * step + step,
        step,
    )
 
    # ------------------------------------------------------------------
    # 6. Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_aspect("equal")

    # Hub
    hub = Circle((0, 0), 0.05, facecolor="white", edgecolor="white",
                 linewidth=1.0, zorder=8)
    ax.add_patch(hub)
 
    # Coloured contour lines (no fill)
    import matplotlib.colors as mcolors
    norm = mcolors.Normalize(vmin=levels[0], vmax=levels[-1])
    cs = ax.contour(XI, YI, ALPHA_GRID, levels=levels,
                    cmap="jet", linewidths=1.5, norm=norm)
    ax.clabel(cs, fmt="%g", fontsize=8, inline=True, inline_spacing=3)
 
    # ScalarMappable for the colorbar
    cf = plt.cm.ScalarMappable(norm=norm, cmap="jet")
    cf.set_array([])
 
    # Reversed-flow region — hatch over filled contour
    ax.contourf(XI, YI, REV_GRID.astype(float), levels=[0.5, 1.5],
                colors=["white"], alpha=0.6, zorder=3)
    ax.contourf(XI, YI, REV_GRID.astype(float), levels=[0.5, 1.5],
                hatches=["////"], colors="none", zorder=4)
    ax.contour(XI, YI, REV_GRID.astype(float), levels=[0.5],
               colors="k", linewidths=1.0, linestyles="--", zorder=5)
    rev_patch = Patch(facecolor="white", edgecolor="black",
                      hatch="////", label="Reversed flow region")
    ax.legend(handles=[rev_patch], loc="lower right", fontsize=8,
              framealpha=0.9, edgecolor="lightgrey")
 
 
    # Rotor disc boundary
    disc = Circle((0, 0), 1.0, fill=False, edgecolor="black",
                  linewidth=1.8, zorder=7)
    ax.add_patch(disc)
 
    # Colourbar
    cbar = fig.colorbar(cf, ax=ax, fraction=0.035, pad=0.04)
    cbar.set_label("Angle of attack α [deg]", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
 
    # Cardinal labels
    off = 1.12
    ax.text( 0,   off, "ψ=180°",      ha="center", va="bottom", fontsize=8.5, color="#333333")
    ax.text( 0,  -off, "ψ=0°",        ha="center", va="top",    fontsize=8.5, color="#333333")
    ax.text( off,  0,  "ψ=90°", ha="left",   va="center", fontsize=8.5, color="#333333")
    ax.text(-off,  0,  "ψ=270°", ha="right", va="center", fontsize=8.5, color="#333333")
    ax.text(-off+0.3,  -off,  "Retreating Side", ha="right", va="center", fontsize=8.5, color="#333333")
    ax.text(off-0.3,  -off,  "Advancing Side", ha="left", va="center", fontsize=8.5, color="#333333")
 
    ax.set_xlim(-1.38, 1.38)
    ax.set_ylim(-1.38, 1.38)
    ax.set_xlabel("Blade x-position (nondimensional)", fontsize=10)
    ax.set_ylabel("Blade y-position (nondimensional)", fontsize=10)
 
    ticks = [-1.0, -0.5, 0.0, 0.5, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.tick_params(labelsize=8)
 
    fig.tight_layout()
    fig.savefig("aoa_disc.png", dpi=150, bbox_inches="tight")
    print("  Plot saved → aoa_disc.png")
    plt.show()
    return fig
 
 
if __name__ == "__main__":
    """plot_flapping_angle(
        V_ff        = 20.0,   # [m/s]
        vi          = 8.72,   # [m/s]
        alpha_c     = 3.0,    # [deg]
        q_deg       = 20.0,   # [deg/s]
        p_deg       = 10.0,   # [deg/s]
        theta0_deg  = 6.0,    # [deg] collective
        theta1c_deg = 2.0,    # [deg] longitudinal cyclic
        theta1s_deg = 1.0,    # [deg] lateral cyclic
    )"""
    
    plot_aoa_disc(
        V_ff        = 20.0,
        vi          = 8.72,
        alpha_c     = 3.0,
        q_deg       = 20.0,
        p_deg       = 10.0,
        theta0_deg  = 6.0,
        theta1c_deg = 2.0,
        theta1s_deg = 1.0,
    )