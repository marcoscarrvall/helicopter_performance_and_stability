import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from data import ApacheAH64

W = ApacheAH64.helicopter["W"]
R = ApacheAH64.main_rotor["R"]
rho = ApacheAH64.flight_condition["rho"]
V_cruise = ApacheAH64.specs["v_cruise"]
f = ApacheAH64.helicopter["flat_plate_area"]

V_max = max(V_cruise * 1.2, 90.0) 

# ==========================================
# 2. ALGORITHM SETUP
# ==========================================
A = np.pi * R**2
vih = np.sqrt(W / (2 * rho * A))

V_array = np.linspace(0, V_max, 200)
vi_array = np.zeros_like(V_array)

def glauert_eq(vi_bar, V_bar, alpha_d):
    # Solves: v_i_bar = 1 / sqrt((V_bar*cos(alpha_d))^2 + (V_bar*sin(alpha_d) + v_i_bar)^2)
    denominator = np.sqrt((V_bar * np.cos(alpha_d))**2 + (V_bar * np.sin(alpha_d) + vi_bar)**2)
    return vi_bar - (1.0 / denominator)

# ==========================================
# 3. NUMERICAL SOLUTION LOOP
# ==========================================
for i, V in enumerate(V_array):
    # Trim condition to find disc angle of attack (alpha_d)
    D_fus = 0.5 * rho * V**2 * f 
    sin_alpha_d = D_fus / W 
    if sin_alpha_d > 1.0: sin_alpha_d = 1.0
    alpha_d = np.arcsin(sin_alpha_d)
    
    # Non-dimensionalize forward velocity
    V_bar = V / vih
    initial_guess = 1.0 if i == 0 else (vi_array[i-1] / vih)
    
    # Solve numerically
    vi_bar_solution = fsolve(glauert_eq, initial_guess, args=(V_bar, alpha_d))[0]
    vi_array[i] = vi_bar_solution * vih

# ==========================================
# 4. CALCULATE CRUISE SPEED POINT EXACTLY
# ==========================================
D_fus_cruise = 0.5 * rho * V_cruise**2 * f
sin_alpha_d_cruise = min(D_fus_cruise / W, 1.0)
alpha_d_cruise = np.arcsin(sin_alpha_d_cruise)

V_bar_cruise = V_cruise / vih
vi_bar_cruise = fsolve(glauert_eq, 0.1, args=(V_bar_cruise, alpha_d_cruise))[0]
vi_cruise_ms = vi_bar_cruise * vih

# ==========================================
# 5. PLOTTING THE RESULTS
# ==========================================
plt.figure(figsize=(10, 6))
plt.plot(V_array, vi_array, label='Induced Velocity ($v_i$)', color='#1f77b4', linewidth=2.5)
plt.axhline(vih, color='gray', linestyle='--', label=f'Hover $v_i$ = {vih:.2f} m/s')

# Plot the specific cruise speed point
plt.plot(V_cruise, vi_cruise_ms, 'ro', markersize=8, label=f'Cruise Speed (285 km/h)')
plt.annotate(f'  {V_cruise:.1f} m/s\n  $v_i$ = {vi_cruise_ms:.2f} m/s', 
             (V_cruise, vi_cruise_ms), textcoords="offset points", xytext=(10, 5), 
             ha='left', va='bottom', fontsize=11, color='red', weight='bold')

plt.title('Variation of Induced Velocity in Forward Flight', fontsize=14)
plt.xlabel('Forward Velocity $V$ (m/s)', fontsize=12)
plt.ylabel('Induced Velocity $v_i$ (m/s)', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=11)
plt.xlim(0, V_max)
plt.ylim(0, vih * 1.2)
plt.tight_layout()

plt.show()