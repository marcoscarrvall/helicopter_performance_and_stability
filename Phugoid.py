import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint

from data import ApacheAH64

def calculate_hover_phugoid_lecture_method():
    # ==========================================
    # 1. HELICOPTER PARAMETERS (Apache A64 Data)
    # ==========================================
    m = ApacheAH64.helicopter["W"] / 9.81  # Mass (kg)
    g = 9.81               # Gravity (m/s^2)
    R = ApacheAH64.main_rotor["R"]              # Rotor radius (m)
    Omega = ApacheAH64.main_rotor["omega"]          # Rotor speed (rad/s)
    c = ApacheAH64.main_rotor["c"]               # Blade chord (m)
    rho = ApacheAH64.flight_condition["rho"]            # Air density (kg/m^3)
    cla = 2 * np.pi      # Lift curve slope (1/rad)
    
    Iy = ApacheAH64.helicopter["I_yy"]          # Pitch moment of inertia (kg*m^2)
    h = ApacheAH64.helicopter["cg_to_hub_z"]         # Mast height above CG (m)
    Ib = ApacheAH64.main_rotor["blade_I"]             # Single blade flapping inertia (kg*m^2)
    theta_0 = np.radians(11.5) # Hover collective pitch (rad)
    
    # Derived Parameters
    V_tip = Omega * R
    A_disc = np.pi * R**2
    
    # ==========================================
    # 2. AERODYNAMIC & FLAPPING PARAMETERS
    # ==========================================
    gamma = (rho * c * cla * R**4) / Ib
    
    C_T = (m * g) / (rho * A_disc * V_tip**2)
    lambda_i = np.sqrt(C_T / 2.0)
    
    A1_du = (1.0 / V_tip) * ((8.0 / 3.0) * theta_0 - 2.0 * lambda_i)
    A2_dq = -16.0 / (gamma * Omega)
    
    # ==========================================
    # 3. ANALYTICAL STABILITY DERIVATIVES
    # ==========================================
    Xu = -g * A1_du
    Mu = (m * g * h / Iy) * A1_du
    
    Mq = (m * g * h / Iy) * A2_dq
    print(f"Derivatives: Xu={Xu:.4f}, Mu={Mu:.4f}, Mq={Mq:.4f}\n")
    
    # ** SLIDE 13 APPROXIMATION: Neglect Xq **
    Xq = 0.0 
    
    # ==========================================
    # 4. STATE MATRIX & EIGENVALUE ANALYSIS
    # ==========================================
    # State vector: [u, q, theta_f]
    A_matrix = np.array([
        [Xu,  Xq, -g],
        [Mu,  Mq,  0],
        [ 0,   1,  0]
    ])
    
    eigenvalues, eigenvectors = np.linalg.eig(A_matrix)
    
    print("=== ALL SYSTEM EIGENVALUES ===")
    for i, eig in enumerate(eigenvalues):
        print(f"Eigenvalue {i+1}: {np.real(eig):.4f} {'+' if np.imag(eig) >= 0 else '-'} {np.abs(np.imag(eig)):.4f}j")
    print("==============================\n")
    
    # Extract the complex conjugate pair (Phugoid Mode)
    complex_roots = [r for r in eigenvalues if np.iscomplex(r)]
    
    if complex_roots:
        phugoid_root = complex_roots[0]
        real_part = np.real(phugoid_root)
        imag_part = np.imag(phugoid_root)
        
        omega_n = np.abs(phugoid_root)                    
        zeta = -real_part / omega_n                       
        period = 2 * np.pi / np.abs(imag_part)            
        
        print("=== HOVER PHUGOID CHARACTERISTICS ===")
        print(f"Eigenvalue                  : {real_part:.4f} +/- {np.abs(imag_part):.4f}j")
        print(f"Natural Frequency (omega_n) : {omega_n:.4f} rad/s")
        print(f"Damping Ratio (zeta)        : {zeta:.4f} (Negative = UNSTABLE)")
        print(f"Period of Oscillation (T)   : {period:.2f} seconds")

    # ==========================================
    # 5. VISUALIZATION (POLE MAP)
    # ==========================================
    plt.figure(figsize=(8, 6))
    real_parts = np.real(eigenvalues)
    imag_parts = np.imag(eigenvalues)
    
    # Calculate plot limits based on the furthest eigenvalue
    limit = max(np.max(np.abs(real_parts)), np.max(np.abs(imag_parts))) * 1.5
    
    plt.xlim(-limit, limit)
    plt.ylim(-limit, limit)
    plt.axvspan(0, limit, color='red', alpha=0.1, label='Unstable Region (RHP)')
    plt.axvspan(-limit, 0, color='green', alpha=0.1, label='Stable Region (LHP)')
    plt.scatter(real_parts, imag_parts, color='red', marker='x', s=100, linewidths=2, label='System Eigenvalues (Poles)')
    plt.axhline(0, color='black', linewidth=1.5)
    plt.axvline(0, color='black', linewidth=1.5)
    plt.title('Hover Phugoid Eigenvalue Plot (Pole Map)', fontsize=14)
    plt.xlabel('Real Part [1/s] (Growth/Decay Rate)', fontsize=12)
    plt.ylabel('Imaginary Part [rad/s] (Oscillation Frequency)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

    # ==========================================
    # 6. TIME HISTORY SIMULATION (SUPERIMPOSED)
    # ==========================================
    # Define the system of differential equations: dx/dt = A * x
    def phugoid_system(x, t):
        return A_matrix.dot(x)
    
    # Initial conditions: [u = 0.1 m/s, q = 0.0 rad/s, theta_f = 0.0 rad]
    x0 = [0.1, 0.0, 0.0]
    
    # Time vector (Simulate for 60 seconds)
    t = np.linspace(0, 60, 1000)
    
    # Integrate the system
    response = odeint(phugoid_system, x0, t)
    
    # Extract states
    u_resp = response[:, 0]
    q_resp = np.degrees(response[:, 1])       # Convert to deg/s
    theta_resp = np.degrees(response[:, 2])   # Convert to deg
    
    # Plotting the Superimposed Time History
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Primary Y-Axis (Left side for Velocity)
    color1 = 'blue'
    ax1.set_xlabel('Time (seconds)', fontsize=12)
    ax1.set_ylabel('Velocity $u$ (m/s)', color=color1, fontsize=12)
    line1 = ax1.plot(t, u_resp, color=color1, linewidth=2, label='Velocity ($u$)')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, linestyle=':', alpha=0.7)

    # Secondary Y-Axis (Right side for Angles)
    ax2 = ax1.twinx()  
    color2 = 'red'
    color3 = 'green'
    ax2.set_ylabel(r'Angles: Pitch $\theta_f$ (deg) & Rate $q$ (deg/s)', color='black', fontsize=12)  
    line2 = ax2.plot(t, q_resp, color=color2, linewidth=2, linestyle='--', label='Pitch Rate ($q$)')
    line3 = ax2.plot(t, theta_resp, color=color3, linewidth=2, label=r'Pitch Attitude ($\theta_f$)')
    ax2.tick_params(axis='y', labelcolor='black')

    # Combine legends from both axes into one box
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left', fontsize=12)

    plt.title('Superimposed Hover Phugoid Response to 0.1 m/s Gust', fontsize=14)
    fig.tight_layout() 
    plt.show()

if __name__ == "__main__":
    calculate_hover_phugoid_lecture_method()