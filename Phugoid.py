import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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
    h = 1.75               # Mast height above CG (m)
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
    print(f"Derivatives: Xu={Xu:.4f}, Mu={Mu:.4f}, Mq={Mq:.4f}")
    
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
        # Formatting to print real and imaginary parts clearly
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
    
    # Plot the eigenvalues
    real_parts = np.real(eigenvalues)
    imag_parts = np.imag(eigenvalues)
    
    # Calculate plot limits based on the furthest eigenvalue
    limit = max(np.max(np.abs(real_parts)), np.max(np.abs(imag_parts))) * 1.5
    
    # Set the limits of the axes FIRST
    plt.xlim(-limit, limit)
    plt.ylim(-limit, limit)
    
    # Use the 'limit' variable to extend the shading all the way to the edges
    plt.axvspan(0, limit, color='red', alpha=0.1, label='Unstable Region (RHP)')
    plt.axvspan(-limit, 0, color='green', alpha=0.1, label='Stable Region (LHP)')
    
    # Use 'x' markers for poles
    plt.scatter(real_parts, imag_parts, color='red', marker='x', s=100, linewidths=2, label='System Eigenvalues (Poles)')
    
    # Add crosshairs for the origin to divide stable (left) and unstable (right) planes
    plt.axhline(0, color='black', linewidth=1.5)
    plt.axvline(0, color='black', linewidth=1.5)
    
    plt.title('Hover Phugoid Eigenvalue Plot (Pole Map)', fontsize=14)
    plt.xlabel('Real Part [1/s] (Growth/Decay Rate)', fontsize=12)
    plt.ylabel('Imaginary Part [rad/s] (Oscillation Frequency)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='best')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    calculate_hover_phugoid_lecture_method()