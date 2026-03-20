import numpy as np
import matplotlib.pyplot as plt

def calculate_trim():
    # ==========================================
    # 1. HELICOPTER PARAMETERS (Apache A64 Data)
    # ==========================================
    m = 10433.0            # Mass (kg)
    g = 9.81               # Gravity (m/s^2)
    W = m * g              # Weight (N)
    
    rho = 1.225            # Air density (kg/m^3)
    Omega = 30.25          # Rotor speed (rad/s)
    R = 7.315              # Rotor radius (m)
    f_area = 1.7           # Equivalent flat plate area for drag (m^2)
    
    Nb = 4                 # Number of blades
    c = 0.53               # Blade chord (m) - UPDATE to your Part 1 value
    cla = 5.7              # Lift curve slope (1/rad) - UPDATE if different
    
    # Derived parameters
    V_tip = Omega * R               # Tip speed (m/s)
    A = np.pi * R**2                # Rotor disc area (m^2)
    sigma = (Nb * c) / (np.pi * R)  # Rotor solidity

    # ==========================================
    # 2. FLIGHT SPEEDS TO EVALUATE
    # ==========================================
    # From hover (0.1 m/s to avoid division by zero) to max speed (approx 80 m/s or 155 knots)
    V_array = np.linspace(0.1, 80, 50) 
    
    # Lists to store results
    theta_0_list = []
    theta_c_list = []
    theta_f_list = []

    # ==========================================
    # 3. TRIM CALCULATION LOOP
    # ==========================================
    for V in V_array:
        # Calculate Parasite Drag and Fuselage Pitch
        D_W = 0.5 * rho * V**2 * f_area
        theta_f = np.arctan(-D_W / W)       # Fuselage pitch angle (rad), negative means nose down
        T = np.sqrt(W**2 + D_W**2)          # Required Thrust (N)
        
        # Non-dimensional parameters
        C_T = T / (rho * A * V_tip**2)      # Thrust coefficient
        mu = V / V_tip                      # Advance ratio
        
        # Initial guess for induced velocity (using hover inflow)
        lambda_i = np.sqrt(C_T / 2.0)
        
        # Iterative solver for lambda_i, theta_0, and a_1 (theta_c)
        tolerance = 1e-6
        max_iter = 500
        
        for iteration in range(max_iter):
            # Matrix components for the 2x2 linear system [A][x] = [B]
            # [A11 A12] [theta_0] = [B1]
            # [A21 A22] [a_1    ] = [B2]
            
            # Equation 1: Longitudinal Flapping Equation
            A11 = -(8.0 / 3.0) * mu
            A12 = 1.0 - 2.5 * mu**2
            B1  = -2.0 * mu**2 * theta_f - 2.0 * mu * lambda_i
            
            # Equation 2: Thrust Equation (BEM)
            A21 = (2.0 / 3.0) * (1.0 + 1.5 * mu**2)
            A22 = mu
            B2  = (4.0 * C_T) / (cla * sigma) + mu * theta_f + lambda_i
            
            # Solve the linear system
            A_mat = np.array([[A11, A12], [A21, A22]])
            B_vec = np.array([B1, B2])
            
            # Extract Collective (theta_0) and Cyclic (theta_c = a_1)
            try:
                solution = np.linalg.solve(A_mat, B_vec)
                theta_0 = solution[0]
                a_1 = solution[1]
            except np.linalg.LinAlgError:
                print(f"Convergence failed at V = {V} m/s")
                break
            
            # Calculate non-dimensional control plane inflow (lambda_c)
            lambda_c = mu * (theta_f - a_1)
            
            # Calculate NEW induced velocity using Glauert's momentum formula
            denominator = 2.0 * np.sqrt(mu**2 + (lambda_i + lambda_c)**2)
            lambda_i_new = C_T / denominator
            
            # Check for convergence
            if abs(lambda_i_new - lambda_i) < tolerance:
                break
            
            # Relaxation factor to ensure steady convergence
            lambda_i = 0.2 * lambda_i_new + 0.8 * lambda_i 

        # Store results (converted to degrees for plotting)
        theta_0_list.append(np.degrees(theta_0))
        theta_c_list.append(np.degrees(a_1))
        theta_f_list.append(np.degrees(theta_f))

    # ==========================================
    # 4. PLOTTING THE RESULTS
    # ==========================================
    plt.figure(figsize=(10, 6))
    
    plt.plot(V_array, theta_0_list, label=r'Collective Pitch ($\theta_0$)', color='blue', linewidth=2)
    plt.plot(V_array, theta_c_list, label=r'Longitudinal Cyclic ($\theta_c$)', color='red', linewidth=2)
    plt.plot(V_array, theta_f_list, label=r'Fuselage Pitch ($\theta_f$)', color='green', linestyle='--', linewidth=2)
    
    plt.title('Helicopter Trim Controls vs. Forward Velocity (Apache A64)', fontsize=14)
    plt.xlabel('Forward Velocity V (m/s)', fontsize=12)
    plt.ylabel('Pitch Angle (Degrees)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.axhline(0, color='black', linewidth=1)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    calculate_trim()