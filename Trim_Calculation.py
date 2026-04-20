import numpy as np
import matplotlib.pyplot as plt

from data import ApacheAH64

def calculate_trim_lecture_method():
    # ==========================================
    # 1. HELICOPTER PARAMETERS (Apache A64 Data)
    # ==========================================
    
    W = ApacheAH64.helicopter["W"]
    R = ApacheAH64.main_rotor["R"]
    rho = ApacheAH64.flight_condition["rho"]

    f_area = ApacheAH64.helicopter["flat_plate_area"]
    if f_area is None:
        f_area = 1.7  

    Omega = ApacheAH64.main_rotor["omega"]         # Rotor speed (rad/s)
    Nb = ApacheAH64.main_rotor["N_blades"]         # Number of blades
    c = ApacheAH64.main_rotor["c"]                 # Blade chord (m)
    cla = 2 * np.pi                                # Lift curve slope (1/rad)
    
    # Derived parameters
    V_tip = Omega * R               # Tip speed (m/s)
    A = np.pi * R**2                # Rotor disc area (m^2)
    sigma = (Nb * c) / (np.pi * R)  # Rotor solidity

    # ==========================================
    # 2. FLIGHT SPEEDS TO EVALUATE
    # ==========================================
    V_array = np.linspace(0.1, 82, 50) # 0.1 m/s to 82 m/s (~155 knots)
    
    theta_0_list = []
    theta_c_list = []
    theta_f_neg_list = [] # List to store negative fuselage pitch (-theta_f)

    # ==========================================
    # 3. TRIM CALCULATION LOOP (Lecture Method)
    # ==========================================
    for V in V_array:
        # Step 1: Fuselage Drag and Required Thrust (Page 5)
        D_fus = 0.5 * rho * V**2 * f_area
        T = np.sqrt(W**2 + D_fus**2)               # Eq. 5
        
        # Calculate actual fuselage pitch (theta_f) using the exact formula
        theta_f = np.arctan(-D_fus / W)
        
        # The drag angle (alpha_d) is the negative of the fuselage pitch
        alpha_d = -theta_f
        
        # Non-dimensional parameters
        C_T = T / (rho * A * V_tip**2)
        mu = V / V_tip
        
        # Step 2: Solve for Induced Velocity (lambda_i) (Page 7, Eq. 10)
        lambda_i = np.sqrt(C_T / 2.0) # Initial guess (hover inflow)
        
        for _ in range(500):
            # Glauert Thrust equation rearranged for fixed-point iteration
            V_term_sin = mu * np.sin(alpha_d)
            V_term_cos = mu * np.cos(alpha_d)
            
            denominator = 2.0 * np.sqrt((V_term_sin + lambda_i)**2 + V_term_cos**2)
            lambda_i_new = C_T / denominator
            
            if abs(lambda_i_new - lambda_i) < 1e-7:
                break
            lambda_i = 0.5 * lambda_i_new + 0.5 * lambda_i # Relaxation
            
        # Step 3: Solve for Collective (theta_0) and Cyclic (theta_c) (Pages 8-9)
        # Matrix A terms
        A11 = (2.0 / 3.0) * (1.0 + 1.5 * mu**2)
        A12 = -mu
        A21 = -(8.0 / 3.0) * mu
        A22 = 1.0 + 1.5 * mu**2  
        
        # Vector B terms
        B1 = (4.0 * C_T) / (cla * sigma) + lambda_i + mu * alpha_d
        B2 = -2.0 * mu * (lambda_i + mu * alpha_d)
        
        A_mat = np.array([[A11, A12], [A21, A22]])
        B_vec = np.array([B1, B2])
        
        # Solve the system
        try:
            solution = np.linalg.solve(A_mat, B_vec)
            theta_0 = solution[0]
            theta_c = solution[1]
        except np.linalg.LinAlgError:
            print(f"Failed to solve at V={V}")
            theta_0, theta_c = 0, 0
            
        # Store results in degrees. We plot -theta_f to keep the curve positive!
        theta_0_list.append(np.degrees(theta_0))
        theta_c_list.append(np.degrees(theta_c))
        theta_f_neg_list.append(np.degrees(-theta_f)) 

    # ==========================================
    # 4. PLOTTING THE RESULTS
    # ==========================================
    plt.figure(figsize=(9, 6))
    
    plt.plot(V_array, theta_0_list, label=r'Collective Pitch ($\theta_0$)', color='blue', linewidth=2)
    plt.plot(V_array, theta_c_list, label=r'Longitudinal Cyclic ($\theta_c$)', color='green', linewidth=2)
    plt.plot(V_array, theta_f_neg_list, label=r'Nose-Down Pitch ($-\theta_f$)', color='red', linestyle='--', linewidth=2)
    
    plt.title('Helicopter Pilot Controls & Pitch vs. Forward Velocity', fontsize=14)
    plt.xlabel('Forward Velocity V (m/s)', fontsize=12)
    plt.ylabel('Angle (Degrees)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    calculate_trim_lecture_method()