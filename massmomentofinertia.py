# ---------------------------------------------------------
# HELICOPTER MASS MOMENT OF INERTIA CALCULATOR (I_yy)
# Includes Centroidal Inertias + Parallel Axis Transfer
# ---------------------------------------------------------

# CONSTANTS: Masses [kg] (From table)
MASS_MF = 6751.0   # Main Fuselage
MASS_LT = 211.68   # Lower Tail
MASS_VT = 296.44   # Vertical Tail
MASS_MR = 745.0    # Main Rotor

# ---------------------------------------------------------
# PARAMETERS: Dimensions [meters]
# TO BE CORRECTED: Replace with actual component dimensions
# L = Length (along x-axis), H = Height (along z-axis)
# ---------------------------------------------------------

# Main Fuselage (Modeled as Ellipsoid)
L_mf = 14.0  
H_mf = 3.5   

# Lower Tail (Modeled as Rectangular Prism)
L_lt = 6.0   
H_lt = 1.0   

# Vertical Tail (Modeled as Rectangular Prism)
L_vt = 1.5   
H_vt = 3.0   

# Main Rotor (Modeled as Rectangular Prism per instructions)
L_mr = 2.0   
H_mr = 0.5   

# ---------------------------------------------------------
# PARAMETERS: Distances from Component Centroid to CG [meters]
# TO BE CORRECTED: Replace with actual moment arms
# ---------------------------------------------------------
d_mf_x = -1.2 
d_mf_z = -0.6 

d_lt_x = 10.5 
d_lt_z = 0.4  

d_vt_x = 9.8  
d_vt_z = 2.9  

d_mr_x = 0.5  
d_mr_z = 3.2  

# ---------------------------------------------------------
# CALCULATION 1: Centroidal Moments of Inertia (I_bar_y)
# ---------------------------------------------------------
# Ellipsoid: (1/20) * m * (L^2 + H^2)
I_bar_mf = (1/20) * MASS_MF * (L_mf**2 + H_mf**2)

# Rectangular Prisms: (1/12) * m * (L^2 + H^2)
I_bar_lt = (1/12) * MASS_LT * (L_lt**2 + H_lt**2)
I_bar_vt = (1/12) * MASS_VT * (L_vt**2 + H_vt**2)
I_bar_mr = (1/12) * MASS_MR * (L_mr**2 + H_mr**2)

# ---------------------------------------------------------
# CALCULATION 2: Transfer Terms (m * d^2)
# d^2 = dx^2 + dz^2
# ---------------------------------------------------------
transfer_mf = MASS_MF * (d_mf_x**2 + d_mf_z**2)
transfer_lt = MASS_LT * (d_lt_x**2 + d_lt_z**2)
transfer_vt = MASS_VT * (d_vt_x**2 + d_vt_z**2)
transfer_mr = MASS_MR * (d_mr_x**2 + d_mr_z**2)

# ---------------------------------------------------------
# CALCULATION 3: Total I_yy (Parallel Axis Theorem)
# I_yy = I_bar + m*d^2
# ---------------------------------------------------------
I_yy_mf = I_bar_mf + transfer_mf
I_yy_lt = I_bar_lt + transfer_lt
I_yy_vt = I_bar_vt + transfer_vt
I_yy_mr = I_bar_mr + transfer_mr

I_yy_total = I_yy_mf + I_yy_lt + I_yy_vt + I_yy_mr

# ---------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------
print("-" * 55)
print("Helicopter Mass Moment of Inertia Calculation (I_yy)")
print("-" * 55)
print(f"{'Component':<18} | {'I_bar (Centroidal)':<18} | {'Total I_yy (About CG)'}")
print("-" * 55)
print(f"{'Main Fuselage:':<18} | {I_bar_mf:>10,.1f} kg*m^2 | {I_yy_mf:>14,.1f} kg*m^2")
print(f"{'Lower Tail:':<18} | {I_bar_lt:>10,.1f} kg*m^2 | {I_yy_lt:>14,.1f} kg*m^2")
print(f"{'Vertical Tail:':<18} | {I_bar_vt:>10,.1f} kg*m^2 | {I_yy_vt:>14,.1f} kg*m^2")
print(f"{'Main Rotor:':<18} | {I_bar_mr:>10,.1f} kg*m^2 | {I_yy_mr:>14,.1f} kg*m^2")
print("-" * 55)
print(f"SYSTEM TOTAL I_yy: {I_yy_total:,.2f} kg*m^2")
print("-" * 55)