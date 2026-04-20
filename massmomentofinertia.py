# ---------------------------------------------------------
# HELICOPTER MASS MOMENT OF INERTIA CALCULATOR (I_yy)
# Includes Centroidal Inertias + Parallel Axis Transfer
# ---------------------------------------------------------

# CONSTANTS: Masses [kg] (From table)
import re
import os

MASS_MF = 6751.0   # Main Fuselage
MASS_LT = 211.68   # Lower Tail
MASS_VT = 296.44   # Vertical Tail
MASS_MR = 745.0    # Main Rotor

# ---------------------------------------------------------
# PARAMETERS: Dimensions [meters]
# TO BE CORRECTED: Replace with actual component dimensions
# ---------------------------------------------------------

# Main Fuselage (Modeled as Solid Ellipsoid)
L_mf = 7.5  # Total Length
H_mf = 2.4   # Total Height

# Convert total dimensions to semi-axes (a and c) for the formula
a_mf = L_mf / 2.0
c_mf = H_mf / 2.0

# Lower Tail (Modeled as Horizontal Solid Cylinder)
L_lt = 5   # Length 
r_lt = 0.5   # Radius (Placeholder - replace with actual radius)

# Vertical Tail (Modeled as Vertical Solid Cylinder)
L_vt = 2.1   # Length/Height 
r_vt = 0.5   # Radius (Placeholder - replace with actual radius)

# Main Rotor (Modeled as Solid Disk)
r_mr = 7.315   # Radius of the rotor disk (Placeholder)

# ---------------------------------------------------------
# PARAMETERS: Distances from Component Centroid to CG [meters]
# ---------------------------------------------------------
d_mf_x = -0.25
d_mf_z = -0.4

d_lt_x = 6.8
d_lt_z = 1.2

d_vt_x = 9
d_vt_z = 0.4 

d_mr_x = 0.25  
d_mr_z = 1.6  

# ---------------------------------------------------------
# CALCULATION 1: Centroidal Moments of Inertia (I_bar_y)
# ---------------------------------------------------------
# Ellipsoid: (1/5) * m * (a^2 + c^2)
I_bar_mf = (1/5) * MASS_MF * (a_mf**2 + c_mf**2)

# Horizontal Cylinder: (1/12)*m*L^2 + (1/4)*m*r^2
I_bar_lt = (1/12) * MASS_LT * (L_lt**2) + (1/4) * MASS_LT * (r_lt**2)

# Vertical Cylinder: (1/12)*m*L^2 + (1/4)*m*r^2
I_bar_vt = (1/12) * MASS_VT * (L_vt**2) + (1/4) * MASS_VT * (r_vt**2)

# Solid Disk: (1/4)*m*r^2
I_bar_mr = (1/4) * MASS_MR * (r_mr**2)

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
print("-" * 65)
print("Helicopter Mass Moment of Inertia Calculation (I_yy)")
print("-" * 65)
print(f"{'Component':<18} | {'I_bar (Centroidal)':<20} | {'Total I_yy (About CG)'}")
print("-" * 65)
print(f"{'Main Fuselage:':<18} | {I_bar_mf:>12,.1f} kg*m^2 | {I_yy_mf:>14,.1f} kg*m^2")
print(f"{'Lower Tail:':<18} | {I_bar_lt:>12,.1f} kg*m^2 | {I_yy_lt:>14,.1f} kg*m^2")
print(f"{'Vertical Tail:':<18} | {I_bar_vt:>12,.1f} kg*m^2 | {I_yy_vt:>14,.1f} kg*m^2")
print(f"{'Main Rotor:':<18} | {I_bar_mr:>12,.1f} kg*m^2 | {I_yy_mr:>14,.1f} kg*m^2")
print("-" * 65)
print(f"SYSTEM TOTAL I_yy: {I_yy_total:,.2f} kg*m^2")
print("-" * 65)

print("\nAttempting to update data.py...")
data_file_path = "data.py"

if os.path.exists(data_file_path):
    # Read the current contents of data.py
    with open(data_file_path, "r", encoding="utf-8") as file:
        content = file.read()
    
    # Regex pattern to find "I_yy": <anything>, and capture the surrounding syntax
    # It looks for "I_yy": followed by any spacing, then anything that isn't a comma, then a comma
    pattern = r'("I_yy"\s*:\s*)[^,]+(,)'
    
    # Replacement string: Keeps the first group ("I_yy": ), inserts the rounded total, keeps the second group (,)
    replacement = rf'\g<1>{I_yy_total:.2f}\g<2>'
    
    # Perform the substitution
    new_content, count = re.subn(pattern, replacement, content)
    
    if count > 0:
        # Write the updated content back to the file
        with open(data_file_path, "w", encoding="utf-8") as file:
            file.write(new_content)
        print(f"✅ Successfully updated 'I_yy' in {data_file_path} to {I_yy_total:.2f}!")
    else:
        print("⚠️ Warning: Could not find the 'I_yy' key in data.py to update. Check the formatting.")
else:
    print(f"❌ Error: {data_file_path} not found in the current directory.")