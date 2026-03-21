class ApacheAH64:
    """
    Data storage class for the AH-64 Apache helicopter parameters 
    needed for performance, stability, and control calculations.
    """
    name = "AH-64"

    # -----------------------------------------------------------------------
    # General Description & Manufacturer Specs
    # -----------------------------------------------------------------------
    specs = {
        "type": "Military",          # [1]
        "payload": 4160,             # [kg] Maximum payload capacity
        "v_max": 81.39,               # [m/s] Maximum speed
        "v_cruise": 79.17,            # [m/s] Cruise speed
        "range": 481,               # [km] Maximum range
        "endurance": 2.733,           # [hrs] Maximum endurance
        "hover_ceiling_ige": 4570,       # [m] Hover ceiling altitude IGE
        "hover_ceiling_oge": 3505,     # [m] Hover ceiling altitude OGE
        "v_best_range_mfg": None,    # [m/s] Manufacturer speed for best range
        "v_best_endurance_mfg": None # [m/s] Manufacturer speed for best endurance
    }

    # -----------------------------------------------------------------------
    # Mass, Inertia & Fuselage Properties
    # -----------------------------------------------------------------------
    helicopter = {
        # Weights converted from lbs [1]
        "mass_empty": 4657.5,        # [kg] Empty weight (10,268 lb)
        "mass_mtow": 8005.9,         # [kg] Max takeoff mass (17,650 lb)
        "fuel_capacity": 1107.7,     # [kg] Fuel capacity (2,442 lb)
        "W": 78537.8,                # [N] Total helicopter weight at MTOW (mass_mtow * 9.81)
        
        "I_yy": None,                # [kg*m^2] Moment of inertia in y-axis
        "cg_to_hub_z": None,         # [m] Vertical distance from rotor hub center to CG
        "fuselage_cd": None,         # [-] Fuselage drag coefficient
        "flat_plate_area": 1.6,     # [m^2] Equivalent flat plate area
        "has_wings": True,           # [bool] Presence of wings (affects drag calculations)        

        "L_tot": 17.729,              # [m] Total length (58 ft 2 in) [1]
    }

    # -----------------------------------------------------------------------
    # Main Rotor & Blade Parameters
    # -----------------------------------------------------------------------
    main_rotor = {
        "airfoil": "HH-02",          # [1]
        "N_blades": 4,               # [-] Number of blades 
        "R": 7.3152,                  # [m] Main rotor radius 
        "c": 0.5334,                  # [m] Blade chord 
        "rpm": 288.8,                # [rev/min] Derived from tip speed & radius
        "omega": 30.25,              # [rad/s] Rotor rotational velocity (726 ft/sec / 24 ft) [1]
        "tip_speed": 221.28,         # [m/s] Tip speed (726 ft/sec) [1]
        "twist": -9,                 # [deg] Blade twist [1]
        "hinge_offset_ratio": 0.038, # [-] e/R ratio [1]
        "blade_mass": None,          # [kg] Mass of a single blade
        "blade_I": 5152.1,           # [kg*m^2] Polar moment of inertia J (3,800 slug-ft^2) [1, 2]
        "sigma": 0.092,              # [-] Rotor solidity [1]
        "cd_p": None,                # [-] Mean profile drag coefficient
        "k": 1.15,                   # [-] Induced drag power factor
        "FM": None                   # [-] Figure of Merit in ACT theory
    }

    # -----------------------------------------------------------------------
    # Tail Rotor Parameters
    # -----------------------------------------------------------------------
    tail_rotor = {
        "airfoil": "NACA 63-414",    # [1]
        "N_blades": 4,               # [-] Number of blades [1]
        "R_tr": 1.3975,               # [m] Tail rotor radius [1]
        "c_tr": 0.253,               # [m] Tail rotor chord (0.83 ft) [1]
        "omega_tr": 147.17,          # [rad/s] TR rotational velocity (677 ft/sec / 4.6 ft) [1]
        "tip_speed_tr": 206.35,      # [m/s] Tail rotor tip speed (677 ft/sec) [1]
        "twist_tr": -8,              # [deg] Tail rotor twist [1]
        "blade_I_tr": 13.56,         # [kg*m^2] Polar moment of inertia J (10 slug-ft^2) [1, 2]
        "sigma_tr": 0.231,           # [-] Tail rotor solidity [1]
        "l_tr": 9.0163,                # [m] Tail rotor distance from main shaft = L_tot - R - R_tr = 17.729 - 7.3152 - 1.3975
        "cd_p_tr": None,             # [-] Tail rotor profile drag coefficient
        "k_tr": 1.4,                 # [-] Tail rotor induced drag factor (1.3 - 1.5) [5, 6]
        "fin_blockage": 1.1          # [-] Vertical fin blockage factor [5, 6]

    }

    # -----------------------------------------------------------------------
    # Engines (Optional additions from specs)
    # -----------------------------------------------------------------------
    engines = {
        "type": "General Electric T700-GE-701", # [1]
        "number": 2,                            # [1]
        "max_to_rating_hp": 3392,               # [hp] Maximum T.O. rating [1]
        "max_usable_power_hp": 2828             # [hp] Maximum usable power [1]
    }

    # -----------------------------------------------------------------------
    # Induced Velocities & Power Calculations (Placeholders)
    # -----------------------------------------------------------------------
    induced_velocity = {
        "hover": None,               
        "forward_flight": None       
    }

    power = {
        "ideal_hover": None,         
        "act_hover": None,           
        "bem_hover": None,           
        "parasite_drag": None,       
        "induced": None,             
        "profile_drag": None,        
        "tail_rotor": None,          
        "total_forward": None        
    }

    flight_condition = {
        "rho": 1.225,                
        "V": None,                   
        "C": 0.0,                    
        "v_best_range_calc": None,   
        "v_best_endurance_calc": None
    }