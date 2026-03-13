import numpy as np
from scipy.optimize import fsolve

# --- 1. System Constants ---
FA0 = 20.0             # Molar feed rate [mol/s]
T0 = 350.0             # Feed temperature [K]
P0 = 580000.0          # Feed pressure [Pa]
R = 8.314              # Ideal gas constant [J/(mol K)]
CpA = 45.0             # Heat capacity [J/(mol K)]
delta_H = -30000.0     # Heat of reaction [J/mol]
V_CSTR = 0.5           # Reactor Volume [m^3]

# --- 2. Heat Exchange Parameters ---
# Convert U from W/(dm^2 K) to W/(m^2 K)
U_dm2 = 0.5            
U = U_dm2 * 100.0      # Now 50.0 W/(m^2 K)
Area = 500.0           # Heat transfer area [m^2]
UA = U * Area          # Total UA = 25000.0 W/K
Ta = 303.0             # Coolant temperature [K]

# --- 3. Kinetics Constants ---
k_ref = 0.2 / 60.0     # Rate constant [s^-1] at 273 K
E_A = 20000.0          # Activation energy [J/mol]
Kc_ref = 30.0          # Equilibrium constant at 400 K

# Initial Volumetric Flow
v0 = (FA0 * R * T0) / P0

# --- 4. Helper Functions ---
def calc_k(T):
    return k_ref * np.exp((-E_A / R) * (1/T - 1/273.0))

def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1/T - 1/400.0))

# --- 5. The Coupled CSTR Equation ---
def cooled_cstr(T):
    # 1. Use the Energy Balance to calculate X as a strict function of T
    numerator = FA0 * CpA * (T - T0) + UA * (T - Ta)
    denominator = FA0 * (-delta_H)
    X = numerator / denominator
    
    # 2. Calculate Local Kinetics
    k = calc_k(T)
    Kc = calc_Kc(T)
    
    # 3. Calculate Local Volumetric Flow & Concentrations 
    v = v0 * (T / T0)
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    
    # 4. Calculate Reaction Rate
    rA = -k * (CA - (CB / Kc))
    
    # 5. CSTR Mole Balance Residual
    # We want (V * -rA) - (FA0 * X) to equal exactly 0
    return (V_CSTR * -rA) - (FA0 * X)

# --- 6. Solve the Equation ---
# We provide an initial guess for the steady-state temperature.
# Since feed is 350 K and coolant is 303 K, a guess around 340 K is safe.
T_guess = 340.0 
T_ss = fsolve(cooled_cstr, T_guess)[0]

# Plug the solved steady-state temperature back into the Energy Balance to get X
X_ss = (FA0 * CpA * (T_ss - T0) + UA * (T_ss - Ta)) / (FA0 * (-delta_H))

# --- 7. Print Results ---
print("-" * 50)
print("Non-Adiabatic CSTR (With Heat Exchange)")
print("-" * 50)
print(f"Reactor Volume             : {V_CSTR} m^3")
print(f"Overall Heat Transfer (UA) : {UA} W/K")
print(f"Steady-State Conversion (X): {X_ss:.4f}")
print(f"Steady-State Temp (T)      : {T_ss:.2f} K")
print("-" * 50)