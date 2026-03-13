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

# Kinetics
k_ref = 0.2 / 60.0     # Rate constant [s^-1] at 273 K
E_A = 20000.0          # Activation energy [J/mol]
Kc_ref = 30.0          # Equilibrium constant at 400 K

# Initial Volumetric Flow
v0 = (FA0 * R * T0) / P0

# --- 2. Helper Functions ---
def calc_k(T):
    return k_ref * np.exp((-E_A / R) * (1/T - 1/273.0))

def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1/T - 1/400.0))

# --- 3. The CSTR Equation to Solve ---
def adiabatic_cstr(X):
    # 1. Adiabatic Energy Balance: T is strictly a function of X
    T = T0 + (-delta_H / CpA) * X
    
    # 2. Local Kinetics
    k = calc_k(T)
    Kc = calc_Kc(T)
    
    # 3. Local Volumetric Flow & Concentrations (Gas phase expands with T)
    v = v0 * (T / T0)
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    
    # 4. Reaction Rate
    rA = -k * (CA - (CB / Kc))
    
    # 5. CSTR Mole Balance
    # Standard form: V = (FA0 * X) / (-rA)
    # We rewrite it as V*(-rA) - FA0*X = 0 to prevent division by zero near equilibrium
    return (V_CSTR * -rA) - (FA0 * X)

# --- 4. Solve the Equation ---
# We provide an initial guess for the solver (e.g., X = 0.4)
X_guess = 0.4
X_ss = fsolve(adiabatic_cstr, X_guess)[0]

# Calculate the final steady-state temperature
T_ss = T0 + (-delta_H / CpA) * X_ss

# --- 5. Print Results ---
print("-" * 50)
print("Adiabatic CSTR (Without Heat Exchange)")
print("-" * 50)
print(f"Reactor Volume             : {V_CSTR} m^3")
print(f"Steady-State Conversion (X): {X_ss:.4f}")
print(f"Steady-State Temp (T)      : {T_ss:.2f} K")
print("-" * 50)