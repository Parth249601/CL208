import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# --- 1. System Constants ---
P0 = 580000.0          # Feed Pressure [Pa]
T0_feed = 350.0        # Initial feed & interstage cooling temperature [K]
FA0 = 20.0             # Molar feed rate of A [mol/s]
R = 8.314              # Ideal gas constant [J/mol*K]
CpA = 45.0             # Heat capacity [J/mol*K]
delta_H = -30000.0     # Heat of reaction [J/mol]

# Kinetic constants
k_ref = 0.2 / 60.0     # Reference rate constant [s^-1] at 273 K
EA = 20000.0           # Activation energy [J/mol]
Kc_ref = 30.0          # Reference equilibrium constant at 400 K

# Initial concentration and volumetric flow
v0 = (FA0 * R * T0_feed) / P0   # Initial volumetric flow rate [m^3/s]
CA0 = FA0 / v0                  # Initial concentration [mol/m^3]

X_target_overall = 0.65

# --- 2. Helper Functions ---
def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1 / T - 1 / 400.0))

def calc_k(T):
    return k_ref * np.exp((-EA / R) * (1 / T - 1 / 273.0))

# --- 3. Stage-by-Stage Calculation ---
current_X = 0.0
stage_count = 0
total_volume = 0.0

print(f"Targeting Overall Conversion > {X_target_overall*100}%\n")
print("=" * 60)
counter = 0
while current_X < X_target_overall:
    stage_count += 1
    X_in = current_X 
    
    # A. Find Equilibrium Conversion for this stage
    def stage_equilibrium(X):
        # Adiabatic energy balance for this specific stage
        T_adiabatic = T0_feed - (delta_H / CpA) * (X - X_in)
        Kc = calc_Kc(T_adiabatic)
        return X - (Kc / (1 + Kc))
    
    # Solve for equilibrium (guess slightly above current conversion)
    X_eq = fsolve(stage_equilibrium, X_in + 0.1)[0]
    
    # B. Apply design rule: 90% of equilibrium conversion
    X_out = 0.9 * X_eq
    
    # C. Calculate exit temperature and concentration
    T_out = T0_feed - (delta_H / CpA) * (X_out - X_in)
    v_out = v0 * (T_out / T0_feed)
    CA_out = (FA0 * (1 - X_out)) / v_out
    
    # D. Calculate PFR Volume for this stage
    def integrand(X):
        # 1. Local T from adiabatic energy balance
        T = T0_feed - (delta_H / CpA) * (X - X_in)
        # 2. Local v
        v = v0 * (T / T0_feed)
        # 3. Local concentrations
        CA = (FA0 * (1 - X)) / v
        CB = (FA0 * X) / v
        # 4. Local kinetics
        k_local = calc_k(T)
        Kc_local = calc_Kc(T)
        # 5. Rate equation
        rA = -k_local * (CA - (CB / Kc_local))
        return FA0 / (-rA)
    
    if(counter == 2):
        V_stage, _ = quad(integrand, X_in, X_target_overall)
        total_volume += V_stage
    # Integrate from X_in to X_out
    else:
        V_stage, _ = quad(integrand, X_in, X_out)
        total_volume += V_stage
    
    # E. Print Results for the Stage
    print(f"STAGE {stage_count}")
    print(f"  Inlet Conversion  : {X_in:.4f}")
    print(f"  Exit Conversion   : {X_out:.4f} (Eq: {X_eq:.4f})")
    print(f"  Exit Temperature  : {T_out:.2f} K")
    print(f"  Exit Conc. (CA)   : {CA_out:.2f} mol/m^3")
    print(f"  Stage Volume      : {V_stage:.4f} m^3")
    print("-" * 60)
    
    counter += 1
    # Update for the next stage
    current_X = X_out

print(f"TOTAL SYSTEM VOLUME : {total_volume:.4f} m^3")
print("=" * 60)