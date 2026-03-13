import numpy as np
from scipy.integrate import quad

# --- 1. Constants ---
P0 = 580000.0          # Feed Pressure [Pa]
FA0 = 20.0             # Molar feed rate of A [mol/s]
R = 8.314              # Ideal gas constant [J/mol*K]
CpA = 45.0             # Heat capacity [J/mol*K]
delta_H = -30000.0     # Heat of reaction [J/mol]

# Kinetic constants
k_ref = 0.2 / 60.0     # Reference rate constant [s^-1] at 273 K
EA = 20000.0           # Activation energy [J/mol]
Kc_ref = 30.0          # Reference equilibrium constant at 400 K

# --- 2. Values from your handwritten derivation ---
T0_feed = 67.536       # Your calculated T0 [K]
v0 = (FA0 * R * T0_feed) / P0  # v0 evaluated at your new T0

# --- 3. Temperature and Kinetics Functions ---
def calc_T(X):
    # Your energy balance: T = T0 - (deltaH/CpA)*X
    return T0_feed - (delta_H / CpA) * X

def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1 / T - 1 / 400.0))

def calc_k(T):
    return k_ref * np.exp((-EA / R) * (1 / T - 1 / 273.0))

# --- 4. Your Exact Integrand ---
def integrand(X):
    T = calc_T(X)
    k = calc_k(T)
    Kc = calc_Kc(T)
    
    # The denominator you derived: (k * T0) / (T * v0) * [1 - X(1 + 1/Kc)]
    denominator = (k * T0_feed) / (T * v0) * (1 - X * (1 + (1 / Kc)))
    
    return 1.0 / denominator

# --- 5. Evaluate the Integral ---
# We integrate up to 0.6499 to avoid dividing by zero at absolute equilibrium
X_target = 0.65

Volume, error = quad(integrand, 0, X_target)

print(f"Calculated v0: {v0:.5f} m³/s")
print(f"Required PFR Volume (up to X={X_target}): {Volume:.4f} m³")