import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import warnings

# Suppress harmless warnings that occur when fsolve guesses extreme values
warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- 1. System Constants ---
FA0 = 20.0             # Molar feed rate of A [mol/s]
T0_feed = 350.0        # Initial reactor temperature [K]
P0 = 580000.0          # Feed Pressure [Pa]
R = 8.314              # Ideal gas constant [J/(mol K)]
CpA = 45.0             # Heat capacity of reacting mixture [J/(mol K)]
delta_H = -30000.0     # Heat of reaction [J/mol]

# Coolant Constants
Ta_inlet = 303.0       # Coolant temperature entering at V = V_final [K]
Fc = 10.0              # Coolant flow rate [mol/s]
Cpc = 75.0             # Coolant heat capacity [J/(mol K)]
Ua = 100.0             # Calculated Ua [W/(m^3 K)]

# Kinetics
k_ref = 0.2 / 60.0     # Rate constant [s^-1] at 273 K
E_A = 20000.0          # Activation energy [J/mol]
Kc_ref = 30.0          # Equilibrium constant at 400 K

v0 = (FA0 * R * T0_feed) / P0

# --- 2. Helper Functions ---
def calc_k(T):
    return k_ref * np.exp((-E_A / R) * (1/T - 1/273.0))

def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1/T - 1/400.0))

# --- 3. The Coupled ODE System ---
def pfr_odes(V, y):
    X = y[0]
    T = y[1]
    Ta = y[2]
    
    # SAFEGUARDS: Prevent solver from breaking math with impossible guesses
    if T <= 1e-5:
        T = 1e-5
        
    k = calc_k(T)
    Kc = calc_Kc(T)
    
    if Kc <= 1e-10:
        Kc = 1e-10

    v = v0 * (T / T0_feed)
    if v <= 1e-10:
        v = 1e-10
        
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    
    rA = -k * (CA - (CB / Kc))
    
    # 1. Mole Balance
    dXdV = -rA / FA0
    
    # 2. Reactor Energy Balance
    dTdV = (Ua * (Ta - T) + rA * delta_H) / (FA0 * CpA)
    
    # 3. Coolant Energy Balance (COUNTER-CURRENT)
    dTadV = - (Ua * (Ta - T)) / (Fc * Cpc)
    
    return [dXdV, dTdV, dTadV]

# --- 4. THE SHOOTING METHOD ---
V_final = 40.0 

def shooting_objective(Ta0_guess):
    y0 = [0.0, T0_feed, Ta0_guess[0]]
    # Using 'Radau' implicit stiff solver. Removed max_step to let it auto-scale.
    sol = solve_ivp(pfr_odes, [0, V_final], y0, method='Radau')
    Ta_end = sol.y[2, -1] 
    return Ta_end - Ta_inlet 

# Initial guess for coolant temp at V=0
initial_guess = [350.0]
Ta0_correct = fsolve(shooting_objective, initial_guess)[0]

print("="*50)
print(f"Calculated Coolant Temp at V=0: {Ta0_correct:.2f} K")

# --- 5. Final Integration with Correct Start ---
y0_final = [0.0, T0_feed, Ta0_correct]
sol = solve_ivp(pfr_odes, [0, V_final], y0_final, method='Radau', dense_output=True)

V_vals = np.linspace(0, V_final, 500)
X_vals, T_vals, Ta_vals = sol.sol(V_vals)
Xe_vals = [calc_Kc(T) / (1 + calc_Kc(T)) for T in T_vals]

print(f"Final Conversion at V={V_final}: {X_vals[-1]:.4f}")
print("="*50)

# --- 6. Plotting Profiles ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Conversions vs Volume
ax1.plot(V_vals, X_vals, 'b-', linewidth=2, label='Actual Conversion ($X_A$)')
ax1.plot(V_vals, Xe_vals, 'g--', linewidth=2, label='Equilibrium Conversion ($X_{Ae}$)')
ax1.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax1.set_ylabel('Conversion', fontsize=12)
ax1.set_title('Conversion Profiles', fontsize=14)
ax1.set_ylim(0, 1.05)
ax1.grid(True, alpha=0.4)
ax1.legend()

# Plot 2: Temperatures vs Volume
ax2.plot(V_vals, T_vals, 'r-', linewidth=2, label='Reactor Temp ($T$)')
ax2.plot(V_vals, Ta_vals, 'c-', linewidth=2, label='Coolant Temp ($T_a$)')
ax2.plot(V_final, Ta_inlet, 'ko', markersize=8, label=f'Boundary: {Ta_inlet} K')
ax2.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax2.set_ylabel('Temperature (K)', fontsize=12)
ax2.set_title('Temperature Profiles (Counter-Current)', fontsize=14)
ax2.grid(True, alpha=0.4)
ax2.legend()

plt.tight_layout()
plt.show()

# --- 7. Plotting Trajectory ---
plt.figure(figsize=(9, 6))
T_background = np.linspace(300, 800, 200)
Xe_background = [calc_Kc(T) / (1 + calc_Kc(T)) for T in T_background]

plt.plot(T_background, Xe_background, 'g--', linewidth=2, label='Equilibrium Curve ($X_e$)')
plt.plot(T_vals, X_vals, 'b-', linewidth=2, label='Reactor Operating Path ($X$)')
plt.plot(T_vals[0], X_vals[0], 'ko', markersize=8, label='Inlet (Feed)')
plt.plot(T_vals[-1], X_vals[-1], 'ro', markersize=8, label='Outlet')

plt.xlabel('Reactor Temperature ($T$) [K]', fontsize=12)
plt.ylabel('Conversion ($X$)', fontsize=12)
plt.title('Conversion vs. Temperature Trajectory', fontsize=14)
plt.xlim(300, 800)
plt.ylim(0, 1.05)
plt.grid(True, alpha=0.4)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()