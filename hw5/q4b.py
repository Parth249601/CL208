import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- 1. System Constants ---
FA0 = 20.0             # Molar feed rate of A [mol/s]
T0_feed = 350.0        # Initial reactor temperature [K]
P0 = 580000.0          # Feed Pressure [Pa]
R = 8.314              # Ideal gas constant [J/(mol K)]
CpA = 45.0             # Heat capacity of reacting mixture [J/(mol K)]
delta_H = -30000.0     # Heat of reaction [J/mol]

# Coolant Constants
Ta0 = 303.0            # Initial coolant temperature [K]
Fc = 10.0              # Coolant flow rate [mol/s]
Cpc = 75.0             # Coolant heat capacity [J/(mol K)]
Ua = 100.0             # Calculated Ua [W/(m^3 K)] or [J/(s m^3 K)]

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

# --- 3. The Coupled ODE System (Now with 3 equations!) ---
def pfr_odes(V, y):
    X = y[0]
    T = y[1]
    Ta = y[2]
    
    k = calc_k(T)
    Kc = calc_Kc(T)
    v = v0 * (T / T0_feed)
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    
    rA = -k * (CA - (CB / Kc))
    
    # 1. Mole Balance
    dXdV = -rA / FA0
    
    # 2. Reactor Energy Balance
    dTdV = (Ua * (Ta - T) + rA * delta_H) / (FA0 * CpA)
    
    # 3. Coolant Energy Balance (Co-current)
    dTadV = (Ua * (T - Ta)) / (Fc * Cpc)
    
    return [dXdV, dTdV, dTadV]

# --- 4. Run the Integrator ---
# We integrate to a very large volume (e.g., 500 m^3) to find the absolute maximum conversion
V_span = (0, 40) 
y0 = [0.0, T0_feed, Ta0] # Initial conditions: X=0, T=350, Ta=303

sol = solve_ivp(pfr_odes, V_span, y0, dense_output=True, max_step=1.0)

# --- 5. Process the Results ---
V_vals = np.linspace(0, sol.t[-1], 1000)
X_vals, T_vals, Ta_vals = sol.sol(V_vals)
Xe_vals = [calc_Kc(T) / (1 + calc_Kc(T)) for T in T_vals]

# The maximum achievable conversion is the final value in the array
X_max = X_vals[-1]
T_final = T_vals[-1]

print("="*50)
print(f"Maximum Achievable Conversion: {X_max:.4f}")
print(f"Final Thermal Equilibrium Temp: {T_final:.2f} K")
print("="*50)

# --- 6. Plotting ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Conversions vs Volume
ax1.plot(V_vals, X_vals, 'b-', linewidth=2, label='Actual Conversion ($X_A$)')
ax1.plot(V_vals, Xe_vals, 'g--', linewidth=2, label='Equilibrium Conversion ($X_{Ae}$)')
ax1.axhline(y=X_max, color='k', linestyle=':', label=f'Max X = {X_max:.3f}')

ax1.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax1.set_ylabel('Conversion', fontsize=12)
ax1.set_title('Conversion Profiles', fontsize=14)
ax1.set_ylim(0, 1.05)
ax1.grid(True, alpha=0.4)
ax1.legend()

# Plot 2: Temperatures vs Volume
ax2.plot(V_vals, T_vals, 'r-', linewidth=2, label='Reactor Temp ($T$)')
ax2.plot(V_vals, Ta_vals, 'c-', linewidth=2, label='Coolant Temp ($T_a$)')

ax2.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax2.set_ylabel('Temperature (K)', fontsize=12)
ax2.set_title('Temperature Profiles', fontsize=14)
ax2.grid(True, alpha=0.4)
ax2.legend()

plt.tight_layout()
plt.show()


# --- 7. Plotting X and X_e vs Temperature ---
plt.figure(figsize=(9, 6))

# Generate a wide, smooth temperature range just for the background Equilibrium curve
T_background = np.linspace(300, 800, 200)
Xe_background = [calc_Kc(T) / (1 + calc_Kc(T)) for T in T_background]

# Plot the theoretical thermodynamic ceiling (Green dashed line)
plt.plot(T_background, Xe_background, 'g--', linewidth=2, label='Equilibrium Curve ($X_e$)')

# Plot the actual path your reactor took (Blue solid line)
plt.plot(T_vals, X_vals, 'b-', linewidth=2, label='Reactor Operating Path ($X$)')

# Mark the starting and ending points for clarity
plt.plot(T_vals[0], X_vals[0], 'ko', markersize=8, label='Inlet (Feed)')
plt.plot(T_vals[-1], X_vals[-1], 'ro', markersize=8, label='Final Thermodynamic Max')

# Formatting
plt.xlabel('Reactor Temperature ($T$) [K]', fontsize=12)
plt.ylabel('Conversion ($X$)', fontsize=12)
plt.title('Conversion vs. Temperature Trajectory', fontsize=14)
plt.xlim(300, 800)
plt.ylim(0, 1.05)
plt.grid(True, alpha=0.4)
plt.legend(fontsize=10)
plt.tight_layout()

plt.show()