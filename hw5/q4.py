import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- 1. System Constants ---
FA0 = 20.0             # Molar feed rate of A [mol/s]
T0_feed = 350.0        # Initial feed temperature [K]
P0 = 580000.0          # Feed Pressure [Pa]
R = 8.314              # Ideal gas constant [J/(mol K)]
CpA = 45.0             # Heat capacity [J/(mol K)]
delta_H = -30000.0     # Heat of reaction [J/mol]
Ta = 303.0             # Coolant temperature [K]

# Kinetics
k_ref = 0.2 / 60.0     # Rate constant [s^-1] at 273 K
E_A = 20000.0          # Activation energy [J/mol]
Kc_ref = 30.0          # Equilibrium constant at 400 K

# Initial volumetric flow rate
v0 = (FA0 * R * T0_feed) / P0

# =====================================================================
# Ua is the overall heat transfer coefficient * heat exchange area.
Ua = 100          # [W/(m^3 K)] or [J/(s m^3 K)]
# =====================================================================

# --- 2. Helper Functions ---
def calc_k(T):
    return k_ref * np.exp((-E_A / R) * (1/T - 1/273.0))

def calc_Kc(T):
    return Kc_ref * np.exp((-delta_H / R) * (1/T - 1/400.0))

# --- 3. The Coupled ODE System ---
def pfr_odes(V, y):
    X = y[0]
    T = y[1]
    
    # Calculate local kinetics
    k = calc_k(T)
    Kc = calc_Kc(T)
    
    # Calculate local volumetric flow and concentrations
    v = v0 * (T / T0_feed)
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    
    # Rate of formation of A (rA is negative because A is consumed)
    rA = -k * (CA - (CB / Kc))
    
    # 1. Mole Balance ODE
    dXdV = -rA / FA0
    
    # 2. Energy Balance ODE (Using your exact handwritten derivation)
    dTdV = (Ua * (Ta - T) + rA * delta_H) / (FA0 * CpA)
    
    return [dXdV, dTdV]

# --- 4. Event tracking to find V at X = 0.75 ---
# solve_ivp can automatically detect when a specific value is reached
def target_conversion(V, y):
    return y[0] - 0.75
# Set to False so the solver continues past 0.75 to complete the plot
target_conversion.terminal = False 

# --- 5. Run the Integrator ---
V_span = (0, 30) # Integrate from 0 to 30 cubic meters (adjust if needed)
y0 = [0.0, T0_feed] # Initial conditions: X=0, T=350 K

# Solve the system
sol = solve_ivp(pfr_odes, V_span, y0, dense_output=True, events=target_conversion, max_step=0.1)

# --- 6. Process the Results ---
V_vals = np.linspace(0, sol.t[-1], 500)
X_vals, T_vals = sol.sol(V_vals)

# Calculate Equilibrium Conversion (Xe) at every point using the local temperature
Xe_vals = [calc_Kc(T) / (1 + calc_Kc(T)) for T in T_vals]

# Check if the target conversion was reached and extract the Volume
if len(sol.t_events[0]) > 0:
    V_target = sol.t_events[0][0]
    print("="*50)
    print(f"SUCCESS: Target conversion of 75% achieved!")
    print(f"Required PFR Volume: {V_target:.4f} m^3")
    print("="*50)
else:
    V_target = None
    print("="*50)
    print("WARNING: Target conversion of 75% was NOT achieved within the volume span.")
    print("Try increasing V_span or check your Ua value. The cooling might be too strong or too weak.")
    print("="*50)

# --- 7. Plotting ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Conversion vs Volume
ax1.plot(V_vals, X_vals, 'b-', linewidth=2, label='Actual Conversion ($X_A$)')
ax1.plot(V_vals, Xe_vals, 'g--', linewidth=2, label='Equilibrium Conversion ($X_{Ae}$)')
if V_target is not None:
    ax1.axvline(x=V_target, color='r', linestyle=':', label=f'V = {V_target:.2f} m³')
    ax1.plot(V_target, 0.75, 'ro') # Red dot at target

ax1.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax1.set_ylabel('Conversion', fontsize=12)
ax1.set_title('Conversion vs. Volume in Cooled PFR', fontsize=14)
ax1.set_ylim(0, 1.05)
ax1.grid(True, alpha=0.4)
ax1.legend()

# Plot 2: Temperature vs Volume
ax2.plot(V_vals, T_vals, 'r-', linewidth=2, label='Reactor Temp ($T$)')
ax2.axhline(y=Ta, color='c', linestyle='--', label=f'Coolant Temp ($T_a = {Ta}$ K)')
if V_target is not None:
    ax2.axvline(x=V_target, color='k', linestyle=':', label=f'Target Volume')

ax2.set_xlabel('Reactor Volume, V ($m^3$)', fontsize=12)
ax2.set_ylabel('Temperature (K)', fontsize=12)
ax2.set_title('Temperature Profile vs. Volume', fontsize=14)
ax2.grid(True, alpha=0.4)
ax2.legend()

plt.tight_layout()
plt.show()