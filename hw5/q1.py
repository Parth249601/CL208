import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import quad

# --- 1. Constants ---
P0 = 580000.0      # Feed Pressure [Pa]
T0 = 350.0         # Feed Temperature [K]
FA0 = 20.0         # Molar feed rate [mol/s]
R = 8.314          # Ideal gas constant [J/mol*K]
CpA = 45.0         # Heat capacity [J/mol*K]
delta_H = -30000.0 # Heat of reaction [J/mol]

# Initial Volumetric Flow Rate (v0) in m³/s
v0 = (FA0 * R * T0) / P0
#? Error: Important I didn't convert this to s^-1
k_ref = 0.2 / 60.0 # Converted to s^-1

# --- 2. Find Equilibrium Conversion ---
def find_equilibrium(X):
    # Adiabatic energy balance
    T = T0 - (delta_H * X) / CpA
    # Equilibrium constant at T
    Kc = 30 * np.exp((-delta_H / R) * (1 / T - 1 / 400))
    # Equilibrium conversion equation
    return X - (Kc / (1 + Kc))

X_eq_guess = 0.5
X_eq = fsolve(find_equilibrium, X_eq_guess)[0]
X_target = 0.9 * X_eq

# --- 3. Calculate PFR Volume ---
def integrand(X):
    # 1. Local temperature
    T = T0 - (delta_H * X) / CpA
    # 2. Local volumetric flow (assuming ideal gas, equimolar A -> B)
    v = v0 * (T / T0) 
    # 3. Local concentrations
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    # 4. Local constants
    k_local = k_ref * np.exp((-20000 / R) * (1 / T - 1 / 273))
    Kc_local = 30 * np.exp((-delta_H / R) * (1 / T - 1 / 400))
    # 5. Exact reversible rate
    #? Never use the formula (1-X/Xe) for reversible reactions, always use the full expression with Kc
    rA = -k_local * (CA - (CB / Kc_local))
    # 6. Levenspiel value
    return FA0 / (-rA)

V_PFR, _ = quad(integrand, 0, X_target)

# --- 4. Calculate CSTR Volume ---
# Exit temperature
T_target = T0 - (delta_H * X_target) / CpA

# Volumetric flow rate at exit
v_target = v0 * (T_target / T0)

# Concentrations at exit
CA_exit = (FA0 * (1 - X_target)) / v_target
CB_exit = (FA0 * X_target) / v_target

# Kinetics at exit temperature
k_target = k_ref * np.exp((-20000 / R) * (1 / T_target - 1 / 273))
Kc_target = 30 * np.exp((-delta_H / R) * (1 / T_target - 1 / 400))

# Reaction rate at exit conditions
rA_exit = -k_target * (CA_exit - (CB_exit / Kc_target))

# CSTR Design Equation
V_CSTR = (FA0 * X_target) / (-rA_exit)

# --- 5. Print Results ---
print("-" * 40)
print(f"Adiabatic Equilibrium Conversion (X_e): {X_eq:.4f}")
print(f"Target Conversion (90% of X_e): {X_target:.4f}")
print(f"Exit Temperature: {T_target:.2f} K")
print("-" * 40)
print(f"Required PFR Volume: {V_PFR:.4f} m³")
print(f"Required CSTR Volume: {V_CSTR:.4f} m³")
print("-" * 40)



def levenspiel_y(X):
    T = T0 - (delta_H * X) / CpA
    v = v0 * (T / T0) 
    CA = (FA0 * (1 - X)) / v
    CB = (FA0 * X) / v
    k_local = k_ref * np.exp((-20000 / R) * (1 / T - 1 / 273))
    Kc_local = 30 * np.exp((-delta_H / R) * (1 / T - 1 / 400))
    rA = -k_local * (CA - (CB / Kc_local))
    return FA0 / (-rA)

# --- 3. Generate Data for Plot ---
# We plot slightly past X_target to show the curve sweeping upward, 
# but stop before X_eq to avoid dividing by zero.
X_vals = np.linspace(0, X_target * 1.05, 100)
Y_vals = [levenspiel_y(x) for x in X_vals]

# Find the specific y-value at the target conversion for the CSTR
Y_target = levenspiel_y(X_target)

# --- 4. Plotting ---
plt.figure(figsize=(10, 6))

# Plot the main Levenspiel curve
plt.plot(X_vals, Y_vals, 'k-', linewidth=2, label=r'Levenspiel Curve ($F_{A0}/-r_A$)')

# Shade PFR Area (Area under the curve)
X_pfr = np.linspace(0, X_target, 100)
Y_pfr = [levenspiel_y(x) for x in X_pfr]
plt.fill_between(X_pfr, 0, Y_pfr, color='blue', alpha=0.3, label='PFR Volume (Integral Area)')

# Draw CSTR Rectangle
# It goes from X=0 to X_target, with a constant height of Y_target
plt.gca().add_patch(plt.Rectangle((0, 0), X_target, Y_target, 
                                  edgecolor='red', facecolor='none', 
                                  linewidth=2, linestyle='--', hatch='//',
                                  label='CSTR Volume (Rectangle Area)'))

# Formatting the plot
plt.axvline(x=X_target, color='gray', linestyle=':')
plt.xlim(0, max(X_vals))
plt.ylim(0, Y_target * 1.2) # Give a little headroom above the CSTR rectangle
plt.xlabel('Conversion ($X$)', fontsize=12)
plt.ylabel(r'$\frac{F_{A0}}{-r_A}$ (m³)', fontsize=14)
plt.title('Levenspiel Plot: PFR vs. CSTR Volume', fontsize=14)

# Add annotation at X=0 to show the starting point
Y_at_zero = levenspiel_y(0)
plt.plot(0, Y_at_zero, 'go', markersize=8, label=f'At X=0: {Y_at_zero:.3f} m³')
plt.annotate(f'Start: X=0\nFA0/-rA={Y_at_zero:.4f}', 
             xy=(0, Y_at_zero), xytext=(0.02, Y_at_zero * 1.3),
             fontsize=9, ha='left',
             arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

plt.legend(loc='upper left', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('levenspiel_plot.png', dpi=300)  # Save the plot as a high-resolution image
plt.show()

print(f"\nAt X=0 (initial): FA0/-rA = {Y_at_zero:.4f} m³")
print(f"At X=X_target (exit): FA0/-rA = {Y_target:.4f} m³")
