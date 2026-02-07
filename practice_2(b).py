import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

#? Constants for simulation
k = 0.1       # Rate constant (L / (mol * hr))
v0 = 10.0     # Initial volumetric flow rate (L/hr)
CA0 = 1.0     # Initial concentration of A (mol/L)
FA0 = v0 * CA0 # Initial molar flow rate of A (mol/hr)
IB0 = 0.5     # Molar injection rate of B (mol / (hr * L))
V_max = 50.0  # Total reactor volume (L)
V_span = (0, V_max)
V_eval = np.linspace(0, V_max, 100)

# --- Part (i): Constant Volumetric Flow Rate ---
def model_const_v(V, y):
    CA, CB = y
    dCAdV = -(k / v0) * CA * CB
    dCBdV = (IB0 / v0) - (k / v0) * CA * CB
    return [dCAdV, dCBdV]

y0_const = [CA0, 0.0]
sol_const = solve_ivp(model_const_v, V_span, y0_const, t_eval=V_eval)

# --- Part (ii): Non-constant Volumetric Flow Rate ---
# Assume B is injected as a solution adding volume (e.g., 0.05 L/hr per L of reactor)
vi = 0.05 # Volumetric injection rate: L_injected / (hr * L_reactor)

def model_var_v(V, y):
    FA, FB, v = y
    # Concentrations are molar flow / volumetric flow
    CA = FA / v
    CB = FB / v
    
    dv_dV = vi
    dFA_dV = -k * CA * CB
    dFB_dV = IB0 - k * CA * CB
    
    return [dFA_dV, dFB_dV, dv_dV]

y0_var = [FA0, 0.0, v0]
sol_var = solve_ivp(model_var_v, V_span, y0_var, t_eval=V_eval)

# Extract concentrations for Part (ii) by dividing molar flow by volume
CA_var = sol_var.y[0] / sol_var.y[2]
CB_var = sol_var.y[1] / sol_var.y[2]

# --- Plotting Results ---
plt.figure(figsize=(12, 6))

# Subplot 1: Concentration of A
plt.subplot(1, 2, 1)
plt.plot(sol_const.t, sol_const.y[0], 'b-', label='Const. $v$ (i)')
plt.plot(sol_var.t, CA_var, 'b--', label='Var. $v$ (ii)')
plt.xlabel('Reactor Volume ($V$) [L]')
plt.ylabel('Concentration of A ($C_A$) [mol/L]')
plt.title('Concentration Profile of A')
plt.legend()
plt.grid(True)

# Subplot 2: Concentration of B
plt.subplot(1, 2, 2)
plt.plot(sol_const.t, sol_const.y[1], 'r-', label='Const. $v$ (i)')
plt.plot(sol_var.t, CB_var, 'r--', label='Var. $v$ (ii)')
plt.xlabel('Reactor Volume ($V$) [L]')
plt.ylabel('Concentration of B ($C_B$) [mol/L]')
plt.title('Concentration Profile of B')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('pfr_simulation_comparison.png')