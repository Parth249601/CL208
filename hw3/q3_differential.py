import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import math


# --- 1. DATA INPUT ---
time = np.array([0, 2, 5, 10, 20, 60, 120, math.inf])
C_R = np.array([10, 13.5, 18.5, 26, 38.5, 69, 85, 91])
C_R_inf = 91  # Concentration at infinity (Equilibrium)

# Initial Concentration of R
C_R0 = C_R[0]

# --- 2. DIFFERENTIAL METHOD CALCULATION ---
# Formula: dC_R/dt = (k1 + k2) * (C_R_inf - C_R)

dC_R_dt = np.diff(C_R) / np.diff(time)  # dC_R/dt using finite difference
# We will use the mid-point time for plotting
time_mid = (time[:-1] + time[1:]) / 2
# Calculate the term (C_R_inf - C_R) at mid-points
C_R_mid = (C_R[:-1] + C_R[1:]) / 2






slope, intercept, r_value, p_value, std_err = linregress((C_R_inf - C_R_mid), dC_R_dt)

# The slope corresponds to (k1 + k2)
k_sum = slope

print(f"Regression Results:")
print(f"Slope (k1 + k2): {slope:.4f} min^-1")
print(f"R-squared: {r_value**2:.4f} (Close to 1.0 indicates a good fit)")

C_total = 100

k1 = k_sum * C_R_inf/C_total
k2 = k_sum - k1

print(f"\nCalculated Rate Constants:")
print(f"Forward Rate Constant (k1): {k1:.4f} min^-1")
print(f"Reverse Rate Constant (k2): {k2:.4f} min^-1")

print(f"\nProposed Rate Equation:")
print(f"-r_A = {k1:.4f} C_A - {k2:.4f} C_R")

plt.plot((C_R_inf - C_R_mid), dC_R_dt, "o", label='dC_R/dt vs (C_R_inf - C_R)', linestyle = "")
plt.plot((C_R_inf - C_R_mid), slope*(C_R_inf - C_R_mid) + intercept, label=f'Best Fit Line (slope={slope:.4f}), intercept={intercept:.4f}', color='red', linewidth=2)

plt.xlabel('(C_R_inf - C_R)')
plt.ylabel('dC_R/dt')
plt.title('Differential Method Plot for Reversible Reaction')
plt.legend()
plt.grid(True)
plt.show()
