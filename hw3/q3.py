import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- 1. DATA INPUT ---
time = np.array([0, 2, 5, 10, 20, 60, 120])
C_R = np.array([10, 13.5, 18.5, 26, 38.5, 69, 85])
C_R_inf = 91  # Concentration at infinity (Equilibrium)

# Initial Concentration of R
C_R0 = C_R[0]

# --- 2. INTEGRAL METHOD CALCULATION ---
# Formula: ln((C_Re - C_R0) / (C_Re - C_R)) = (k1 + k2) * t
# We calculate the term inside the log first
numerator = C_R_inf - C_R0
denominator = C_R_inf - C_R

# Avoid division by zero if C_R reaches C_R_inf exactly (though typically it approaches asymptotically)
# Filter out the infinite point or any point where C_R >= C_R_inf for the log calculation
valid_indices = denominator > 0
time_fit = time[valid_indices]
y_values = np.log(numerator / denominator[valid_indices])

# --- 3. LINEAR REGRESSION ---
# Fit y = mx + c
slope, intercept, r_value, p_value, std_err = linregress(time_fit, y_values)

# The slope corresponds to (k1 + k2)
k_sum = slope

print(f"Regression Results:")
print(f"Slope (k1 + k2): {slope:.4f} min^-1")
print(f"R-squared: {r_value**2:.4f} (Close to 1.0 indicates a good fit)")

# --- 4. CALCULATE INDIVIDUAL CONSTANTS (k1, k2) ---
# We assume the total concentration is C_total = C_A + C_R.
# Based on problem: Feed is 90% n-pentane, 10% i-pentane.
# If C_R0 = 10 represents 10%, then C_total = 100 mol/L.
C_total = 100 
C_A_inf = C_total - C_R_inf

# Equilibrium Constant K_eq = k1 / k2 = [R]e / [A]e
K_eq = C_R_inf / C_A_inf

# System of equations:
# 1) k1 + k2 = slope
# 2) k1 / k2 = K_eq  -> k1 = K_eq * k2

# Substitute (2) into (1):
# K_eq * k2 + k2 = slope -> k2 * (1 + K_eq) = slope
k2 = k_sum / (1 + K_eq)
k1 = K_eq * k2

print(f"\nCalculated Rate Constants:")
print(f"Equilibrium Constant (Kc): {K_eq:.4f}")
print(f"Forward Rate Constant (k1): {k1:.4f} min^-1")
print(f"Reverse Rate Constant (k2): {k2:.4f} min^-1")

print(f"\nProposed Rate Equation:")
print(f"-r_A = {k1:.4f} C_A - {k2:.4f} C_R")

# --- 5. PLOT ---
plt.figure(figsize=(8, 5))
plt.scatter(time_fit, y_values, label='Experimental Data', color='blue')
plt.plot(time_fit, slope * time_fit + intercept, label=f'Fit: y={slope:.3f}x + {intercept:.3f}', color='red')
plt.xlabel('Time (min)')
plt.ylabel('ln((C_Re - C_R0) / (C_Re - C_R))')
plt.title('Integral Method Plot for Reversible Reaction')
plt.legend()
plt.grid(True)
plt.show()