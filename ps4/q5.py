import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data
c_a = np.array([1, 0.84, 0.68, 0.53, 0.38, 0.27, 0.16, 0.09, 0.04, 0.018, 0.006, 0.0025])
times = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
c_e_0 = 0.01  # millimol/L

# --- PART 1: DIFFERENTIAL METHOD 
print("--- Differential Method ---")
dca_dt = np.diff(c_a) / np.diff(times)
avg_c_a = (c_a[:-1] + c_a[1:]) / 2

def rate_law(c_a, k, c_m):
    return k * c_a * c_e_0 / (c_a + c_m)

popt, pcov = curve_fit(rate_law, avg_c_a, -dca_dt)
k_diff, cm_diff = popt
print(f"Estimated k: {k_diff:.4f}, Estimated c_m: {cm_diff:.4f}")

# Plotting Differential
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
c_a_fit = np.linspace(min(c_a), max(c_a), 100)
fitted_rate = rate_law(c_a_fit, k_diff, cm_diff)
plt.scatter(avg_c_a, -dca_dt, label='Data Points')
plt.plot(c_a_fit, fitted_rate, label='Differential Fit', color='red')
plt.xlabel('Conc (mmol/L)')
plt.ylabel('Rate (-r_a)')
plt.title('Differential Analysis')
plt.legend()
plt.grid(True)

# --- PART 2: INTEGRAL METHOD (New Code) ---
print("\n--- Integral Method ---")

# 1. Define the integrated function: t = f(C_A)
# We treat C_A as the independent variable 'x' and Time as dependent 'y' for fitting
def integrated_model(ca_val, k_val, cm_val):
    ca0 = c_a[0]
    # Equation: t = (1/(k*Ce0)) * [ (CA0 - CA) + Cm * ln(CA0/CA) ]
    term1 = (ca0 - ca_val)
    term2 = cm_val * np.log(ca0 / ca_val)
    return (1 / (k_val * c_e_0)) * (term1 + term2)

# 2. Fit the curve
# We use the differential results as initial guesses (p0) to ensure convergence
popt_int, pcov_int = curve_fit(integrated_model, c_a, times, p0=[k_diff, cm_diff])
k_int, cm_int = popt_int

print(f"Estimated k: {k_int:.4f}, Estimated c_m: {cm_int:.4f}")

# 3. Plotting Integral
plt.subplot(1, 2, 2)
# Generate smooth points for the curve
# Note: We generate smooth C_A values, then calculate the corresponding Time
ca_smooth = np.linspace(c_a[0], c_a[-1], 100)
t_predicted = integrated_model(ca_smooth, k_int, cm_int)

plt.scatter(times, c_a, label='Experimental Data', color='blue')
# Plot calculated time vs smooth concentration
plt.plot(t_predicted, ca_smooth, label='Integral Fit', color='green', linestyle='--')

plt.xlabel('Time (min)')
plt.ylabel('Concentration ($C_A$)')
plt.title('Integral Analysis ($C_A$ vs t)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


#? --- Differential Method ---
#? Estimated k: 20.2490, Estimated c_m: 0.2113

#? --- Integral Method ---
#? Estimated k: 19.9572, Estimated c_m: 0.1990

