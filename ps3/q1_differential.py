# import numpy as np
# import matplotlib.pyplot as plt


# #? Constants given in the problem

# P_initial_cold = 1.0
# T_cold = 298.15
# T_0 = 373.15  #Intial temperature in K
# R = 0.0821    # Ideal gas constant in L·atm/(K·mol)
# P_0 = P_initial_cold * (T_0 / T_cold) # roughly 1.2515



# C_A_0 = P_0/(R * T_0)  # Initial concentration of A in mol/L


# #! Reaction 2A -> B
# times = np.array([1,2,3,4,5,6,7,8,9,10,15,20])

# pressure_total = np.array([1.14,1.04,0.982,0.94,0.905,0.87,0.85,0.832,0.815,0.8,0.754,0.728])  # in atm

# pressure_A = 2*pressure_total - P_0  # Partial pressure of A assuming ideal behavior

# concentration_A = pressure_A / (R * T_0)  # Convert pressure to concentration using ideal gas law

# avg_concentrations = np.zeros(len(concentration_A) - 1)
# for i in range(len(concentration_A)-1):
#     avg_concentrations[i] = (concentration_A[i] + concentration_A[i+1]) / 2


# dCA_dt = np.diff(concentration_A) / np.diff(times)
# # dCA_dt = np.zeros(len(concentration_A) - 1)
# # for i in range(len(concentration_A) - 1):
# #     dCA_dt[i] = (concentration_A[i+1] - concentration_A[i]) / (times[i+1] - times[i])

# magnitude_dCA_dt = np.abs(dCA_dt)

# ln_avg_concentrations = np.log(avg_concentrations)
# ln_magnitude_dCA_dt = np.log(magnitude_dCA_dt)

# plt.plot(ln_avg_concentrations, ln_magnitude_dCA_dt, 'o-')
# plt.xlabel('ln(Average Concentration of A (mol/L))')
# plt.ylabel('ln(Magnitude of dC_A/dt (mol/L·hr))')
# plt.title('ln(dC_A/dt) vs ln(Average C_A)')
# plt.grid(True)
# plt.show()
# #to show slope of the plot
# coefficients = np.polyfit(ln_avg_concentrations, ln_magnitude_dCA_dt, 1)
# slope, intercept = coefficients[0], coefficients[1]

# print(f"Slope of the plot (Reaction Order): {slope}, Intercept: {intercept}")


import numpy as np
import matplotlib.pyplot as plt

# --- 1. Constants ---
T_cold = 25 + 273.15   # 298.15 K
T_rxn = 100 + 273.15   # 373.15 K
R = 0.08206            # L·atm/(K·mol)

# Calculate P_0 at reaction temperature (Gay-Lussac's Law)
P_0 = 1.0 * (T_rxn / T_cold) 
print(f"Calculated P_A0: {P_0:.4f} atm")

# --- 2. Data ---
times = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20])
pressure_total = np.array([1.14, 1.04, 0.982, 0.94, 0.905, 0.87, 0.85, 0.832, 0.815, 0.8, 0.754, 0.728])

# --- 3. Stoichiometry & Conversion ---
# p_A = 2*P_total - P_A0
pressure_A = 2 * pressure_total - P_0
concentration_A = pressure_A / (R * T_rxn)

# --- 4. Differential Analysis (Interval Method) ---
# Calculate -dCA/dt for each interval
dCA = -np.diff(concentration_A) # Negative because concentration decreases
dt = np.diff(times)
rate = dCA / dt

# Calculate Average Concentration for each interval to match the rate
C_avg = (concentration_A[:-1] + concentration_A[1:]) / 2

# --- 5. Log-Log Transformation ---
ln_rate = np.log(rate)
ln_C = np.log(C_avg)

# --- 6. Plotting & Fitting ---
coefficients = np.polyfit(ln_C, ln_rate, 1)
order = coefficients[0]
ln_k = coefficients[1]
k = np.exp(ln_k)

print(f"Order of Reaction (n): {order:.4f}")
print(f"Rate Constant (k): {k:.4f}")

plt.plot(ln_C, ln_rate, 'bo', label='Experimental Data')
plt.plot(ln_C, np.polyval(coefficients, ln_C), 'r-', label=f'Fit: n={order:.2f}')
plt.xlabel('ln(Concentration $C_A$)')
plt.ylabel('ln(Rate $-r_A$)')
plt.title('Differential Method Analysis')
plt.legend()
plt.grid(True)
plt.show()