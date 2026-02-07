import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
#? Constants given in the problem

P_initial_cold = 1.0
T_cold = 298.15
T_0 = 373.15  #Intial temperature in K
R = 0.0821    # Ideal gas constant in L·atm/(K·mol)
P_0 = P_initial_cold * (T_0 / T_cold) # roughly 1.2515


C_A_0 = P_0/(R * T_0)  # Initial concentration of A in mol/L


#! Reaction 2A -> B
times = np.array([1,2,3,4,5,6,7,8,9,10,15,20])

pressure_total = np.array([1.14,1.04,0.982,0.94,0.905,0.87,0.85,0.832,0.815,0.8,0.754,0.728])  # in atm

pressure_A = 2*pressure_total - P_0  # Partial pressure of A assuming ideal behavior
concentration_A = pressure_A / (R * T_0)  # Convert pressure to concentration using ideal gas law


#? Assuming the reaction is first order

x_values = np.log(C_A_0 / concentration_A)
y_values = times

# Calculate Pearson correlation
corr, _ = pearsonr(x_values, y_values)
print(f'Pearson Correlation: {corr}')

#Lets try to fit a linear regression model
x_values_reshaped = x_values.reshape(-1, 1)
model = LinearRegression()
model.fit(x_values_reshaped, y_values)
r_squared = model.score(x_values_reshaped, y_values)
print(f"R-squared: {r_squared}")

# plt.xlabel('ln(C_A0 / C_A)')
# plt.ylabel('Time (hr)')
# plt.plot(x_values, y_values, 'o', label='Data Points')
# plt.show()


#? Assuming the reaction is second order 

y_values_2 = ((1 / concentration_A) - (1 / C_A_0))
x_values_2 = times

plt.ylabel('1/C_A - 1/C_A0')
plt.xlabel('Time (min)')
plt.plot(x_values_2, y_values_2, 'o', label='Data Points')
plt.show()