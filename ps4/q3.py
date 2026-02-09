import numpy as np
import matplotlib.pyplot as plt

positions = np.array([0,5,10,15,20,25,30])
X_values = np.array([0,1.93/100,3.82/100, 5.68/100, 7.48/100, 9.25/100, 11.0/100])


#? Solving the question using differential method

c_a_0 = 2.33 * 10**(-6) # mol/cm^3
v_0 = 19.6 #cm^3/s
dia = 0.158 #cm
area = np.pi * (dia/2)**2 # cm^2

#? Calculating dX/dV values using finite difference method
dX_dV = np.diff(X_values) /(area* np.diff(positions))

#? Midpoint conversions for better accuracy
X_mid = (X_values[:-1] + X_values[1:]) / 2

#? Plotting ln(dX/dV) vs ln(1-X) to find the order of the reaction
coefficients = np.polyfit(np.log(1-X_mid), np.log(dX_dV), 1)
slope, intercept = coefficients[0], coefficients[1]

n = slope
k = np.exp(intercept) * (v_0 / c_a_0**(n-1))

plt.scatter(np.log(1-X_mid), np.log(dX_dV), label='Data Points', color='blue')
plt.plot(np.log(1-X_mid), slope*np.log(1-X_mid) + intercept, label='Fitted Line', color='red')
plt.xlabel('ln(1-X)')
plt.ylabel('ln(dX/dV)')
plt.title('Differential Method: ln(dX/dV) vs ln(1-X)')
plt.legend()
plt.text(0.98, 0.02, f'Slope (n) = {n:.6f}\nIntercept = {intercept:.6f}', transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.savefig('q3_differential_method.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"Order of the reaction (n) = {n}")
print(f"Rate constant (k) = {k} units depend on the order of the reaction")



#? Solving using regression analysis (integral method) to find the order of the reaction


x_values = area * positions
y_values = -np.log(1-X_values)

coefficients_integral = np.polyfit(x_values, y_values, 1)
slope_integral, intercept_integral = coefficients_integral[0], coefficients_integral[1]

k_integral = v_0 * slope_integral  

plt.scatter(x_values, y_values, label='Data Points', color='blue')
plt.plot(x_values, slope_integral*x_values + intercept_integral, label='Fitted Line', color='red')
plt.xlabel('V')
plt.ylabel('-ln(1-X)')
plt.title('Integral Method: -ln(1-X) vs V')
plt.legend()
plt.text(0.98, 0.02, f'Slope (n) = {slope_integral:.6f}\nIntercept = {intercept_integral:.6f}', transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.savefig('q3_integral_method.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"Slope of the line = {slope_integral}")
print(f"Rate constant (k) = {k_integral} 1/s")

#! Fits the data perfectly 
#!Slope of the line = 0.19801710680603773
#!Rate constant (k) = 3.88113529339834 1/s