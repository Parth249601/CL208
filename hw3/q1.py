import numpy as np
import matplotlib.pyplot as plt

#? Given data points
v_0 = np.array([10,3,1.2,0.5]) #L/hr
c_a = np.array([85.7,66.7,50,33.4]) #mol/L
V = 1 #L
c_a_0 = 100 # Initial concentration

X = 2*(c_a_0 - c_a)/(2*c_a_0 - c_a)

rate = -(v_0*c_a_0 - v_0*(1-(X/2))*c_a)/V

abs_rate = np.abs(rate)

slope, intercept = np.polyfit(np.log(c_a), np.log(abs_rate), 1)

plt.figure(figsize=(8,6))
plt.plot(np.log(c_a), slope*np.log(c_a) + intercept, label=f'Best Fit Line (slope={slope:.2f}), intercept={intercept:.2f}', color='red', linewidth=2)
plt.plot(np.log(c_a), np.log(abs_rate), marker='o', label='Original Data Points', linestyle='', markersize=8)
plt.xlabel('Log Concentration of A (mol/L)')
plt.ylabel('Log Rate (mol/L/hr)')
plt.title('Rate vs Concentration of A')
plt.legend()
plt.grid(True)
plt.show()