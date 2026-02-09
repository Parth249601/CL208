import numpy as np
import matplotlib.pyplot as plt

times = np.array([0,5,9,15,22,30,40,60])
c_a = np.array([2,1.6,1.35,1.1,0.87,0.70,0.53,0.35])

#?  Using Differential method to find the order of the reaction
#? dc_a/dt = -k * c_a^n
#? We can find the order by plotting ln(dc_a/dt) vs ln(c_a)


#? Calculate dc_a/dt using finite difference method
dc_a_dt = np.diff(c_a) / np.diff(times)
c_a_mid = (c_a[:-1] + c_a[1:]) / 2  # Midpoint concentrations for better accuracy

#? Plot ln(dc_a/dt) vs ln(c_a)
coefficients = np.polyfit(np.log(c_a_mid), np.log(-dc_a_dt), 1)
slope, intercept = coefficients[0], coefficients[1]
k = np.exp(intercept)
n = slope
plt.scatter(np.log(c_a_mid), np.log(-dc_a_dt), label='Data Points', color='blue')
plt.plot(np.log(c_a_mid), slope*np.log(c_a_mid) + intercept, label='Fitted Line', color='red')
plt.xlabel('ln(c_a)')
plt.ylabel('ln(-dc_a/dt)')
plt.title('Differential Method: ln(-dc_a/dt) vs ln(c_a)')
plt.legend()
plt.text(0.98, 0.02, f'Slope (k) = {n:.6f}', transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.savefig('q2_differential_method.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"Rate constant (k) = {k}")
print(f"Order of the reaction (n) = {n}")


