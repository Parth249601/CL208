import numpy as np
import matplotlib.pyplot as plt

times = np.array([0,100,200,300,400])
c_a = np.array([1000,500,333,250,200])

#? Linearize the data

y = np.log(c_a/c_a[0])

x = times

#? Fit a line to the linearized data

coefficients = np.polyfit(x, y, 1)
slope, intercept = coefficients[0], coefficients[1]

k = slope

#? Plot the original data and the fitted line

plt.scatter(x,y, label='Data Points', color='blue')
plt.plot(x, slope*x + intercept, label='Fitted Line', color='red')
plt.xlabel('Time (s)')
plt.ylabel('ln(c_a/c_a[0])')
plt.title('Linearized Data and Fitted Line')
plt.text(0.98, 0.02, f'Slope (k) = {slope:.6f}', transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.legend()
plt.savefig('q1_linearized_first_order.png', dpi=300, bbox_inches='tight')
plt.show()


#! Not a first order reaction as can be seen in the plot q1_integral_order=1.png

#? Trying the second order fit 

y_second_order = 1/c_a - 1/c_a[0]
coefficients_second_order = np.polyfit(x, y_second_order, 1)
slope_second_order, intercept_second_order = coefficients_second_order[0], coefficients_second_order[1]
k_second_order = slope_second_order

#? plot the data
plt.scatter(x,y_second_order, label='Data Points', color='blue')
plt.plot(x, slope_second_order*x + intercept_second_order, label='Fitted Line', color='red')
plt.xlabel('Time (s)')
plt.ylabel('1/c_a - 1/c_a[0]')
plt.title('Second Order Linearized Data and Fitted Line')
plt.text(0.98, 0.02, f'Slope (k) = {slope_second_order:.6f}', transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.legend()
plt.savefig('q1_linearized_second_order.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Rate constant for first order reaction: {k}")
print(f"Rate constant for second order reaction: {k_second_order}")

