import numpy as np
import matplotlib.pyplot as plt

#? Defining Constants

CA_0 = 1000 # mol/m^3

times = np.array([0,100,200,300,400])
CA = np.array([1000,500,333,250,200])

#! trying to check if its first order
# plt.plot(times, np.log(CA/CA_0), "-o")
# plt.xlabel("Time")
# plt.ylabel("ln(CA/CA_0)")
# plt.show()

#? not a first order reaction

#! trying to check if its second order

plt.plot(times, (1/CA - 1/CA_0), "-o")
plt.xlabel("Time")
plt.ylabel("1/CA - 1/CA_0")
plt.show()

coeffs = np.polyfit(times, 1/CA - 1/CA_0, 1)
slope = coeffs[0] #~10^-5
print(f"slope is {slope}")


