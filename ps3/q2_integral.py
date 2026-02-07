import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#? Defining constants

CA_0 = 2 #mol/dm^3

space_times = np.array([15,38,100,300,1200])
CA = np.array([1.5,1.25,1.0,0.75,0.5])

rate_of_reaction = (CA_0 - CA)/space_times

#! Trying First order fit

# plt.plot(CA, rate_of_reaction, "-o")
# plt.xlabel("CA")
# plt.ylabel("Rate of reaction")
# plt.title("Rate of Reaction vs Concentrations")
# plt.show()

#! Trying Second order fit
# plt.plot(CA**2, rate_of_reaction, "-o")
# plt.xlabel("CA^2")
# plt.ylabel("Rate of reaction")
# plt.title("Rate of Reaction vs Concentrations")
# plt.show()

#! Trying Third order fit

plt.plot(np.power(CA, 3), rate_of_reaction, "-o")
plt.xlabel("CA^3")
plt.ylabel("Rate of reaction")
plt.title("Rate of Reaction vs Concentrations")
plt.show()

coefficients = np.polyfit(np.power(CA,3), rate_of_reaction, 1)
slope, intercept = coefficients[0], coefficients[1]
print(f"The slope of RoR vs CA^3 is {slope} and intercept is {intercept}")
