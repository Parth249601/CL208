import matplotlib.pyplot as plt
import numpy as np


times = np.array([1,2,3,4,5,6,7,8,9,10,15,20])

pressures = np.array([1.14,1.04,0.982,0.94,0.905,0.87,0.85,0.832,0.815,0.8,0.754,0.728])
p_0 = 1.25
extent_of_reactions = 2*(1- pressures/p_0)

# log_extent_of_reactions = np.log(extent_of_reactions)
# exp_extent_of_reactions = np.exp(extent_of_reactions)


# plt.scatter(times,exp_extent_of_reactions)
# plt.plot(times,extent_of_reactions)

#inverse_extent_of_reactions = (1/(1- extent_of_reactions)) - 1


inverse_extent_of_reactions = (1/extent_of_reactions) - 1

inverse_times = 1/times

plt.plot(inverse_times,inverse_extent_of_reactions)
# plt.xlabel("Time")
# plt.ylabel("Extent of reaction")





plt.title("Inverse Function vs Time")

# Calculate slope and intercept
# slope, intercept = np.polyfit(times, inverse_extent_of_reactions, 1)
# print(f"Slope: {slope}")
# print(f"Intercept: {intercept}")

plt.show()  