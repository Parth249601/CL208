import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

# 1. Define the system of ODEs
def reversible_reaction(t,y,k1,k2,k1_prime,k2_prime):
    Ca, Cb, Cc = y
    dCadt = -k1*Ca + k1_prime*Cb
    dCbdt = k1*Ca - (k1_prime + k2)*Cb + k2_prime*Cc
    dCcdt = k2*Cb - k2_prime*Cc
    return [dCadt, dCbdt, dCcdt]

#? Constants of the raction

# 2. Parameters (Rate Constants)

k1, k1_p = 0.1, 0.01
k2, k2_p = 0.003, 0.003

# 3. Initial Conditions [Ca0, Cb0, Cc0]
# Example: Starting with 1.0 units of A and nothing of B or C
y0 = [1.0, 0.0, 0.0]

t_span = (0, 1000)  # Time span for the simulation
t_eval = np.linspace(0,1000,10000)    #Time points for evaluation

sol = solve_ivp(
    reversible_reaction,
    t_span,
    y0,
    args=(k1, k2, k1_p, k2_p),
    t_eval=t_eval
)





plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='$C_A$ (Reactant)')
plt.plot(sol.t, sol.y[1], label='$C_B$ (Intermediate)')
plt.plot(sol.t, sol.y[2], label='$C_C$ (Product)')

plt.xlabel('Time (t)')
plt.ylabel('Concentration')
plt.title(r'Concentration Profiles for Reversible Reactions $A \leftrightarrow B \leftrightarrow C$')
plt.legend()
plt.grid(True)
plt.show()


