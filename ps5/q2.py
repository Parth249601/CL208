import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def reaction_system(t,y, k1, k2,v_0, C_b_0, V):
    C_a , C_b, C_c, C_d = y
    dCadt = -k1 * C_a*C_b - k2 * C_a* (C_b**2)
    dCbdt = v_0*C_b_0/V -k1 * C_a*C_b - 2*k2 * C_a* (C_b**2)
    dCcdt = k1 * C_a*C_b
    dCddt = k2 * C_a* (C_b**2)

    return [dCadt, dCbdt, dCcdt, dCddt]

#? Constants of the raction

k1 = 0.1
k2 = 0.002
V = 100
y_0 = [10.0, 20.0, 0, 0]

t_span = (0, 40)  # Time span for the simulation
t_eval = np.linspace(0,40,1000)    #Time points for evaluation

v_0 = 1.0

sol = solve_ivp(
    reaction_system,
    t_span,
    y_0,
    args=(k1, k2,v_0, y_0[1], V),
    t_eval=t_eval
)

selectivity_c = k1/(k2*sol.y[1])
selectivity_d = 1/selectivity_c

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left subplot: Concentration profiles
ax1.plot(sol.t, sol.y[0], label='$C_A$ (Reactant)', linewidth=2)
ax1.plot(sol.t, sol.y[1], label='$C_B$ (Reactant)', linewidth=2)
ax1.plot(sol.t, sol.y[2], label='$C_C$ (Product)', linewidth=2)
ax1.plot(sol.t, sol.y[3], label='$C_D$ (Product)', linewidth=2)
ax1.set_xlabel('Time', fontsize=11)
ax1.set_ylabel('Concentration (mol/L)', fontsize=11)
ax1.set_title('Concentration Profiles', fontsize=12, fontweight='bold')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# Right subplot: Selectivity profiles
ax2.plot(sol.t, selectivity_c, label='Selectivity of C', linewidth=2, color='tab:green')
ax2.plot(sol.t, selectivity_d, label='Selectivity of D', linewidth=2, color='tab:red')
ax2.set_xlabel('Time', fontsize=11)
ax2.set_ylabel('Selectivity (mol C or D / mol C+D)', fontsize=11)
ax2.set_title('Product Selectivity', fontsize=12, fontweight='bold')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

fig.suptitle(r'Competing Reactions: $A + B \rightarrow C$ and $A + 2B \rightarrow D$', 
             fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()