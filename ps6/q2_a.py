import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp



def reaction_system(V,y,k1,k2,v_0, I_b_0):
    F_a, F_b, F_d, F_u = y

    C_a = F_a/v_0
    C_b = F_b/v_0
    C_d = F_d/v_0
    C_u = F_u/v_0

    dFadV = -k2*C_a*(C_b**2) - k1*(C_a**2)*(C_b)
    dFbdV = I_b_0 - k2*C_a*(C_b**2) - k1*(C_a**2)*(C_b)
    dFddV = k1*(C_a**2)*C_b
    dFudV = k2*C_a*(C_b**2)

    return [dFadV, dFbdV, dFddV, dFudV]

#! Constants of the reaction

k1 = 2 #? L^2/(mol^2*sec)
k2 = 3 #? L^2/(mol^2*sec)
V_final = 50 #? L
V_span = (0, 50)  #? volume span for the simulation
V_eval = np.linspace(0, 50,1000)    #?volume points for evaluation
v_0 = 5 #? L/sec
I_b_0 = 4/50 #? mol/(s * L)

y_0 = [4.0, 0, 0, 0] #? Initial molar flow rates of A, B, D, and U


sol = solve_ivp(
    reaction_system,
    V_span,
    y_0,
    args=(k1, k2, v_0, I_b_0),
    t_eval=V_eval
)


selectivity_d = (sol.y[0][1:])/(sol.y[1][1:]) * (k1/k2)  #? Selectivity of D
selectivity_u = 1 / selectivity_d  #? Selectivity of U

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
# Left subplot: Molar flow rates
ax1.plot(sol.t, sol.y[0], label='$F_A$ (Reactant)', linewidth=2)
ax1.plot(sol.t, sol.y[1], label='$F_B$ (Reactant)', linewidth=2)
ax1.plot(sol.t, sol.y[2], label='$F_D$ (Product)', linewidth=2)
ax1.plot(sol.t, sol.y[3], label='$F_U$ (Product)', linewidth=2)
ax1.set_xlabel('Volume (L)', fontsize=11)
ax1.set_ylabel('Molar Flow Rate (mol/s)', fontsize=11)
ax1.set_title('Molar Flow Rates vs Volume', fontsize=12, fontweight='bold')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)
# Right subplot: Selectivity profiles
ax2.plot(sol.t[1:], selectivity_d, label='Selectivity of D', linewidth=2, color='tab:green')
ax2.plot(sol.t[1:], selectivity_u, label='Selectivity of U', linewidth=2, color='tab:red')
ax2.set_xlabel('Volume (L)', fontsize=11)
ax2.set_ylabel('Selectivity', fontsize=11)
ax2.set_title('Product Selectivity', fontsize=12, fontweight='bold')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)
fig.suptitle(r'Competing Reactions in a MR: $A + B \rightarrow D$ and $A + B \rightarrow U$', 
             fontsize=13, fontweight='bold')
plt.show()