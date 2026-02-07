import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def pbr_system(W, F, k1, k2, Cto, alpha):
    """
    Defines the system of ODEs for the PBR.
    W: Independent variable (Catalyst Weight)
    F: List of dependent variables [Fa, Fb, Fc, Fd]
    """
    Fa, Fb, Fc, Fd = F
    
    # 1. Calculate Total Molar Flow Rate
    Ft = Fa + Fb + Fc + Fd
    
    # Avoid division by zero if flow stops
    if Ft <= 0:
        return [0, 0, 0, 0]

    # 2. Calculate Pressure Drop Term (1 - alpha * W)
    # The equation in your image suggests: C = Cto * (Fi/Ft) * (1 - alpha*W)
    pressure_term = (1 - alpha * W)
    
    # physical constraint: pressure cannot be negative
    if pressure_term < 0:
        pressure_term = 0 

    # 3. Calculate Concentrations (Ci)
    # Using the formula from the handwritten note
    Ca = Cto * (Fa / Ft) * pressure_term
    Cb = Cto * (Fb / Ft) * pressure_term
    Cc = Cto * (Fc / Ft) * pressure_term
    Cd = Cto * (Fd / Ft) * pressure_term
    
    # 4. Define Reaction Rates (based on handwritten eq 2 and 3)
    # Rate 1 term: k1 * Ca * Cb^2
    r1 = k1 * Ca * (Cb**2)
    
    # Rate 2 term: k2 * Ca^2 * Cc^3
    r2 = k2 * (Ca**2) * (Cc**3)
    
    # 5. Define Differential Equations (dFi/dW)
    #Equation 1: dFa/dW = -r1 - r2 (assuming A is consumed in both reactions)

    #dFa_dW = -r1 - 1.5*r2  # Adjust coefficients based on your exact stoichiometry from the image

    # Equation 2: dFb/dW = -2 * r1
    dFb_dW = -2 * r1
    
    # Equation 3: dFc/dW = r1 - r2
    dFc_dW = r1 - r2
    
    # Equation 4: dFd/dW = (1/3) * r2
    dFd_dW = (1.0/3.0) * r2
    
    # Equation 1: dFa/dW
    # The image is slightly messy for A. Based on the logic of chemical stoichiometry:
    # If B is -2*r1 (consuming 2 B), A is likely consuming 1 A.
    # If D is +1/3*r2, and C is -1*r2, A is likely involved in r2 as well.
    # You can edit the coefficients below to match your exact derivation:
    dFa_dW = -1 * r1 - 1.5 * r2  
    
    return [dFa_dW, dFb_dW, dFc_dW, dFd_dW]

# --- PARAMETERS (You must replace these with your problem's actual values) ---
k1 = 100          # Rate constant 1
k2 = 1500        # Rate constant 2
Cto = 0.2       # Initial Total Concentration (mol/L)
alpha = 0.0019    # Pressure drop parameter (kg^-1) or similar units
W_final = 1000.0    # Final weight of catalyst (kg)

# Initial Molar Flow Rates (mol/s) - Replace with your values
Fa0 = 100.0
Fb0 = 100.0
Fc0 = 0.0
Fd0 = 0.0
F_init = [Fa0, Fb0, Fc0, Fd0]

# --- SOLVING THE ODEs ---
# Integration range: from W=0 to W_final
w_span = (0, W_final)

# Solve using 'Radau' or 'BDF' method as these equations can be stiff
solution = solve_ivp(
    fun=lambda w, y: pbr_system(w, y, k1, k2, Cto, alpha),
    t_span=w_span,
    y0=F_init,
    method='Radau', 
    dense_output=True
)

# --- PLOTTING RESULTS ---
W_plot = np.linspace(0, W_final, 1000)
F_sol = solution.sol(W_plot)

plt.figure(figsize=(10, 6))
plt.plot(W_plot, F_sol[0], label='$F_A$')
plt.plot(W_plot, F_sol[1], label='$F_B$')
plt.plot(W_plot, F_sol[2], label='$F_C$')
plt.plot(W_plot, F_sol[3], label='$F_D$')

plt.title('Molar Flow Rates vs Catalyst Weight')
plt.xlabel('Catalyst Weight (W)')
plt.ylabel('Molar Flow Rate (F)')
plt.legend()
plt.grid(True)
plt.show()

# Print final values
print(f"Final Flow Rates at W={W_final}:")
print(f"Fa: {F_sol[0][-1]:.4f}")
print(f"Fb: {F_sol[1][-1]:.4f}")
print(f"Fc: {F_sol[2][-1]:.4f}")
print(f"Fd: {F_sol[3][-1]:.4f}")