import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- 1. Define Problem Parameters ---
UA = 3600         # cal/(min * K)
Cpa = 40          # cal/(mol * K)
dH_rx = -80000    # cal/mol A
Keq_400 = 100     # Equilibrium constant at 400 K
k_400 = 1         # Rate constant at 400 K (min^-1)
E_R = 20000       # Activation energy / Gas constant (K)
V = 10            # Volume (dm^3)
v0 = 1            # Volumetric flow rate (dm^3/min)
Fa0 = 10          # Molar flow rate of A (mol/min)
T0 = 310          # Feed temperature (K) -> 37 C
Ta = 310          # Ambient temperature (K) -> 37 C
R_cal = 1.987     # Ideal gas constant (cal/(mol * K))

# Calculate space time
tau = V / v0      # min

# --- 2. Define Core Functions ---
def k(T):
    """Rate constant as a function of temperature."""
    return k_400 * np.exp(E_R * (1/400 - 1/T))

def Keq(T):
    """Equilibrium constant as a function of temperature."""
    return Keq_400 * np.exp((-dH_rx / R_cal) * (1/T- 1/400))

def X_MB(T):
    """Conversion from the Material Balance."""
    return (tau * k(T)) / (1 + tau * k(T) * (1 + 1/Keq(T)))

def G(T):
    """Heat Generated curve."""
    return (-dH_rx) * X_MB(T)

def R(T):
    """Heat Removed line."""
    return Cpa * (T - T0) + (UA / Fa0) * (T - Ta)

def energy_balance(T):
    """Objective function to find steady states (G(T) = R(T))."""
    return G(T) - R(T)

# --- 3. Solve for Steady States ---
# Provide initial guesses near the expected intersections
initial_guesses = [310, 380, 420]
steady_states = []

for guess in initial_guesses:
    root, info, ier, mesg = fsolve(energy_balance, guess, full_output=True)
    if ier == 1:
        T_ss = root[0]
        # Avoid appending duplicate roots due to convergence proximity
        if not any(np.isclose(T_ss, ss, atol=1.0) for ss in steady_states):
            steady_states.append(T_ss)

steady_states.sort()

# --- 4. Plotting ---
# Generate temperature data points for smooth curves
T_range = np.linspace(290, 450, 1000)
G_vals = G(T_range)
R_vals = R(T_range)

# Initialize a professional-looking figure
plt.figure(figsize=(10, 6))

# Plot G(T) and R(T) with distinct colors and thicker lines
plt.plot(T_range, G_vals, label='$G(T)$ (Heat Generated)', color='#D32F2F', linewidth=2.5)
plt.plot(T_range, R_vals, label='$R(T)$ (Heat Removed)', color='#1976D2', linewidth=2.5)

# Highlight and annotate the steady-state intersection points
for i, T_ss in enumerate(steady_states):
    # Plot the point
    plt.plot(T_ss, G(T_ss), 'ko', markersize=8, zorder=5)
    # Draw a vertical dashed line down to the x-axis
    plt.vlines(T_ss, 0, G(T_ss), linestyles='dashed', colors='black', alpha=0.5)
    
    # Adjust annotation text positioning to avoid overlapping the curves
    y_offset = -5000 if i == 0 else 3000
    plt.text(T_ss + 2, G(T_ss) + y_offset, f'{T_ss:.1f} K', fontsize=12, fontweight='bold', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.8))

# Formatting the plot
plt.title('Energy Balance: Heat Generation and Removal vs. Temperature', fontsize=16, fontweight='bold', pad=15)
plt.xlabel('Reactor Temperature, $T$ (K)', fontsize=14)
plt.ylabel('Heat (cal/mol A)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12, loc='upper left')

# Set axis limits for a clean view
plt.xlim(290, 450)
plt.ylim(0, 90000)

plt.tight_layout()
plt.show()

print("="*50)
print("Steady-State Temperatures (K):")
for i, T_ss in enumerate(steady_states):
    print(f"Steady State {i+1}: {T_ss:.2f} K")
print("="*50)
print("Finding the conversion at each steady state:")
for i, T_ss in enumerate(steady_states):
    X_ss = X_MB(T_ss)
    print(f"Steady State {i+1}: Conversion (X) = {X_ss:.4f}")