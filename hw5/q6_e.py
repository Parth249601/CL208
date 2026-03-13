import numpy as np
import matplotlib.pyplot as plt

# --- 1. Define Problem Parameters ---
UA = 3600         
Cpa = 40          
dH_rx = -80000    
Keq_400 = 100     
k_400 = 1         
E_R = 20000       
V = 10            
v0 = 1            
Fa0 = 10          
T0 = 310          
R_cal = 1.987     

tau = V / v0      

# --- 2. Core Functions ---
def k(T):
    return k_400 * np.exp(E_R * (1/400 - 1/T))

def Keq(T):
    # Retaining the corrected van 't Hoff equation
    return Keq_400 * np.exp((dH_rx / R_cal) * (1/400 - 1/T))

def X_MB(T):
    return (tau * k(T)) / (1 + tau * k(T) * (1 + 1/Keq(T)))

def G(T):
    return (-dH_rx) * X_MB(T)

def calc_Ta(T):
    """Calculates ambient temperature (Ta) directly from reactor temperature (T)."""
    return T + (Fa0 / UA) * (Cpa * (T - T0) - G(T))

# --- 3. Generate Data ---
# Use a finely spaced array for a smooth curve
T_reactor_range = np.linspace(300, 450, 5000)
Ta_range = calc_Ta(T_reactor_range)

# --- 4. Identify Ignition and Extinction Points ---
# Find the indices where the slope changes sign (local min/max of Ta)
dTa_dT = np.diff(Ta_range)

# Ignition point corresponds to the local maximum of Ta before the curve drops
idx_ignite = np.where((dTa_dT[:-1] > 0) & (dTa_dT[1:] < 0))[0][0] + 1
Ta_ignition = Ta_range[idx_ignite]
T_ignition = T_reactor_range[idx_ignite]

# Extinction point corresponds to the local minimum of Ta after the drop
idx_extinct = np.where((dTa_dT[:-1] < 0) & (dTa_dT[1:] > 0))[0][0] + 1
Ta_extinction = Ta_range[idx_extinct]
T_extinction = T_reactor_range[idx_extinct]

print(f"Ignition Temperature (Ta): {Ta_ignition:.1f} K")
print(f"Extinction Temperature (Ta): {Ta_extinction:.1f} K")

# --- 5. Plotting ---
plt.figure(figsize=(10, 7))

# Plot the main S-curve
plt.plot(Ta_range, T_reactor_range, color='#2E7D32', linewidth=2.5, label='Steady States')

# Mark Ignition Point
plt.plot(Ta_ignition, T_ignition, 'ro', markersize=8)
plt.annotate(f'Ignition\n$T_a$ = {Ta_ignition:.1f} K', 
             xy=(Ta_ignition, T_ignition), xytext=(Ta_ignition + 5, T_ignition - 15),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6),
             fontsize=11, fontweight='bold', ha='center')

# Mark Extinction Point
plt.plot(Ta_extinction, T_extinction, 'bo', markersize=8)
plt.annotate(f'Extinction\n$T_a$ = {Ta_extinction:.1f} K', 
             xy=(Ta_extinction, T_extinction), xytext=(Ta_extinction - 10, T_extinction + 15),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=6),
             fontsize=11, fontweight='bold', ha='center')

# Formatting
plt.title('Reactor Temperature vs. Ambient Temperature ($T_a$)', fontsize=16, fontweight='bold', pad=15)
plt.xlabel('Ambient Temperature, $T_a$ (K)', fontsize=14)
plt.ylabel('Reactor Temperature, $T$ (K)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12, loc='upper left')

# Highlight the multiple steady state region
plt.axvspan(Ta_extinction, Ta_ignition, color='gray', alpha=0.15, label='Multiple Steady States Region')

plt.tight_layout()
plt.show()