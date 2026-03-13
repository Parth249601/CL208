import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- 1. System Constants ---
T0_feed = 350.0        # Initial feed temperature [K]
T_cool = 350.0         # Temperature to cool down to between stages [K]
R = 8.314              # Ideal gas constant [J/mol*K]
CpA = 45.0             # Heat capacity [J/mol*K]
delta_H = -30000.0     # Heat of reaction [J/mol]

# Set your final desired conversion to determine how many stages you need
X_target_overall = 0.65

# --- 2. Equilibrium Curve Function ---
def calc_Xeq(T):
    Kc = 30 * np.exp((-delta_H / R) * (1 / T - 1 / 400))
    return Kc / (1 + Kc)

# --- 3. Stage-by-Stage Calculation ---
stages_T = [] # To store temperatures for plotting
stages_X = [] # To store conversions for plotting

current_X = 0.0
current_T = T0_feed
stage_count = 0

print(f"Targeting Overall Conversion of {X_target_overall}\n")
print("-" * 50)

while current_X < X_target_overall:
    stage_count += 1
    
    # 1. Define the operating line equation for THIS specific stage
    # T_adiabatic = T_in + (-delta_H / CpA) * (X - X_in)
    def find_stage_equilibrium(X):
        T_adiabatic = current_T - (delta_H / CpA) * (X - current_X)
        return X - calc_Xeq(T_adiabatic)
    
    # 2. Find intersection with the equilibrium curve
    X_eq_stage = fsolve(find_stage_equilibrium, current_X + 0.1)[0]
    
    # 3. PFR achieves 90% of this stage's adiabatic equilibrium conversion
    X_out = 0.9 * X_eq_stage
    
    # 4. Temperature at the exit of this PFR
    T_out = current_T - (delta_H / CpA) * (X_out - current_X)
    
    # Store reaction path for plotting (Slanted line up and right)
    stages_T.append([current_T, T_out])
    stages_X.append([current_X, X_out])
    
    print(f"Stage {stage_count} (Reaction): X goes from {current_X:.4f} to {X_out:.4f}")
    print(f"                   T goes from {current_T:.1f} K to {T_out:.1f} K")
    
    # 5. Update for next stage (Interstage Cooling)
    current_X = X_out
    
    # If we haven't reached the overall target, cool it down for the next stage
    if current_X < X_target_overall:
        stages_T.append([T_out, T_cool])        # Cooling path (Horizontal line left)
        stages_X.append([current_X, current_X]) # Constant X during cooling
        current_T = T_cool
        print(f"Stage {stage_count} (Cooling) : T cooled back to {T_cool:.1f} K at X = {current_X:.4f}")
        print("-" * 50)

print("-" * 50)
print(f"Total Stages Required: {stage_count}")

# --- 4. Plotting ---
T_range = np.linspace(300, 800, 200)
X_eq_range = calc_Xeq(T_range)

plt.figure(figsize=(10, 6))

# Plot the equilibrium curve (Blue)
plt.plot(T_range, X_eq_range, 'b-', linewidth=2, label='Equilibrium Conversion ($X_e$)')

# Plot the sawtooth stages
for i in range(len(stages_T)):
    if i % 2 == 0:
        # Reaction path (Green dashed)
        label = 'Adiabatic Reaction Path' if i == 0 else ""
        plt.plot(stages_T[i], stages_X[i], 'g--', linewidth=2, label=label)
    else:
        # Cooling path (Red solid)
        label = 'Interstage Cooling' if i == 1 else ""
        plt.plot(stages_T[i], stages_X[i], 'r-', linewidth=2, label=label)

# Formatting
plt.axvline(x=T0_feed, color='gray', linestyle=':', label=f'Feed Temp ({T0_feed} K)')
plt.axhline(y=X_target_overall, color='purple', linestyle=':', label=f'Overall Target X ({X_target_overall})')

plt.xlabel('Temperature (K)', fontsize=12)
plt.ylabel('Conversion ($X$)', fontsize=12)
plt.title('Multistage PFR with Interstage Cooling', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.5)
plt.xlim(300, 800)
plt.ylim(0, 1.05)
plt.tight_layout()
plt.show()