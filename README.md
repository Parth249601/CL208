# CL208 — Chemical Reaction Engineering

> **Semester 4 | Department of Chemical Engineering**

Computational solutions for CL208 (Chemical Reaction Engineering) coursework — encompassing problem sets, homework assignments, and in-class exercises. All simulations are implemented in **Python 3** using standard scientific computing libraries.

---

## Table of Contents

- [Dependencies](#dependencies)
- [Repository Structure](#repository-structure)
- [extra_class/](#extra_class)
- [hw2/](#hw2)
- [hw3/](#hw3)
- [ps3/](#ps3)
- [ps5/](#ps5)
- [Usage](#usage)

---

## Dependencies

| Package        | Purpose                                                             |
| -------------- | ------------------------------------------------------------------- |
| `numpy`        | Numerical array operations and linear algebra                       |
| `scipy`        | ODE solvers (`solve_ivp`) and statistical regression (`linregress`) |
| `matplotlib`   | Plotting and visualization                                          |
| `scikit-learn` | Linear regression fitting (used in PS3)                             |
| `pandas`       | Data handling (imported in PS3)                                     |

Install all dependencies:

```bash
pip install numpy scipy matplotlib scikit-learn pandas
```

---

## Repository Structure

```
CRE/
├── extra_class/
│   └── complex_reactions.py
├── hw2/
│   ├── 24B0395-HW2.pdf
│   ├── CRE-HW2-Part1.pdf
│   ├── q2.pdf
│   └── q5.pdf
├── hw3/
│   ├── q1.py
│   ├── q3.py
│   ├── q3_differential.py
│   ├── q1.png
│   ├── q3_differential.png
│   └── q3_integral.png
├── ps3/
│   ├── q1_differential.py
│   ├── q1_integral.py
│   ├── q2_integral.py
│   ├── q3_integral.py
│   └── plots/
│       ├── q1_diff_method_plot.png
│       └── q1_integral_method_plot.png
├── ps5/
│   ├── q1.py
│   ├── q2.py
│   ├── q1.png
│   ├── q2.png
│   └── q2_selectivity.png
└── README.md
```

---

## extra_class/

### `complex_reactions.py` — Packed Bed Reactor (PBR) with Multiple Reactions

Solves a system of coupled ODEs describing molar flow rates in a **Packed Bed Reactor (PBR)** with two competing reactions and pressure drop.

**Reactions modeled:**

```math
A + 2B \rightarrow C \qquad (r_1 = k_1 \, C_A \, C_B^2)
```

```math
2A + 3C \rightarrow D \qquad (r_2 = k_2 \, C_A^2 \, C_C^3)
```

**Key features:**

- Accounts for **pressure drop** via the Ergun equation parameter α, modifying concentrations as Cᵢ = C_T0 · (Fᵢ/F_T) · (1 − αW).
- Tracks four species: F_A, F_B, F_C, F_D as functions of catalyst weight W.
- Solved using `scipy.integrate.solve_ivp` with the **Radau** (implicit Runge–Kutta) method for stiff systems.
- Outputs molar flow rate profiles vs. catalyst weight.

**Parameters:**

| Symbol | Value       | Description                  |
| ------ | ----------- | ---------------------------- |
| k₁     | 100         | Rate constant for reaction 1 |
| k₂     | 1500        | Rate constant for reaction 2 |
| C_T0   | 0.2 mol/L   | Total inlet concentration    |
| α      | 0.0019 kg⁻¹ | Pressure drop parameter      |
| W      | 0 – 1000 kg | Catalyst weight range        |

---

## hw2/

Contains the **submitted PDF solutions** for Homework 2:

| File                | Description                  |
| ------------------- | ---------------------------- |
| `24B0395-HW2.pdf`   | Complete homework submission |
| `CRE-HW2-Part1.pdf` | Part 1 solutions             |
| `q2.pdf`            | Question 2 detailed solution |
| `q5.pdf`            | Question 5 detailed solution |

---

## hw3/

### `q1.py` — Rate Law Determination from CSTR Data (Differential Analysis)

Determines the **reaction order** and **rate constant** for a liquid-phase reaction using CSTR data.

**Method:** Log-log analysis of rate vs. concentration.

- Converts exit concentration data from a CSTR into reaction rates using:

```math
-r_A = \frac{v_0 \, C_{A0} - v_0 \left(1 - \frac{X}{2}\right) C_A}{V}
```

- Performs linear regression on ln(−r_A) vs. ln(C_A) to extract the reaction order (slope) and rate constant (intercept).
- Plots the log-log fit alongside experimental data points.

**Given Data:**

| v₀ (L/hr) | C_A (mol/L) |
| --------- | ----------- |
| 10        | 85.7        |
| 3         | 66.7        |
| 1.2       | 50          |
| 0.5       | 33.4        |

---

### `q3.py` — Reversible Isomerization Kinetics (Integral Method)

Analyzes the **reversible first-order isomerization** of _n_-pentane to _i_-pentane using the integral method.

```math
A \underset{k_2}{\overset{k_1}{\rightleftharpoons}} R
```

**Method:**

- Linearizes the integrated rate law:

```math
\ln\!\left(\frac{C_{R_e} - C_{R_0}}{C_{R_e} - C_R}\right) = (k_1 + k_2) \, t
```

- Uses `scipy.stats.linregress` to extract (k₁ + k₂) from the slope.
- Decomposes into individual constants using the equilibrium relation K_eq = k₁/k₂ = C_Re/C_Ae.

**Results reported:** k₁, k₂, K_eq, and the proposed rate equation −r_A = k₁·C_A − k₂·C_R.

---

### `q3_differential.py` — Reversible Isomerization Kinetics (Differential Method)

Same reaction system as `q3.py`, solved via the **differential method**.

**Method:**

- Computes dC_R/dt by finite differences at interval midpoints.
- Plots dC_R/dt vs. (C_Re − C_R); the slope yields (k₁ + k₂).
- Decomposes rate constants identically to the integral approach.

**Output plots:**

- `q3_differential.png` — Differential method fit.
- `q3_integral.png` — Integral method fit (from `q3.py`).

---

## ps3/

### Problem Set 3 — Rate Law Determination from Experimental Data

---

### `q1_differential.py` — Gas-Phase Decomposition (Differential Method)

Determines the kinetics of a **gas-phase decomposition** reaction:

```math
2A \rightarrow B
```

conducted at 100 °C in a constant-volume batch reactor with total pressure measured over time.

**Method:**

- Corrects initial pressure for temperature using Gay-Lussac's Law: P₀ = P_cold · (T_rxn / T_cold).
- Extracts partial pressure of A from total pressure: P_A = 2·P_T − P_A0.
- Converts to concentration: C_A = P_A / (RT).
- Applies the differential method: finite-difference −dC_A/dt vs. C_A on a log-log scale.
- Linear regression yields the **reaction order** n and **rate constant** k.

---

### `q1_integral.py` — Gas-Phase Decomposition (Integral Method)

Same reaction as `q1_differential.py`, analyzed via the **integral method**.

**Method:**

- Tests candidate rate laws (first-order and second-order) by plotting the appropriate linearized forms:
  - First order: ln(C_A0 / C_A) vs. t
  - Second order: (1/C_A − 1/C_A0) vs. t
- Evaluates goodness of fit via Pearson correlation (r) and R² from `sklearn.linear_model.LinearRegression`.

**Output plots:** `plots/q1_integral_method_plot.png`, `plots/q1_diff_method_plot.png`

---

### `q2_integral.py` — CSTR Rate Law from Space-Time Data

Determines the rate law for a **liquid-phase reaction in a CSTR** given space-time (τ) and exit concentration data.

**Given Data:**

| τ (s) | C_A (mol/dm³) |
| ----- | ------------- |
| 15    | 1.5           |
| 38    | 1.25          |
| 100   | 1.0           |
| 300   | 0.75          |
| 1200  | 0.5           |

**Method:**

- Computes rates from the CSTR design equation: −r_A = (C_A0 − C_A) / τ.
- Tests first, second, and third-order fits by plotting rate vs. C_Aⁿ.
- Best fit found for **third-order kinetics**: −r_A = k·C_A³.

---

### `q3_integral.py` — Batch Reactor Order Determination

Determines the reaction order for a **batch liquid-phase reaction** with concentration-time data.

**Given Data:**

| t (s) | C_A (mol/m³) |
| ----- | ------------ |
| 0     | 1000         |
| 100   | 500          |
| 200   | 333          |
| 300   | 250          |
| 400   | 200          |

**Method:**

- Tests linearized integral forms for first and second order.
- Confirms **second-order kinetics** from the linearity of (1/C_A − 1/C_A0) vs. t.
- Extracts the rate constant k from the slope.

---

## ps5/

### Problem Set 5 — Multiple Reactions and Selectivity

---

### `q1.py` — Series-Parallel Reversible Reactions (A ⇌ B ⇌ C)

Simulates a system of **consecutive reversible reactions** in a batch reactor:

```math
A \underset{k_1'}{\overset{k_1}{\rightleftharpoons}} B \underset{k_2'}{\overset{k_2}{\rightleftharpoons}} C
```

**Key features:**

- Solves three coupled ODEs for C_A, C_B, C_C using `solve_ivp`.
- Demonstrates the characteristic rise-and-fall behavior of the intermediate B.
- Plots concentration profiles of all three species over time.

**Parameters:**

| Constant | Value     |
| -------- | --------- |
| k₁       | 0.1       |
| k₁'      | 0.01      |
| k₂       | 0.003     |
| k₂'      | 0.003     |
| C_A0     | 1.0 mol/L |

---

### `q2.py` — Competing Reactions with Semi-Batch Feed and Selectivity Analysis

Models **parallel competing reactions** in a semi-batch reactor where species B is continuously fed:

```math
A + B \xrightarrow{k_1} C
```

```math
A + 2B \xrightarrow{k_2} D
```

**Key features:**

- Semi-batch operation: B is fed at volumetric flow rate v₀ into a reactor of volume V.
- Solves four coupled ODEs for C_A, C_B, C_C, C_D.
- Computes **instantaneous selectivity**:
  - S_C/D = k₁ / (k₂·C_B), favoring product C at low C_B.
  - S_D/C = (k₂·C_B) / k₁, favoring product D at high C_B.
- Side-by-side subplots display concentration profiles and selectivity evolution.

**Parameters:**

| Symbol | Value      | Description                |
| ------ | ---------- | -------------------------- |
| k₁     | 0.1        | Rate constant for A+B → C  |
| k₂     | 0.002      | Rate constant for A+2B → D |
| V      | 100 L      | Reactor volume             |
| v₀     | 1.0 L/time | Volumetric feed rate of B  |
| C_A0   | 10.0 mol/L | Initial concentration of A |
| C_B0   | 20.0 mol/L | Initial concentration of B |

**Output plots:** `q2.png`, `q2_selectivity.png`

---

## Usage

Run any script from the repository root:

```bash
python <folder>/<script>.py
```

For example:

```bash
python ps5/q2.py
python hw3/q3.py
python extra_class/complex_reactions.py
```

> All scripts generate `matplotlib` plots that display interactively. Close the plot window to allow the script to continue execution (if multiple plots are generated sequentially).

---

## License

This repository contains academic coursework and is intended for educational purposes only.
