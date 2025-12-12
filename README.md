# Effect of Wall Suction/Injection on the Laminar Boundary Layer Over a Flat Plate  
### A Numerical Study Using RK4, Finite Difference Method, and the Shooting Technique

This repository contains the complete C++ numerical implementation for solving the **modified Blasius boundary-layer equation** under the influence of **uniform wall suction or injection**. The work is based on the term paper *“Effect of Wall Suction/Injection on the Laminar Boundary Layer over a Flat Plate: A Numerical Study”* submitted to IIT Delhi.

The project investigates how wall mass flux affects velocity profiles, boundary-layer thickness, wall shear, and skin-friction coefficient. The numerical solution is obtained using:

- **Shooting Method** (Secant-based)  
- **4th-order Runge–Kutta (RK4) method**  
- **Finite Difference Method (FDM)** for verification  

---

## 📘 Problem Summary

The classical Blasius similarity equation is modified to include wall mass flux:

\[
f''' + \frac{1}{2} f f'' = 0
\]

Boundary conditions:

\[
f(0) = S, \quad f'(0) = 0, \quad f'(\infty) = 1
\]

Where:  
- **S > 0 → Suction** (flow drawn into the wall)  
- **S < 0 → Injection** (fluid blown out of the wall)

Suction reduces boundary-layer thickness and increases wall shear, while injection produces the opposite effect.

---

## 🧮 Numerical Methods Implemented

### **1. Shooting Method (Secant Approach)**
Solves the boundary value problem by guessing the unknown curvature \( f''(0) \) and iteratively refining it until:

\[
f'(\eta_{\text{max}}) \approx 1
\]

### **2. Runge–Kutta 4 (RK4)**
The third-order ODE is converted into a first-order system:

\[
f' = g, \quad g' = h, \quad h' = -\frac{1}{2} f h
\]

RK4 advances the system over the similarity domain.

### **3. Finite Difference Method (FDM)**
Used for cross-validation of the RK4 solver.

---

## 📊 Key Results

From the simulations:

### **Velocity Profiles**
- Suction → boundary layer becomes thinner  
- Injection → boundary layer thickens  

### **Wall Shear \( f''(0) \)**
- Increases with suction  
- Decreases with injection  
- Matches Blasius benchmark:  
  \[
  f''(0) \approx 0.33206 \quad \text{for } S = 0
  \]

### **Skin-Friction Coefficient**
\[
C_f = \frac{2 f''(0)}{\sqrt{Re_x}}
\]
- Highest for suction  
- Lowest for injection  

### **Boundary-Layer Thickness**
Defined by:
\[
f'(\eta_\delta) = 0.99
\]
- Thickest for injection  
- Thinnest for suction  

---

## 🛠️ Code Structure
project/
│── src/
│   ├── blasius_rk4.cpp          # Shooting method + RK4 solver
│   ├── blasius_fdm.cpp          # FDM implementation
│   └── utils.h / utils.cpp       # Derivative functions
│
│── data/
│   ├── velocity_profiles/
│   ├── shear_values/
│   └── boundary_layer_thickness/
│
│── plots/
│   ├── velocity_profiles.png
│   ├── fpp_vs_S.png
│   ├── Cf_vs_S.png
│   └── eta_delta_vs_S.png
│
└── term_paper.pdf

---

## 🚀 Compilation & Execution

### **Compile**
```bash
g++ blasius_rk4.cpp -o rk4
g++ blasius_fdm.cpp -o fdm
./rk4
./fdm
