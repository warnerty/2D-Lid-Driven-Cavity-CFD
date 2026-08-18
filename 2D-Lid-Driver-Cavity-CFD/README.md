# 2D Lid-Driven Cavity CFD Solver

**Author:** Tyler Warner

A parallel C++ Computational Fluid Dynamics (CFD) solver for two-dimensional incompressible lid-driven cavity flow using finite-difference discretization, a pressure-projection method, and OpenMP.

## Overview

This project was developed to build a stronger understanding of the numerical methods underlying Computational Fluid Dynamics rather than relying solely on commercial CFD software.

The solver models the classic **2D lid-driven cavity problem**. A square cavity is filled with an incompressible fluid while the upper wall moves horizontally at a constant velocity. Momentum transferred from the moving lid produces a dominant recirculating vortex within the cavity.

The current case is solved at:

$$Re = 100$$

using a **101 × 101 computational grid**.

The project combines:

- C++
- Incompressible Navier–Stokes equations
- Finite-difference methods
- Pressure–velocity coupling
- OpenMP parallel programming
- Numerical convergence monitoring
- Python post-processing
- Benchmark validation

---

## Physical Problem

The computational domain is a unit square:

$$L_x = L_y = 1$$

The upper wall moves in the positive x-direction:

$$U_{\text{lid}} = 1$$

while the remaining walls are stationary.

The fluid properties used are:

$$\rho = 1$$

$$\nu = 0.01$$

Therefore, the Reynolds number is:

$$Re = \frac{U_{\text{lid}}L}{\nu} = \frac{(1)(1)}{0.01} = 100$$

The moving upper wall transfers momentum to the fluid, producing the characteristic clockwise recirculation of the lid-driven cavity problem.

---

## Governing Equations

The solver approximates the two-dimensional incompressible Navier–Stokes equations.

### X-Momentum Equation

$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\right)$$

### Y-Momentum Equation

$$\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial y} + \nu\left(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\right)$$

where:

| Symbol | Description |
|:---:|---|
| $u$ | Horizontal velocity |
| $v$ | Vertical velocity |
| $p$ | Pressure |
| $\rho$ | Fluid density |
| $\nu$ | Kinematic viscosity |
| $t$ | Time |

### Continuity Equation

For incompressible flow, conservation of mass requires:

$$\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$$

or equivalently:

$$\nabla\cdot\mathbf{u}=0$$

The velocity field must therefore remain divergence-free.

---

## Finite-Difference Discretization

The continuous governing equations are converted into algebraic equations on a **101 × 101 grid**.

For the unit-square domain:

$$\Delta x = \Delta y = \frac{1}{100} = 0.01$$

Spatial derivatives are approximated using neighboring grid points.

For example, a central-difference approximation for the first derivative is:

$$\frac{\partial u}{\partial x}\approx\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}$$

and the second derivative is approximated as:

$$\frac{\partial^2u}{\partial x^2}\approx\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{\Delta x^2}$$

Similar approximations are used for the remaining derivatives appearing in the momentum and pressure equations.

---

## Pressure-Projection Method

The solver uses a predictor-corrector / projection approach rather than directly solving velocity and pressure simultaneously.

### 1. Velocity Prediction

Intermediate velocities $u^\ast$ and $v^\ast$ are first calculated from the convection and viscous-diffusion terms.

Conceptually:

```text
Current Velocity
       ↓
Convection + Diffusion
       ↓
Predicted Velocity (u*, v*)
```

### 2. Pressure Poisson Equation

The predicted velocity field does not necessarily satisfy incompressible continuity.

Its divergence is therefore used to construct a pressure Poisson equation:

$$\nabla^2p=\frac{\rho}{\Delta t}\left(\frac{\partial u^\ast}{\partial x}+\frac{\partial v^\ast}{\partial y}\right)$$

The pressure equation is solved iteratively.

### 3. Velocity Correction

The resulting pressure gradient is used to correct the predicted velocity:

$$u^{n+1}=u^\ast-\frac{\Delta t}{\rho}\frac{\partial p}{\partial x}$$

$$v^{n+1}=v^\ast-\frac{\Delta t}{\rho}\frac{\partial p}{\partial y}$$

This correction drives the numerical velocity field toward satisfying:

$$\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$$

The overall numerical procedure is therefore:

```text
u, v
 ↓
Velocity Prediction
 ↓
u*, v*
 ↓
Calculate Divergence
 ↓
Pressure Poisson Equation
 ↓
Pressure Field
 ↓
Velocity Correction
 ↓
Apply Boundary Conditions
 ↓
Check Convergence
 ↓
Repeat
```

---

## Boundary Conditions

No-slip boundary conditions are applied to the stationary walls.

### Bottom Wall

$$u = v = 0$$

### Left and Right Walls

$$u = v = 0$$

### Moving Lid

$$u = U_{\text{lid}} = 1,\qquad v = 0$$

The upper-wall motion drives the circulation throughout the cavity.

---

## Parallel Programming

OpenMP is used to parallelize computationally intensive grid operations.

Parallel regions use directives such as:

```cpp
#pragma omp parallel for collapse(2)
```

This allows independent grid calculations to be distributed among multiple CPU threads rather than being performed entirely sequentially.

The program can be compiled using GCC with:

```bash
g++ -O2 -fopenmp 2D-liddriven-cavityflow.cpp -o cavity
```

where:

- `-O2` enables compiler optimization.
- `-fopenmp` enables OpenMP support.

---

## Convergence Monitoring

The solver monitors changes in the numerical solution rather than relying only on a predetermined iteration count.

The velocity residual is based on the maximum change in velocity between successive iterations:

$$R_u=\max\left(\left|u^{n+1}-u^n\right|,\left|v^{n+1}-v^n\right|\right)$$

The solver also evaluates continuity error using:

$$\left|\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right|$$

which provides a numerical measure of how closely the velocity field satisfies incompressible mass conservation.

A pressure-equation residual is also monitored as an additional convergence diagnostic.

---

## Post-Processing

The C++ solver exports the computed solution to CSV format:

```text
x, y, u, v, p
```

Python is then used separately for post-processing and visualization using:

- NumPy
- pandas
- Matplotlib

The resulting visualizations include:

- Pressure field
- Velocity magnitude
- Streamlines
- Velocity vectors
- Horizontal centerline velocity
- Vertical centerline velocity

This keeps the numerical solver and visualization workflow separate:

```text
C++ CFD Solver
      ↓
   CSV Data
      ↓
Python Post-Processing
      ↓
Visualization + Validation
```

---

## Validation

A numerical solution should not be considered accurate simply because the resulting flow field appears physically reasonable.

For this reason, the solver was compared against the published lid-driven cavity benchmark data of **Ghia, Ghia, and Shin (1982)**.

At $Re=100$, two centerline velocity profiles were examined:

$$u(0.5,y)$$

and

$$v(x,0.5)$$

The calculated C++ profiles were plotted against the corresponding benchmark values.

The resulting profiles show reasonable agreement with the published data, providing an initial validation that the solver reproduces the expected lid-driven cavity behavior.

---

## Current Simulation Parameters

| Parameter | Value |
|---|---:|
| Domain | 1 × 1 |
| Grid | 101 × 101 |
| Lid Velocity | 1.0 |
| Density | 1.0 |
| Kinematic Viscosity | 0.01 |
| Reynolds Number | 100 |
| Time Step | 0.001 |

---

## Future Work

Potential extensions of the project include:

- Quantitative error analysis against benchmark data
- Grid-convergence studies
- OpenMP speedup and parallel-efficiency testing
- Higher Reynolds number simulations
- Improved spatial discretization
- Improved pressure solvers
- More advanced pressure–velocity coupling methods

---

## Reference

Benchmark velocity data used for validation:

U. Ghia, K. N. Ghia, and C. T. Shin, “High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method,” *Journal of Computational Physics*, vol. 48, no. 3, pp. 387–411, 1982. DOI: 10.1016/0021-9991(82)90058-4.

---

## Disclaimer

This project is an educational CFD implementation developed to study numerical methods, C++, fluid dynamics, and parallel computing. It is not intended to replace validated production CFD software.
