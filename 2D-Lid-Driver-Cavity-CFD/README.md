# 2D Lid-Driven Cavity CFD Solver

**Author:** Tyler Warner

## Overview

This project implements a two-dimensional incompressible Computational Fluid Dynamics (CFD) solver in C++ for the classic lid-driven cavity problem.

The solver uses finite-difference discretization and a pressure-projection method to solve the incompressible Navier–Stokes equations. OpenMP is used to parallelize computationally intensive grid operations.

The current simulation uses a square cavity with a moving upper lid and operates at a Reynolds number of **Re = 100**.

---

## Governing Equations

The simulation solves the two-dimensional incompressible Navier–Stokes equations for horizontal velocity \(u\), vertical velocity \(v\), and pressure \(p\).

### X-Momentum Equation

$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\right)$$

### Y-Momentum Equation

$$\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial y} + \nu\left(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\right)$$

The terms in these equations represent:

- **Transient term** — change in velocity with time
- **Convective terms** — transport of momentum by the moving fluid
- **Pressure-gradient term** — acceleration caused by pressure gradients
- **Viscous terms** — diffusion of momentum due to fluid viscosity

### Continuity Equation

For incompressible flow, conservation of mass requires:

$$\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$$

This requires the velocity field to be divergence-free.

---

## Reynolds Number

The Reynolds number is defined as:

$$Re=\frac{UL}{\nu}$$

For this simulation:

$$Re=\frac{(1.0)(1.0)}{0.01}=100$$

where:

- \(U\) = lid velocity
- \(L\) = cavity length
- \(\nu\) = kinematic viscosity

---

## Numerical Method

The continuous governing equations are converted into discrete equations using finite-difference approximations on a **101 × 101 computational grid**.

For example, a first derivative is approximated using a central difference:

$$\frac{\partial u}{\partial x}\approx\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}$$

A second derivative is approximated by:

$$\frac{\partial^2u}{\partial x^2}\approx\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{\Delta x^2}$$

Similar approximations are used for the remaining derivatives in the momentum and pressure equations.

---

## Pressure Projection Method

The solver first calculates intermediate velocity fields \(u^*\) and \(v^*\) using the convection and viscous terms.

Because the predicted velocity field does not necessarily satisfy continuity, its divergence is used to construct the pressure Poisson equation:

$$\nabla^2 p = \frac{\rho}{\Delta t}\left(\frac{\partial u^{*}}{\partial x} + \frac{\partial v^{*}}{\partial y}\right)$$

After solving for pressure, the velocity field is corrected using:

$$u^{n+1}=u^*-\frac{\Delta t}{\rho}\frac{\partial p}{\partial x}$$

$$v^{n+1}=v^*-\frac{\Delta t}{\rho}\frac{\partial p}{\partial y}$$

This process drives the velocity field toward satisfying the incompressible continuity equation.

---

## Boundary Conditions

The upper wall moves in the positive x-direction with a velocity of:

$$U_{\text{lid}}=1.0$$

The remaining walls are stationary and use no-slip boundary conditions:

$$u=v=0$$

The motion of the upper wall transfers momentum to the fluid and produces the characteristic clockwise recirculation within the cavity.

---

## Parallel Computing

OpenMP is used to parallelize several computationally intensive grid operations.

The solver can be compiled with GCC using:

```bash
g++ -O2 -fopenmp 2D-liddriven-cavityflow.cpp -o cavity
