# 2D Lid-Driven Cavity CFD Solver

## Governing Equations

This simulation models two-dimensional, incompressible flow using the Navier–Stokes equations. The solution consists of the horizontal velocity component $u$, vertical velocity component $v$, and pressure $p$.

### Conservation of Momentum

The momentum equation in the **x-direction** is

$$
\frac{\partial u}{\partial t}
+ u\frac{\partial u}{\partial x}
+ v\frac{\partial u}{\partial y}
=
-\frac{1}{\rho}\frac{\partial p}{\partial x}
+ \nu
\left(
\frac{\partial^2 u}{\partial x^2}
+
\frac{\partial^2 u}{\partial y^2}
\right).
$$

The momentum equation in the **y-direction** is

$$
\frac{\partial v}{\partial t}
+ u\frac{\partial v}{\partial x}
+ v\frac{\partial v}{\partial y}
=
-\frac{1}{\rho}\frac{\partial p}{\partial y}
+ \nu
\left(
\frac{\partial^2 v}{\partial x^2}
+
\frac{\partial^2 v}{\partial y^2}
\right).
$$

The terms in these equations represent:

- **Transient term** — change in velocity with time
- **Convective terms** — transport of momentum by the moving fluid
- **Pressure-gradient term** — acceleration caused by pressure gradients
- **Viscous terms** — diffusion of momentum due to fluid viscosity

### Conservation of Mass

Because the fluid is incompressible, the velocity field must satisfy the continuity equation:

$$
\frac{\partial u}{\partial x}
+
\frac{\partial v}{\partial y}
= 0.
$$

Equivalently,

$$
\nabla \cdot \mathbf{u} = 0,
$$

meaning the velocity field is divergence-free and conserves mass throughout the cavity.

### Reynolds Number

The flow regime is characterized by the Reynolds number:

$$
Re = \frac{UL}{\nu},
$$

where $U$ is the lid velocity, $L$ is the cavity length, and $\nu$ is the kinematic viscosity.

For the current simulation,

$$
Re
=
\frac{(1.0)(1.0)}{0.01}
=
100.
$$

### Finite-Difference Discretization

The continuous Navier–Stokes equations are converted into algebraic equations that can be evaluated on the computational grid.

For example, the first derivative of horizontal velocity in the $x$-direction is approximated using a central difference:

$$
\left.\frac{\partial u}{\partial x}\right|_{i,j}
\approx
\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}.
$$

The second derivative is approximated by

$$
\left.\frac{\partial^2 u}{\partial x^2}\right|_{i,j}
\approx
\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{\Delta x^2}.
$$

Similar finite-difference expressions are used for the remaining spatial derivatives in the momentum and pressure equations.

### Pressure Projection

An intermediate velocity field, $u^*$ and $v^*$, is first calculated from the convection and viscous terms.

Because this intermediate field does not necessarily satisfy continuity, its divergence is used to construct the pressure Poisson equation:

$$
\nabla^2 p
=
\frac{\rho}{\Delta t}
\left(
\frac{\partial u^*}{\partial x}
+
\frac{\partial v^*}{\partial y}
\right).
$$

After solving for pressure, the intermediate velocities are corrected using the pressure gradients:

$$
u^{n+1}
=
u^*
-
\frac{\Delta t}{\rho}
\frac{\partial p}{\partial x},
$$

$$
v^{n+1}
=
v^*
-
\frac{\Delta t}{\rho}
\frac{\partial p}{\partial y}.
$$

This pressure-correction step drives the velocity field toward satisfying the incompressible continuity equation.

### Variables

| Symbol | Description |
|:---:|---|
| $u$ | Velocity in the $x$-direction |
| $v$ | Velocity in the $y$-direction |
| $p$ | Pressure |
| $\rho$ | Fluid density |
| $\nu$ | Kinematic viscosity |
| $U$ | Moving-lid velocity |
| $L$ | Cavity length |
| $Re$ | Reynolds number |
| $\Delta t$ | Numerical time step |
| $\Delta x,\Delta y$ | Grid spacing |
