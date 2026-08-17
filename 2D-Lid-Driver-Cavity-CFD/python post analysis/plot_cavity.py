"""
2D Lid-Driven Cavity CFD Post-Processing
Author: Tyler Warner

Visualization and validation tools for the C++ lid-driven cavity CFD solver.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# read cfd results

data = pd.read_csv("cavity_results.csv")

x_values = np.sort(data["x"].unique())
y_values = np.sort(data["y"].unique())

Nx = len(x_values)
Ny = len(y_values)

X = data["x"].values.reshape(Ny, Nx)
Y = data["y"].values.reshape(Ny, Nx)

U = data["u"].values.reshape(Ny, Nx)
V = data["v"].values.reshape(Ny, Nx)
P = data["p"].values.reshape(Ny, Nx)

speed = np.sqrt(U**2 + V**2)


# velocity vector plot
plt.figure(figsize=(7, 7))

skip = 5

plt.quiver(
    X[::skip, ::skip],
    Y[::skip, ::skip],
    U[::skip, ::skip],
    V[::skip, ::skip]
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Lid-Driven Cavity Velocity Vectors")

plt.axis("equal")
plt.tight_layout()

plt.savefig("velocity_vectors.png", dpi=300)


# streamline plot
plt.figure(figsize=(7, 7))

plt.streamplot(
    X,
    Y,
    U,
    V,
    density=1.5
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Lid-Driven Cavity Streamlines")

plt.axis("equal")
plt.tight_layout()

plt.savefig("streamlines.png", dpi=300)


# velocity magnitude contour
plt.figure(figsize=(7, 7))

velocity_contour = plt.contourf(
    X,
    Y,
    speed,
    levels=30
)

plt.colorbar(
    velocity_contour,
    label="Velocity Magnitude"
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Velocity Magnitude")

plt.axis("equal")
plt.tight_layout()

plt.savefig("velocity_magnitude.png", dpi=300)


# pressure contour 
plt.figure(figsize=(7, 7))

pressure_contour = plt.contourf(
    X,
    Y,
    P,
    levels=30
)

plt.colorbar(
    pressure_contour,
    label="Pressure"
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Pressure Field")

plt.axis("equal")
plt.tight_layout()

plt.savefig("pressure.png", dpi=300)

# display all plots
plt.show()

print("Plots created:")
print("velocity_vectors.png")
print("streamlines.png")
print("velocity_magnitude.png")
print("pressure.png")