"""
2D Lid-Driven Cavity CFD Post-Processing
Author: Tyler Warner

Visualization and validation tools for the C++ lid-driven cavity CFD solver.
"""

import pandas as pd
import matplotlib.pyplot as plt

# read cfd results
u_data = pd.read_csv("u_centerline.csv")
v_data = pd.read_csv("v_centerline.csv")

# ghia et al. (1982) benchmark data for comparison
# Re = 100
# u velocity along x = 0.5

ghia_y = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 
         0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]

ghia_u = [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, 
         -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000]

# ghia et al. (1982) benchmark data for comparison
# v velocity along y = 0.5

ghia_x = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 
         0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]

ghia_v = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, 
         -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000]

# plot u centerline
plt.figure(figsize=(7, 7))
plt.plot(u_data["u"], u_data["y"], label="C++ Solver")
plt.scatter(ghia_u, ghia_y, label="Ghia et al. (1982)")
plt.xlabel("u velocity")
plt.ylabel("y")
plt.title("Vertical Centerline Velocity\n"
          "Lid-Driven Cavity Re = 100")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("validation_u_centerline.png", dpi = 300)

# plot v centerline
plt.figure(figsize=(7, 7))
plt.plot(v_data["x"], v_data["v"], label="C++ Solver")
plt.scatter(ghia_x, ghia_v, label="Ghia et al. (1982)")
plt.xlabel("x")
plt.ylabel("v velocity")
plt.title(
    "Horizontal Centerline Velocity\n"
    "Lid-Driven Cavity, Re = 100"
)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(
    "validation_v_centerline.png",
    dpi=300
)
plt.show()