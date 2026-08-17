/*
===============================================================================
2D Lid-Driven Cavity Flow CFD Solver
===============================================================================

Author: Tyler Warner
Description:
    A two-dimensional incompressible CFD solver for the lid-driven cavity
    problem using the Navier-Stokes equations and a finite-difference
    projection method.

    The solver uses a pressure Poisson equation to enforce incompressibility
    and OpenMP to parallelize computationally intensive grid operations.

Features:
    - 2D incompressible Navier-Stokes equations
    - Finite-difference discretization
    - Pressure Poisson solver
    - Velocity and pressure residual monitoring
    - Continuity error monitoring
    - OpenMP parallelization
    - CSV output for Python post-processing
    - Validation against lid-driven cavity benchmark data

Case:
    Reynolds Number: Re = 100
    Grid: 101 x 101
    Moving Lid Velocity: U = 1.0

===============================================================================
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <algorithm>

using namespace std;

int main() {
    // grid size
    const int Nx = 101;
    const int Ny = 101;
    // physical dimensions
    const double Lx = 1.0;
    const double Ly = 1.0;
    // grid spacing
    const double dx = Lx / (Nx - 1);
    const double dy = Ly / (Ny - 1);
    // fluid properties
    const double rho = 1.0; // density
    const double nu = 0.01; // kinematic viscosity
    const double U_lid = 1.0; // lid velocity
    const double Re = (U_lid * Lx) / nu; // Reynolds number
    // simulation settings
    const double dt = 0.001; // time step
    const int maxTimeSteps = 10000; // number of iterations
    const int pressureIterations = 200; // Jacobi iterations used for each time step to solve for pressure
    const double velocityTolerance = 1e-6; // tolerance for velocity convergence

    // arrays
    vector<vector<double>> u(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> v(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> uStar(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> vStar(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> uOld(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> vOld(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> p(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> pNew(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> rhs(Ny, vector<double>(Nx, 0.0));

    // starting information
    cout << "2D Lid-Driven Cavity Flow CFD Solver\n";
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "OpenMP threads available: " << omp_get_max_threads() << endl;

    // main time-stepping loop
    for (int step = 0; step < maxTimeSteps; step++) {
        // save previous velocity fields for convergence check
        uOld = u;
        vOld = v;

        // predictor step: compute intermediate velocities uStar and vStar
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // first derivatives
                double du_dx = (u[j][i + 1] - u[j][i - 1]) / (2 * dx);
                double du_dy = (u[j + 1][i] - u[j - 1][i]) / (2 * dy);
                double dv_dx = (v[j][i + 1] - v[j][i - 1]) / (2 * dx);
                double dv_dy = (v[j + 1][i] - v[j - 1][i]) / (2 * dy);
                // second derivatives
                double d2u_dx2 = (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / (dx * dx);
                double d2u_dy2 = (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / (dy * dy);
                double d2v_dx2 = (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / (dx * dx);
                double d2v_dy2 = (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / (dy * dy);
                // temporary velocities (uStar, vStar) without pressure correction
                uStar[j][i] = u[j][i] + dt * (-u[j][i] * du_dx - v[j][i] * du_dy + nu * (d2u_dx2 + d2u_dy2));
                vStar[j][i] = v[j][i] + dt * (-u[j][i] * dv_dx - v[j][i] * dv_dy + nu * (d2v_dx2 + d2v_dy2));
            }
    }

    // apply velocity boundary conditions for uStar and vStar
    for (int i = 0; i < Nx; i++) {
        // bottom wall (stationary)
        uStar[0][i] = 0.0;
        vStar[0][i] = 0.0;
    }
    for (int i = 0; i < Nx; i++) {
        // top wall (moving lid)
        uStar[Ny - 1][i] = U_lid;
        vStar[Ny - 1][i] = 0.0;
    }
    for (int j = 0; j < Ny; j++) {
        // left wall
        uStar[j][0] = 0.0;
        vStar[j][0] = 0.0;
        // right wall
        uStar[j][Nx - 1] = 0.0;
        vStar[j][Nx - 1] = 0.0;
    }

    // build pressure Poisson equation rhs
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            double duStar_dx = (uStar[j][i + 1] - uStar[j][i - 1]) / (2 * dx);
            double dvStar_dy = (vStar[j + 1][i] - vStar[j - 1][i]) / (2 * dy);
            rhs[j][i] = (rho / dt) * (duStar_dx + dvStar_dy);
        }
    }

    // pressure poisson solver
    for (int pressureStep = 0; pressureStep < pressureIterations; pressureStep++) {
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                pNew[j][i] = ((p[j][i + 1] + p[j][i - 1]) * dy * dy + (p[j + 1][i] + p[j - 1][i]) * dx * dx - rhs[j][i] * dx * dx * dy * dy) / (2.0 * (dx * dx + dy * dy));
            }
        }
        // apply pressure boundary conditions
        // left and right: dp/dx = 0
        for (int j = 0; j < Ny; j++) {
            pNew[j][0] = pNew[j][1]; // left wall
            pNew[j][Nx - 1] = pNew[j][Nx - 2]; // right wall
        }
        // bottom: dp/dy = 0
        for (int i = 0; i < Nx; i++) {
            pNew[0][i] = pNew[1][i]; // bottom wall
        }
        // top boundary as pressure reference
        for (int i = 0; i < Nx; i++) {
            pNew[Ny - 1][i] = 0.0;
        }

        p.swap(pNew); // update pressure field

        }

    // pressure correction
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            double dp_dx = (p[j][i + 1] - p[j][i - 1]) / (2 * dx);
            double dp_dy = (p[j + 1][i] - p[j - 1][i]) / (2 * dy);
            u[j][i] = uStar[j][i] - (dt / rho) * dp_dx;
            v[j][i] = vStar[j][i] - (dt / rho) * dp_dy;
        }
    }

    // final velocity boundary conditions
    // bottom wall
    for (int i = 0; i < Nx; i++) {
        u[0][i] = 0.0;
        v[0][i] = 0.0;
    }
    // top moving lid
    for (int i = 1; i < Nx - 1; i++) {
        u[Ny - 1][i] = U_lid;
        v[Ny - 1][i] = 0.0;
    }

    // left and right walls
    for (int j = 0; j < Ny; j++) {
        u[j][0] = 0.0;
        v[j][0] = 0.0;
        u[j][Nx - 1] = 0.0;
        v[j][Nx - 1] = 0.0;
    }

    // velocity residual
    double velocityResidual = 0.0;
    #pragma omp parallel for collapse(2) reduction(max:velocityResidual)
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            double du = fabs(u[j][i] - uOld[j][i]);
            double dv = fabs(v[j][i] - vOld[j][i]);
            double localResidual = max(du, dv);
            if (localResidual > velocityResidual) {
                velocityResidual = localResidual;
            }
        }
    }

    // continuity error
    double continuityError = 0.0;
    #pragma omp parallel for collapse(2) reduction(max:continuityError)
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            double du_dx = (u[j][i + 1] - u[j][i - 1]) / (2 * dx);
            double dv_dy = (v[j + 1][i] - v[j - 1][i]) / (2 * dy);
            double divergence = fabs(du_dx + dv_dy);
            if (divergence > continuityError) {
                continuityError = divergence;
            }
        }
    }

    // pressure residual
    double pressureResidual = 0.0;
    #pragma omp parallel for collapse(2) reduction(max:pressureResidual)
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            double laplacianP = (p[j][i + 1] - 2.0 * p[j][i] + p[j][i - 1]) / (dx * dx) + (p[j + 1][i] - 2.0 * p[j][i] + p[j - 1][i]) / (dy * dy);
            double localResidual = fabs(laplacianP - rhs[j][i]);

            if (localResidual > pressureResidual) {
                pressureResidual = localResidual;
            }
        }
    }

    // print convergence information
    if (step % 500 == 0) {
        cout << "\nIteration: " << step << endl;
        cout << "Velocity residual = " << velocityResidual << endl;
        cout << "Pressure Residual = " << pressureResidual << endl;
        cout << "Continuity error = " << continuityError << endl;
    }

    // optional convergence stop
    if (step > 1000 && velocityResidual < velocityTolerance) {
        cout << "\nVelocity field converged at iteration " << step << endl;
        break;
    }
}

int midX = Nx / 2;
int midY = Ny / 2;

cout << "Simulation completed.";
cout << "\nCenter:\n";
cout << "u = " << u[midY][midX] << endl;
cout << "v = " << v[midY][midX] << endl;
cout << "p = " << p[midY][midX] << endl;

cout << "\nUpper Center:\n";
cout << "u = " << u[3 * Ny / 4][midX] << endl;
cout << "v = " << v[3 * Ny / 4][midX] << endl;

cout << "\nLower center:\n";
cout << "u = " << u[Ny / 4][midX] << endl;
cout << "v = " << v[Ny / 4][midX] << endl;

cout << "\nLeft center:\n";
cout << "u = " << u[midY][Nx / 4] << endl;
cout << "v = " << v[midY][Nx / 4] << endl;

cout << "\nRight center:\n";
cout << "u = " << u[midY][3 * Nx / 4] << endl;
cout << "v = " << v[midY][3 * Nx / 4] << endl;

// results to csv
ofstream outputFile("cavity_results.csv");
outputFile << "x,y,u,v,p\n";
for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        double y = j * dy;
        outputFile
            << x << ","
            << y << ","
            << u[j][i] << ","
            << v[j][i] << ","
            << p[j][i] << "\n";
    }
}

outputFile.close();
cout << "\nResults writen to cavity_results.csv\n";

// write centerline data for validation
ofstream uCenterlineFile("u_centerline.csv");
uCenterlineFile << "y,u\n";
for (int j = 0; j < Ny; j++) {
    double y = j * dy;
    uCenterlineFile << y << "," << u[j][midX] << "\n";
}

uCenterlineFile.close();

// horizontal centerline: v(x) at y = 0.5
ofstream vCenterlineFile("v_centerline.csv");
vCenterlineFile << "x,v\n";
for (int i = 0; i < Nx; i++) {
    double x = i * dx;
    vCenterlineFile << x << "," << v[midY][i] << "\n";
}

vCenterlineFile.close();
cout << "\nCenterline data written to u_centerline.csv and v_centerline.csv\n";
return 0;
}

