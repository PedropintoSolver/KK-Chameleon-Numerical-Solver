#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solve_chameleon_bvp.py

Nonlinear spherical BVP solver for a KK-Chameleon radion field using scipy.solve_bvp.
Outputs:
  r_m.txt        -> radius grid (meters)
  phi.txt        -> scalar field profile
  a_phi.txt      -> acceleration profile a_phi = -beta * dphi/dr (proxy units)
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import root_scalar
import os
import sys

# ----------------------- Model parameters (edit as needed) -----------------------

beta = 1.0                      # dimensionless matter coupling
alpha = 4.0 / np.sqrt(6)        # geometric exponential factor
Lambda = 1e-3                   # template energy scale (arbitrary units)
V0 = 0.5 * Lambda**4            # vacuum Casimir amplitude (arbitrary units)

rho_wall = 1.0                  # effective density in the wall
rho_vac = 1e-12                 # residual gas density inside cavity

# Geometry
R = 0.10                        # chamber radius in meters
r_min = 1e-6                    # avoid singularity at r=0

# Solver settings
rtol = 1e-8
atol = 1e-10
max_nodes = 8000

OUT_DIR = "bvp_output"


# ----------------------- Potential derivative components -----------------------

def dV_geom_dphi(phi):
    return -alpha * (Lambda**4) * np.exp(-alpha * phi)

def dV_void_dphi(phi):
    return alpha * V0 * np.exp(-alpha * phi)

def dV_baryon_dphi(phi, rho):
    return beta * rho * np.exp(beta * phi)

def dVeff_dphi(phi, rho):
    return dV_geom_dphi(phi) + dV_void_dphi(phi) + dV_baryon_dphi(phi, rho)


# ----------------------- Utility: find phi that solves dV_eff/dphi = 0 -----------------------

def find_phi_min(rho, phi_guess=1.0):
    def f(phi):
        return dVeff_dphi(phi, rho)

    brackets = [(-50, 50), (-20, 20), (-10, 10), (-5, 5), (0, 10)]
    for a, b in brackets:
        try:
            sol = root_scalar(f, bracket=[a, b], method='brentq', maxiter=200)
            if sol.converged:
                return sol.root
        except Exception:
            pass

    try:
        sol = root_scalar(f, x0=phi_guess, x1=phi_guess + 1.0,
                          method='secant', maxiter=200)
        if sol.converged:
            return sol.root
    except Exception:
        pass

    return phi_guess


# ----------------------- BVP differential system -----------------------

def fun(r, y):
    phi = y[0]
    dphidr = y[1]
    d2phidr2 = - (2.0 / r) * dphidr + dVeff_dphi(phi, rho_vac)
    return np.vstack((dphidr, d2phidr2))


# ----------------------- Boundary conditions -----------------------

def bc(ya, yb):
    return np.array([
        ya[1],          # dphi/dr(r_min) = 0
        yb[0] - phi_wall  # phi(R) = phi_wall
    ])


# ----------------------- Main script -----------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    global phi_wall, phi_vac
    phi_wall = find_phi_min(rho_wall, phi_guess=0.0)
    phi_vac = find_phi_min(rho_vac, phi_guess=2.0)

    print(f"phi_wall = {phi_wall:.6e}")
    print(f"phi_vac  = {phi_vac:.6e}")

    N = 800
    r = np.linspace(r_min, R, N)

    phi_init = phi_vac + (phi_wall - phi_vac) * (r / R)**2
    dphi_init = np.gradient(phi_init, r)
    y_init = np.vstack((phi_init, dphi_init))

    print("Starting BVP solver...")
    sol = solve_bvp(fun, bc, r, y_init, tol=rtol,
                    max_nodes=max_nodes, verbose=2)

    if not sol.success:
        print("Initial solve did not converge. Retrying with linear guess...")
        phi_init2 = phi_vac + (phi_wall - phi_vac) * (r / R)
        dphi_init2 = np.gradient(phi_init2, r)
        y_init2 = np.vstack((phi_init2, dphi_init2))

        sol = solve_bvp(fun, bc, r, y_init2, tol=rtol,
                        max_nodes=2 * max_nodes, verbose=2)

    if not sol.success:
        print("WARNING: Solver still did not converge.")
        print(sol.message)

    r_fine = np.linspace(r_min, R, 2000)
    y_fine = sol.sol(r_fine)

    phi_sol = y_fine[0]
    dphidr_sol = y_fine[1]
    a_phi = -beta * dphidr_sol

    np.savetxt(os.path.join(OUT_DIR, "r_m.txt"), r_fine)
    np.savetxt(os.path.join(OUT_DIR, "phi.txt"), phi_sol)
    np.savetxt(os.path.join(OUT_DIR, "a_phi.txt"), a_phi)

    print(f"Saved outputs to folder: {OUT_DIR}")
    print(f"a_phi range: {a_phi.min():.3e} ... {a_phi.max():.3e}")

    # Optional plotting
    try:
        import matplotlib.pyplot as plt
        plt.plot(r_fine, phi_sol)
        plt.xlabel("r (m)")
        plt.ylabel("phi(r)")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "phi_profile.png"), dpi=300)
        plt.close()

        plt.plot(r_fine, a_phi)
        plt.xlabel("r (m)")
        plt.ylabel("a_phi (proxy units)")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "a_phi_profile.png"), dpi=300)
        plt.close()
    except Exception:
        pass


if __name__ == "__main__":
    main()
