#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solve_chameleon_bvp.py

Solve a nonlinear spherical boundary-value problem (BVP) for a chameleon-like
Kaluza-Klein radion using scipy.integrate.solve_bvp.

Outputs (text files):
  /mnt/data/bvp_output/r_m.txt    -> radius grid (m)
  /mnt/data/bvp_output/phi.txt    -> phi(r) profile (dimensionless)
  /mnt/data/bvp_output/a_phi.txt  -> a_phi(r) = -beta * dphi/dr (proxy units)

Notes:
- The potential used here is a flexible template. If you need strict SI
  calibration, adjust Lambda, V0, rho_wall, rho_vac to consistent units and
  provide the physical mapping of phi -> phi_phys / Mpl.
- The code enforces phi(R) = phi_wall as an outer boundary (thin-shell).
- Regularity at r=0 is approximated by starting at r_min > 0 and enforcing dphi/dr(r_min)=0.
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import root_scalar
import os
import sys

# ----------------------------- User parameters (edit as needed) -----------------------------
beta = 1.0                # dimensionless coupling to matter
alpha = 4.0 / np.sqrt(6) # geometric exponent factor (KK reduction)
Lambda = 1e-3            # energy scale (template units)
V0 = 0.5 * Lambda**4     # void Casimir amplitude (same units as Lambda^4)
# Densities (in template / dimensionless units consistent with potential)
rho_wall = 1.0           # wall effective density (template)
rho_vac = 1e-12          # vacuum residual gas density inside chamber (template)

# Geometry
R = 0.10                 # chamber radius in meters (e.g., 0.10 m = 10 cm)
r_min = 1e-6             # small cutoff to avoid singularity at r=0

# Solver tolerances and settings
rtol = 1e-8
atol = 1e-10
max_nodes = 8000

# Output directory
OUT_DIR = "/mnt/data/bvp_output"

# ----------------------------- Potential derivative functions -----------------------------
def dV_geom_dphi(phi):
    """Derivative of geometric exponential term: d/dphi [ Lambda^4 e^{-alpha phi} ]"""
    return -alpha * (Lambda**4) * np.exp(-alpha * phi)

def dV_void_dphi(phi):
    """Derivative of void (fermionic Casimir) term: -V0 e^{-alpha phi} -> derivative = +alpha V0 e^{-alpha phi}"""
    return alpha * V0 * np.exp(-alpha * phi)

def dV_baryon_dphi(phi, rho):
    """Derivative of baryonic coupling term: rho * exp(beta phi) -> beta * rho * exp(beta phi)"""
    return beta * rho * np.exp(beta * phi)

def dVeff_dphi(phi, rho):
    """Total derivative of effective potential dV_eff/dphi"""
    return dV_geom_dphi(phi) + dV_void_dphi(phi) + dV_baryon_dphi(phi, rho)

# ----------------------------- Utility: find local minimum phi_c for given density -----------------------------
def find_phi_min(rho, phi_guess=1.0):
    """
    Find phi such that dVeff_dphi(phi, rho) = 0 using bracketing.
    Returns phi_min (float). If root finding fails, returns phi_guess.
    """
    def f(phi):
        return dVeff_dphi(phi, rho)

    # Try several brackets for robustness
    brackets = [(-50, 50), (-20, 20), (-10, 10), (-5, 5), (0, 10)]
    for a, b in brackets:
        try:
            sol = root_scalar(f, bracket=[a, b], method='brentq', maxiter=200)
            if sol.converged:
                return sol.root
        except Exception:
            continue

    # Try secant with guesses
    try:
        sol = root_scalar(f, x0=phi_guess, x1=phi_guess + 1.0, method='secant', maxiter=200)
        if sol.converged:
            return sol.root
    except Exception:
        pass

    # If all fails, return the guess (user should check this)
    return phi_guess

# ----------------------------- BVP ODE system -----------------------------
def fun(r, y):
    """
    r: array of radii
    y: array shape (2, len(r)) with y[0]=phi, y[1]=dphi/dr
    returns dy/dr as (2, len(r))
    Spherical equation:
      d2phi/dr2 + (2/r) dphi/dr = dVeff/dphi(phi, rho(r))
    -> d2phi/dr2 = - (2/r) dphi/dr + dVeff_dphi
    Here we solve inside the chamber (rho = rho_vac).
    """
    phi = y[0]
    dphidr = y[1]
    # Use rho_vac inside chamber; if you need exterior region, extend domain
    d2phidr2 = - (2.0 / r) * dphidr + dVeff_dphi(phi, rho_vac)
    return np.vstack((dphidr, d2phidr2))

# ----------------------------- Boundary conditions -----------------------------
def bc(ya, yb):
    """
    ya = y(r_min), yb = y(R)
    Enforce:
      - regularity at center: dphi/dr (r_min) = 0 (approx center)
      - thin-shell: phi(R) = phi_wall
    Returns residual array of length 2.
    """
    return np.array([ya[1], yb[0] - phi_wall])

# ----------------------------- Main execution -----------------------------
def main():
    # Create output folder
    os.makedirs(OUT_DIR, exist_ok=True)

    # Find minima for wall and vacuum
    global phi_wall, phi_vac
    phi_wall = find_phi_min(rho_wall, phi_guess=0.0)
    phi_vac = find_phi_min(rho_vac, phi_guess=3.0)
    print(f"Computed minima: phi_wall = {phi_wall:.12e}, phi_vac = {phi_vac:.12e}")

    # Initial mesh & guess
    N = 800
    r = np.linspace(r_min, R, N)

    # Initial guess interpolation between center value phi_vac and boundary phi_wall
    # Use a smooth quadratic profile as initial guess
    phi_init = phi_vac + (phi_wall - phi_vac) * (r / R)**2
    dphi_init = np.gradient(phi_init, r)
    y_init = np.vstack((phi_init, dphi_init))

    # Solve BVP
    print("Starting solve_bvp...")
    sol = solve_bvp(fun, bc, r, y_init, tol=rtol, max_nodes=max_nodes, verbose=2)

    # Fallback attempt with linear initial guess if solver didn't converge
    if not sol.success:
        print("Initial solve_bvp did not converge. Retrying with relaxed initial guess...")
        r2 = np.linspace(r_min, R, N)
        phi_init2 = phi_vac + (phi_wall - phi_vac) * (r2 / R)
        dphi_init2 = np.gradient(phi_init2, r2)
        y_init2 = np.vstack((phi_init2, dphi_init2))
        sol = solve_bvp(fun, bc, r2, y_init2, tol=rtol, max_nodes=2*max_nodes, verbose=2)
        if not sol.success:
            print("Warning: solver did not converge after retry. Message:", sol.message, file=sys.stderr)

    # Evaluate solution on fine grid
    r_fine = np.linspace(r_min, R, 2000)
    y_fine = sol.sol(r_fine)
    phi_sol = y_fine[0]
    dphidr_sol = y_fine[1]

    # Compute acceleration profile (proxy): a_phi = -beta * dphi/dr
    a_phi = -beta * dphidr_sol

    # Save outputs
    np.savetxt(os.path.join(OUT_DIR, "r_m.txt"), r_fine)
    np.savetxt(os.path.join(OUT_DIR, "phi.txt"), phi_sol)
    np.savetxt(os.path.join(OUT_DIR, "a_phi.txt"), a_phi)
    print(f"Saved outputs to {OUT_DIR}")

    # Print short summary
    print(f"r range: {r_fine[0]:.3e} - {r_fine[-1]:.3e} m")
    print(f"phi range: {phi_sol.min():.6e} .. {phi_sol.max():.6e}")
    print(f"a_phi range (proxy units): {a_phi.min():.6e} .. {a_phi.max():.6e}")

    # Optional: quick plotting (uncomment if you have a display or want files)
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6,4))
        plt.plot(r_fine*100, phi_sol, label=r'$\phi(r)$')
        plt.xlabel('r (cm)')
        plt.ylabel('phi (dimensionless)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "phi_profile.png"), dpi=300)
        plt.close()

        plt.figure(figsize=(6,4))
        plt.plot(r_fine*100, a_phi, label=r'$a_\\phi(r)$')
        plt.xlabel('r (cm)')
        plt.ylabel('a_phi (proxy units)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "a_phi_profile.png"), dpi=300)
        plt.close()

        print(f"Saved quick-plot PNGs to {OUT_DIR}")
    except Exception:
        # plotting optional; ignore if headless environment has no display
        pass

if __name__ == "__main__":
    main()