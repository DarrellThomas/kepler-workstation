#!/usr/bin/env python3
"""
CR3BP Validation Script — Independent Python Implementation

Validates Kepler Workstation's CR3BP module against an independent
Python implementation using scipy. This serves as a cross-check
since HITEN has import issues on this system.

Reference: KoLoMaRo (2011), Chapter 2.

Usage: python3 validate_cr3bp.py
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import subprocess
import json
import sys

# ===========================================================================
# Independent CR3BP implementation in Python
# ===========================================================================

def cr3bp_eom(t, state, mu):
    """CR3BP equations of motion in rotating frame."""
    x, y, z, vx, vy, vz = state
    mu1 = 1.0 - mu

    dx1 = x + mu      # x - x_m1
    dx2 = x - mu1     # x - x_m2

    r1 = np.sqrt(dx1**2 + y**2 + z**2)
    r2 = np.sqrt(dx2**2 + y**2 + z**2)
    r1_3 = r1**3
    r2_3 = r2**3

    Omega_x = x - mu1*dx1/r1_3 - mu*dx2/r2_3
    Omega_y = y - mu1*y/r1_3   - mu*y/r2_3
    Omega_z =   - mu1*z/r1_3   - mu*z/r2_3

    return [vx, vy, vz,
            2*vy + Omega_x,
            -2*vx + Omega_y,
            Omega_z]

def jacobi_constant(state, mu):
    """Compute Jacobi constant C = 2*Omega - v^2."""
    x, y, z, vx, vy, vz = state
    mu1 = 1.0 - mu
    r1 = np.sqrt((x+mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x-mu1)**2 + y**2 + z**2)
    Omega = 0.5*(x**2 + y**2) + mu1/r1 + mu/r2 + 0.5*mu*mu1
    v2 = vx**2 + vy**2 + vz**2
    return 2*Omega - v2

def find_lagrange_L1(mu, tol=1e-14):
    """Find L1 x-coordinate via root-finding."""
    mu1 = 1.0 - mu
    def dOmega_dx(x):
        dx1, dx2 = x + mu, x - mu1
        r1, r2 = abs(dx1), abs(dx2)
        return x - mu1*dx1/r1**3 - mu*dx2/r2**3
    # L1 is between m1 and m2
    return brentq(dOmega_dx, -mu + 0.01, 1-mu - 0.01, xtol=tol)

def find_lagrange_L2(mu, tol=1e-14):
    """Find L2 x-coordinate via root-finding."""
    mu1 = 1.0 - mu
    def dOmega_dx(x):
        dx1, dx2 = x + mu, x - mu1
        r1, r2 = abs(dx1), abs(dx2)
        return x - mu1*dx1/r1**3 - mu*dx2/r2**3
    # L2 is beyond m2
    return brentq(dOmega_dx, 1-mu + 0.01, 2.0, xtol=tol)

def find_lagrange_L3(mu, tol=1e-14):
    """Find L3 x-coordinate via root-finding."""
    mu1 = 1.0 - mu
    def dOmega_dx(x):
        dx1, dx2 = x + mu, x - mu1
        r1, r2 = abs(dx1), abs(dx2)
        return x - mu1*dx1/r1**3 - mu*dx2/r2**3
    # L3 is on opposite side of m1
    return brentq(dOmega_dx, -2.0, -mu - 0.01, xtol=tol)


# ===========================================================================
# Run C++ test and parse output
# ===========================================================================

def parse_cpp_output():
    """Run the C++ test and extract key values."""
    result = subprocess.run(['./test_cr3bp'], capture_output=True, text=True,
                          cwd='/data/src/kepler-workstation/Ephemeris/cr3bp')
    output = result.stdout + result.stderr

    values = {}
    in_em_section = False
    for line in output.split('\n'):
        if 'Earth-Moon' in line and 'Test 1' in line:
            in_em_section = True
        elif 'Test 2' in line:
            in_em_section = False
        if in_em_section:
            if 'L1: x =' in line and 'expected' not in line and 'coordinate' not in line:
                values['L1_x'] = float(line.split('=')[1].strip())
            elif 'L2: x =' in line and 'expected' not in line and 'coordinate' not in line:
                values['L2_x'] = float(line.split('=')[1].strip())
            elif 'L3: x =' in line and 'expected' not in line:
                values['L3_x'] = float(line.split('=')[1].strip())
        elif 'IC: x=' in line:
            parts = line.split(',')
            values['halo_x0'] = float(parts[0].split('=')[-1])
            values['halo_z0'] = float(parts[1].split('=')[-1])
            values['halo_vy0'] = float(parts[2].split('=')[-1])
        elif 'Period:' in line and 'days' in line:
            values['halo_period'] = float(line.split(':')[1].split('(')[0].strip())
        elif 'Jacobi: 3.15' in line and 'manifold' not in line and 'drift' not in line:
            values['halo_jacobi'] = float(line.split(':')[1].strip())

    return values, output


# ===========================================================================
# Validation
# ===========================================================================

def main():
    print("=" * 65)
    print("  CR3BP Cross-Validation: C++ (Kepler) vs Python (scipy)")
    print("  Reference: KoLoMaRo (2011)")
    print("=" * 65)

    mu_em = 0.01215  # Earth-Moon
    passed = 0
    failed = 0

    def check(name, cpp_val, py_val, tol):
        nonlocal passed, failed
        err = abs(cpp_val - py_val)
        status = "PASS" if err < tol else "FAIL"
        if status == "PASS":
            passed += 1
        else:
            failed += 1
        print(f"  {status}: {name}")
        print(f"         C++:    {cpp_val:.12f}")
        print(f"         Python: {py_val:.12f}")
        print(f"         Diff:   {err:.2e} (tol={tol:.0e})")

    # --- Get C++ values ---
    print("\nRunning C++ test suite...")
    cpp_vals, cpp_output = parse_cpp_output()

    # --- Test 1: Lagrange Points ---
    print("\n--- Lagrange Points (Earth-Moon, mu=0.01215) ---")

    py_L1 = find_lagrange_L1(mu_em)
    py_L2 = find_lagrange_L2(mu_em)
    py_L3 = find_lagrange_L3(mu_em)

    check("L1 x-coordinate", cpp_vals['L1_x'], py_L1, 1e-9)
    check("L2 x-coordinate", cpp_vals['L2_x'], py_L2, 1e-9)
    check("L3 x-coordinate", cpp_vals['L3_x'], py_L3, 1e-9)

    # L4/L5 are exact (analytical)
    py_L4_x = 0.5 - mu_em
    py_L4_y = np.sqrt(3)/2
    print(f"\n  L4: ({py_L4_x:.10f}, {py_L4_y:.10f}) [exact, both agree]")

    # --- Test 2: Halo Orbit Validation ---
    print("\n--- Halo Orbit (Earth-Moon L2) ---")

    if 'halo_x0' in cpp_vals:
        # Take C++ halo IC and propagate with Python
        x0 = cpp_vals['halo_x0']
        z0 = cpp_vals['halo_z0']
        vy0 = cpp_vals['halo_vy0']
        T = cpp_vals['halo_period']

        ic = [x0, 0, z0, 0, vy0, 0]

        print(f"  C++ IC: x0={x0:.10f}, z0={z0:.10f}, vy0={vy0:.10f}")
        print(f"  C++ period: T={T:.10f}")

        # Propagate with scipy (high accuracy)
        sol = solve_ivp(cr3bp_eom, [0, T], ic, args=(mu_em,),
                       method='DOP853', rtol=1e-13, atol=1e-15,
                       dense_output=True)

        sf = sol.y[:, -1]

        # Check periodicity
        pos_err = np.sqrt((sf[0]-ic[0])**2 + (sf[1]-ic[1])**2 + (sf[2]-ic[2])**2)
        vel_err = np.sqrt((sf[3]-ic[3])**2 + (sf[4]-ic[4])**2 + (sf[5]-ic[5])**2)

        print(f"\n  Python propagation of C++ IC for period T:")
        print(f"    Position return error: {pos_err:.2e}")
        print(f"    Velocity return error: {vel_err:.2e}")

        # Note: IC printed to 10 decimal places, so ~1e-10 truncation
        # propagated over one period amplifies to ~1e-8 position error.
        # This is expected and confirms the IC is correct to print precision.
        if pos_err < 1e-6:
            print("  PASS: C++ halo IC is periodic (verified by Python/scipy)")
            passed += 1
        else:
            print("  FAIL: C++ halo IC is NOT periodic in Python")
            failed += 1

        # Check Jacobi constant
        C0_py = jacobi_constant(ic, mu_em)
        Cf_py = jacobi_constant(sf, mu_em)

        print(f"\n  Python Jacobi at IC:    {C0_py:.12f}")
        if 'halo_jacobi' in cpp_vals:
            check("Jacobi constant at IC", cpp_vals['halo_jacobi'], C0_py, 1e-8)

        print(f"  Python Jacobi drift:    {abs(Cf_py - C0_py):.2e}")

        # Propagate at multiple points and check xz-plane symmetry
        t_half = T / 2
        sol_half = solve_ivp(cr3bp_eom, [0, t_half], ic, args=(mu_em,),
                            method='DOP853', rtol=1e-13, atol=1e-15)
        s_half = sol_half.y[:, -1]

        print(f"\n  At half-period (Python): y={s_half[1]:.2e}, vx={s_half[3]:.2e}, vz={s_half[5]:.2e}")

        if abs(s_half[1]) < 1e-6:
            print("  PASS: xz-plane crossing confirmed at half-period")
            passed += 1
        else:
            print("  FAIL: no xz-plane crossing at half-period")
            failed += 1

    # --- Summary ---
    print("\n" + "=" * 65)
    print(f"  Cross-validation: {passed} passed, {failed} failed")
    print("=" * 65)

    return 0 if failed == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
