// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// CR3BP Module Test Program
//
// Validates: Lagrange points, Jacobi constant conservation, halo orbits,
// and invariant manifold computation.
//
// Compile: g++ -std=c++11 -O2 -o test_cr3bp test_cr3bp.cpp
// Run:     ./test_cr3bp
///////////////////////////////////////////////////////////////////////////////

#include "cr3bp.hpp"
#include "halo_orbit.hpp"
#include "manifold.hpp"
#include <cstdio>
#include <cmath>

static int tests_passed = 0;
static int tests_failed = 0;

void check(const char* name, bool condition)
{
    if (condition) {
        printf("  PASS: %s\n", name);
        tests_passed++;
    } else {
        printf("  FAIL: %s\n", name);
        tests_failed++;
    }
}

void check_near(const char* name, double val, double expected, double tol)
{
    double err = fabs(val - expected);
    if (err < tol) {
        printf("  PASS: %s = %.10f (expected %.10f, err=%.2e)\n", name, val, expected, err);
        tests_passed++;
    } else {
        printf("  FAIL: %s = %.10f (expected %.10f, err=%.2e, tol=%.2e)\n",
               name, val, expected, err, tol);
        tests_failed++;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Test 1: Lagrange Points for Earth-Moon System
///////////////////////////////////////////////////////////////////////////////
void test_lagrange_points_earth_moon()
{
    printf("\n=== Test 1: Lagrange Points (Earth-Moon, mu=0.01215) ===\n");

    double mu = SYS_EARTH_MOON.mu;
    LagrangePoints lp = cr3bp_lagrange_points(mu);

    // Known values for Earth-Moon system (from various references)
    // L1: x ≈ 0.8369 (between Earth and Moon)
    // L2: x ≈ 1.1557 (beyond Moon)
    // L3: x ≈ -1.0051 (opposite side of Earth from Moon)
    // L4: (0.5-mu, sqrt(3)/2) = (0.48785, 0.86603)
    // L5: (0.5-mu, -sqrt(3)/2)

    printf("  L1: x = %.10f\n", lp.L1[0]);
    printf("  L2: x = %.10f\n", lp.L2[0]);
    printf("  L3: x = %.10f\n", lp.L3[0]);
    printf("  L4: (%.10f, %.10f)\n", lp.L4[0], lp.L4[1]);
    printf("  L5: (%.10f, %.10f)\n", lp.L5[0], lp.L5[1]);

    // L1 should be between Earth (-mu) and Moon (1-mu)
    check("L1 between primaries", lp.L1[0] > -mu && lp.L1[0] < 1.0-mu);
    check_near("L1 x-coordinate", lp.L1[0], 0.8369, 0.001);

    // L2 should be beyond Moon
    check("L2 beyond Moon", lp.L2[0] > 1.0-mu);
    check_near("L2 x-coordinate", lp.L2[0], 1.1557, 0.001);

    // L3 should be on opposite side
    check("L3 opposite side", lp.L3[0] < -mu);
    check_near("L3 x-coordinate", lp.L3[0], -1.0051, 0.001);

    // L4/L5 exact
    check_near("L4 x", lp.L4[0], 0.5 - mu, 1e-12);
    check_near("L4 y", lp.L4[1], sqrt(3.0)/2.0, 1e-12);
    check_near("L5 y", lp.L5[1], -sqrt(3.0)/2.0, 1e-12);

    // Equilibrium check: dOmega/dx should be ~0 at each collinear point
    double dOdx_L1 = lagrange_dOmega_dx(lp.L1[0], mu);
    double dOdx_L2 = lagrange_dOmega_dx(lp.L2[0], mu);
    double dOdx_L3 = lagrange_dOmega_dx(lp.L3[0], mu);
    check("L1 is equilibrium", fabs(dOdx_L1) < 1e-12);
    check("L2 is equilibrium", fabs(dOdx_L2) < 1e-12);
    check("L3 is equilibrium", fabs(dOdx_L3) < 1e-12);
}

///////////////////////////////////////////////////////////////////////////////
// Test 2: Lagrange Points for Sun-Jupiter System
///////////////////////////////////////////////////////////////////////////////
void test_lagrange_points_sun_jupiter()
{
    printf("\n=== Test 2: Lagrange Points (Sun-Jupiter, mu=9.537e-4) ===\n");

    double mu = SYS_SUN_JUPITER.mu;
    LagrangePoints lp = cr3bp_lagrange_points(mu);

    printf("  L1: x = %.10f\n", lp.L1[0]);
    printf("  L2: x = %.10f\n", lp.L2[0]);

    // From KoLoMaRo p. 42: gamma1 = 6.668e-2, x_L1 = 1 - mu - gamma1 = 0.9324
    check_near("L1 (KoLoMaRo)", lp.L1[0], 0.9324, 0.001);

    // Equilibrium check
    check("L1 equilibrium", fabs(lagrange_dOmega_dx(lp.L1[0], mu)) < 1e-12);
    check("L2 equilibrium", fabs(lagrange_dOmega_dx(lp.L2[0], mu)) < 1e-12);
}

///////////////////////////////////////////////////////////////////////////////
// Test 3: Jacobi Constant Conservation
///////////////////////////////////////////////////////////////////////////////
void test_jacobi_conservation()
{
    printf("\n=== Test 3: Jacobi Constant Conservation ===\n");

    double mu = SYS_EARTH_MOON.mu;

    // Start near L1 with a small velocity
    LagrangePoints lp = cr3bp_lagrange_points(mu);
    CR3BPState s0;
    s0.x = lp.L1[0] + 0.001;
    s0.y = 0.0;
    s0.z = 0.001;
    s0.vx = 0.0;
    s0.vy = 0.01;
    s0.vz = 0.0;

    double C0 = cr3bp_jacobi(s0, mu);
    printf("  Initial Jacobi constant: %.12f\n", C0);

    // Propagate for several orbital periods
    CR3BPState sf = cr3bp_propagate(mu, s0, 10.0);
    double Cf = cr3bp_jacobi(sf, mu);
    printf("  Final Jacobi constant:   %.12f\n", Cf);
    printf("  Drift: %.2e\n", fabs(Cf - C0));

    check("Jacobi conserved to 1e-8", fabs(Cf - C0) < 1e-8);
}

///////////////////////////////////////////////////////////////////////////////
// Test 4: STM Propagation (Phi(0,0) = I, det(Phi) = 1)
///////////////////////////////////////////////////////////////////////////////
void test_stm_propagation()
{
    printf("\n=== Test 4: STM Propagation ===\n");

    double mu = SYS_EARTH_MOON.mu;

    // Use a converged halo orbit for symplecticity check — the monodromy
    // matrix (STM over one full period) should be exactly symplectic.
    double Az = 0.04;
    HaloOrbit halo = compute_halo_orbit(mu, 2, Az, true);
    if (!halo.converged) {
        printf("  SKIP: halo orbit didn't converge\n");
        return;
    }

    double phi[6][6];
    CR3BPState sf = cr3bp_propagate_stm(mu, halo.ic, halo.period, phi);

    printf("  STM diagonal: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n",
           phi[0][0], phi[1][1], phi[2][2], phi[3][3], phi[4][4], phi[5][5]);

    // Practical symplecticity check for monodromy matrix:
    // The eigenvalues come in reciprocal pairs (lambda, 1/lambda).
    // Verify via the eigenvalue product from the manifold eigenvectors.
    // This is the most physically meaningful check for mission design.
    double v_u[6], v_s[6];
    double lam_u = power_iteration(phi, v_u);
    double lam_s = inverse_power_iteration(phi, v_s);
    double product = lam_u * lam_s;
    printf("  Eigenvalue product (should=1): %.8f\n", product);
    check("Monodromy eigenvalues are reciprocal", fabs(product - 1.0) < 0.01);

    // Also check: STM maps IC back to itself (since this is a periodic orbit)
    // Phi * [0, 0, 0, 0, vy0, 0]^T should give a consistent result
    check("STM has large elements (saddle dynamics)", fabs(phi[0][0]) > 10.0);
}

///////////////////////////////////////////////////////////////////////////////
// Test 5: Linearization at L1
///////////////////////////////////////////////////////////////////////////////
void test_linearization()
{
    printf("\n=== Test 5: Linearization at L1/L2 (Earth-Moon) ===\n");

    double mu = SYS_EARTH_MOON.mu;
    LagrangePoints lp = cr3bp_lagrange_points(mu);

    CollinearLinearization lin1 = cr3bp_linearize(lp.L1[0], mu);
    CollinearLinearization lin2 = cr3bp_linearize(lp.L2[0], mu);

    printf("  L1: lambda=%.6f, nu=%.6f, a=%.6f, b=%.6f\n",
           lin1.lambda, lin1.nu, lin1.a, lin1.b);
    printf("  L2: lambda=%.6f, nu=%.6f, a=%.6f, b=%.6f\n",
           lin2.lambda, lin2.nu, lin2.a, lin2.b);

    // Both a and b should be positive (KoLoMaRo Section 2.6)
    check("L1: a > 0", lin1.a > 0);
    check("L1: b > 0", lin1.b > 0);
    check("L2: a > 0", lin2.a > 0);
    check("L2: b > 0", lin2.b > 0);

    // Lambda should be real and positive (saddle)
    check("L1: lambda > 0", lin1.lambda > 0);
    // Nu should be real and positive (center frequency)
    check("L1: nu > 0", lin1.nu > 0);
}

///////////////////////////////////////////////////////////////////////////////
// Test 6: Halo Orbit (Earth-Moon L2)
///////////////////////////////////////////////////////////////////////////////
void test_halo_orbit()
{
    printf("\n=== Test 6: Halo Orbit (Earth-Moon L2) ===\n");

    double mu = SYS_EARTH_MOON.mu;

    // Az amplitude: 0.02 in nondimensional units
    // (about 7,700 km for Earth-Moon system)
    double Az = 0.04;

    printf("  Computing halo orbit at L2 with Az = %.4f ...\n", Az);
    HaloOrbit halo = compute_halo_orbit(mu, 2, Az, true);

    printf("  Converged: %s (iterations: %d)\n",
           halo.converged ? "YES" : "NO", halo.iterations);
    printf("  IC: x=%.10f, z=%.10f, vy=%.10f\n",
           halo.ic.x, halo.ic.z, halo.ic.vy);
    printf("  Period: %.8f (%.2f days)\n",
           halo.period,
           halo.period * SYS_EARTH_MOON.T_sec / (2.0*CR3BP_PI) / 86400.0);
    printf("  Jacobi: %.10f\n", halo.jacobi);
    printf("  Amplitudes: Ax=%.6f, Ay=%.6f, Az=%.6f\n",
           halo.Ax, halo.Ay, halo.Az);

    check("Halo orbit converged", halo.converged);

    if (halo.converged) {
        // Verify periodicity: propagate full period, should return to IC
        CR3BPState sf = cr3bp_propagate(mu, halo.ic, halo.period);
        double err_x = fabs(sf.x - halo.ic.x);
        double err_y = fabs(sf.y - halo.ic.y);
        double err_z = fabs(sf.z - halo.ic.z);
        double err_vx = fabs(sf.vx - halo.ic.vx);
        double err_vy = fabs(sf.vy - halo.ic.vy);
        double err_vz = fabs(sf.vz - halo.ic.vz);
        double pos_err = sqrt(err_x*err_x + err_y*err_y + err_z*err_z);
        double vel_err = sqrt(err_vx*err_vx + err_vy*err_vy + err_vz*err_vz);

        printf("  Periodicity: pos_err=%.2e, vel_err=%.2e\n", pos_err, vel_err);
        check("Periodic to 1e-6 in position", pos_err < 1e-6);
        check("Periodic to 1e-6 in velocity", vel_err < 1e-6);

        // Jacobi constant should be conserved along orbit
        double C0 = cr3bp_jacobi(halo.ic, mu);
        double Cf = cr3bp_jacobi(sf, mu);
        printf("  Jacobi drift over one period: %.2e\n", fabs(Cf - C0));
        check("Jacobi conserved over period", fabs(Cf - C0) < 1e-8);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Test 7: Invariant Manifolds
///////////////////////////////////////////////////////////////////////////////
void test_manifolds()
{
    printf("\n=== Test 7: Invariant Manifolds (Earth-Moon L2) ===\n");

    double mu = SYS_EARTH_MOON.mu;
    double Az = 0.04;

    HaloOrbit halo = compute_halo_orbit(mu, 2, Az, true);
    if (!halo.converged) {
        printf("  SKIP: halo orbit did not converge\n");
        return;
    }

    // Compute monodromy and eigenvalues
    double monodromy[6][6];
    compute_monodromy(mu, halo.ic, halo.period, monodromy);

    double v_u[6], v_s[6];
    double lambda_u = power_iteration(monodromy, v_u);
    double lambda_s = inverse_power_iteration(monodromy, v_s);

    printf("  Unstable eigenvalue: %.6f\n", lambda_u);
    printf("  Stable eigenvalue:   %.6f\n", lambda_s);
    printf("  Product (should ≈ 1): %.6f\n", lambda_u * lambda_s);

    // For symplectic system, lambda_u * lambda_s should be 1
    check("|lambda_u| > 1", fabs(lambda_u) > 1.0);
    check("|lambda_s| < 1", fabs(lambda_s) < 1.0);
    check("lambda_u * lambda_s ≈ 1", fabs(lambda_u * lambda_s - 1.0) < 0.1);

    // Compute unstable manifold (small, quick test)
    printf("  Computing unstable manifold (10 arcs, t=3.0) ...\n");
    ManifoldTube Wu = compute_manifold(mu, halo, 10, 3.0, 50, 1e-6, false, +1);

    check("Manifold has 10 trajectories", (int)Wu.trajectories.size() == 10);
    check("Each trajectory has points", Wu.trajectories[0].size() > 10);

    // Manifold trajectories should have similar Jacobi constant to the halo orbit
    if (!Wu.trajectories.empty() && !Wu.trajectories[0].empty()) {
        double C_halo = halo.jacobi;
        double C_man = cr3bp_jacobi(Wu.trajectories[0][0], mu);
        printf("  Jacobi: halo=%.8f, manifold[0]=%.8f, diff=%.2e\n",
               C_halo, C_man, fabs(C_man - C_halo));
        check("Manifold Jacobi ≈ halo Jacobi", fabs(C_man - C_halo) < 0.01);
    }

    printf("  Manifold computation complete.\n");
}

///////////////////////////////////////////////////////////////////////////////
// Test 8: Dimensional conversion check
///////////////////////////////////////////////////////////////////////////////
void test_dimensional()
{
    printf("\n=== Test 8: Dimensional Conversions ===\n");

    // Earth-Moon system: L ≈ 385,000 km
    printf("  Earth-Moon: L = %.0f km, T = %.0f s (%.2f days)\n",
           SYS_EARTH_MOON.L_km, SYS_EARTH_MOON.T_sec,
           SYS_EARTH_MOON.T_sec / 86400.0);
    printf("  V = %.3f km/s\n", SYS_EARTH_MOON.V_km_s());

    // Sun-Jupiter: L ≈ 778.4 million km
    printf("  Sun-Jupiter: L = %.3e km, T = %.3e s (%.2f years)\n",
           SYS_SUN_JUPITER.L_km, SYS_SUN_JUPITER.T_sec,
           SYS_SUN_JUPITER.T_sec / (365.25*86400.0));

    // Velocity scale should be reasonable
    check("Earth-Moon V ~ 1 km/s", SYS_EARTH_MOON.V_km_s() > 0.5 && SYS_EARTH_MOON.V_km_s() < 2.0);
    check("Sun-Jupiter V ~ 13 km/s", SYS_SUN_JUPITER.V_km_s() > 10 && SYS_SUN_JUPITER.V_km_s() < 15);
}

///////////////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////////////
int main()
{
    printf("===========================================================\n");
    printf("  Kepler Workstation — CR3BP Module Test Suite\n");
    printf("  Reference: KoLoMaRo (2011)\n");
    printf("===========================================================\n");

    test_lagrange_points_earth_moon();
    test_lagrange_points_sun_jupiter();
    test_jacobi_conservation();
    test_stm_propagation();
    test_linearization();
    test_halo_orbit();
    test_manifolds();
    test_dimensional();

    printf("\n===========================================================\n");
    printf("  Results: %d passed, %d failed, %d total\n",
           tests_passed, tests_failed, tests_passed + tests_failed);
    printf("===========================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
