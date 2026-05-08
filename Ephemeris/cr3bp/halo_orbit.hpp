// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Halo Orbit Computation for the CR3BP
//
// Two-stage approach following KoLoMaRo Chapter 6:
//   1. Richardson 3rd-order analytical approximation (initial guess)
//   2. Single-shooting differential correction using STM (converge to
//      true periodic orbit)
//
// Reference:
//   Richardson, D.L. (1980). "Analytic Construction of Periodic Orbits
//     About the Collinear Points." Celestial Mechanics, 22, 241-253.
//   Koon, Lo, Marsden, Ross (2011). Chapter 6.
//   Howell, K.C. (1984). "Three-Dimensional, Periodic, 'Halo' Orbits."
//     Celestial Mechanics, 32, 53-71.
//
// Conventions:
//   - Coordinates are in the normalized coordinate system centered at
//     L1 or L2, scaled by gamma (distance from Lagrange point to m2).
//     See KoLoMaRo Eq. 6.3 ff.
//   - Northern halo (class I): Az > 0, m = 1
//   - Southern halo (class II): Az < 0 (mirror across xy-plane), m = 3
//
// Compiles standalone: g++ -std=c++11 -O2
//
// 20260508 Created for Kepler Workstation CR3BP module
///////////////////////////////////////////////////////////////////////////////

#ifndef HALO_ORBIT_HPP
#define HALO_ORBIT_HPP

#include "cr3bp.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// Richardson 3rd-Order Halo Orbit Approximation
//
// Generates initial conditions for a halo orbit near L1 or L2 given
// the out-of-plane amplitude Az (nondimensional, in CR3BP units).
//
// KoLoMaRo Section 6.6, pp. 154-158
///////////////////////////////////////////////////////////////////////////////

struct RichardsonCoeffs
{
    // System parameters
    double mu;
    double gamma;    // distance from Li to m2 (nondimensional)
    int    point;    // 1 for L1, 2 for L2
    double x_Li;    // x-coordinate of Lagrange point

    // Legendre expansion coefficients c_n (Eq. below 6.3.1)
    double c2, c3, c4;

    // Eigenvalues
    double lambda;   // saddle eigenvalue
    double omega_p;  // in-plane frequency
    double omega_v;  // out-of-plane frequency = sqrt(c2)
    double kappa;    // y/x amplitude ratio

    // Richardson coefficients
    double a21, a22, a23, a24;
    double b21, b22;
    double d21;
    double a31, a32;
    double b31, b32;
    double d31, d32;
    double s1, s2;
    double l1, l2;
    double Delta;    // correction term Omega_p^2 - c2
};

// Compute Richardson coefficients for L1 or L2
// point: 1 for L1, 2 for L2
inline RichardsonCoeffs richardson_coefficients(double mu, int point)
{
    RichardsonCoeffs R;
    R.mu = mu;
    R.point = point;

    // Compute Lagrange points
    LagrangePoints lp = cr3bp_lagrange_points(mu);
    if (point == 1) {
        R.x_Li = lp.L1[0];
        R.gamma = (1.0 - mu) - R.x_Li; // distance from L1 to m2
    } else {
        R.x_Li = lp.L2[0];
        R.gamma = R.x_Li - (1.0 - mu); // distance from L2 to m2
    }

    double g = R.gamma;
    double g3 = g * g * g;

    // Sign convention: upper sign for L1, lower for L2
    // c_n = (1/gamma^3) * (±^n * mu + (-1)^n * (1-mu)*gamma^{n+1} / (1∓gamma)^{n+1})
    double mu1 = 1.0 - mu;
    double sign1 = (point == 1) ? 1.0 : -1.0;  // ± for L1/L2

    // c_n formula from KoLoMaRo below Eq. 6.3.1:
    // c_n = (1/gamma^3) * ((±1)^n * mu + (-1)^n * (1-mu) * gamma^{n+1} / (1 ∓ gamma)^{n+1})
    for (int n = 2; n <= 4; n++) {
        double sign_n = (n % 2 == 0) ? 1.0 : -1.0; // (-1)^n
        double pm_n;
        if (point == 1)
            pm_n = pow(1.0, n);  // (+1)^n = 1 always for L1
        else
            pm_n = pow(-1.0, n); // (-1)^n for L2

        double denom = (point == 1) ? (1.0 - g) : (1.0 + g);
        double cn = (1.0 / g3) * (pm_n * mu + sign_n * mu1 * pow(g, n+1) / pow(denom, n+1));

        if (n == 2) R.c2 = cn;
        else if (n == 3) R.c3 = cn;
        else R.c4 = cn;
    }

    // Eigenvalues from linearized equations (KoLoMaRo Eq. below 6.4.1)
    // lambda^2 = (c2 - 2 + sqrt(9*c2^2 - 8*c2)) / 2
    // omega_p^2 = (2 - c2 + sqrt(9*c2^2 - 8*c2)) / 2
    // omega_v^2 = c2
    double disc = sqrt(9.0*R.c2*R.c2 - 8.0*R.c2);
    R.lambda = sqrt((R.c2 - 2.0 + disc) / 2.0);
    R.omega_p = sqrt((2.0 - R.c2 + disc) / 2.0);
    R.omega_v = sqrt(R.c2);

    // kappa = (omega_p^2 + 1 + 2*c2) / (2*omega_p)  (Eq. 6.4.3)
    R.kappa = (R.omega_p*R.omega_p + 1.0 + 2.0*R.c2) / (2.0*R.omega_p);

    // Delta = omega_p^2 - c2  (correction term)
    R.Delta = R.omega_p*R.omega_p - R.c2;

    double wp = R.omega_p;
    double k = R.kappa;
    double c2 = R.c2, c3 = R.c3, c4 = R.c4;

    // d1, d2 denominators used in several coefficients
    double d1 = 3.0*wp*wp / k * (k*(6.0*wp*wp - 1.0) - 2.0*wp);
    double d2 = 8.0*wp*wp / k * (k*(11.0*wp*wp - 1.0) - 2.0*wp);

    // Richardson coefficients (KoLoMaRo pp. 156-157)
    R.a21 = 3.0*c3*(k*k - 2.0) / (4.0*(1.0 + 2.0*c2));
    R.a22 = 3.0*c3 / (4.0*(1.0 + 2.0*c2));
    R.a23 = -3.0*c3*wp / (4.0*k*d1) * (3.0*k*k*k*wp - 6.0*k*(k - wp) + 4.0);
    R.a24 = -3.0*c3*wp / (4.0*k*d1) * (2.0 + 3.0*k*wp);

    R.b21 = -3.0*c3*wp / (2.0*d1) * (3.0*k*wp - 4.0);
    R.b22 = -3.0*c3*wp / d1;

    R.d21 = -c3 / (2.0*wp*wp);

    R.a31 = -9.0*wp / (4.0*d2) * (4.0*c3*(k*R.a23 - R.b21) + k*c4*(4.0 + k*k))
             + (9.0*wp*wp + 1.0 - c2) / (2.0*d2) *
               (3.0*c3*(2.0*R.a23 - k*R.b21) + c4*(2.0 + 3.0*k*k));

    R.a32 = -9.0*wp / (4.0*d2) * (4.0*c3*(k*R.a24 - R.b22) + k*c4)
             - 3.0/(2.0*d2) * (9.0*wp*wp + 1.0 - c2) *
               (c3*(k*R.b22 + R.d21 - 2.0*R.a24) - c4);

    R.b31 = 3.0/(8.0*d2) * 8.0*wp * (3.0*c3*(k*R.b21 - 2.0*R.a23) - c4*(2.0 + 3.0*k*k))
             + 3.0/(8.0*d2) * (9.0*wp*wp + 1.0 + 2.0*c2) *
               (4.0*c3*(k*R.a23 - R.b21) + k*c4*(4.0 + k*k));

    R.b32 = 9.0*wp/d2 * (c3*(k*R.b22 + R.d21 - 2.0*R.a24) - c4)
             + 3.0*(9.0*wp*wp + 1.0 + 2.0*c2) / (8.0*d2) *
               (4.0*c3*(k*R.a24 - R.b22) + k*c4);

    R.d31 = 3.0 / (64.0*wp*wp) * (4.0*c3*R.a24 + c4);
    R.d32 = 3.0 / (64.0*wp*wp) * (4.0*c3*(R.a23 - R.d21) + c4*(4.0 + k*k));

    // s1, s2: frequency correction coefficients (Eq. 6.6.3)
    R.s1 = 1.0 / (2.0*wp*(wp*(1.0 + k*k) - 2.0*k))
            * (1.5*c3*(2.0*R.a21*(k*k - 2.0) - R.a23*(k*k + 2.0) - 2.0*k*R.b21)
               - 3.0/8.0*c4*(3.0*k*k*k*k - 8.0*k*k + 8.0));

    R.s2 = 1.0 / (2.0*wp*(wp*(1.0 + k*k) - 2.0*k))
            * (1.5*c3*(2.0*R.a22*(k*k - 2.0) + R.a24*(k*k + 2.0) + 2.0*k*R.b22 + 5.0*R.d21)
               + 3.0/8.0*c4*(12.0 - k*k));

    // l1, l2: amplitude constraint coefficients (Eq. 6.6.4)
    R.l1 = -1.5*c3*(2.0*R.a21 + R.a23 + 5.0*R.d21)
            - 3.0/8.0*c4*(12.0 - k*k) + 2.0*wp*wp*R.s1;

    R.l2 = 1.5*c3*(R.a24 - 2.0*R.a22)
           + 9.0/8.0*c4 + 2.0*wp*wp*R.s2;

    return R;
}

///////////////////////////////////////////////////////////////////////////////
// Generate initial conditions from Richardson approximation
//
// Inputs:
//   R    — Richardson coefficients (from richardson_coefficients)
//   Az   — out-of-plane amplitude (nondimensional, in CR3BP units)
//          Positive for northern halo (class I)
//   phi  — phase angle (typically 0 for standard halo)
//
// Output: CR3BPState at t=0 (the xz-plane crossing)
///////////////////////////////////////////////////////////////////////////////

inline CR3BPState richardson_initial_guess(const RichardsonCoeffs& R,
                                           double Az, double phi = 0.0)
{
    double g = R.gamma;
    double wp = R.omega_p;
    double k = R.kappa;

    // Determine halo class
    // m=1: northern (delta_m = 2-1 = 1, so delta_m > 0)
    // m=3: southern (delta_m = 2-3 = -1, so delta_m < 0)
    int m = (Az >= 0) ? 1 : 3;
    double delta_m = 2.0 - m;  // +1 for northern, -1 for southern
    double Az_abs = fabs(Az);

    // Convert Az from CR3BP units to normalized (L-point centered, gamma-scaled)
    double Az_n = Az_abs / g;

    // Ax from amplitude constraint (Eq. 6.6.4):
    // l1*Ax^2 + l2*Az^2 + Delta = 0
    // Ax = sqrt((-l2*Az_n^2 - Delta) / l1)
    double Ax_sq = (-R.l2*Az_n*Az_n - R.Delta) / R.l1;
    if (Ax_sq < 0) {
        // Az too small for halo orbit — below bifurcation
        // Use minimum amplitude
        Ax_sq = fabs(R.Delta / R.l1);
        Az_n = 0.0;
    }
    double Ax = sqrt(Ax_sq);

    // Frequency correction (Eq. 6.6.3): nu = 1 + s1*Ax^2 + s2*Az^2
    double nu = 1.0 + R.s1*Ax*Ax + R.s2*Az_n*Az_n;

    // tau1 = omega_p * tau + phi, evaluated at tau = 0
    double tau1 = phi;

    // 3rd-order solution in normalized coords (KoLoMaRo p. 156)
    // At tau = 0 (initial time, xz-plane crossing for phi=0):
    double cos_t1 = cos(tau1);
    double sin_t1 = sin(tau1);
    double cos_2t1 = cos(2.0*tau1);
    double sin_2t1 = sin(2.0*tau1);
    double cos_3t1 = cos(3.0*tau1);
    double sin_3t1 = sin(3.0*tau1);

    // Position in normalized (gamma-scaled, L-point centered) coords
    double x_n = R.a21*Ax*Ax + R.a22*Az_n*Az_n - Ax*cos_t1
                 + (R.a23*Ax*Ax - R.a24*Az_n*Az_n)*cos_2t1
                 + (R.a31*Ax*Ax*Ax - R.a32*Ax*Az_n*Az_n)*cos_3t1;

    double y_n = k*Ax*sin_t1
                 + (R.b21*Ax*Ax - R.b22*Az_n*Az_n)*sin_2t1
                 + (R.b31*Ax*Ax*Ax - R.b32*Ax*Az_n*Az_n)*sin_3t1;

    double z_n = delta_m*Az_n*cos_t1
                 + delta_m*R.d21*Ax*Az_n*(cos_2t1 - 3.0)
                 + delta_m*(R.d32*Az_n*Ax*Ax - R.d31*Az_n*Az_n*Az_n)*cos_3t1;

    // Velocity in normalized coords (d/dtau derivatives, times omega_p*nu)
    double wpnu = wp * nu;
    double vx_n = wpnu * (Ax*sin_t1
                  - 2.0*(R.a23*Ax*Ax - R.a24*Az_n*Az_n)*sin_2t1
                  - 3.0*(R.a31*Ax*Ax*Ax - R.a32*Ax*Az_n*Az_n)*sin_3t1);

    double vy_n = wpnu * (k*Ax*cos_t1
                  + 2.0*(R.b21*Ax*Ax - R.b22*Az_n*Az_n)*cos_2t1
                  + 3.0*(R.b31*Ax*Ax*Ax - R.b32*Ax*Az_n*Az_n)*cos_3t1);

    double vz_n = wpnu * (-delta_m*Az_n*sin_t1
                  - 2.0*delta_m*R.d21*Ax*Az_n*sin_2t1
                  - 3.0*delta_m*(R.d32*Az_n*Ax*Ax - R.d31*Az_n*Az_n*Az_n)*sin_3t1);

    // Convert from normalized coords to CR3BP rotating frame
    // x_cr3bp = x_Li ± gamma * x_n  (+ for L2, - for L1)
    // y_cr3bp = gamma * y_n
    // z_cr3bp = gamma * z_n
    double sign = (R.point == 1) ? -1.0 : 1.0;

    CR3BPState s;
    s.x  = R.x_Li + sign * g * x_n;
    s.y  = g * y_n;
    s.z  = g * z_n;
    s.vx = sign * g * vx_n;
    s.vy = g * vy_n;
    s.vz = g * vz_n;

    return s;
}

///////////////////////////////////////////////////////////////////////////////
// Differential Correction for Halo Orbits
//
// Single-shooting method: exploit xz-plane symmetry to find periodic orbits.
//
// Initial state: (x0, 0, z0, 0, vy0, 0)  — on xz-plane, perpendicular crossing
// Target: propagate to next y=0 crossing, require vx_f = vz_f = 0
//
// Use STM to compute corrections to (x0, z0, vy0).
//
// KoLoMaRo Section 6.7, pp. 159-160
///////////////////////////////////////////////////////////////////////////////

struct HaloOrbit
{
    CR3BPState ic;       // initial conditions (at y=0 crossing)
    double period;       // full orbital period (nondimensional)
    double jacobi;       // Jacobi constant
    double Az;           // out-of-plane amplitude
    double Ax;           // in-plane amplitude (x)
    double Ay;           // in-plane amplitude (y)
    bool converged;
    int iterations;
};

// Event detection: find y=0 crossing during propagation
// Returns time of crossing (refined by bisection)
inline double find_y_crossing(double mu, const double s0[42], double t_guess,
                              double* s_out)
{
    // Propagate slightly past expected half-period
    double s[42];
    memcpy(s, s0, 42*sizeof(double));

    CR3BPIntegrator integ;
    double dt_step = 0.001;

    // First: propagate with small steps until y changes sign
    double t = 0.0;
    double y_prev = s[1];
    double t_cross = t_guess;
    bool found = false;

    // Coarse propagation
    int n_steps = (int)(t_guess * 1.5 / dt_step) + 100;
    double s_prev[42];
    memcpy(s_prev, s, 42*sizeof(double));

    for (int i = 0; i < n_steps; i++) {
        memcpy(s_prev, s, 42*sizeof(double));
        double t_prev = t;
        t = integ.propagate(cr3bp_eom_stm, mu, s, 42, t, t + dt_step, dt_step*0.1);

        // Check for y sign change (skip first few steps)
        if (i > 5 && y_prev * s[1] < 0) {
            // Bisection refinement
            double ta = t_prev, tb = t;
            double sa[42], sb[42];
            memcpy(sa, s_prev, 42*sizeof(double));
            memcpy(sb, s, 42*sizeof(double));

            for (int j = 0; j < 60; j++) {
                double tm = (ta + tb) / 2.0;
                double sm[42];
                memcpy(sm, sa, 42*sizeof(double));
                integ.propagate(cr3bp_eom_stm, mu, sm, 42, ta, tm, dt_step*0.01);

                if (sm[1] * sa[1] < 0) {
                    tb = tm;
                    memcpy(sb, sm, 42*sizeof(double));
                } else {
                    ta = tm;
                    memcpy(sa, sm, 42*sizeof(double));
                }

                if (fabs(tb - ta) < 1e-14) break;
            }

            t_cross = (ta + tb) / 2.0;
            memcpy(s_out, sa, 42*sizeof(double));

            // Final refinement to exact y=0
            double sm[42];
            memcpy(sm, s0, 42*sizeof(double));
            integ.propagate(cr3bp_eom_stm, mu, sm, 42, 0.0, t_cross, dt_step*0.01);
            memcpy(s_out, sm, 42*sizeof(double));

            found = true;
            break;
        }
        y_prev = s[1];
    }

    if (!found) {
        // Fallback: just propagate to t_guess
        memcpy(s_out, s0, 42*sizeof(double));
        integ.propagate(cr3bp_eom_stm, mu, s_out, 42, 0.0, t_guess, dt_step*0.01);
        t_cross = t_guess;
    }

    return t_cross;
}

// Differential correction iteration
// Corrects (x0, z0, vy0) to satisfy (vx_f, vz_f) = 0 at y=0 crossing
inline HaloOrbit halo_differential_correction(double mu, CR3BPState ic_guess,
                                              double period_guess,
                                              int max_iter = 30,
                                              double tol = 1e-11)
{
    HaloOrbit halo;
    halo.converged = false;
    halo.iterations = 0;

    // Enforce symmetry: y=0, vx=0, vz=0 at initial state
    CR3BPState ic = ic_guess;
    ic.y = 0.0;
    ic.vx = 0.0;
    ic.vz = 0.0;

    double half_period = period_guess / 2.0;

    for (int iter = 0; iter < max_iter; iter++) {
        halo.iterations = iter + 1;

        // Set up 42-state with STM = I
        double s0[42];
        cr3bp_stm_identity(s0, ic);

        // Propagate to y=0 crossing (half period)
        double sf[42];
        double t_half = find_y_crossing(mu, s0, half_period, sf);

        // Check convergence: vx_f and vz_f should be zero
        double vx_f = sf[3];
        double vz_f = sf[5];

        if (fabs(vx_f) < tol && fabs(vz_f) < tol) {
            halo.converged = true;
            halo.ic = ic;
            halo.period = 2.0 * t_half;
            halo.jacobi = cr3bp_jacobi(ic, mu);

            // Compute amplitudes by propagating full orbit
            std::vector<CR3BPState> traj = cr3bp_propagate_record(mu, ic, halo.period, 500);
            double xmin = 1e20, xmax = -1e20;
            double ymin = 1e20, ymax = -1e20;
            double zmin = 1e20, zmax = -1e20;
            for (const auto& st : traj) {
                if (st.x < xmin) xmin = st.x;
                if (st.x > xmax) xmax = st.x;
                if (st.y < ymin) ymin = st.y;
                if (st.y > ymax) ymax = st.y;
                if (st.z < zmin) zmin = st.z;
                if (st.z > zmax) zmax = st.z;
            }
            halo.Ax = (xmax - xmin) / 2.0;
            halo.Ay = (ymax - ymin) / 2.0;
            halo.Az = (zmax - zmin) / 2.0;

            return halo;
        }

        // Extract STM at half-period
        double phi[6][6];
        cr3bp_stm_extract(sf, phi);

        // Acceleration at final state (needed for correction)
        double dsdt[6];
        cr3bp_eom(mu, sf, dsdt);
        double ax_f = dsdt[3]; // vx_dot at final state
        double az_f = dsdt[5]; // vz_dot at final state
        double vy_f = sf[4];   // vy at final state (nonzero, used for y-crossing)

        // We have 3 free variables: (x0, z0, vy0) = (ic.x, ic.z, ic.vy)
        // 2 constraints: vx_f = 0, vz_f = 0
        // Since we have more unknowns than constraints, we fix one variable.
        // Standard approach: fix z0, correct (x0, vy0)
        //
        // From STM and chain rule at y=0 crossing:
        // [dvx_f]   [phi(4,1) - ax_f*phi(2,1)/vy_f, phi(4,5) - ax_f*phi(2,5)/vy_f] [dx0 ]
        // [dvz_f] = [phi(6,1) - az_f*phi(2,1)/vy_f, phi(6,5) - az_f*phi(2,5)/vy_f] [dvy0]

        double M[2][2];
        M[0][0] = phi[3][0] - ax_f * phi[1][0] / vy_f;
        M[0][1] = phi[3][4] - ax_f * phi[1][4] / vy_f;
        M[1][0] = phi[5][0] - az_f * phi[1][0] / vy_f;
        M[1][1] = phi[5][4] - az_f * phi[1][4] / vy_f;

        // Solve 2x2: M * [dx0; dvy0] = -[vx_f; vz_f]
        double det = M[0][0]*M[1][1] - M[0][1]*M[1][0];
        if (fabs(det) < 1e-20) {
            // Singular — try fixing x0 instead, correct (z0, vy0)
            M[0][0] = phi[3][2] - ax_f * phi[1][2] / vy_f;
            M[0][1] = phi[3][4] - ax_f * phi[1][4] / vy_f;
            M[1][0] = phi[5][2] - az_f * phi[1][2] / vy_f;
            M[1][1] = phi[5][4] - az_f * phi[1][4] / vy_f;

            det = M[0][0]*M[1][1] - M[0][1]*M[1][0];
            if (fabs(det) < 1e-20) break; // truly singular

            double dz0  = (-vx_f * M[1][1] + vz_f * M[0][1]) / det;
            double dvy0 = (-vz_f * M[0][0] + vx_f * M[1][0]) / det;

            ic.z  += dz0;
            ic.vy += dvy0;
        } else {
            double dx0  = (-vx_f * M[1][1] + vz_f * M[0][1]) / det;
            double dvy0 = (-vz_f * M[0][0] + vx_f * M[1][0]) / det;

            ic.x  += dx0;
            ic.vy += dvy0;
        }

        // Update half-period estimate
        half_period = t_half;
    }

    // Did not converge — return best guess
    halo.ic = ic;
    halo.period = 2.0 * half_period;
    halo.jacobi = cr3bp_jacobi(ic, mu);
    halo.Az = 0;
    halo.Ax = 0;
    halo.Ay = 0;
    return halo;
}

///////////////////////////////////////////////////////////////////////////////
// High-level: compute a halo orbit given system and amplitude
//
// Inputs:
//   mu    — mass parameter
//   point — 1 for L1, 2 for L2
//   Az    — out-of-plane amplitude (nondimensional, CR3BP units)
//           Use system.L_km to convert: Az_nd = Az_km / L_km
//   northern — true for class I (northern), false for class II (southern)
//
// Returns: HaloOrbit with initial conditions, period, and amplitudes
///////////////////////////////////////////////////////////////////////////////

inline HaloOrbit compute_halo_orbit(double mu, int point, double Az,
                                    bool northern = true,
                                    int max_iter = 30, double tol = 1e-11)
{
    // Get Richardson coefficients
    RichardsonCoeffs R = richardson_coefficients(mu, point);

    // Richardson initial guess
    double Az_signed = northern ? fabs(Az) : -fabs(Az);
    CR3BPState ic_guess = richardson_initial_guess(R, Az_signed);

    // Estimated period: T = 2*pi / (omega_p * nu)
    double Az_n = fabs(Az) / R.gamma;
    double Ax_sq = (-R.l2*Az_n*Az_n - R.Delta) / R.l1;
    if (Ax_sq < 0) Ax_sq = fabs(R.Delta / R.l1);
    double Ax = sqrt(Ax_sq);
    double nu = 1.0 + R.s1*Ax*Ax + R.s2*Az_n*Az_n;
    double period_guess = 2.0 * CR3BP_PI / (R.omega_p * nu);

    // Differential correction
    return halo_differential_correction(mu, ic_guess, period_guess, max_iter, tol);
}

///////////////////////////////////////////////////////////////////////////////
// Halo orbit family continuation
//
// Generate a family of halo orbits by stepping in Az amplitude.
// Each orbit uses the previous converged orbit as seed for the next.
///////////////////////////////////////////////////////////////////////////////

inline std::vector<HaloOrbit> compute_halo_family(
    double mu, int point, double Az_start, double Az_end,
    int n_orbits, bool northern = true)
{
    std::vector<HaloOrbit> family;
    family.reserve(n_orbits);

    double dAz = (Az_end - Az_start) / (n_orbits - 1);

    HaloOrbit prev;
    prev.converged = false;

    for (int i = 0; i < n_orbits; i++) {
        double Az = Az_start + i * dAz;
        HaloOrbit halo;

        if (prev.converged && i > 0) {
            // Use previous orbit as seed with slight perturbation
            CR3BPState ic = prev.ic;
            // Scale z proportionally to new Az
            double scale = Az / (Az - dAz);
            ic.z *= scale;
            halo = halo_differential_correction(mu, ic, prev.period);
        } else {
            // Use Richardson from scratch
            halo = compute_halo_orbit(mu, point, Az, northern);
        }

        family.push_back(halo);
        if (halo.converged) prev = halo;
    }

    return family;
}

#endif // HALO_ORBIT_HPP
