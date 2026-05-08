// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Invariant Manifold Computation for CR3BP Periodic Orbits
//
// Given a periodic orbit (e.g., halo orbit), computes the stable and
// unstable invariant manifolds by:
//   1. Computing the monodromy matrix (STM over one full period)
//   2. Finding the stable/unstable eigenvectors (real eigenvalues)
//   3. Perturbing the orbit state along the eigenvectors at N points
//   4. Propagating forward (unstable) or backward (stable) in time
//
// The resulting manifold "tubes" are the dynamical highways that
// control transport between realms in the CR3BP.
//
// Reference: KoLoMaRo (2011), Chapters 7 and 9.
//
// Compiles standalone: g++ -std=c++11 -O2
//
// 20260508 Created for Kepler Workstation CR3BP module
///////////////////////////////////////////////////////////////////////////////

#ifndef MANIFOLD_HPP
#define MANIFOLD_HPP

#include "cr3bp.hpp"
#include "halo_orbit.hpp"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
// Monodromy Matrix and Eigendecomposition
///////////////////////////////////////////////////////////////////////////////

// Compute monodromy matrix: STM over one full period of a periodic orbit
inline void compute_monodromy(double mu, const CR3BPState& ic, double period,
                              double monodromy[6][6])
{
    double phi[6][6];
    cr3bp_propagate_stm(mu, ic, period, phi, 0.0005);
    memcpy(monodromy, phi, sizeof(double)*36);
}

// Power iteration to find dominant eigenvalue/eigenvector of a 6x6 matrix
// Returns the eigenvalue (real) and the corresponding eigenvector
inline double power_iteration(const double M[6][6], double evec[6],
                              int max_iter = 1000, double tol = 1e-12)
{
    // Initial guess
    for (int i = 0; i < 6; i++) evec[i] = 1.0;
    double eigenvalue = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        // y = M * x
        double y[6] = {0};
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                y[i] += M[i][j] * evec[j];

        // Find component with largest magnitude
        double max_val = 0;
        for (int i = 0; i < 6; i++)
            if (fabs(y[i]) > max_val) max_val = fabs(y[i]);

        if (max_val < 1e-30) break;

        // Normalize
        double ev_new = max_val;
        for (int i = 0; i < 6; i++) evec[i] = y[i] / max_val;

        // Check convergence
        if (fabs(ev_new - eigenvalue) < tol * fabs(ev_new)) {
            eigenvalue = ev_new;
            break;
        }
        eigenvalue = ev_new;
    }

    // Determine sign: check if M*v ≈ lambda*v
    double Mv[6] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            Mv[i] += M[i][j] * evec[j];

    // Check if eigenvalue should be negative
    double dot_pos = 0, dot_neg = 0;
    for (int i = 0; i < 6; i++) {
        dot_pos += (Mv[i] - eigenvalue * evec[i]) * (Mv[i] - eigenvalue * evec[i]);
        dot_neg += (Mv[i] + eigenvalue * evec[i]) * (Mv[i] + eigenvalue * evec[i]);
    }
    if (dot_neg < dot_pos) eigenvalue = -eigenvalue;

    return eigenvalue;
}

// Inverse power iteration to find smallest eigenvalue
// (stable eigenvalue = 1/unstable eigenvalue for symplectic matrices)
inline double inverse_power_iteration(const double M[6][6], double evec[6],
                                      int max_iter = 1000, double tol = 1e-12)
{
    // For a symplectic matrix (like the monodromy), if lambda is an eigenvalue,
    // then 1/lambda is also an eigenvalue. The stable eigenvector corresponds
    // to the eigenvalue with |lambda| < 1.
    //
    // We use the fact that for the monodromy matrix of a halo orbit,
    // the eigenvalues come in pairs: (lambda, 1/lambda) and two pairs on
    // the unit circle. We can find the stable direction by solving
    // M * v = (1/lambda_u) * v, i.e., M^{-1} * v = lambda_u * v
    //
    // Rather than inverting M, we use inverse iteration: solve M*y = x

    // Simple approach: solve (M - sigma*I)y = x with sigma near 0
    // But for simplicity, since monodromy is symplectic, if we know lambda_u,
    // we know lambda_s = 1/lambda_u and can use shifted power iteration.

    // For now: use power iteration on M^T (transpose) to find left eigenvector,
    // which gives us the stable direction since M^T and M^{-T} share eigenvectors
    // for symplectic matrices.

    // Actually, simplest robust approach: get unstable eigenvector from power
    // iteration, then get stable eigenvector from power iteration on the
    // inverse. For symplectic M, M^{-1} = J * M^T * J^{-1} where J is the
    // symplectic form. But let's just do power iteration on M^T to get the
    // left eigenvector of the smallest eigenvalue.

    // More practical: just use power iteration to find largest eigenvalue of
    // M, deflate, repeat. Or: compute M^{-1} via symplectic property.

    // Symplectic inverse: for a 2n×2n symplectic matrix M:
    // M^{-1} = J_{2n}^T * M^T * J_{2n}
    // where J_{2n} = [0 I_n; -I_n 0]

    // Build M^{-1} using symplectic property
    double Minv[6][6];

    // J = [0 I; -I 0] for 3+3 system
    // M^{-1} = J^T * M^T * J
    // J^T = -J for the standard symplectic form

    // Direct computation:
    // If M = [A B; C D] (3x3 blocks), then
    // M^{-1} = [D^T  -B^T; -C^T  A^T]
    // This is a well-known property of symplectic matrices.

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Minv[i][j]     =  M[j+3][i+3]; // D^T
            Minv[i][j+3]   = -M[j][i+3];   // -B^T
            Minv[i+3][j]   = -M[j+3][i];   // -C^T
            Minv[i+3][j+3] =  M[j][i];     // A^T
        }
    }

    // Power iteration on M^{-1} gives eigenvector for largest eigenvalue
    // of M^{-1} = 1/lambda_min of M = lambda_stable
    double lambda_inv = power_iteration(Minv, evec, max_iter, tol);

    // The eigenvalue of M corresponding to this eigenvector is 1/lambda_inv
    return 1.0 / lambda_inv;
}

///////////////////////////////////////////////////////////////////////////////
// Manifold Tube Computation
///////////////////////////////////////////////////////////////////////////////

struct ManifoldTube
{
    std::vector<std::vector<CR3BPState>> trajectories;
    // trajectories[i] = the i-th manifold trajectory (time series)
    // i indexes the departure point along the periodic orbit
    // Each trajectory is propagated for the specified duration

    bool is_stable;   // true = stable (W^s), false = unstable (W^u)
    int  branch;      // +1 or -1 (interior/exterior branch)
};

// Compute invariant manifold tubes for a periodic orbit
//
// Inputs:
//   mu          — mass parameter
//   halo        — converged halo orbit
//   n_points    — number of points along the orbit to launch manifold arcs
//   t_manifold  — propagation time for each manifold arc (nondimensional)
//   n_record    — number of points to record per arc
//   epsilon     — perturbation magnitude along eigenvector
//   stable      — true for stable manifold (W^s), false for unstable (W^u)
//   branch      — +1 or -1 (two branches of each manifold)
//
// Returns: ManifoldTube with n_points trajectories

inline ManifoldTube compute_manifold(
    double mu, const HaloOrbit& halo,
    int n_points, double t_manifold, int n_record = 200,
    double epsilon = 1e-6, bool stable = false, int branch = 1)
{
    ManifoldTube tube;
    tube.is_stable = stable;
    tube.branch = branch;
    tube.trajectories.resize(n_points);

    // Step 1: Compute monodromy matrix
    double monodromy[6][6];
    compute_monodromy(mu, halo.ic, halo.period, monodromy);

    // Step 2: Find unstable eigenvector (largest eigenvalue)
    double v_unstable[6], v_stable[6];
    double lambda_u = power_iteration(monodromy, v_unstable);
    double lambda_s = inverse_power_iteration(monodromy, v_stable);

    // Choose which eigenvector to use
    double* v_eig = stable ? v_stable : v_unstable;

    // Normalize eigenvector (position part only, for consistent perturbation size)
    double pos_norm = sqrt(v_eig[0]*v_eig[0] + v_eig[1]*v_eig[1] + v_eig[2]*v_eig[2]);
    if (pos_norm > 1e-20) {
        for (int i = 0; i < 6; i++) v_eig[i] /= pos_norm;
    }

    // Step 3: For each point along the orbit, compute local eigenvector
    // via STM transport: v(t) = Phi(t,0) * v(0)
    // Then perturb and propagate

    double dt_orbit = halo.period / n_points;

    for (int i = 0; i < n_points; i++) {
        double t_depart = i * dt_orbit;

        // Propagate orbit + STM from t=0 to t_depart
        double s42[42];
        cr3bp_stm_identity(s42, halo.ic);

        CR3BPIntegrator integ;
        if (t_depart > 1e-14) {
            integ.propagate(cr3bp_eom_stm, mu, s42, 42, 0.0, t_depart, 0.001);
        }

        // State on the periodic orbit at this point
        CR3BPState orbit_state;
        orbit_state.from_array(s42);

        // Transport eigenvector: v(t) = Phi(t,0) * v(0)
        double phi[6][6];
        cr3bp_stm_extract(s42, phi);

        double v_local[6] = {0};
        for (int r = 0; r < 6; r++)
            for (int c = 0; c < 6; c++)
                v_local[r] += phi[r][c] * v_eig[c];

        // Re-normalize position part
        double pn = sqrt(v_local[0]*v_local[0] + v_local[1]*v_local[1] + v_local[2]*v_local[2]);
        if (pn > 1e-20)
            for (int r = 0; r < 6; r++) v_local[r] /= pn;

        // Perturb state along eigenvector
        CR3BPState perturbed;
        double sign = (double)branch;
        perturbed.x  = orbit_state.x  + sign * epsilon * v_local[0];
        perturbed.y  = orbit_state.y  + sign * epsilon * v_local[1];
        perturbed.z  = orbit_state.z  + sign * epsilon * v_local[2];
        perturbed.vx = orbit_state.vx + sign * epsilon * v_local[3];
        perturbed.vy = orbit_state.vy + sign * epsilon * v_local[4];
        perturbed.vz = orbit_state.vz + sign * epsilon * v_local[5];

        // Propagate: forward for unstable manifold, backward for stable
        double t_prop = stable ? -t_manifold : t_manifold;
        tube.trajectories[i] = cr3bp_propagate_record(mu, perturbed, t_prop, n_record, 0.001);
    }

    return tube;
}

///////////////////////////////////////////////////////////////////////////////
// Convenience: compute all four manifold branches
// (stable interior, stable exterior, unstable interior, unstable exterior)
///////////////////////////////////////////////////////////////////////////////

struct ManifoldSet
{
    ManifoldTube W_s_plus;   // stable, interior branch
    ManifoldTube W_s_minus;  // stable, exterior branch
    ManifoldTube W_u_plus;   // unstable, interior branch
    ManifoldTube W_u_minus;  // unstable, exterior branch
};

inline ManifoldSet compute_all_manifolds(
    double mu, const HaloOrbit& halo,
    int n_points = 50, double t_manifold = 6.0, int n_record = 200,
    double epsilon = 1e-6)
{
    ManifoldSet ms;
    ms.W_s_plus  = compute_manifold(mu, halo, n_points, t_manifold, n_record, epsilon, true,  +1);
    ms.W_s_minus = compute_manifold(mu, halo, n_points, t_manifold, n_record, epsilon, true,  -1);
    ms.W_u_plus  = compute_manifold(mu, halo, n_points, t_manifold, n_record, epsilon, false, +1);
    ms.W_u_minus = compute_manifold(mu, halo, n_points, t_manifold, n_record, epsilon, false, -1);
    return ms;
}

#endif // MANIFOLD_HPP
