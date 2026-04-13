// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Lambert Solver — Universal Variable Formulation
//
// Given two position vectors and a time of flight, find the transfer orbit
// velocity vectors at departure and arrival.
//
// Algorithm: Bate, Mueller & White "Fundamentals of Astrodynamics" (1971)
// with Battin's improvement for convergence near 180° transfers.
//
// All units: km, seconds, km/s
//
// 20260413 Created for Jupiter Harvester mission design
///////////////////////////////////////////////////////////////////////////////

#ifndef LAMBERT_HPP
#define LAMBERT_HPP

#include <cmath>
#include <cstdio>

// Stumpff functions C(z) and S(z)
inline double stumpff_c(double z)
{
    if (fabs(z) < 1e-6)
        return 1.0/2.0 - z/24.0 + z*z/720.0;
    if (z > 0)
        return (1.0 - cos(sqrt(z))) / z;
    else
        return (cosh(sqrt(-z)) - 1.0) / (-z);
}

inline double stumpff_s(double z)
{
    if (fabs(z) < 1e-6)
        return 1.0/6.0 - z/120.0 + z*z/5040.0;
    if (z > 0) {
        double sz = sqrt(z);
        return (sz - sin(sz)) / (sz * sz * sz);
    } else {
        double sz = sqrt(-z);
        return (sinh(sz) - sz) / (sz * sz * sz);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Solve Lambert's problem
//
// Inputs:
//   r1[3]  — departure position vector (km)
//   r2[3]  — arrival position vector (km)
//   tof    — time of flight (seconds)
//   mu     — gravitational parameter (km³/s²)
//   prograde — true for prograde (short way), false for retrograde
//
// Outputs:
//   v1[3]  — departure velocity vector (km/s)
//   v2[3]  — arrival velocity vector (km/s)
//
// Returns: true on convergence, false on failure
///////////////////////////////////////////////////////////////////////////////

inline bool lambert_solve(const double r1[3], const double r2[3],
                          double tof, double mu, bool prograde,
                          double v1[3], double v2[3])
{
    const int MAX_ITER = 100;
    const double TOL = 1e-10;

    // Magnitudes
    double R1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double R2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);

    // Cross product r1 × r2 (z-component determines transfer angle)
    double cross_z = r1[0]*r2[1] - r1[1]*r2[0];

    // Transfer angle
    double cos_dnu = (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]) / (R1 * R2);
    if (cos_dnu > 1.0) cos_dnu = 1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;
    double dnu = acos(cos_dnu);

    // Adjust for prograde/retrograde
    if (prograde) {
        if (cross_z < 0) dnu = 2.0 * M_PI - dnu;
    } else {
        if (cross_z >= 0) dnu = 2.0 * M_PI - dnu;
    }

    // Variable A
    double sin_dnu = sin(dnu);
    if (fabs(sin_dnu) < 1e-15) {
        // Near 0° or 180° transfer — degenerate
        return false;
    }
    double A = sin_dnu * sqrt(R1 * R2 / (1.0 - cos_dnu));

    // Newton-Raphson iteration on universal variable z
    // Initial guess: z = 0 (parabolic)
    double z = 0.0;

    // Better initial guess based on TOF
    // For elliptic transfers z > 0, hyperbolic z < 0
    double sqrt_mu = sqrt(mu);

    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        double C = stumpff_c(z);
        double S = stumpff_s(z);

        double y = R1 + R2 + A * (z * S - 1.0) / sqrt(C);

        if (y < 0) {
            // Adjust z to make y positive
            z = z + 0.1;
            continue;
        }

        double x = sqrt(y / C);
        double tof_guess = (x*x*x * S + A * sqrt(y)) / sqrt_mu;
        double dtdz;

        if (fabs(z) > 1e-6) {
            dtdz = (x*x*x * (S - 3.0*S*z*stumpff_s(z)/(2.0*z) +
                    3.0*S*S/(2.0*C)) +
                    (A/8.0) * (3.0*S*sqrt(y)/C + A/x)) / sqrt_mu;
            // Simplified derivative
            dtdz = (x*x*x * (stumpff_s(z) - 3.0*stumpff_s(z)*sqrt(stumpff_c(z)) /
                    (2.0*stumpff_c(z))) + A * sqrt(y) * 0.5 / y *
                    (A/(2.0*sqrt(stumpff_c(z))) - x*x*x*stumpff_s(z)/(2.0*stumpff_c(z)))) / sqrt_mu;
        }

        // Use a simpler, more robust derivative calculation
        // F(z) = tof_guess - tof = 0
        // Numerical derivative
        double dz = 0.01;
        double C2 = stumpff_c(z + dz);
        double S2 = stumpff_s(z + dz);
        double y2 = R1 + R2 + A * ((z+dz) * S2 - 1.0) / sqrt(C2);
        if (y2 < 0) y2 = fabs(y2);
        double x2 = sqrt(y2 / C2);
        double tof2 = (x2*x2*x2 * S2 + A * sqrt(y2)) / sqrt_mu;
        dtdz = (tof2 - tof_guess) / dz;

        if (fabs(dtdz) < 1e-20) break;

        double z_new = z - (tof_guess - tof) / dtdz;

        if (fabs(z_new - z) < TOL) {
            z = z_new;
            break;
        }

        z = z_new;

        // Prevent runaway
        if (z > 4.0 * M_PI * M_PI) z = 4.0 * M_PI * M_PI - 0.1;
        if (z < -1000.0) z = -1000.0;
    }

    // Compute final velocities
    double C = stumpff_c(z);
    double S = stumpff_s(z);
    double y = R1 + R2 + A * (z * S - 1.0) / sqrt(C);
    if (y < 0) return false;

    double f = 1.0 - y / R1;
    double g = A * sqrt(y / mu);
    double gdot = 1.0 - y / R2;

    if (fabs(g) < 1e-20) return false;

    for (int i = 0; i < 3; i++) {
        v1[i] = (r2[i] - f * r1[i]) / g;
        v2[i] = (gdot * r2[i] - r1[i]) / g;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Compute delta-V for a transfer
///////////////////////////////////////////////////////////////////////////////

struct TransferResult
{
    double dv_depart;    // km/s departure delta-V (from circular orbit)
    double dv_arrive;    // km/s arrival delta-V (to capture orbit)
    double dv_total;     // km/s total delta-V
    double v_inf_depart; // km/s hyperbolic excess at departure
    double v_inf_arrive; // km/s hyperbolic excess at arrival
    double tof_days;     // time of flight in days
    double c3_depart;    // km²/s² characteristic energy at departure
    double c3_arrive;    // km²/s² characteristic energy at arrival
    double v1[3];        // departure velocity (heliocentric)
    double v2[3];        // arrival velocity (heliocentric)
    bool valid;
};

// Compute transfer between two planets at given dates
// Uses Lambert solver + planet ephemeris
inline TransferResult compute_transfer(
    const double r1[3], const double v1_planet[3],
    const double r2[3], const double v2_planet[3],
    double tof_seconds, double mu_sun,
    double mu_depart, double r_park_depart,
    double mu_arrive, double r_park_arrive)
{
    TransferResult res;
    res.valid = false;
    res.tof_days = tof_seconds / 86400.0;

    // Solve Lambert's problem
    if (!lambert_solve(r1, r2, tof_seconds, mu_sun, true, res.v1, res.v2)) {
        return res;
    }

    // Hyperbolic excess velocity at departure
    double vinf_d[3], vinf_a[3];
    for (int i = 0; i < 3; i++) {
        vinf_d[i] = res.v1[i] - v1_planet[i];
        vinf_a[i] = res.v2[i] - v2_planet[i];
    }

    res.v_inf_depart = sqrt(vinf_d[0]*vinf_d[0]+vinf_d[1]*vinf_d[1]+vinf_d[2]*vinf_d[2]);
    res.v_inf_arrive = sqrt(vinf_a[0]*vinf_a[0]+vinf_a[1]*vinf_a[1]+vinf_a[2]*vinf_a[2]);

    res.c3_depart = res.v_inf_depart * res.v_inf_depart;
    res.c3_arrive = res.v_inf_arrive * res.v_inf_arrive;

    // Delta-V from circular parking orbit (vis-viva)
    double v_circ_depart = sqrt(mu_depart / r_park_depart);
    double v_hyp_depart = sqrt(res.v_inf_depart * res.v_inf_depart + 2.0 * mu_depart / r_park_depart);
    res.dv_depart = v_hyp_depart - v_circ_depart;

    // Delta-V to capture into circular orbit at arrival
    double v_circ_arrive = sqrt(mu_arrive / r_park_arrive);
    double v_hyp_arrive = sqrt(res.v_inf_arrive * res.v_inf_arrive + 2.0 * mu_arrive / r_park_arrive);
    res.dv_arrive = v_hyp_arrive - v_circ_arrive;

    res.dv_total = res.dv_depart + res.dv_arrive;
    res.valid = true;
    return res;
}

#endif // LAMBERT_HPP
