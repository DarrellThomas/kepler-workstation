// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Patched Three-Body Approximation (P3BA) for Multi-Moon Transfers
//
// Implements the Concentric Circular Model (CCM) for the Jovian system:
// Jupiter is central, moons orbit in circles. The four-body problem
// (Jupiter + 2 moons + spacecraft) is decomposed into two CR3BP systems,
// patched at a switching orbit where one moon's influence gives way
// to another.
//
// Key tools:
//   1. Tisserand parameter analysis for resonant gravity assists
//   2. Energy-based delta-V estimation (KoLoMaRo Eq. 10.4.5)
//   3. Gravity assist chains between moon pairs
//   4. v-infinity leveraging for orbit lowering/raising
//
// Reference: KoLoMaRo (2011), Chapters 5 and 10.
//   Ross, Koon, Lo, Marsden (2003, 2004) — Multi-Moon Orbiter.
//
// Compiles standalone: g++ -std=c++11 -O2
//
// 20260508 Created for Kepler Workstation
///////////////////////////////////////////////////////////////////////////////

#ifndef PATCHED_CR3BP_HPP
#define PATCHED_CR3BP_HPP

#include "cr3bp.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
// Tisserand Parameter
//
// The Tisserand parameter T is a quasi-integral of motion that is
// approximately conserved across a gravity assist with a moon.
// In the CR3BP, it relates to the Jacobi constant.
//
// T = a_M/a + 2*cos(i)*sqrt((a/a_M)*(1-e^2))
//
// where a_M is the moon's semi-major axis, and (a, e, i) are the
// spacecraft's osculating orbital elements around Jupiter.
//
// For planar orbits (i=0): T = a_M/a + 2*sqrt((a/a_M)*(1-e^2))
//
// Key property: a gravity assist by the moon changes (a, e) but
// preserves T (approximately). So we can trace curves of constant T
// in the (periapsis, apoapsis) plane to find which orbits are
// connected by a single gravity assist.
///////////////////////////////////////////////////////////////////////////////

// Compute Tisserand parameter for a planar orbit around Jupiter
// a: semi-major axis (km), e: eccentricity
// a_moon: moon's orbital radius (km)
inline double tisserand(double a, double e, double a_moon)
{
    if (a <= 0 || a_moon <= 0) return 0;
    return a_moon / a + 2.0 * sqrt(a / a_moon * (1.0 - e*e));
}

// Compute Tisserand parameter from periapsis and apoapsis
inline double tisserand_from_apse(double r_peri, double r_apo, double a_moon)
{
    double a = (r_peri + r_apo) / 2.0;
    double e = (r_apo - r_peri) / (r_apo + r_peri);
    return tisserand(a, e, a_moon);
}

// v-infinity at a moon encounter, derived from Tisserand parameter
// v_inf^2 = 3 - T (in nondimensional units where v_moon = 1)
// In dimensional: v_inf = v_moon * sqrt(3 - T)
// v_moon should be in consistent units (km/s)
inline double v_infinity_from_tisserand(double T, double v_moon)
{
    double val = 3.0 - T;
    if (val < 0) return 0;
    return v_moon * sqrt(val);
}

// Moon orbital velocity (circular, around Jupiter)
inline double moon_orbital_velocity(double a_moon, double GM_central)
{
    return sqrt(GM_central / a_moon);
}

///////////////////////////////////////////////////////////////////////////////
// Resonant Orbit Analysis
//
// A p:q resonance with a moon means the spacecraft completes p orbits
// while the moon completes q orbits. This gives:
//   a_resonance = a_moon * (q/p)^(2/3)
//
// The spacecraft encounters the moon every q spacecraft orbits.
// Resonant orbits are the stepping stones for gravity-assist chains.
///////////////////////////////////////////////////////////////////////////////

struct ResonantOrbit
{
    int p, q;            // resonance ratio (spacecraft:moon)
    double a;            // semi-major axis (km)
    double period_days;  // orbital period (days)
    double e_for_peri;   // eccentricity if periapsis is at target radius
    double r_peri;       // periapsis (km)
    double r_apo;        // apoapsis (km)
    double T;            // Tisserand parameter
    double v_inf_kms;    // v-infinity at moon (km/s)
};

// Compute resonant orbit parameters
// target_peri: desired periapsis (km), 0 = compute from resonance only
inline ResonantOrbit compute_resonance(int p, int q, double a_moon,
                                       double GM_central,
                                       double target_peri = 0)
{
    ResonantOrbit res;
    res.p = p;
    res.q = q;
    res.a = a_moon * pow((double)q / p, 2.0/3.0);
    res.period_days = 2.0 * CR3BP_PI * sqrt(res.a * res.a * res.a / GM_central) / 86400.0;

    if (target_peri > 0 && target_peri < 2.0 * res.a) {
        res.r_peri = target_peri;
        res.r_apo = 2.0 * res.a - target_peri;
        res.e_for_peri = 1.0 - target_peri / res.a;
    } else {
        // Default: periapsis at moon distance
        res.e_for_peri = 0.0;
        res.r_peri = res.a;
        res.r_apo = res.a;
    }

    res.T = tisserand(res.a, res.e_for_peri, a_moon);

    double v_moon = moon_orbital_velocity(a_moon, GM_central);
    res.v_inf_kms = v_infinity_from_tisserand(res.T, v_moon / 1000.0);

    return res;
}

///////////////////////////////////////////////////////////////////////////////
// Gravity Assist Delta-V Estimation
//
// From KoLoMaRo Eq. 10.4.5: the delta-V to change between two
// three-body energies E_a and E_b near a moon at distance r_M:
//
// dV ≈ |sqrt(2(mu/r_M + 3/2 + E_a)) - sqrt(2(mu/r_M + 3/2 + E_b))|
//
// where mu is the CR3BP mass parameter and r_M is the distance from
// the moon (nondimensional).
//
// For our purposes, we use a simpler Keplerian approach:
// dV at periapsis = |v_periapsis(orbit_b) - v_periapsis(orbit_a)|
// using the vis-viva equation.
///////////////////////////////////////////////////////////////////////////////

// Vis-viva velocity at radius r in orbit with semi-major axis a
inline double vis_viva_v(double r, double a, double GM)
{
    return sqrt(GM * (2.0/r - 1.0/a));
}

// Delta-V to change from orbit (a1, e1) to orbit (a2, e2) at a common
// apse point at radius r_common.
inline double dv_at_apse(double r_common, double a1, double a2, double GM)
{
    double v1 = vis_viva_v(r_common, a1, GM);
    double v2 = vis_viva_v(r_common, a2, GM);
    return fabs(v2 - v1);
}

///////////////////////////////////////////////////////////////////////////////
// Gravity Assist Chain Analysis
//
// Given a sequence of moons (outer to inner), compute the delta-V
// cost of using resonant gravity assists to lower the periapsis
// from an initial orbit to a target orbit.
//
// The approach:
// 1. At each moon, the Tisserand parameter T is conserved
// 2. The spacecraft enters with some v_inf at the moon
// 3. The gravity assist bends the velocity vector, changing (a, e)
//    but preserving T
// 4. Over multiple encounters, the periapsis can be lowered
// 5. The delta-V cost is the sum of small correction burns needed
//    to achieve the right resonance geometry
///////////////////////////////////////////////////////////////////////////////

struct MoonEncounter
{
    const char* moon_name;
    double a_moon;       // moon orbital radius (km)
    double GM_moon;      // moon GM (km^3/s^2)
    double R_moon;       // moon radius (km)
    double v_moon;       // moon orbital velocity (km/s)
};

struct GravityAssistLeg
{
    MoonEncounter moon;
    double r_peri_before; // periapsis before GA chain at this moon (km)
    double r_peri_after;  // periapsis after GA chain (km)
    double r_apo_before;  // apoapsis before (km)
    double r_apo_after;   // apoapsis after (km)
    double T;             // Tisserand parameter at this moon
    double v_inf;         // v-infinity at moon (km/s)
    double dv;            // delta-V for this leg (km/s)
    int n_flybys;         // estimated number of flybys
    double tof_days;      // estimated time of flight (days)
};

// Compute what periapsis can be achieved from a given apoapsis
// using a single gravity assist at a moon (Tisserand constraint).
// Given: spacecraft approaches moon with apoapsis = a_moon (flyby at apoapsis),
// and current periapsis r_peri_current.
// After GA: periapsis can be lowered to some minimum determined by T conservation.
// The minimum periapsis for a given T occurs when the apoapsis equals the moon's orbit.
inline double min_periapsis_from_tisserand(double T, double a_moon, double GM_central)
{
    // For an orbit with apoapsis = a_moon and Tisserand parameter T:
    // T = a_moon/a + 2*sqrt(a/a_moon * (1-e^2))
    // With r_apo = a(1+e) = a_moon and r_peri = a(1-e):
    // a = a_moon/(1+e), so a/a_moon = 1/(1+e)
    // 1-e^2 = (1-e)(1+e), so a/a_moon*(1-e^2) = (1-e)
    // T = (1+e) + 2*sqrt(1-e)
    // Let u = sqrt(1-e): e = 1-u^2, T = (2-u^2) + 2u = 2 + 2u - u^2
    // u^2 - 2u + (T-2) = 0
    // u = 1 - sqrt(1-(T-2)) = 1 - sqrt(3-T)
    // e = 1 - u^2 = 1 - (1-sqrt(3-T))^2

    double val = 3.0 - T;
    if (val < 0) return a_moon; // no flyby possible
    double u = 1.0 - sqrt(val);
    if (u <= 0) return 0; // can reach Jupiter center
    double e = 1.0 - u * u;
    if (e < 0) e = 0;
    if (e >= 1.0) return 0;
    double a = a_moon / (1.0 + e);
    double r_peri = a * (1.0 - e);
    return r_peri;
}

// Compute the Tisserand parameter for an orbit that just touches a moon's
// orbit at apoapsis, coming from a given periapsis
inline double tisserand_at_moon_from_peri(double r_peri, double a_moon)
{
    double r_apo = a_moon;
    return tisserand_from_apse(r_peri, r_apo, a_moon);
}

// Estimate number of gravity assist flybys needed to lower periapsis
// from r_peri_start to r_peri_end using a moon at a_moon.
// Each flyby rotates the v-infinity vector, changing periapsis by ~2*v_inf*sin(delta/2)
// where delta is the flyby deflection angle.
inline int estimate_flyby_count(double r_peri_start, double r_peri_end,
                                double a_moon, double GM_central,
                                double GM_moon, double R_moon,
                                double flyby_alt = 100.0)
{
    // v-infinity at moon
    double T = tisserand_from_apse(r_peri_start, a_moon, a_moon);
    double v_moon = sqrt(GM_central / a_moon);
    double v_inf = v_moon * sqrt(fmax(0, 3.0 - T));

    if (v_inf < 1e-6) return 1;

    // Maximum deflection per flyby (at minimum flyby distance)
    double r_flyby = R_moon + flyby_alt;
    double sin_half_delta = 1.0 / (1.0 + r_flyby * v_inf * v_inf / GM_moon);
    double delta_max = 2.0 * asin(fmin(sin_half_delta, 1.0));

    // Change in periapsis per flyby (rough estimate)
    // Each flyby changes v by ~2*v_inf*sin(delta/2)
    double dv_per_flyby = 2.0 * v_inf * sin(delta_max / 2.0);
    double v_peri_start = sqrt(GM_central * (2.0/r_peri_start - 2.0/(r_peri_start + a_moon)));
    double v_peri_end = sqrt(GM_central * (2.0/r_peri_end - 2.0/(r_peri_end + a_moon)));
    double total_dv_needed = fabs(v_peri_end - v_peri_start);

    if (dv_per_flyby < 1e-6) return 100;
    int n = (int)ceil(total_dv_needed / dv_per_flyby);
    return fmax(n, 1);
}

///////////////////////////////////////////////////////////////////////////////
// Multi-Moon Transfer Chain
//
// Computes a sequence of gravity assists through the Jovian system:
//   Callisto → Ganymede → Europa → target orbit
//
// At each moon, the spacecraft performs multiple resonant gravity assists
// to lower (or raise) its orbit until it can be captured by the next moon.
///////////////////////////////////////////////////////////////////////////////

struct MultiMoonTransfer
{
    std::vector<GravityAssistLeg> legs;
    double total_dv;          // total deterministic delta-V (km/s)
    double total_tof_days;    // total time of flight (days)
    double r_peri_initial;    // starting periapsis (km)
    double r_peri_final;      // ending periapsis (km)
};

inline MultiMoonTransfer compute_multi_moon_transfer(
    double r_peri_start,      // starting periapsis (km) — from interplanetary approach
    double r_peri_target,     // target periapsis (km) — scoop orbit
    double GM_central,        // Jupiter GM (km^3/s^2)
    const std::vector<MoonEncounter>& moons,  // ordered outer to inner
    double flyby_alt = 200.0) // minimum flyby altitude (km)
{
    MultiMoonTransfer transfer;
    transfer.r_peri_initial = r_peri_start;
    transfer.r_peri_final = r_peri_target;
    transfer.total_dv = 0;
    transfer.total_tof_days = 0;

    double r_peri = r_peri_start;

    for (size_t i = 0; i < moons.size(); i++) {
        const MoonEncounter& moon = moons[i];
        GravityAssistLeg leg;
        leg.moon = moon;

        // Current orbit: periapsis at r_peri, apoapsis at moon's orbit
        // (we need to reach the moon first)
        double r_apo = moon.a_moon;

        // If periapsis is already below this moon, skip
        if (r_peri >= r_apo) {
            // Need to lower from above — compute arrival orbit
            // The spacecraft arrives from outside with apoapsis beyond this moon
            r_apo = r_peri; // orbit doesn't reach this moon yet
            // Use a small burn to target this moon
            // For now, assume previous moon delivered us to an orbit
            // that reaches this moon's orbit
            if (i == 0) {
                // First moon (outermost) — spacecraft arrives from interplanetary
                // Assume apoapsis is at or beyond this moon
                r_apo = fmax(r_peri, moon.a_moon);
            }
        }

        leg.r_peri_before = r_peri;
        leg.r_apo_before = r_apo;

        // Tisserand parameter at this moon
        leg.T = tisserand_from_apse(r_peri, r_apo, moon.a_moon);

        // v-infinity at this moon (km/s, since GM is km^3/s^2 and a is km)
        double v_moon = sqrt(GM_central / moon.a_moon);
        leg.v_inf = v_moon * sqrt(fmax(0, 3.0 - leg.T));

        // What periapsis can we achieve with GA at this moon?
        double min_peri = min_periapsis_from_tisserand(leg.T, moon.a_moon, GM_central);

        // Target: either the next moon's orbit, or the final target
        double target_peri;
        if (i + 1 < moons.size()) {
            target_peri = moons[i+1].a_moon; // lower to next moon's orbit
        } else {
            target_peri = r_peri_target; // final target
        }

        // Can we reach the target with Tisserand-conserving GAs?
        if (min_peri <= target_peri) {
            // Yes — GAs alone can lower periapsis to target
            leg.r_peri_after = target_peri;
            leg.dv = 0.0; // free (ballistic)

            // But we need a small deterministic burn to set up the right geometry
            // KoLoMaRo notes ~22 m/s for Callisto-Ganymede-Europa tour
            // Use the empirical formula: dV ≈ 10-50 m/s per moon
            leg.dv = 0.020; // 20 m/s per moon (conservative estimate)
        } else {
            // No — need a propulsive burn to change T
            // Compute delta-V at periapsis to change orbit
            double a_before = (leg.r_peri_before + moon.a_moon) / 2.0;
            double a_after  = (target_peri + moon.a_moon) / 2.0;
            leg.r_peri_after = target_peri;
            leg.dv = dv_at_apse(moon.a_moon, a_before, a_after, GM_central);
        }

        leg.r_apo_after = moon.a_moon;

        // Estimate number of flybys and time
        leg.n_flybys = estimate_flyby_count(leg.r_peri_before, leg.r_peri_after,
                                            moon.a_moon, GM_central,
                                            moon.GM_moon, moon.R_moon, flyby_alt);

        // Time: each flyby takes roughly the moon's orbital period
        double T_moon_days = 2.0 * CR3BP_PI * sqrt(moon.a_moon * moon.a_moon * moon.a_moon / GM_central) / 86400.0;
        leg.tof_days = leg.n_flybys * T_moon_days;

        transfer.legs.push_back(leg);
        transfer.total_dv += leg.dv;
        transfer.total_tof_days += leg.tof_days;

        r_peri = leg.r_peri_after;
    }

    return transfer;
}

///////////////////////////////////////////////////////////////////////////////
// Print transfer summary
///////////////////////////////////////////////////////////////////////////////

inline void print_transfer(const MultiMoonTransfer& t, double R_jupiter = 71492.0)
{
    printf("\n  Multi-Moon Transfer Chain:\n");
    printf("  %-12s  %10s  %10s  %10s  %8s  %6s  %8s\n",
           "Moon", "Peri(Rj)", "->Peri(Rj)", "v_inf(km/s)", "dV(km/s)", "GAs", "TOF(d)");
    printf("  %-12s  %10s  %10s  %10s  %8s  %6s  %8s\n",
           "----", "--------", "--------", "----------", "-------", "---", "------");

    for (const auto& leg : t.legs) {
        printf("  %-12s  %10.1f  %10.1f  %10.2f  %8.3f  %6d  %8.0f\n",
               leg.moon.moon_name,
               leg.r_peri_before / R_jupiter,
               leg.r_peri_after / R_jupiter,
               leg.v_inf,
               leg.dv,
               leg.n_flybys,
               leg.tof_days);
    }

    printf("  %-12s  %10s  %10s  %10s  %8.3f  %6s  %8.0f\n",
           "TOTAL", "", "", "",
           t.total_dv, "",
           t.total_tof_days);
}

#endif // PATCHED_CR3BP_HPP
