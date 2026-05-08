// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Jupiter Harvester — Multi-Moon Gravity Assist Analysis
//
// Evaluates using Callisto/Ganymede/Europa resonant gravity assists
// to lower the Bucket's orbit from interplanetary capture to the
// deep scoop orbit (350 km above Jupiter's 1-bar level).
//
// This is the "last mile" problem that Sun-Jupiter manifold tubes
// cannot solve — the moons operate at the right distance scale.
//
// Reference: KoLoMaRo (2011), Chapter 10.
//   Ross, Koon, Lo, Marsden (2003) — Multi-Moon Orbiter.
//   Figure 1.2.10: Tour of Jupiter's moons with dV ~ 22 m/s.
//
// Compile: g++ -std=c++11 -O2 -I../../Ephemeris/cr3bp -o multi_moon_analysis multi_moon_analysis.cpp
// Run:     ./multi_moon_analysis
///////////////////////////////////////////////////////////////////////////////

#include "cr3bp.hpp"
#include "patched_cr3bp.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// Bucket Mission Constants
///////////////////////////////////////////////////////////////////////////////

const double DRY_MASS      = 8000.0;
const double TANK_CAP      = 33000.0;
const double ISP           = 3000.0;
const double VE            = ISP * 9.80665e-3;
const double H2_FRAC       = 0.75;
const double HE_FRAC       = 0.25;
const double H2_TANK       = TANK_CAP * H2_FRAC;
const double HE_TANK       = TANK_CAP * HE_FRAC;

const double DV_EARTH_DEPART = 3.5;
const double DV_CAPTURE      = 0.5;
const double DV_JUP_DEPART   = 5.0;
const double DV_MARGIN       = 0.5;
const double DV_TOTAL_BASE   = DV_EARTH_DEPART + DV_CAPTURE + DV_JUP_DEPART + DV_MARGIN;

const double R_JUPITER     = 71492.0;
const double SCOOP_ALT     = 350.0;
const double R_SCOOP_PERI  = R_JUPITER + SCOOP_ALT;
const double R_SCOOP_APO   = R_JUPITER * 5.0;
const double A_SCOOP       = (R_SCOOP_PERI + R_SCOOP_APO) / 2.0;

double compute_payload(double dv_total)
{
    double mass_ratio = exp(dv_total / VE);
    double h2_payload = (H2_TANK - (DRY_MASS + HE_TANK) * (mass_ratio - 1.0)) / mass_ratio;
    return (h2_payload > 0) ? h2_payload : 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////////////
int main()
{
    printf("===================================================================\n");
    printf("  JUPITER HARVESTER — Multi-Moon Gravity Assist Analysis\n");
    printf("  Kepler Workstation | Patched Three-Body Approximation\n");
    printf("  Reference: KoLoMaRo (2011), Chapter 10\n");
    printf("===================================================================\n");

    double GM_J = BODY_JUPITER.GM;

    // =====================================================================
    // Section 1: Jovian Moon System Overview
    // =====================================================================
    printf("\n--- Section 1: Jovian Moon System ---\n\n");
    printf("  %-12s  %10s  %8s  %10s  %8s  %10s\n",
           "Moon", "a(km)", "a(Rj)", "Period(d)", "R(km)", "v_orb(km/s)");
    printf("  %-12s  %10s  %8s  %10s  %8s  %10s\n",
           "----", "------", "-----", "--------", "-----", "----------");

    const CelestialBody* moons[] = {&BODY_IO, &BODY_EUROPA, &BODY_GANYMEDE, &BODY_CALLISTO};
    for (int i = 0; i < 4; i++) {
        double v_orb = sqrt(GM_J / moons[i]->a);
        printf("  %-12s  %10.0f  %8.1f  %10.2f  %8.0f  %10.2f\n",
               moons[i]->name, moons[i]->a, moons[i]->a / R_JUPITER,
               moons[i]->period_days, moons[i]->radius, v_orb);
    }
    printf("\n  Scoop orbit: %.0f x %.0f km (%.1f x %.1f Rj)\n",
           R_SCOOP_PERI, R_SCOOP_APO, R_SCOOP_PERI/R_JUPITER, R_SCOOP_APO/R_JUPITER);

    // 4:2:1 resonance of Io:Europa:Ganymede
    printf("\n  Laplace resonance: Io:Europa:Ganymede = 4:2:1\n");
    printf("    Io period:       %.3f days\n", BODY_IO.period_days);
    printf("    Europa period:   %.3f days (2x Io: %.3f)\n",
           BODY_EUROPA.period_days, 2*BODY_IO.period_days);
    printf("    Ganymede period: %.3f days (4x Io: %.3f)\n",
           BODY_GANYMEDE.period_days, 4*BODY_IO.period_days);

    // =====================================================================
    // Section 2: Tisserand Analysis — What's Reachable?
    // =====================================================================
    printf("\n--- Section 2: Tisserand Reachability Analysis ---\n\n");

    // For each moon, compute the Tisserand parameter for an orbit
    // with apoapsis at the moon and various periapsis values
    printf("  Can gravity assists at each moon reach the scoop orbit?\n\n");
    printf("  %-12s  %10s  %12s  %10s  %10s\n",
           "Moon", "a_moon(Rj)", "T(r_peri=1Rj)", "min_peri(Rj)", "Reaches?");

    for (int i = 3; i >= 0; i--) {
        double a_m = moons[i]->a;
        double T_at_scoop = tisserand_from_apse(R_SCOOP_PERI, a_m, a_m);
        double min_peri = min_periapsis_from_tisserand(T_at_scoop, a_m, GM_J);

        printf("  %-12s  %10.1f  %12.4f  %10.2f  %10s\n",
               moons[i]->name, a_m/R_JUPITER, T_at_scoop,
               min_peri/R_JUPITER,
               (min_peri <= R_SCOOP_PERI * 1.5) ? "YES" : "no");
    }

    printf("\n  Interpretation: if min_peri <= 1.0 Rj, a single gravity-assist\n");
    printf("  chain at that moon can reach the scoop orbit from the moon's orbit.\n");

    // =====================================================================
    // Section 3: Multi-Moon Transfer Chain
    // =====================================================================
    printf("\n--- Section 3: Multi-Moon Gravity Assist Chain ---\n");

    // Define moon encounters (outer to inner)
    std::vector<MoonEncounter> encounter_chain;

    MoonEncounter callisto = {"Callisto", BODY_CALLISTO.a, BODY_CALLISTO.GM,
                              BODY_CALLISTO.radius,
                              sqrt(GM_J / BODY_CALLISTO.a)};
    MoonEncounter ganymede = {"Ganymede", BODY_GANYMEDE.a, BODY_GANYMEDE.GM,
                              BODY_GANYMEDE.radius,
                              sqrt(GM_J / BODY_GANYMEDE.a)};
    MoonEncounter europa   = {"Europa",   BODY_EUROPA.a,   BODY_EUROPA.GM,
                              BODY_EUROPA.radius,
                              sqrt(GM_J / BODY_EUROPA.a)};
    MoonEncounter io       = {"Io",       BODY_IO.a,       BODY_IO.GM,
                              BODY_IO.radius,
                              sqrt(GM_J / BODY_IO.a)};

    // Scenario A: Start from Callisto orbit (arriving from interplanetary)
    printf("\n  Scenario A: Callisto -> Ganymede -> Europa -> Scoop Orbit\n");
    encounter_chain = {callisto, ganymede, europa};
    MultiMoonTransfer transferA = compute_multi_moon_transfer(
        BODY_CALLISTO.a, R_SCOOP_PERI, GM_J, encounter_chain, 200.0);
    print_transfer(transferA);

    // Scenario B: Start from Ganymede orbit
    printf("\n  Scenario B: Ganymede -> Europa -> Scoop Orbit\n");
    encounter_chain = {ganymede, europa};
    MultiMoonTransfer transferB = compute_multi_moon_transfer(
        BODY_GANYMEDE.a, R_SCOOP_PERI, GM_J, encounter_chain, 200.0);
    print_transfer(transferB);

    // Scenario C: Full chain with Io
    printf("\n  Scenario C: Callisto -> Ganymede -> Europa -> Io -> Scoop Orbit\n");
    encounter_chain = {callisto, ganymede, europa, io};
    MultiMoonTransfer transferC = compute_multi_moon_transfer(
        BODY_CALLISTO.a, R_SCOOP_PERI, GM_J, encounter_chain, 200.0);
    print_transfer(transferC);

    // =====================================================================
    // Section 4: Bucket Delta-V Comparison
    // =====================================================================
    printf("\n--- Section 4: Bucket Delta-V Comparison ---\n");

    // Best scenario for the Bucket: use manifold to reach Callisto's orbit,
    // then multi-moon chain to scoop orbit
    // Capture: manifold + moon chain replaces the 0.5 km/s capture
    // Departure: reverse the chain (scoop -> Europa -> Ganymede -> Callisto -> manifold)

    MultiMoonTransfer& best = transferC; // use full chain

    // Capture: interplanetary to scoop orbit via moons
    double dv_capture_moons = best.total_dv;

    // Departure: scoop orbit back out via moons (roughly same cost, reversed)
    double dv_depart_moons = best.total_dv;

    // However: the current 5.0 km/s departure includes the Oberth burn at
    // scoop periapsis. The Oberth effect makes deep burns very efficient.
    // From scoop periapsis (71,842 km, v = 54.2 km/s), the Oberth effect
    // gives a huge velocity boost. The multi-moon chain avoids this by using
    // free gravity assists, but takes much longer.

    // Vis-viva: velocity at scoop periapsis in the scoop orbit
    double v_scoop_peri = sqrt(GM_J * (2.0/R_SCOOP_PERI - 1.0/A_SCOOP));

    // To escape Jupiter from scoop periapsis (parabolic):
    double v_escape_peri = sqrt(2.0 * GM_J / R_SCOOP_PERI);
    double dv_oberth_escape = v_escape_peri - v_scoop_peri;

    printf("\n  Key velocities at scoop periapsis (%.0f km):\n", R_SCOOP_PERI);
    printf("    v_scoop (in orbit):  %.2f km/s\n", v_scoop_peri);
    printf("    v_escape (parabolic): %.2f km/s\n", v_escape_peri);
    printf("    Oberth escape dV:    %.2f km/s\n", dv_oberth_escape);
    printf("    Current budget:      %.2f km/s (includes v_inf for Earth return)\n", DV_JUP_DEPART);

    double dv_manifold_total = DV_EARTH_DEPART + dv_capture_moons + dv_depart_moons + DV_MARGIN;

    printf("\n  %-28s %10s %10s %10s\n", "", "Baseline", "Manifold+", "Savings");
    printf("  %-28s %10s %10s %10s\n",   "", "",         "Moons",     "");
    printf("  %-28s %10s %10s %10s\n",   "", "--------", "--------",  "-------");
    printf("  %-28s %9.2f  %9.2f  %+8.2f km/s\n",
           "Earth departure", DV_EARTH_DEPART, DV_EARTH_DEPART, 0.0);
    printf("  %-28s %9.2f  %9.3f  %+8.2f km/s\n",
           "Jupiter capture (moon GAs)", DV_CAPTURE, dv_capture_moons,
           DV_CAPTURE - dv_capture_moons);
    printf("  %-28s %9.2f  %9.3f  %+8.2f km/s\n",
           "Jupiter departure (moon GAs)", DV_JUP_DEPART, dv_depart_moons,
           DV_JUP_DEPART - dv_depart_moons);
    printf("  %-28s %9.2f  %9.2f  %9.2f km/s\n",
           "Margin", DV_MARGIN, DV_MARGIN, 0.0);
    printf("  %-28s %9.2f  %9.3f  %+8.2f km/s\n",
           "TOTAL", DV_TOTAL_BASE, dv_manifold_total,
           DV_TOTAL_BASE - dv_manifold_total);

    printf("\n  Time cost: %.0f days (%.1f years) per capture+departure chain\n",
           2 * best.total_tof_days, 2 * best.total_tof_days / 365.25);

    // Mass budget
    double h2_base = compute_payload(DV_TOTAL_BASE);
    double h2_moon = compute_payload(dv_manifold_total);
    double total_base = h2_base + HE_TANK;
    double total_moon = h2_moon + HE_TANK;

    printf("\n  Mass budget:\n");
    printf("    %-28s %8.1f t  %8.1f t  %+7.1f t\n",
           "H2 delivered", h2_base/1000, h2_moon/1000, (h2_moon-h2_base)/1000);
    printf("    %-28s %8.1f t  %8.1f t  %7.1f t\n",
           "He delivered", HE_TANK/1000, HE_TANK/1000, 0.0);
    printf("    %-28s %8.1f t  %8.1f t  %+7.1f t\n",
           "TOTAL delivered", total_base/1000, total_moon/1000, (total_moon-total_base)/1000);
    printf("    %-28s %8.2fx  %8.2fx\n",
           "Delivery ratio (x dry)", total_base/DRY_MASS, total_moon/DRY_MASS);

    // Fleet impact
    double loops_per_year = 10.0;
    printf("\n  Fleet impact (100 buckets, ~10 deliveries/year):\n");
    printf("    Baseline: %.0f tonnes/year\n", total_base/1000 * loops_per_year);
    printf("    Moon GAs: %.0f tonnes/year\n", total_moon/1000 * loops_per_year);
    printf("    Gain:     %+.0f tonnes/year\n", (total_moon-total_base)/1000 * loops_per_year);

    // =====================================================================
    // Section 5: The Tradeoff — Delta-V vs Time
    // =====================================================================
    printf("\n--- Section 5: The Tradeoff ---\n");
    printf("\n  The multi-moon gravity assist chain trades TIME for FUEL.\n");
    printf("  KoLoMaRo Eq. 10.4.6: dV_tot ~ 340 - 0.60 * TOF (m/s, days)\n\n");

    printf("  %-20s %8s %8s %10s\n", "Strategy", "dV(km/s)", "TOF(yr)", "Payload(t)");
    printf("  %-20s %8s %8s %10s\n", "--------", "-------", "------", "---------");

    // Current baseline
    printf("  %-20s %8.2f %8.1f %10.1f\n",
           "Baseline (Oberth)", DV_TOTAL_BASE,
           0.0, total_base/1000);

    // Moon GA chain
    double tof_moon_yr = 2 * best.total_tof_days / 365.25;
    printf("  %-20s %8.3f %8.1f %10.1f\n",
           "Full moon chain", dv_manifold_total,
           tof_moon_yr, total_moon/1000);

    // Hybrid: use Oberth for departure but moons for capture
    double dv_hybrid = DV_EARTH_DEPART + dv_capture_moons + DV_JUP_DEPART + DV_MARGIN;
    double h2_hybrid = compute_payload(dv_hybrid);
    double total_hybrid = h2_hybrid + HE_TANK;
    printf("  %-20s %8.3f %8.1f %10.1f\n",
           "Hybrid (moon+Oberth)", dv_hybrid,
           best.total_tof_days / 365.25, total_hybrid/1000);

    printf("\n===================================================================\n");
    printf("  CONCLUSION\n");
    printf("===================================================================\n");
    printf("\n  THE OBERTH EFFECT WINS.\n\n");
    printf("  Multi-moon gravity assists at 6-26 Rj cannot compete with a\n");
    printf("  single Oberth burn at scoop periapsis (1.0 Rj, v=54 km/s).\n");
    printf("  The deep gravity well is the Bucket's ADVANTAGE, not a cost:\n\n");
    printf("    - Oberth escape from 1 Rj: %.1f km/s for parabolic escape\n", dv_oberth_escape);
    printf("    - Moon GA chain to reach 1 Rj: %.1f km/s (propulsive burns)\n", best.total_dv);
    printf("    - Current baseline (VEEJ + Oberth): %.1f km/s total\n\n", DV_TOTAL_BASE);
    printf("  The KoLoMaRo 22 m/s result (Ch 10, Fig 1.2.10) applies to\n");
    printf("  transfers BETWEEN moon orbits (Ganymede L1 -> Europa L2),\n");
    printf("  not to pumping periapsis from 15 Rj down to 1 Rj. The Bucket\n");
    printf("  needs to dive deep into Jupiter's atmosphere — a fundamentally\n");
    printf("  different problem than the Multi-Moon Orbiter science mission.\n\n");
    printf("  DESIGN VALIDATION: The current VEEJ + magnetic-drag capture\n");
    printf("  + Oberth departure strategy is near-optimal for the Bucket's\n");
    printf("  unique mission profile. The Hall thruster (Isp=3000s) burning\n");
    printf("  scooped H2 at scoop periapsis is the most mass-efficient way\n");
    printf("  to escape Jupiter's deep gravity well.\n");
    printf("===================================================================\n");

    return 0;
}
