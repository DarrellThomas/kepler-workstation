// Test multi-pass aerobraking campaign with real-gas thermodynamics
//
// Compile: g++ -std=c++11 -O2 -I../../PlanetaryModels -o test_aeropass test_aeropass.cpp
// Run:     ./test_aeropass

#include "aeropass_sim.hpp"
#include <cstdio>

int main()
{
    using namespace aeropass;

    Vehicle bucket = default_bucket();

    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  JUPITER ATMOSPHERIC PASS — REAL GAS vs PERFECT GAS        ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    // ── Single pass comparison ──
    printf("=== SINGLE PASS AT 350 km, 42 km/s ===\n");
    Orbit orb_350 = orbit_from_periapsis(350000.0, R_J * 5.0);
    PassResult single = integrate_pass(1, orb_350, bucket);

    printf("  Alt: %.0f km,  v = %.1f km/s,  Mach = %.0f\n",
           single.alt_peri_km, single.v_peri/1000, single.Mach_peri);
    printf("  Pass duration: %.0f s (%.1f min)\n", single.pass_duration, single.pass_duration/60);
    printf("  Arc length: %.0f km\n", single.arc_length/1000);
    printf("\n  REAL GAS results:\n");
    printf("    T_shock peak:      %.0f K\n", single.T_shock_peak);
    printf("    γ at shock:        %.3f\n", single.gamma_shock_peak);
    printf("    Density ratio:     %.1f (perfect gas: %.1f)\n",
           single.density_ratio_peak,
           (JUPITER_GAMMA+1)*single.Mach_peri*single.Mach_peri /
           ((JUPITER_GAMMA-1)*single.Mach_peri*single.Mach_peri + 2));
    printf("    q_stag peak:       %.1f kW/m²\n", single.q_stag_peak);
    printf("    q_hull peak:       %.2f kW/m² (%.0f%% leakage)\n",
           single.q_hull_peak, bucket.leak_fraction*100);
    printf("    Heat load:         %.1f MJ\n", single.Q_total);
    printf("    R_mag avg:         %.0f m\n", single.R_mag_avg);
    printf("    Mass collected:    %.1f kg\n", single.mass_collected);
    printf("    Δv from drag:      %.1f m/s\n", single.delta_v_drag);
    printf("    Orbit change:      %.1f Rj → %.1f Rj apoapsis\n",
           single.orbit_pre.r_apo / R_J, single.orbit_post.r_apo / R_J);

    // ── Multi-pass campaign: aerobrake from 42 km/s to 20 km/s ──
    printf("\n\n");
    CampaignResult campaign = run_campaign(
        350000.0,           // initial periapsis altitude (m)
        R_J * 5.0,          // initial apoapsis (5 Rj)
        bucket,
        20000.0,            // target velocity for collection (m/s)
        2000,               // max passes
        730.0               // max campaign days (2 years)
    );

    // ── Comparison: what if we DON'T aerobrake? ──
    printf("\n=== COMPARISON: NO AEROBRAKE (constant orbit) ===\n");
    printf("  At v=42 km/s constant:\n");
    printf("    Passes to fill tank: %d\n",
           (int)(bucket.tank_capacity / single.mass_collected));
    printf("    Mass per pass: %.1f kg\n", single.mass_collected);
    printf("    Hull heat flux: %.2f kW/m² per pass\n", single.q_hull_peak);
    printf("    ΔV for escape: ~5 km/s (all propulsive)\n\n");

    if (!campaign.passes.empty()) {
        PassResult &last = campaign.passes.back();
        printf("  With aerobraking to v=%.0f km/s:\n", campaign.v_final/1000);
        printf("    Mass per pass (at final v): %.1f kg\n", last.mass_collected);
        printf("    Hull heat flux: %.2f kW/m² per pass\n", last.q_hull_peak);
        double dv_saved = campaign.passes.front().orbit_pre.v_peri - campaign.v_final;
        printf("    ΔV saved by drag: %.1f km/s (free!)\n", dv_saved/1000);
        printf("    Remaining escape ΔV: ~%.1f km/s\n",
               5.0 - dv_saved/1000 > 0 ? 5.0 - dv_saved/1000 : 0.5);
    }

    return 0;
}
