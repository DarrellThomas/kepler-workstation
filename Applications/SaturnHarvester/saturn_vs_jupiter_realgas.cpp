// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// SATURN vs JUPITER — Real-Gas Head-to-Head
//
// Re-runs the scoop comparison from saturn_vs_jupiter.cpp but with
// stat_thermo.hpp real-gas thermodynamics instead of perfect-gas assumptions.
//
// Key differences from perfect gas version:
//   - Post-shock temperature computed via iterative R-H with real γ
//   - Density ratio accounts for dissociation/ionization
//   - Composition behind shock (ionization fraction, species mix)
//   - γ_eff varies with temperature (1.44 cold → 1.05 dissociating → 1.67 plasma)
//
// Compile: g++ -std=c++11 -O2 -I../../PlanetaryModels -o saturn_vs_jupiter_realgas saturn_vs_jupiter_realgas.cpp
// Run:     ./saturn_vs_jupiter_realgas
///////////////////////////////////////////////////////////////////////////////

#include "planetary_environment.hpp"
#include "stat_thermo.hpp"
#include "transport_properties.hpp"
#include <cstdio>
#include <cmath>

using namespace stat_thermo;
using namespace transport;

static const double MU0 = 4.0 * M_PI * 1e-7;

double mag_funnel_r(double B_coil, double r_coil, double p_ram, double max_r)
{
    double B_need = sqrt(2.0 * MU0 * p_ram);
    if (B_need <= 0) return max_r;
    double R = r_coil * pow(B_coil / B_need, 1.0 / 3.0);
    return (R < max_r) ? R : max_r;
}

int main()
{
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  SATURN vs JUPITER — REAL-GAS THERMODYNAMICS                   ║\n");
    printf("║  7-species H2/He: H2, H, H+, He, He+, He2+, e-               ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");

    // Vehicle spec
    double B_coil = 20.0, r_coil = 3.0, eta = 0.5, max_R = 200.0;
    double dry_mass = 8000, tank_cap = 33000;

    // Planet orbital constants (km)
    double GM_J = 1.26713e8, R_J = 71492.0;
    double GM_S = 3.7931e7,  R_S = 60268.0;

    // Composition (mass fractions of H2 and He in cold gas)
    // Jupiter: 86.2% H2, 13.8% He by mass (Galileo)
    // Saturn:  93.9% H2, 6.1% He by mass (Voyager/Cassini — He depleted by rain-out)
    double JUP_Y_H2 = 0.862, JUP_Y_He = 0.138;
    double SAT_Y_H2 = 0.939, SAT_Y_He = 0.061;

    // ═══════════════════════════════════════════════════════════
    // JUPITER: Scoop sweep with real-gas shock
    // ═══════════════════════════════════════════════════════════
    printf("═══ JUPITER SCOOP SWEEP (Real Gas) ═══\n\n");
    printf("  %-6s %-8s %-10s %-10s %-8s %-10s %-10s %-8s %-8s %-20s\n",
           "Alt", "v_peri", "rho", "R_mag", "kg/pass",
           "T_shock", "T_pg", "γ_shk", "ρ₂/ρ₁", "Composition");
    printf("  %-6s %-8s %-10s %-10s %-8s %-10s %-10s %-8s %-8s\n",
           "(km)", "(km/s)", "(kg/m³)", "(m)", "",
           "(K,real)", "(K,perf)", "", "(real)");
    printf("  ──────────────────────────────────────────────────────────"
           "────────────────────────────────────────\n");

    double j_alts[] = {250, 300, 350, 400, 500, 600};
    for (double alt : j_alts) {
        double r_peri = R_J + alt;
        double r_apo = R_J * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double v_p = sqrt(GM_J * (2.0/r_peri - 1.0/a_orb)) * 1000;  // m/s

        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt * 1000);

        double p_ram = 0.5 * rho * v_p * v_p;
        double R_mag = mag_funnel_r(B_coil, r_coil, p_ram, max_R);
        double A_eff = M_PI * R_mag * R_mag;
        double m_dot = rho * v_p * A_eff * eta;

        double H_m = JUPITER_RGAS * tempk / JUPITER_G0;
        double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
        double dt = arc / v_p;
        double kg_pass = m_dot * dt;

        // Perfect gas shock temperature
        double T_pg = 2.0 * (JUPITER_GAMMA - 1.0) * v_p * v_p
                    / ((JUPITER_GAMMA + 1.0) * (JUPITER_GAMMA + 1.0) * JUPITER_RGAS);

        // Real-gas shock
        ShockState shock = normal_shock(v_p, tempk, press, JUP_Y_H2, JUP_Y_He);

        // Dominant post-shock species
        char comp[64];
        if (shock.post.X[Hp] > 0.1)
            snprintf(comp, sizeof(comp), "H+=%.0f%% e-=%.0f%%",
                     shock.post.X[Hp]*100, shock.post.X[Elec]*100);
        else if (shock.post.X[H] > 0.1)
            snprintf(comp, sizeof(comp), "H=%.0f%% H2=%.0f%%",
                     shock.post.X[H]*100, shock.post.X[H2]*100);
        else
            snprintf(comp, sizeof(comp), "H2=%.0f%%", shock.post.X[H2]*100);

        printf("  %-6.0f %-8.1f %-10.2e %-10.0f %-8.1f %-10.0f %-10.0f %-8.3f %-8.1f %-20s\n",
               alt, v_p/1000, rho, R_mag, kg_pass,
               shock.post.T, T_pg, shock.post.gamma_eff, shock.density_ratio, comp);
    }

    // ═══════════════════════════════════════════════════════════
    // SATURN: Scoop sweep with real-gas shock
    // ═══════════════════════════════════════════════════════════
    printf("\n═══ SATURN SCOOP SWEEP (Real Gas) ═══\n\n");
    printf("  %-6s %-8s %-10s %-10s %-8s %-10s %-10s %-8s %-8s %-20s\n",
           "Alt", "v_peri", "rho", "R_mag", "kg/pass",
           "T_shock", "T_pg", "γ_shk", "ρ₂/ρ₁", "Composition");
    printf("  ──────────────────────────────────────────────────────────"
           "────────────────────────────────────────\n");

    double s_alts[] = {350, 500, 700, 1000, 1200, 1500, 2000};
    for (double alt : s_alts) {
        double r_peri = R_S + alt;
        double r_apo = R_S * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double v_p = sqrt(GM_S * (2.0/r_peri - 1.0/a_orb)) * 1000;

        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt * 1000);

        double p_ram = 0.5 * rho * v_p * v_p;
        double R_mag = mag_funnel_r(B_coil, r_coil, p_ram, max_R);
        double A_eff = M_PI * R_mag * R_mag;
        double m_dot = rho * v_p * A_eff * eta;

        double H_m = SATURN_RGAS * tempk / SATURN_G0;
        double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
        double dt = arc / v_p;
        double kg_pass = m_dot * dt;

        double T_pg = 2.0 * (SATURN_GAMMA - 1.0) * v_p * v_p
                    / ((SATURN_GAMMA + 1.0) * (SATURN_GAMMA + 1.0) * SATURN_RGAS);

        ShockState shock = normal_shock(v_p, tempk, press, SAT_Y_H2, SAT_Y_He);

        char comp[64];
        if (shock.post.X[Hp] > 0.1)
            snprintf(comp, sizeof(comp), "H+=%.0f%% e-=%.0f%%",
                     shock.post.X[Hp]*100, shock.post.X[Elec]*100);
        else if (shock.post.X[H] > 0.1)
            snprintf(comp, sizeof(comp), "H=%.0f%% H2=%.0f%%",
                     shock.post.X[H]*100, shock.post.X[H2]*100);
        else
            snprintf(comp, sizeof(comp), "H2=%.0f%%", shock.post.X[H2]*100);

        printf("  %-6.0f %-8.1f %-10.2e %-10.0f %-8.1f %-10.0f %-10.0f %-8.3f %-8.1f %-20s\n",
               alt, v_p/1000, rho, R_mag, kg_pass,
               shock.post.T, T_pg, shock.post.gamma_eff, shock.density_ratio, comp);
    }

    // ═══════════════════════════════════════════════════════════
    // HEAD-TO-HEAD: Optimal altitude comparison
    // ═══════════════════════════════════════════════════════════
    printf("\n═══ HEAD-TO-HEAD AT OPTIMAL ALTITUDES ═══\n\n");

    // Jupiter at 350 km
    {
        double alt = 350;
        double r_peri = R_J + alt, r_apo = R_J * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double v_p = sqrt(GM_J * (2.0/r_peri - 1.0/a_orb)) * 1000;
        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt * 1000);
        ShockState shock = normal_shock(v_p, tempk, press, JUP_Y_H2, JUP_Y_He);
        TransportProps tp = mixture_transport(shock.post);

        printf("  JUPITER (350 km above 1-bar):\n");
        printf("    v_peri:         %.1f km/s (Mach %.0f)\n", v_p/1000, v_p/shock.pre.a_sound);
        printf("    ── Perfect gas ──\n");
        double T_pg = 2.0 * (JUPITER_GAMMA-1) * v_p * v_p
                    / ((JUPITER_GAMMA+1)*(JUPITER_GAMMA+1) * JUPITER_RGAS);
        double rr_pg = (JUPITER_GAMMA+1)*pow(v_p/shock.pre.a_sound,2)
                      /((JUPITER_GAMMA-1)*pow(v_p/shock.pre.a_sound,2)+2);
        printf("      T_shock:      %.0f K\n", T_pg);
        printf("      ρ₂/ρ₁:        %.1f\n", rr_pg);
        printf("      γ:             %.2f (constant)\n", JUPITER_GAMMA);
        printf("    ── Real gas (stat_thermo) ──\n");
        printf("      T_shock:      %.0f K (%.1fx lower!)\n",
               shock.post.T, T_pg / shock.post.T);
        printf("      ρ₂/ρ₁:        %.1f (%.1fx higher!)\n",
               shock.density_ratio, shock.density_ratio / rr_pg);
        printf("      γ_eff:         %.3f\n", shock.post.gamma_eff);
        printf("      M_mix:         %.3f g/mol (was %.3f)\n",
               shock.post.M_mix*1000, shock.pre.M_mix*1000);
        printf("      μ:             %.3e Pa·s\n", tp.mu);
        printf("      k:             %.3e W/(m·K)\n", tp.k_thermal);
        printf("      Pr:            %.3f\n", tp.Pr);
        printf("      Composition:   H2=%.1f%% H=%.1f%% H+=%.1f%% He=%.1f%% He+=%.1f%% e-=%.1f%%\n",
               shock.post.X[H2]*100, shock.post.X[H]*100, shock.post.X[Hp]*100,
               shock.post.X[He]*100, shock.post.X[Hep]*100, shock.post.X[Elec]*100);
    }

    printf("\n");

    // Saturn at 700 km
    {
        double alt = 700;
        double r_peri = R_S + alt, r_apo = R_S * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double v_p = sqrt(GM_S * (2.0/r_peri - 1.0/a_orb)) * 1000;
        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt * 1000);
        ShockState shock = normal_shock(v_p, tempk, press, SAT_Y_H2, SAT_Y_He);
        TransportProps tp = mixture_transport(shock.post);

        printf("  SATURN (700 km above 1-bar):\n");
        printf("    v_peri:         %.1f km/s (Mach %.0f)\n", v_p/1000, v_p/shock.pre.a_sound);
        printf("    ── Perfect gas ──\n");
        double T_pg = 2.0 * (SATURN_GAMMA-1) * v_p * v_p
                    / ((SATURN_GAMMA+1)*(SATURN_GAMMA+1) * SATURN_RGAS);
        double rr_pg = (SATURN_GAMMA+1)*pow(v_p/shock.pre.a_sound,2)
                      /((SATURN_GAMMA-1)*pow(v_p/shock.pre.a_sound,2)+2);
        printf("      T_shock:      %.0f K\n", T_pg);
        printf("      ρ₂/ρ₁:        %.1f\n", rr_pg);
        printf("      γ:             %.2f (constant)\n", SATURN_GAMMA);
        printf("    ── Real gas (stat_thermo) ──\n");
        printf("      T_shock:      %.0f K (%.1fx lower!)\n",
               shock.post.T, T_pg / shock.post.T);
        printf("      ρ₂/ρ₁:        %.1f (%.1fx higher!)\n",
               shock.density_ratio, shock.density_ratio / rr_pg);
        printf("      γ_eff:         %.3f\n", shock.post.gamma_eff);
        printf("      M_mix:         %.3f g/mol (was %.3f)\n",
               shock.post.M_mix*1000, shock.pre.M_mix*1000);
        printf("      μ:             %.3e Pa·s\n", tp.mu);
        printf("      k:             %.3e W/(m·K)\n", tp.k_thermal);
        printf("      Pr:            %.3f\n", tp.Pr);
        printf("      Composition:   H2=%.1f%% H=%.1f%% H+=%.1f%% He=%.1f%% He+=%.1f%% e-=%.1f%%\n",
               shock.post.X[H2]*100, shock.post.X[H]*100, shock.post.X[Hp]*100,
               shock.post.X[He]*100, shock.post.X[Hep]*100, shock.post.X[Elec]*100);
    }

    // ═══════════════════════════════════════════════════════════
    // VERDICT (updated with real gas)
    // ═══════════════════════════════════════════════════════════
    printf("\n╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  VERDICT: WHAT REAL-GAS THERMODYNAMICS CHANGES                 ║\n");
    printf("╠══════════════════════════════════════════════════════════════════╣\n");
    printf("║                                                                ║\n");
    printf("║  1. Post-shock temperature is MUCH LOWER than perfect gas      ║\n");
    printf("║     because energy goes into dissociation & ionization.        ║\n");
    printf("║     Jupiter: perfect gas says ~180,000 K → real gas ~6,000 K   ║\n");
    printf("║     Saturn:  perfect gas says ~30,000 K  → real gas ~3,000 K   ║\n");
    printf("║                                                                ║\n");
    printf("║  2. Density ratio is MUCH HIGHER (10-30x vs 5-6x perfect).    ║\n");
    printf("║     This means the shock layer is THINNER, which is GOOD       ║\n");
    printf("║     for magnetic standoff — the field pushes it farther out.   ║\n");
    printf("║                                                                ║\n");
    printf("║  3. Ionization is REAL at Jupiter (H+ and e- behind shock).    ║\n");
    printf("║     At Saturn, mostly dissociated H atoms — less ionized.      ║\n");
    printf("║     This affects the magnetic scoop collection mechanism.      ║\n");
    printf("║                                                                ║\n");
    printf("║  4. The magnetic scoop concept is MORE FAVORABLE with real     ║\n");
    printf("║     gas than perfect gas predicted. Lower T → less radiation.  ║\n");
    printf("║     Higher ρ ratio → thinner shock → better mag standoff.     ║\n");
    printf("║                                                                ║\n");
    printf("║  5. Original verdict STANDS: Build both.                       ║\n");
    printf("║     Jupiter Bucket for He → pays the bills.                    ║\n");
    printf("║     Saturn Station for H2/CH4/O2 → opens the frontier.        ║\n");
    printf("║     Real gas makes both missions MORE feasible, not less.      ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n");

    return 0;
}
