// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// SATURN vs JUPITER — Head-to-Head Harvester Analysis
//
// Runs the numbers on both planets for the magnetic scoop concept:
//   1. Atmosphere at scoop altitudes (density, pressure, scale height)
//   2. Scoop physics (velocity, ram pressure, funnel, mass per pass)
//   3. Self-ionization check (bow shock temperature)
//   4. 2-year campaign simulation at each planet
//   5. Delta-V budget (Earth departure, capture, escape, return)
//   6. Resource mix (H2, He, CH4)
//   7. Moon resources (Titan, Enceladus vs Europa)
//   8. Radiation environment
//   9. Mass budget and payload delivered
//  10. Verdict
//
// References:
//   Jupiter: Galileo Probe ASI (Seiff et al. 1998), Juno (Iess et al. 2018)
//   Saturn: Voyager RSS (Lindal 1985), Cassini CIRS/INMS (Flasar 2005)
//   Titan: Huygens HASI (Fulchignoni et al. 2005)
//   Enceladus: Cassini INMS (Waite et al. 2009)
//
// Compile: g++ -std=c++11 -O2 -I../../PlanetaryModels -o saturn_vs_jupiter saturn_vs_jupiter.cpp
// Run:     ./saturn_vs_jupiter
///////////////////////////////////////////////////////////////////////////////

#include "planetary_environment.hpp"
#include <cstdio>
#include <cmath>

static const double MU0 = 4.0 * M_PI * 1e-7;

// Magnetic funnel radius: dipole B = B0*(r0/r)^3, set B = B_need
double mag_funnel_radius(double B_coil, double r_coil, double p_ram, double max_r)
{
    double B_need = sqrt(2.0 * MU0 * p_ram);
    if (B_need <= 0) return max_r;
    double R = r_coil * pow(B_coil / B_need, 1.0/3.0);
    if (R > max_r) R = max_r;
    return R;
}

// Bow shock post-shock temperature (Rankine-Hugoniot, strong shock limit M>>1)
// T_shock = 2*(gamma-1)*v^2 / ((gamma+1)^2 * R_specific)
double shock_temperature(double v_ms, double R_specific, double gamma)
{
    return 2.0 * (gamma - 1.0) * v_ms * v_ms
           / ((gamma + 1.0) * (gamma + 1.0) * R_specific);
}

int main()
{
    printf("===================================================================\n");
    printf("  SATURN vs JUPITER — Harvester Head-to-Head Analysis\n");
    printf("  Kepler Workstation\n");
    printf("===================================================================\n");

    // =====================================================================
    // Vehicle Spec (identical Bucket for fair comparison)
    // =====================================================================
    double dry_mass   = 8000;     // kg
    double tank_cap   = 33000;    // kg
    double isp        = 3000;     // s (Hall thruster)
    double ve         = isp * 9.80665e-3; // km/s
    double B_coil     = 20.0;     // Tesla
    double r_coil     = 3.0;      // meters
    double eta_collect = 0.5;     // 50% collection efficiency
    double max_R_mag   = 200.0;   // m practical funnel limit

    // Planet constants (km for orbital mechanics)
    double GM_J = 1.26713e8, R_J = 71492.0;   // km
    double GM_S = 3.7931e7,  R_S = 60268.0;   // km

    // Atmosphere composition (mass fractions)
    // Jupiter: ~86% H2 / 14% He by volume → 75% H2 / 25% He by mass
    // Saturn:  96.3% H2 / 3.25% He / 0.45% CH4 by vol → 90.6/6.1/3.4 by mass
    double JUP_H2_FRAC = 0.75, JUP_HE_FRAC = 0.25, JUP_CH4_FRAC = 0.0;
    double SAT_H2_FRAC = 0.906, SAT_HE_FRAC = 0.061, SAT_CH4_FRAC = 0.034;

    printf("\n--- Vehicle Spec (same Bucket for both planets) ---\n");
    printf("  Dry: %g kg | Tank: %g kg | Isp: %d s | Magnet: %gT, %gm\n\n",
           dry_mass, tank_cap, (int)isp, B_coil, r_coil);

    // =====================================================================
    // Section 1: Physical Constants Comparison
    // =====================================================================
    printf("--- Section 1: Physical Constants ---\n\n");
    printf("  %-20s %14s %14s %10s\n", "", "Jupiter", "Saturn", "Ratio S/J");
    printf("  %-20s %14s %14s %10s\n", "", "-------", "------", "--------");
    printf("  %-20s %14.0f %14.0f %10.2f\n", "Radius (km)", R_J, R_S, R_S/R_J);
    printf("  %-20s %14.3e %14.3e %10.3f\n", "GM (km3/s2)", GM_J, GM_S, GM_S/GM_J);
    printf("  %-20s %14.2f %14.2f %10.2f\n", "Surface g (m/s2)", JUPITER_G0, SATURN_G0, SATURN_G0/JUPITER_G0);
    printf("  %-20s %14.0f %14.0f %10.2f\n", "R_gas (J/kg/K)", JUPITER_RGAS, SATURN_RGAS, SATURN_RGAS/JUPITER_RGAS);
    printf("  %-20s %14.2f %14.2f %10.2f\n", "v_esc surface (km/s)",
           sqrt(2*GM_J/R_J), sqrt(2*GM_S/R_S), sqrt(2*GM_S/R_S)/sqrt(2*GM_J/R_J));

    // Scale heights at 1-bar
    double H_J = JUPITER_RGAS * 166.0 / JUPITER_G0; // m
    double H_S = SATURN_RGAS * 134.0 / SATURN_G0;
    printf("  %-20s %14.0f %14.0f %10.2f\n", "Scale height (m)", H_J, H_S, H_S/H_J);
    printf("  %-20s %14s %14s\n", "He mass frac", "25%", "6.1%");

    // Radiation
    printf("\n  Radiation environment:\n");
    printf("    Jupiter: ~500 rem/day at Io orbit, ~36 rem/day at Europa\n");
    printf("             Lethal to unshielded electronics in weeks\n");
    printf("    Saturn:  ~0.05 rem/day (1000x less than Jupiter)\n");
    printf("             Benign for electronics and long-duration ops\n");

    // =====================================================================
    // Section 2: Atmosphere at Scoop Altitudes
    // =====================================================================
    printf("\n--- Section 2: Atmosphere Comparison ---\n\n");

    printf("  JUPITER (above 1-bar):\n");
    printf("  %8s %10s %12s %10s %10s\n", "Alt(km)", "T(K)", "P(Pa)", "rho(kg/m3)", "H(km)");
    double j_alts[] = {200, 300, 350, 400, 500, 700, 1000};
    for (double alt : j_alts) {
        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt*1000);
        double H = JUPITER_RGAS * tempk / JUPITER_G0 / 1000;
        printf("  %8.0f %10.1f %12.3e %10.3e %10.1f\n", alt, tempk, press, rho, H);
    }

    printf("\n  SATURN (above 1-bar):\n");
    printf("  %8s %10s %12s %10s %10s\n", "Alt(km)", "T(K)", "P(Pa)", "rho(kg/m3)", "H(km)");
    double s_alts[] = {200, 300, 350, 400, 500, 700, 1000, 1500, 2000};
    for (double alt : s_alts) {
        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt*1000);
        double H = SATURN_RGAS * tempk / SATURN_G0 / 1000;
        printf("  %8.0f %10.1f %12.3e %10.3e %10.1f\n", alt, tempk, press, rho, H);
    }

    printf("\n  KEY INSIGHT: Saturn's scale height is ~2.5x Jupiter's (lower gravity,\n");
    printf("  similar R_gas). The atmosphere extends MUCH higher. At 350 km altitude,\n");
    printf("  Saturn is orders of magnitude denser than Jupiter.\n");

    // =====================================================================
    // Section 3: Scoop Physics Sweep
    // =====================================================================
    printf("\n--- Section 3: Scoop Physics at Optimal Altitude ---\n\n");

    // Sweep altitudes to find optimal scoop altitude for each planet
    printf("  JUPITER scoop sweep (orbit: alt_peri x 5Rj):\n");
    printf("  %6s %10s %10s %10s %10s %10s %10s\n",
           "Alt", "v_peri", "rho", "p_ram", "R_mag", "kg/pass", "T_shock");
    printf("  %6s %10s %10s %10s %10s %10s %10s\n",
           "(km)", "(km/s)", "(kg/m3)", "(Pa)", "(m)", "", "(K)");

    double best_j_kg = 0, best_j_alt = 0;
    double j_sweep[] = {250, 300, 350, 400, 500, 600};
    for (double alt : j_sweep) {
        double r_peri = R_J + alt;
        double r_apo = R_J * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double T_orb = 2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_J);
        double v_p = sqrt(GM_J*(2.0/r_peri - 1.0/a_orb)) * 1000; // m/s

        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt*1000);

        double p_ram = 0.5 * rho * v_p * v_p;
        double R_mag = mag_funnel_radius(B_coil, r_coil, p_ram, max_R_mag);
        double A_eff = M_PI * R_mag * R_mag;
        double m_dot = rho * v_p * A_eff * eta_collect;

        // Pass duration
        double H_m = JUPITER_RGAS * tempk / JUPITER_G0;
        double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
        double dt = arc / v_p;
        double kg_pass = m_dot * dt;

        double T_shock = shock_temperature(v_p, JUPITER_RGAS, JUPITER_GAMMA);

        printf("  %6.0f %10.1f %10.2e %10.3f %10.0f %10.2f %10.0f\n",
               alt, v_p/1000, rho, p_ram, R_mag, kg_pass, T_shock);

        if (kg_pass > best_j_kg) { best_j_kg = kg_pass; best_j_alt = alt; }
    }
    printf("  → Best: %.0f km alt, %.2f kg/pass\n", best_j_alt, best_j_kg);

    printf("\n  SATURN scoop sweep (orbit: alt_peri x 5Rs):\n");
    printf("  %6s %10s %10s %10s %10s %10s %10s\n",
           "Alt", "v_peri", "rho", "p_ram", "R_mag", "kg/pass", "T_shock");
    printf("  %6s %10s %10s %10s %10s %10s %10s\n",
           "(km)", "(km/s)", "(kg/m3)", "(Pa)", "(m)", "", "(K)");

    double best_s_kg = 0, best_s_alt = 0;
    double s_sweep[] = {350, 500, 700, 1000, 1200, 1500, 1800, 2000};
    for (double alt : s_sweep) {
        double r_peri = R_S + alt;
        double r_apo = R_S * 5;
        double a_orb = (r_peri + r_apo) / 2;
        double T_orb = 2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_S);
        double v_p = sqrt(GM_S*(2.0/r_peri - 1.0/a_orb)) * 1000; // m/s

        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt*1000);

        double p_ram = 0.5 * rho * v_p * v_p;
        double R_mag = mag_funnel_radius(B_coil, r_coil, p_ram, max_R_mag);
        double A_eff = M_PI * R_mag * R_mag;
        double m_dot = rho * v_p * A_eff * eta_collect;

        double H_m = SATURN_RGAS * tempk / SATURN_G0;
        double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
        double dt = arc / v_p;
        double kg_pass = m_dot * dt;

        double T_shock = shock_temperature(v_p, SATURN_RGAS, SATURN_GAMMA);

        printf("  %6.0f %10.1f %10.2e %10.3f %10.0f %10.2f %10.0f\n",
               alt, v_p/1000, rho, p_ram, R_mag, kg_pass, T_shock);

        if (kg_pass > best_s_kg) { best_s_kg = kg_pass; best_s_alt = alt; }
    }
    printf("  → Best: %.0f km alt, %.2f kg/pass\n", best_s_alt, best_s_kg);

    printf("\n  IONIZATION CHECK:\n");
    // At Jupiter best alt
    {
        double r_peri = R_J + best_j_alt;
        double v_p = sqrt(GM_J*(2.0/r_peri - 1.0/(r_peri + R_J*5)*2)) * 1000;
        // Simpler: recompute
        double a_orb = (r_peri + R_J*5) / 2;
        v_p = sqrt(GM_J*(2.0/r_peri - 1.0/a_orb)) * 1000;
        double m_h2 = 2 * 1.67e-27;
        double KE = 0.5 * m_h2 * v_p * v_p;
        double eV = KE / 1.6e-19;
        double T_shock = shock_temperature(v_p, JUPITER_RGAS, JUPITER_GAMMA);
        printf("    Jupiter (%.0f km): v=%.1f km/s, H2 KE=%.1f eV (need 15.4), "
               "T_shock=%.0f K → %s\n",
               best_j_alt, v_p/1000, eV, T_shock,
               T_shock > 10000 ? "IONIZED" : "MARGINAL");
    }
    {
        double r_peri = R_S + best_s_alt;
        double a_orb = (r_peri + R_S*5) / 2;
        double v_p = sqrt(GM_S*(2.0/r_peri - 1.0/a_orb)) * 1000;
        double m_h2 = 2 * 1.67e-27;
        double KE = 0.5 * m_h2 * v_p * v_p;
        double eV = KE / 1.6e-19;
        double T_shock = shock_temperature(v_p, SATURN_RGAS, SATURN_GAMMA);
        printf("    Saturn  (%.0f km): v=%.1f km/s, H2 KE=%.1f eV (need 15.4), "
               "T_shock=%.0f K → %s\n",
               best_s_alt, v_p/1000, eV, T_shock,
               T_shock > 10000 ? "IONIZED" : "MARGINAL");
    }
    printf("\n    Per-molecule KE is below ionization threshold at Saturn,\n");
    printf("    BUT the bow shock thermalizes kinetic energy collectively.\n");
    printf("    Post-shock temperature >> 10,000 K at both planets → full ionization.\n");
    printf("    The magnetic scoop works at Saturn.\n");

    // =====================================================================
    // Section 4: Campaign Simulation (both planets, best altitude)
    // =====================================================================
    printf("\n--- Section 4: 2-Year Scoop Campaign ---\n");

    // Jupiter campaign
    double j_peri_alt = best_j_alt;
    double j_r_peri = R_J + j_peri_alt;
    double j_r_apo = R_J * 5;
    double j_a = (j_r_peri + j_r_apo) / 2;
    double j_T_orb = 2*M_PI*sqrt(j_a*j_a*j_a/GM_J);
    double campaign_sec = 2 * 365.25 * 86400;
    int j_max_passes = (int)(campaign_sec / j_T_orb);

    double j_total = 0;
    int j_passes_used = 0;
    for (int p = 1; p <= j_max_passes; p++) {
        double progress = (double)p / j_max_passes;
        double alt = j_peri_alt + 50 - 100*progress; // vary altitude
        if (alt < 250) alt = 250;
        double rp = R_J + alt;
        double vp = sqrt(GM_J*(2.0/rp - 1.0/j_a)) * 1000;
        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt*1000);
        double p_ram = 0.5*rho*vp*vp;
        double Rm = mag_funnel_radius(B_coil, r_coil, p_ram, max_R_mag);
        double Ae = M_PI*Rm*Rm;
        double H_m = JUPITER_RGAS*tempk/JUPITER_G0;
        double arc = 2.0*sqrt(2.0*rp*1000*H_m);
        double dt = arc/vp;
        double harvested = rho*vp*Ae*eta_collect*dt;
        j_total += harvested;
        j_passes_used = p;
        if (j_total >= tank_cap) break;
    }

    printf("\n  JUPITER (%.0f km periapsis, %g km x %.0f Rj orbit):\n",
           j_peri_alt, j_peri_alt, j_r_apo/R_J);
    printf("    Orbit period:   %.1f hours\n", j_T_orb/3600);
    printf("    Max passes/2yr: %d\n", j_max_passes);
    printf("    Passes used:    %d\n", j_passes_used);
    printf("    Total scooped:  %.1f tonnes (%.0f%% of tank)\n",
           j_total/1000, j_total/tank_cap*100);
    printf("    H2: %.1ft, He: %.1ft\n",
           j_total*JUP_H2_FRAC/1000, j_total*JUP_HE_FRAC/1000);

    // Saturn campaign
    double s_peri_alt = best_s_alt;
    double s_r_peri = R_S + s_peri_alt;
    double s_r_apo = R_S * 5;
    double s_a = (s_r_peri + s_r_apo) / 2;
    double s_T_orb = 2*M_PI*sqrt(s_a*s_a*s_a/GM_S);
    int s_max_passes = (int)(campaign_sec / s_T_orb);

    double s_total = 0;
    int s_passes_used = 0;
    for (int p = 1; p <= s_max_passes; p++) {
        double progress = (double)p / s_max_passes;
        double alt = s_peri_alt + 100 - 200*progress;
        if (alt < s_peri_alt - 200) alt = s_peri_alt - 200;
        double rp = R_S + alt;
        double vp = sqrt(GM_S*(2.0/rp - 1.0/s_a)) * 1000;
        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt*1000);
        double p_ram = 0.5*rho*vp*vp;
        double Rm = mag_funnel_radius(B_coil, r_coil, p_ram, max_R_mag);
        double Ae = M_PI*Rm*Rm;
        double H_m = SATURN_RGAS*tempk/SATURN_G0;
        double arc = 2.0*sqrt(2.0*rp*1000*H_m);
        double dt = arc/vp;
        double harvested = rho*vp*Ae*eta_collect*dt;
        s_total += harvested;
        s_passes_used = p;
        if (s_total >= tank_cap) break;
    }

    printf("\n  SATURN (%.0f km periapsis, %g km x %.0f Rs orbit):\n",
           s_peri_alt, s_peri_alt, s_r_apo/R_S);
    printf("    Orbit period:   %.1f hours\n", s_T_orb/3600);
    printf("    Max passes/2yr: %d\n", s_max_passes);
    printf("    Passes used:    %d\n", s_passes_used);
    printf("    Total scooped:  %.1f tonnes (%.0f%% of tank)\n",
           s_total/1000, s_total/tank_cap*100);
    printf("    H2: %.1ft, He: %.1ft, CH4: %.1ft\n",
           s_total*SAT_H2_FRAC/1000, s_total*SAT_HE_FRAC/1000,
           s_total*SAT_CH4_FRAC/1000);

    // Clamp to tank capacity
    double j_scooped = (j_total > tank_cap) ? tank_cap : j_total;
    double s_scooped = (s_total > tank_cap) ? tank_cap : s_total;

    // =====================================================================
    // Section 5: Delta-V Budget
    // =====================================================================
    printf("\n--- Section 5: Delta-V Budget ---\n\n");

    // Jupiter Oberth escape from scoop orbit
    double j_v_peri = sqrt(GM_J * (2.0/j_r_peri - 1.0/j_a));
    double j_v_esc = sqrt(2.0 * GM_J / j_r_peri);
    double j_dv_oberth = j_v_esc - j_v_peri;

    // Saturn Oberth escape from scoop orbit
    double s_v_peri = sqrt(GM_S * (2.0/s_r_peri - 1.0/s_a));
    double s_v_esc = sqrt(2.0 * GM_S / s_r_peri);
    double s_dv_oberth = s_v_esc - s_v_peri;

    // Full budget
    double j_dv_depart = 3.5;   // from LEO via VEEJ
    double j_dv_capture = 0.5;  // magnetic drag + small burn
    double j_dv_escape = 5.0;   // Oberth + v_inf for Earth return
    double j_dv_margin = 0.5;
    double j_dv_total = j_dv_depart + j_dv_capture + j_dv_escape + j_dv_margin;

    // Saturn: longer transfer, use VVEJGA (Cassini-like)
    // v_inf at Saturn ~5.5 km/s via gravity assists
    double s_dv_depart = 3.8;   // from LEO (slightly higher C3 for Saturn)
    double s_dv_capture = 0.5;  // similar magnetic drag concept
    // Oberth escape is cheaper at Saturn (shallower well)
    double s_dv_escape = j_dv_escape * (s_dv_oberth / j_dv_oberth); // scale by Oberth ratio
    if (s_dv_escape < s_dv_oberth + 1.0) s_dv_escape = s_dv_oberth + 1.0; // need v_inf too
    double s_dv_margin = 0.5;
    double s_dv_total = s_dv_depart + s_dv_capture + s_dv_escape + s_dv_margin;

    printf("  %-24s %10s %10s %10s\n", "", "Jupiter", "Saturn", "Savings");
    printf("  %-24s %10s %10s %10s\n", "", "-------", "------", "-------");
    printf("  %-24s %9.1f  %9.1f  %+8.1f km/s\n", "Earth departure",
           j_dv_depart, s_dv_depart, j_dv_depart - s_dv_depart);
    printf("  %-24s %9.1f  %9.1f  %+8.1f km/s\n", "Planet capture",
           j_dv_capture, s_dv_capture, j_dv_capture - s_dv_capture);
    printf("  %-24s %9.1f  %9.1f  %+8.1f km/s\n", "Planet escape",
           j_dv_escape, s_dv_escape, j_dv_escape - s_dv_escape);
    printf("  %-24s %9.1f  %9.1f  %+8.1f km/s\n", "Margin",
           j_dv_margin, s_dv_margin, 0.0);
    printf("  %-24s %9.1f  %9.1f  %+8.1f km/s\n", "TOTAL",
           j_dv_total, s_dv_total, j_dv_total - s_dv_total);

    printf("\n  Oberth effect at periapsis:\n");
    printf("    Jupiter: v_orbit=%.1f, v_escape=%.1f, dv_Oberth=%.1f km/s\n",
           j_v_peri, j_v_esc, j_dv_oberth);
    printf("    Saturn:  v_orbit=%.1f, v_escape=%.1f, dv_Oberth=%.1f km/s\n",
           s_v_peri, s_v_esc, s_dv_oberth);
    printf("    Saturn Oberth is %.0f%% of Jupiter's → less fuel to escape\n",
           s_dv_oberth/j_dv_oberth*100);

    // =====================================================================
    // Section 6: Mass Budget & Payload
    // =====================================================================
    printf("\n--- Section 6: Mass Budget ---\n\n");

    // Tsiolkovsky: m0/mf = exp(dv/ve)
    double j_mr = exp(j_dv_total / ve);
    double s_mr = exp(s_dv_total / ve);

    double j_h2_tank = j_scooped * JUP_H2_FRAC;
    double j_he_tank = j_scooped * JUP_HE_FRAC;
    double j_ch4_tank = 0;
    double j_h2_payload = (j_h2_tank - (dry_mass + j_he_tank) * (j_mr - 1.0)) / j_mr;
    if (j_h2_payload < 0) j_h2_payload = 0;
    double j_delivered = j_h2_payload + j_he_tank;

    double s_h2_tank = s_scooped * SAT_H2_FRAC;
    double s_he_tank = s_scooped * SAT_HE_FRAC;
    double s_ch4_tank = s_scooped * SAT_CH4_FRAC;
    double s_h2_payload = (s_h2_tank - (dry_mass + s_he_tank + s_ch4_tank) * (s_mr - 1.0)) / s_mr;
    if (s_h2_payload < 0) s_h2_payload = 0;
    double s_delivered = s_h2_payload + s_he_tank + s_ch4_tank;

    printf("  %-28s %10s %10s\n", "", "Jupiter", "Saturn");
    printf("  %-28s %10s %10s\n", "", "-------", "------");
    printf("  %-28s %9.1f  %9.1f km/s\n", "Total dV", j_dv_total, s_dv_total);
    printf("  %-28s %9.3f  %9.3f\n", "Mass ratio", j_mr, s_mr);
    printf("  %-28s %9.1f  %9.1f t\n", "Scooped total", j_scooped/1000, s_scooped/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "H2 scooped", j_h2_tank/1000, s_h2_tank/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "He scooped", j_he_tank/1000, s_he_tank/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "CH4 scooped", j_ch4_tank/1000, s_ch4_tank/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "H2 as fuel", (j_h2_tank-j_h2_payload)/1000,
           (s_h2_tank-s_h2_payload)/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "H2 delivered", j_h2_payload/1000, s_h2_payload/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "He delivered", j_he_tank/1000, s_he_tank/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "CH4 delivered", j_ch4_tank/1000, s_ch4_tank/1000);
    printf("  %-28s %9.1f  %9.1f t\n", "TOTAL DELIVERED", j_delivered/1000, s_delivered/1000);
    printf("  %-28s %9.1fx %9.1fx\n", "Delivery ratio (x dry)",
           j_delivered/dry_mass, s_delivered/dry_mass);

    // =====================================================================
    // Section 7: Timeline & Annual Delivery
    // =====================================================================
    printf("\n--- Section 7: Timeline ---\n\n");

    // Jupiter loop: ~5yr transit out, 2yr scoop, ~3yr return = ~10yr
    double j_transit_out = 5.5;  // years (VEEJ)
    double j_scoop_time = 2.0;
    double j_transit_ret = 3.0;
    double j_loop = j_transit_out + j_scoop_time + j_transit_ret;

    // Saturn loop: ~7yr transit out, 2yr scoop, ~5yr return = ~14yr
    double s_transit_out = 7.0;  // years (VVEJGA, like Cassini)
    double s_scoop_time = 2.0;
    double s_transit_ret = 5.0;
    double s_loop = s_transit_out + s_scoop_time + s_transit_ret;

    printf("  %-28s %10s %10s\n", "", "Jupiter", "Saturn");
    printf("  %-28s %10s %10s\n", "", "-------", "------");
    printf("  %-28s %9.1f  %9.1f yr\n", "Outbound transit", j_transit_out, s_transit_out);
    printf("  %-28s %9.1f  %9.1f yr\n", "Scoop campaign", j_scoop_time, s_scoop_time);
    printf("  %-28s %9.1f  %9.1f yr\n", "Return transit", j_transit_ret, s_transit_ret);
    printf("  %-28s %9.1f  %9.1f yr\n", "LOOP TIME", j_loop, s_loop);

    double j_synodic = 398.88 / 365.25;  // years
    double s_synodic = 378.09 / 365.25;

    printf("  %-28s %9.2f  %9.2f yr\n", "Synodic period", j_synodic, s_synodic);

    // Fleet: 10 buckets, annual launches
    int fleet = 10;
    double j_annual_t = fleet * (1.0/j_loop) * j_delivered / 1000;
    double s_annual_t = fleet * (1.0/s_loop) * s_delivered / 1000;

    printf("\n  Fleet (10 buckets, continuous ops):\n");
    printf("    Jupiter: %.1f tonnes/year delivered\n", j_annual_t);
    printf("    Saturn:  %.1f tonnes/year delivered\n", s_annual_t);

    // =====================================================================
    // Section 8: Moon Resources
    // =====================================================================
    printf("\n--- Section 8: Moon Resources ---\n\n");

    printf("  JUPITER MOONS:\n");
    printf("    Europa:   H2O ice crust, subsurface ocean\n");
    printf("              Mine ice → electrolyze → O2 + H2\n");
    printf("              Surface g=1.31 m/s2, escape=2.03 km/s\n");
    printf("              Radiation: ~5.4 Sv/day (lethal, robotic only)\n\n");
    printf("    Ganymede: H2O ice + rock, weak magnetosphere\n");
    printf("              Surface g=1.43 m/s2, escape=2.74 km/s\n");
    printf("              Some radiation shielding from intrinsic field\n\n");

    printf("  SATURN MOONS:\n");
    printf("    Titan:    1.47 bar N2/CH4 atmosphere, 94 K surface\n");
    printf("              METHANE LAKES (Ligeia Mare: ~9,000 km2)\n");
    printf("              Pump liquid CH4 directly — no mining needed\n");
    printf("              Surface g=1.35 m/s2, escape=2.64 km/s\n");
    printf("              Radiation: negligible (Saturn's field is weak)\n");

    // Titan atmosphere data
    double t_rho, t_press, t_tempk;
    atmosphere_titan(t_rho, t_press, t_tempk, 0);
    printf("              Atmosphere: %.0f K, %.0f Pa (%.2f bar), %.2f kg/m3\n",
           t_tempk, t_press, t_press/1e5, t_rho);
    printf("              → Aerobrake for FREE at Titan (1.47 bar!)\n");
    printf("              → Parachute landing, surface ops, balloon science\n\n");

    printf("    Enceladus: WATER GEYSERS from south pole (Cassini INMS)\n");
    printf("              Plume composition: 98%% H2O, 1%% CO2, 0.8%% H2\n");
    printf("              Fly through plume → collect H2O directly\n");
    printf("              Surface g=0.11 m/s2, escape=0.24 km/s\n");
    printf("              Radiation: negligible\n");
    printf("              No mining equipment needed — just fly through!\n\n");

    // Compare O2 sourcing
    double o2_needed = 20000; // kg per mission
    printf("  O2 sourcing comparison (%.0f kg needed per mission):\n", o2_needed);
    printf("    Europa: Mine ice, haul to surface, electrolyze\n");
    printf("            Land/launch dV: ~4 km/s round trip\n");
    printf("            Radiation: 5.4 Sv/day (robotic, limited ops window)\n");
    printf("    Enceladus: Fly through geysers at 0.2 km/s\n");
    printf("            Plume intercept dV: ~0.1 km/s\n");
    printf("            Radiation: safe for crewed ops\n");
    printf("    Winner: Enceladus (not even close)\n");

    // =====================================================================
    // Section 9: Economics
    // =====================================================================
    printf("\n--- Section 9: Value Per Trip ---\n\n");

    double he_price = 1500;  // $/kg (future shortage)
    double h2_price = 500;   // $/kg (propellant depot value)
    double ch4_price = 200;  // $/kg (propellant/life support)

    double j_value = j_h2_payload * h2_price / 1000
                   + j_he_tank * he_price / 1000;
    double s_value = s_h2_payload * h2_price / 1000
                   + s_he_tank * he_price / 1000
                   + s_ch4_tank * ch4_price / 1000;

    printf("  %-28s %10s %10s\n", "", "Jupiter", "Saturn");
    printf("  %-28s %10s %10s\n", "", "-------", "------");
    printf("  %-28s  $%7.0fM  $%7.0fM\n", "H2 value ($500/kg)",
           j_h2_payload*h2_price/1e6, s_h2_payload*h2_price/1e6);
    printf("  %-28s  $%7.0fM  $%7.0fM\n", "He value ($1500/kg)",
           j_he_tank*he_price/1e6, s_he_tank*he_price/1e6);
    printf("  %-28s  $%7.0fM  $%7.0fM\n", "CH4 value ($200/kg)",
           j_ch4_tank*ch4_price/1e6, s_ch4_tank*ch4_price/1e6);
    printf("  %-28s  $%7.0fM  $%7.0fM\n", "TOTAL PER TRIP",
           j_value/1000, s_value/1000);
    printf("  %-28s  $%7.0fM  $%7.0fM\n", "Per year (10 buckets)",
           j_value/1000 * fleet/j_loop, s_value/1000 * fleet/s_loop);

    // =====================================================================
    // Section 10: Verdict
    // =====================================================================
    printf("\n===================================================================\n");
    printf("  VERDICT\n");
    printf("===================================================================\n");

    printf("\n  JUPITER wins on:\n");
    printf("    - He fraction: 25%% vs 6.1%% by mass → 4x more He per trip\n");
    printf("    - Shorter loop: %.0f yr vs %.0f yr → faster turnaround\n", j_loop, s_loop);
    printf("    - Higher orbital velocity → better Oberth effect\n");
    printf("    - He economics: He is THE high-value product\n");
    printf("    - Proven trajectory design (VEEJ, well-studied)\n");

    printf("\n  SATURN wins on:\n");
    printf("    - Radiation: 1000x less → electronics last decades, crewable\n");
    printf("    - Escape dV: %.1f km/s vs %.1f km/s → %.0f%% less fuel to leave\n",
           s_dv_oberth, j_dv_oberth, (1-s_dv_oberth/j_dv_oberth)*100);
    printf("    - Payload ratio: %.1fx vs %.1fx dry mass\n",
           s_delivered/dry_mass, j_delivered/dry_mass);
    printf("    - Titan: liquid methane on the surface, 1.47 bar atmo\n");
    printf("    - Enceladus: free water from geysers → O2 with zero mining\n");
    printf("    - Atmospheric density: extended atmosphere, high scale height\n");
    printf("    - Human-rated: could support crewed operations\n");

    printf("\n  THE REAL QUESTION:\n");
    printf("    What are you harvesting FOR?\n\n");
    printf("    If HELIUM for Earth: Jupiter. No contest. 4x the He per kg.\n");
    printf("    Earth's helium shortage is a Jupiter problem.\n\n");
    printf("    If PROPELLANT DEPOT for outer solar system: Saturn.\n");
    printf("    H2 from atmosphere + CH4 from Titan + O2 from Enceladus\n");
    printf("    = self-sustaining fuel depot at Saturn.\n");
    printf("    Crewed. Reusable. Low radiation. Accessible moons.\n");
    printf("    The gas station on the highway to the outer planets.\n\n");
    printf("    If BOTH: Start at Jupiter (He for Earth economics),\n");
    printf("    then build the Saturn depot for expansion.\n");

    printf("\n===================================================================\n");
    printf("  RECOMMENDATION: Build both. Different missions.\n");
    printf("  Jupiter Bucket for He → pays the bills.\n");
    printf("  Saturn Station for H2/CH4/O2 → opens the frontier.\n");
    printf("===================================================================\n");

    return 0;
}
