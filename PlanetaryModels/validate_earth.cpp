// Validation: Earth atmosphere conditions from airline cruise to re-entry
//
// Test 1: FL390 at 420 KTAS — verify SAT and TAT match ISA + real-world
// Test 2: Earth re-entry at 7.8 km/s (LEO) — compare to published data
// Test 3: Earth re-entry at 11 km/s (lunar return) — Apollo-class heating
// Test 4: Earth re-entry at 14.7 km/s (Jupiter return) — Bucket aerocapture
//
// SAT = Static Air Temperature (freestream, = ISA at altitude)
// TAT = Total Air Temperature (stagnation, = SAT * (1 + (γ-1)/2 * M²))
// Recovery factor r ≈ 0.98 for turbulent BL: TAT_probe = SAT*(1 + r*(γ-1)/2*M²)
//
// Compile: g++ -std=c++11 -O2 -o validate_earth validate_earth.cpp
// Run:     ./validate_earth

#include "stat_thermo.hpp"
#include "transport_properties.hpp"
#include "planetary_environment.hpp"
#include <cstdio>
#include <cmath>

using namespace stat_thermo;
using namespace transport;

// ISA (International Standard Atmosphere) for Earth
// Simplified: troposphere + stratosphere
void isa_earth(double alt_ft, double &T_K, double &p_Pa, double &rho)
{
    double alt_m = alt_ft * 0.3048;

    if (alt_m < 11000) {
        // Troposphere: lapse rate -6.5 K/km
        T_K = 288.15 - 0.0065 * alt_m;
        p_Pa = 101325.0 * pow(T_K / 288.15, 5.2561);
    } else if (alt_m < 20000) {
        // Lower stratosphere: isothermal at 216.65 K
        T_K = 216.65;
        double p_11 = 101325.0 * pow(216.65 / 288.15, 5.2561);
        p_Pa = p_11 * exp(-9.80665 * (alt_m - 11000) / (287.053 * 216.65));
    } else if (alt_m < 32000) {
        // Upper stratosphere: lapse +1.0 K/km
        T_K = 216.65 + 0.001 * (alt_m - 20000);
        double p_20 = 101325.0 * pow(216.65 / 288.15, 5.2561)
                     * exp(-9.80665 * 9000 / (287.053 * 216.65));
        p_Pa = p_20 * pow(T_K / 216.65, -34.1632);
    } else {
        // Simplified upper atmosphere (50-100 km for re-entry)
        // Exponential decay
        T_K = 270.0 - 0.003 * (alt_m - 32000);
        if (T_K < 180) T_K = 180;
        double p_32 = 868.0;  // Pa at 32 km
        double H = 287.053 * T_K / 9.80665;
        p_Pa = p_32 * exp(-(alt_m - 32000) / H);
    }

    // Air properties
    double R_air = 287.053;  // J/(kg·K) for dry air
    rho = p_Pa / (R_air * T_K);
}

// Speed of sound in air
double a_air(double T_K)
{
    return sqrt(1.4 * 287.053 * T_K);
}

// Knots true airspeed to m/s
double ktas_to_ms(double ktas) { return ktas * 0.514444; }

// Sutton-Graves for Earth (N2/O2)
double sutton_graves_earth(double rho, double v, double R_nose)
{
    // k = 1.7415e-4 for air (Sutton-Graves)
    return 1.7415e-4 * sqrt(rho / R_nose) * pow(v, 3) / 1e3;  // kW/m²
}


int main()
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  VALIDATION: EARTH ATMOSPHERE — ISA + REENTRY              ║\n");
    printf("║  Compare model against known physics at every speed regime ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    // ═══════════════════════════════════════════════════════════════
    // TEST 1: FL390 at 420 KTAS — Airline Cruise
    // ═══════════════════════════════════════════════════════════════
    printf("═══ TEST 1: AIRLINE CRUISE — FL390, 420 KTAS ═══\n\n");

    double alt_ft = 39000;
    double v_ktas = 420;
    double v_ms = ktas_to_ms(v_ktas);

    double T_isa, p_isa, rho_isa;
    isa_earth(alt_ft, T_isa, p_isa, rho_isa);

    double a = a_air(T_isa);
    double Mach = v_ms / a;

    // SAT = freestream static temperature
    double SAT_K = T_isa;
    double SAT_C = SAT_K - 273.15;

    // TAT = total (stagnation) temperature
    // TAT = SAT * (1 + (γ-1)/2 * M²)
    double gamma_air = 1.4;
    double TAT_K = SAT_K * (1.0 + (gamma_air - 1.0) / 2.0 * Mach * Mach);
    double TAT_C = TAT_K - 273.15;

    // Probe TAT (with recovery factor r ≈ 0.98 for turbulent BL)
    double r_recovery = 0.98;
    double TAT_probe_K = SAT_K * (1.0 + r_recovery * (gamma_air - 1.0) / 2.0 * Mach * Mach);
    double TAT_probe_C = TAT_probe_K - 273.15;

    // Ram rise
    double ram_rise = TAT_C - SAT_C;

    printf("  Altitude:          FL%03.0f (%.0f ft / %.0f m)\n",
           alt_ft / 100, alt_ft, alt_ft * 0.3048);
    printf("  True airspeed:     %g KTAS (%.1f m/s)\n", v_ktas, v_ms);
    printf("  ISA temperature:   %.2f K (%.2f °C)\n", T_isa, T_isa - 273.15);
    printf("  ISA pressure:      %.0f Pa (%.1f hPa / %.2f inHg)\n",
           p_isa, p_isa / 100, p_isa / 3386.39);
    printf("  ISA density:       %.4f kg/m³\n", rho_isa);
    printf("  Speed of sound:    %.1f m/s (%.0f kts)\n", a, a / 0.514444);
    printf("  Mach number:       %.3f\n", Mach);
    printf("\n");
    printf("  ┌─────────────────────────────────────────────────┐\n");
    printf("  │  SAT (Static):    %+7.2f °C  (%6.2f K)         │\n", SAT_C, SAT_K);
    printf("  │  TAT (Total):     %+7.2f °C  (%6.2f K)         │\n", TAT_C, TAT_K);
    printf("  │  TAT probe (r=.98): %+6.2f °C  (%6.2f K)       │\n", TAT_probe_C, TAT_probe_K);
    printf("  │  Ram rise:        %+7.2f °C                     │\n", ram_rise);
    printf("  └─────────────────────────────────────────────────┘\n");
    printf("\n");
    printf("  Expected (ISA @ FL390):\n");
    printf("    SAT: -56.5 °C (ISA standard tropopause)\n");
    printf("    TAT: ~-24 to -28 °C (typical for M0.78-0.82)\n");
    printf("    Ram rise: ~28-30 °C\n");
    printf("\n");

    // Also compute at other common cruise conditions
    printf("  ── Cruise condition sweep ──\n");
    printf("  %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n",
           "FL", "KTAS", "Mach", "SAT(°C)", "TAT(°C)", "Ram(°C)", "q(Pa)");
    printf("  ──────────────────────────────────────────────────────────\n");

    double fls[] = {350, 370, 390, 410, 430, 450};
    double speeds[] = {380, 400, 420, 440, 460, 480};
    for (int i = 0; i < 6; i++) {
        double T, p, rho;
        isa_earth(fls[i] * 100, T, p, rho);
        double v = ktas_to_ms(speeds[i]);
        double M = v / a_air(T);
        double sat = T - 273.15;
        double tat = T * (1.0 + 0.2 * M * M) - 273.15;
        double q = 0.5 * rho * v * v;
        printf("  %-8.0f %-8.0f %-8.3f %-8.1f %-8.1f %-8.1f %-8.0f\n",
               fls[i], speeds[i], M, sat, tat, tat - sat, q);
    }

    // ═══════════════════════════════════════════════════════════════
    // TEST 2: Earth Re-entry from LEO (7.8 km/s)
    // ═══════════════════════════════════════════════════════════════
    printf("\n\n═══ TEST 2: EARTH REENTRY — LEO (7.8 km/s, 70 km) ═══\n\n");

    double v_leo = 7800;  // m/s
    double alt_entry_m = 70000;  // 70 km
    double T_70, p_70, rho_70;
    isa_earth(alt_entry_m / 0.3048, T_70, p_70, rho_70);

    double M_leo = v_leo / a_air(T_70);
    double R_nose_capsule = 1.5;  // m (Apollo-class)
    double q_leo = sutton_graves_earth(rho_70, v_leo, R_nose_capsule);

    // Perfect gas post-shock
    double T_shock_pg = T_70 * (2.0 * gamma_air * M_leo * M_leo - (gamma_air - 1.0))
                       * ((gamma_air - 1.0) * M_leo * M_leo + 2.0)
                       / ((gamma_air + 1.0) * (gamma_air + 1.0) * M_leo * M_leo);
    double rho_ratio_pg = (gamma_air + 1.0) * M_leo * M_leo
                        / ((gamma_air - 1.0) * M_leo * M_leo + 2.0);

    printf("  Altitude:          70 km\n");
    printf("  Velocity:          %.1f km/s (Mach %.0f)\n", v_leo / 1000, M_leo);
    printf("  Freestream:        T=%.1f K, p=%.2f Pa, ρ=%.2e kg/m³\n", T_70, p_70, rho_70);
    printf("  Stagnation heat:   %.0f kW/m² (Sutton-Graves, R_nose=%.1fm)\n", q_leo, R_nose_capsule);
    printf("  Perfect gas shock: T=%.0f K, ρ₂/ρ₁=%.1f\n", T_shock_pg, rho_ratio_pg);
    printf("\n");
    printf("  Expected (published, Apollo-class):\n");
    printf("    Peak heating:    ~1000-3000 kW/m² at peak dynamic pressure\n");
    printf("    Post-shock T:    ~7000-11000 K (real gas, with dissociation)\n");
    printf("    ρ₂/ρ₁:          ~12-15 (real gas) vs %.1f (perfect gas)\n", rho_ratio_pg);
    printf("    Real gas absorbs energy in N₂/O₂ dissociation → lower T, higher ρ\n");

    // ═══════════════════════════════════════════════════════════════
    // TEST 3: Lunar return (11 km/s)
    // ═══════════════════════════════════════════════════════════════
    printf("\n\n═══ TEST 3: EARTH REENTRY — LUNAR RETURN (11 km/s, 70 km) ═══\n\n");

    double v_lunar = 11000;
    double M_lunar = v_lunar / a_air(T_70);
    double q_lunar = sutton_graves_earth(rho_70, v_lunar, R_nose_capsule);
    double T_shock_lunar_pg = T_70 * (2.0 * gamma_air * M_lunar * M_lunar - (gamma_air - 1.0))
                             * ((gamma_air - 1.0) * M_lunar * M_lunar + 2.0)
                             / ((gamma_air + 1.0) * (gamma_air + 1.0) * M_lunar * M_lunar);

    printf("  Velocity:          11 km/s (Mach %.0f)\n", M_lunar);
    printf("  Stagnation heat:   %.0f kW/m² (Sutton-Graves)\n", q_lunar);
    printf("  Perfect gas shock: T=%.0f K\n", T_shock_lunar_pg);
    printf("\n");
    printf("  Expected (Apollo 4 flight data):\n");
    printf("    Peak heating:    ~5000 kW/m²\n");
    printf("    Post-shock T:    ~11000-15000 K (heavily dissociated + partially ionized)\n");
    printf("    This is where real-gas corrections matter most\n");

    // ═══════════════════════════════════════════════════════════════
    // TEST 4: Jupiter return aerocapture (14.7 km/s)
    // ═══════════════════════════════════════════════════════════════
    printf("\n\n═══ TEST 4: EARTH REENTRY — JUPITER RETURN (14.7 km/s, 70 km) ═══\n\n");

    double v_bucket = 14700;
    double M_bucket = v_bucket / a_air(T_70);
    double q_bucket = sutton_graves_earth(rho_70, v_bucket, R_nose_capsule);
    double T_shock_bucket_pg = T_70 * (2.0 * gamma_air * M_bucket * M_bucket - (gamma_air - 1.0))
                              * ((gamma_air - 1.0) * M_bucket * M_bucket + 2.0)
                              / ((gamma_air + 1.0) * (gamma_air + 1.0) * M_bucket * M_bucket);

    printf("  Velocity:          14.7 km/s (Mach %.0f) — Bucket return from Jupiter\n", M_bucket);
    printf("  Stagnation heat:   %.0f kW/m² (Sutton-Graves)\n", q_bucket);
    printf("  Perfect gas shock: T=%.0f K\n", T_shock_bucket_pg);
    printf("  h_stagnation:      %.1f MJ/kg (v²/2)\n", 0.5 * v_bucket * v_bucket / 1e6);
    printf("\n");
    printf("  At this speed, air behind the shock is FULLY ionized plasma.\n");
    printf("  Perfect gas overpredicts T by ~3x, underpredicts ρ ratio by ~2x.\n");
    printf("  Real gas: T_shock ~ 15000-20000 K (energy goes into dissociation+ionization)\n");
    printf("  This is the regime where stat_thermo is essential.\n");

    // ═══════════════════════════════════════════════════════════════
    // TEST 5: Use our H2/He stat_thermo at Jupiter conditions for comparison
    // ═══════════════════════════════════════════════════════════════
    printf("\n\n═══ TEST 5: JUPITER ENTRY — REAL GAS H2/He (42 km/s) ═══\n\n");

    double rho_j, p_j, T_j;
    atmosphere_jupiter(rho_j, p_j, T_j, 350000.0);  // 350 km above 1-bar

    double v_jupiter = 42000;
    GasMixture pre = equilibrium_state(T_j, p_j, 0.862, 0.138);
    double M_jup = v_jupiter / pre.a_sound;

    printf("  Freestream:        T=%.0f K, ρ=%.2e kg/m³, a=%.0f m/s\n",
           T_j, rho_j, pre.a_sound);
    printf("  Velocity:          42 km/s (Mach %.0f)\n", M_jup);
    printf("  h_stagnation:      %.0f MJ/kg\n", (pre.h + 0.5*v_jupiter*v_jupiter)/1e6);

    ShockState shock = normal_shock(v_jupiter, T_j, p_j, 0.862, 0.138);
    printf("  Post-shock (real gas):\n");
    printf("    T = %.0f K\n", shock.post.T);
    printf("    ρ₂/ρ₁ = %.1f (vs perfect gas %.1f)\n", shock.density_ratio,
           (JUPITER_GAMMA+1)*M_jup*M_jup/((JUPITER_GAMMA-1)*M_jup*M_jup+2));
    printf("    γ_eff = %.3f\n", shock.post.gamma_eff);
    printf("    Composition: H₂=%.1f%% H=%.1f%% H⁺=%.1f%% He=%.1f%% He⁺=%.1f%% e⁻=%.1f%%\n",
           shock.post.X[H2]*100, shock.post.X[H]*100, shock.post.X[Hp]*100,
           shock.post.X[He]*100, shock.post.X[Hep]*100, shock.post.X[Elec]*100);

    // ═══════════════════════════════════════════════════════════════
    // SUMMARY TABLE
    // ═══════════════════════════════════════════════════════════════
    printf("\n\n╔══════════════════════════════════════════════════════════════════════╗\n");
    printf("║  COMPARISON: AIRLINE CRUISE → JUPITER ENTRY                        ║\n");
    printf("╠══════════════════════════════════════════════════════════════════════╣\n");
    printf("║  %-18s %-10s %-10s %-10s %-12s     ║\n",
           "Scenario", "v(km/s)", "Mach", "T_shock(K)", "q(kW/m²)");
    printf("║  ──────────────────────────────────────────────────────────         ║\n");
    printf("║  %-18s %-10.2f %-10.2f %-10.0f %-12s     ║\n",
           "FL390 cruise", v_ms/1000, Mach, TAT_K, "N/A (subsonic)");
    printf("║  %-18s %-10.1f %-10.0f %-10.0f %-12.0f     ║\n",
           "LEO reentry", v_leo/1000.0, M_leo, T_shock_pg, q_leo);
    printf("║  %-18s %-10.1f %-10.0f %-10.0f %-12.0f     ║\n",
           "Lunar return", v_lunar/1000.0, M_lunar, T_shock_lunar_pg, q_lunar);
    printf("║  %-18s %-10.1f %-10.0f %-10.0f %-12.0f     ║\n",
           "Bucket return", v_bucket/1000.0, M_bucket, T_shock_bucket_pg, q_bucket);
    printf("║  %-18s %-10.1f %-10.0f %-10.0f %-12s     ║\n",
           "Jupiter scoop", v_jupiter/1000.0, M_jup, shock.post.T, "mag shield");
    printf("║                                                                    ║\n");
    printf("║  NOTE: T_shock values for Earth entries are PERFECT GAS estimates. ║\n");
    printf("║  Real gas T is 2-3x lower (energy absorbed by dissociation).       ║\n");
    printf("║  Jupiter entry uses our stat_thermo real-gas model.                ║\n");
    printf("╚══════════════════════════════════════════════════════════════════════╝\n");

    return 0;
}
