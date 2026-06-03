// Test and validation for stat_thermo.hpp
// Sweeps temperature from 200 K to 200,000 K and prints composition/properties.
// Validates cold-gas limit against perfect-gas constants.
//
// Compile: g++ -std=c++11 -O2 -o test_stat_thermo test_stat_thermo.cpp
// Run:     ./test_stat_thermo

#include "stat_thermo.hpp"
#include <cstdio>
#include <cmath>

using namespace stat_thermo;

int main()
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  STATISTICAL THERMODYNAMICS — 7-SPECIES H2/He VALIDATION   ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    // Jupiter composition: 86.2% H2, 13.8% He by mass
    double Y_H2 = 0.862;
    double Y_He = 0.138;
    double p = 1.0e5;  // 1 bar

    // ── Test 1: Cold gas limit ──
    printf("=== TEST 1: Cold Gas Limit (T=166 K, p=1 bar) ===\n");
    printf("Expected: γ ≈ 1.44, M ≈ 2.22 g/mol, R ≈ 3745 J/(kg·K)\n\n");

    GasMixture cold = equilibrium_state(166.0, p, Y_H2, Y_He);
    print_mixture(cold, "T = 166 K (Jupiter 1-bar)");

    double gamma_err = fabs(cold.gamma_eff - 1.44) / 1.44 * 100;
    double M_err = fabs(cold.M_mix * 1000 - 2.22) / 2.22 * 100;
    double R_err = fabs(cold.R_mix - 3745) / 3745 * 100;
    printf("\n  γ error: %.2f%%  M error: %.2f%%  R error: %.2f%%\n", gamma_err, M_err, R_err);
    printf("  %s\n\n", (gamma_err < 5 && M_err < 5 && R_err < 5) ?
           "PASS — matches perfect gas" : "CHECK — deviation from perfect gas");

    // ── Test 2: Temperature sweep ──
    printf("=== TEST 2: Temperature Sweep (200 K — 200,000 K at 1 bar) ===\n\n");
    printf("%-10s %-8s %-8s %-10s %-10s %-8s %-8s %-8s %-8s %-8s %-8s\n",
           "T(K)", "γ_eff", "M(g/mol)", "Cp(J/kgK)", "a(m/s)",
           "X_H2", "X_H", "X_H+", "X_He", "X_He+", "X_e-");
    printf("────────────────────────────────────────────────────────────────"
           "──────────────────────────────────────────────\n");

    double temps[] = {200, 300, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000,
                      6000, 8000, 10000, 15000, 20000, 30000, 50000, 100000, 200000};
    int n_temps = sizeof(temps) / sizeof(temps[0]);

    for (int i = 0; i < n_temps; i++) {
        GasMixture m = equilibrium_state(temps[i], p, Y_H2, Y_He);
        printf("%-10.0f %-8.4f %-8.4f %-10.1f %-10.1f %-8.2e %-8.2e %-8.2e %-8.2e %-8.2e %-8.2e\n",
               m.T, m.gamma_eff, m.M_mix*1000, m.Cp, m.a_sound,
               m.X[H2], m.X[H], m.X[Hp], m.X[He], m.X[Hep], m.X[Elec]);
    }

    // ── Test 3: Dissociation zone detail ──
    printf("\n=== TEST 3: Dissociation Zone (1000-6000 K) ===\n\n");
    printf("%-8s %-10s %-10s %-10s %-8s %-8s\n",
           "T(K)", "α_dissoc", "X_H2", "X_H", "γ_eff", "M(g/mol)");
    printf("──────────────────────────────────────────────────────────\n");

    for (double T = 1000; T <= 6000; T += 250) {
        GasMixture m = equilibrium_state(T, p, Y_H2, Y_He);
        printf("%-8.0f %-10.6f %-10.6f %-10.6f %-8.4f %-8.4f\n",
               T, m.alpha_d, m.X[H2], m.X[H], m.gamma_eff, m.M_mix*1000);
    }

    // ── Test 4: Normal shock at Jupiter conditions ──
    printf("\n=== TEST 4: Normal Shock — Jupiter Entry ===\n");
    printf("Freestream: 350 km above 1-bar, v = 42 km/s\n");

    // Get freestream conditions at 350 km
    // (use approximate values from planetary_environment.hpp)
    double T_inf = 120.0;   // K (stratosphere)
    double p_inf = 0.1;     // Pa (very thin, 350 km up)
    double v_inf = 42000.0; // m/s

    ShockState shock = normal_shock(v_inf, T_inf, p_inf, Y_H2, Y_He);
    print_shock(shock);

    // ── Test 5: Lower velocity shock (Mach reduction scenario) ──
    printf("\n=== TEST 5: Reduced Velocity Shock (v = 15 km/s, ~aerobraked) ===\n");
    double v_low = 15000.0;
    ShockState shock_low = normal_shock(v_low, T_inf, p_inf, Y_H2, Y_He);
    print_shock(shock_low);

    // ── Test 6: Pressure dependence ──
    printf("\n=== TEST 6: Pressure Dependence at T=5000 K ===\n");
    printf("(Dissociation shifts with pressure — Le Chatelier's principle)\n\n");
    printf("%-12s %-10s %-8s %-8s %-8s\n",
           "p(Pa)", "α_dissoc", "X_H2", "X_H", "γ_eff");
    printf("──────────────────────────────────────────────\n");

    double pressures[] = {1.0, 10.0, 100.0, 1000.0, 1e4, 1e5, 1e6};
    for (int i = 0; i < 7; i++) {
        GasMixture m = equilibrium_state(5000.0, pressures[i], Y_H2, Y_He);
        printf("%-12.1e %-10.6f %-8.4f %-8.4f %-8.4f\n",
               pressures[i], m.alpha_d, m.X[H2], m.X[H], m.gamma_eff);
    }

    // ── Summary ──
    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  KEY OBSERVATIONS                                          ║\n");
    printf("╠══════════════════════════════════════════════════════════════╣\n");

    GasMixture m300 = equilibrium_state(300, p, Y_H2, Y_He);
    GasMixture m5k  = equilibrium_state(5000, p, Y_H2, Y_He);
    GasMixture m20k = equilibrium_state(20000, p, Y_H2, Y_He);
    GasMixture m100k = equilibrium_state(100000, p, Y_H2, Y_He);

    printf("║  300 K:    γ=%.3f  M=%.2f g/mol  (cold H2/He)          ║\n",
           m300.gamma_eff, m300.M_mix*1000);
    printf("║  5000 K:   γ=%.3f  M=%.2f g/mol  (dissociating)        ║\n",
           m5k.gamma_eff, m5k.M_mix*1000);
    printf("║  20000 K:  γ=%.3f  M=%.2f g/mol  (ionizing)            ║\n",
           m20k.gamma_eff, m20k.M_mix*1000);
    printf("║  100000 K: γ=%.3f  M=%.2f g/mol  (fully ionized)       ║\n",
           m100k.gamma_eff, m100k.M_mix*1000);
    printf("║                                                            ║\n");
    printf("║  Perfect gas assumption (constant γ=1.44) is WRONG above   ║\n");
    printf("║  2000 K — which is exactly where the scoop operates.       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}
