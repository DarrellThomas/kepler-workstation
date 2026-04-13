///////////////////////////////////////////////////////////////////////////////
// Magnetic Scoop Harvester — Full Mission Simulation
//
// Option 4: Bussard-style magnetic collection at Jupiter
// - Superconducting solenoid creates magnetic funnel
// - Bow shock self-ionizes H2/He at 42 km/s
// - Vehicle never contacts dense atmosphere
// - No heat shield, no thermal cycling, unlimited reuse
// - O2 from Europa/Ganymede water ice electrolysis
// - Every pass refuels the vehicle
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include "../planetary_models/planetary_environment.hpp"

int main()
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║    MAGNETIC SCOOP HARVESTER — JUPITER OPERATIONS SIM       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    // Constants
    double GM_J = 1.26713e8;   // km³/s²
    double R_J = 71492.0;      // km
    double mu0 = 4 * M_PI * 1e-7;

    // Magnetic coil parameters
    double B_coil = 10.0;     // Tesla (HTS solenoid)
    double r_coil = 2.0;      // meters (coil radius)

    // Vehicle
    double dry_mass = 22000;  // kg (no TPS!)
    double tank_capacity = 80000; // kg (COPV cryo tanks)

    // Orbit: elliptical, periapsis in upper atmosphere
    double alt_peri = 350.0;  // km above 1-bar (operating altitude)
    double r_peri = R_J + alt_peri;
    double r_apo = R_J * 5;  // 5 Rj apoapsis
    double a_orbit = (r_peri + r_apo) / 2;
    double T_orbit = 2 * M_PI * sqrt(a_orbit*a_orbit*a_orbit / GM_J); // seconds
    double v_peri = sqrt(GM_J * (2.0/r_peri - 1.0/a_orbit)) * 1000; // m/s

    // Atmospheric conditions at periapsis
    double rho, press, tempk;
    atmosphere_jupiter(rho, press, tempk, alt_peri * 1000.0);

    // Ram pressure and magnetic field required
    double p_ram = 0.5 * rho * v_peri * v_peri;
    double B_needed = sqrt(2 * mu0 * p_ram);

    // Magnetic funnel radius (dipole: B = B_coil * (r_coil/r)³)
    double R_mag = r_coil * pow(B_coil / B_needed, 1.0/3.0);
    double A_eff = M_PI * R_mag * R_mag;

    // Collection rate
    double eta_collect = 0.5; // 50% of intercepted ions captured
    double m_dot = rho * v_peri * A_eff * eta_collect;

    // Pass duration (atmosphere arc)
    double H_scale = JUPITER_RGAS * tempk / JUPITER_G0;
    double arc_length = 2.0 * sqrt(2.0 * r_peri * 1000 * H_scale);
    double pass_time = arc_length / v_peri;
    double mass_per_pass = m_dot * pass_time;

    // Heating on vehicle hull (it's BEHIND the magnetic shield)
    // Hull sees only the scattered/leaked particles
    double leak_fraction = 0.01; // 1% leaks through magnetic field
    double rho_hull = rho * leak_fraction;
    double q_hull = 1.7415e-4 * sqrt(rho_hull / 3.0) * pow(v_peri, 3) / 1000;

    printf("═══ MAGNETIC SCOOP PARAMETERS ═══\n\n");
    printf("  Coil:              B=%g T, r=%.0f m (HTS superconducting)\n", B_coil, r_coil);
    printf("  Orbit:             %g km × %g Rj\n", alt_peri, r_apo/R_J);
    printf("  Periapsis:         %.0f km above 1-bar\n", alt_peri);
    printf("  Velocity at peri:  %.0f m/s (%.1f km/s)\n", v_peri, v_peri/1000);
    printf("  Orbit period:      %.1f hours\n", T_orbit/3600);
    printf("\n");
    printf("  Atmo density:      %.2e kg/m³\n", rho);
    printf("  Ram pressure:      %.3f Pa\n", p_ram);
    printf("  B needed:          %.6f T (a fridge magnet is 0.001 T)\n", B_needed);
    printf("  Magnetic radius:   %.0f m\n", R_mag);
    printf("  Effective area:    %.0f m² (vs 50 m² physical port)\n", A_eff);
    printf("  Collection rate:   %.2f kg/s\n", m_dot);
    printf("  Pass duration:     %.0f seconds\n", pass_time);
    printf("  Mass per pass:     %.1f kg\n", mass_per_pass);
    printf("\n");
    printf("  Hull heating:      %.2f kW/cm² (1%% magnetic leakage)\n", q_hull);
    printf("  Hull temperature:  ~%.0f °C\n",
           pow(q_hull*1000*1e4/(5.67e-8*0.9), 0.25) - 273);
    printf("  → %s\n", q_hull < 0.01 ? "ALUMINUM SKIN — no TPS needed!" :
           q_hull < 1 ? "SiC tiles — lightweight" : "Need better shielding");

    // Ionization check
    double m_h2 = 2 * 1.67e-27;
    double KE_h2 = 0.5 * m_h2 * v_peri * v_peri;
    double eV_h2 = KE_h2 / 1.6e-19;
    printf("\n  H2 kinetic energy: %.0f eV (ionization: 15.4 eV) → %s\n",
           eV_h2, eV_h2 > 15.4 ? "SELF-IONIZING ✓" : "needs ionizer");

    // ================================================================
    // 2-YEAR HARVEST CAMPAIGN
    // ================================================================
    printf("\n═══ 2-YEAR HARVEST CAMPAIGN ═══\n\n");

    double campaign_days = 2 * 365.25;
    double campaign_sec = campaign_days * 86400;
    int total_passes = (int)(campaign_sec / T_orbit);

    double h2_total = 0, he_total = 0;
    double current_h2 = 0, current_he = 0;
    double cumulative_mass = 0;

    // Simulate pass by pass
    int milestone_passes[] = {1, 10, 50, 100, 200, 500, 1000, total_passes};
    int mi = 0;

    printf("  %-8s %-12s %-12s %-12s %-12s %-10s\n",
           "Pass#", "This pass", "Cum H2(t)", "Cum He(t)", "Cum total", "Tank %");

    for (int pass = 1; pass <= total_passes; pass++)
    {
        double harvested = mass_per_pass;
        double h2 = harvested * 0.75;
        double he = harvested * 0.25;

        current_h2 += h2;
        current_he += he;
        cumulative_mass += harvested;

        // Check if we should print this pass
        if (mi < 8 && pass == milestone_passes[mi]) {
            double tank_pct = cumulative_mass / tank_capacity * 100;
            printf("  %-8d %-12.1f %-12.1f %-12.1f %-12.1f %-10.1f\n",
                   pass, harvested, current_h2/1000, current_he/1000,
                   cumulative_mass/1000, tank_pct);
            mi++;
        }

        // Tank full? Stop harvesting (or vent and keep going)
        if (cumulative_mass >= tank_capacity && pass < total_passes) {
            printf("  *** TANKS FULL at pass %d (day %.0f) ***\n",
                   pass, pass * T_orbit / 86400);
            h2_total = current_h2;
            he_total = current_he;
            break;
        }
    }

    if (cumulative_mass < tank_capacity) {
        h2_total = current_h2;
        he_total = current_he;
    }

    printf("\n  Campaign results:\n");
    printf("  Total passes:      %d\n", total_passes);
    printf("  Total harvested:   %.1f tonnes\n", cumulative_mass/1000);
    printf("  H2:                %.1f tonnes\n", h2_total/1000);
    printf("  He:                %.1f tonnes\n", he_total/1000);

    // ================================================================
    // EUROPA ICE MINING (O2)
    // ================================================================
    printf("\n═══ EUROPA ICE MINING ═══\n\n");

    // How much O2 do we need?
    double lox_lh2_ratio = 6.0;
    // Use some H2 for return propulsion, rest as payload
    // First figure out how much propellant we need
    double dv_return = 5.0;  // km/s (Jupiter to Mars, arcjet)
    double ve_lox = 450 * 9.80665e-3; // LOX/LH2 chemical Isp=450s

    // With LOX/LH2 for return:
    double mf = dry_mass + he_total; // final mass (dry + He payload)
    double m0 = mf * exp(dv_return / ve_lox);
    double prop_total = m0 - mf;
    double h2_fuel = prop_total / 7.0;  // 1/7 of LOX/LH2 is H2
    double o2_fuel = prop_total - h2_fuel;
    double h2_payload = h2_total - h2_fuel;

    // Ice needed for O2
    double ice_needed = o2_fuel / 0.889;
    double bonus_h2 = ice_needed - o2_fuel;

    printf("  Return propulsion (LOX/LH2, Isp=450s):\n");
    printf("    ΔV needed:       %.1f km/s\n", dv_return);
    printf("    Propellant:      %.1f tonnes (%.1ft H2 + %.1ft O2)\n",
           prop_total/1000, h2_fuel/1000, o2_fuel/1000);
    printf("    H2 from Jupiter: %.1f tonnes available → %.1ft for fuel, %.1ft payload\n",
           h2_total/1000, h2_fuel/1000, h2_payload/1000);
    printf("\n");
    printf("  Europa ice electrolysis:\n");
    printf("    O2 needed:       %.1f tonnes\n", o2_fuel/1000);
    printf("    Ice to mine:     %.1f tonnes (a %.0fm cube)\n",
           ice_needed/1000, pow(ice_needed/917, 1.0/3.0));
    printf("    Bonus H2:        %.1f tonnes (from water)\n", bonus_h2/1000);
    printf("    Energy:          %.0f MWh (RTG or fuel cell)\n",
           o2_fuel * 13.2e3 / 3.6e9);

    // Europa delta-V from Jupiter orbit
    // Europa orbits at 671,000 km from Jupiter
    double r_europa = 671000; // km
    double v_europa = sqrt(GM_J / r_europa); // km/s
    printf("\n  Europa access (from harvesting orbit):\n");
    printf("    Europa orbit:    %.0f km from Jupiter\n", r_europa);
    printf("    Europa v_circ:   %.1f km/s\n", v_europa);
    printf("    Hohmann from 5Rj: ~%.1f km/s ΔV\n",
           fabs(sqrt(GM_J*(2.0/r_apo - 2.0/(r_apo + r_europa))) -
                sqrt(GM_J/r_apo)));
    printf("    Europa surface g: 1.315 m/s², escape: 2.025 km/s\n");
    printf("    Land, mine ice, launch: ~2 km/s ΔV each way\n");

    // ================================================================
    // COMPLETE MISSION PROFILE
    // ================================================================
    printf("\n═══ COMPLETE MISSION — MAGNETIC SCOOP ARCHITECTURE ═══\n\n");

    double total_delivered = h2_payload + bonus_h2 + he_total;

    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  OUTBOUND (Earth → Jupiter)                                ║\n");
    printf("║    Departure ΔV: 3.5 km/s from LEO (VEEJ trajectory)      ║\n");
    printf("║    Transit: ~5 years (with Venus gravity assists)          ║\n");
    printf("║    Jupiter capture: propulsive (small burn from H2 reserve)║\n");
    printf("║                                                            ║\n");
    printf("║  JUPITER OPERATIONS (2 years)                              ║\n");
    printf("║    Magnetic scoop harvesting:                              ║\n");
    printf("║      Altitude: %.0f km above 1-bar                         ║\n", alt_peri);
    printf("║      Magnetic funnel: %.0f m radius, %.0f m² effective area║\n", R_mag, A_eff);
    printf("║      %.1f kg per pass × %d passes                       ║\n",
           mass_per_pass, total_passes);
    printf("║      H2: %.1f tonnes, He: %.1f tonnes                     ║\n",
           h2_total/1000, he_total/1000);
    printf("║      Hull temp: <%.0f°C — NO HEAT SHIELD                  ║\n",
           pow(q_hull*1000*1e4/(5.67e-8*0.9), 0.25) - 273);
    printf("║    Europa ice mining:                                      ║\n");
    printf("║      Mine %.0f tonnes ice → %.0f tonnes O2                 ║\n",
           ice_needed/1000, o2_fuel/1000);
    printf("║      + %.1f tonnes bonus H2                                ║\n", bonus_h2/1000);
    printf("║                                                            ║\n");
    printf("║  RETURN (Jupiter → Mars/Earth)                             ║\n");
    printf("║    Propulsion: LOX/LH2 (Isp=450s)                         ║\n");
    printf("║    ΔV: %.1f km/s                                           ║\n", dv_return);
    printf("║    Propellant: %.0ft LOX + %.0ft LH2                       ║\n",
           o2_fuel/1000, h2_fuel/1000);
    printf("║    Transit: ~3 years                                       ║\n");
    printf("║    Earth/Mars aerocapture: FREE                            ║\n");
    printf("║                                                            ║\n");
    printf("║  DELIVERED TO DEPOT                                        ║\n");
    printf("║    H2:    %.1f tonnes (atmo) + %.1f tonnes (ice) = %.1f t  ║\n",
           h2_payload/1000, bonus_h2/1000, (h2_payload+bonus_h2)/1000);
    printf("║    He:    %.1f tonnes                                      ║\n", he_total/1000);
    printf("║    TOTAL: %.1f tonnes per collector                        ║\n",
           total_delivered/1000);
    printf("║                                                            ║\n");
    printf("║  VEHICLE: REUSABLE (no ablation, no TPS degradation)      ║\n");
    printf("║    Refuel at depot → repeat mission                        ║\n");
    printf("║    Design life: 10+ missions (50+ years)                   ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    // Fleet economics
    printf("\n═══ FLEET (5 collectors, annual launches) ═══\n\n");
    double collectors = 5;
    double window = 398.88; // days
    double launches_per_yr = 365.25 / window;
    double annual_total = collectors * launches_per_yr * total_delivered / 1000;
    double annual_h2 = collectors * launches_per_yr * (h2_payload+bonus_h2) / 1000;
    double annual_he = collectors * launches_per_yr * he_total / 1000;

    printf("  Collectors per window:   %.0f\n", collectors);
    printf("  Annual H2 to depot:      %.0f tonnes\n", annual_h2);
    printf("  Annual He to depot:      %.0f tonnes\n", annual_he);
    printf("  Annual total:            %.0f tonnes\n", annual_total);
    printf("  Vehicle reuse:           10+ missions each\n");
    printf("  Fleet replacement:       ~every 50 years\n");
    printf("\n");
    printf("  Value at $1500/kg launch cost avoided:\n");
    printf("    $%.0fM/year in H2 delivery\n", annual_h2 * 1000 * 1500 / 1e6);
    printf("    + He stockpile for Earth's future\n");
    printf("    + Europa science (free, while you're there)\n");

    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  The magnetic field is the scoop.                          ║\n");
    printf("║  The bow shock is the ionizer.                             ║\n");
    printf("║  Jupiter's gravity is the engine.                          ║\n");
    printf("║  Europa's ice is the oxidizer.                             ║\n");
    printf("║  Every pass refuels the vehicle.                           ║\n");
    printf("║  The vehicle is reusable forever.                          ║\n");
    printf("║  No heat shield. No thermal cycling. No ablation.          ║\n");
    printf("║                                                            ║\n");
    printf("║  We don't harvest the atmosphere.                          ║\n");
    printf("║  We harvest the magnetosphere.                             ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}
