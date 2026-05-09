// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// SATURN STATION — Crewed Titan Base + Atmospheric Mining Fleet
//
// Concept: Permanent crewed base on Titan with a fleet of simple scoop
// drones doing continuous atmospheric passes. Process and stockpile He/H2
// at the station, ship bulk cargo back to Earth.
//
// The Bucket is a truck. This is a refinery.
//
// Jupiter is hostile — 500 rem/day at Io, 36 rem/day at Europa.
// Electronics die in weeks. No crew. Every vehicle is on a clock.
//
// Saturn is benign — 0.05 rem/day. Safe for decades. Titan has a
// 1.47 bar N2/CH4 atmosphere you can land in with a parachute.
// Enceladus has free water. You can BUILD here.
//
// The 6.1% He fraction is fine when you're mining continuously.
//
// Compile: g++ -std=c++11 -O2 -I../../PlanetaryModels -o saturn_station saturn_station.cpp
// Run:     ./saturn_station
///////////////////////////////////////////////////////////////////////////////

#include "planetary_environment.hpp"
#include <cstdio>
#include <cmath>

static const double MU0 = 4.0 * M_PI * 1e-7;

double mag_funnel_radius(double B, double r0, double p_ram, double max_r)
{
    double B_need = sqrt(2.0 * MU0 * p_ram);
    if (B_need <= 0) return max_r;
    double R = r0 * pow(B / B_need, 1.0/3.0);
    return (R > max_r) ? max_r : R;
}

int main()
{
    printf("===================================================================\n");
    printf("  SATURN STATION — Titan Base + Atmospheric Mining Fleet\n");
    printf("  Kepler Workstation\n");
    printf("===================================================================\n");

    // Saturn orbital mechanics (km)
    double GM_S = 3.7931e7;
    double R_S  = 60268.0;

    // Atmosphere composition (mass fractions)
    double H2_FRAC  = 0.906;
    double HE_FRAC  = 0.061;
    double CH4_FRAC = 0.034;

    // =====================================================================
    // Section 1: The Architecture
    // =====================================================================
    printf("\n--- Section 1: Saturn Station Architecture ---\n\n");

    printf("  TITAN BASE (surface, 94 K, 1.47 bar)\n");
    printf("    - Crew habitat: 6-12 personnel, ISRU life support\n");
    printf("    - He/H2/CH4 separation & cryogenic storage\n");
    printf("    - Power: 500 kW fission reactor (Kilopower x5)\n");
    printf("    - Landing/launch via Titan atmosphere (parachute/propulsive)\n");
    printf("    - Titan surface g: 1.35 m/s2, escape: 2.64 km/s\n");
    printf("    - Titan aerobrake: FREE (1.47 bar atmosphere)\n\n");

    printf("  SCOOP DRONES (fleet of simple atmospheric miners)\n");
    printf("    - No interplanetary capability needed\n");
    printf("    - Just: magnet + cryo tank + Hall thruster + avionics\n");
    printf("    - Self-fueling: burn scooped H2 for orbit maintenance\n\n");

    printf("  TENDER (shuttle between scoop orbit and Titan)\n");
    printf("    - Picks up full tanks at scoop orbit apoapsis\n");
    printf("    - Delivers to Titan depot for processing\n\n");

    printf("  BULK CARGO (Earth return)\n");
    printf("    - Simple cryo tanker, no crew\n");
    printf("    - Departs Saturn, aerocaptures at Earth\n");

    // =====================================================================
    // Section 2: Scoop Drone Design
    // =====================================================================
    printf("\n--- Section 2: Scoop Drone Design ---\n\n");

    double drone_dry  = 2000;    // kg (magnet, structure, avionics, thruster)
    double drone_tank = 10000;   // kg capacity
    double drone_B    = 20.0;    // Tesla
    double drone_r    = 3.0;     // m coil radius
    double drone_eta  = 0.5;     // collection efficiency
    double drone_isp  = 3000;    // Hall thruster
    double max_R      = 200.0;

    printf("  Dry mass:     %g kg\n", drone_dry);
    printf("  Tank:         %g kg\n", drone_tank);
    printf("  Magnet:       %gT, %gm coil\n", drone_B, drone_r);
    printf("  Propulsion:   Hall thruster, Isp=%d s (burns scooped H2)\n", (int)drone_isp);
    printf("  No heat shield, no interplanetary systems, no VEEJ trajectory\n");
    printf("  Just a magnet and a bucket.\n");

    // =====================================================================
    // Section 3: Altitude Optimization for Continuous Ops
    // =====================================================================
    printf("\n--- Section 3: Scoop Altitude Trade Study ---\n\n");

    printf("  %-6s %8s %10s %8s %8s %8s %10s %8s\n",
           "Alt", "v_peri", "rho", "R_mag", "kg/pass", "passes", "days/fill", "t/yr");
    printf("  %-6s %8s %10s %8s %8s %8s %10s %8s\n",
           "(km)", "(km/s)", "(kg/m3)", "(m)", "", "to fill", "", "");

    struct AltResult {
        double alt, v_peri, rho, R_mag, kg_pass, T_orb;
        double passes_fill, days_fill, tonnes_yr;
    };
    AltResult results[10];
    int n_results = 0;

    double best_tpy = 0;
    double best_alt = 0;

    double alts[] = {500, 700, 1000, 1200, 1500, 1800, 2000, 2500};
    for (double alt : alts) {
        double r_peri = R_S + alt;
        double r_apo  = R_S * 5;
        double a_orb  = (r_peri + r_apo) / 2;
        double T_orb  = 2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_S);
        double v_p    = sqrt(GM_S*(2.0/r_peri - 1.0/a_orb)) * 1000; // m/s

        double rho, press, tempk;
        atmosphere_saturn(rho, press, tempk, alt*1000);

        double p_ram = 0.5 * rho * v_p * v_p;
        double Rm = mag_funnel_radius(drone_B, drone_r, p_ram, max_R);
        double Ae = M_PI * Rm * Rm;
        double m_dot = rho * v_p * Ae * drone_eta;

        double H_m = SATURN_RGAS * tempk / SATURN_G0;
        double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
        double dt = arc / v_p;
        double kg_pass = m_dot * dt;

        double passes = drone_tank / kg_pass;
        if (passes < 1) passes = 1;
        double days_fill = passes * T_orb / 86400;
        double fills_yr = 365.25 / days_fill;
        double tonnes_yr = fills_yr * drone_tank / 1000;

        printf("  %-6.0f %8.1f %10.2e %8.0f %8.1f %8.0f %10.1f %8.0f\n",
               alt, v_p/1000, rho, Rm, kg_pass, ceil(passes), days_fill, tonnes_yr);

        if (n_results < 10) {
            results[n_results] = {alt, v_p/1000, rho, Rm, kg_pass, T_orb,
                                  passes, days_fill, tonnes_yr};
            n_results++;
        }

        // Constraint: kg/pass must be < 10% of drone mass for manageable drag
        // At low altitudes, scooping >> vehicle mass is unphysical
        double dm_ratio = kg_pass / (drone_dry + drone_tank);
        if (tonnes_yr > best_tpy && dm_ratio < 0.10) {
            best_tpy = tonnes_yr; best_alt = alt;
        }
    }

    printf("\n  → Best continuous mining altitude: %.0f km (%.0f tonnes/yr/drone)\n",
           best_alt, best_tpy);

    // Pick operating altitude
    // Use best_alt but also consider vehicle lifetime (lower alt = more drag)
    // For crewed ops, prefer slightly higher alt for safety
    double op_alt = best_alt;
    // Find the result for this altitude
    AltResult op = results[0];
    for (int i = 0; i < n_results; i++) {
        if (results[i].alt == op_alt) op = results[i];
    }

    // =====================================================================
    // Section 4: Fleet Operations
    // =====================================================================
    printf("\n--- Section 4: Fleet Mining Operations ---\n\n");

    int n_drones_options[] = {5, 10, 20, 50};

    printf("  Operating altitude: %.0f km | %.1f kg/pass | %.0f days/fill\n",
           op.alt, op.kg_pass, op.days_fill);
    printf("  Orbit period: %.1f hours\n\n", op.T_orb/3600);

    printf("  %-8s %10s %10s %10s %10s %10s\n",
           "Drones", "t_raw/yr", "He(t/yr)", "H2(t/yr)", "CH4(t/yr)", "He value");
    printf("  %-8s %10s %10s %10s %10s %10s\n",
           "------", "--------", "--------", "--------", "---------", "--------");

    for (int nd : n_drones_options) {
        double raw = nd * best_tpy;
        double he  = raw * HE_FRAC;
        double h2  = raw * H2_FRAC;
        double ch4 = raw * CH4_FRAC;
        double he_val = he * 1000 * 1500 / 1e9; // $B at $1500/kg
        printf("  %-8d %10.0f %10.1f %10.0f %10.0f    $%.1fB/yr\n",
               nd, raw, he, h2, ch4, he_val);
    }

    // Baseline fleet: 20 drones
    int fleet_size = 20;
    double annual_raw = fleet_size * best_tpy;
    double annual_he  = annual_raw * HE_FRAC;
    double annual_h2  = annual_raw * H2_FRAC;
    double annual_ch4 = annual_raw * CH4_FRAC;

    printf("\n  BASELINE: %d drones at %.0f km altitude\n", fleet_size, op.alt);
    printf("    Raw material:   %.0f tonnes/year\n", annual_raw);
    printf("    Helium-3/4:     %.0f tonnes/year\n", annual_he);
    printf("    Hydrogen:       %.0f tonnes/year\n", annual_h2);
    printf("    Methane:        %.0f tonnes/year (also from Titan directly)\n", annual_ch4);

    // Orbit maintenance: atmospheric drag during each scoop pass lowers
    // apoapsis. Re-boost burn at apoapsis (Oberth advantage) restores orbit.
    //
    // Drag impulse per pass: dm * v_relative (momentum of captured gas)
    // where v_rel = v_orbital - v_corotation (~22 km/s at Saturn)
    // But: reboost at apoapsis costs only dv * (r_peri/r_apo) due to
    // the orbit mechanics — a small burn high up restores periapsis speed.
    double dm_pass = op.kg_pass;
    double M_drone = drone_dry + drone_tank/2;
    double omega_sat = 2*M_PI / (10.7 * 3600); // Saturn rotation rate
    double v_corot = omega_sat * (R_S + op.alt); // co-rotation velocity (km/s)
    double v_rel = op.v_peri - v_corot;
    double dv_drag = dm_pass * v_rel / (M_drone + dm_pass); // raw momentum loss
    double rp_ra = (R_S + op.alt) / (R_S * 5); // periapsis/apoapsis ratio
    double dv_reboost = dv_drag * rp_ra; // apoapsis reboost (Oberth saves fuel)

    double ve_hall = drone_isp * 9.80665e-3;
    double h2_per_reboost = M_drone * (1 - exp(-dv_reboost / ve_hall));
    double passes_per_fill = ceil(op.passes_fill);
    double h2_fuel_per_fill = h2_per_reboost * passes_per_fill;
    double h2_scooped = drone_tank * H2_FRAC;
    double net_h2 = h2_scooped - h2_fuel_per_fill;

    printf("\n  Orbit maintenance (self-fueling):\n");
    printf("    Drag per pass:    %.2f km/s (momentum transfer)\n", dv_drag);
    printf("    Reboost at apo:   %.3f km/s (r_peri/r_apo = %.2f)\n",
           dv_reboost, rp_ra);
    printf("    H2 fuel per fill: %.0f kg (%.0f passes)\n",
           h2_fuel_per_fill, passes_per_fill);
    printf("    H2 scooped:       %.0f kg\n", h2_scooped);
    printf("    Net H2 per fill:  %.0f kg (%.0f%% retained)\n",
           net_h2, net_h2/h2_scooped*100);
    printf("    → Self-sustaining: burns scooped H2, never needs refueling\n");

    // =====================================================================
    // Section 5: Processing & Staging
    // =====================================================================
    printf("\n--- Section 5: Titan Processing Station ---\n\n");

    printf("  SEPARATION PLANT\n");
    printf("    Input:  %.0f tonnes/year raw (%.0f%% H2, %.0f%% He, %.0f%% CH4)\n",
           annual_raw, H2_FRAC*100, HE_FRAC*100, CH4_FRAC*100);
    printf("    Process: cryogenic distillation (H2 boils at 20K, He at 4K)\n");
    printf("    Output:\n");
    printf("      He:   %.0f tonnes/year → cryogenic storage\n", annual_he);
    printf("      H2:   %.0f tonnes/year → propellant + fuel cells\n", annual_h2);
    printf("      CH4:  %.0f tonnes/year → propellant (+ unlimited from Titan lakes)\n",
           annual_ch4);

    // Titan surface resources
    printf("\n  TITAN SURFACE RESOURCES\n");
    printf("    Methane: Ligeia Mare alone = 9,000 km2 of liquid CH4\n");
    printf("             Estimated volume: ~2,000 km3 = 0.9 trillion tonnes\n");
    printf("             Just pump it. No mining, no drilling.\n");
    printf("    Nitrogen: 95%% of atmosphere (1.47 bar) → life support\n");
    printf("    Water ice: subsurface, available for mining\n");

    // Enceladus water
    printf("\n  ENCELADUS WATER\n");
    printf("    Plume output: ~200 kg/s continuously (Cassini measurement)\n");
    printf("    = 6,300 tonnes/year flowing into space FOR FREE\n");
    printf("    Electrolyze: H2O → H2 + O2\n");
    printf("    O2 for crew life support + LOX/LH2 propellant\n");
    printf("    Collection: fly through plume, passive scoop\n");
    printf("    dV to reach Enceladus from Titan: ~2 km/s\n");

    // =====================================================================
    // Section 6: Earth Return
    // =====================================================================
    printf("\n--- Section 6: Earth Return ---\n\n");

    // Bulk cargo tanker: giant cryo ship, never lands anywhere.
    // Sized to match transit time — 5 years each way, 10 year round trip.
    // Like oil tankers: huge because the voyage is long.
    // H2 propellant is free (scooped at Saturn).

    double dv_needed = 5.0; // km/s (Saturn escape + Earth transfer)
    double cargo_isp = 3000.0;
    double cargo_ve = cargo_isp * 9.80665e-3; // 29.4 km/s
    double cargo_mr = exp(dv_needed / cargo_ve); // 1.185

    // Size the tanker to carry 1000 tonnes He per trip
    // (enough to justify the 10-year round trip)
    double cargo_he = 1000e3;   // 1000 tonnes He payload
    double cargo_dry = 50e3;    // 50 tonnes structure + cryo + avionics
    double cargo_mf = cargo_dry + cargo_he;
    double cargo_m0 = cargo_mf * cargo_mr;
    double cargo_h2 = cargo_m0 - cargo_mf;

    printf("  BULK CARGO TANKER (never lands, sized for 10yr round trip)\n");
    printf("    Dry mass:     %.0f tonnes (structure, cryo insulation, avionics)\n", cargo_dry/1000);
    printf("    He payload:   %.0f tonnes\n", cargo_he/1000);
    printf("    H2 propellant: %.0f tonnes (scooped at Saturn — free)\n", cargo_h2/1000);
    printf("    Total mass:   %.0f tonnes\n", cargo_m0/1000);
    printf("    dV budget:    %.1f km/s (Isp=%.0fs, mass ratio=%.3f)\n",
           dv_needed, cargo_isp, cargo_mr);
    printf("    Propulsion:   Ion/Hall cluster, burns scooped H2\n");
    printf("    Earth arrival: aerocapture (free)\n");
    printf("    Transit:      ~5 years each way, 10yr round trip\n");
    printf("    Supplies:     carries crew supplies TO Saturn on return leg\n\n");

    // Pipeline: ships in transit at any time
    double transit_yr = 5.0;
    double round_trip = 2 * transit_yr;
    double ships_per_yr = annual_he / (cargo_he / 1000);
    int pipeline_depth = (int)ceil(ships_per_yr * round_trip);
    printf("  Ships needed:   %.1f departures/year\n", ships_per_yr);
    printf("  Pipeline depth: %d ships in transit at any time\n", pipeline_depth);
    printf("  Total fleet:    %d tankers (pipeline + 2 at Saturn loading)\n",
           pipeline_depth + 2);

    double he_to_earth = ships_per_yr * cargo_he / 1000;
    printf("  He to Earth:    %.0f tonnes/year\n", he_to_earth);

    // =====================================================================
    // Section 7: vs Jupiter Bucket Fleet
    // =====================================================================
    printf("\n--- Section 7: Saturn Station vs Jupiter Bucket Fleet ---\n\n");

    // Jupiter fleet: 10 Buckets, 10.5 year loop, 8.2t He per trip
    double j_fleet = 10;
    double j_loop = 10.5;
    double j_he_per_trip = 8.25; // tonnes
    double j_annual_he = j_fleet / j_loop * j_he_per_trip;

    printf("  %-32s %12s %12s\n", "", "Jupiter", "Saturn");
    printf("  %-32s %12s %12s\n", "", "Bucket x10", "Station");
    printf("  %-32s %12s %12s\n", "", "----------", "--------");
    printf("  %-32s %12s %12s\n", "Architecture", "Unmanned", "Crewed base");
    printf("  %-32s %12s %12s\n", "Location", "Flyby", "Titan surface");
    printf("  %-32s %12s %12s\n", "Operations", "Get in/out", "Permanent");
    printf("  %-32s %11.0f  %11d\n", "Mining vehicles", j_fleet, fleet_size);
    printf("  %-32s %11.1f  %11.0f t/yr\n", "He delivered to Earth", j_annual_he, he_to_earth);
    printf("  %-32s %12s %12s\n", "Crewed?", "No", "Yes (6-12)");
    printf("  %-32s %12s %12s\n", "Radiation", "Lethal", "Safe");
    printf("  %-32s %12s %12s\n", "O2 source", "Europa ice", "Enceladus");
    printf("  %-32s %12s %12s\n", "CH4 source", "None", "Titan lakes");
    printf("  %-32s %12s %12s\n", "Vehicle reuse", "Maybe", "Indefinite");
    printf("  %-32s %12s %12s\n", "Expandable?", "Linear", "Add drones");

    double ratio = he_to_earth / j_annual_he;
    printf("\n  He delivery ratio: Saturn / Jupiter = %.0fx\n", ratio);
    printf("  Saturn Station delivers %.0fx more He per year.\n\n", ratio);

    printf("  Why? Because:\n");
    printf("    Jupiter Bucket: 1 trip per vehicle per 10.5 years\n");
    printf("    Saturn Drone:   %.0f fills per year, continuous\n",
           365.25 / op.days_fill);
    printf("    The 6.1%% He fraction doesn't matter when you mine\n");
    printf("    %.0f tonnes of atmosphere per year per drone.\n", best_tpy);

    // =====================================================================
    // Section 8: Economics
    // =====================================================================
    printf("\n--- Section 8: Economics ---\n\n");

    // Upfront costs (rough)
    double cost_titan_base = 50e9;     // $50B (includes transit)
    double cost_per_drone  = 200e6;    // $200M each
    double cost_cargo_ship = 500e6;    // $500M each
    double cost_enceladus  = 5e9;      // $5B for water collection setup
    int cargo_ships_needed = pipeline_depth + 2;

    double total_upfront = cost_titan_base
                         + fleet_size * cost_per_drone
                         + cargo_ships_needed * cost_cargo_ship
                         + cost_enceladus;

    double annual_he_revenue = he_to_earth * 1000 * 1500; // $1500/kg
    double annual_h2_revenue = annual_h2 * 1000 * 500;    // local use credit
    double payback_years = total_upfront / annual_he_revenue;

    printf("  Upfront investment:\n");
    printf("    Titan base + transit:      $%.0fB\n", cost_titan_base/1e9);
    printf("    %d scoop drones:            $%.0fB\n", fleet_size, fleet_size*cost_per_drone/1e9);
    printf("    %d cargo ships:             $%.0fB\n", cargo_ships_needed, cargo_ships_needed*cost_cargo_ship/1e9);
    printf("    Enceladus water ops:       $%.0fB\n", cost_enceladus/1e9);
    printf("    TOTAL:                     $%.0fB\n", total_upfront/1e9);
    printf("\n");
    printf("  Annual revenue (He at $1500/kg):\n");
    printf("    He to Earth: %.0f tonnes x $1.5M/t = $%.0fB/year\n",
           he_to_earth, annual_he_revenue/1e9);
    printf("    Payback:     %.1f years (from first delivery)\n", payback_years);

    // Compare Jupiter
    double j_upfront = j_fleet * 245e6; // $245M per bucket
    double j_annual_rev = j_annual_he * 1000 * 1500;
    double j_payback = j_upfront / j_annual_rev;

    printf("\n  Compare Jupiter Bucket fleet:\n");
    printf("    Upfront: $%.1fB (10 buckets at $245M)\n", j_upfront/1e9);
    printf("    Revenue: $%.0fM/year (%.1f t He)\n", j_annual_rev/1e6, j_annual_he);
    printf("    Payback: %.1f years\n", j_payback);

    printf("\n  Jupiter is CHEAPER to start.\n");
    printf("  Saturn is BIGGER when running.\n");
    printf("  Build Jupiter first. Fund Saturn with Jupiter He revenue.\n");

    // =====================================================================
    // Section 9: The 30-Year Plan
    // =====================================================================
    printf("\n--- Section 9: The 30-Year Plan ---\n\n");

    printf("  Phase 1 (Years 1-10): JUPITER BUCKETS\n");
    printf("    Launch 10 Buckets via VEEJ\n");
    printf("    First He delivery: Year 10.5\n");
    printf("    Revenue: ~$%.0fM/year starting Year 11\n", j_annual_rev/1e6);
    printf("    Fund Saturn Station development\n\n");

    printf("  Phase 2 (Years 10-17): BUILD SATURN STATION\n");
    printf("    Year 10: Launch Titan base (VVEJGA, 7yr transit)\n");
    printf("    Year 12: Launch drone fleet (20 drones)\n");
    printf("    Year 14: Launch cargo ships\n");
    printf("    Year 17: Titan base operational, first drones arrive\n\n");

    printf("  Phase 3 (Years 17-22): RAMP UP\n");
    printf("    Drones begin continuous mining\n");
    printf("    Process and stockpile He at Titan\n");
    printf("    First bulk cargo departs Saturn: Year 19\n");
    printf("    First Saturn He reaches Earth: Year 24\n\n");

    printf("  Phase 4 (Years 24+): FULL OPERATIONS\n");
    printf("    Jupiter: ~%.0f tonnes He/year (steady state)\n", j_annual_he);
    printf("    Saturn:  ~%.0f tonnes He/year (scaling up)\n", he_to_earth);
    printf("    Combined: ~%.0f tonnes He/year\n", j_annual_he + he_to_earth);
    printf("    Earth He consumption: ~30,000 tonnes/year\n");
    printf("    Fraction replaced: %.1f%%\n",
           (j_annual_he + he_to_earth) / 30000 * 100);
    printf("    At 100 Saturn drones: %.0f%% of Earth consumption\n",
           (j_annual_he + 100.0/fleet_size * he_to_earth) / 30000 * 100);

    // =====================================================================
    // Conclusion
    // =====================================================================
    printf("\n===================================================================\n");
    printf("  CONCLUSION\n");
    printf("===================================================================\n");
    printf("\n");
    printf("  Jupiter is where you START. Unmanned. Proven. Fast payback.\n");
    printf("  Saturn is where you SCALE. Crewed. Permanent. Industrial.\n");
    printf("\n");
    printf("  The 6.1%% He fraction is irrelevant at scale.\n");
    printf("  When you mine %.0f tonnes/year with %d drones, you get\n",
           annual_raw, fleet_size);
    printf("  %.0f tonnes of He — more than %d Jupiter Buckets.\n",
           annual_he, (int)(annual_he / j_he_per_trip));
    printf("\n");
    printf("  And Saturn gives you everything else:\n");
    printf("    - Titan: CH4 lakes, N2 atmosphere, landing site, crew base\n");
    printf("    - Enceladus: free water geysers → O2 + H2\n");
    printf("    - Low radiation: safe for humans and electronics\n");
    printf("    - Expandable: add drones, not missions\n");
    printf("\n");
    printf("  Jupiter is hostile. You raid it.\n");
    printf("  Saturn is habitable. You settle it.\n");
    printf("===================================================================\n");

    return 0;
}
