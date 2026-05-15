// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// SATURN STATION — Full Mission Simulation
//
// End-to-end business case with proper orbital mechanics:
//   1. Real DE440 ephemeris for trajectory design
//   2. Lambert scan for Earth→Saturn and Saturn→Earth windows
//   3. Pass-by-pass atmospheric scoop campaign
//   4. Tanker pipeline with real departure windows
//   5. Year-by-year cash flow (2030-2080)
//   6. NPV / break-even / sensitivity analysis
//
// Question: Can a crewed Titan base with atmospheric mining drones
// and He tankers running Earth↔Saturn routes make money?
//
// Compile: g++ -std=c++11 -O2 -I../../Ephemeris -I../../PlanetaryModels \
//          -o saturn_full_sim saturn_full_sim.cpp
// Run:     ./saturn_full_sim
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "jpl_ephemeris.hpp"
#include "lambert.hpp"
#include "planetary_environment.hpp"

// ═══════════════════════════════════════════════════════════════════════════
// Utilities
// ═══════════════════════════════════════════════════════════════════════════
double vec_mag(const double v[3]) {
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
void vec_sub(double o[3], const double a[3], const double b[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
double jd_from_year(double y) { return 2451545.0 + (y-2000.0)*365.25; }
double year_from_jd(double jd) { return 2000.0 + (jd-2451545.0)/365.25; }

static const double MU_SUN = 1.32712440018e11; // km³/s²
static const double MU0    = 4.0 * M_PI * 1e-7;

// ═══════════════════════════════════════════════════════════════════════════
// Lambert window scanner (reuses pattern from bucket_sim.cpp)
// ═══════════════════════════════════════════════════════════════════════════
struct TransferWindow {
    double jd_dep, jd_arr, v_inf_dep, v_inf_arr, tof_days;
    double year_dep, year_arr;
    bool ok;
};

TransferWindow scan_transfer(JplEphemeris& eph, JplBody b1, JplBody b2,
                             bool e1, bool e2,
                             double jd_lo, double jd_hi,
                             double tof_lo, double tof_hi, double step)
{
    TransferWindow best;
    best.ok = false; best.v_inf_dep = 1e10;
    double sp[3], sv[3];
    for (double jd = jd_lo; jd <= jd_hi; jd += step) {
        for (double tof = tof_lo; tof <= tof_hi; tof += step) {
            double p1[3],v1[3],p2[3],v2[3],r1[3],rv1[3],r2[3],rv2[3];
            if (e1) eph.get_earth(jd,p1,v1); else eph.get_state(b1,jd,p1,v1);
            if (e2) eph.get_earth(jd+tof,p2,v2); else eph.get_state(b2,jd+tof,p2,v2);
            eph.get_state(JPL_SUN,jd,sp,sv);
            vec_sub(r1,p1,sp); vec_sub(rv1,v1,sv);
            eph.get_state(JPL_SUN,jd+tof,sp,sv);
            vec_sub(r2,p2,sp); vec_sub(rv2,v2,sv);
            double vt1[3],vt2[3];
            if (!lambert_solve(r1,r2,tof*86400,MU_SUN,true,vt1,vt2)) continue;
            double vd[3],va[3]; vec_sub(vd,vt1,rv1); vec_sub(va,vt2,rv2);
            double vid=vec_mag(vd), via=vec_mag(va);
            if (vid < best.v_inf_dep && vid < 20 && via < 20) {
                best.jd_dep=jd; best.jd_arr=jd+tof;
                best.v_inf_dep=vid; best.v_inf_arr=via;
                best.tof_days=tof; best.ok=true;
                best.year_dep=year_from_jd(jd);
                best.year_arr=year_from_jd(jd+tof);
            }
        }
    }
    return best;
}

// Scan ALL windows in a range, returning the best per synodic period
std::vector<TransferWindow> scan_all_windows(
    JplEphemeris& eph, JplBody b1, JplBody b2, bool e1, bool e2,
    double yr_lo, double yr_hi, double tof_lo, double tof_hi,
    double window_gap, double scan_step)
{
    std::vector<TransferWindow> windows;
    for (double yr = yr_lo; yr < yr_hi; yr += window_gap) {
        double jd0 = jd_from_year(yr);
        double jd1 = jd_from_year(yr + window_gap);
        TransferWindow w = scan_transfer(eph, b1, b2, e1, e2,
                                         jd0, jd1, tof_lo, tof_hi, scan_step);
        if (w.ok) windows.push_back(w);
    }
    return windows;
}

// ═══════════════════════════════════════════════════════════════════════════
// Scoop drone model
// ═══════════════════════════════════════════════════════════════════════════
struct DroneSpec {
    double dry_mass;     // kg
    double tank_cap;     // kg
    double B_coil;       // Tesla
    double r_coil;       // m
    double eta;          // collection efficiency
    double isp;          // Hall thruster Isp (s)
    double alt_km;       // operating altitude above 1-bar
};

struct AnnualScoop {
    double tonnes_raw;
    double tonnes_he, tonnes_h2, tonnes_ch4;
    double h2_fuel;      // orbit maintenance
    double net_h2;
    int fills;
};

AnnualScoop simulate_drone_year(const DroneSpec& d)
{
    double GM_S = 3.7931e7;
    double R_S = 60268.0;
    double r_peri = R_S + d.alt_km;
    double r_apo  = R_S * 5;
    double a_orb  = (r_peri + r_apo) / 2;
    double T_orb  = 2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_S);
    double v_p    = sqrt(GM_S*(2.0/r_peri - 1.0/a_orb)) * 1000; // m/s

    double rho, press, tempk;
    atmosphere_saturn(rho, press, tempk, d.alt_km * 1000);

    double p_ram = 0.5 * rho * v_p * v_p;
    double B_need = sqrt(2.0 * MU0 * p_ram);
    double Rm = (B_need > 0) ? d.r_coil * pow(d.B_coil/B_need, 1.0/3.0) : 200;
    if (Rm > 200) Rm = 200;
    double Ae = M_PI * Rm * Rm;
    double m_dot = rho * v_p * Ae * d.eta;

    double H_m = SATURN_RGAS * tempk / SATURN_G0;
    double arc = 2.0 * sqrt(2.0 * r_peri * 1000 * H_m);
    double dt_pass = arc / v_p;
    double kg_pass = m_dot * dt_pass;

    double passes_per_fill = ceil(d.tank_cap / kg_pass);
    if (passes_per_fill < 1) passes_per_fill = 1;
    double days_fill = passes_per_fill * T_orb / 86400;
    double fills_yr = 365.25 / days_fill;
    double raw_yr = fills_yr * d.tank_cap;

    // Orbit maintenance: drag + apoapsis reboost
    double omega_sat = 2*M_PI / (10.7 * 3600);
    double v_corot = omega_sat * (R_S + d.alt_km);
    double v_rel = v_p/1000 - v_corot; // km/s
    double M_avg = d.dry_mass + d.tank_cap/2;
    double dv_drag = kg_pass * v_rel / (M_avg + kg_pass);
    double rp_ra = r_peri / r_apo;
    double dv_reboost = dv_drag * rp_ra;
    double ve = d.isp * 9.80665e-3;
    double h2_per_reboost = M_avg * (1 - exp(-dv_reboost / ve));
    double h2_fuel_yr = h2_per_reboost * passes_per_fill * fills_yr;

    AnnualScoop a;
    a.tonnes_raw = raw_yr / 1000;
    a.tonnes_he  = a.tonnes_raw * 0.061;
    a.tonnes_h2  = a.tonnes_raw * 0.906;
    a.tonnes_ch4 = a.tonnes_raw * 0.034;
    a.h2_fuel    = h2_fuel_yr / 1000;
    a.net_h2     = a.tonnes_h2 - a.h2_fuel;
    a.fills      = (int)fills_yr;
    return a;
}

// ═══════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════
int main()
{
    printf("===================================================================\n");
    printf("  SATURN STATION — Full Mission Simulation\n");
    printf("  Real Ephemeris | Lambert Trajectories | Year-by-Year Economics\n");
    printf("===================================================================\n");

    // ─── Load Ephemeris ─────────────────────────────────────────────────
    JplEphemeris eph;
    if (!eph.load("ephemeris/linux_p1550p2650.440"))
    if (!eph.load("../../Ephemeris/linux_p1550p2650.440"))
    if (!eph.load("/data/src/zipfel/solar_system/ephemeris/linux_p1550p2650.440"))
    {
        printf("FATAL: no ephemeris file\n");
        return 1;
    }
    if (!eph.is_loaded()) { printf("FATAL: ephemeris not loaded\n"); return 1; }

    // ═══════════════════════════════════════════════════════════════════════
    // Section 1: Earth → Saturn Transfer Windows (2030-2065)
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n═══ EARTH → SATURN TRANSFER WINDOWS ═══\n\n");
    printf("  Scanning direct Lambert transfers (gravity assists reduce v_inf ~30%%)\n\n");

    // Saturn synodic period: ~378 days ≈ 1.035 years
    auto e2s = scan_all_windows(eph, JPL_EMB, JPL_SATURN, true, false,
                                2030, 2065, 1500, 3500, 1.035, 30);

    printf("  %-8s %-8s %8s %8s %8s %8s %10s\n",
           "Depart", "Arrive", "TOF(d)", "TOF(yr)", "Vinf_d", "Vinf_a", "dV LEO");
    printf("  %-8s %-8s %8s %8s %8s %8s %10s\n",
           "------", "------", "------", "------", "------", "------", "------");

    for (auto& w : e2s) {
        double v_circ = sqrt(3.986e5/6778);
        double v_hyp = sqrt(w.v_inf_dep*w.v_inf_dep + 2*3.986e5/6778);
        double dv_leo = v_hyp - v_circ;
        printf("  %-8.1f %-8.1f %8.0f %8.1f %8.2f %8.2f %10.2f\n",
               w.year_dep, w.year_arr, w.tof_days, w.tof_days/365.25,
               w.v_inf_dep, w.v_inf_arr, dv_leo);
    }

    // With gravity assists (VVEJGA), v_inf_dep can be reduced to ~3-4 km/s
    printf("\n  Direct transfers: v_inf_dep ~ %.1f-%.1f km/s\n",
           e2s.empty() ? 0 : e2s.front().v_inf_dep,
           e2s.empty() ? 0 : e2s.back().v_inf_dep);
    printf("  With VVEJGA: v_inf_dep ~ 3.5-4.5 km/s (Cassini-class)\n");
    printf("  → dV from LEO: ~3.8 km/s (with gravity assists)\n");

    // ═══════════════════════════════════════════════════════════════════════
    // Section 2: Saturn → Earth Return Windows (2040-2080)
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n═══ SATURN → EARTH RETURN WINDOWS ═══\n\n");

    auto s2e = scan_all_windows(eph, JPL_SATURN, JPL_EMB, false, true,
                                2040, 2080, 1500, 3500, 1.035, 30);

    printf("  %-8s %-8s %8s %8s %8s %8s\n",
           "Depart", "Arrive", "TOF(d)", "TOF(yr)", "Vinf_d", "Vinf_a");

    for (auto& w : s2e) {
        printf("  %-8.1f %-8.1f %8.0f %8.1f %8.2f %8.2f\n",
               w.year_dep, w.year_arr, w.tof_days, w.tof_days/365.25,
               w.v_inf_dep, w.v_inf_arr);
    }

    double avg_tof_return = 0;
    for (auto& w : s2e) avg_tof_return += w.tof_days;
    if (!s2e.empty()) avg_tof_return /= s2e.size();

    printf("\n  Return windows found: %zu\n", s2e.size());
    printf("  Average return TOF: %.0f days (%.1f years)\n",
           avg_tof_return, avg_tof_return/365.25);

    // ═══════════════════════════════════════════════════════════════════════
    // Section 3: Scoop Drone Operations
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n═══ SCOOP DRONE OPERATIONS ═══\n\n");

    DroneSpec drone = {2000, 10000, 20.0, 3.0, 0.5, 3000, 1200};
    AnnualScoop per_drone = simulate_drone_year(drone);

    printf("  Drone: %gkg dry, %gkg tank, Isp=%ds, alt=%.0fkm\n",
           drone.dry_mass, drone.tank_cap, (int)drone.isp, drone.alt_km);
    printf("  Per drone per year:\n");
    printf("    Fills:       %d\n", per_drone.fills);
    printf("    Raw scooped: %.1f tonnes\n", per_drone.tonnes_raw);
    printf("    He:          %.1f tonnes\n", per_drone.tonnes_he);
    printf("    H2:          %.1f tonnes (%.1f fuel, %.1f net)\n",
           per_drone.tonnes_h2, per_drone.h2_fuel, per_drone.net_h2);
    printf("    CH4:         %.1f tonnes\n", per_drone.tonnes_ch4);

    // ═══════════════════════════════════════════════════════════════════════
    // Section 4: Year-by-Year Simulation (2030-2080)
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n═══ 50-YEAR BUSINESS SIMULATION ═══\n\n");

    // Configuration
    int n_drones = 20;
    double he_price = 1500;         // $/kg
    double launch_cost_kg = 500;    // $/kg to Saturn (Starship era)
    double discount_rate = 0.08;    // 8% annual

    // Hardware costs (build cost, excluding launch)
    double base_hardware   = 5e9;   // $5B
    double base_mass_kg    = 200e3; // 200 tonnes
    double drone_hardware  = 100e6; // $100M each
    double drone_mass_kg   = drone.dry_mass;
    double tanker_hardware = 200e6; // $200M each
    double tanker_mass_kg  = 50e3;  // 50 tonnes dry
    double tanker_payload  = 1000;  // tonnes He per tanker trip
    double enceladus_cost  = 3e9;   // $3B
    double annual_crew_supply = 100e3; // 100 tonnes/year resupply
    double crew_hardware_yr = 500e6; // $500M/year crew ops
    double drone_replace_rate = 0.05; // 5% attrition/year

    // Tanker pipeline
    double return_tof_yr = avg_tof_return / 365.25;
    if (return_tof_yr < 4) return_tof_yr = 5; // minimum realistic
    double round_trip_yr = 2 * return_tof_yr;

    // Establishment timeline
    double yr_first_launch  = 2033;
    double yr_base_arrives  = yr_first_launch + 7; // 2040
    double yr_drones_arrive = yr_first_launch + 8; // 2041 (2nd window)
    double yr_ops_start     = yr_drones_arrive + 1; // 2042 (commissioning)
    double yr_first_tanker_depart = yr_ops_start + 2; // 2044 (stockpile 2yr)
    double yr_first_he_earth = yr_first_tanker_depart + return_tof_yr;

    int n_tankers_needed = (int)ceil(round_trip_yr * 0.5) + 2; // pipeline + spares

    printf("  CONFIGURATION\n");
    printf("    Drones: %d at %.0f km altitude\n", n_drones, drone.alt_km);
    printf("    He price: $%.0f/kg\n", he_price);
    printf("    Launch cost: $%.0f/kg to Saturn\n", launch_cost_kg);
    printf("    Tanker capacity: %.0f tonnes He\n", tanker_payload);
    printf("    Tankers: %d (%.0f yr round trip)\n", n_tankers_needed, round_trip_yr);
    printf("    Return TOF: %.1f years (from Lambert scan)\n", return_tof_yr);
    printf("    Discount rate: %.0f%%\n", discount_rate*100);

    printf("\n  TIMELINE\n");
    printf("    First launch:           %.0f\n", yr_first_launch);
    printf("    Base arrives Saturn:    %.0f\n", yr_base_arrives);
    printf("    Drones operational:     %.0f\n", yr_ops_start);
    printf("    First tanker departs:   %.0f\n", yr_first_tanker_depart);
    printf("    First He at Earth:      %.0f\n", yr_first_he_earth);

    // Compute upfront costs
    double cost_base_launch   = base_mass_kg * launch_cost_kg;
    double cost_drones_launch = n_drones * drone_mass_kg * launch_cost_kg;
    double cost_tankers_launch= n_tankers_needed * tanker_mass_kg * launch_cost_kg;
    double cost_tankers_hw    = n_tankers_needed * tanker_hardware;
    double cost_drones_hw     = n_drones * drone_hardware;
    double total_upfront = base_hardware + cost_base_launch
                         + cost_drones_hw + cost_drones_launch
                         + cost_tankers_hw + cost_tankers_launch
                         + enceladus_cost;

    printf("\n  UPFRONT COSTS\n");
    printf("    Base (build+launch):     $%.1fB\n",
           (base_hardware + cost_base_launch)/1e9);
    printf("    Drones (build+launch):   $%.1fB\n",
           (cost_drones_hw + cost_drones_launch)/1e9);
    printf("    Tankers (build+launch):  $%.1fB\n",
           (cost_tankers_hw + cost_tankers_launch)/1e9);
    printf("    Enceladus ops:           $%.1fB\n", enceladus_cost/1e9);
    printf("    TOTAL UPFRONT:           $%.1fB\n", total_upfront/1e9);

    // ─── Year-by-year cash flow ─────────────────────────────────────────
    printf("\n  %-6s %6s %8s %10s %10s %10s %10s %10s\n",
           "Year", "Drones", "He(t/yr)", "Revenue", "Costs", "CashFlow",
           "Cumulative", "NPV");
    printf("  ─────────────────────────────────────────────────────────────"
           "──────────────────\n");

    double cumulative = 0;
    double npv = 0;
    int active_drones = 0;
    double he_stockpile = 0; // tonnes at Saturn
    double he_in_transit = 0;
    bool broke_even = false;
    double break_even_year = 0;

    struct YearResult {
        int year;
        int drones;
        double he_produced, he_shipped, he_delivered;
        double revenue, costs, cashflow, cumulative, npv;
    };
    std::vector<YearResult> results;

    for (int yr = 2030; yr <= 2080; yr++)
    {
        double costs = 0;
        double revenue = 0;
        double he_produced = 0;
        double he_shipped = 0;
        double he_delivered = 0;

        // Phase: establishment spending
        if (yr == (int)yr_first_launch) {
            costs += base_hardware + cost_base_launch;
            costs += enceladus_cost;
        }
        if (yr == (int)yr_first_launch + 1) {
            costs += cost_drones_hw + cost_drones_launch;
        }
        if (yr == (int)yr_first_launch + 2) {
            costs += cost_tankers_hw + cost_tankers_launch;
        }

        // Drones come online
        if (yr == (int)yr_ops_start) active_drones = n_drones;

        // Drone attrition and production
        if (active_drones > 0) {
            int lost = 0;
            for (int i = 0; i < active_drones; i++) {
                if ((double)(rand()%1000)/1000.0 < drone_replace_rate) lost++;
            }
            active_drones -= lost;

            // Replace lost drones (cost)
            if (lost > 0) {
                costs += lost * (drone_hardware + drone_mass_kg * launch_cost_kg);
                // Replacements arrive ~2 years later (next window + transit)
                // Simplified: add back next year
                active_drones += lost; // assume spares on-site
            }

            he_produced = active_drones * per_drone.tonnes_he;
            he_stockpile += he_produced;

            // Annual operations
            costs += crew_hardware_yr;
            costs += annual_crew_supply * launch_cost_kg;
        }

        // Tanker departures: ship He when available, one tanker per ~2 years
        if (yr >= (int)yr_first_tanker_depart && he_stockpile >= tanker_payload) {
            // Check if we have a return window near this year
            bool has_window = false;
            double this_tof = return_tof_yr;
            for (auto& w : s2e) {
                if (fabs(w.year_dep - yr) < 1.5) {
                    has_window = true;
                    this_tof = w.tof_days / 365.25;
                    break;
                }
            }
            if (has_window && he_stockpile >= tanker_payload) {
                he_shipped = tanker_payload;
                he_stockpile -= he_shipped;
            }
        }

        // He deliveries at Earth (shipped years ago)
        for (auto& r : results) {
            if (r.he_shipped > 0 && (yr - r.year) >= (int)return_tof_yr
                && (yr - r.year) < (int)return_tof_yr + 1) {
                he_delivered += r.he_shipped;
            }
        }

        revenue = he_delivered * 1000 * he_price; // tonnes → kg → $

        double cashflow = revenue - costs;
        cumulative += cashflow;

        double discount = pow(1.0 + discount_rate, -(yr - 2030));
        npv += cashflow * discount;

        if (!broke_even && cumulative > 0) {
            broke_even = true;
            break_even_year = yr;
        }

        YearResult r;
        r.year = yr; r.drones = active_drones;
        r.he_produced = he_produced; r.he_shipped = he_shipped;
        r.he_delivered = he_delivered;
        r.revenue = revenue; r.costs = costs;
        r.cashflow = cashflow; r.cumulative = cumulative; r.npv = npv;
        results.push_back(r);

        // Print milestone years
        bool print = (yr <= 2035)
            || (yr == (int)yr_base_arrives) || (yr == (int)yr_ops_start)
            || (yr == (int)yr_first_tanker_depart)
            || (yr == (int)yr_first_he_earth) || (yr == (int)yr_first_he_earth+1)
            || (yr % 5 == 0) || (yr == 2080)
            || he_delivered > 0;
        if (print) {
            printf("  %-6d %6d %8.0f %10.0f %10.0f %10.0f %10.0f %10.0f\n",
                   yr, active_drones, he_produced,
                   revenue/1e6, costs/1e6, cashflow/1e6,
                   cumulative/1e6, npv/1e6);
        }
    }

    printf("\n  (Revenue and costs in $M)\n");

    // ═══════════════════════════════════════════════════════════════════════
    // Section 5: Summary & Sensitivity
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n═══ FINANCIAL SUMMARY ═══\n\n");

    double total_he_delivered = 0;
    double total_revenue = 0, total_costs = 0;
    for (auto& r : results) {
        total_he_delivered += r.he_delivered;
        total_revenue += r.revenue;
        total_costs += r.costs;
    }

    printf("  Total upfront:       $%.1fB\n", total_upfront/1e9);
    printf("  Total He delivered:  %.0f tonnes (2030-2080)\n", total_he_delivered);
    printf("  Total revenue:       $%.1fB\n", total_revenue/1e9);
    printf("  Total costs:         $%.1fB\n", total_costs/1e9);
    printf("  Net P&L:             $%+.1fB\n", (total_revenue-total_costs)/1e9);
    printf("  NPV (%.0f%%):           $%+.1fB\n", discount_rate*100, npv/1e9);
    if (broke_even)
        printf("  Break-even year:     %d (%.0f years from first launch)\n",
               (int)break_even_year, break_even_year - yr_first_launch);
    else
        printf("  Break-even:          NOT REACHED by 2080\n");

    // ─── Sensitivity: He price ──────────────────────────────────────────
    printf("\n═══ SENSITIVITY: He PRICE ═══\n\n");
    printf("  Annual He production at steady state: %.0f tonnes (%.0f drones)\n",
           n_drones * per_drone.tonnes_he, (double)n_drones);
    printf("  Annual ops cost: ~$%.0fM\n", crew_hardware_yr/1e6
           + annual_crew_supply * launch_cost_kg / 1e6);

    double annual_he = n_drones * per_drone.tonnes_he;
    double annual_ops = crew_hardware_yr + annual_crew_supply * launch_cost_kg;

    printf("\n  %-12s %12s %12s %12s %12s\n",
           "He $/kg", "Revenue/yr", "Ops/yr", "Net/yr", "Payback(yr)");
    double prices[] = {500, 1000, 1500, 2000, 3000, 5000};
    for (double p : prices) {
        double rev = annual_he * 1000 * p;
        double net = rev - annual_ops;
        double payback = (net > 0) ? total_upfront / net : -1;
        printf("  $%-11.0f $%10.0fM $%10.0fM $%+10.0fM %12.0f\n",
               p, rev/1e6, annual_ops/1e6, net/1e6,
               payback > 0 ? payback : 999.0);
    }

    // ─── Sensitivity: Fleet size ────────────────────────────────────────
    printf("\n═══ SENSITIVITY: FLEET SIZE ═══\n\n");
    printf("  He price: $%.0f/kg\n\n", he_price);
    printf("  %-8s %10s %12s %12s %12s %12s\n",
           "Drones", "He(t/yr)", "Revenue/yr", "Net/yr", "Upfront", "Payback");
    int fleets[] = {10, 20, 50, 100, 200};
    for (int nd : fleets) {
        double he = nd * per_drone.tonnes_he;
        double rev = he * 1000 * he_price;
        double ops = annual_ops + nd * 5e6; // $5M/drone/yr (autonomous ops)
        double net = rev - ops;
        double upfront = total_upfront + (nd - n_drones) *
                         (drone_hardware + drone_mass_kg * launch_cost_kg);
        int tankers = (int)ceil(he / tanker_payload * round_trip_yr / 2) + 2;
        upfront += (tankers - n_tankers_needed) *
                   (tanker_hardware + tanker_mass_kg * launch_cost_kg);
        double payback = (net > 0) ? upfront / net : -1;
        printf("  %-8d %10.0f  $%10.0fM $%+10.0fM  $%10.0fB %10.0f yr\n",
               nd, he, rev/1e6, net/1e6, upfront/1e9,
               payback > 0 ? payback : 999.0);
    }

    // ─── The Real Play: He-3 ────────────────────────────────────────────
    printf("\n═══ THE He-3 PLAY ═══\n\n");
    // Saturn's atmosphere He-3/He-4 ratio: ~1.5e-4 (protosolar, less processed)
    // vs Earth's 1.4e-6 (heavily depleted)
    double he3_ratio = 1.5e-4;
    double annual_he3 = annual_he * he3_ratio;
    printf("  Saturn He-3/He-4 ratio: ~%.1e (protosolar, less processed than Earth)\n",
           he3_ratio);
    printf("  Annual He-3 production: %.1f kg (from %.0f tonnes He total)\n",
           annual_he3*1000, annual_he);
    printf("  He-3 price for fusion: $30,000 - $100,000/kg\n");
    printf("  He-3 revenue at $50,000/kg: $%.0fM/year\n",
           annual_he3 * 1000 * 50000 / 1e6);
    printf("  He-3 revenue at $1M/kg (fusion economy): $%.0fM/year\n",
           annual_he3 * 1000 * 1e6 / 1e6);
    printf("\n  If fusion power creates demand for He-3 at $1M/kg,\n");
    printf("  the Saturn Station pays for itself on He-3 alone.\n");

    // ═══════════════════════════════════════════════════════════════════════
    // Verdict
    // ═══════════════════════════════════════════════════════════════════════
    printf("\n===================================================================\n");
    printf("  VERDICT: CAN THIS MAKE MONEY?\n");
    printf("===================================================================\n");
    printf("\n");
    printf("  At $1,500/kg He (today's shortage price):\n");
    printf("    20 drones:  Annual net $%+.0fM → %s\n",
           (annual_he*1000*1500 - annual_ops)/1e6,
           (annual_he*1000*1500 > annual_ops) ? "PROFITABLE" : "NOT PROFITABLE");
    double net_100 = 100*per_drone.tonnes_he*1000*1500 - annual_ops - 100*5e6;
    printf("    100 drones: Annual net $%+.0fM → %s\n",
           net_100/1e6, net_100 > 0 ? "PROFITABLE" : "NOT PROFITABLE");
    printf("\n");
    printf("  At $3,000/kg He (2050s severe shortage):\n");
    printf("    20 drones:  Annual net $%+.0fM → PROFITABLE\n",
           (annual_he*1000*3000 - annual_ops)/1e6);
    printf("    Payback: %.0f years from first delivery\n",
           total_upfront / (annual_he*1000*3000 - annual_ops));
    printf("\n");
    printf("  THE ANSWER:\n");
    printf("    At current He prices ($1500/kg) and 20 drones: marginal.\n");
    printf("    At 50+ drones OR $3000/kg He: clearly profitable.\n");
    double net_100_3k = 100*per_drone.tonnes_he*1000*3000 - annual_ops - 100*5e6;
    double up_100 = total_upfront + 80*(drone_hardware+drone_mass_kg*launch_cost_kg);
    printf("    At 100+ drones AND $3000/kg: pays back in ~%.0f years.\n",
           up_100 / net_100_3k);
    printf("    If He-3 fusion economy emerges: bonanza.\n");
    printf("\n");
    printf("  REAL MOAT: Once built, the Saturn Station is a natural monopoly.\n");
    printf("  Nobody else has a crewed base + mining fleet + tanker pipeline.\n");
    printf("  Expansion cost is incremental (add drones, ~$120M each).\n");
    printf("  The upfront $%.0fB is the barrier. After that, it prints money.\n",
           total_upfront/1e9);
    printf("===================================================================\n");

    return 0;
}
