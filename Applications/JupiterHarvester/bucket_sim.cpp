///////////////////////////////////////////////////////////////////////////////
// THE BUCKET — Full Loop Simulation
//
// End-to-end simulation of one Bucket doing one complete loop:
//   1. Earth departure (VEEJ trajectory via Lambert)
//   2. Venus flyby (gravity assist, no harvest — keep it simple)
//   3. Earth flyby (gravity assist)
//   4. Jupiter approach and capture
//   5. Magnetic scoop campaign (pass-by-pass)
//   6. Jupiter departure (Oberth + Hall thruster)
//   7. Earth return and aerocapture
//
// Uses real ephemeris (DE440) for trajectory design.
// Magnetic scoop physics modeled per-pass.
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <vector>
#include "jpl_ephemeris.hpp"
#include "lambert.hpp"
#include "../planetary_models/planetary_environment.hpp"

// Utilities
double vec_mag(const double v[3]) {
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
void vec_sub(double o[3], const double a[3], const double b[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
double jd_from_year(double y) { return 2451545.0 + (y-2000.0)*365.25; }
double year_from_jd(double jd) { return 2000.0 + (jd-2451545.0)/365.25; }

// Quick Lambert scan — returns best v_inf_depart
struct Xfer { double jd_dep, jd_arr, v_inf_d, v_inf_a, tof_days; bool ok; };

Xfer scan_transfer(JplEphemeris& eph, JplBody b1, JplBody b2,
                   bool e1, bool e2,
                   double jd_lo, double jd_hi,
                   double tof_lo, double tof_hi, double step)
{
    Xfer best; best.ok = false; best.v_inf_d = 1e10;
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
            if (!lambert_solve(r1,r2,tof*86400,1.32712440018e11,true,vt1,vt2)) continue;
            double vd[3],va[3]; vec_sub(vd,vt1,rv1); vec_sub(va,vt2,rv2);
            double vid = vec_mag(vd), via = vec_mag(va);
            if (vid < best.v_inf_d && vid < 15 && via < 15) {
                best.jd_dep=jd; best.jd_arr=jd+tof; best.v_inf_d=vid;
                best.v_inf_a=via; best.tof_days=tof; best.ok=true;
            }
        }
    }
    return best;
}

int main()
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║              THE BUCKET — FULL LOOP SIMULATION             ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    JplEphemeris eph;
    if (!eph.load("ephemeris/linux_p1550p2650.440"))
        eph.load("/data/src/zipfel/solar_system/ephemeris/linux_p1550p2650.440");
    if (!eph.is_loaded()) { printf("FATAL: no ephemeris\n"); return 1; }

    // ═══ VEHICLE SPEC ═══
    double dry_mass = 8000;     // kg
    double tank_cap = 33000;    // kg max H2+He storage
    double isp = 3000;          // Hall thruster
    double ve = isp * 9.80665e-3; // km/s

    // Magnetic scoop parameters
    double B_coil = 20.0;      // Tesla
    double r_coil = 3.0;       // meters
    double mu0 = 4*M_PI*1e-7;
    double eta_collect = 0.5;

    printf("═══ VEHICLE SPEC ═══\n");
    printf("  Dry mass: %g kg | Tank: %g kg | Isp: %d s\n", dry_mass, tank_cap, (int)isp);
    printf("  Magnet: %gT, %gm radius | Collection eff: %.0f%%\n\n", B_coil, r_coil, eta_collect*100);

    // ═══ LEG 1: Earth → Venus ═══
    printf("═══ LEG 1: Earth → Venus ═══\n");
    double jd0 = jd_from_year(2029.0);
    Xfer leg1 = scan_transfer(eph, JPL_EMB, JPL_VENUS, true, false,
                              jd0, jd_from_year(2031), 90, 250, 10);
    if (!leg1.ok) { printf("No solution\n"); return 1; }

    double v_circ_leo = sqrt(3.986e5/6778.0);
    double v_hyp = sqrt(leg1.v_inf_d*leg1.v_inf_d + 2*3.986e5/6778.0);
    double dv_depart = v_hyp - v_circ_leo;

    printf("  Depart %.2f → Arrive %.2f (%.0f days)\n",
           year_from_jd(leg1.jd_dep), year_from_jd(leg1.jd_arr), leg1.tof_days);
    printf("  V∞=%.2f km/s, ΔV from LEO=%.2f km/s\n\n", leg1.v_inf_d, dv_depart);

    // ═══ LEG 2: Venus → Earth ═══
    printf("═══ LEG 2: Venus → Earth (gravity assist) ═══\n");
    Xfer leg2 = scan_transfer(eph, JPL_VENUS, JPL_EMB, false, true,
                              leg1.jd_arr, leg1.jd_arr+500, 100, 400, 10);
    if (!leg2.ok) { printf("No solution\n"); return 1; }
    printf("  Depart %.2f → Arrive %.2f (%.0f days)\n",
           year_from_jd(leg2.jd_dep), year_from_jd(leg2.jd_arr), leg2.tof_days);
    printf("  V∞ depart=%.2f, V∞ arrive=%.2f km/s\n\n", leg2.v_inf_d, leg2.v_inf_a);

    // ═══ LEG 3: Earth → Jupiter ═══
    printf("═══ LEG 3: Earth → Jupiter ═══\n");
    Xfer leg3 = scan_transfer(eph, JPL_EMB, JPL_JUPITER, true, false,
                              leg2.jd_arr, leg2.jd_arr+400, 500, 1200, 20);
    if (!leg3.ok) { printf("No solution\n"); return 1; }

    double v_entry_jup = sqrt(leg3.v_inf_a*leg3.v_inf_a + 2*1.26713e8/(71492+350));
    printf("  Depart %.2f → Arrive %.2f (%.0f days, %.1f yr)\n",
           year_from_jd(leg3.jd_dep), year_from_jd(leg3.jd_arr),
           leg3.tof_days, leg3.tof_days/365.25);
    printf("  V∞=%.2f km/s, entry at 350km=%.1f km/s\n\n", leg3.v_inf_a, v_entry_jup);

    // ═══ JUPITER CAPTURE ═══
    printf("═══ JUPITER CAPTURE ═══\n");
    double dv_capture = 0.5; // small burn + magnetic drag over multiple passes
    printf("  Capture burn: %.1f km/s (+ magnetic atmospheric drag)\n", dv_capture);
    printf("  Settling into 350km × 5Rj ellipse over ~10 orbits\n\n");

    // ═══ SCOOP CAMPAIGN ═══
    printf("═══ SCOOP CAMPAIGN (pass-by-pass) ═══\n\n");

    double GM_J = 1.26713e8;
    double R_J = 71492.0;
    double alt_peri = 350.0;
    double r_peri = R_J + alt_peri;
    double r_apo = R_J * 5;
    double a_orb = (r_peri + r_apo)/2;
    double T_orb = 2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_J);
    double v_peri = sqrt(GM_J*(2.0/r_peri - 1.0/a_orb))*1000; // m/s

    double campaign_sec = 2 * 365.25 * 86400;
    int max_passes = (int)(campaign_sec / T_orb);

    printf("  Orbit: %.0f km × %.0f Rj, period=%.1f hrs\n",
           alt_peri, r_apo/R_J, T_orb/3600);
    printf("  V_peri: %.1f km/s, max passes in 2yr: %d\n\n", v_peri/1000, max_passes);

    printf("  %-6s %-8s %-10s %-10s %-10s %-10s %-8s\n",
           "Pass", "Alt(km)", "rho", "R_mag(m)", "kg/pass", "Total(t)", "Tank%");

    double h2_tank = 0, he_tank = 0, total_tank = 0;
    int actual_passes = 0;

    // Simulate varying altitude: start high, gradually dip lower as confidence grows
    for (int pass = 1; pass <= max_passes; pass++)
    {
        // Altitude schedule: start at 400km, drop to 300km over campaign
        double progress = (double)pass / max_passes;
        double alt = 400.0 - 100.0 * progress; // 400→300 km
        double r_p = R_J + alt;

        // Recalc velocity for this periapsis (vis-viva)
        double v_p = sqrt(GM_J*(2.0/r_p - 1.0/a_orb))*1000;

        // Atmosphere
        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt * 1000.0);

        // Magnetic funnel
        double p_ram = 0.5 * rho * v_p * v_p;
        double B_need = sqrt(2*mu0*p_ram);
        double R_mag = (B_need > 0) ? r_coil * pow(B_coil/B_need, 1.0/3.0) : 1000;
        if (R_mag > 200) R_mag = 200; // practical limit
        double A_eff = M_PI * R_mag * R_mag;

        // Mass flow
        double m_dot = rho * v_p * A_eff * eta_collect;

        // Pass duration
        double H = JUPITER_RGAS * tempk / JUPITER_G0;
        double arc = 2.0 * sqrt(2.0 * r_p * 1000 * H);
        double dt_pass = arc / v_p;

        double harvested = m_dot * dt_pass;
        double h2 = harvested * 0.75;
        double he = harvested * 0.25;

        h2_tank += h2;
        he_tank += he;
        total_tank += harvested;
        actual_passes = pass;

        // Print milestones
        if (pass == 1 || pass == 10 || pass == 50 || pass == 100 ||
            pass == 200 || pass == 500 || pass == max_passes ||
            total_tank >= tank_cap)
        {
            printf("  %-6d %-8.0f %-10.2e %-10.0f %-10.2f %-10.1f %-8.1f\n",
                   pass, alt, rho, R_mag, harvested, total_tank/1000,
                   total_tank/tank_cap*100);
        }

        if (total_tank >= tank_cap) {
            printf("  *** TANK FULL at pass %d (day %.0f of campaign) ***\n",
                   pass, pass*T_orb/86400);
            break;
        }
    }

    printf("\n  Campaign summary:\n");
    printf("    Passes used:     %d of %d available\n", actual_passes, max_passes);
    printf("    H2 in tank:      %.1f tonnes\n", h2_tank/1000);
    printf("    He in tank:      %.1f tonnes\n", he_tank/1000);
    printf("    Total:           %.1f tonnes (%.0f%% of capacity)\n",
           total_tank/1000, total_tank/tank_cap*100);

    // ═══ JUPITER DEPARTURE ═══
    printf("\n═══ JUPITER DEPARTURE ═══\n\n");

    // Oberth burn at periapsis + Hall thruster spiral
    double dv_jup_depart = 5.0;
    printf("  Oberth burn at periapsis + spiral: %.1f km/s\n", dv_jup_depart);

    // ═══ RETURN LEG ═══
    printf("\n═══ RETURN: Jupiter → Earth ═══\n");
    double jd_jup_depart = leg3.jd_arr + campaign_sec/86400 + 30; // +30 days margin
    Xfer leg_ret = scan_transfer(eph, JPL_JUPITER, JPL_EMB, false, true,
                                 jd_jup_depart, jd_jup_depart+800, 500, 1200, 20);
    if (!leg_ret.ok) { printf("  No return solution found\n"); return 1; }

    printf("  Depart %.2f → Arrive %.2f (%.0f days, %.1f yr)\n",
           year_from_jd(leg_ret.jd_dep), year_from_jd(leg_ret.jd_arr),
           leg_ret.tof_days, leg_ret.tof_days/365.25);
    printf("  V∞ depart=%.2f, V∞ arrive Earth=%.2f km/s\n", leg_ret.v_inf_d, leg_ret.v_inf_a);

    double v_earth_entry = sqrt(leg_ret.v_inf_a*leg_ret.v_inf_a + 2*3.986e5/6578);
    printf("  Earth aerocapture at %.1f km/s: FREE\n\n", v_earth_entry);

    // ═══ MASS BUDGET ═══
    printf("═══ MASS BUDGET ═══\n\n");

    double dv_total = dv_depart + dv_capture + dv_jup_depart + 0.5; // +margin
    double mass_ratio = exp(dv_total / ve);

    double h2_payload = (h2_tank - (dry_mass + he_tank)*(mass_ratio - 1)) / mass_ratio;
    double h2_fuel = h2_tank - h2_payload;
    double delivered = h2_payload + he_tank;

    printf("  ΔV total loop:     %.2f km/s\n", dv_total);
    printf("  Mass ratio:        %.4f\n", mass_ratio);
    printf("  H2 scooped:        %.1f tonnes\n", h2_tank/1000);
    printf("  H2 as fuel:        %.1f tonnes\n", h2_fuel/1000);
    printf("  H2 delivered:      %.1f tonnes\n", h2_payload/1000);
    printf("  He delivered:      %.1f tonnes\n", he_tank/1000);
    printf("  ─────────────────────────────\n");
    printf("  TOTAL DELIVERED:   %.1f tonnes\n", delivered/1000);
    printf("  Delivery ratio:    %.1f× dry mass\n", delivered/dry_mass);

    // ═══ TIMELINE ═══
    double jd_start = leg1.jd_dep;
    double jd_end = leg_ret.jd_arr;
    double total_years = (jd_end - jd_start) / 365.25;

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  MISSION TIMELINE                                          ║\n");
    printf("╠══════════════════════════════════════════════════════════════╣\n");
    printf("║  %.2f  Earth departure (%.1f km/s from LEO)            ║\n",
           year_from_jd(leg1.jd_dep), dv_depart);
    printf("║  %.2f  Venus flyby (gravity assist)                    ║\n",
           year_from_jd(leg1.jd_arr));
    printf("║  %.2f  Earth flyby (gravity assist)                    ║\n",
           year_from_jd(leg2.jd_arr));
    printf("║  %.2f  Jupiter arrival (capture %.1f km/s)              ║\n",
           year_from_jd(leg3.jd_arr), dv_capture);
    printf("║  +2 years   Magnetic scoop: %d passes, %.0ft captured       ║\n",
           actual_passes, total_tank/1000);
    printf("║  %.2f  Jupiter departure (%.1f km/s Oberth+Hall)       ║\n",
           year_from_jd(leg_ret.jd_dep), dv_jup_depart);
    printf("║  %.2f  Earth arrival (aerocapture at %.0f km/s)         ║\n",
           year_from_jd(leg_ret.jd_arr), v_earth_entry);
    printf("║                                                            ║\n");
    printf("║  LOOP TIME:    %.1f years                                  ║\n", total_years);
    printf("║  DELIVERED:    %.1f tonnes (%.1ft H2 + %.1ft He)          ║\n",
           delivered/1000, h2_payload/1000, he_tank/1000);
    printf("║  VEHICLE:      %.0ft dry → reuse on next loop              ║\n", dry_mass/1000.0);
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    // ═══ FLEET PROJECTION ═══
    printf("\n═══ FLEET: 10 BUCKETS PER STARSHIP ═══\n\n");

    int buckets = 10;
    double window = 1.09; // years
    double attrition = 0.10;
    int pipeline = (int)ceil(total_years / window);
    int fleet = buckets * pipeline;
    double eff = buckets * (1-attrition) / window;
    double annual_tonnes = eff * delivered / 1000;
    double cost_per = 245e6;

    printf("  Loop time:         %.1f years\n", total_years);
    printf("  Pipeline depth:    %d windows\n", pipeline);
    printf("  Total fleet:       %d buckets\n", fleet);
    printf("  Annual delivery:   %.0f tonnes\n", annual_tonnes);
    printf("  Fleet cost:        $%.1fB\n", fleet*cost_per/1e9);
    printf("  Annual value:      $%.0fM (at $1500/kg)\n", annual_tonnes*1e3*1500/1e6);

    printf("\n");
    printf("  First 5 loops:\n");
    for (int i = 0; i < 5; i++) {
        double launch = year_from_jd(jd_start) + i*window;
        double arrive_j = launch + (year_from_jd(leg3.jd_arr)-year_from_jd(jd_start));
        double arrive_e = launch + total_years;
        printf("    Batch %d: launch %.1f → Jupiter %.1f → Earth %.1f → %.0ft delivered\n",
               i+1, launch, arrive_j, arrive_e, delivered/1000);
    }

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  8 tonnes. 3 valves. 1 magnet. 1 thruster.                ║\n");
    printf("║  It flies to Jupiter. It fills up. It comes home.          ║\n");
    printf("║  It does this until it breaks.                             ║\n");
    printf("║  Then we send another one.                                 ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}
