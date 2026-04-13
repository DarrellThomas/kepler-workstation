///////////////////////////////////////////////////////////////////////////////
// THE BUCKET CONSTELLATION — Full Fleet Simulation
//
// Simulate 10 years of operations: launching batches of Buckets every
// synodic window, tracking each through the complete loop, accounting
// for attrition, deliveries, fuel consumption, and economics.
//
// Uses real DE440 ephemeris for each launch window's trajectory.
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "jpl_ephemeris.hpp"
#include "lambert.hpp"
#include "../planetary_models/planetary_environment.hpp"

double vec_mag(const double v[3]) { return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
void vec_sub(double o[3], const double a[3], const double b[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
double jd_from_year(double y) { return 2451545.0+(y-2000.0)*365.25; }
double year_from_jd(double jd) { return 2000.0+(jd-2451545.0)/365.25; }

struct LambertResult { double v_inf_d, v_inf_a, tof_days; double jd_d, jd_a; bool ok; };

LambertResult quick_lambert(JplEphemeris& eph, JplBody b1, JplBody b2,
                            bool e1, bool e2,
                            double jd_lo, double jd_hi,
                            double tof_lo, double tof_hi, double step)
{
    LambertResult best; best.ok=false; best.v_inf_d=1e10;
    double sp[3],sv[3];
    for (double jd=jd_lo;jd<=jd_hi;jd+=step) {
        for (double tof=tof_lo;tof<=tof_hi;tof+=step) {
            double p1[3],v1[3],p2[3],v2[3],r1[3],rv1[3],r2[3],rv2[3];
            if(e1) eph.get_earth(jd,p1,v1); else eph.get_state(b1,jd,p1,v1);
            if(e2) eph.get_earth(jd+tof,p2,v2); else eph.get_state(b2,jd+tof,p2,v2);
            eph.get_state(JPL_SUN,jd,sp,sv);
            vec_sub(r1,p1,sp); vec_sub(rv1,v1,sv);
            eph.get_state(JPL_SUN,jd+tof,sp,sv);
            vec_sub(r2,p2,sp); vec_sub(rv2,v2,sv);
            double vt1[3],vt2[3];
            if(!lambert_solve(r1,r2,tof*86400,1.32712440018e11,true,vt1,vt2)) continue;
            double vd[3],va[3]; vec_sub(vd,vt1,rv1); vec_sub(va,vt2,rv2);
            double vid=vec_mag(vd), via=vec_mag(va);
            if(vid<best.v_inf_d && vid<15 && via<15) {
                best.jd_d=jd; best.jd_a=jd+tof; best.v_inf_d=vid;
                best.v_inf_a=via; best.tof_days=tof; best.ok=true;
            }
        }
    }
    return best;
}

// Bucket state machine
enum BucketPhase {
    PHASE_OUTBOUND_VENUS,    // Earth → Venus
    PHASE_OUTBOUND_EARTH,    // Venus → Earth flyby
    PHASE_OUTBOUND_JUPITER,  // Earth → Jupiter
    PHASE_CAPTURE,           // Jupiter capture
    PHASE_SCOOPING,          // Magnetic scoop operations
    PHASE_DEPART_JUPITER,    // Jupiter departure
    PHASE_RETURN,            // Jupiter → Earth
    PHASE_DELIVERED,         // At Earth depot, payload offloaded
    PHASE_LOST               // Failed (attrition)
};

struct Bucket {
    int id;
    int batch;
    BucketPhase phase;
    double jd_launch;
    double jd_phase_end;     // when current phase completes
    double h2_tank;          // kg
    double he_tank;          // kg
    double h2_fuel_used;     // kg total fuel burned this loop
    int scoop_passes;
    bool alive;
    int loop_count;          // how many complete loops
    double total_h2_delivered;
    double total_he_delivered;
};

int main()
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║        THE BUCKET CONSTELLATION — FLEET SIMULATION         ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    JplEphemeris eph;
    if (!eph.load("ephemeris/linux_p1550p2650.440"))
        eph.load("/data/src/zipfel/solar_system/ephemeris/linux_p1550p2650.440");
    if (!eph.is_loaded()) { printf("FATAL: no ephemeris\n"); return 1; }

    // Fleet parameters
    const int BUCKETS_PER_BATCH = 10;
    const double SYNODIC_WINDOW = 398.88;  // days
    const double ATTRITION_RATE = 0.10;    // 10% loss per loop
    const double DRY_MASS = 8000;          // kg
    const double TANK_CAP = 33000;         // kg
    const int ISP = 3000;
    const double VE = ISP * 9.80665e-3;    // km/s

    // Simulation time
    double sim_start_year = 2029.0;
    double sim_end_year = 2065.0;
    double sim_start_jd = jd_from_year(sim_start_year);
    double sim_end_jd = jd_from_year(sim_end_year);

    printf("Simulating %.0f - %.0f (%.0f years)\n",
           sim_start_year, sim_end_year, sim_end_year - sim_start_year);
    printf("Buckets per batch: %d, Window: %.0f days\n\n", BUCKETS_PER_BATCH, SYNODIC_WINDOW);

    // Pre-compute trajectories for each launch window
    printf("Computing trajectories for each launch window...\n");

    struct WindowTrajectory {
        double jd_launch;
        double jd_venus;     // Venus flyby
        double jd_earth;     // Earth flyby
        double jd_jupiter;   // Jupiter arrival
        double jd_depart;    // Jupiter departure (after scoop)
        double jd_home;      // Earth return
        double dv_depart_leo;// km/s from LEO
        double v_inf_jupiter;// km/s at Jupiter
        double v_inf_earth;  // km/s at Earth return
        bool valid;
        double fill_days;    // days to fill tanks
    };

    std::vector<WindowTrajectory> windows;
    int window_count = 0;

    for (double jd = sim_start_jd; jd < sim_end_jd; jd += SYNODIC_WINDOW) {
        WindowTrajectory w;
        w.jd_launch = jd;
        w.valid = false;
        double yr = year_from_jd(jd);

        // Leg 1: Earth → Venus
        LambertResult l1 = quick_lambert(eph, JPL_EMB, JPL_VENUS, true, false,
                                         jd, jd+365, 90, 250, 15);
        if (!l1.ok) { windows.push_back(w); continue; }
        w.jd_venus = l1.jd_a;

        // Leg 2: Venus → Earth
        LambertResult l2 = quick_lambert(eph, JPL_VENUS, JPL_EMB, false, true,
                                         l1.jd_a, l1.jd_a+500, 100, 400, 15);
        if (!l2.ok) { windows.push_back(w); continue; }
        w.jd_earth = l2.jd_a;

        // Leg 3: Earth → Jupiter
        LambertResult l3 = quick_lambert(eph, JPL_EMB, JPL_JUPITER, true, false,
                                         l2.jd_a, l2.jd_a+400, 500, 1200, 30);
        if (!l3.ok) { windows.push_back(w); continue; }
        w.jd_jupiter = l3.jd_a;
        w.v_inf_jupiter = l3.v_inf_a;

        // LEO departure delta-V
        double v_hyp = sqrt(l1.v_inf_d*l1.v_inf_d + 2*3.986e5/6778);
        w.dv_depart_leo = v_hyp - sqrt(3.986e5/6778);

        // Scoop campaign: simulate fill time
        double GM_J=1.26713e8, R_J=71492;
        double r_peri=R_J+350, r_apo=R_J*5;
        double a_orb=(r_peri+r_apo)/2;
        double T_orb=2*M_PI*sqrt(a_orb*a_orb*a_orb/GM_J);
        double mu0=4*M_PI*1e-7, B_coil=20, r_coil=3;

        double cumulative = 0;
        int pass = 0;
        while (cumulative < TANK_CAP && pass < 2000) {
            pass++;
            double progress = (double)pass / 1137.0;
            double alt = 400 - 100*fmin(progress, 1.0);
            double rho, press, tempk;
            atmosphere_jupiter(rho, press, tempk, alt*1000);
            double v_p = sqrt(GM_J*(2/(R_J+alt)-1/a_orb))*1000;
            double p_ram = 0.5*rho*v_p*v_p;
            double B_need = sqrt(2*mu0*p_ram);
            double R_mag = (B_need>0) ? r_coil*pow(B_coil/B_need,1.0/3.0) : 200;
            if (R_mag>200) R_mag=200;
            double A_eff = M_PI*R_mag*R_mag;
            double H = JUPITER_RGAS*tempk/JUPITER_G0;
            double arc = 2*sqrt(2*(R_J+alt)*1000*H);
            double dt = arc/v_p;
            cumulative += rho*v_p*A_eff*0.5*dt;
        }
        w.fill_days = pass * T_orb / 86400;

        // Departure from Jupiter
        double jd_filled = w.jd_jupiter + w.fill_days + 30; // +30 margin
        LambertResult l4 = quick_lambert(eph, JPL_JUPITER, JPL_EMB, false, true,
                                         jd_filled, jd_filled+800, 500, 1200, 30);
        if (!l4.ok) { windows.push_back(w); continue; }
        w.jd_depart = l4.jd_d;
        w.jd_home = l4.jd_a;
        w.v_inf_earth = l4.v_inf_a;
        w.valid = true;

        window_count++;
        windows.push_back(w);
    }

    printf("Found %d valid launch windows out of %d\n\n", window_count, (int)windows.size());

    // Create bucket fleet
    std::vector<Bucket> fleet;
    int next_id = 1;

    // Simulation loop: step through time day by day
    double cumulative_he_delivered = 0;
    double cumulative_h2_delivered = 0;
    int total_launched = 0;
    int total_lost = 0;
    int total_deliveries = 0;

    // Annual tracking
    struct AnnualStats {
        double year;
        int launched, active, lost, delivered_count;
        double he_tonnes, h2_tonnes;
    };
    std::vector<AnnualStats> annual;

    double current_year = sim_start_year;
    int win_idx = 0;

    printf("═══ CONSTELLATION OPERATIONS LOG ═══\n\n");
    printf("%-8s %-6s %-6s %-6s %-8s %-10s %-10s %-10s\n",
           "Year", "Launch", "Active", "Lost", "Deliver", "He(t)", "H2(t)", "CumHe(t)");
    printf("─────────────────────────────────────────────────────────────────────\n");

    for (double yr = sim_start_year; yr < sim_end_year; yr += 1.0) {
        AnnualStats stats;
        stats.year = yr;
        stats.launched = 0;
        stats.lost = 0;
        stats.delivered_count = 0;
        stats.he_tonnes = 0;
        stats.h2_tonnes = 0;

        // Launch batches that fall in this year
        for (int wi = 0; wi < (int)windows.size(); wi++) {
            double launch_yr = year_from_jd(windows[wi].jd_launch);
            if (launch_yr >= yr && launch_yr < yr+1 && windows[wi].valid) {
                WindowTrajectory& w = windows[wi];
                for (int b = 0; b < BUCKETS_PER_BATCH; b++) {
                    Bucket bk;
                    bk.id = next_id++;
                    bk.batch = wi;
                    bk.phase = PHASE_OUTBOUND_VENUS;
                    bk.jd_launch = w.jd_launch;
                    bk.jd_phase_end = w.jd_home; // simplify: track completion time
                    bk.h2_tank = 0;
                    bk.he_tank = 0;
                    bk.h2_fuel_used = 0;
                    bk.scoop_passes = 0;
                    bk.alive = true;
                    bk.loop_count = 0;
                    bk.total_h2_delivered = 0;
                    bk.total_he_delivered = 0;
                    fleet.push_back(bk);
                    stats.launched++;
                    total_launched++;
                }
            }
        }

        // Check for deliveries and attrition this year
        for (auto& bk : fleet) {
            if (!bk.alive) continue;
            double arrive_yr = year_from_jd(bk.jd_phase_end);
            if (arrive_yr >= yr && arrive_yr < yr+1) {
                // Bucket completes its loop this year
                // Attrition check
                double roll = (double)(rand() % 1000) / 1000.0;
                if (roll < ATTRITION_RATE) {
                    bk.alive = false;
                    bk.phase = PHASE_LOST;
                    stats.lost++;
                    total_lost++;
                    continue;
                }

                // Successful delivery!
                bk.phase = PHASE_DELIVERED;
                bk.loop_count++;

                // Calculate payload using Ganymede-assist ΔV budget
                double dv_return = 3.5;  // with Ganymede assist
                double dv_outbound = 4.0;
                double dv_total = dv_return + dv_outbound + 0.5;
                double mr = exp(dv_total / VE);

                double scooped_h2 = TANK_CAP * 0.75;
                double scooped_he = TANK_CAP * 0.25;
                double h2_payload = (scooped_h2-(DRY_MASS+scooped_he)*(mr-1))/mr;
                double he_payload = scooped_he;

                if (h2_payload < 0) h2_payload = 0;

                bk.total_h2_delivered += h2_payload;
                bk.total_he_delivered += he_payload;
                cumulative_h2_delivered += h2_payload;
                cumulative_he_delivered += he_payload;
                stats.h2_tonnes += h2_payload / 1000;
                stats.he_tonnes += he_payload / 1000;
                stats.delivered_count++;
                total_deliveries++;

                // Relaunch on next available window
                for (int wi = 0; wi < (int)windows.size(); wi++) {
                    if (windows[wi].valid && year_from_jd(windows[wi].jd_launch) > yr) {
                        bk.jd_launch = windows[wi].jd_launch;
                        bk.jd_phase_end = windows[wi].jd_home;
                        bk.phase = PHASE_OUTBOUND_VENUS;
                        break;
                    }
                }
            }
        }

        int active = 0;
        for (auto& bk : fleet) if (bk.alive) active++;
        stats.active = active;

        if (stats.launched > 0 || stats.delivered_count > 0 || ((int)yr % 5 == 0)) {
            printf("%-8.0f %-6d %-6d %-6d %-8d %-10.1f %-10.1f %-10.1f\n",
                   yr, stats.launched, active, stats.lost, stats.delivered_count,
                   stats.he_tonnes, stats.h2_tonnes, cumulative_he_delivered/1000);
        }
        annual.push_back(stats);
    }

    // Summary
    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  CONSTELLATION SUMMARY: %.0f - %.0f                          ║\n",
           sim_start_year, sim_end_year);
    printf("╠══════════════════════════════════════════════════════════════╣\n");
    printf("║  Total launched:     %d buckets                           ║\n", total_launched);
    printf("║  Total lost:         %d (%.0f%% attrition)                  ║\n",
           total_lost, (double)total_lost/total_launched*100);
    printf("║  Total deliveries:   %d                                   ║\n", total_deliveries);
    printf("║  He delivered:       %.0f tonnes                           ║\n",
           cumulative_he_delivered/1000);
    printf("║  H2 delivered:       %.0f tonnes                           ║\n",
           cumulative_h2_delivered/1000);
    printf("║  Total delivered:    %.0f tonnes                           ║\n",
           (cumulative_he_delivered+cumulative_h2_delivered)/1000);
    printf("╠══════════════════════════════════════════════════════════════╣\n");

    // Economics
    double he_price_now = 220;
    double he_price_future = 1100; // 5× by 2060
    double build_cost = total_launched * 245e6;
    double revenue_now = cumulative_he_delivered * he_price_now;
    double revenue_future = cumulative_he_delivered * he_price_future;

    printf("║  ECONOMICS (at today\\'s He price \$220/kg)                  ║\n");
    printf("║    He revenue:       \$%.1fB                                ║\n",
           revenue_now/1e9);
    printf("║    Fleet cost:       \$%.1fB                                ║\n",
           build_cost/1e9);
    printf("║    ROI:              %.0f%%                                  ║\n",
           (revenue_now-build_cost)/build_cost*100);
    printf("║                                                            ║\n");
    printf("║  ECONOMICS (at 5× shortage \$1100/kg by 2060)              ║\n");
    printf("║    He revenue:       \$%.1fB                                ║\n",
           revenue_future/1e9);
    printf("║    Fleet cost:       \$%.1fB                                ║\n",
           build_cost/1e9);
    printf("║    ROI:              %.0f%%                                  ║\n",
           (revenue_future-build_cost)/build_cost*100);
    printf("╠══════════════════════════════════════════════════════════════╣\n");

    // Helium impact
    double he_earth_reserve = 8.6e6; // tonnes
    double he_consumption = 30000;    // tonnes/year
    double years_sim = sim_end_year - sim_start_year;
    double he_consumed_period = he_consumption * years_sim;
    double he_fraction_replaced = cumulative_he_delivered/1000 / he_consumed_period * 100;

    printf("║  HELIUM IMPACT                                             ║\n");
    printf("║    Earth consumed (%.0fyr): %.0f tonnes                      ║\n",
           years_sim, he_consumed_period);
    printf("║    Jupiter delivered:     %.0f tonnes                       ║\n",
           cumulative_he_delivered/1000);
    printf("║    Replaced:              %.1f%% of consumption             ║\n",
           he_fraction_replaced);
    printf("║    At 10× fleet:          %.0f%% of consumption             ║\n",
           he_fraction_replaced*10);
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    // Per-bucket lifetime stats
    printf("\n═══ BUCKET LIFETIME EXAMPLES ═══\n\n");
    int shown = 0;
    for (auto& bk : fleet) {
        if (bk.loop_count > 0 && shown < 5) {
            printf("  Bucket #%d: %d loops, delivered %.1ft He + %.1ft H2, %s\n",
                   bk.id, bk.loop_count,
                   bk.total_he_delivered/1000, bk.total_h2_delivered/1000,
                   bk.alive ? "still active" : "LOST");
            shown++;
        }
    }

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  %d buckets launched. %d still flying.                     ║\n",
           total_launched, total_launched-total_lost);
    printf("║  %.0f tonnes of helium delivered to Earth.                  ║\n",
           cumulative_he_delivered/1000);
    printf("║  Every kg extends civilization's helium supply.            ║\n");
    printf("║                                                            ║\n");
    printf("║  The Bucket constellation doesn't stop.                    ║\n");
    printf("║  It just keeps flying.                                     ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}
