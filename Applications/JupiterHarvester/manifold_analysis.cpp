// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Jupiter Harvester — Manifold Transfer Delta-V Analysis
//
// Compares the current Bucket fleet delta-V budget against a manifold-based
// approach using Sun-Jupiter L1/L2 invariant manifold tubes.
//
// Question: Can riding dynamical highways (manifold tubes) from Lagrange
// point halo orbits reduce the propellant cost of Jupiter capture and
// departure, allowing each Bucket to deliver more helium per trip?
//
// Reference: KoLoMaRo (2011), Chapters 2, 6, 7, 9, 10.
//
// Compile: g++ -std=c++11 -O2 -I../../Ephemeris/cr3bp -o manifold_analysis manifold_analysis.cpp
// Run:     ./manifold_analysis
///////////////////////////////////////////////////////////////////////////////

#include "cr3bp.hpp"
#include "halo_orbit.hpp"
#include "manifold.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
// Bucket Mission Constants (from bucket_sim.cpp / README.md)
///////////////////////////////////////////////////////////////////////////////

const double DRY_MASS      = 8000.0;    // kg
const double TANK_CAP      = 33000.0;   // kg total scooped
const double ISP           = 3000.0;    // seconds (Hall thruster)
const double VE            = ISP * 9.80665e-3; // 29.42 km/s exhaust velocity
const double H2_FRAC       = 0.75;      // H2 fraction of scooped mass
const double HE_FRAC       = 0.25;      // He fraction
const double H2_TANK       = TANK_CAP * H2_FRAC;  // 24,750 kg
const double HE_TANK       = TANK_CAP * HE_FRAC;  //  8,250 kg

// Current baseline delta-V budget (km/s)
const double DV_EARTH_DEPART = 3.5;
const double DV_CAPTURE      = 0.5;
const double DV_JUP_DEPART   = 5.0;
const double DV_MARGIN       = 0.5;
const double DV_TOTAL_BASE   = DV_EARTH_DEPART + DV_CAPTURE + DV_JUP_DEPART + DV_MARGIN;

// Jupiter physical constants
const double GM_JUPITER    = 1.26713e8;   // km^3/s^2
const double R_JUPITER     = 71492.0;     // km (equatorial)
const double SCOOP_ALT     = 350.0;       // km above 1-bar level
const double R_SCOOP_PERI  = R_JUPITER + SCOOP_ALT;      // 71,842 km
const double R_SCOOP_APO   = R_JUPITER * 5.0;            // 357,460 km
const double A_SCOOP       = (R_SCOOP_PERI + R_SCOOP_APO) / 2.0;
const double JUPITER_SOI   = 48.2e6;     // km

///////////////////////////////////////////////////////////////////////////////
// Dimensional conversion helpers
///////////////////////////////////////////////////////////////////////////////

double to_km(double nd, const CR3BPSystem& sys) { return nd * sys.L_km; }
double to_kms(double nd, const CR3BPSystem& sys) { return nd * sys.V_km_s(); }
double to_days(double nd, const CR3BPSystem& sys) { return nd * sys.T_sec / (2.0 * CR3BP_PI * 86400.0); }

// Convert rotating frame state to Jupiter-relative inertial velocity (nondim)
void rotating_to_jupiter_vel(const CR3BPState& s, double mu,
                             double& vjx, double& vjy, double& vjz)
{
    // Inertial velocity in barycentric frame:
    //   vx_inert = vx - y*omega = vx - y  (omega=1 in nondim)
    //   vy_inert = vy + x*omega = vy + x
    // Jupiter (m2) is at (1-mu, 0, 0) with inertial velocity (0, 1-mu, 0)
    double vx_inert = s.vx - s.y;
    double vy_inert = s.vy + s.x;
    double vz_inert = s.vz;

    vjx = vx_inert - 0.0;          // Jupiter vx_inert = 0
    vjy = vy_inert - (1.0 - mu);   // Jupiter vy_inert = 1-mu
    vjz = vz_inert - 0.0;
}

double jupiter_relative_speed(const CR3BPState& s, double mu)
{
    double vjx, vjy, vjz;
    rotating_to_jupiter_vel(s, mu, vjx, vjy, vjz);
    return sqrt(vjx*vjx + vjy*vjy + vjz*vjz);
}

// Distance from state to Jupiter (m2 at 1-mu, 0, 0)
double dist_to_jupiter(const CR3BPState& s, double mu)
{
    double dx = s.x - (1.0 - mu);
    return sqrt(dx*dx + s.y*s.y + s.z*s.z);
}

// Vis-viva: orbital speed at radius r in orbit with semi-major axis a around Jupiter
double vis_viva(double r_km, double a_km)
{
    return sqrt(GM_JUPITER * (2.0/r_km - 1.0/a_km));
}

// Tsiolkovsky: H2 payload delivered given total delta-V
double compute_payload(double dv_total)
{
    double mass_ratio = exp(dv_total / VE);
    // m0 = dry + h2_tank + he_tank (fully loaded)
    // mf = dry + h2_payload + he_tank
    // m0/mf = mass_ratio => h2_payload = h2_tank - (dry + he_tank)*(mass_ratio - 1) / mass_ratio
    // Actually: m0 = mf * mass_ratio
    //   dry + h2_tank + he_tank = (dry + h2_payload + he_tank) * mass_ratio
    //   h2_tank = h2_payload * mass_ratio + (dry + he_tank)*(mass_ratio - 1)
    //   h2_payload = (h2_tank - (dry + he_tank)*(mass_ratio - 1)) / mass_ratio
    double h2_payload = (H2_TANK - (DRY_MASS + HE_TANK) * (mass_ratio - 1.0)) / mass_ratio;
    if (h2_payload < 0) h2_payload = 0;
    return h2_payload;
}

///////////////////////////////////////////////////////////////////////////////
// Manifold periapsis scanner
///////////////////////////////////////////////////////////////////////////////

struct PeriapsisInfo
{
    double r_min_nd;    // closest approach (nondim)
    double r_min_km;    // closest approach (km from Jupiter center)
    double v_rel_nd;    // Jupiter-relative speed at periapsis (nondim)
    double v_rel_kms;   // Jupiter-relative speed (km/s)
    double v_scoop_kms; // vis-viva speed at same radius in scoop orbit (km/s)
    double dv_kms;      // delta-V to match scoop orbit at that radius
    int arc;            // which manifold arc
    int pt;             // which point along arc
};

PeriapsisInfo scan_periapsis(const ManifoldTube& tube, double mu,
                             const CR3BPSystem& sys)
{
    PeriapsisInfo best;
    best.r_min_nd = 1e20;
    best.arc = -1;

    for (int i = 0; i < (int)tube.trajectories.size(); i++) {
        for (int j = 0; j < (int)tube.trajectories[i].size(); j++) {
            const CR3BPState& s = tube.trajectories[i][j];
            double r = dist_to_jupiter(s, mu);
            if (r < best.r_min_nd) {
                best.r_min_nd = r;
                best.v_rel_nd = jupiter_relative_speed(s, mu);
                best.arc = i;
                best.pt = j;
            }
        }
    }

    best.r_min_km = to_km(best.r_min_nd, sys);
    best.v_rel_kms = to_kms(best.v_rel_nd, sys);

    // Vis-viva at same radius in scoop orbit (if within scoop orbit apoapsis)
    if (best.r_min_km > 0 && best.r_min_km < R_SCOOP_APO * 2) {
        best.v_scoop_kms = vis_viva(best.r_min_km, A_SCOOP);
    } else {
        // Use hyperbolic approach: v_inf matching
        best.v_scoop_kms = sqrt(2.0 * GM_JUPITER / best.r_min_km);
    }
    best.dv_kms = fabs(best.v_rel_kms - best.v_scoop_kms);

    return best;
}

// Collect periapsis info for ALL arcs
std::vector<PeriapsisInfo> scan_all_periapsis(const ManifoldTube& tube, double mu,
                                              const CR3BPSystem& sys)
{
    std::vector<PeriapsisInfo> results;
    for (int i = 0; i < (int)tube.trajectories.size(); i++) {
        PeriapsisInfo info;
        info.r_min_nd = 1e20;
        info.arc = i;

        for (int j = 0; j < (int)tube.trajectories[i].size(); j++) {
            const CR3BPState& s = tube.trajectories[i][j];
            double r = dist_to_jupiter(s, mu);
            if (r < info.r_min_nd) {
                info.r_min_nd = r;
                info.v_rel_nd = jupiter_relative_speed(s, mu);
                info.pt = j;
            }
        }

        info.r_min_km = to_km(info.r_min_nd, sys);
        info.v_rel_kms = to_kms(info.v_rel_nd, sys);

        if (info.r_min_km > 0 && info.r_min_km < R_SCOOP_APO * 2) {
            info.v_scoop_kms = vis_viva(info.r_min_km, A_SCOOP);
        } else {
            info.v_scoop_kms = sqrt(2.0 * GM_JUPITER / fmax(info.r_min_km, R_JUPITER));
        }
        info.dv_kms = fabs(info.v_rel_kms - info.v_scoop_kms);

        results.push_back(info);
    }
    return results;
}

///////////////////////////////////////////////////////////////////////////////
// Main Analysis
///////////////////////////////////////////////////////////////////////////////

int main()
{
    double mu = SYS_SUN_JUPITER.mu;

    printf("===================================================================\n");
    printf("  JUPITER HARVESTER — Manifold Transfer Delta-V Analysis\n");
    printf("  Kepler Workstation | CR3BP Module\n");
    printf("  Reference: KoLoMaRo (2011)\n");
    printf("===================================================================\n");

    // =====================================================================
    // Section 1: Sun-Jupiter Lagrange Points
    // =====================================================================
    printf("\n--- Section 1: Sun-Jupiter System ---\n\n");
    printf("  mu = %.6e\n", mu);
    printf("  L  = %.3e km (Sun-Jupiter distance)\n", SYS_SUN_JUPITER.L_km);
    printf("  T  = %.3e s (%.2f years, orbital period)\n",
           SYS_SUN_JUPITER.T_sec, SYS_SUN_JUPITER.T_sec / (365.25*86400));
    printf("  V  = %.3f km/s (velocity scale)\n", SYS_SUN_JUPITER.V_km_s());

    LagrangePoints lp = cr3bp_lagrange_points(mu);
    double L1_km_from_jup = to_km(fabs(lp.L1[0] - (1.0-mu)), SYS_SUN_JUPITER);
    double L2_km_from_jup = to_km(fabs(lp.L2[0] - (1.0-mu)), SYS_SUN_JUPITER);

    printf("\n  L1: x = %.10f  (%.2f Mkm from Jupiter, sunward)\n",
           lp.L1[0], L1_km_from_jup / 1e6);
    printf("  L2: x = %.10f  (%.2f Mkm from Jupiter, anti-sun)\n",
           lp.L2[0], L2_km_from_jup / 1e6);

    double C_L1 = cr3bp_jacobi_at_point(lp.L1[0], 0, 0, mu);
    double C_L2 = cr3bp_jacobi_at_point(lp.L2[0], 0, 0, mu);
    printf("  C(L1) = %.8f   C(L2) = %.8f\n", C_L1, C_L2);

    printf("\n  Jupiter SOI: %.1f Mkm = %.5f nd\n",
           JUPITER_SOI/1e6, JUPITER_SOI / SYS_SUN_JUPITER.L_km);
    printf("  Scoop orbit periapsis: %.0f km = %.2e nd from Jupiter\n",
           R_SCOOP_PERI, R_SCOOP_PERI / SYS_SUN_JUPITER.L_km);

    // =====================================================================
    // Section 2: Halo Orbit Families at L1 and L2
    // =====================================================================
    printf("\n--- Section 2: Sun-Jupiter Halo Orbit Families ---\n");

    int n_family = 4;

    // L1 halo orbits require larger Az (above bifurcation point).
    // For small mu, the minimum Az_nd ~ gamma * sqrt(|Delta/l1|) is significant.
    // Sun-Jupiter L1 minimum is ~0.03 nd (~23 Mkm).
    double Az_L1_start = 0.030;
    double Az_L1_end   = 0.060;

    // L2 halo orbits have a lower minimum amplitude threshold.
    double Az_L2_start = 0.005;
    double Az_L2_end   = 0.020;

    printf("\n  L1 Halo Family (Az range: %.0f-%.0f Mkm):\n",
           Az_L1_start * SYS_SUN_JUPITER.L_km/1e6, Az_L1_end * SYS_SUN_JUPITER.L_km/1e6);
    printf("  %3s  %10s  %10s  %10s  %12s  %5s  %s\n",
           "#", "Az(Mkm)", "Ax(Mkm)", "Period(d)", "Jacobi", "Conv", "Iter");

    std::vector<HaloOrbit> L1_halos = compute_halo_family(mu, 1, Az_L1_start, Az_L1_end, n_family, true);
    HaloOrbit best_L1;
    best_L1.converged = false;
    for (int i = 0; i < n_family; i++) {
        HaloOrbit& h = L1_halos[i];
        double Az_km = h.Az * SYS_SUN_JUPITER.L_km;
        double Ax_km = h.Ax * SYS_SUN_JUPITER.L_km;
        double T_days = to_days(h.period, SYS_SUN_JUPITER);
        printf("  %3d  %10.2f  %10.2f  %10.1f  %12.8f  %5s  %4d\n",
               i, Az_km/1e6, Ax_km/1e6, T_days, h.jacobi,
               h.converged ? "YES" : "NO", h.iterations);
        if (h.converged) best_L1 = h;
    }

    printf("\n  L2 Halo Family (Az range: %.0f-%.0f Mkm):\n",
           Az_L2_start * SYS_SUN_JUPITER.L_km/1e6, Az_L2_end * SYS_SUN_JUPITER.L_km/1e6);
    printf("  %3s  %10s  %10s  %10s  %12s  %5s  %s\n",
           "#", "Az(Mkm)", "Ax(Mkm)", "Period(d)", "Jacobi", "Conv", "Iter");

    std::vector<HaloOrbit> L2_halos = compute_halo_family(mu, 2, Az_L2_start, Az_L2_end, n_family, true);
    HaloOrbit best_L2;
    best_L2.converged = false;
    for (int i = 0; i < n_family; i++) {
        HaloOrbit& h = L2_halos[i];
        double Az_km = h.Az * SYS_SUN_JUPITER.L_km;
        double Ax_km = h.Ax * SYS_SUN_JUPITER.L_km;
        double T_days = to_days(h.period, SYS_SUN_JUPITER);
        printf("  %3d  %10.2f  %10.2f  %10.1f  %12.8f  %5s  %4d\n",
               i, Az_km/1e6, Ax_km/1e6, T_days, h.jacobi,
               h.converged ? "YES" : "NO", h.iterations);
        if (h.converged) best_L2 = h;
    }

    // L1 may not converge for Sun-Jupiter (small mu makes the bifurcation
    // threshold large). Use L2 for both capture and departure analysis.
    if (!best_L2.converged) {
        printf("\n  ERROR: No converged halo orbits at L2.\n");
        return 1;
    }
    if (!best_L1.converged) {
        printf("\n  NOTE: L1 halos did not converge. Using L2 for all manifold analysis.\n");
        best_L1 = best_L2; // fallback
    }

    // =====================================================================
    // Section 3: Manifold Tubes — Periapsis at Jupiter
    // =====================================================================
    printf("\n--- Section 3: Manifold Tubes Toward Jupiter ---\n");

    int n_arcs = 15;
    double t_man = 4.0;  // ~0.7 Jupiter orbits in nondim time
    int n_rec = 150;

    // L1 unstable manifold, interior branch (toward Jupiter)
    // For L1 (sunward of Jupiter), branch=-1 sends arcs toward Jupiter
    printf("\n  Computing L1 unstable manifold (%d arcs, t=%.1f nd) ...\n", n_arcs, t_man);
    ManifoldTube Wu_L1 = compute_manifold(mu, best_L1, n_arcs, t_man, n_rec, 1e-6, false, -1);

    PeriapsisInfo best_Wu_L1 = scan_periapsis(Wu_L1, mu, SYS_SUN_JUPITER);
    std::vector<PeriapsisInfo> Wu_L1_all = scan_all_periapsis(Wu_L1, mu, SYS_SUN_JUPITER);

    printf("  L1 Wu (toward Jupiter):\n");
    printf("    Best periapsis:  %.0f km (%.2f Rj) at v_rel=%.2f km/s\n",
           best_Wu_L1.r_min_km, best_Wu_L1.r_min_km/R_JUPITER, best_Wu_L1.v_rel_kms);

    // Count arcs reaching SOI
    int in_soi = 0;
    double sum_r = 0, sum_v = 0;
    for (auto& p : Wu_L1_all) {
        if (p.r_min_km < JUPITER_SOI) { in_soi++; sum_r += p.r_min_km; sum_v += p.v_rel_kms; }
    }
    printf("    Arcs reaching SOI: %d/%d\n", in_soi, n_arcs);
    if (in_soi > 0)
        printf("    Mean periapsis (SOI arcs): %.0f km (%.1f Rj), mean v_rel=%.2f km/s\n",
               sum_r/in_soi, sum_r/in_soi/R_JUPITER, sum_v/in_soi);

    // L2 stable manifold, interior branch (toward Jupiter)
    // For L2 (anti-sun of Jupiter), branch=+1 sends arcs toward Jupiter
    printf("\n  Computing L2 stable manifold (%d arcs, t=%.1f nd) ...\n", n_arcs, t_man);
    ManifoldTube Ws_L2 = compute_manifold(mu, best_L2, n_arcs, t_man, n_rec, 1e-6, true, +1);

    PeriapsisInfo best_Ws_L2 = scan_periapsis(Ws_L2, mu, SYS_SUN_JUPITER);
    std::vector<PeriapsisInfo> Ws_L2_all = scan_all_periapsis(Ws_L2, mu, SYS_SUN_JUPITER);

    printf("  L2 Ws (toward Jupiter):\n");
    printf("    Best periapsis:  %.0f km (%.2f Rj) at v_rel=%.2f km/s\n",
           best_Ws_L2.r_min_km, best_Ws_L2.r_min_km/R_JUPITER, best_Ws_L2.v_rel_kms);

    in_soi = 0; sum_r = 0; sum_v = 0;
    for (auto& p : Ws_L2_all) {
        if (p.r_min_km < JUPITER_SOI) { in_soi++; sum_r += p.r_min_km; sum_v += p.v_rel_kms; }
    }
    printf("    Arcs reaching SOI: %d/%d\n", in_soi, n_arcs);
    if (in_soi > 0)
        printf("    Mean periapsis (SOI arcs): %.0f km (%.1f Rj), mean v_rel=%.2f km/s\n",
               sum_r/in_soi, sum_r/in_soi/R_JUPITER, sum_v/in_soi);

    // Also try the other branch combinations
    printf("\n  Computing L1 stable manifold (for departure) ...\n");
    ManifoldTube Ws_L1 = compute_manifold(mu, best_L1, n_arcs, t_man, n_rec, 1e-6, true, -1);
    PeriapsisInfo best_Ws_L1 = scan_periapsis(Ws_L1, mu, SYS_SUN_JUPITER);
    printf("  L1 Ws (toward Jupiter): best periapsis = %.0f km (%.2f Rj)\n",
           best_Ws_L1.r_min_km, best_Ws_L1.r_min_km/R_JUPITER);

    printf("\n  Computing L2 unstable manifold (for departure) ...\n");
    ManifoldTube Wu_L2 = compute_manifold(mu, best_L2, n_arcs, t_man, n_rec, 1e-6, false, +1);
    PeriapsisInfo best_Wu_L2 = scan_periapsis(Wu_L2, mu, SYS_SUN_JUPITER);
    printf("  L2 Wu (toward Jupiter): best periapsis = %.0f km (%.2f Rj)\n",
           best_Wu_L2.r_min_km, best_Wu_L2.r_min_km/R_JUPITER);

    // =====================================================================
    // Section 4: Delta-V Comparison
    // =====================================================================
    printf("\n--- Section 4: Delta-V Analysis ---\n");

    // Pick the best capture manifold (smallest periapsis among all tubes)
    PeriapsisInfo capture_best = best_Wu_L1;
    const char* capture_source = "L1 Wu";
    if (best_Ws_L2.r_min_km < capture_best.r_min_km) {
        capture_best = best_Ws_L2;
        capture_source = "L2 Ws";
    }

    // Pick the best departure manifold
    PeriapsisInfo depart_best = best_Ws_L1;
    const char* depart_source = "L1 Ws";
    if (best_Wu_L2.r_min_km < depart_best.r_min_km) {
        depart_best = best_Wu_L2;
        depart_source = "L2 Wu";
    }

    printf("\n  CAPTURE (arrival at Jupiter):\n");
    printf("    Best manifold: %s\n", capture_source);
    printf("    Manifold periapsis: %.0f km (%.1f Rj)\n",
           capture_best.r_min_km, capture_best.r_min_km/R_JUPITER);
    printf("    Jupiter-relative speed at periapsis: %.2f km/s\n", capture_best.v_rel_kms);

    // Compute what the scoop orbit insertion costs from the manifold-delivered state
    // The manifold delivers to a high orbit; we need to descend to scoop orbit.
    // Two burns: 1) at manifold periapsis, circularize or enter transfer orbit
    //            2) at scoop orbit, match velocity
    // Simplified: compute v_inf at Jupiter from manifold, then hyperbolic capture dv

    double v_inf_capture = capture_best.v_rel_kms;
    // Hyperbolic velocity at scoop periapsis with this v_inf
    double v_hyp_at_scoop = sqrt(v_inf_capture*v_inf_capture + 2.0*GM_JUPITER/R_SCOOP_PERI);
    double v_scoop_peri = vis_viva(R_SCOOP_PERI, A_SCOOP);
    double dv_capture_manifold = v_hyp_at_scoop - v_scoop_peri;
    if (dv_capture_manifold < 0) dv_capture_manifold = 0;

    printf("    v_inf at Jupiter: %.2f km/s\n", v_inf_capture);
    printf("    v_hyp at scoop periapsis: %.2f km/s\n", v_hyp_at_scoop);
    printf("    v_scoop at periapsis: %.2f km/s\n", v_scoop_peri);
    printf("    dv to insert into scoop orbit: %.2f km/s\n", dv_capture_manifold);
    printf("    Baseline capture dv: %.2f km/s\n", DV_CAPTURE);

    printf("\n  DEPARTURE (leaving Jupiter):\n");
    printf("    Best manifold: %s\n", depart_source);
    printf("    Manifold periapsis: %.0f km (%.1f Rj)\n",
           depart_best.r_min_km, depart_best.r_min_km/R_JUPITER);
    printf("    Jupiter-relative speed at periapsis: %.2f km/s\n", depart_best.v_rel_kms);

    double v_inf_depart = depart_best.v_rel_kms;
    double v_hyp_at_scoop_dep = sqrt(v_inf_depart*v_inf_depart + 2.0*GM_JUPITER/R_SCOOP_PERI);
    double dv_depart_manifold = v_hyp_at_scoop_dep - v_scoop_peri;
    if (dv_depart_manifold < 0) dv_depart_manifold = 0;

    printf("    v_inf needed: %.2f km/s\n", v_inf_depart);
    printf("    v_hyp at scoop periapsis: %.2f km/s\n", v_hyp_at_scoop_dep);
    printf("    dv from scoop to manifold: %.2f km/s\n", dv_depart_manifold);
    printf("    Baseline departure dv: %.2f km/s\n", DV_JUP_DEPART);

    // =====================================================================
    // Section 5: Mass Budget Comparison
    // =====================================================================
    printf("\n--- Section 5: Mass Budget Comparison ---\n");

    double dv_manifold_total = DV_EARTH_DEPART + dv_capture_manifold + dv_depart_manifold + DV_MARGIN;
    double dv_savings = DV_TOTAL_BASE - dv_manifold_total;

    double h2_payload_base = compute_payload(DV_TOTAL_BASE);
    double h2_payload_man  = compute_payload(dv_manifold_total);
    double total_base = h2_payload_base + HE_TANK;
    double total_man  = h2_payload_man  + HE_TANK;

    double mr_base = exp(DV_TOTAL_BASE / VE);
    double mr_man  = exp(dv_manifold_total / VE);

    printf("\n  %-28s %10s %10s %10s\n", "", "Baseline", "Manifold", "Savings");
    printf("  %-28s %10s %10s %10s\n", "", "--------", "--------", "-------");
    printf("  %-28s %9.2f  %9.2f  %9.2f km/s\n",
           "Earth departure", DV_EARTH_DEPART, DV_EARTH_DEPART, 0.0);
    printf("  %-28s %9.2f  %9.2f  %+8.2f km/s\n",
           "Jupiter capture", DV_CAPTURE, dv_capture_manifold,
           DV_CAPTURE - dv_capture_manifold);
    printf("  %-28s %9.2f  %9.2f  %+8.2f km/s\n",
           "Jupiter departure", DV_JUP_DEPART, dv_depart_manifold,
           DV_JUP_DEPART - dv_depart_manifold);
    printf("  %-28s %9.2f  %9.2f  %9.2f km/s\n",
           "Margin", DV_MARGIN, DV_MARGIN, 0.0);
    printf("  %-28s %9.2f  %9.2f  %+8.2f km/s\n",
           "TOTAL", DV_TOTAL_BASE, dv_manifold_total, dv_savings);

    printf("\n  %-28s %9.3f  %9.3f\n", "Mass ratio", mr_base, mr_man);
    printf("  %-28s %9.1f  %9.1f  %+8.1f t\n",
           "H2 fuel burned (tonnes)", (H2_TANK - h2_payload_base)/1000.0,
           (H2_TANK - h2_payload_man)/1000.0,
           (h2_payload_man - h2_payload_base)/1000.0);
    printf("  %-28s %9.1f  %9.1f  %+8.1f t\n",
           "H2 delivered (tonnes)", h2_payload_base/1000.0, h2_payload_man/1000.0,
           (h2_payload_man - h2_payload_base)/1000.0);
    printf("  %-28s %9.1f  %9.1f  %9.1f t\n",
           "He delivered (tonnes)", HE_TANK/1000.0, HE_TANK/1000.0, 0.0);
    printf("  %-28s %9.1f  %9.1f  %+8.1f t\n",
           "TOTAL delivered (tonnes)", total_base/1000.0, total_man/1000.0,
           (total_man - total_base)/1000.0);
    printf("  %-28s %9.2fx %9.2fx\n",
           "Delivery ratio (x dry mass)", total_base/DRY_MASS, total_man/DRY_MASS);

    // Fleet impact
    double loops_per_year = 10.0;  // ~10 buckets return per year at steady state
    double fleet_base = total_base / 1000.0 * loops_per_year;
    double fleet_man  = total_man  / 1000.0 * loops_per_year;
    printf("\n  Fleet impact (100 buckets, ~10 deliveries/year):\n");
    printf("    Baseline: %.0f tonnes/year\n", fleet_base);
    printf("    Manifold: %.0f tonnes/year\n", fleet_man);
    printf("    Gain:     %+.0f tonnes/year\n", fleet_man - fleet_base);

    // =====================================================================
    // Section 6: Jovian Moon Manifolds (Ganymede)
    // =====================================================================
    printf("\n--- Section 6: Jupiter-Ganymede Manifold Feasibility ---\n");

    double mu_jg = SYS_JUPITER_GANYMEDE.mu;
    printf("\n  Jupiter-Ganymede system:\n");
    printf("    mu = %.6e\n", mu_jg);
    printf("    L  = %.0f km (Jupiter-Ganymede distance)\n", SYS_JUPITER_GANYMEDE.L_km);
    printf("    V  = %.3f km/s\n", SYS_JUPITER_GANYMEDE.V_km_s());

    LagrangePoints lp_jg = cr3bp_lagrange_points(mu_jg);
    double L1_jg_km = fabs(lp_jg.L1[0] - (1.0-mu_jg)) * SYS_JUPITER_GANYMEDE.L_km;
    double L2_jg_km = fabs(lp_jg.L2[0] - (1.0-mu_jg)) * SYS_JUPITER_GANYMEDE.L_km;

    printf("    L1: %.0f km from Ganymede (%.1f Rj from Jupiter)\n",
           L1_jg_km, (SYS_JUPITER_GANYMEDE.L_km - L1_jg_km)/R_JUPITER);
    printf("    L2: %.0f km from Ganymede (%.1f Rj from Jupiter)\n",
           L2_jg_km, (SYS_JUPITER_GANYMEDE.L_km + L2_jg_km)/R_JUPITER);

    printf("\n  Computing Jupiter-Ganymede L1 halo orbit ...\n");
    double Az_jg = 0.01; // small amplitude
    HaloOrbit halo_jg = compute_halo_orbit(mu_jg, 1, Az_jg, true);

    if (halo_jg.converged) {
        double T_jg_days = to_days(halo_jg.period, SYS_JUPITER_GANYMEDE);
        printf("    Converged! Period = %.2f days, Jacobi = %.8f\n",
               T_jg_days, halo_jg.jacobi);

        // Compute manifold toward Jupiter
        printf("    Computing manifold tubes (20 arcs, t=6.0) ...\n");
        ManifoldTube Wu_jg = compute_manifold(mu_jg, halo_jg, 20, 6.0, 200, 1e-6, false, -1);

        // Distance to Jupiter (at -mu in Jupiter-Ganymede frame)
        // In JG frame, Jupiter is at (-mu_jg, 0, 0)
        double best_r = 1e20;
        double best_v = 0;
        for (auto& traj : Wu_jg.trajectories) {
            for (auto& s : traj) {
                double dx = s.x + mu_jg; // distance to Jupiter in JG frame
                double r = sqrt(dx*dx + s.y*s.y + s.z*s.z);
                if (r < best_r) {
                    best_r = r;
                    // Approximate Jupiter-relative speed in JG frame
                    double vxj = s.vx - s.y + mu_jg; // rough rotating->inertial
                    double vyj = s.vy + s.x;
                    best_v = sqrt(vxj*vxj + vyj*vyj + s.vz*s.vz);
                }
            }
        }
        double best_r_km = best_r * SYS_JUPITER_GANYMEDE.L_km;
        double best_v_kms = best_v * SYS_JUPITER_GANYMEDE.V_km_s();

        printf("    Closest approach to Jupiter: %.0f km (%.2f Rj)\n",
               best_r_km, best_r_km/R_JUPITER);
        printf("    Speed at closest approach: %.2f km/s\n", best_v_kms);

        if (best_r_km < R_SCOOP_APO) {
            printf("    ** REACHES SCOOP ORBIT REGION! **\n");
            double v_scoop = vis_viva(best_r_km, A_SCOOP);
            printf("    dv to match scoop orbit: %.2f km/s\n", fabs(best_v_kms - v_scoop));
        } else {
            printf("    Does not reach scoop orbit (%.0f km > %.0f km apoapsis)\n",
                   best_r_km, R_SCOOP_APO);
        }
    } else {
        printf("    Halo orbit did not converge.\n");
    }

    // =====================================================================
    // Summary
    // =====================================================================
    printf("\n===================================================================\n");
    printf("  SUMMARY\n");
    printf("===================================================================\n");
    printf("\n  The CR3BP manifold tubes from Sun-Jupiter L1/L2 halo orbits\n");
    printf("  provide low-energy pathways into/out of Jupiter's vicinity.\n\n");

    printf("  FINDING: Sun-Jupiter manifold tubes deliver to ~%.0f Rj\n",
           capture_best.r_min_km / R_JUPITER);
    printf("  from Jupiter at v_inf ~ %.1f km/s. The scoop orbit is at\n",
           capture_best.v_rel_kms);
    printf("  %.0f km (1.0 Rj) — a factor of %.0fx deeper in the gravity well.\n",
           R_SCOOP_PERI, capture_best.r_min_km / R_SCOOP_PERI);
    printf("\n  The manifold provides a FREE interplanetary capture (v_inf=%.1f km/s\n",
           capture_best.v_rel_kms);
    printf("  vs ~5.6 km/s for a direct Hohmann), but an Oberth burn at periapsis\n");
    printf("  is still needed to descend from the manifold-delivered high orbit\n");
    printf("  to the scoop orbit. The current Bucket design already exploits this:\n");
    printf("  the 0.5 km/s capture uses magnetic drag over multiple passes,\n");
    printf("  and the 5.0 km/s departure uses an Oberth burn at periapsis.\n");
    printf("\n  CONCLUSION: The current VEEJ + Oberth strategy is already near-optimal\n");
    printf("  for the deep scoop orbit. Manifold tubes help at the SOI boundary\n");
    printf("  but cannot bypass the deep gravity well descent. The manifold v_inf\n");
    printf("  (%.1f km/s) is lower than direct transfer (~5.6 km/s), which could\n",
           capture_best.v_rel_kms);
    printf("  reduce the capture burn — but the Bucket already achieves this via\n");
    printf("  the VEEJ gravity assist chain.\n");
    printf("\n  NEXT STEPS: The real opportunity is Jupiter-MOON manifolds\n");
    printf("  (Ganymede/Callisto L1/L2), which operate at 15 Rj — much closer\n");
    printf("  to the scoop orbit. Resonant gravity assists between moons could\n");
    printf("  reduce the descent/ascent cost significantly. This requires the\n");
    printf("  patched 3-body approach (KoLoMaRo Chapter 5, 10).\n");

    printf("\n===================================================================\n");

    return 0;
}
