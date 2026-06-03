///////////////////////////////////////////////////////////////////////////////
// PORKCHOP PLOT GENERATOR — Earth→Mars Launch Windows
//
// Sweeps departure dates and times of flight to compute the full
// "porkchop plot" of interplanetary transfer performance.
//
// Uses JPL DE440 ephemeris + Lambert solver for real planetary positions.
// Outputs CSV data suitable for plotting with Python/matplotlib or gnuplot.
//
// Compile:
//   g++ -std=c++11 -O2 -o porkchop_plot porkchop_plot.cpp
//   (requires ../../Ephemeris/linux_p1550p2650.440)
//
// Output:
//   porkchop_earth_mars.csv  — full sweep data
//   stdout                   — top 5 windows + ASCII summary
//
// 20260603 Created for Kepler Workstation OrbitalTransfer application
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
// CONSTANTS
///////////////////////////////////////////////////////////////////////////////

static const double GM_SUN      = 1.32712440018e11; // km³/s²
static const double AU_KM       = 149597870.7;      // km/AU
static const double SEC_PER_DAY = 86400.0;

static const double R_EARTH_KM  = 6371.0;
static const double R_MARS_KM   = 3390.0;
static const double PARK_ALT_KM = 300.0;
static const double MU_EARTH    = 3.986004418e5;
static const double MU_MARS     = 4.282837e4;

// Search parameters
static const double DEPART_START_YEAR = 2026.0;
static const double DEPART_END_YEAR   = 2036.0;
static const double DEPART_STEP_DAYS  = 20.0;   // departure grid spacing
static const double TOF_MIN_DAYS      = 50.0;
static const double TOF_MAX_DAYS      = 450.0;
static const double TOF_STEP_DAYS     = 10.0;   // TOF grid spacing

///////////////////////////////////////////////////////////////////////////////
// VECTOR HELPERS
///////////////////////////////////////////////////////////////////////////////

static inline double vmag(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline void vsub(double o[3], const double a[3], const double b[3]) {
    o[0] = a[0] - b[0]; o[1] = a[1] - b[1]; o[2] = a[2] - b[2];
}

static double jd_from_year(double y) {
    return 2451545.0 + (y - 2000.0) * 365.25;
}

static double year_from_jd(double jd) {
    return 2000.0 + (jd - 2451545.0) / 365.25;
}

// Stumpff functions (inline from lambert.hpp — needed for standalone compile)
static inline double stumpff_c(double z) {
    if (fabs(z) < 1e-6) return 1.0/2.0 - z/24.0 + z*z/720.0;
    if (z > 0) return (1.0 - cos(sqrt(z))) / z;
    return (cosh(sqrt(-z)) - 1.0) / (-z);
}

static inline double stumpff_s(double z) {
    if (fabs(z) < 1e-6) return 1.0/6.0 - z/120.0 + z*z/5040.0;
    if (z > 0) {
        double sz = sqrt(z);
        return (sz - sin(sz)) / (sz * sz * sz);
    }
    double sz = sqrt(-z);
    return (sinh(sz) - sz) / (sz * sz * sz);
}


///////////////////////////////////////////////////////////////////////////////
// LAMBERT SOLVER (bundled for standalone compile — mirrors lambert.hpp)
///////////////////////////////////////////////////////////////////////////////

static bool lambert_solve_bundled(const double r1[3], const double r2[3],
                                   double tof, double mu, bool prograde,
                                   double v1[3], double v2[3])
{
    const int MAX_ITER = 100;
    const double TOL = 1e-10;

    double R1 = vmag(r1);
    double R2 = vmag(r2);
    double cross_z = r1[0]*r2[1] - r1[1]*r2[0];
    double cos_dnu = (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]) / (R1 * R2);
    if (cos_dnu > 1.0) cos_dnu = 1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;
    double dnu = acos(cos_dnu);

    if (prograde) { if (cross_z < 0) dnu = 2.0 * M_PI - dnu; }
    else          { if (cross_z >= 0) dnu = 2.0 * M_PI - dnu; }

    double sin_dnu = sin(dnu);
    if (fabs(sin_dnu) < 1e-15) return false;
    double A = sin_dnu * sqrt(R1 * R2 / (1.0 - cos_dnu));

    double z = 0.0;
    double sqrt_mu = sqrt(mu);

    for (int iter = 0; iter < MAX_ITER; iter++) {
        double C = stumpff_c(z);
        double S = stumpff_s(z);
        double y = R1 + R2 + A * (z * S - 1.0) / sqrt(C);
        if (y < 0) { z += 0.1; continue; }

        double x = sqrt(y / C);
        double tof_guess = (x*x*x * S + A * sqrt(y)) / sqrt_mu;

        // Numerical derivative
        double dz = 0.01;
        double C2 = stumpff_c(z + dz);
        double S2 = stumpff_s(z + dz);
        double y2 = R1 + R2 + A * ((z+dz) * S2 - 1.0) / sqrt(C2);
        if (y2 < 0) y2 = fabs(y2);
        double x2 = sqrt(y2 / C2);
        double tof2 = (x2*x2*x2 * S2 + A * sqrt(y2)) / sqrt_mu;
        double dtdz = (tof2 - tof_guess) / dz;

        if (fabs(dtdz) < 1e-20) break;
        double z_new = z - (tof_guess - tof) / dtdz;
        if (fabs(z_new - z) < TOL) { z = z_new; break; }
        z = z_new;
        if (z > 4.0*M_PI*M_PI) z = 4.0*M_PI*M_PI - 0.1;
        if (z < -1000.0) z = -1000.0;
    }

    double C = stumpff_c(z);
    double S = stumpff_s(z);
    double y = R1 + R2 + A * (z * S - 1.0) / sqrt(C);
    if (y < 0) return false;

    double f = 1.0 - y / R1;
    double g = A * sqrt(y / mu);
    double gdot = 1.0 - y / R2;
    if (fabs(g) < 1e-20) return false;

    for (int i = 0; i < 3; i++) {
        v1[i] = (r2[i] - f * r1[i]) / g;
        v2[i] = (gdot * r2[i] - r1[i]) / g;
    }
    return true;
}


///////////////////////////////////////////////////////////////////////////////
// TRANSFER RESULT
///////////////////////////////////////////////////////////////////////////////

struct PorkchopPoint {
    double jd_depart;
    double jd_arrive;
    double tof_days;
    double c3_depart;      // km²/s²
    double v_inf_depart;   // km/s
    double v_inf_arrive;   // km/s
    double dv_total;       // km/s (approx, from parking orbits)
    bool   valid;
};

PorkchopPoint compute_one(const double r1[3], const double rv1[3],
                           const double r2[3], const double rv2[3],
                           double jd_d, double tof_days)
{
    PorkchopPoint p;
    p.jd_depart = jd_d;
    p.jd_arrive = jd_d + tof_days;
    p.tof_days  = tof_days;
    p.valid     = false;

    double vt1[3], vt2[3];
    if (!lambert_solve_bundled(r1, r2, tof_days * SEC_PER_DAY, GM_SUN, true, vt1, vt2))
        return p;

    double vd[3], va[3];
    vsub(vd, vt1, rv1);
    vsub(va, vt2, rv2);
    double vid = vmag(vd);
    double via = vmag(va);

    // Filter unreasonable transfers (C3 < 2500 km²/s²)
    if (vid > 50.0 || via > 50.0) return p;

    p.v_inf_depart = vid;
    p.v_inf_arrive = via;
    p.c3_depart    = vid * vid;

    // Approximate total ΔV from parking orbits
    double r_park_e = R_EARTH_KM + PARK_ALT_KM;
    double r_park_m = R_MARS_KM  + PARK_ALT_KM;
    double v_circ_e = sqrt(MU_EARTH / r_park_e);
    double v_circ_m = sqrt(MU_MARS  / r_park_m);
    double v_hyp_d  = sqrt(vid*vid + 2.0*MU_EARTH/r_park_e);
    double v_hyp_a  = sqrt(via*via + 2.0*MU_MARS /r_park_m);
    p.dv_total      = (v_hyp_d - v_circ_e) + (v_hyp_a - v_circ_m);
    p.valid         = true;

    return p;
}


///////////////////////////////////////////////////////////////////////////////
// CIRCULAR-COPLANAR FALLBACK (no ephemeris needed)
///////////////////////////////////////////////////////////////////////////////

void compute_circular_sweep(std::vector<PorkchopPoint>& points)
{
    // Circular, coplanar Earth and Mars orbits
    double r_earth = 1.00000011 * AU_KM;
    double r_mars  = 1.523679   * AU_KM;
    double v_earth = sqrt(GM_SUN / r_earth);  // km/s
    double v_mars  = sqrt(GM_SUN / r_mars);

    // Earth and Mars velocity vectors (circular, in xy-plane)
    // Position angle changes with time; for the circular model we
    // parameterize by the phase angle between Earth and Mars.

    double jd_start = jd_from_year(DEPART_START_YEAR);
    double jd_end   = jd_from_year(DEPART_END_YEAR);
    double omega_earth = (v_earth / r_earth) * SEC_PER_DAY;  // rad/day
    double omega_mars  = (v_mars  / r_mars)  * SEC_PER_DAY;

    // Mars starts ~180° ahead of Earth (approximate opposition),
    // producing periodic favorable transfer windows as Earth laps Mars.
    double mars_phase_offset = M_PI;

    for (double jd_d = jd_start; jd_d <= jd_end; jd_d += DEPART_STEP_DAYS) {
        for (double tof = TOF_MIN_DAYS; tof <= TOF_MAX_DAYS; tof += TOF_STEP_DAYS) {
            double jd_a = jd_d + tof;

            // Earth position at departure (angle from reference)
            double theta_e_d = omega_earth * (jd_d - jd_start);
            double theta_m_a = mars_phase_offset + omega_mars * (jd_a - jd_start);

            double r1[3] = { r_earth * cos(theta_e_d), r_earth * sin(theta_e_d), 0 };
            double r2[3] = { r_mars  * cos(theta_m_a), r_mars  * sin(theta_m_a),  0 };

            // Planet velocities (circular, perpendicular to position)
            double rv1[3] = { -v_earth * sin(theta_e_d), v_earth * cos(theta_e_d), 0 };
            double rv2[3] = { -v_mars  * sin(theta_m_a), v_mars  * cos(theta_m_a),  0 };

            PorkchopPoint p = compute_one(r1, rv1, r2, rv2, jd_d, tof);
            if (p.valid) points.push_back(p);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// EPHEMERIS-BASED SWEEP
///////////////////////////////////////////////////////////////////////////////

#ifdef HAS_EPHEMERIS
#include "../../Ephemeris/jpl_ephemeris.hpp"

void compute_ephemeris_sweep(JplEphemeris& eph, std::vector<PorkchopPoint>& points)
{
    double jd_start = jd_from_year(DEPART_START_YEAR);
    double jd_end   = jd_from_year(DEPART_END_YEAR);
    double sp[3], sv[3];

    printf("  Computing %.0f departure dates × %.0f TOF values...\n",
           (jd_end - jd_start) / DEPART_STEP_DAYS,
           (TOF_MAX_DAYS - TOF_MIN_DAYS) / TOF_STEP_DAYS);

    int count = 0;
    for (double jd_d = jd_start; jd_d <= jd_end; jd_d += DEPART_STEP_DAYS) {
        double ep[3], ev[3];
        if (!eph.get_earth(jd_d, ep, ev)) continue;
        eph.get_state(JPL_SUN, jd_d, sp, sv);
        double r1[3], rv1[3];
        vsub(r1, ep, sp);
        vsub(rv1, ev, sv);

        for (double tof = TOF_MIN_DAYS; tof <= TOF_MAX_DAYS; tof += TOF_STEP_DAYS) {
            double jd_a = jd_d + tof;
            double mp[3], mv[3];
            if (!eph.get_state(JPL_MARS, jd_a, mp, mv)) continue;
            eph.get_state(JPL_SUN, jd_a, sp, sv);
            double r2[3], rv2[3];
            vsub(r2, mp, sp);
            vsub(rv2, mv, sv);

            PorkchopPoint p = compute_one(r1, rv1, r2, rv2, jd_d, tof);
            if (p.valid) points.push_back(p);
            count++;
        }
    }
    printf("  Evaluated %d transfer windows, %zu converged.\n",
           count, points.size());
}
#endif


///////////////////////////////////////////////////////////////////////////////
// OUTPUT
///////////////////////////////////////////////////////////////////////////////

void write_csv(const std::vector<PorkchopPoint>& points, const char* path)
{
    FILE* f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "Cannot write %s\n", path);
        return;
    }
    fprintf(f, "departure_jd,arrival_jd,departure_year,arrival_year,tof_days,"
               "c3_depart_km2s2,v_inf_depart_km_s,v_inf_arrive_km_s,"
               "dv_total_km_s\n");
    for (const auto& p : points) {
        fprintf(f, "%.3f,%.3f,%.4f,%.4f,%.1f,%.3f,%.3f,%.3f,%.3f\n",
                p.jd_depart, p.jd_arrive,
                year_from_jd(p.jd_depart), year_from_jd(p.jd_arrive),
                p.tof_days, p.c3_depart, p.v_inf_depart, p.v_inf_arrive,
                p.dv_total);
    }
    fclose(f);
    printf("  Wrote %zu points to %s\n", points.size(), path);
}

void print_top_windows(const std::vector<PorkchopPoint>& points, int n)
{
    // Sort by C3 (ascending)
    std::vector<PorkchopPoint> sorted = points;
    std::sort(sorted.begin(), sorted.end(),
              [](const PorkchopPoint& a, const PorkchopPoint& b) {
                  return a.c3_depart < b.c3_depart;
              });

    printf("\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  TOP %d EARTH→MARS TRANSFER WINDOWS\n", n);
    printf("  (sorted by minimum C3 at departure)\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  %-6s %-10s %-8s %-10s %-8s %-8s\n",
           "Rank", "Departure", "TOF(d)", "Arrival", "C3", "ΔV_tot");
    printf("  %-6s %-10s %-8s %-10s %-8s %-8s\n",
           "------", "----------", "--------", "----------", "--------", "--------");

    int shown = 0;
    double last_c3 = -1;
    for (const auto& p : sorted) {
        // Skip near-duplicates (within 0.1 C3 of previous)
        if (shown > 0 && fabs(p.c3_depart - last_c3) < 0.01) continue;
        printf("  %-6d %10.2f %8.1f %10.2f %8.3f %8.3f\n",
               shown + 1,
               year_from_jd(p.jd_depart), p.tof_days,
               year_from_jd(p.jd_arrive),
               p.c3_depart, p.dv_total);
        last_c3 = p.c3_depart;
        shown++;
        if (shown >= n) break;
    }
    printf("\n");
}

void print_summary(const std::vector<PorkchopPoint>& points)
{
    if (points.empty()) { printf("  No valid points.\n"); return; }

    double min_c3 = 1e10, max_c3 = 0;
    double min_tof = 1e10, max_tof = 0;
    for (const auto& p : points) {
        if (p.c3_depart < min_c3) min_c3 = p.c3_depart;
        if (p.c3_depart > max_c3) max_c3 = p.c3_depart;
        if (p.tof_days  < min_tof) min_tof = p.tof_days;
        if (p.tof_days  > max_tof) max_tof = p.tof_days;
    }

    printf("\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  PORKCHOP SUMMARY — Earth→Mars %.0f–%.0f\n",
           DEPART_START_YEAR, DEPART_END_YEAR);
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  Valid transfer windows: %zu\n", points.size());
    printf("  C3 range:               %.3f – %.3f km²/s²\n", min_c3, max_c3);
    printf("  TOF range:              %.0f – %.0f days\n", min_tof, max_tof);
    printf("  Departure grid:         %.0f-day steps\n", DEPART_STEP_DAYS);
    printf("  TOF grid:               %.0f-day steps\n", TOF_STEP_DAYS);
    printf("\n");
    printf("  Plot with Python:\n");
    printf("    import pandas as pd, matplotlib.pyplot as plt\n");
    printf("    df = pd.read_csv('porkchop_earth_mars.csv')\n");
    printf("    plt.tricontourf(df['departure_year'], df['tof_days'], df['c3_depart'])\n");
    printf("    plt.colorbar(label='C3 (km²/s²)')\n");
    printf("    plt.xlabel('Departure Year')\n");
    printf("    plt.ylabel('Time of Flight (days)')\n");
    printf("    plt.title('Earth→Mars Porkchop Plot')\n");
    printf("    plt.show()\n");
    printf("\n");
}


///////////////////////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║  KEPLER WORKSTATION — Porkchop Plot Generator        ║\n");
    printf("║  Earth → Mars Transfer Windows %.0f–%.0f               ║\n",
           DEPART_START_YEAR, DEPART_END_YEAR);
    printf("╚══════════════════════════════════════════════════════╝\n");
    printf("\n");

    std::vector<PorkchopPoint> points;

#ifdef HAS_EPHEMERIS
    const char* eph_path = "../../Ephemeris/linux_p1550p2650.440";
    JplEphemeris eph;
    if (eph.load(eph_path)) {
        printf("  Using JPL DE440 ephemeris for real planetary positions.\n\n");
        compute_ephemeris_sweep(eph, points);
        eph.close();
    } else {
        printf("  Ephemeris not found: %s\n", eph_path);
        printf("  Falling back to circular-coplanar model.\n\n");
        compute_circular_sweep(points);
    }
#else
    printf("  Compiled without ephemeris support.\n");
    printf("  Using circular-coplanar model (approximate).\n");
    printf("  Recompile with -DHAS_EPHEMERIS for real positions.\n\n");
    compute_circular_sweep(points);
#endif

    if (points.empty()) {
        printf("  No valid transfer windows found.\n");
        return 1;
    }

    // Output CSV
    write_csv(points, "porkchop_earth_mars.csv");

    // Print summary and top windows
    print_summary(points);
    print_top_windows(points, 5);

    printf("══════════════════════════════════════════════════════════════\n");
    printf("  CSV saved to porkchop_earth_mars.csv\n");
    printf("  Ready for plotting or import into mission planning tools.\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");

    return 0;
}
