///////////////////////////////////////////////////////////////////////////////
// EARTH-MARS HOHMANN TRANSFER DEMO
//
// Computes the Hohmann transfer orbit from Earth to Mars:
//   1. Ideal Hohmann (circular, coplanar) — the textbook answer
//   2. Real transfer using JPL DE440 ephemeris + Lambert solver
//
// Compares theory to reality. Educational tool for orbital mechanics.
//
// Compile:
//   g++ -std=c++11 -O2 -o hohmann_demo hohmann_demo.cpp
//   (requires ../../Ephemeris/linux_p1550p2650.440 for real transfer)
//
// 20260603 Created for Kepler Workstation OrbitalTransfer application
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////
// IDEAL HOHMANN TRANSFER (circular, coplanar)
///////////////////////////////////////////////////////////////////////////////

// Physical constants
static const double AU_KM       = 149597870.7;    // km per AU
static const double GM_SUN      = 1.32712440018e11; // km³/s²
static const double SEC_PER_DAY = 86400.0;

// Earth and Mars orbit data
static const double EARTH_SMA_AU = 1.00000011;   // semi-major axis (AU)
static const double MARS_SMA_AU  = 1.523679;     // semi-major axis (AU)

// Parking orbits
static const double R_EARTH_KM   = 6371.0;
static const double R_MARS_KM    = 3390.0;
static const double PARK_ALT_KM  = 300.0;
static const double MU_EARTH     = 3.986004418e5; // km³/s²
static const double MU_MARS      = 4.282837e4;    // km³/s²

struct IdealHohmann {
    double r_depart_km;    // departure orbit radius
    double r_arrive_km;    // arrival orbit radius
    double a_transfer_km;  // transfer orbit semi-major axis
    double v_depart_km_s;  // departure velocity from parking orbit
    double v_arrive_km_s;  // arrival velocity at target
    double dv_depart;      // departure ΔV (km/s)
    double dv_arrive;      // arrival ΔV (km/s)
    double dv_total;       // total ΔV (km/s)
    double tof_days;       // time of flight (days)
    double c3_depart;      // C3 at departure (km²/s²)
};

IdealHohmann compute_ideal_hohmann()
{
    IdealHohmann h;
    h.r_depart_km = EARTH_SMA_AU * AU_KM;
    h.r_arrive_km = MARS_SMA_AU * AU_KM;
    h.a_transfer_km = (h.r_depart_km + h.r_arrive_km) / 2.0;

    // Vis-viva: v = sqrt(mu * (2/r - 1/a))
    double v_earth_orbit = sqrt(GM_SUN / h.r_depart_km);  // Earth orbital velocity
    double v_mars_orbit  = sqrt(GM_SUN / h.r_arrive_km);   // Mars orbital velocity

    // Transfer orbit velocities at perihelion (departure) and aphelion (arrival)
    double v_trans_depart = sqrt(GM_SUN * (2.0/h.r_depart_km - 1.0/h.a_transfer_km));
    double v_trans_arrive = sqrt(GM_SUN * (2.0/h.r_arrive_km - 1.0/h.a_transfer_km));

    // Hyperbolic excess (v_infinity) at each end
    h.c3_depart     = (v_trans_depart - v_earth_orbit) * (v_trans_depart - v_earth_orbit);
    double v_inf_d  = v_trans_depart - v_earth_orbit;
    double v_inf_a  = v_mars_orbit - v_trans_arrive;  // positive if braking

    // ΔV from parking orbit (vis-viva for hyperbolic departure)
    double r_park_earth = R_EARTH_KM + PARK_ALT_KM;
    double r_park_mars  = R_MARS_KM  + PARK_ALT_KM;
    double v_circ_earth = sqrt(MU_EARTH / r_park_earth);
    double v_circ_mars  = sqrt(MU_MARS  / r_park_mars);

    double v_hyp_depart = sqrt(v_inf_d * v_inf_d + 2.0 * MU_EARTH / r_park_earth);
    double v_hyp_arrive = sqrt(v_inf_a * v_inf_a + 2.0 * MU_MARS  / r_park_mars);

    h.dv_depart = v_hyp_depart - v_circ_earth;
    h.dv_arrive = v_hyp_arrive - v_circ_mars;
    h.dv_total  = h.dv_depart + h.dv_arrive;
    h.v_depart_km_s = v_trans_depart;
    h.v_arrive_km_s = v_trans_arrive;

    // Time of flight = half the transfer orbit period
    double period_transfer = 2.0 * M_PI * sqrt(h.a_transfer_km*h.a_transfer_km*h.a_transfer_km / GM_SUN);
    h.tof_days = period_transfer / (2.0 * SEC_PER_DAY);

    return h;
}

void print_ideal_hohmann(const IdealHohmann& h)
{
    printf("\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  IDEAL HOHMANN TRANSFER — Earth to Mars\n");
    printf("  (Circular, coplanar — the textbook answer)\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  Assumptions:\n");
    printf("    Earth orbit: circular at %.6f AU\n", EARTH_SMA_AU);
    printf("    Mars  orbit: circular at %.6f AU\n", MARS_SMA_AU);
    printf("    Parking orbits: %d km altitude\n", (int)PARK_ALT_KM);
    printf("\n");
    printf("  Transfer orbit:\n");
    printf("    Perihelion (Earth):  %.3f million km\n", h.r_depart_km / 1e6);
    printf("    Aphelion (Mars):     %.3f million km\n", h.r_arrive_km / 1e6);
    printf("    Semi-major axis:     %.3f million km\n", h.a_transfer_km / 1e6);
    printf("    Eccentricity:        %.4f\n",
           (h.r_arrive_km - h.r_depart_km) / (h.r_arrive_km + h.r_depart_km));
    printf("\n");
    printf("  Velocities (heliocentric):\n");
    printf("    Earth orbital:       %.3f km/s\n", sqrt(GM_SUN / h.r_depart_km));
    printf("    Transfer at depart:  %.3f km/s\n", h.v_depart_km_s);
    printf("    Transfer at arrival: %.3f km/s\n", h.v_arrive_km_s);
    printf("    Mars orbital:        %.3f km/s\n", sqrt(GM_SUN / h.r_arrive_km));
    printf("\n");
    printf("  Performance:\n");
    printf("    C3 departure:        %.3f km²/s²\n", h.c3_depart);
    printf("    ΔV departure:        %.3f km/s\n", h.dv_depart);
    printf("    ΔV arrival:          %.3f km/s\n", h.dv_arrive);
    printf("    ΔV total:            %.3f km/s\n", h.dv_total);
    printf("    Time of flight:      %.1f days (%.2f years)\n",
           h.tof_days, h.tof_days / 365.25);
    printf("\n");
}


///////////////////////////////////////////////////////////////////////////////
// REAL TRANSFER — uses JPL DE440 ephemeris + Lambert solver
///////////////////////////////////////////////////////////////////////////////

#ifdef HAS_EPHEMERIS

#include "../../Ephemeris/jpl_ephemeris.hpp"
#include "../../Ephemeris/lambert.hpp"

// Helper: vector magnitude
static inline double vmag(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Helper: vector subtraction
static inline void vsub(double o[3], const double a[3], const double b[3]) {
    o[0] = a[0] - b[0]; o[1] = a[1] - b[1]; o[2] = a[2] - b[2];
}

// Helper: Julian date from calendar year
static double jd_from_year(double year) {
    return 2451545.0 + (year - 2000.0) * 365.25;
}

// Helper: calendar year from JD
static double year_from_jd(double jd) {
    return 2000.0 + (jd - 2451545.0) / 365.25;
}

// Best transfer window result
struct BestWindow {
    double jd_depart;
    double jd_arrive;
    double tof_days;
    double v_inf_depart;
    double v_inf_arrive;
    double c3_depart;
    double dv_total;
    bool   valid;
};

BestWindow find_best_window(JplEphemeris& eph,
                             double jd_start, double jd_end,
                             double tof_min, double tof_max, double tof_step)
{
    BestWindow best;
    best.valid = false;
    best.c3_depart = 1e10;

    double sp[3], sv[3];  // Sun position/velocity

    for (double jd_d = jd_start; jd_d <= jd_end; jd_d += 30.0) {  // 30-day departure steps
        double ep[3], ev[3];  // Earth
        double mp[3], mv[3];  // Mars

        if (!eph.get_earth(jd_d, ep, ev)) continue;

        for (double tof = tof_min; tof <= tof_max; tof += tof_step) {
            double jd_a = jd_d + tof;

            if (!eph.get_state(JPL_MARS, jd_a, mp, mv)) continue;

            // Get Sun-relative positions
            eph.get_state(JPL_SUN, jd_d, sp, sv);
            double r1[3], rv1[3];
            vsub(r1, ep, sp);
            vsub(rv1, ev, sv);

            eph.get_state(JPL_SUN, jd_a, sp, sv);
            double r2[3], rv2[3];
            vsub(r2, mp, sp);
            vsub(rv2, mv, sv);

            // Solve Lambert
            double vt1[3], vt2[3];
            if (!lambert_solve(r1, r2, tof * SEC_PER_DAY, GM_SUN, true, vt1, vt2))
                continue;

            // Compute v_inf
            double vd[3], va[3];
            vsub(vd, vt1, rv1);
            vsub(va, vt2, rv2);
            double vid = vmag(vd);
            double via = vmag(va);

            double c3 = vid * vid;

            // Only consider reasonable transfers (C3 < 50 km²/s²)
            if (c3 < 50.0 && c3 < best.c3_depart) {
                best.c3_depart   = c3;
                best.v_inf_depart = vid;
                best.v_inf_arrive = via;
                best.tof_days     = tof;
                best.jd_depart    = jd_d;
                best.jd_arrive    = jd_a;

                // Approximate total ΔV (from parking orbits)
                double r_park_earth = R_EARTH_KM + PARK_ALT_KM;
                double r_park_mars  = R_MARS_KM  + PARK_ALT_KM;
                double v_circ_e = sqrt(MU_EARTH / r_park_earth);
                double v_circ_m = sqrt(MU_MARS  / r_park_mars);
                double v_hyp_d  = sqrt(vid*vid + 2.0*MU_EARTH/r_park_earth);
                double v_hyp_a  = sqrt(via*via + 2.0*MU_MARS /r_park_mars);
                best.dv_total    = (v_hyp_d - v_circ_e) + (v_hyp_a - v_circ_m);
                best.valid       = true;
            }
        }
    }
    return best;
}

void print_real_transfer(JplEphemeris& eph)
{
    printf("\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  REAL TRANSFER — JPL DE440 Ephemeris + Lambert Solver\n");
    printf("  Searching %.0f–%.0f for best Earth→Mars window\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");

    // Search 2026–2037 (spans several synodic periods)
    double jd_start = jd_from_year(2026.0);
    double jd_end   = jd_from_year(2037.0);

    printf("  Searching %.0f departure dates × %.0f TOF values...\n",
           (jd_end - jd_start) / 30.0, (400.0 - 50.0) / 10.0);

    BestWindow best = find_best_window(eph, jd_start, jd_end, 50.0, 400.0, 10.0);

    if (!best.valid) {
        printf("  No valid transfer found in search range.\n");
        return;
    }

    double dep_year = year_from_jd(best.jd_depart);
    double arr_year = year_from_jd(best.jd_arrive);

    printf("\n");
    printf("  Best transfer window found:\n");
    printf("    Departure:  JD %.1f  (~%.2f)\n", best.jd_depart, dep_year);
    printf("    Arrival:    JD %.1f  (~%.2f)\n", best.jd_arrive, arr_year);
    printf("    TOF:        %.1f days (%.2f years)\n",
           best.tof_days, best.tof_days / 365.25);
    printf("\n");
    printf("  Performance:\n");
    printf("    v_inf departure:  %.3f km/s\n", best.v_inf_depart);
    printf("    v_inf arrival:    %.3f km/s\n", best.v_inf_arrive);
    printf("    C3 departure:     %.3f km²/s²\n", best.c3_depart);
    printf("    ΔV total:         %.3f km/s\n", best.dv_total);
    printf("\n");
}

#endif // HAS_EPHEMERIS


///////////////////////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║  KEPLER WORKSTATION — Orbital Transfer Demo          ║\n");
    printf("║  Earth → Mars Hohmann Transfer                       ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");

    // 1. Ideal Hohmann (always available)
    IdealHohmann ideal = compute_ideal_hohmann();
    print_ideal_hohmann(ideal);

    // 2. Real transfer (requires ephemeris file)
#ifdef HAS_EPHEMERIS
    const char* eph_path = "../../Ephemeris/linux_p1550p2650.440";

    JplEphemeris eph;
    if (eph.load(eph_path)) {
        print_real_transfer(eph);

        // Comparison
        printf("══════════════════════════════════════════════════════════════\n");
        printf("  COMPARISON: Ideal vs Real\n");
        printf("══════════════════════════════════════════════════════════════\n");
        printf("\n");
        printf("  %-25s %12s %12s\n", "Parameter", "Ideal", "Real");
        printf("  %-25s %12s %12s\n", "-------------------------", "------------", "------------");
        printf("  %-25s %10.3f   %10s\n", "C3 departure (km²/s²)",
               ideal.c3_depart, "(varies)");
        printf("  %-25s %10.3f   %10s\n", "ΔV total (km/s)",
               ideal.dv_total, "(varies)");
        printf("  %-25s %10.1f   %10s\n", "TOF (days)",
               ideal.tof_days, "(varies)");
        printf("\n");
        printf("  The ideal Hohmann assumes circular, coplanar orbits.\n");
        printf("  Real transfers vary with planetary positions — use\n");
        printf("  porkchop_plot for the full window analysis.\n");

        eph.close();
    } else {
        printf("\n");
        printf("  [Ephemeris file not found: %s]\n", eph_path);
        printf("  Download from: https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/\n");
        printf("  Recompile with: g++ -std=c++11 -O2 -DHAS_EPHEMERIS -o hohmann_demo hohmann_demo.cpp\n");
    }
#else
    printf("\n");
    printf("  [Real transfer not compiled — add -DHAS_EPHEMERIS and download DE440]\n");
    printf("  [Showing ideal Hohmann only.]\n");
#endif

    printf("\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  Next: Run porkchop_plot to see Earth→Mars launch windows\n");
    printf("  over the next decade with full C3 contour data.\n");
    printf("══════════════════════════════════════════════════════════════\n");
    printf("\n");

    return 0;
}
