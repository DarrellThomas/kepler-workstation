// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Test: JPL Ephemeris Reader + N-Body Propagator
//
// 1. Load DE440 and print planet positions at J2000 epoch
// 2. Validate against known positions
// 3. Propagate a spacecraft from Earth orbit toward Jupiter
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include "jpl_ephemeris.hpp"
#include "nbody_propagator.hpp"

// Julian Date for J2000 epoch: 2000-Jan-1.5 TDB
static const double JD_J2000 = 2451545.0;

// Julian Date for a recent date: 2025-Jan-1.0
static const double JD_2025 = 2460676.5;

void print_vec(const char* label, const double v[3])
{
    printf("  %-12s %14.3f %14.3f %14.3f  |r|=%14.3f km\n",
           label, v[0], v[1], v[2],
           sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

int main()
{
    printf("=== JPL DE440 Ephemeris Test ===\n\n");

    JplEphemeris eph;
    if (!eph.load("ephemeris/linux_p1550p2650.440")) {
        printf("Failed to load ephemeris. Trying alternate path...\n");
        if (!eph.load("/data/src/zipfel/solar_system/ephemeris/linux_p1550p2650.440")) {
            printf("FATAL: Cannot load DE440\n");
            return 1;
        }
    }

    printf("\n--- Planet positions at J2000 epoch (JD %.1f) ---\n", JD_J2000);
    printf("  %-12s %14s %14s %14s  %14s\n", "Body", "X(km)", "Y(km)", "Z(km)", "|r|");

    double pos[3], vel[3];
    const char* names[] = {"Mercury","Venus","EMB","Mars","Jupiter",
                           "Saturn","Uranus","Neptune","Pluto","Moon","Sun"};

    for (int b = 1; b <= 11; b++) {
        if (eph.get_state(static_cast<JplBody>(b), JD_J2000, pos, vel)) {
            print_vec(names[b-1], pos);
        }
    }

    // Get Earth specifically
    printf("\n--- Earth (derived from EMB and Moon) ---\n");
    eph.get_earth(JD_J2000, pos, vel);
    print_vec("Earth", pos);
    printf("  Earth vel:   %10.6f %10.6f %10.6f  |v|=%.3f km/s\n",
           vel[0], vel[1], vel[2],
           sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]));

    // Known value: Earth-Sun distance at J2000 ~1 AU = 149,597,870.7 km
    double earth_pos[3], earth_vel[3];
    double sun_pos[3], sun_vel[3];
    eph.get_earth(JD_J2000, earth_pos, earth_vel);
    eph.get_state(JPL_SUN, JD_J2000, sun_pos, sun_vel);
    double dx = earth_pos[0]-sun_pos[0], dy = earth_pos[1]-sun_pos[1], dz = earth_pos[2]-sun_pos[2];
    double earth_sun_dist = sqrt(dx*dx+dy*dy+dz*dz);
    printf("\n  Earth-Sun distance: %.1f km (expected ~149,597,871 km = 1 AU)\n", earth_sun_dist);
    printf("  Ratio to AU: %.6f\n", earth_sun_dist / 149597870.7);

    // Jupiter distance from Sun
    double jup_pos[3], jup_vel[3];
    eph.get_state(JPL_JUPITER, JD_J2000, jup_pos, jup_vel);
    dx = jup_pos[0]-sun_pos[0]; dy = jup_pos[1]-sun_pos[1]; dz = jup_pos[2]-sun_pos[2];
    double jup_sun_dist = sqrt(dx*dx+dy*dy+dz*dz);
    printf("  Jupiter-Sun distance: %.1f km (expected ~5.2 AU = ~778M km)\n", jup_sun_dist);
    printf("  In AU: %.3f\n", jup_sun_dist / 149597870.7);

    // Test at 2025
    printf("\n--- Planet positions at 2025-Jan-1 (JD %.1f) ---\n", JD_2025);
    for (int b = 1; b <= 11; b++) {
        if (eph.get_state(static_cast<JplBody>(b), JD_2025, pos, vel)) {
            double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
            printf("  %-12s |r|=%14.1f km\n", names[b-1], r);
        }
    }

    printf("\n=== N-Body Propagator Test ===\n\n");

    NBodyPropagator prop;
    prop.set_ephemeris(&eph);
    prop.set_tolerances(1e-10, 1e-10);
    prop.set_step_limits(1.0, 3600.0); // 1s to 1h steps

    // Start spacecraft in Earth orbit (LEO, ~400 km altitude)
    // Position: Earth + 6778 km in X direction
    // Velocity: Earth velocity + 7.67 km/s in Y direction (circular orbit)
    SpacecraftState sc;
    sc.jd = JD_2025;
    sc.mass = 10000.0; // 10 tonnes

    eph.get_earth(JD_2025, earth_pos, earth_vel);
    sc.pos[0] = earth_pos[0] + 6778.0; // 400 km above surface
    sc.pos[1] = earth_pos[1];
    sc.pos[2] = earth_pos[2];
    sc.vel[0] = earth_vel[0];
    sc.vel[1] = earth_vel[1] + 7.67;   // circular orbit velocity
    sc.vel[2] = earth_vel[2];

    printf("Spacecraft initial state (LEO, 400 km):\n");
    print_vec("Position", sc.pos);
    printf("  Velocity:    %10.3f %10.3f %10.3f  |v|=%.3f km/s\n",
           sc.vel[0], sc.vel[1], sc.vel[2],
           sqrt(sc.vel[0]*sc.vel[0]+sc.vel[1]*sc.vel[1]+sc.vel[2]*sc.vel[2]));

    // Propagate for 1 orbit (~92 minutes = 5520 seconds)
    double jd_1orbit = sc.jd + 5520.0 / 86400.0;
    printf("\nPropagating 1 orbit (92 min)...\n");
    int steps = prop.propagate(sc, jd_1orbit);

    // Check: should be back near starting position relative to Earth
    eph.get_earth(sc.jd, earth_pos, earth_vel);
    double rel_pos[3], rel_vel[3];
    for (int i = 0; i < 3; i++) {
        rel_pos[i] = sc.pos[i] - earth_pos[i];
        rel_vel[i] = sc.vel[i] - earth_vel[i];
    }
    double alt = sqrt(rel_pos[0]*rel_pos[0]+rel_pos[1]*rel_pos[1]+rel_pos[2]*rel_pos[2]) - 6378.0;

    printf("After 1 orbit (%d steps):\n", steps);
    printf("  Altitude above Earth: %.1f km (started at ~400 km)\n", alt);
    printf("  Relative velocity: %.3f km/s (started at 7.67 km/s)\n",
           sqrt(rel_vel[0]*rel_vel[0]+rel_vel[1]*rel_vel[1]+rel_vel[2]*rel_vel[2]));

    // Propagate for 1 day
    printf("\nPropagating 1 day from LEO...\n");
    sc.jd = JD_2025;
    eph.get_earth(JD_2025, earth_pos, earth_vel);
    sc.pos[0] = earth_pos[0] + 6778.0;
    sc.pos[1] = earth_pos[1];
    sc.pos[2] = earth_pos[2];
    sc.vel[0] = earth_vel[0];
    sc.vel[1] = earth_vel[1] + 7.67;
    sc.vel[2] = earth_vel[2];

    double jd_1day = sc.jd + 1.0;
    steps = prop.propagate(sc, jd_1day);

    eph.get_earth(sc.jd, earth_pos, earth_vel);
    for (int i = 0; i < 3; i++) rel_pos[i] = sc.pos[i] - earth_pos[i];
    alt = sqrt(rel_pos[0]*rel_pos[0]+rel_pos[1]*rel_pos[1]+rel_pos[2]*rel_pos[2]) - 6378.0;

    printf("After 1 day (%d steps):\n", steps);
    printf("  Altitude above Earth: %.1f km (should still be ~400 km for stable orbit)\n", alt);
    printf("  Nearest body: %d, dist: %.1f km, in_soi: %d\n",
           sc.nearest_body, sc.dist_nearest, sc.in_soi);

    printf("\n=== Test Complete ===\n");
    return 0;
}
