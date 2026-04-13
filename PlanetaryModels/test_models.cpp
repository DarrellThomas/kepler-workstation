///////////////////////////////////////////////////////////////////////////////
// Validation test for planetary environment models
// Compares outputs against known reference values
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include "planetary_environment.hpp"

int main()
{
    double rho, press, tempk;

    printf("=== VENUS ATMOSPHERE ===\n");
    printf("%-12s %12s %12s %12s %12s\n", "Alt(km)", "T(K)", "P(Pa)", "rho(kg/m3)", "Vsound(m/s)");
    double venus_alts[] = {0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 100000};
    for (double alt : venus_alts)
    {
        atmosphere_venus(rho, press, tempk, alt);
        double vs = speed_of_sound_venus(tempk);
        printf("%-12.0f %12.2f %12.0f %12.4f %12.2f\n", alt/1000, tempk, press, rho, vs);
    }
    printf("\nVenus surface validation:\n");
    atmosphere_venus(rho, press, tempk, 0);
    printf("  T=%.1f K (expected 737 K)\n", tempk);
    printf("  P=%.0f Pa = %.1f bar (expected 92 bar)\n", press, press/1e5);
    printf("  rho=%.2f kg/m3 (expected ~65)\n", rho);

    printf("\n=== MARS ATMOSPHERE ===\n");
    printf("%-12s %12s %12s %12s %12s\n", "Alt(km)", "T(K)", "P(Pa)", "rho(kg/m3)", "Vsound(m/s)");

    double mars_alts[] = {0, 1000, 5000, 7000, 10000, 20000, 40000, 60000, 80000};
    for (double alt : mars_alts)
    {
        atmosphere_mars(rho, press, tempk, alt);
        double vs = speed_of_sound_mars(tempk);
        printf("%-12.0f %12.2f %12.2f %12.6f %12.2f\n", alt/1000, tempk, press, rho, vs);
    }

    // Known reference: Mars surface ~T=210K, P=636Pa, rho=0.020 kg/m3
    atmosphere_mars(rho, press, tempk, 0);
    printf("\nMars surface validation:\n");
    printf("  T=%.1f K (expected ~210 K)\n", tempk);
    printf("  P=%.1f Pa (expected ~636 Pa)\n", press);
    printf("  rho=%.4f kg/m3 (expected ~0.020)\n", rho);

    printf("\n=== JUPITER ATMOSPHERE ===\n");
    printf("%-12s %12s %12s %12s %12s\n", "Alt(km)", "T(K)", "P(Pa)", "rho(kg/m3)", "Vsound(m/s)");

    // Altitudes above 1-bar level (can be negative for deep atmosphere)
    double jup_alts[] = {-130000, -100000, -60000, -30000, 0, 10000, 30000, 50000, 100000, 200000};
    for (double alt : jup_alts)
    {
        atmosphere_jupiter(rho, press, tempk, alt);
        double vs = speed_of_sound_jupiter(tempk);
        printf("%-12.0f %12.2f %12.2f %12.6f %12.2f\n", alt/1000, tempk, press, rho, vs);
    }

    // Galileo probe reference values:
    // At 1 bar (h=0): T~166 K, P=1e5 Pa
    // At 24 bars (h~-133 km): T~428 K, P~24e5 Pa
    atmosphere_jupiter(rho, press, tempk, 0);
    printf("\nJupiter 1-bar level validation:\n");
    printf("  T=%.1f K (Galileo: 166 K)\n", tempk);
    printf("  P=%.0f Pa (expected: 100000 Pa)\n", press);

    atmosphere_jupiter(rho, press, tempk, -133000);
    printf("\nJupiter ~24-bar level (-133 km) validation:\n");
    printf("  T=%.1f K (Galileo: ~428 K)\n", tempk);

    printf("\n=== EARTH GRAVITY (upgraded J2-J6) ===\n");
    double earth_g[3];

    // Equator surface
    double earth_eq[] = {EARTH_RADIUS, 0, 0};
    gravity_earth(earth_g, earth_eq);
    double g_eq = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    printf("Equator surface: g=%.6f m/s2 (expected ~9.7803)\n", g_eq);

    // Pole surface (polar radius ~6356752 m)
    double earth_pole[] = {0, 0, 6356752.0};
    gravity_earth(earth_g, earth_pole);
    double g_pole = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    printf("Pole surface:    g=%.6f m/s2 (expected ~9.8322)\n", g_pole);

    // LEO (400 km)
    double earth_leo[] = {EARTH_RADIUS + 400000, 0, 0};
    gravity_earth(earth_g, earth_leo);
    double g_leo = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    printf("400 km LEO:      g=%.6f m/s2 (expected ~8.69)\n", g_leo);

    // GEO (35786 km)
    double earth_geo[] = {EARTH_RADIUS + 35786000, 0, 0};
    gravity_earth(earth_g, earth_geo);
    double g_geo = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    printf("GEO (35786 km):  g=%.6f m/s2 (expected ~0.224)\n", g_geo);

    // J3 asymmetry test: compare 45°N vs 45°S
    double r45 = EARTH_RADIUS;
    double lat45 = 45.0 * M_PI / 180.0;
    double earth_n45[] = {r45*cos(lat45), 0, r45*sin(lat45)};
    double earth_s45[] = {r45*cos(lat45), 0, -r45*sin(lat45)};
    gravity_earth(earth_g, earth_n45);
    double g_n45 = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    gravity_earth(earth_g, earth_s45);
    double g_s45 = sqrt(earth_g[0]*earth_g[0]+earth_g[1]*earth_g[1]+earth_g[2]*earth_g[2]);
    printf("45N vs 45S:      g_N=%.6f g_S=%.6f diff=%.3e m/s2 (J3 pear shape)\n",
           g_n45, g_s45, g_n45-g_s45);

    printf("\n=== VENUS GRAVITY ===\n");
    double venus_g[3];
    double venus_eq[] = {VENUS_RADIUS, 0, 0};
    gravity_venus(venus_g, venus_eq);
    printf("Equator surface: g=%.4f m/s2 (expected ~8.87)\n",
           sqrt(venus_g[0]*venus_g[0]+venus_g[1]*venus_g[1]+venus_g[2]*venus_g[2]));

    printf("\n=== MOON GRAVITY (GRAIL) ===\n");
    double moon_g[3];
    double moon_eq[] = {MOON_RADIUS, 0, 0};
    gravity_moon(moon_g, moon_eq);
    printf("Equator surface: g=%.4f m/s2 (expected ~1.62)\n",
           sqrt(moon_g[0]*moon_g[0]+moon_g[1]*moon_g[1]+moon_g[2]*moon_g[2]));
    double moon_pole[] = {0, 0, MOON_RADIUS};
    gravity_moon(moon_g, moon_pole);
    printf("Pole surface:    g=%.4f m/s2\n",
           sqrt(moon_g[0]*moon_g[0]+moon_g[1]*moon_g[1]+moon_g[2]*moon_g[2]));

    printf("\n=== MARS GRAVITY ===\n");
    double mars_g[3];

    // Surface at equator
    double mars_eq[] = {MARS_RADIUS, 0, 0};
    gravity_mars(mars_g, mars_eq);
    printf("Equator surface: g=%.4f m/s2 (expected ~3.71)\n",
           sqrt(mars_g[0]*mars_g[0]+mars_g[1]*mars_g[1]+mars_g[2]*mars_g[2]));

    // Surface at pole
    double mars_pole[] = {0, 0, MARS_RADIUS};
    gravity_mars(mars_g, mars_pole);
    printf("Pole surface:    g=%.4f m/s2 (expected ~3.76, stronger due to J2)\n",
           sqrt(mars_g[0]*mars_g[0]+mars_g[1]*mars_g[1]+mars_g[2]*mars_g[2]));

    // At 100 km altitude
    double mars_100km[] = {MARS_RADIUS + 100000, 0, 0};
    gravity_mars(mars_g, mars_100km);
    printf("100 km equator:  g=%.4f m/s2\n",
           sqrt(mars_g[0]*mars_g[0]+mars_g[1]*mars_g[1]+mars_g[2]*mars_g[2]));

    printf("\n=== JUPITER GRAVITY ===\n");
    double jup_g[3];

    // "Surface" at 1-bar level (equatorial radius)
    double jup_eq[] = {JUPITER_RADIUS, 0, 0};
    gravity_jupiter(jup_g, jup_eq);
    printf("Equator 1-bar:   g=%.4f m/s2 (expected ~24.79)\n",
           sqrt(jup_g[0]*jup_g[0]+jup_g[1]*jup_g[1]+jup_g[2]*jup_g[2]));

    // Pole
    double jup_pole[] = {0, 0, JUPITER_RADIUS * 0.9353}; // polar radius ~66854 km
    gravity_jupiter(jup_g, jup_pole);
    printf("Pole surface:    g=%.4f m/s2 (expected ~28.4, stronger at pole)\n",
           sqrt(jup_g[0]*jup_g[0]+jup_g[1]*jup_g[1]+jup_g[2]*jup_g[2]));

    // At 1000 km altitude
    double jup_1000[] = {JUPITER_RADIUS + 1000000, 0, 0};
    gravity_jupiter(jup_g, jup_1000);
    printf("1000 km equator: g=%.4f m/s2\n",
           sqrt(jup_g[0]*jup_g[0]+jup_g[1]*jup_g[1]+jup_g[2]*jup_g[2]));

    printf("\n=== COMPARISON: ALL BODIES (surface/reference level) ===\n");
    printf("%-12s %10s %10s %10s %10s %10s\n", "Property", "Venus", "Earth", "Moon", "Mars", "Jupiter");

    double e_rho=1.225, e_press=101325, e_tempk=288.15;
    double v_rho, v_press, v_tempk;
    double m_rho, m_press, m_tempk;
    double j_rho, j_press, j_tempk;
    atmosphere_venus(v_rho, v_press, v_tempk, 0);
    atmosphere_mars(m_rho, m_press, m_tempk, 0);
    atmosphere_jupiter(j_rho, j_press, j_tempk, 0);

    printf("%-12s %10.1f %10.1f %10s %10.1f %10.1f\n", "T(K)", v_tempk, e_tempk, "N/A", m_tempk, j_tempk);
    printf("%-12s %10.0f %10.0f %10s %10.0f %10.0f\n", "P(Pa)", v_press, e_press, "0", m_press, j_press);
    printf("%-12s %10.2f %10.4f %10s %10.4f %10.4f\n", "rho(kg/m3)", v_rho, e_rho, "0", m_rho, j_rho);
    printf("%-12s %10.2f %10.2f %10.2f %10.2f %10.2f\n", "g(m/s2)", VENUS_G0, 9.81, MOON_G0, MARS_G0, JUPITER_G0);
    printf("%-12s %10.1f %10.1f %10.1f %10.1f %10.1f\n", "R(km)", VENUS_RADIUS/1e3, EARTH_RADIUS/1e3, MOON_RADIUS/1e3, MARS_RADIUS/1e3, JUPITER_RADIUS/1e3);
    printf("%-12s %10.3e %10.3e %10.3e %10.3e %10.3e\n", "J2", VENUS_J2, EARTH_J2, MOON_J2, MARS_J2, JUPITER_J2);

    return 0;
}
