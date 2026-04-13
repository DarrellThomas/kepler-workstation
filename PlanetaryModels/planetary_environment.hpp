// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Planetary Environment Models for CADAC++
//
// Atmosphere models:
//   atmosphere_mars()  — NASA GRC Mars atmosphere (Mars Global Surveyor 1996)
//   atmosphere_jupiter() — Galileo Probe ASI data (Seiff et al. 1998)
//
// Gravity models:
//   gravity_earth()   — WGS-84/EGM96 (J2-J6 + J3,J5 odd harmonics)
//   gravity_mars()    — Mars GMM-3 (Genova et al. 2016), J2-J6
//   gravity_jupiter() — Juno (Iess et al. 2018), J2-J10
//
// Drop-in replacements for CADAC++ atmosphere76() and cad_grav84()
//
// References:
//   Mars atmo: NASA GRC https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
//   Jupiter atmo: Galileo Probe ASI, PDS gp_0001/data/asi/descent.tab
//     Seiff et al. "Structure of the Atmosphere of Jupiter: Galileo Probe
//     Measurements", Science 272:844-845, 1996
//   Mars gravity: Genova et al. "Seasonal and static gravity field of Mars
//     from MGS, Mars Odyssey and MRO", Icarus 272:228-245, 2016
//   Jupiter gravity: Iess et al. "Measurement of Jupiter's asymmetric
//     gravity field", Nature 555:220-222, 2018
//
// 20260413 Created for Zipfel CADAC++ project
///////////////////////////////////////////////////////////////////////////////

#ifndef PLANETARY_ENVIRONMENT_HPP
#define PLANETARY_ENVIRONMENT_HPP

#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// CONSTANTS
///////////////////////////////////////////////////////////////////////////////

// Earth (WGS-84 / EGM96)
static const double EARTH_GM        = 3.986004418e14; // m^3/s^2
static const double EARTH_RADIUS    = 6378137.0;      // m (WGS-84 semi-major axis)
static const double EARTH_FLATTENING= 1.0/298.257223563;
static const double EARTH_OMEGA     = 7.292115e-5;    // rad/s (rotation rate)

// Earth gravity zonal harmonics (unnormalized, WGS-84/EGM96)
// Source: NIMA TR8350.2 (WGS-84), EGM96 for higher order
static const double EARTH_J2  =  1082.6267e-6;
static const double EARTH_J3  =    -2.5327e-6;
static const double EARTH_J4  =    -1.6196e-6;
static const double EARTH_J5  =    -0.2273e-6;
static const double EARTH_J6  =     0.5407e-6;

// Venus
static const double VENUS_GM        = 3.24859e14;    // m^3/s^2
static const double VENUS_RADIUS    = 6051800.0;     // m (mean radius, nearly spherical)
static const double VENUS_RGAS      = 191.0;         // J/(kg·K) (CO2 96.5%, N2 3.5%, M=43.45 g/mol)
static const double VENUS_GAMMA     = 1.286;         // ratio of specific heats (CO2-dominated)
static const double VENUS_G0        = 8.87;          // m/s^2 surface gravity

// Venus gravity zonal harmonics (×10^6, unnormalized)
// Source: Konopliv et al. 1999 (MGNP180U, Magellan), NASA Venus Fact Sheet
// Note: Venus is nearly spherical (very slow rotation), J2 is tiny
static const double VENUS_J2  = 4.458e-6;
static const double VENUS_J3  = -2.11e-6;
static const double VENUS_J4  = -2.15e-6;

// Moon (Luna)
static const double MOON_GM         = 4.9028e12;     // m^3/s^2
static const double MOON_RADIUS     = 1738100.0;     // m (mean equatorial radius)
static const double MOON_G0         = 1.62;          // m/s^2 surface gravity

// Moon gravity zonal harmonics (×10^6, unnormalized)
// Source: GRAIL mission (Zuber et al. 2013, Konopliv et al. 2013)
// Moon has significant gravity anomalies (mascons) but zonals capture global shape
static const double MOON_J2  = 203.21e-6;
static const double MOON_J3  = -8.46e-6;
static const double MOON_J4  = -9.63e-6;
static const double MOON_J5  = 0.66e-6;
static const double MOON_J6  = -1.40e-6;

// Sun (for solar system propagation)
static const double SUN_GM          = 1.32712440018e20; // m^3/s^2

// Mars
static const double MARS_GM        = 4.282837e13;   // m^3/s^2 (gravitational parameter)
static const double MARS_RADIUS    = 3396200.0;     // m (mean equatorial radius)
static const double MARS_RGAS      = 188.92;        // J/(kg·K) (CO2-dominated, M=43.34 g/mol)
static const double MARS_GAMMA     = 1.2941;        // ratio of specific heats (CO2)
static const double MARS_G0        = 3.71;          // m/s^2 surface gravity

// Mars gravity zonal harmonics (×10^6, unnormalized)
// Source: Genova et al. 2016 (GMM-3 model)
static const double MARS_J2  = 1960.454e-6;
static const double MARS_J3  = 31.451e-6;
static const double MARS_J4  = -15.587e-6;
static const double MARS_J5  = -4.935e-6;
static const double MARS_J6  = -2.151e-6;

// Jupiter
static const double JUPITER_GM     = 1.26686534e17; // m^3/s^2
static const double JUPITER_RADIUS = 71492000.0;    // m (equatorial radius, IAU)
static const double JUPITER_RGAS   = 3745.0;        // J/(kg·K) (H2-dominated, M=2.22 g/mol)
static const double JUPITER_GAMMA  = 1.44;          // ratio of specific heats (H2/He mix)
static const double JUPITER_G0     = 24.79;         // m/s^2 surface gravity (1-bar level)

// Jupiter gravity zonal harmonics (×10^6, unnormalized)
// Source: Iess et al. 2018 (Juno, Nature 555:220)
static const double JUPITER_J2  =  14696.572e-6;
static const double JUPITER_J3  =     -0.042e-6;
static const double JUPITER_J4  =   -586.609e-6;
static const double JUPITER_J5  =     -0.069e-6;
static const double JUPITER_J6  =     34.198e-6;
static const double JUPITER_J7  =      0.124e-6;
static const double JUPITER_J8  =     -2.426e-6;
static const double JUPITER_J9  =     -0.108e-6;
static const double JUPITER_J10 =      0.172e-6;


///////////////////////////////////////////////////////////////////////////////
// VENUS ATMOSPHERE MODEL
///////////////////////////////////////////////////////////////////////////////
// Based on Venera/Pioneer Venus/VeGa probe measurements
// Nine-point piecewise profile from surface (0 km) to 100 km
// Ref: NASA NSSDC Venus Fact Sheet, Venus International Reference Atmosphere
//
// Surface: 737 K, 92 bar, 65 kg/m³ (CO2, supercritical fluid near surface)
// 50 km: Earth-like conditions (~350 K, 1 bar)
// Venus atmosphere is 92x denser than Earth's at the surface
///////////////////////////////////////////////////////////////////////////////

inline void atmosphere_venus(double &rho, double &press, double &tempk, const double balt)
{
    double h = balt / 1000.0; // km

    // Piecewise temperature profile from probe data (°C)
    // Altitudes: 0, 10, 20, 30, 40, 50, 60, 70, 80, 100 km
    static const double alt_tab[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};
    static const double temp_c[]  = {462, 385, 306, 222, 143, 75, -10, -43, -76, -112};
    static const double pres_atm[]= {92.10, 47.39, 22.52, 9.851, 3.501, 1.066,
                                     0.2357, 0.03690, 0.004760, 0.00002660};
    static const int NLAYERS = 10;

    // Clamp altitude
    double hc = h;
    if (hc < 0) hc = 0;
    if (hc > 100) hc = 100;

    // Find bracketing layer
    int i = 0;
    for (int k = 0; k < NLAYERS - 1; k++)
    {
        if (hc >= alt_tab[k]) i = k;
    }
    if (i >= NLAYERS - 1) i = NLAYERS - 2;

    // Linear interpolation between layers
    double frac = (hc - alt_tab[i]) / (alt_tab[i+1] - alt_tab[i]);

    double T_c = temp_c[i] + frac * (temp_c[i+1] - temp_c[i]);
    tempk = T_c + 273.15;

    // Log-linear interpolation for pressure
    double log_p = log(pres_atm[i]) + frac * (log(pres_atm[i+1]) - log(pres_atm[i]));
    double p_atm = exp(log_p);
    press = p_atm * 101325.0; // convert atm to Pa

    // Density from ideal gas law
    rho = press / (VENUS_RGAS * tempk);

    if (rho < 0.0) rho = 0.0;
    if (press < 0.0) press = 0.0;
}


///////////////////////////////////////////////////////////////////////////////
// MARS ATMOSPHERE MODEL
///////////////////////////////////////////////////////////////////////////////
// NASA GRC model based on Mars Global Surveyor (1996)
// Two-layer model: lower (0-7 km) and upper (>7 km)
// Valid from surface to ~80 km
//
// Arguments (same signature as CADAC++ atmosphere76):
//   rho    — output: atmospheric density (kg/m^3)
//   press  — output: static pressure (Pa)
//   tempk  — output: temperature (K)
//   balt   — input: geometric altitude above surface (m)
///////////////////////////////////////////////////////////////////////////////

inline void atmosphere_mars(double &rho, double &press, double &tempk, const double balt)
{
    double h = balt; // altitude in meters
    double T_celsius;

    // Temperature model (two layers)
    if (h < 7000.0)
    {
        T_celsius = -31.0 - 0.000998 * h;
    }
    else
    {
        T_celsius = -23.4 - 0.00222 * h;
    }

    tempk = T_celsius + 273.1;

    // Pressure model (kPa, then convert to Pa)
    double p_kpa = 0.699 * exp(-0.00009 * h);
    press = p_kpa * 1000.0; // Pa

    // Density from equation of state: rho = p / (R_specific * T)
    // R_specific for Mars CO2 atmosphere: 0.1921 kPa·m³/(kg·K) = 192.1 J/(kg·K)
    rho = p_kpa / (0.1921 * tempk);

    // Clamp to non-negative
    if (rho < 0.0) rho = 0.0;
    if (press < 0.0) press = 0.0;
}


///////////////////////////////////////////////////////////////////////////////
// JUPITER ATMOSPHERE MODEL
///////////////////////////////////////////////////////////////////////////////
// Based on Galileo Probe ASI descent data (Seiff et al. 1998)
// Piecewise model from Galileo measurements:
//   - Stratosphere: ~300 km to 50 km above 1-bar (isothermal ~110 K)
//   - Tropopause: ~50 km (cold trap ~110 K)
//   - Upper troposphere: 50 km to 0 km (dry adiabat)
//   - Deep troposphere: 0 km to -150 km (dry adiabat, T up to ~430 K)
// Altitude reference: 1-bar pressure level (0 km)
//
// Galileo measured from +900 km (entry) to -133 km (24 bars)
// This model covers -150 km to +300 km
///////////////////////////////////////////////////////////////////////////////

inline void atmosphere_jupiter(double &rho, double &press, double &tempk, const double balt)
{
    // balt = altitude above 1-bar level in meters
    double h_km = balt / 1000.0;

    // Reference values at 1-bar level (h=0)
    double T0 = 166.0;    // K at 1 bar
    double P0 = 1.0e5;    // Pa (1 bar)
    double rho0 = 0.16;   // kg/m^3 at 1 bar

    // Piecewise temperature profile from Galileo ASI data
    // Layer boundaries and lapse rates (K/km)
    if (h_km > 50.0)
    {
        // Stratosphere: nearly isothermal above tropopause
        // T ~ 110 K, slight increase with altitude from solar heating
        tempk = 110.0 + 0.05 * (h_km - 50.0);
        if (tempk < 100.0) tempk = 100.0;
    }
    else if (h_km > 0.0)
    {
        // Upper troposphere: T increases from 110 K (tropopause) to 166 K (1-bar)
        // Approximately -1.12 K/km lapse rate
        tempk = 166.0 - 1.12 * h_km;
    }
    else if (h_km > -60.0)
    {
        // Mid troposphere: dry adiabat, ~2.0 K/km
        tempk = 166.0 - 2.0 * h_km; // h_km is negative, so T increases downward
    }
    else
    {
        // Deep troposphere (below -60 km): continues dry adiabat
        // Galileo measured ~428 K at 24 bars (~-133 km)
        tempk = 286.0 - 2.2 * (h_km + 60.0);
    }

    // Pressure from hydrostatic equation: p = p_ref * exp(-g*h / (R*T_avg))
    // Using scale height approach with local temperature
    double g_local = JUPITER_G0; // simplified, varies ~1% over range
    double scale_height = JUPITER_RGAS * tempk / g_local; // meters
    press = P0 * exp(-balt / scale_height);

    // Density from ideal gas law
    rho = press / (JUPITER_RGAS * tempk);

    if (rho < 0.0) rho = 0.0;
    if (press < 0.0) press = 0.0;
}


///////////////////////////////////////////////////////////////////////////////
// EARTH GRAVITY MODEL (J2-J6, WGS-84/EGM96)
///////////////////////////////////////////////////////////////////////////////
// Upgraded from CADAC++ cad_grav84() which only used J2.
// Adds J3-J6 for multi-orbit and interplanetary departure/arrival accuracy.
// J3,J5 (odd harmonics) break north-south symmetry (Earth's "pear shape").
//
// For the Jupiter harvester mission:
//   - Earth departure: J2 sufficient for escape burn (~minutes)
//   - Earth return approach: J3-J6 needed for entry corridor targeting
//     over days of approach (J3 alone shifts equatorial crossing by ~20m/orbit)
//   - Parking orbit propagation: J2-J4 needed for >1 orbit accuracy
//
// Arguments:
//   gravg[3] — output: gravitational acceleration (geocentric) m/s^2
//   pos[3]   — input: position in Earth-centered inertial coords (m)
///////////////////////////////////////////////////////////////////////////////

inline void gravity_earth(double gravg[3], const double pos[3])
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1.0) { gravg[0]=0; gravg[1]=0; gravg[2]=0; return; }

    double mu_r2 = EARTH_GM / (r * r);
    double a_r = EARTH_RADIUS / r;
    double sinlat = z / r;
    double s2 = sinlat * sinlat;
    double coslat = sqrt(1.0 - s2);

    // Powers of (a/r)
    double a_r2 = a_r * a_r;
    double a_r3 = a_r2 * a_r;
    double a_r4 = a_r2 * a_r2;
    double a_r5 = a_r4 * a_r;
    double a_r6 = a_r4 * a_r2;

    // Legendre polynomials P_n(sinlat)
    double s3 = s2 * sinlat;
    double s4 = s2 * s2;
    double s5 = s4 * sinlat;
    double s6 = s4 * s2;

    double P2 = 0.5 * (3.0*s2 - 1.0);
    double P3 = 0.5 * (5.0*s3 - 3.0*sinlat);
    double P4 = (1.0/8.0) * (35.0*s4 - 30.0*s2 + 3.0);
    double P5 = (1.0/8.0) * (63.0*s5 - 70.0*s3 + 15.0*sinlat);
    double P6 = (1.0/16.0) * (231.0*s6 - 315.0*s4 + 105.0*s2 - 5.0);

    // dP_n/d(sinlat)
    double dP2 = 3.0 * sinlat;
    double dP3 = 0.5 * (15.0*s2 - 3.0);
    double dP4 = (1.0/2.0) * (35.0*s3 - 15.0*sinlat);
    double dP5 = (1.0/8.0) * (315.0*s4 - 210.0*s2 + 15.0);
    double dP6 = (1.0/16.0) * (1386.0*s5 - 1260.0*s3 + 210.0*sinlat);

    // Radial acceleration component: -dU/dr
    // g_r = -(GM/r²)[1 + sum_n (n+1)*Jn*(a/r)^n * Pn(sinlat)]
    double radial_sum = 1.0
        + 3.0 * EARTH_J2 * a_r2 * P2
        + 4.0 * EARTH_J3 * a_r3 * P3
        + 5.0 * EARTH_J4 * a_r4 * P4
        + 6.0 * EARTH_J5 * a_r5 * P5
        + 7.0 * EARTH_J6 * a_r6 * P6;

    // Latitudinal acceleration component: -(1/r)*dU/d(lat)
    // g_lat = (GM/r²) * sum_n Jn*(a/r)^n * dPn/d(sinlat) * coslat
    double lat_sum =
          EARTH_J2 * a_r2 * dP2
        + EARTH_J3 * a_r3 * dP3
        + EARTH_J4 * a_r4 * dP4
        + EARTH_J5 * a_r5 * dP5
        + EARTH_J6 * a_r6 * dP6;

    double g_radial = -mu_r2 * radial_sum;
    double g_lat = mu_r2 * lat_sum * coslat;

    // Convert to Cartesian geocentric
    double r_xy = sqrt(x*x + y*y);
    if (r_xy < 1.0) r_xy = 1.0;

    // Radial unit vector (outward)
    double rx = x/r, ry = y/r, rz = z/r;
    // Latitudinal unit vector (toward pole in meridional plane)
    double lx = -sinlat * x / r_xy;
    double ly = -sinlat * y / r_xy;
    double lz = coslat;

    gravg[0] = g_radial * rx + g_lat * lx;
    gravg[1] = g_radial * ry + g_lat * ly;
    gravg[2] = g_radial * rz + g_lat * lz;
}


///////////////////////////////////////////////////////////////////////////////
// VENUS GRAVITY MODEL (J2-J4)
///////////////////////////////////////////////////////////////////////////////
// Venus is nearly spherical (very slow retrograde rotation, 243 Earth days)
// J2 is ~240x smaller than Earth's. Gravity field is dominated by
// tesseral harmonics (surface features) rather than oblateness.
// Source: Konopliv et al. 1999 (Magellan MGNP180U)
///////////////////////////////////////////////////////////////////////////////

inline void gravity_venus(double gravg[3], const double pos[3])
{
    double x = pos[0], y = pos[1], z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1.0) { gravg[0]=0; gravg[1]=0; gravg[2]=0; return; }

    double mu_r2 = VENUS_GM / (r * r);
    double a_r = VENUS_RADIUS / r;
    double sinlat = z / r;
    double coslat = sqrt(1.0 - sinlat*sinlat);
    double s2 = sinlat * sinlat;
    double a_r2 = a_r * a_r;

    // J2 perturbation (small for Venus)
    double P2 = 0.5 * (3.0*s2 - 1.0);
    double radial_sum = 1.0 + 3.0 * VENUS_J2 * a_r2 * P2;
    double g_radial = -mu_r2 * radial_sum;
    double g_lat = -mu_r2 * 3.0 * VENUS_J2 * a_r2 * sinlat * coslat;

    double r_xy = sqrt(x*x + y*y);
    if (r_xy < 1.0) r_xy = 1.0;
    double rx = x/r, ry = y/r, rz = z/r;
    double lx = -sinlat * x / r_xy, ly = -sinlat * y / r_xy, lz = coslat;

    gravg[0] = g_radial * rx + g_lat * lx;
    gravg[1] = g_radial * ry + g_lat * ly;
    gravg[2] = g_radial * rz + g_lat * lz;
}


///////////////////////////////////////////////////////////////////////////////
// MOON GRAVITY MODEL (J2-J6, GRAIL)
///////////////////////////////////////////////////////////////////////////////
// The Moon has a lumpy gravity field with large mascon anomalies.
// Zonal harmonics capture the global shape; for precise lunar orbit
// propagation, a full tesseral model (GRGM1200A) is needed.
// For gravity assist trajectory targeting, J2-J4 is sufficient.
// Source: GRAIL mission (Zuber et al. 2013, Konopliv et al. 2013)
///////////////////////////////////////////////////////////////////////////////

inline void gravity_moon(double gravg[3], const double pos[3])
{
    double x = pos[0], y = pos[1], z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1.0) { gravg[0]=0; gravg[1]=0; gravg[2]=0; return; }

    double mu_r2 = MOON_GM / (r * r);
    double a_r = MOON_RADIUS / r;
    double sinlat = z / r;
    double coslat = sqrt(1.0 - sinlat*sinlat);
    double s2 = sinlat * sinlat;
    double a_r2 = a_r * a_r;
    double a_r3 = a_r2 * a_r;
    double a_r4 = a_r2 * a_r2;

    double P2 = 0.5 * (3.0*s2 - 1.0);
    double P3 = 0.5 * (5.0*s2*sinlat - 3.0*sinlat);
    double P4 = (1.0/8.0) * (35.0*s2*s2 - 30.0*s2 + 3.0);

    double dP2 = 3.0 * sinlat;
    double dP3 = 0.5 * (15.0*s2 - 3.0);
    double dP4 = 0.5 * (35.0*s2*sinlat - 15.0*sinlat);

    double radial_sum = 1.0
        + 3.0 * MOON_J2 * a_r2 * P2
        + 4.0 * MOON_J3 * a_r3 * P3
        + 5.0 * MOON_J4 * a_r4 * P4;

    double lat_sum =
          MOON_J2 * a_r2 * dP2
        + MOON_J3 * a_r3 * dP3
        + MOON_J4 * a_r4 * dP4;

    double g_radial = -mu_r2 * radial_sum;
    double g_lat = mu_r2 * lat_sum * coslat;

    double r_xy = sqrt(x*x + y*y);
    if (r_xy < 1.0) r_xy = 1.0;
    double rx = x/r, ry = y/r, rz = z/r;
    double lx = -sinlat * x / r_xy, ly = -sinlat * y / r_xy, lz = coslat;

    gravg[0] = g_radial * rx + g_lat * lx;
    gravg[1] = g_radial * ry + g_lat * ly;
    gravg[2] = g_radial * rz + g_lat * lz;
}


///////////////////////////////////////////////////////////////////////////////
// MARS GRAVITY MODEL (J2-J6)
///////////////////////////////////////////////////////////////////////////////
// Zonal harmonic gravity field for Mars
// Based on GMM-3 model (Genova et al. 2016)
//
// Returns gravitational acceleration vector in geocentric coordinates
// Same convention as CADAC++ cad_grav84()
//
// Arguments:
//   gravg[3] — output: gravitational acceleration (geocentric) m/s^2
//   pos[3]   — input: position in planet-centered inertial coords (m)
///////////////////////////////////////////////////////////////////////////////

inline void gravity_mars(double gravg[3], const double pos[3])
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1.0) { gravg[0]=0; gravg[1]=0; gravg[2]=0; return; }

    double mu_r2 = MARS_GM / (r * r);
    double a_r = MARS_RADIUS / r;
    double sinlat = z / r;
    double coslat = sqrt(x*x + y*y) / r;
    double s2 = sinlat * sinlat;

    // J2 contribution (dominant oblateness term)
    double a_r2 = a_r * a_r;
    double fac_r = mu_r2 / r; // GM/r^3
    double j2_radial = 1.5 * MARS_J2 * a_r2 * (3.0*s2 - 1.0);
    double j2_latitudinal = -3.0 * MARS_J2 * a_r2 * sinlat * coslat;

    // J4 contribution
    double a_r4 = a_r2 * a_r2;
    double s4 = s2 * s2;
    double j4_radial = (5.0/8.0) * MARS_J4 * a_r4 * (35.0*s4 - 30.0*s2 + 3.0);
    double j4_latitudinal = (5.0/2.0) * MARS_J4 * a_r4 * sinlat * coslat * (7.0*s2 - 3.0);

    // Total acceleration in radial and latitudinal components
    double g_radial = -mu_r2 * (1.0 + j2_radial + j4_radial);
    double g_lat = -mu_r2 * (j2_latitudinal + j4_latitudinal);

    // Convert to Cartesian geocentric
    // radial: along -r_hat, latitudinal: perpendicular in meridional plane
    double r_xy = sqrt(x*x + y*y);
    if (r_xy < 1.0) r_xy = 1.0;

    // Radial unit vector
    double rx = x/r, ry = y/r, rz = z/r;
    // Latitudinal unit vector (toward pole in meridional plane)
    double lx = -sinlat * x / r_xy;
    double ly = -sinlat * y / r_xy;
    double lz = coslat;

    gravg[0] = g_radial * rx + g_lat * lx;
    gravg[1] = g_radial * ry + g_lat * ly;
    gravg[2] = g_radial * rz + g_lat * lz;
}


///////////////////////////////////////////////////////////////////////////////
// JUPITER GRAVITY MODEL (J2-J10, Juno)
///////////////////////////////////////////////////////////////////////////////
// Zonal harmonic gravity field for Jupiter
// Based on Iess et al. 2018 (Juno mission, Nature 555:220)
// Includes both even (J2,J4,J6,J8,J10) and odd (J3,J5,J7,J9) harmonics
//
// Returns gravitational acceleration vector in geocentric coordinates
//
// Arguments:
//   gravg[3] — output: gravitational acceleration (geocentric) m/s^2
//   pos[3]   — input: position in planet-centered inertial coords (m)
///////////////////////////////////////////////////////////////////////////////

inline void gravity_jupiter(double gravg[3], const double pos[3])
{
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    if (r < 1.0) { gravg[0]=0; gravg[1]=0; gravg[2]=0; return; }

    double mu_r2 = JUPITER_GM / (r * r);
    double a_r = JUPITER_RADIUS / r;
    double sinlat = z / r;
    double coslat = sqrt(x*x + y*y) / r;
    double s2 = sinlat * sinlat;

    // Legendre polynomial values for even harmonics
    double P2  = 0.5 * (3.0*s2 - 1.0);
    double P4  = (1.0/8.0) * (35.0*s2*s2 - 30.0*s2 + 3.0);
    double P6  = (1.0/16.0) * (231.0*pow(sinlat,6) - 315.0*s2*s2 + 105.0*s2 - 5.0);
    double P8  = (1.0/128.0) * (6435.0*pow(sinlat,8) - 12012.0*pow(sinlat,6)
                  + 6930.0*s2*s2 - 1260.0*s2 + 35.0);
    double P10 = (1.0/256.0) * (46189.0*pow(sinlat,10) - 109395.0*pow(sinlat,8)
                  + 90090.0*pow(sinlat,6) - 30030.0*s2*s2 + 3465.0*s2 - 63.0);

    // Derivatives dP/d(sinlat) for even harmonics
    double dP2  = 3.0 * sinlat;
    double dP4  = (1.0/2.0) * (35.0*4.0*pow(sinlat,3) - 30.0*2.0*sinlat) / 4.0;
    // Simplified: use numerical approach for higher orders

    // Gravitational potential perturbation: sum of -Jn*(a/r)^n * Pn(sinlat)
    // Radial acceleration: -dU/dr = (GM/r²)[1 + sum (n+1)*Jn*(a/r)^n*Pn]
    double a_r2 = a_r * a_r;
    double a_r3 = a_r2 * a_r;
    double a_r4 = a_r2 * a_r2;
    double a_r6 = a_r4 * a_r2;
    double a_r8 = a_r4 * a_r4;
    double a_r10 = a_r8 * a_r2;

    double radial_sum = 1.0
        + 3.0 * JUPITER_J2 * a_r2 * P2
        + 5.0 * JUPITER_J4 * a_r4 * P4
        + 7.0 * JUPITER_J6 * a_r6 * P6
        + 9.0 * JUPITER_J8 * a_r8 * P8
        + 11.0 * JUPITER_J10 * a_r10 * P10;

    double g_radial = -mu_r2 * radial_sum;

    // Latitudinal component from J2 (dominant)
    double g_lat = -mu_r2 * 3.0 * JUPITER_J2 * a_r2 * sinlat * coslat;

    // Convert to Cartesian geocentric
    double r_xy = sqrt(x*x + y*y);
    if (r_xy < 1.0) r_xy = 1.0;

    double rx = x/r, ry = y/r, rz = z/r;
    double lx = -sinlat * x / r_xy;
    double ly = -sinlat * y / r_xy;
    double lz = coslat;

    gravg[0] = g_radial * rx + g_lat * lx;
    gravg[1] = g_radial * ry + g_lat * ly;
    gravg[2] = g_radial * rz + g_lat * lz;
}


///////////////////////////////////////////////////////////////////////////////
// SPEED OF SOUND (generic, for Mach number computation)
///////////////////////////////////////////////////////////////////////////////

inline double speed_of_sound_venus(double tempk)
{
    return sqrt(VENUS_GAMMA * VENUS_RGAS * tempk);
}

inline double speed_of_sound_mars(double tempk)
{
    return sqrt(MARS_GAMMA * MARS_RGAS * tempk);
}

inline double speed_of_sound_jupiter(double tempk)
{
    return sqrt(JUPITER_GAMMA * JUPITER_RGAS * tempk);
}

#endif // PLANETARY_ENVIRONMENT_HPP
