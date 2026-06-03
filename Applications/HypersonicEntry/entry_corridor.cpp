///////////////////////////////////////////////////////////////////////////////
// ATMOSPHERIC ENTRY CORRIDOR ANALYSIS
//
// Computes viable entry corridors for Earth, Mars, Jupiter, Venus, Titan,
// Neptune, and Uranus using real atmospheric models.
//
// For each body, sweeps entry flight path angles and simulates the
// trajectory to find:
//   - Skip-out boundary (too shallow — vehicle exits atmosphere)
//   - G-limit boundary (too steep — exceeds structural limits)
//   - Heat-rate boundary (exceeds TPS limits)
//
// Physics: 3-DOF ballistic entry with real atmosphere profiles.
// Integration: RK4 with adaptive altitude step.
//
// Compile:
//   g++ -std=c++11 -O2 -o entry_corridor entry_corridor.cpp
//
// 20260603 Created for Kepler Workstation HypersonicEntry application
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <cstdlib>

// Planet IDs for the sweep
enum Planet {
    PL_EARTH, PL_MARS, PL_JUPITER, PL_VENUS,
    PL_TITAN, PL_NEPTUNE, PL_URANUS,
    PL_COUNT
};

static const char* planet_name[] = {
    "Earth", "Mars", "Jupiter", "Venus",
    "Titan", "Neptune", "Uranus"
};

// Planet physical parameters
struct PlanetParams {
    double radius_m;     // reference radius (m)
    double gm;           // gravitational parameter (m³/s²)
    double g0;           // surface gravity (m/s²)
    double rgas;         // specific gas constant (J/kg·K)
};

static const PlanetParams planet_params[PL_COUNT] = {
    // Earth
    { 6378137.0,     3.986004418e14, 9.80665,  287.058 },
    // Mars
    { 3396200.0,     4.282837e13,    3.71,     188.92  },
    // Jupiter
    { 71492000.0,    1.26686534e17,  24.79,    3745.0  },
    // Venus
    { 6051800.0,     3.24859e14,     8.87,     191.0   },
    // Titan
    { 2574700.0,     8.978138e12,    1.352,    296.8   },
    // Neptune
    { 24764000.0,    6.835100e15,    11.15,    3615.0  },
    // Uranus
    { 25559000.0,    5.793951e15,    8.87,     3641.0  },
};


///////////////////////////////////////////////////////////////////////////////
// ATMOSPHERE MODELS (bundled — mirrors planetary_environment.hpp)
///////////////////////////////////////////////////////////////////////////////

static inline void get_atmosphere(Planet pl, double alt_m,
                                   double &rho, double &press, double &tempk)
{
    double h_km = alt_m / 1000.0;

    switch (pl) {
    case PL_EARTH: {
        // US Standard Atmosphere 1976, 7-layer (0-84.852 km)
        if (h_km < 0.0) h_km = 0.0;
        if (h_km > 84.852) h_km = 84.852;
        static const double H[8]  = {0,11,20,32,47,51,71,84.852};
        static const double LR[7] = {-6.5,0,1.0,2.8,0,-2.8,-2.0};
        static const double Tb[7] = {288.15,216.65,216.65,228.65,270.65,270.65,214.65};
        const double R=287.058, g0=9.80665;
        int i=0; for(int k=1;k<7;k++) if(h_km>=H[k]) i=k;
        double lapse=LR[i], T_base=Tb[i];
        tempk = (fabs(lapse)<1e-9) ? T_base : T_base+lapse*(h_km-H[i]);
        double Pb=101325.0;
        for(int k=0;k<i;k++){
            double dh=H[k+1]-H[k];
            if(fabs(LR[k])<1e-9) Pb*=exp(-g0*dh*1000/(R*Tb[k]));
            else { double Tt=Tb[k]+LR[k]*dh; Pb*=pow(Tt/Tb[k],-g0/(R*LR[k]*0.001)); }
        }
        double dh=h_km-H[i];
        if(fabs(lapse)<1e-9) press=Pb*exp(-g0*dh*1000/(R*T_base));
        else { double Tt=T_base+lapse*dh; press=Pb*pow(Tt/T_base,-g0/(R*lapse*0.001)); }
        rho=press/(R*tempk);
        break;
    }
    case PL_MARS: {
        double T_c = (alt_m < 7000.0) ? -31.0-0.000998*alt_m : -23.4-0.00222*alt_m;
        tempk = T_c + 273.1;
        double p_kpa = 0.699*exp(-0.00009*alt_m);
        press = p_kpa*1000.0;
        rho = p_kpa/(0.1921*tempk);
        break;
    }
    case PL_JUPITER: {
        if (h_km > 50.0)       tempk = 110.0 + 0.05*(h_km-50.0);
        else if (h_km > 0.0)   tempk = 166.0 - 1.12*h_km;
        else if (h_km > -60.0) tempk = 166.0 - 2.0*h_km;
        else                    tempk = 286.0 - 2.2*(h_km+60.0);
        if (tempk<100) tempk=100;
        double sh = 3745.0*tempk/24.79;
        press = 1.0e5*exp(-alt_m/sh);
        rho = press/(3745.0*tempk);
        break;
    }
    case PL_VENUS: {
        static const double at[]={0,10,20,30,40,50,60,70,80,100};
        static const double tc[]={462,385,306,222,143,75,-10,-43,-76,-112};
        static const double pa[]={92.1,47.39,22.52,9.851,3.501,1.066,0.2357,0.0369,0.00476,0.0000266};
        if(h_km<0)h_km=0; if(h_km>100)h_km=100;
        int i=0; for(int k=0;k<9;k++) if(h_km>=at[k]) i=k;
        if(i>=9)i=8;
        double frac=(h_km-at[i])/(at[i+1]-at[i]);
        tempk=(tc[i]+frac*(tc[i+1]-tc[i]))+273.15;
        double lp=log(pa[i])+frac*(log(pa[i+1])-log(pa[i]));
        press=exp(lp)*101325.0;
        rho=press/(191.0*tempk);
        break;
    }
    case PL_TITAN: {
        if (h_km>250)           tempk=170.0+0.1*(h_km-250);
        else if(h_km>40)        tempk=70.0+(170-70)/(250-40)*(h_km-40);
        else if(h_km>0)         tempk=94.0-0.6*h_km;
        else                    tempk=94.0;
        double sh=296.8*tempk/1.352;
        press=1.47e5*exp(-alt_m/sh);
        rho=press/(296.8*tempk);
        break;
    }
    case PL_NEPTUNE: {
        if (h_km>50)            tempk=55.0+0.4*(h_km-50);
        else if(h_km>0)         tempk=72.0-0.34*h_km;
        else                    tempk=72.0-0.9*h_km;
        double sh=3615.0*tempk/11.15;
        press=1.0e5*exp(-alt_m/sh);
        rho=press/(3615.0*tempk);
        break;
    }
    case PL_URANUS: {
        if (h_km>50)            tempk=53.0+0.3*(h_km-50);
        else if(h_km>0)         tempk=76.0-0.46*h_km;
        else                    tempk=76.0-0.8*h_km;
        double sh=3641.0*tempk/8.87;
        press=1.0e5*exp(-alt_m/sh);
        rho=press/(3641.0*tempk);
        break;
    }
    default: rho=0; press=0; tempk=0;
    }
    if(rho<0)rho=0; if(press<0)press=0;
}


///////////////////////////////////////////////////////////////////////////////
// VEHICLE MODEL
///////////////////////////////////////////////////////////////////////////////

struct Vehicle {
    double mass_kg;        // vehicle mass
    double area_m2;        // reference area (drag)
    double cd;             // drag coefficient
    double rn_m;           // nose radius (for heating)
    double g_limit;        // max g-load
    double q_limit_kw_m2;  // max stagnation heat rate (kW/m²)
};

// Default: blunt-body entry capsule (Apollo/Orion class)
static const Vehicle DEFAULT_VEHICLE = {
    5500.0,    // kg
    12.57,     // m² (~4m diameter)
    1.3,       // CD (blunt body)
    2.0,       // m nose radius
    10.0,      // max g
    1500.0     // kW/m² (carbon phenolic TPS)
};


///////////////////////////////////////////////////////////////////////////////
// TRAJECTORY SIMULATION (3-DOF ballistic entry, RK4)
///////////////////////////////////////////////////////////////////////////////

struct EntryState {
    double alt;    // altitude above reference (m)
    double vel;    // velocity (m/s)
    double fpa;    // flight path angle (rad, negative = descending)
    double range;  // downrange (m)
};

struct EntryResult {
    bool   skip_out;     // vehicle exited atmosphere
    bool   impact;       // vehicle hit surface
    bool   g_exceed;     // exceeded g-limit
    bool   q_exceed;     // exceeded heat-rate limit
    double alt_min;      // minimum altitude reached
    double g_max;        // peak deceleration (g)
    double q_max_kw_m2;  // peak heat rate (kW/m²)
    double vel_impact;   // velocity at impact (m/s)
};

// Sutton-Graves stagnation heat flux (kW/m²)
static inline double heat_rate_kw_m2(double rho, double vel, double rn_m)
{
    // q = 1.83e-4 * sqrt(rho / rn) * vel^3  (kW/m²)
    // Sutton & Graves, NASA TR R-382 (1971)
    return 1.83e-4 * sqrt(rho / rn_m) * vel * vel * vel;
}

EntryResult simulate_entry(Planet pl, const Vehicle& v,
                            double alt_entry, double vel_entry,
                            double fpa_entry_deg)
{
    EntryResult res;
    res.skip_out  = false;
    res.impact    = false;
    res.g_exceed  = false;
    res.q_exceed  = false;
    res.alt_min   = alt_entry;
    res.g_max     = 0.0;
    res.q_max_kw_m2 = 0.0;
    res.vel_impact = 0.0;

    const PlanetParams& pp = planet_params[pl];
    double beta = v.mass_kg / (v.cd * v.area_m2);  // ballistic coefficient

    EntryState s;
    s.alt   = alt_entry;
    s.vel   = vel_entry;
    s.fpa   = fpa_entry_deg * M_PI / 180.0;  // rad
    s.range = 0.0;

    double dt = 0.1;  // time step
    const double DT_MIN = 0.001;
    const double DT_MAX = 5.0;
    const int    MAX_STEPS = 500000;
    const double ATMOS_CEILING = (pl == PL_JUPITER || pl == PL_NEPTUNE || pl == PL_URANUS)
                                  ? 500000.0 : 200000.0;  // m

    for (int step = 0; step < MAX_STEPS; step++) {
        if (s.alt > ATMOS_CEILING && s.fpa > 0.0) {
            res.skip_out = true;
            break;
        }
        if (s.alt <= 0.0) {
            res.impact = true;
            res.vel_impact = s.vel;
            break;
        }
        if (s.vel < 100.0) {
            res.impact = true;
            res.vel_impact = s.vel;
            break;
        }

        // Atmosphere at current altitude
        double rho, press, tempk;
        get_atmosphere(pl, s.alt, rho, press, tempk);

        // Gravity
        double r = pp.radius_m + s.alt;
        double g = pp.gm / (r * r);

        // Drag deceleration
        double drag = 0.5 * rho * s.vel * s.vel / beta;
        double g_load = drag / 9.80665;

        // Heat rate
        double q_kw = (rho > 1e-15) ? heat_rate_kw_m2(rho, s.vel, v.rn_m) : 0.0;

        // Check limits
        if (g_load > res.g_max)        res.g_max = g_load;
        if (q_kw  > res.q_max_kw_m2)   res.q_max_kw_m2 = q_kw;
        if (s.alt < res.alt_min)        res.alt_min = s.alt;
        if (g_load > v.g_limit)         { res.g_exceed = true; break; }
        if (q_kw   > v.q_limit_kw_m2)   { res.q_exceed = true; break; }

        // Adaptive timestep based on density
        if (rho > 1e-6) {
            dt = 0.1 / sqrt(rho * 1000.0);  // smaller steps in denser air
            if (dt < DT_MIN) dt = DT_MIN;
            if (dt > DT_MAX) dt = DT_MAX;
        } else {
            dt = DT_MAX;
        }

        // RK4 integration
        // Derivatives: dh/dt = v*sin(γ), dv/dt = -drag - g*sin(γ), dγ/dt = (v/r - g/v)*cos(γ)
        auto derivs = [&](const EntryState& st, EntryState& ds, double& h_dot) {
            (void)h_dot;
            double rho_l, press_l, tempk_l;
            get_atmosphere(pl, st.alt, rho_l, press_l, tempk_l);
            double r_l = pp.radius_m + st.alt;
            double g_l = pp.gm / (r_l * r_l);
            double drag_l = 0.5 * rho_l * st.vel * st.vel / beta;
            ds.vel   = -drag_l - g_l * sin(st.fpa);
            ds.fpa   = (st.vel / r_l - g_l / st.vel) * cos(st.fpa);
            ds.alt   = st.vel * sin(st.fpa);
            ds.range = st.vel * cos(st.fpa);
        };

        // RK4 step
        EntryState k1, k2, k3, k4, stmp;
        double dummy;
        derivs(s, k1, dummy);
        stmp.alt = s.alt + 0.5*dt*k1.alt;
        stmp.vel = s.vel + 0.5*dt*k1.vel;
        stmp.fpa = s.fpa + 0.5*dt*k1.fpa;
        stmp.range = s.range + 0.5*dt*k1.range;
        derivs(stmp, k2, dummy);
        stmp.alt = s.alt + 0.5*dt*k2.alt;
        stmp.vel = s.vel + 0.5*dt*k2.vel;
        stmp.fpa = s.fpa + 0.5*dt*k2.fpa;
        stmp.range = s.range + 0.5*dt*k2.range;
        derivs(stmp, k3, dummy);
        stmp.alt = s.alt + dt*k3.alt;
        stmp.vel = s.vel + dt*k3.vel;
        stmp.fpa = s.fpa + dt*k3.fpa;
        stmp.range = s.range + dt*k3.range;
        derivs(stmp, k4, dummy);

        s.alt   += dt/6.0 * (k1.alt   + 2*k2.alt   + 2*k3.alt   + k4.alt);
        s.vel   += dt/6.0 * (k1.vel   + 2*k2.vel   + 2*k3.vel   + k4.vel);
        s.fpa   += dt/6.0 * (k1.fpa   + 2*k2.fpa   + 2*k3.fpa   + k4.fpa);
        s.range += dt/6.0 * (k1.range + 2*k2.range + 2*k3.range + k4.range);
    }

    return res;
}


///////////////////////////////////////////////////////////////////////////////
// ENTRY CORRIDOR SWEEP
///////////////////////////////////////////////////////////////////////////////

struct CorridorResult {
    double skip_out_deg;   // shallowest angle that doesn't skip out
    double g_limit_deg;    // steepest angle within g-limit
    double q_limit_deg;    // steepest angle within heat-rate limit
    double corridor_width; // viable entry corridor (deg)
    double g_peak;         // peak g at corridor center
    double q_peak_kw_m2;   // peak heat rate at corridor center
};

CorridorResult find_corridor(Planet pl, const Vehicle& v,
                              double alt_entry, double vel_entry)
{
    CorridorResult cr;
    cr.skip_out_deg   = -99;
    cr.g_limit_deg    = -99;
    cr.q_limit_deg    = -99;
    cr.corridor_width = 0;
    cr.g_peak         = 0;
    cr.q_peak_kw_m2   = 0;

    // Sweep from very shallow (-1°) to very steep (-60°)
    // Step: 0.1° for shallow, 0.5° for steep
    double best_fpa = -10.0;
    double best_margin = -1e10;

    for (double fpa = -0.5; fpa >= -60.0; ) {
        EntryResult r = simulate_entry(pl, v, alt_entry, vel_entry, fpa);

        // First non-skip-out angle
        if (cr.skip_out_deg < -90 && !r.skip_out && !r.g_exceed && !r.q_exceed) {
            cr.skip_out_deg = fpa;
        }

        // Last angle before g-limit exceeded
        if (!r.g_exceed) {
            cr.g_limit_deg = fpa;
            if (r.g_max > cr.g_peak) cr.g_peak = r.g_max;
        }

        // Last angle before heat-rate limit exceeded
        if (!r.q_exceed) {
            cr.q_limit_deg = fpa;
            if (r.q_max_kw_m2 > cr.q_peak_kw_m2) cr.q_peak_kw_m2 = r.q_max_kw_m2;
        }

        // Best trajectory (most margin to both limits)
        double margin = (v.g_limit - r.g_max) + 0.01*(v.q_limit_kw_m2 - r.q_max_kw_m2);
        if (!r.skip_out && !r.g_exceed && !r.q_exceed && margin > best_margin) {
            best_margin = margin;
            best_fpa = fpa;
        }

        // Adaptive step
        if (fpa > -3.0)       fpa -= 0.1;
        else if (fpa > -10.0)  fpa -= 0.25;
        else if (fpa > -20.0)  fpa -= 0.5;
        else                    fpa -= 1.0;
    }

    // Corridor width (most restrictive bounds)
    double upper = cr.skip_out_deg;
    double lower = (cr.g_limit_deg < cr.q_limit_deg) ? cr.g_limit_deg : cr.q_limit_deg;
    if (upper < -90 || lower < -90) {
        cr.corridor_width = 0;
    } else {
        cr.corridor_width = fabs(lower - upper);
    }

    return cr;
}


///////////////////////////////////////////////////////////////////////////////
// ENTRY VELOCITIES BY BODY
///////////////////////////////////////////////////////////////////////////////

struct EntryConfig {
    Planet     planet;
    double     alt_entry_m;    // entry interface altitude
    double     vel_entry_kms;  // entry velocity (km/s)
};

static const EntryConfig entry_configs[] = {
    { PL_EARTH,   120000, 11.0  },  // LEO return / lunar return
    { PL_MARS,    125000, 6.0   },  // interplanetary arrival
    { PL_JUPITER, 500000, 48.0  },  // Galileo-class entry probe
    { PL_VENUS,   150000, 11.5  },  // interplanetary arrival
    { PL_TITAN,   500000, 6.5   },  // Saturn system arrival
    { PL_NEPTUNE, 500000, 23.0  },  // interplanetary arrival (~similar to Jupiter but less)
    { PL_URANUS,  500000, 22.0  },  // interplanetary arrival
};


///////////////////////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////////////////////

int main()
{
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║  KEPLER WORKSTATION — Atmospheric Entry Corridors    ║\n");
    printf("║  7 bodies · 3-DOF ballistic · real atmospheres       ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");
    printf("\n");

    printf("  Vehicle: %.0f kg blunt-body capsule\n", DEFAULT_VEHICLE.mass_kg);
    printf("    Ballistic coefficient: %.1f kg/m²\n",
           DEFAULT_VEHICLE.mass_kg / (DEFAULT_VEHICLE.cd * DEFAULT_VEHICLE.area_m2));
    printf("    Limits: %.0f g, %.0f kW/m² peak heat rate\n",
           DEFAULT_VEHICLE.g_limit, DEFAULT_VEHICLE.q_limit_kw_m2);
    printf("\n");

    printf("══════════════════════════════════════════════════════════════════════\n");
    printf("  %-12s %6s %8s %10s %10s %8s %8s\n",
           "Body", "v_inf", "Skip-out", "G-limit", "Q-limit", "Width", "Peak g");
    printf("  %-12s %6s %8s %10s %10s %8s %8s\n",
           "------------", "------", "--------", "----------", "----------", "--------", "--------");
    printf("══════════════════════════════════════════════════════════════════════\n");

    for (int i = 0; i < 7; i++) {
        const EntryConfig& ec = entry_configs[i];

        printf("  Simulating %s entry (%.1f km/s, %.0f km interface)...\n",
               planet_name[ec.planet], ec.vel_entry_kms, ec.alt_entry_m/1000.0);

        CorridorResult cr = find_corridor(ec.planet, DEFAULT_VEHICLE,
                                           ec.alt_entry_m, ec.vel_entry_kms * 1000.0);

        printf("  %-12s %5.1f  %7.1f° %9.1f° %9.1f° %7.1f° %7.1f\n",
               planet_name[ec.planet],
               ec.vel_entry_kms,
               cr.skip_out_deg,
               cr.g_limit_deg,
               cr.q_limit_deg,
               cr.corridor_width,
               cr.g_peak);

        if (cr.corridor_width > 0) {
            printf("    Corridor: %.1f° to %.1f° (total %.1f°)  |  "
                   "Peak: %.1f g, %.0f kW/m²\n",
                   cr.skip_out_deg,
                   (cr.g_limit_deg < cr.q_limit_deg) ? cr.g_limit_deg : cr.q_limit_deg,
                   cr.corridor_width,
                   cr.g_peak, cr.q_peak_kw_m2);
        } else {
            printf("    No viable corridor found with current vehicle limits.\n");
        }
        printf("\n");
    }

    printf("══════════════════════════════════════════════════════════════════════\n");
    printf("  Interpretation:\n");
    printf("    Skip-out:  shallowest angle that captures (more negative = easier)\n");
    printf("    G-limit:   steepest angle within structural g-limit\n");
    printf("    Q-limit:   steepest angle within TPS heat-rate limit\n");
    printf("    Width:     viable entry corridor (total degrees)\n");
    printf("    Wider corridor = easier targeting = more margin\n");
    printf("══════════════════════════════════════════════════════════════════════\n");
    printf("\n");

    printf("  References:\n");
    printf("    Allen & Eggers (1958) — analytical entry mechanics\n");
    printf("    Sutton & Graves (1971) — stagnation heat flux correlation\n");
    printf("    US Standard Atmosphere 1976 — Earth\n");
    printf("    Galileo Probe ASI (Seiff et al. 1998) — Jupiter\n");
    printf("    Voyager 2 RSS (Lindal 1987, 1990) — Uranus, Neptune\n");
    printf("    Cassini-Huygens (Niemann et al. 2005) — Titan\n");
    printf("    Venera/Pioneer Venus — Venus\n");
    printf("    Mars Global Surveyor — Mars\n");
    printf("\n");

    return 0;
}
