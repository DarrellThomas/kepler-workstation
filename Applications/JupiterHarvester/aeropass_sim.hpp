// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// TRAJECTORY-COUPLED ATMOSPHERIC PASS SIMULATION
//
// Integrates real-gas thermodynamics along the atmospheric pass arc:
//   - Per-timestep: atmosphere query → real-gas shock → heating → drag → orbit update
//   - Multi-pass campaign with aerobraking (Mach reduction) or collection modes
//   - Tracks cumulative heat load, mass collected, orbit decay
//
// Couples stat_thermo.hpp (7-species H2/He equilibrium) with orbital mechanics
// to answer: what does a multi-pass aerobraking + collection campaign actually
// look like with real physics?
//
// Key insight: at Mach 100, perfect gas says density ratio = 5.8.
// Real gas says 10+. That changes EVERYTHING — shock standoff, drag, heating.
//
// References:
//   Sutton-Graves stagnation heating: Sutton & Graves, NASA TR R-195 (1971)
//   Fay-Riddell (equilibrium): Fay & Riddell, J. Aero. Sci. 25:73 (1958)
//   Drag in continuum regime: Cd ≈ 2.0 for bluff body behind mag. shield
///////////////////////////////////////////////////////////////////////////////

#ifndef AEROPASS_SIM_HPP
#define AEROPASS_SIM_HPP

#include "../../PlanetaryModels/stat_thermo.hpp"
#include "../../PlanetaryModels/planetary_environment.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

namespace aeropass {

using namespace stat_thermo;

///////////////////////////////////////////////////////////////////////////////
// CONSTANTS
///////////////////////////////////////////////////////////////////////////////

static const double MU0     = 4.0 * M_PI * 1e-7;
static const double GM_J    = 1.26686534e17;    // m³/s² (Jupiter)
static const double R_J     = 71492000.0;       // m (Jupiter equatorial radius)
static const double SIGMA_SB = 5.670374419e-8;  // Stefan-Boltzmann (W/(m²·K⁴))

///////////////////////////////////////////////////////////////////////////////
// VEHICLE CONFIGURATION
///////////////////////////////////////////////////////////////////////////////

struct Vehicle {
    double dry_mass;        // kg
    double tank_capacity;   // kg (max H2+He storage)
    double B_coil;          // Tesla (superconducting solenoid)
    double r_coil;          // m (coil radius)
    double eta_collect;     // collection efficiency (0-1)
    double Cd;              // drag coefficient (bluff body behind mag shield)
    double A_ref;           // reference area for drag (m²) — physical cross-section
    double emissivity;      // hull surface emissivity
    double max_heat_flux;   // kW/m² — hull material limit
    double leak_fraction;   // fraction of particles leaking through mag field
};

inline Vehicle default_bucket()
{
    Vehicle v;
    v.dry_mass      = 8000;     // kg
    v.tank_capacity = 33000;    // kg
    v.B_coil        = 20.0;    // Tesla
    v.r_coil        = 3.0;     // m
    v.eta_collect   = 0.5;
    v.Cd            = 2.0;     // bluff body
    v.A_ref         = M_PI * 3.0 * 3.0; // physical cross-section ~28 m²
    v.emissivity    = 0.9;
    v.max_heat_flux = 50.0;    // kW/m² (SiC limit, no TPS needed below this)
    v.leak_fraction = 0.01;    // 1% leaks through magnetic field
    return v;
}


///////////////////////////////////////////////////////////////////////////////
// ORBIT STATE
///////////////////////////////////////////////////////////////////////////////

struct Orbit {
    double a;       // semi-major axis (m)
    double e;       // eccentricity
    double r_peri;  // periapsis radius (m)
    double r_apo;   // apoapsis radius (m)
    double v_peri;  // velocity at periapsis (m/s)
    double T_period;// orbital period (s)
    double energy;  // specific orbital energy (J/kg)

    // Altitude above 1-bar level at periapsis (m)
    double alt_peri() const { return r_peri - R_J; }
    double alt_peri_km() const { return (r_peri - R_J) / 1000.0; }
};

inline Orbit orbit_from_periapsis(double alt_peri_m, double r_apo_m)
{
    Orbit o;
    o.r_peri = R_J + alt_peri_m;
    o.r_apo = r_apo_m;
    o.a = (o.r_peri + o.r_apo) / 2.0;
    o.e = (o.r_apo - o.r_peri) / (o.r_apo + o.r_peri);
    o.v_peri = sqrt(GM_J * (2.0 / o.r_peri - 1.0 / o.a));
    o.T_period = 2.0 * M_PI * sqrt(o.a * o.a * o.a / GM_J);
    o.energy = -GM_J / (2.0 * o.a);
    return o;
}

// Update orbit after energy change (from drag)
inline Orbit orbit_after_energy_loss(const Orbit &old, double delta_E_per_kg)
{
    // delta_E_per_kg is negative (energy lost to drag)
    Orbit o;
    o.energy = old.energy + delta_E_per_kg;
    o.a = -GM_J / (2.0 * o.energy);

    // Keep periapsis approximately the same (drag at periapsis lowers apoapsis)
    o.r_peri = old.r_peri;  // first approximation
    o.r_apo = 2.0 * o.a - o.r_peri;

    // If apoapsis went below periapsis, we've circularized
    if (o.r_apo < o.r_peri) {
        o.r_apo = o.r_peri;
        o.a = o.r_peri;
        o.e = 0.0;
    } else {
        o.e = (o.r_apo - o.r_peri) / (o.r_apo + o.r_peri);
    }

    o.v_peri = sqrt(GM_J * (2.0 / o.r_peri - 1.0 / o.a));
    o.T_period = 2.0 * M_PI * sqrt(o.a * o.a * o.a / GM_J);
    return o;
}


///////////////////////////////////////////////////////////////////////////////
// MAGNETIC FUNNEL PHYSICS
///////////////////////////////////////////////////////////////////////////////

// Magnetic funnel radius: dipole field B(r) = B_coil * (r_coil/r)³
// Set equal to ram pressure: B²/(2μ₀) = ½ρv²
inline double mag_funnel_radius(double B_coil, double r_coil,
                                 double rho, double v, double max_r = 500.0)
{
    double p_ram = 0.5 * rho * v * v;
    if (p_ram <= 0) return max_r;
    double B_need = sqrt(2.0 * MU0 * p_ram);
    double R = r_coil * pow(B_coil / B_need, 1.0 / 3.0);
    return (R < max_r) ? R : max_r;
}


///////////////////////////////////////////////////////////////////////////////
// PER-PASS RESULTS
///////////////////////////////////////////////////////////////////////////////

struct PassResult {
    int pass_number;
    Orbit orbit_pre;            // orbit before this pass
    Orbit orbit_post;           // orbit after drag
    double alt_peri_km;         // periapsis altitude (km)
    double v_peri;              // periapsis velocity (m/s)
    double Mach_peri;           // Mach number at periapsis
    double pass_duration;       // seconds in atmosphere
    double arc_length;          // meters of atmosphere traversed

    // Thermodynamics
    double T_shock_peak;        // peak post-shock temperature (K)
    double gamma_shock_peak;    // γ at peak shock conditions
    double density_ratio_peak;  // ρ₂/ρ₁ at peak

    // Heating
    double q_stag_peak;         // peak stagnation-point heat flux (kW/m²)
    double q_hull_peak;         // peak hull heat flux (kW/m², after mag shield)
    double Q_total;             // integrated heat load over pass (MJ)

    // Collection
    double R_mag_avg;           // average magnetic funnel radius (m)
    double mass_collected;      // kg collected this pass
    double h2_collected;        // kg H2
    double he_collected;        // kg He

    // Drag
    double delta_v_drag;        // velocity change from drag (m/s)
    double delta_E_drag;        // specific energy change from drag (J/kg)

    // Cumulative
    double total_collected;     // cumulative kg in tank
    double tank_fraction;       // fraction of tank full
};


///////////////////////////////////////////////////////////////////////////////
// SINGLE-PASS ATMOSPHERIC INTEGRATION
//
// Integrates along the pass arc with real-gas thermodynamics.
// The trajectory is approximated as a straight line through a
// spherically-symmetric atmosphere at the periapsis altitude.
///////////////////////////////////////////////////////////////////////////////

inline PassResult integrate_pass(int pass_num, const Orbit &orb,
                                  const Vehicle &veh,
                                  double Y_H2 = 0.862, double Y_He = 0.138)
{
    PassResult res;
    res.pass_number = pass_num;
    res.orbit_pre = orb;
    res.alt_peri_km = orb.alt_peri_km();
    res.v_peri = orb.v_peri;

    // Atmosphere at periapsis
    double rho_peri, p_peri, T_peri;
    atmosphere_jupiter(rho_peri, p_peri, T_peri, orb.alt_peri());

    // Speed of sound (freestream, perfect gas is fine here — cold gas)
    double a_inf = sqrt(JUPITER_GAMMA * JUPITER_RGAS * T_peri);
    res.Mach_peri = orb.v_peri / a_inf;

    // Scale height at periapsis
    double H = JUPITER_RGAS * T_peri / JUPITER_G0;

    // Arc length: chord through atmosphere ≈ 2√(2Rh)
    res.arc_length = 2.0 * sqrt(2.0 * orb.r_peri * H);
    res.pass_duration = res.arc_length / orb.v_peri;

    // Integration along the arc
    // Parameterize by distance s from periapsis, s ∈ [-arc/2, +arc/2]
    // Altitude varies as: h(s) ≈ h_peri + s²/(2R)
    int N_steps = 200;
    double ds = res.arc_length / N_steps;
    double s_start = -res.arc_length / 2.0;

    // Accumulators
    double Q_total = 0.0;              // integrated heat (J)
    double drag_impulse = 0.0;         // integrated drag force·dt (N·s)
    double mass_collected = 0.0;       // kg collected
    double q_stag_peak = 0.0;
    double q_hull_peak = 0.0;
    double T_shock_peak = 0.0;
    double gamma_shock_peak = 0.0;
    double density_ratio_peak = 0.0;
    double R_mag_sum = 0.0;

    double v_local = orb.v_peri;  // velocity (updated with drag)
    double vehicle_mass = veh.dry_mass + mass_collected;  // updated each step

    for (int i = 0; i < N_steps; i++) {
        double s = s_start + (i + 0.5) * ds;
        double dt = ds / v_local;

        // Local altitude: h(s) = h_peri + s²/(2R_peri)
        double alt_local = orb.alt_peri() + s * s / (2.0 * orb.r_peri);
        if (alt_local < 0) alt_local = 0;

        // Local atmosphere
        double rho, press, tempk;
        atmosphere_jupiter(rho, press, tempk, alt_local);

        if (rho < 1e-15) continue;  // above atmosphere

        // ── Real-gas shock ──
        // Quick Rankine-Hugoniot with real γ
        GasMixture freestream = equilibrium_state(tempk, press, Y_H2, Y_He);
        double M_local = v_local / freestream.a_sound;

        // Post-shock temperature estimate (iterative with real gas)
        // Use energy conservation: h₂ + v₂²/2 = h₁ + v₁²/2
        double h_stag = freestream.h + 0.5 * v_local * v_local;

        // Estimate post-shock T by bisection on h(T) = h_stag
        double T_lo = tempk, T_hi = 300000.0;
        double T_shock = 0.5 * (T_lo + T_hi);
        for (int bis = 0; bis < 40; bis++) {
            GasMixture trial = equilibrium_state(T_shock, press * 10.0, Y_H2, Y_He);
            if (trial.h < h_stag)
                T_lo = T_shock;
            else
                T_hi = T_shock;
            T_shock = 0.5 * (T_lo + T_hi);
        }

        GasMixture post_shock = equilibrium_state(T_shock, press * 10.0, Y_H2, Y_He);
        double rho_ratio = post_shock.rho / rho;

        if (T_shock > T_shock_peak) {
            T_shock_peak = T_shock;
            gamma_shock_peak = post_shock.gamma_eff;
            density_ratio_peak = rho_ratio;
        }

        // ── Stagnation-point heating ──
        // Sutton-Graves with real-gas correction
        // q_sg = k * sqrt(rho/R_nose) * v³
        // Real-gas correction factor: (h_stag - h_wall) / (h_stag_perfect)
        // Approximate: multiply by sqrt(rho_ratio / rho_ratio_perfect)
        double R_nose = 3.0;  // effective nose radius (m)
        double k_sg = 1.7415e-4;  // Sutton-Graves constant for H2/He
        double q_stag = k_sg * sqrt(rho / R_nose) * pow(v_local, 3.0) / 1e3;  // kW/m²

        // Real-gas correction: density ratio increases → thinner shock → more heating
        double perfect_rho_ratio = (JUPITER_GAMMA + 1.0) * M_local * M_local /
                                   ((JUPITER_GAMMA - 1.0) * M_local * M_local + 2.0);
        double rg_correction = sqrt(rho_ratio / perfect_rho_ratio);
        q_stag *= rg_correction;

        // Hull heating (behind magnetic shield)
        double q_hull = q_stag * veh.leak_fraction;

        if (q_stag > q_stag_peak) q_stag_peak = q_stag;
        if (q_hull > q_hull_peak) q_hull_peak = q_hull;

        // Accumulate heat load (on hull)
        Q_total += q_hull * 1e3 * veh.A_ref * dt;  // J

        // ── Drag force ──
        // Effective drag area includes magnetic funnel interaction
        double R_mag = mag_funnel_radius(veh.B_coil, veh.r_coil, rho, v_local);
        R_mag_sum += R_mag;
        double A_drag = M_PI * R_mag * R_mag * 0.1;  // drag area ≈ 10% of funnel area
        // Plus physical body drag
        A_drag += veh.A_ref * veh.Cd;

        double F_drag = 0.5 * rho * v_local * v_local * A_drag;
        drag_impulse += F_drag * dt;

        // Velocity loss from drag (this step)
        vehicle_mass = veh.dry_mass + mass_collected;
        double dv = F_drag * dt / vehicle_mass;
        v_local -= dv;
        if (v_local < 100.0) v_local = 100.0;  // minimum — we're not stopping

        // ── Mass collection ──
        double A_collect = M_PI * R_mag * R_mag;
        double m_dot = rho * v_local * A_collect * veh.eta_collect;
        mass_collected += m_dot * dt;
    }

    // Results
    res.T_shock_peak = T_shock_peak;
    res.gamma_shock_peak = gamma_shock_peak;
    res.density_ratio_peak = density_ratio_peak;
    res.q_stag_peak = q_stag_peak;
    res.q_hull_peak = q_hull_peak;
    res.Q_total = Q_total / 1e6;  // MJ
    res.R_mag_avg = R_mag_sum / N_steps;
    res.mass_collected = mass_collected;
    res.h2_collected = mass_collected * 0.75;
    res.he_collected = mass_collected * 0.25;

    // Drag effects on orbit
    res.delta_v_drag = drag_impulse / (veh.dry_mass + mass_collected);
    res.delta_E_drag = -res.delta_v_drag * orb.v_peri;  // approximate

    // Update orbit
    res.orbit_post = orbit_after_energy_loss(orb, res.delta_E_drag);

    return res;
}


///////////////////////////////////////////////////////////////////////////////
// MULTI-PASS CAMPAIGN
///////////////////////////////////////////////////////////////////////////////

enum CampaignMode {
    AEROBRAKE,      // intentional drag to reduce velocity
    COLLECT,        // collection mode at target velocity
    COMPLETE        // tank full or campaign ended
};

struct CampaignResult {
    std::vector<PassResult> passes;
    int total_passes;
    double total_collected;     // kg
    double total_h2;            // kg
    double total_he;            // kg
    double total_heat_load;     // MJ
    double max_heat_flux;       // kW/m² (hull, peak across all passes)
    double v_final;             // final periapsis velocity (m/s)
    double Mach_final;
    Orbit final_orbit;
    double campaign_days;
};

inline CampaignResult run_campaign(
    double alt_peri_init_m,     // initial periapsis altitude (m above 1-bar)
    double r_apo_init_m,        // initial apoapsis radius (m)
    const Vehicle &veh,
    double v_target = 20000.0,  // m/s — switch from aerobrake to collect
    int max_passes = 2000,
    double max_campaign_days = 730.0,  // 2 years default
    double Y_H2 = 0.862,
    double Y_He = 0.138)
{
    CampaignResult campaign;
    campaign.total_collected = 0;
    campaign.total_h2 = 0;
    campaign.total_he = 0;
    campaign.total_heat_load = 0;
    campaign.max_heat_flux = 0;
    campaign.campaign_days = 0;

    Orbit orb = orbit_from_periapsis(alt_peri_init_m, r_apo_init_m);

    printf("\n╔══════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║  MULTI-PASS ATMOSPHERIC CAMPAIGN — REAL GAS THERMODYNAMICS                 ║\n");
    printf("╚══════════════════════════════════════════════════════════════════════════════╝\n\n");
    printf("  Initial orbit: %.0f km × %.1f Rj, v_peri = %.1f km/s (Mach %.0f)\n",
           orb.alt_peri() / 1000, orb.r_apo / R_J, orb.v_peri / 1000,
           orb.v_peri / sqrt(JUPITER_GAMMA * JUPITER_RGAS * 166.0));
    printf("  Target velocity for collection: %.0f m/s (%.1f km/s)\n", v_target, v_target/1000);
    printf("  Vehicle: %.0f kg dry, %.0f kg tank, B=%.0f T\n\n", veh.dry_mass, veh.tank_capacity, veh.B_coil);

    printf("  %-5s %-6s %-8s %-8s %-8s %-8s %-8s %-10s %-10s %-8s %-8s\n",
           "Pass", "Mode", "Alt(km)", "v(km/s)", "Mach", "T_shk(K)", "γ_shk",
           "q_hull", "kg/pass", "Tank(t)", "Apo(Rj)");
    printf("  ──────────────────────────────────────────────────────────────"
           "──────────────────────────────────────────\n");

    CampaignMode mode = AEROBRAKE;

    for (int pass = 1; pass <= max_passes; pass++) {
        // Check time limit
        campaign.campaign_days += orb.T_period / 86400.0;
        if (campaign.campaign_days > max_campaign_days) {
            mode = COMPLETE;
            break;
        }

        // Determine mode
        if (campaign.total_collected >= veh.tank_capacity) {
            mode = COMPLETE;
            break;
        }

        if (orb.v_peri <= v_target && mode == AEROBRAKE) {
            mode = COLLECT;
            printf("\n  *** SWITCHED TO COLLECTION MODE at pass %d "
                   "(v=%.1f km/s, Mach %.0f) ***\n\n", pass, orb.v_peri/1000,
                   orb.v_peri / sqrt(JUPITER_GAMMA * JUPITER_RGAS * 166.0));
        }

        // Integrate this pass
        PassResult res = integrate_pass(pass, orb, veh, Y_H2, Y_He);

        // Accumulate
        campaign.total_collected += res.mass_collected;
        campaign.total_h2 += res.h2_collected;
        campaign.total_he += res.he_collected;
        campaign.total_heat_load += res.Q_total;
        if (res.q_hull_peak > campaign.max_heat_flux)
            campaign.max_heat_flux = res.q_hull_peak;

        res.total_collected = campaign.total_collected;
        res.tank_fraction = campaign.total_collected / veh.tank_capacity;

        campaign.passes.push_back(res);

        // Print milestones
        bool milestone = (pass <= 5) || (pass <= 20 && pass % 5 == 0) ||
                         (pass <= 100 && pass % 20 == 0) ||
                         (pass % 100 == 0) || (mode == COLLECT && pass % 50 == 0) ||
                         (campaign.total_collected >= veh.tank_capacity);

        if (milestone) {
            const char *mode_str = (mode == AEROBRAKE) ? "BRAKE" :
                                   (mode == COLLECT)   ? "COLL"  : "DONE";
            printf("  %-5d %-6s %-8.0f %-8.1f %-8.0f %-8.0f %-8.3f %-10.2f %-10.1f %-8.1f %-8.1f\n",
                   pass, mode_str,
                   res.alt_peri_km, res.v_peri/1000, res.Mach_peri,
                   res.T_shock_peak, res.gamma_shock_peak,
                   res.q_hull_peak, res.mass_collected,
                   campaign.total_collected/1000, res.orbit_post.r_apo / R_J);
        }

        // Update orbit for next pass
        orb = res.orbit_post;

        // Safety: if orbit decays too much, stop
        if (orb.alt_peri() < 50000.0) {  // 50 km minimum
            printf("\n  *** WARNING: Periapsis dropped below 50 km — aborting ***\n");
            mode = COMPLETE;
            break;
        }
    }

    campaign.total_passes = (int)campaign.passes.size();
    campaign.v_final = orb.v_peri;
    campaign.Mach_final = orb.v_peri / sqrt(JUPITER_GAMMA * JUPITER_RGAS * 166.0);
    campaign.final_orbit = orb;

    // Summary
    printf("\n  ╔════════════════════════════════════════════════════════╗\n");
    printf("  ║  CAMPAIGN SUMMARY                                    ║\n");
    printf("  ╠════════════════════════════════════════════════════════╣\n");
    printf("  ║  Passes:        %-6d (%.0f days)                   ║\n",
           campaign.total_passes, campaign.campaign_days);
    printf("  ║  Final orbit:   %.0f km × %.1f Rj                   ║\n",
           orb.alt_peri()/1000, orb.r_apo/R_J);
    printf("  ║  Final v_peri:  %.1f km/s (Mach %.0f)                ║\n",
           campaign.v_final/1000, campaign.Mach_final);
    printf("  ║  Collected:     %.1f tonnes (%.0f%% of tank)          ║\n",
           campaign.total_collected/1000, campaign.total_collected/veh.tank_capacity*100);
    printf("  ║    H2:          %.1f tonnes                          ║\n", campaign.total_h2/1000);
    printf("  ║    He:          %.1f tonnes                          ║\n", campaign.total_he/1000);
    printf("  ║  Heat load:     %.1f MJ total                        ║\n", campaign.total_heat_load);
    printf("  ║  Peak hull q:   %.2f kW/m²                          ║\n", campaign.max_heat_flux);
    printf("  ║  Hull status:   %s    ║\n",
           campaign.max_heat_flux < 1.0 ? "ALUMINUM SKIN — no TPS needed" :
           campaign.max_heat_flux < 50.0 ? "SiC tiles — lightweight      " :
           "NEEDS TPS — too hot               ");

    // Delta-V accounting
    double dv_saved = campaign.passes.front().orbit_pre.v_peri - campaign.v_final;
    printf("  ║  ΔV from drag:  %.1f km/s (FREE deceleration)       ║\n", dv_saved/1000);
    printf("  ╚════════════════════════════════════════════════════════╝\n");

    return campaign;
}

} // namespace aeropass

#endif // AEROPASS_SIM_HPP
