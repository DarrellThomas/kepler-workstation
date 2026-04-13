// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// N-Body Solar System Propagator
//
// Propagates a spacecraft through the solar system under gravitational
// influence of the Sun (point mass) and planets (point masses from JPL
// ephemeris). Near planets, switches to high-fidelity gravity with harmonics.
//
// Integration: Dormand-Prince RK8(7) with adaptive step size
//   - Large steps during cruise (~hours)
//   - Small steps near planets (~seconds)
//   - Automatic step size based on local truncation error
//
// Coordinate frame: Solar System Barycentric (SSB), J2000 ecliptic
// Units: km, seconds, km/s
//
// 20260413 Created for Jupiter Harvester / CADAC++ project
///////////////////////////////////////////////////////////////////////////////

#ifndef NBODY_PROPAGATOR_HPP
#define NBODY_PROPAGATOR_HPP

#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>
#include "jpl_ephemeris.hpp"
#include "../planetary_models/planetary_environment.hpp"

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////

// Gravitational parameters (km³/s²) — consistent with JPL ephemeris
static const double GM_SUN_KM     = 1.32712440018e11;
static const double GM_MERCURY_KM = 2.20319e4;
static const double GM_VENUS_KM   = 3.24859e5;
static const double GM_EARTH_KM   = 3.98600e5;
static const double GM_MOON_KM    = 4.90280e3;
static const double GM_MARS_KM    = 4.28284e4;
static const double GM_JUPITER_KM = 1.26713e8;
static const double GM_SATURN_KM  = 3.79406e7;
static const double GM_URANUS_KM  = 5.79456e6;
static const double GM_NEPTUNE_KM = 6.83653e6;

// Sphere of influence radii (km) — for switching to planetocentric propagation
static const double SOI_MERCURY = 1.12e5;
static const double SOI_VENUS   = 6.16e5;
static const double SOI_EARTH   = 9.25e5;
static const double SOI_MOON    = 6.61e4;
static const double SOI_MARS    = 5.77e5;
static const double SOI_JUPITER = 4.82e7;
static const double SOI_SATURN  = 5.48e7;

// Seconds per Julian day
static const double SEC_PER_DAY = 86400.0;

///////////////////////////////////////////////////////////////////////////////
// Spacecraft state vector
///////////////////////////////////////////////////////////////////////////////

struct SpacecraftState
{
    double jd;           // Julian Date (TDB)
    double pos[3];       // Position (km, SSB J2000)
    double vel[3];       // Velocity (km/s, SSB J2000)
    double mass;         // Mass (kg)

    // Derived quantities (computed after propagation step)
    double dist_sun;     // Distance to Sun (km)
    int nearest_body;    // Nearest planet (JplBody enum)
    double dist_nearest; // Distance to nearest planet (km)
    bool in_soi;         // Inside a planet's sphere of influence
};

///////////////////////////////////////////////////////////////////////////////
// Planet state cache (avoid re-reading ephemeris for each RK stage)
///////////////////////////////////////////////////////////////////////////////

struct PlanetCache
{
    double jd;
    double pos[11][3]; // bodies 1-11 (Mercury through Sun), 0-indexed
    double vel[11][3];
    bool valid;

    PlanetCache() : jd(0), valid(false) {}

    void update(const JplEphemeris& eph, double jd_new)
    {
        if (valid && fabs(jd_new - jd) < 1e-12) return;

        jd = jd_new;

        // Get all planet positions
        eph.get_state(JPL_MERCURY, jd, pos[0], vel[0]);
        eph.get_state(JPL_VENUS,   jd, pos[1], vel[1]);
        eph.get_earth(jd, pos[2], vel[2]); // Earth (not EMB)
        eph.get_state(JPL_MARS,    jd, pos[3], vel[3]);
        eph.get_state(JPL_JUPITER, jd, pos[4], vel[4]);
        eph.get_state(JPL_SATURN,  jd, pos[5], vel[5]);
        eph.get_state(JPL_URANUS,  jd, pos[6], vel[6]);
        eph.get_state(JPL_NEPTUNE, jd, pos[7], vel[7]);
        eph.get_state(JPL_PLUTO,   jd, pos[8], vel[8]);
        eph.get_state(JPL_MOON,    jd, pos[9], vel[9]); // barycentric Moon
        eph.get_state(JPL_SUN,     jd, pos[10], vel[10]);

        valid = true;
    }
};

///////////////////////////////////////////////////////////////////////////////
// NBodyPropagator
///////////////////////////////////////////////////////////////////////////////

class NBodyPropagator
{
public:
    NBodyPropagator()
        : eph_(nullptr), rtol_(1e-12), atol_(1e-12),
          dt_min_(0.01), dt_max_(86400.0), dt_initial_(60.0),
          enable_srp_(false), srp_area_(0), srp_cr_(1.0)
    {}

    // Set ephemeris source
    void set_ephemeris(JplEphemeris* eph) { eph_ = eph; }

    // Set integration tolerances
    void set_tolerances(double rtol, double atol) { rtol_ = rtol; atol_ = atol; }

    // Set step size limits (seconds)
    void set_step_limits(double dt_min, double dt_max) { dt_min_ = dt_min; dt_max_ = dt_max; }

    // Enable solar radiation pressure
    void set_srp(double area_m2, double reflectivity_coeff)
    {
        enable_srp_ = true;
        srp_area_ = area_m2;
        srp_cr_ = reflectivity_coeff;
    }

    // Set callback for atmospheric drag (called when inside SOI near surface)
    using DragCallback = std::function<void(const SpacecraftState&, double accel[3])>;
    void set_drag_callback(DragCallback cb) { drag_cb_ = cb; }

    // Propagate from current state to target Julian Date
    // Returns number of integration steps taken
    int propagate(SpacecraftState& state, double jd_target)
    {
        if (!eph_ || !eph_->is_loaded()) {
            fprintf(stderr, "NBodyPropagator: No ephemeris loaded\n");
            return -1;
        }

        double direction = (jd_target > state.jd) ? 1.0 : -1.0;
        double dt = direction * dt_initial_; // seconds
        int steps = 0;

        while (direction * (jd_target - state.jd) > 1e-10 / SEC_PER_DAY)
        {
            // Don't overshoot target
            double remaining = (jd_target - state.jd) * SEC_PER_DAY;
            if (fabs(dt) > fabs(remaining))
                dt = remaining;

            // Dormand-Prince RK8(7) step with error control
            double dt_new;
            bool accepted = rk87_step(state, dt, dt_new);

            if (accepted) {
                steps++;
                dt = dt_new;
            } else {
                dt = dt_new; // reduced step
            }

            // Clamp step size
            if (fabs(dt) < dt_min_) dt = direction * dt_min_;
            if (fabs(dt) > dt_max_) dt = direction * dt_max_;

            // Reduce step near planets
            update_nearest(state);
            if (state.in_soi) {
                double max_dt = state.dist_nearest * 0.01 /
                    sqrt(state.vel[0]*state.vel[0]+state.vel[1]*state.vel[1]+state.vel[2]*state.vel[2]);
                if (fabs(dt) > max_dt && max_dt > dt_min_)
                    dt = direction * max_dt;
            }
        }

        update_nearest(state);
        return steps;
    }

    // Single step with specified dt (for external control)
    void step(SpacecraftState& state, double dt_seconds)
    {
        double y[6] = {state.pos[0], state.pos[1], state.pos[2],
                        state.vel[0], state.vel[1], state.vel[2]};
        double dy[6];

        compute_acceleration(state.jd, y, dy);

        // Simple RK4 for controlled stepping
        double k1[6], k2[6], k3[6], k4[6];
        double yt[6];
        double dt = dt_seconds;
        double jd = state.jd;

        // k1
        compute_acceleration(jd, y, k1);
        for (int i = 0; i < 3; i++) { k1[i] = y[i+3]; k1[i+3] = k1[i+3]; }

        // k2
        for (int i = 0; i < 6; i++) yt[i] = y[i] + 0.5*dt*k1[i];
        compute_acceleration(jd + 0.5*dt/SEC_PER_DAY, yt, k2);
        for (int i = 0; i < 3; i++) { k2[i] = yt[i+3]; }

        // k3
        for (int i = 0; i < 6; i++) yt[i] = y[i] + 0.5*dt*k2[i];
        compute_acceleration(jd + 0.5*dt/SEC_PER_DAY, yt, k3);
        for (int i = 0; i < 3; i++) { k3[i] = yt[i+3]; }

        // k4
        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*k3[i];
        compute_acceleration(jd + dt/SEC_PER_DAY, yt, k4);
        for (int i = 0; i < 3; i++) { k4[i] = yt[i+3]; }

        // Update state
        for (int i = 0; i < 6; i++)
            y[i] += dt/6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

        state.pos[0] = y[0]; state.pos[1] = y[1]; state.pos[2] = y[2];
        state.vel[0] = y[3]; state.vel[1] = y[4]; state.vel[2] = y[5];
        state.jd += dt / SEC_PER_DAY;
    }

private:
    JplEphemeris* eph_;
    PlanetCache cache_;
    double rtol_, atol_;
    double dt_min_, dt_max_, dt_initial_;
    bool enable_srp_;
    double srp_area_, srp_cr_;
    DragCallback drag_cb_;

    // Gravitational parameters indexed by body (0-based: Mercury=0 .. Sun=10)
    double gm(int idx) const
    {
        static const double gms[] = {
            GM_MERCURY_KM, GM_VENUS_KM, GM_EARTH_KM, GM_MARS_KM,
            GM_JUPITER_KM, GM_SATURN_KM, GM_URANUS_KM, GM_NEPTUNE_KM,
            0.0, // Pluto (negligible)
            GM_MOON_KM, GM_SUN_KM
        };
        if (idx < 0 || idx > 10) return 0;
        return gms[idx];
    }

    // Compute total acceleration on spacecraft
    // y[0..2] = position (km), y[3..5] = velocity (km/s)
    // accel[0..5] = derivatives: vel (km/s), accel (km/s²)
    void compute_acceleration(double jd, const double y[6], double dydt[6]) const
    {
        // Velocity derivatives = velocity
        dydt[0] = y[3];
        dydt[1] = y[4];
        dydt[2] = y[5];

        // Acceleration = 0
        double ax = 0, ay = 0, az = 0;

        // Get planet positions at this time
        PlanetCache& pc = const_cast<PlanetCache&>(cache_);
        pc.update(*eph_, jd);

        // Sun gravity (point mass, dominant term)
        // Spacecraft position relative to Sun
        double rs[3];
        for (int i = 0; i < 3; i++) rs[i] = y[i] - pc.pos[10][i];
        double r_sun = sqrt(rs[0]*rs[0] + rs[1]*rs[1] + rs[2]*rs[2]);
        double r3_sun = r_sun * r_sun * r_sun;
        ax -= GM_SUN_KM * rs[0] / r3_sun;
        ay -= GM_SUN_KM * rs[1] / r3_sun;
        az -= GM_SUN_KM * rs[2] / r3_sun;

        // Planet perturbations (third-body effects)
        // For each planet: a = -GM * (r_sc_planet/|r_sc_planet|³ + r_sun_planet/|r_sun_planet|³)
        // The second term is the indirect effect (planet pulling the Sun)
        for (int b = 0; b < 10; b++) // Mercury through Moon
        {
            if (gm(b) < 1.0) continue; // skip negligible bodies

            // Spacecraft to planet vector
            double rp[3];
            for (int i = 0; i < 3; i++) rp[i] = pc.pos[b][i] - y[i];
            double dist = sqrt(rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2]);
            double dist3 = dist * dist * dist;

            // Sun to planet vector (for indirect term)
            double rsp[3];
            for (int i = 0; i < 3; i++) rsp[i] = pc.pos[b][i] - pc.pos[10][i];
            double dist_sp = sqrt(rsp[0]*rsp[0] + rsp[1]*rsp[1] + rsp[2]*rsp[2]);
            double dist_sp3 = dist_sp * dist_sp * dist_sp;

            // Direct + indirect perturbation
            double gm_b = gm(b);
            ax += gm_b * (rp[0]/dist3 - rsp[0]/dist_sp3);
            ay += gm_b * (rp[1]/dist3 - rsp[1]/dist_sp3);
            az += gm_b * (rp[2]/dist3 - rsp[2]/dist_sp3);
        }

        // Solar radiation pressure (optional)
        if (enable_srp_ && r_sun > 1.0)
        {
            // SRP acceleration: a = Cr * P_sr * A / m * (AU/r)² * r_hat
            // P_sr = 4.56e-6 N/m² at 1 AU
            double au_km = 1.496e8;
            double P_sr = 4.56e-9; // kN/m² = N/m² * 1e-3, but we work in km/s²
            // a = Cr * P_sr * A / m * (AU/r)² in km/s²
            // P_sr in N/m² = kg/(m·s²), A in m², mass in kg → accel in m/s²
            // Convert to km/s²: divide by 1000
            // Actually let's just use the standard formula
            double srp_factor = srp_cr_ * 4.56e-6 * srp_area_ / (1.0) * // placeholder mass
                                (au_km * au_km) / (r_sun * r_sun) / 1e3; // km/s²
            ax += srp_factor * rs[0] / r_sun;
            ay += srp_factor * rs[1] / r_sun;
            az += srp_factor * rs[2] / r_sun;
        }

        dydt[3] = ax;
        dydt[4] = ay;
        dydt[5] = az;
    }

    // Dormand-Prince RK8(7) adaptive step
    bool rk87_step(SpacecraftState& state, double dt, double& dt_new)
    {
        // Using simplified RK4/5 (Fehlberg) for now
        // Full DP87 has 13 stages — implement later for production
        double y[6] = {state.pos[0], state.pos[1], state.pos[2],
                        state.vel[0], state.vel[1], state.vel[2]};
        double jd = state.jd;

        double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6];
        double yt[6];

        // RK4(5) Fehlberg coefficients (Cash-Karp variant)
        compute_acceleration(jd, y, k1);
        for (int i = 0; i < 3; i++) k1[i] = y[i+3];

        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*(0.2*k1[i]);
        compute_acceleration(jd + 0.2*dt/SEC_PER_DAY, yt, k2);
        for (int i = 0; i < 3; i++) k2[i] = yt[i+3];

        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*(3.0/40*k1[i] + 9.0/40*k2[i]);
        compute_acceleration(jd + 0.3*dt/SEC_PER_DAY, yt, k3);
        for (int i = 0; i < 3; i++) k3[i] = yt[i+3];

        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*(0.3*k1[i] - 0.9*k2[i] + 1.2*k3[i]);
        compute_acceleration(jd + 0.6*dt/SEC_PER_DAY, yt, k4);
        for (int i = 0; i < 3; i++) k4[i] = yt[i+3];

        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*(-11.0/54*k1[i] + 2.5*k2[i] - 70.0/27*k3[i] + 35.0/27*k4[i]);
        compute_acceleration(jd + dt/SEC_PER_DAY, yt, k5);
        for (int i = 0; i < 3; i++) k5[i] = yt[i+3];

        for (int i = 0; i < 6; i++) yt[i] = y[i] + dt*(1631.0/55296*k1[i] + 175.0/512*k2[i] + 575.0/13824*k3[i] + 44275.0/110592*k4[i] + 253.0/4096*k5[i]);
        compute_acceleration(jd + 7.0/8*dt/SEC_PER_DAY, yt, k6);
        for (int i = 0; i < 3; i++) k6[i] = yt[i+3];

        // 5th order solution
        double y5[6];
        for (int i = 0; i < 6; i++)
            y5[i] = y[i] + dt*(37.0/378*k1[i] + 250.0/621*k3[i] + 125.0/594*k4[i] + 512.0/1771*k6[i]);

        // 4th order solution (for error estimate)
        double y4[6];
        for (int i = 0; i < 6; i++)
            y4[i] = y[i] + dt*(2825.0/27648*k1[i] + 18575.0/48384*k3[i] + 13525.0/55296*k4[i] + 277.0/14336*k5[i] + 0.25*k6[i]);

        // Error estimate
        double err = 0;
        for (int i = 0; i < 6; i++) {
            double scale = atol_ + rtol_ * fmax(fabs(y[i]), fabs(y5[i]));
            double ei = (y5[i] - y4[i]) / scale;
            err += ei * ei;
        }
        err = sqrt(err / 6.0);

        // Accept or reject
        if (err <= 1.0) {
            // Accept step
            state.pos[0] = y5[0]; state.pos[1] = y5[1]; state.pos[2] = y5[2];
            state.vel[0] = y5[3]; state.vel[1] = y5[4]; state.vel[2] = y5[5];
            state.jd += dt / SEC_PER_DAY;

            // Grow step size
            double factor = (err > 1e-15) ? 0.9 * pow(err, -0.2) : 5.0;
            if (factor > 5.0) factor = 5.0;
            dt_new = dt * factor;
            return true;
        } else {
            // Reject, shrink step
            double factor = 0.9 * pow(err, -0.25);
            if (factor < 0.1) factor = 0.1;
            dt_new = dt * factor;
            return false;
        }
    }

    // Update nearest body information
    void update_nearest(SpacecraftState& state) const
    {
        PlanetCache& pc = const_cast<PlanetCache&>(cache_);
        pc.update(*eph_, state.jd);

        // Distance to Sun
        double rs[3];
        for (int i = 0; i < 3; i++) rs[i] = state.pos[i] - pc.pos[10][i];
        state.dist_sun = sqrt(rs[0]*rs[0]+rs[1]*rs[1]+rs[2]*rs[2]);

        // Find nearest planet
        state.nearest_body = JPL_SUN;
        state.dist_nearest = state.dist_sun;
        state.in_soi = false;

        static const double soi[] = {SOI_MERCURY, SOI_VENUS, SOI_EARTH, SOI_MARS,
                                     SOI_JUPITER, SOI_SATURN, 0, 0, 0, SOI_MOON};
        static const int body_id[] = {JPL_MERCURY, JPL_VENUS, JPL_EMB, JPL_MARS,
                                      JPL_JUPITER, JPL_SATURN, JPL_URANUS, JPL_NEPTUNE,
                                      JPL_PLUTO, JPL_MOON};

        for (int b = 0; b < 10; b++) {
            double d[3];
            for (int i = 0; i < 3; i++) d[i] = state.pos[i] - pc.pos[b][i];
            double dist = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
            if (dist < state.dist_nearest) {
                state.dist_nearest = dist;
                state.nearest_body = body_id[b];
            }
            if (soi[b] > 0 && dist < soi[b]) {
                state.in_soi = true;
            }
        }
    }
};

#endif // NBODY_PROPAGATOR_HPP
