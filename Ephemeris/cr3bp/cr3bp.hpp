// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Circular Restricted Three-Body Problem (CR3BP)
//
// Core dynamics for the CR3BP in the synodic (rotating) reference frame
// with nondimensional units.
//
// Reference: Koon, Lo, Marsden, Ross — "Dynamical Systems, the Three-Body
//   Problem and Space Mission Design" (2011), Chapter 2.
//   Szebehely — "Theory of Orbits" (1967).
//
// Conventions:
//   - m1 (primary, larger mass) at (-mu, 0, 0)
//   - m2 (secondary, smaller mass) at (1-mu, 0, 0)
//   - mu = m2 / (m1 + m2), mu in [0, 0.5]
//   - Unit of length: distance between m1 and m2
//   - Unit of time: chosen so orbital period = 2*pi
//   - Unit of mass: m1 + m2 = 1
//   - State vector: [x, y, z, vx, vy, vz]
//
// Compiles standalone: g++ -std=c++11 -O2
//
// 20260508 Created for Kepler Workstation CR3BP module
///////////////////////////////////////////////////////////////////////////////

#ifndef CR3BP_HPP
#define CR3BP_HPP

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////

static const double CR3BP_PI = 3.14159265358979323846;

///////////////////////////////////////////////////////////////////////////////
// Known system mass parameters (mu = m2/(m1+m2))
// Source: KoLoMaRo Table 2.2.1 and JPL SSD
///////////////////////////////////////////////////////////////////////////////

struct CR3BPSystem
{
    const char* name;
    double mu;          // mass parameter
    double L_km;        // characteristic length (km)
    double T_sec;       // characteristic time (seconds)

    // Derived: velocity scale V = L / (T / 2*pi) = 2*pi*L / T
    double V_km_s() const { return 2.0 * CR3BP_PI * L_km / T_sec; }
};

// Standard systems from KoLoMaRo Table 2.2.1
static const CR3BPSystem SYS_SUN_JUPITER    = {"Sun-Jupiter",       9.537e-4,  7.784e8,  3.733e8};
static const CR3BPSystem SYS_SUN_EARTH      = {"Sun-(Earth+Moon)",  3.036e-6,  1.496e8,  3.147e7};
static const CR3BPSystem SYS_EARTH_MOON     = {"Earth-Moon",        1.215e-2,  3.850e5,  2.361e6};
static const CR3BPSystem SYS_MARS_PHOBOS    = {"Mars-Phobos",       1.667e-8,  9.380e3,  2.749e4};
static const CR3BPSystem SYS_JUPITER_IO     = {"Jupiter-Io",        4.704e-5,  4.218e5,  1.524e5};
static const CR3BPSystem SYS_JUPITER_EUROPA = {"Jupiter-Europa",    2.528e-5,  6.711e5,  3.060e5};
static const CR3BPSystem SYS_JUPITER_GANYMEDE={"Jupiter-Ganymede",  7.804e-5,  1.070e6,  6.165e5};
static const CR3BPSystem SYS_JUPITER_CALLISTO={"Jupiter-Callisto",  5.667e-5,  1.883e6,  1.438e6};
static const CR3BPSystem SYS_SATURN_TITAN    = {"Saturn-Titan",      2.366e-4,  1.222e6,  1.374e6};
static const CR3BPSystem SYS_SATURN_ENCELADUS={"Saturn-Enceladus",  1.901e-7,  2.380e5,  1.184e5};
static const CR3BPSystem SYS_SATURN_RHEA    = {"Saturn-Rhea",       4.058e-6,  5.271e5,  3.903e5};
static const CR3BPSystem SYS_SATURN_DIONE   = {"Saturn-Dione",      1.927e-6,  3.774e5,  2.364e5};

///////////////////////////////////////////////////////////////////////////////
// Celestial Body Physical Parameters
// Source: JPL SSD, JUP365 ephemeris (Jacobson 2021)
// https://ssd.jpl.nasa.gov/sats/phys_par/
// https://ssd.jpl.nasa.gov/sats/elem/
///////////////////////////////////////////////////////////////////////////////

struct CelestialBody
{
    const char* name;
    double GM;           // gravitational parameter (km^3/s^2)
    double radius;       // mean radius (km)
    double a;            // semi-major axis of orbit (km)
    double ecc;          // orbital eccentricity
    double period_days;  // orbital period (days)
};

// Jupiter system
static const CelestialBody BODY_JUPITER   = {"Jupiter",   1.26713e8, 71492.0, 7.784e8,   0.0489, 4332.59};
static const CelestialBody BODY_IO        = {"Io",        5959.916,  1821.49, 421800.0,  0.004,  1.769138};
static const CelestialBody BODY_EUROPA    = {"Europa",    3202.712,  1560.80, 671100.0,  0.009,  3.551181};
static const CelestialBody BODY_GANYMEDE  = {"Ganymede",  9887.833,  2631.20, 1070400.0, 0.001,  7.154553};
static const CelestialBody BODY_CALLISTO  = {"Callisto",  7179.283,  2410.30, 1882700.0, 0.007,  16.689018};

// Earth system
static const CelestialBody BODY_EARTH     = {"Earth",     3.98600e5, 6371.0,  1.496e8,   0.0167, 365.256};
static const CelestialBody BODY_MOON      = {"Moon",      4902.800,  1737.4,  384400.0,  0.0549, 27.322};

// Saturn system
// Source: JPL SSD, Cassini mission (Jacobson et al. 2006, 2008)
static const CelestialBody BODY_SATURN    = {"Saturn",    3.7931e7,  60268.0,  1.4267e9,  0.0565, 10759.22};
static const CelestialBody BODY_MIMAS     = {"Mimas",     2.504,     198.2,   185539.0,  0.020,  0.942};
static const CelestialBody BODY_ENCELADUS = {"Enceladus", 7.211,     252.1,   238042.0,  0.005,  1.370};
static const CelestialBody BODY_TETHYS    = {"Tethys",    41.21,     533.0,   294672.0,  0.001,  1.888};
static const CelestialBody BODY_DIONE     = {"Dione",     73.11,     561.4,   377415.0,  0.002,  2.737};
static const CelestialBody BODY_RHEA      = {"Rhea",      153.94,    764.3,   527068.0,  0.001,  4.518};
static const CelestialBody BODY_TITAN     = {"Titan",     8978.14,   2574.7,  1221870.0, 0.029,  15.945};
static const CelestialBody BODY_IAPETUS   = {"Iapetus",   120.51,    735.6,   3560854.0, 0.029,  79.322};

// Sun
static const CelestialBody BODY_SUN       = {"Sun",       1.32712e11, 696000.0, 0.0, 0.0, 0.0};

///////////////////////////////////////////////////////////////////////////////
// CR3BP State
///////////////////////////////////////////////////////////////////////////////

struct CR3BPState
{
    double x, y, z;     // position (nondimensional, rotating frame)
    double vx, vy, vz;  // velocity (nondimensional, rotating frame)

    // Pack/unpack to array
    void to_array(double s[6]) const
    {
        s[0] = x; s[1] = y; s[2] = z;
        s[3] = vx; s[4] = vy; s[5] = vz;
    }

    void from_array(const double s[6])
    {
        x = s[0]; y = s[1]; z = s[2];
        vx = s[3]; vy = s[4]; vz = s[5];
    }
};

///////////////////////////////////////////////////////////////////////////////
// CR3BP Core Functions
///////////////////////////////////////////////////////////////////////////////

// Distance from state to m1 (primary, at -mu, 0, 0)
inline double cr3bp_r1(double x, double y, double z, double mu)
{
    double dx = x + mu;
    return sqrt(dx*dx + y*y + z*z);
}

// Distance from state to m2 (secondary, at 1-mu, 0, 0)
inline double cr3bp_r2(double x, double y, double z, double mu)
{
    double dx = x - 1.0 + mu;
    return sqrt(dx*dx + y*y + z*z);
}

// Pseudo-potential Omega(x,y,z)
// Omega = 1/2(x^2 + y^2) + (1-mu)/r1 + mu/r2 + 1/2 mu(1-mu)
// Note: -Omega = U_bar from KoLoMaRo Eq. 2.3.9
inline double cr3bp_omega(double x, double y, double z, double mu)
{
    double r1 = cr3bp_r1(x, y, z, mu);
    double r2 = cr3bp_r2(x, y, z, mu);
    return 0.5*(x*x + y*y) + (1.0-mu)/r1 + mu/r2 + 0.5*mu*(1.0-mu);
}

// Jacobi constant: C = 2*Omega - v^2
// Conserved quantity (energy integral) of the CR3BP
inline double cr3bp_jacobi(const CR3BPState& s, double mu)
{
    double v2 = s.vx*s.vx + s.vy*s.vy + s.vz*s.vz;
    return 2.0 * cr3bp_omega(s.x, s.y, s.z, mu) - v2;
}

///////////////////////////////////////////////////////////////////////////////
// CR3BP Equations of Motion
//
// State: s[0..5] = [x, y, z, vx, vy, vz]
// Derivatives: dsdt[0..5]
//
// From KoLoMaRo Eq. 2.3.8:
//   x'' - 2y' = Omega_x
//   y'' + 2x' = Omega_y
//   z''       = Omega_z
///////////////////////////////////////////////////////////////////////////////

inline void cr3bp_eom(double mu, const double s[6], double dsdt[6])
{
    double x  = s[0], y  = s[1], z  = s[2];
    double vx = s[3], vy = s[4], vz = s[5];

    double mu1 = 1.0 - mu;  // mass of primary
    double dx1 = x + mu;    // x - x_m1
    double dx2 = x - mu1;   // x - x_m2

    double r1 = sqrt(dx1*dx1 + y*y + z*z);
    double r2 = sqrt(dx2*dx2 + y*y + z*z);
    double r1_3 = r1 * r1 * r1;
    double r2_3 = r2 * r2 * r2;

    // Partial derivatives of pseudo-potential Omega
    double Omega_x = x - mu1*dx1/r1_3 - mu*dx2/r2_3;
    double Omega_y = y - mu1*y/r1_3    - mu*y/r2_3;
    double Omega_z =   - mu1*z/r1_3    - mu*z/r2_3;

    // Equations of motion
    dsdt[0] = vx;
    dsdt[1] = vy;
    dsdt[2] = vz;
    dsdt[3] = 2.0*vy + Omega_x;
    dsdt[4] = -2.0*vx + Omega_y;
    dsdt[5] = Omega_z;
}

///////////////////////////////////////////////////////////////////////////////
// Second partial derivatives of Omega (gravity gradient tensor)
// Needed for STM propagation and linearization
//
// From KoLoMaRo Section 2.6
///////////////////////////////////////////////////////////////////////////////

struct OmegaPartials
{
    double Uxx, Uyy, Uzz;
    double Uxy, Uxz, Uyz;
};

inline OmegaPartials cr3bp_partials(double x, double y, double z, double mu)
{
    double mu1 = 1.0 - mu;
    double dx1 = x + mu;
    double dx2 = x - 1.0 + mu;

    double r1 = sqrt(dx1*dx1 + y*y + z*z);
    double r2 = sqrt(dx2*dx2 + y*y + z*z);
    double r1_3 = r1*r1*r1;
    double r2_3 = r2*r2*r2;
    double r1_5 = r1_3*r1*r1;
    double r2_5 = r2_3*r2*r2;

    OmegaPartials p;
    p.Uxx = 1.0 - mu1/r1_3 + 3.0*mu1*dx1*dx1/r1_5 - mu/r2_3 + 3.0*mu*dx2*dx2/r2_5;
    p.Uyy = 1.0 - mu1/r1_3 + 3.0*mu1*y*y/r1_5     - mu/r2_3 + 3.0*mu*y*y/r2_5;
    p.Uzz =     - mu1/r1_3 + 3.0*mu1*z*z/r1_5      - mu/r2_3 + 3.0*mu*z*z/r2_5;
    p.Uxy = 3.0*mu1*dx1*y/r1_5 + 3.0*mu*dx2*y/r2_5;
    p.Uxz = 3.0*mu1*dx1*z/r1_5 + 3.0*mu*dx2*z/r2_5;
    p.Uyz = 3.0*mu1*y*z/r1_5   + 3.0*mu*y*z/r2_5;

    return p;
}

///////////////////////////////////////////////////////////////////////////////
// CR3BP + STM Equations of Motion (42-state system)
//
// State: s[0..5]  = [x, y, z, vx, vy, vz]  (trajectory)
//        s[6..41] = Phi (6x6 STM, row-major)
//
// dPhi/dt = A(t) * Phi
// where A is the Jacobian of the CR3BP EOM
///////////////////////////////////////////////////////////////////////////////

inline void cr3bp_eom_stm(double mu, const double s[42], double dsdt[42])
{
    // Trajectory derivatives (same as cr3bp_eom)
    cr3bp_eom(mu, s, dsdt);

    // Build Jacobian A = [0 I; G C] where:
    //   C = [0 2 0; -2 0 0; 0 0 0]  (Coriolis)
    //   G = Hessian of Omega
    OmegaPartials p = cr3bp_partials(s[0], s[1], s[2], mu);

    // A matrix (6x6), stored as A[i][j]
    double A[6][6];
    memset(A, 0, sizeof(A));

    // Upper-right block: identity
    A[0][3] = 1.0;  A[1][4] = 1.0;  A[2][5] = 1.0;

    // Lower-left block: gravity gradient
    A[3][0] = p.Uxx;  A[3][1] = p.Uxy;  A[3][2] = p.Uxz;
    A[4][0] = p.Uxy;  A[4][1] = p.Uyy;  A[4][2] = p.Uyz;
    A[5][0] = p.Uxz;  A[5][1] = p.Uyz;  A[5][2] = p.Uzz;

    // Lower-right block: Coriolis
    A[3][4] = 2.0;
    A[4][3] = -2.0;

    // STM derivative: dPhi/dt = A * Phi
    // Phi is s[6..41], stored row-major: Phi(i,j) = s[6 + i*6 + j]
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double sum = 0.0;
            for (int k = 0; k < 6; k++) {
                sum += A[i][k] * s[6 + k*6 + j];
            }
            dsdt[6 + i*6 + j] = sum;
        }
    }
}

// Initialize STM to identity
inline void cr3bp_stm_identity(double s[42], const CR3BPState& state)
{
    state.to_array(s);
    for (int i = 0; i < 36; i++) s[6+i] = 0.0;
    for (int i = 0; i < 6; i++) s[6 + i*6 + i] = 1.0; // identity
}

// Extract STM from 42-state vector
inline void cr3bp_stm_extract(const double s[42], double phi[6][6])
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            phi[i][j] = s[6 + i*6 + j];
}

///////////////////////////////////////////////////////////////////////////////
// RK4/5 Integrator (Cash-Karp)
//
// Adaptive step-size Runge-Kutta for CR3BP propagation.
// Matches the integrator style in nbody_propagator.hpp.
///////////////////////////////////////////////////////////////////////////////

// Function pointer type for EOM
typedef void (*CR3BPEomFunc)(double mu, const double* s, double* dsdt);

struct CR3BPIntegrator
{
    double rtol;
    double atol;
    double dt_min;
    double dt_max;

    CR3BPIntegrator()
        : rtol(1e-12), atol(1e-14), dt_min(1e-12), dt_max(0.5) {}

    // Propagate state from t to t_target using adaptive RK4/5
    // Returns final time reached
    double propagate(CR3BPEomFunc eom, double mu, double* s, int n,
                     double t, double t_target, double dt_init = 0.001)
    {
        double direction = (t_target > t) ? 1.0 : -1.0;
        double dt = direction * fabs(dt_init);

        while (direction * (t_target - t) > 1e-15)
        {
            // Don't overshoot
            double remaining = t_target - t;
            if (fabs(dt) > fabs(remaining))
                dt = remaining;

            double dt_new;
            bool accepted = rk45_step(eom, mu, s, n, t, dt, dt_new);

            if (accepted) {
                t += dt;
                dt = dt_new;
            } else {
                dt = dt_new;
            }

            // Clamp
            if (fabs(dt) < dt_min) dt = direction * dt_min;
            if (fabs(dt) > dt_max) dt = direction * dt_max;
        }
        return t;
    }

    // Single RK4 step (fixed step, for simple use)
    static void rk4_step(CR3BPEomFunc eom, double mu, double* s, int n,
                         double t, double dt)
    {
        std::vector<double> k1(n), k2(n), k3(n), k4(n), tmp(n);

        eom(mu, s, k1.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + 0.5*dt*k1[i];
        eom(mu, tmp.data(), k2.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + 0.5*dt*k2[i];
        eom(mu, tmp.data(), k3.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*k3[i];
        eom(mu, tmp.data(), k4.data());

        for (int i = 0; i < n; i++)
            s[i] += dt/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

private:
    // Cash-Karp RK4(5) adaptive step
    bool rk45_step(CR3BPEomFunc eom, double mu, double* s, int n,
                   double t, double dt, double& dt_new)
    {
        std::vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n);
        std::vector<double> tmp(n), y5(n), y4(n);

        eom(mu, s, k1.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*(0.2*k1[i]);
        eom(mu, tmp.data(), k2.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*(3.0/40*k1[i] + 9.0/40*k2[i]);
        eom(mu, tmp.data(), k3.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*(0.3*k1[i] - 0.9*k2[i] + 1.2*k3[i]);
        eom(mu, tmp.data(), k4.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*(-11.0/54*k1[i] + 2.5*k2[i] - 70.0/27*k3[i] + 35.0/27*k4[i]);
        eom(mu, tmp.data(), k5.data());

        for (int i = 0; i < n; i++) tmp[i] = s[i] + dt*(1631.0/55296*k1[i] + 175.0/512*k2[i] + 575.0/13824*k3[i] + 44275.0/110592*k4[i] + 253.0/4096*k5[i]);
        eom(mu, tmp.data(), k6.data());

        // 5th order solution
        for (int i = 0; i < n; i++)
            y5[i] = s[i] + dt*(37.0/378*k1[i] + 250.0/621*k3[i] + 125.0/594*k4[i] + 512.0/1771*k6[i]);

        // 4th order solution
        for (int i = 0; i < n; i++)
            y4[i] = s[i] + dt*(2825.0/27648*k1[i] + 18575.0/48384*k3[i] + 13525.0/55296*k4[i] + 277.0/14336*k5[i] + 0.25*k6[i]);

        // Error estimate
        // For 42-state (trajectory + STM), use all states to control STM accuracy
        int n_err = n;
        double err = 0;
        for (int i = 0; i < n_err; i++) {
            double scale = atol + rtol * fmax(fabs(s[i]), fabs(y5[i]));
            double ei = (y5[i] - y4[i]) / scale;
            err += ei * ei;
        }
        err = sqrt(err / n_err);

        if (err <= 1.0) {
            // Accept
            for (int i = 0; i < n; i++) s[i] = y5[i];
            double factor = (err > 1e-15) ? 0.9 * pow(err, -0.2) : 5.0;
            if (factor > 5.0) factor = 5.0;
            dt_new = dt * factor;
            return true;
        } else {
            // Reject
            double factor = 0.9 * pow(err, -0.25);
            if (factor < 0.1) factor = 0.1;
            dt_new = dt * factor;
            return false;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////
// Convenience propagation functions
///////////////////////////////////////////////////////////////////////////////

// Propagate CR3BPState for time dt (nondimensional)
inline CR3BPState cr3bp_propagate(double mu, const CR3BPState& s0, double dt,
                                  double step = 0.01)
{
    double s[6];
    s0.to_array(s);
    CR3BPIntegrator integ;
    integ.propagate(cr3bp_eom, mu, s, 6, 0.0, dt, step);
    CR3BPState sf;
    sf.from_array(s);
    return sf;
}

// Propagate with STM, returns final state and STM
inline CR3BPState cr3bp_propagate_stm(double mu, const CR3BPState& s0,
                                      double dt, double phi[6][6],
                                      double step = 0.01)
{
    double s[42];
    cr3bp_stm_identity(s, s0);
    CR3BPIntegrator integ;
    integ.propagate(cr3bp_eom_stm, mu, s, 42, 0.0, dt, step);
    CR3BPState sf;
    sf.from_array(s);
    cr3bp_stm_extract(s, phi);
    return sf;
}

// Propagate and record trajectory at fixed intervals
inline std::vector<CR3BPState> cr3bp_propagate_record(
    double mu, const CR3BPState& s0, double dt_total,
    int n_points, double step = 0.001)
{
    std::vector<CR3BPState> traj;
    traj.reserve(n_points + 1);

    double s[6];
    s0.to_array(s);
    CR3BPIntegrator integ;

    double dt_out = dt_total / n_points;
    double t = 0.0;

    CR3BPState st;
    st.from_array(s);
    traj.push_back(st);

    for (int i = 0; i < n_points; i++) {
        t = integ.propagate(cr3bp_eom, mu, s, 6, t, t + dt_out, fabs(step));
        st.from_array(s);
        traj.push_back(st);
    }
    return traj;
}

///////////////////////////////////////////////////////////////////////////////
// Lagrange Point Solver
//
// L1, L2, L3: collinear points on x-axis (y = z = 0)
//   Solved via Newton's method on dOmega/dx = 0 along y=0.
//   Initial guess from Hill sphere approximation (KoLoMaRo Eq. 2.5.3)
//
// L4, L5: equilateral triangle points
//   L4 = (1/2 - mu, +sqrt(3)/2, 0)
//   L5 = (1/2 - mu, -sqrt(3)/2, 0)
///////////////////////////////////////////////////////////////////////////////

struct LagrangePoints
{
    double L1[3]; // {x, y, z}
    double L2[3];
    double L3[3];
    double L4[3];
    double L5[3];
};

// dOmega/dx along y=z=0
// Omega(x,0,0) = 1/2 x^2 + (1-mu)/|x+mu| + mu/|x-1+mu| + 1/2 mu(1-mu)
// dOmega/dx = x - (1-mu)(x+mu)/|x+mu|^3 - mu(x-1+mu)/|x-1+mu|^3
inline double lagrange_dOmega_dx(double x, double mu)
{
    double mu1 = 1.0 - mu;
    double dx1 = x + mu;
    double dx2 = x - 1.0 + mu;
    double r1 = fabs(dx1);
    double r2 = fabs(dx2);
    double r1_3 = r1*r1*r1;
    double r2_3 = r2*r2*r2;
    return x - mu1*dx1/r1_3 - mu*dx2/r2_3;
}

// d2Omega/dx2 along y=z=0
inline double lagrange_d2Omega_dx2(double x, double mu)
{
    double mu1 = 1.0 - mu;
    double dx1 = x + mu;
    double dx2 = x - 1.0 + mu;
    double r1 = fabs(dx1);
    double r2 = fabs(dx2);
    double r1_3 = r1*r1*r1;
    double r2_3 = r2*r2*r2;
    double r1_5 = r1_3*r1*r1;
    double r2_5 = r2_3*r2*r2;
    return 1.0 - mu1/r1_3 + 3.0*mu1*dx1*dx1/r1_5 - mu/r2_3 + 3.0*mu*dx2*dx2/r2_5;
}

inline LagrangePoints cr3bp_lagrange_points(double mu)
{
    LagrangePoints lp;
    const int MAX_ITER = 100;
    const double TOL = 1e-14;
    double mu1 = 1.0 - mu;

    // Hill sphere radius
    double rh = pow(mu / 3.0, 1.0/3.0);

    // --- L1: between m1 and m2, x in (-mu, 1-mu) ---
    // Initial guess: x_L1 = 1 - mu - gamma1
    // gamma1 ≈ rh * (1 - rh/3 - rh^2/9)
    double gamma1 = rh * (1.0 - rh/3.0 - rh*rh/9.0);
    double x = 1.0 - mu - gamma1;
    for (int i = 0; i < MAX_ITER; i++) {
        double f = lagrange_dOmega_dx(x, mu);
        double df = lagrange_d2Omega_dx2(x, mu);
        double dx = -f / df;
        x += dx;
        if (fabs(dx) < TOL) break;
    }
    lp.L1[0] = x; lp.L1[1] = 0; lp.L1[2] = 0;

    // --- L2: beyond m2 (away from m1), x > 1-mu ---
    // gamma2 ≈ rh * (1 + rh/3 - rh^2/9)
    double gamma2 = rh * (1.0 + rh/3.0 - rh*rh/9.0);
    x = 1.0 - mu + gamma2;
    for (int i = 0; i < MAX_ITER; i++) {
        double f = lagrange_dOmega_dx(x, mu);
        double df = lagrange_d2Omega_dx2(x, mu);
        double dx = -f / df;
        x += dx;
        if (fabs(dx) < TOL) break;
    }
    lp.L2[0] = x; lp.L2[1] = 0; lp.L2[2] = 0;

    // --- L3: opposite side of m1, x < -mu ---
    // Initial guess from Szebehely: x_L3 ≈ -1 - 5*mu/12
    x = -1.0 - 5.0*mu/12.0;
    for (int i = 0; i < MAX_ITER; i++) {
        double f = lagrange_dOmega_dx(x, mu);
        double df = lagrange_d2Omega_dx2(x, mu);
        double dx = -f / df;
        x += dx;
        if (fabs(dx) < TOL) break;
    }
    lp.L3[0] = x; lp.L3[1] = 0; lp.L3[2] = 0;

    // --- L4: equilateral triangle, upper (y > 0) ---
    lp.L4[0] = 0.5 - mu;
    lp.L4[1] = sqrt(3.0) / 2.0;
    lp.L4[2] = 0;

    // --- L5: equilateral triangle, lower (y < 0) ---
    lp.L5[0] = 0.5 - mu;
    lp.L5[1] = -sqrt(3.0) / 2.0;
    lp.L5[2] = 0;

    return lp;
}

///////////////////////////////////////////////////////////////////////////////
// Linearization at a collinear Lagrange point
//
// Eigenvalues: ±lambda (real, saddle), ±i*nu (imaginary, center)
// From KoLoMaRo Section 2.7
///////////////////////////////////////////////////////////////////////////////

struct CollinearLinearization
{
    double x_eq;       // equilibrium x-coordinate
    double mu_bar;     // mu|x_e-1+mu|^{-3} + (1-mu)|x_e+mu|^{-3}
    double a, b;       // a = 2*mu_bar + 1, b = mu_bar - 1
    double lambda;     // real eigenvalue (saddle)
    double nu;         // imaginary eigenvalue frequency (center)
};

inline CollinearLinearization cr3bp_linearize(double x_eq, double mu)
{
    CollinearLinearization lin;
    lin.x_eq = x_eq;

    double dx1 = fabs(x_eq + mu);
    double dx2 = fabs(x_eq - 1.0 + mu);
    double r1_3 = dx1*dx1*dx1;
    double r2_3 = dx2*dx2*dx2;

    lin.mu_bar = mu / r2_3 + (1.0-mu) / r1_3;
    lin.a = 2.0 * lin.mu_bar + 1.0;
    lin.b = lin.mu_bar - 1.0;

    // Characteristic polynomial: beta^4 + (2-mu_bar)*beta^2 + (1+mu_bar-2*mu_bar^2) = 0
    // Let alpha = beta^2:
    // alpha = (mu_bar - 2 ± sqrt(9*mu_bar^2 - 8*mu_bar)) / 2
    double disc = 9.0*lin.mu_bar*lin.mu_bar - 8.0*lin.mu_bar;
    double alpha1 = (lin.mu_bar - 2.0 + sqrt(disc)) / 2.0; // positive root
    double alpha2 = (lin.mu_bar - 2.0 - sqrt(disc)) / 2.0; // negative root

    lin.lambda = sqrt(alpha1);   // real eigenvalue (saddle direction)
    lin.nu = sqrt(-alpha2);      // imaginary eigenvalue frequency

    return lin;
}

///////////////////////////////////////////////////////////////////////////////
// Utility: Jacobi constant at a Lagrange point (energy level)
///////////////////////////////////////////////////////////////////////////////

inline double cr3bp_jacobi_at_point(double px, double py, double pz, double mu)
{
    // Zero velocity at equilibrium
    return 2.0 * cr3bp_omega(px, py, pz, mu);
}

#endif // CR3BP_HPP
