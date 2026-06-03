// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// Statistical Thermodynamics for H2/He Mixtures
//
// 7-species real-gas model from first principles:
//   H2, H, H+, He, He+, He2+, e-
//
// Computes equilibrium composition, thermodynamic properties, and mixture
// state for temperatures from 200 K to 200,000 K — covering cold gas through
// fully ionized plasma.
//
// Physics:
//   - Partition functions (translational, rotational, vibrational, electronic)
//   - H2 dissociation equilibrium (H2 ⇌ 2H)
//   - Saha ionization (H ⇌ H+ + e-, He ⇌ He+ + e-, He+ ⇌ He2+ + e-)
//   - Mixture-averaged h(T,p), Cp(T,p), γ_eff(T,p), speed of sound
//
// Data sources:
//   Molecular constants: NIST Chemistry WebBook
//   Ionization energies: NIST Atomic Spectra Database
//   Dissociation energy: Herzberg (1970), D0(H2) = 4.478 eV
//   Electronic degeneracies: atomic physics (ground state + first excited)
//
// References:
//   Vincenti & Kruger, "Introduction to Physical Gas Dynamics" (1965)
//   Park, "Nonequilibrium Hypersonic Aerothermodynamics" (1990)
//   Anderson, "Hypersonic and High-Temperature Gas Dynamics" (2006)
//
// Compile test: g++ -std=c++11 -O2 -o test_stat_thermo test_stat_thermo.cpp
///////////////////////////////////////////////////////////////////////////////

#ifndef STAT_THERMO_HPP
#define STAT_THERMO_HPP

#include <cmath>
#include <cstdio>
#include <algorithm>

namespace stat_thermo {

///////////////////////////////////////////////////////////////////////////////
// PHYSICAL CONSTANTS (SI)
///////////////////////////////////////////////////////////////////////////////

static const double KB      = 1.380649e-23;     // Boltzmann constant (J/K)
static const double NA      = 6.02214076e23;    // Avogadro's number (1/mol)
static const double RGAS    = 8.314462618;      // Universal gas constant (J/(mol·K))
static const double HPLANCK = 6.62607015e-34;   // Planck constant (J·s)
static const double ME      = 9.1093837015e-31; // Electron mass (kg)
static const double EV_TO_J = 1.602176634e-19;  // eV to Joules
static const double PI      = 3.14159265358979323846;

///////////////////////////////////////////////////////////////////////////////
// SPECIES DATA
///////////////////////////////////////////////////////////////////////////////

// Species indices
enum Species { H2 = 0, H = 1, Hp = 2, He = 3, Hep = 4, He2p = 5, Elec = 6 };
static const int NSPECIES = 7;

// Molecular weights (kg/mol)
static const double MW[NSPECIES] = {
    2.01588e-3,     // H2
    1.00794e-3,     // H
    1.00739e-3,     // H+ (proton, minus electron mass)
    4.00260e-3,     // He
    4.00205e-3,     // He+
    4.00150e-3,     // He2+
    5.48580e-7      // e-
};

// Particle masses (kg)
static const double MASS[NSPECIES] = {
    3.3474e-27,     // H2
    1.6737e-27,     // H
    1.6726e-27,     // H+
    6.6465e-27,     // He
    6.6454e-27,     // He+
    6.6443e-27,     // He2+
    9.1094e-31      // e-
};

// Electronic degeneracies (ground state statistical weight)
static const double GELEC[NSPECIES] = {
    1.0,    // H2 (Σ ground state)
    2.0,    // H  (2S1/2, g = 2)
    1.0,    // H+ (bare proton)
    1.0,    // He (1S0)
    2.0,    // He+ (2S1/2, hydrogen-like)
    1.0,    // He2+ (bare nucleus)
    2.0     // e- (spin degeneracy)
};

// H2 molecular constants
static const double THETA_ROT_H2  = 87.6;      // Rotational temperature (K)
static const double THETA_VIB_H2  = 6332.0;    // Vibrational temperature (K)

// Dissociation and ionization energies (J)
static const double D0_H2     = 4.478  * EV_TO_J;  // H2 → 2H
static const double I_H       = 13.598 * EV_TO_J;  // H  → H+ + e-
static const double I_HE1     = 24.587 * EV_TO_J;  // He → He+ + e-
static const double I_HE2     = 54.418 * EV_TO_J;  // He+ → He2+ + e-

// Dissociation/ionization energies per unit mass (J/kg) for enthalpy
static const double D0_H2_MASS = D0_H2 * NA / MW[H2];   // J/kg
static const double I_H_MASS   = I_H   * NA / MW[H];     // J/kg
static const double I_HE1_MASS = I_HE1 * NA / MW[He];    // J/kg
static const double I_HE2_MASS = I_HE2 * NA / MW[Hep];   // J/kg


///////////////////////////////////////////////////////////////////////////////
// PER-SPECIES THERMODYNAMIC PROPERTIES
///////////////////////////////////////////////////////////////////////////////

// Translational Cv contribution: (3/2) R/M  (all species, per unit mass)
inline double cv_trans(int sp)
{
    return 1.5 * RGAS / MW[sp];
}

// Rotational Cv for H2 (fully excited: R/M for diatomic, with high-T limit)
// At T >> θ_rot (87.6 K), Q_rot ≈ T/(2θ_rot) → Cv_rot = R/M
// At T < θ_rot, quantum corrections needed (but we're always above 200 K)
inline double cv_rot_h2(double T)
{
    if (T < 50.0) return 0.0;  // frozen out
    // Full quantum sum is expensive; use Einstein function approximation
    // For H2, θ_rot = 87.6 K. At T > 200K, essentially fully excited.
    double x = THETA_ROT_H2 / T;
    if (x < 0.1) return RGAS / MW[H2];  // classical limit
    // Smooth transition
    double f = x * x * exp(x) / ((exp(x) - 1.0) * (exp(x) - 1.0));
    return RGAS / MW[H2] * f;
}

// Vibrational Cv for H2 (Einstein oscillator)
// Cv_vib = R/M * (θ_vib/T)^2 * exp(θ_vib/T) / (exp(θ_vib/T) - 1)^2
inline double cv_vib_h2(double T)
{
    double x = THETA_VIB_H2 / T;
    if (x > 50.0) return 0.0;  // frozen out (T << 6332 K)
    double ex = exp(x);
    return RGAS / MW[H2] * x * x * ex / ((ex - 1.0) * (ex - 1.0));
}

// Per-species Cv (J/(kg·K)) — internal modes only (add translational separately)
inline double cv_internal(int sp, double T)
{
    if (sp == H2)
        return cv_rot_h2(T) + cv_vib_h2(T);
    // All other species are monatomic (or electron): no internal modes
    return 0.0;
}

// Per-species Cp (J/(kg·K)) = Cv_trans + Cv_internal + R/M
inline double cp_species(int sp, double T)
{
    return cv_trans(sp) + cv_internal(sp, T) + RGAS / MW[sp];
}

// Per-species enthalpy h(T) (J/kg) relative to elements at 0 K
// h = ∫Cp dT + formation enthalpy
// Formation enthalpy accounts for dissociation/ionization energy
inline double enthalpy_species(int sp, double T)
{
    double h_trans = 2.5 * RGAS * T / MW[sp];  // (3/2 + 1) kT per particle

    switch (sp) {
        case H2: {
            // H2: translational + rotational + vibrational
            // Reference: H2 molecule, h_f = 0 at 0 K
            double h_rot = RGAS * T / MW[H2];  // classical limit (T >> 87.6 K)
            double x = THETA_VIB_H2 / T;
            double h_vib = (x > 50.0) ? 0.0 :
                RGAS * THETA_VIB_H2 / MW[H2] / (exp(x) - 1.0);
            // Zero-point vibration energy (included in D0 definition)
            return h_trans + h_rot + h_vib;
        }
        case H:
            // H atom: monatomic, h_f = D0(H2)/2 per atom
            return h_trans + 0.5 * D0_H2_MASS * MW[H2] / (2.0 * MW[H]);
        case Hp:
            // H+: monatomic, h_f = D0(H2)/2 + I_H per proton
            return h_trans + 0.5 * D0_H2_MASS * MW[H2] / (2.0 * MW[Hp])
                           + I_H_MASS * MW[H] / MW[Hp];
        case He:
            // He: monatomic noble gas, h_f = 0
            return h_trans;
        case Hep:
            // He+: h_f = I_He1
            return h_trans + I_HE1_MASS;
        case He2p:
            // He2+: h_f = I_He1 + I_He2
            return h_trans + I_HE1_MASS * MW[He] / MW[He2p]
                           + I_HE2_MASS * MW[Hep] / MW[He2p];
        case Elec:
            // electron: translational only
            return h_trans;
        default:
            return h_trans;
    }
}


///////////////////////////////////////////////////////////////////////////////
// EQUILIBRIUM COMPOSITION SOLVER
///////////////////////////////////////////////////////////////////////////////

// Saha function: n_ion * n_e / n_atom = S(T)
// S(T) = (2 * g_ion / g_atom) * (2π m_e kT / h²)^(3/2) * exp(-I / kT)
inline double saha_rhs(double T, double I_energy, double g_ion, double g_atom)
{
    double lambda_e = HPLANCK / sqrt(2.0 * PI * ME * KB * T);
    double ne_factor = 2.0 / (lambda_e * lambda_e * lambda_e);
    return (g_ion / g_atom) * ne_factor * exp(-I_energy / (KB * T));
}

// H2 dissociation equilibrium constant Kp(T) in Pa
// For H2 ⇌ 2H: Kp = p_H^2 / p_H2
// From partition functions:
//   Kp(T) = (kT) * (Q_H / V)^2 / (Q_H2 / V) * exp(-D0/kT)
inline double dissociation_kp(double T)
{
    // Translational partition function ratio (per unit volume)
    // Q_trans/V = (2π m kT / h²)^(3/2)
    double lambda_H2 = HPLANCK / sqrt(2.0 * PI * MASS[H2] * KB * T);
    double lambda_H  = HPLANCK / sqrt(2.0 * PI * MASS[H]  * KB * T);
    double qt_H2 = 1.0 / (lambda_H2 * lambda_H2 * lambda_H2);
    double qt_H  = 1.0 / (lambda_H  * lambda_H  * lambda_H);

    // Internal partition functions
    // H2: Q_rot * Q_vib * Q_elec
    double Q_rot_H2 = T / (2.0 * THETA_ROT_H2);  // high-T limit
    double x = THETA_VIB_H2 / T;
    double Q_vib_H2 = (x > 50.0) ? 1.0 : 1.0 / (1.0 - exp(-x));
    double Q_int_H2 = Q_rot_H2 * Q_vib_H2 * GELEC[H2];

    // H: Q_elec only (monatomic)
    double Q_int_H = GELEC[H];  // ground state only at T < 50000 K

    // Kp in pressure units (Pa)
    // Kp = kT * (qt_H * Q_int_H)^2 / (qt_H2 * Q_int_H2) * exp(-D0/kT)
    double Kp = KB * T * (qt_H * Q_int_H) * (qt_H * Q_int_H)
                / (qt_H2 * Q_int_H2) * exp(-D0_H2 / (KB * T));
    return Kp;
}


///////////////////////////////////////////////////////////////////////////////
// GAS MIXTURE STATE
///////////////////////////////////////////////////////////////////////////////

struct GasMixture {
    double T;               // Temperature (K)
    double p;               // Pressure (Pa)
    double rho;             // Density (kg/m³)
    double X[NSPECIES];     // Mole fractions
    double Y[NSPECIES];     // Mass fractions
    double h;               // Mixture specific enthalpy (J/kg)
    double Cp;              // Mixture specific heat at constant pressure (J/(kg·K))
    double Cv;              // Mixture specific heat at constant volume (J/(kg·K))
    double gamma_eff;       // Effective ratio of specific heats
    double a_sound;         // Equilibrium speed of sound (m/s)
    double M_mix;           // Mixture molecular weight (kg/mol)
    double R_mix;           // Mixture gas constant (J/(kg·K))
    double alpha_d;         // Dissociation fraction of H2
    double alpha_iH;        // Ionization fraction of H
    double alpha_iHe1;      // First ionization fraction of He
    double alpha_iHe2;      // Second ionization fraction of He+
    double n_e;             // Electron number density (1/m³)
};


///////////////////////////////////////////////////////////////////////////////
// EQUILIBRIUM STATE SOLVER
//
// Given T, p, and initial elemental composition (mass fractions of H2 and He
// in the cold gas), compute the equilibrium 7-species mixture state.
//
// Y_H2_cold: mass fraction of H2 in cold gas (e.g., 0.862 for Jupiter)
// Y_He_cold: mass fraction of He in cold gas (e.g., 0.138 for Jupiter)
///////////////////////////////////////////////////////////////////////////////

inline GasMixture equilibrium_state(double T, double p,
                                     double Y_H2_cold = 0.862,
                                     double Y_He_cold = 0.138)
{
    GasMixture mix;
    mix.T = T;
    mix.p = p;

    // Convert cold-gas mass fractions to mole fractions
    double n_H2_0 = Y_H2_cold / MW[H2];  // moles of H2 per kg of mixture
    double n_He_0 = Y_He_cold / MW[He];   // moles of He per kg of mixture

    // Total hydrogen nuclei (per kg cold gas): 2 * n_H2_0 (in mol/kg)
    double N_H_total = 2.0 * n_H2_0;     // mol H nuclei per kg
    double N_He_total = n_He_0;           // mol He nuclei per kg

    // ── Step 1: H2 Dissociation ──
    // H2 ⇌ 2H,  Kp = (p_H)^2 / p_H2 = Kp(T)
    double Kp_d = dissociation_kp(T);

    // Solve for dissociation fraction α_d
    // Let α_d = fraction of H2 that dissociates
    // In terms of partial pressures (ideal mixture):
    //   n_H2 = n_H2_0 * (1 - α_d)
    //   n_H  = 2 * n_H2_0 * α_d
    //   Kp = [x_H^2 / x_H2] * p = [4 α_d^2 / (1 - α_d)] * p / (1 + α_d + n_He_0/n_H2_0)
    //
    // For robustness, iterate:
    double alpha_d = 0.0;

    if (T > 500.0 && Kp_d > 1e-30) {
        // Solve: Kp = 4α²p / [(1-α)(1+α+r)]  for α ∈ [0,1)
        // where r = n_He_0 / n_H2_0
        double r = N_He_total / n_H2_0;

        // Good initial guess from low-dissociation limit:
        // When α << 1: Kp ≈ 4α²p/(1+r), so α ≈ sqrt(Kp*(1+r)/(4p))
        double alpha_guess = sqrt(Kp_d * (1.0 + r) / (4.0 * p));
        if (alpha_guess > 0.99) alpha_guess = 0.5;
        if (alpha_guess < 1e-15) alpha_guess = 1e-15;
        alpha_d = alpha_guess;

        for (int iter = 0; iter < 100; iter++) {
            double one_ma = 1.0 - alpha_d;
            if (one_ma < 1e-15) one_ma = 1e-15;
            double denom = one_ma * (1.0 + alpha_d + r);
            double f = 4.0 * alpha_d * alpha_d * p / denom - Kp_d;

            // df/dα via quotient rule on 4α²/(( 1-α)(1+α+r))
            double num = 4.0 * alpha_d * alpha_d;
            double dnum = 8.0 * alpha_d;
            double ddenom = -(1.0 + alpha_d + r) + one_ma;  // = -2α - r
            double df = p * (dnum * denom - num * ddenom) / (denom * denom);

            if (fabs(df) < 1e-50) break;
            double da = -f / df;

            // Damped Newton to stay in bounds
            if (alpha_d + da < 0) da = -alpha_d * 0.5;
            if (alpha_d + da > 1.0 - 1e-12) da = (1.0 - 1e-12 - alpha_d) * 0.5;

            alpha_d += da;
            alpha_d = std::max(1e-30, std::min(1.0 - 1e-12, alpha_d));
            if (fabs(da) / (alpha_d + 1e-30) < 1e-12) break;
        }

        // If Kp is huge, α should be near 1
        if (Kp_d / p > 1e6) alpha_d = 1.0 - 1e-10;
    }

    // Species amounts after dissociation (mol per kg cold gas)
    double n_H2 = n_H2_0 * (1.0 - alpha_d);
    double n_H  = 2.0 * n_H2_0 * alpha_d;
    double n_He = N_He_total;

    // ── Step 2: Ionization (Saha) ──
    // Solve for ionization fractions given T and total number density

    // Total moles per kg (before ionization)
    double n_total_pre = n_H2 + n_H + n_He;
    double M_mix_pre = 1.0 / n_total_pre;  // kg/mol (mixture)
    double rho_est = p * M_mix_pre / (RGAS * T);
    double n_density = rho_est * NA / M_mix_pre;  // particles per m³

    // Saha RHS values (number density form: n_ion * n_e / n_neutral)
    double S_H   = saha_rhs(T, I_H,   GELEC[Hp],   GELEC[H]);
    double S_He1 = saha_rhs(T, I_HE1, GELEC[Hep],  GELEC[He]);
    double S_He2 = saha_rhs(T, I_HE2, GELEC[He2p], GELEC[Hep]);

    // Ionization fractions
    double alpha_iH = 0.0, alpha_iHe1 = 0.0, alpha_iHe2 = 0.0;
    double n_e_density = 0.0;

    if (T > 3000.0) {
        // Iterative solve: coupled Saha equations with charge neutrality
        // n_e = n_H+ + n_He+ + 2*n_He2+
        //
        // For each species: α_i = S_i / (n_e + S_i)  [from Saha with n_e]

        // Number densities of neutral species (per m³)
        double nd_H  = n_H  * NA * rho_est;  // approximate
        double nd_He = n_He * NA * rho_est;

        // Initial guess for n_e
        n_e_density = 1.0;  // start small

        for (int iter = 0; iter < 100; iter++) {
            double ne_old = n_e_density;

            // H ionization: n_Hp * n_e / n_H = S_H
            // n_H_total_atoms = nd_H (from dissociation step)
            // n_H_neutral = nd_H / (1 + S_H/n_e)
            // n_Hp = nd_H * S_H / (n_e + S_H)
            if (nd_H > 0 && S_H > 0)
                alpha_iH = S_H / (n_e_density + S_H);
            else
                alpha_iH = 0.0;

            // He first ionization
            // n_He_total = nd_He
            // n_He_neutral = nd_He / (1 + S_He1/n_e * (1 + S_He2/n_e))
            // n_Hep = n_He_neutral * S_He1 / n_e
            // n_He2p = n_Hep * S_He2 / n_e
            if (nd_He > 0 && S_He1 > 0) {
                double ratio1 = S_He1 / (n_e_density > 0 ? n_e_density : 1.0);
                double ratio2 = S_He2 / (n_e_density > 0 ? n_e_density : 1.0);
                double denom_He = 1.0 + ratio1 * (1.0 + ratio2);
                double f_He0 = 1.0 / denom_He;
                double f_Hep = ratio1 / denom_He;
                double f_He2p = ratio1 * ratio2 / denom_He;
                alpha_iHe1 = f_Hep + f_He2p;  // fraction of He that's ionized at all
                alpha_iHe2 = (f_Hep > 0) ? f_He2p / (f_Hep + f_He2p + 1e-40) : 0.0;

                // Electron contribution from He
                double ne_He = nd_He * (f_Hep + 2.0 * f_He2p);
                double ne_H = nd_H * alpha_iH;
                n_e_density = ne_H + ne_He;
            } else {
                n_e_density = nd_H * alpha_iH;
            }

            if (n_e_density < 1.0) n_e_density = 1.0;

            // Check convergence
            if (fabs(n_e_density - ne_old) / (ne_old + 1.0) < 1e-10)
                break;
        }
    }

    // ── Step 3: Final composition ──
    // Moles per kg cold gas
    double n_Hp = n_H * alpha_iH;
    double n_H_final = n_H * (1.0 - alpha_iH);

    // He ionization fractions (of total He)
    double f_He0, f_Hep_mol, f_He2p_mol;
    if (T > 3000.0 && S_He1 > 0) {
        double ne = (n_e_density > 0) ? n_e_density : 1.0;
        double r1 = S_He1 / ne;
        double r2 = S_He2 / ne;
        double dd = 1.0 + r1 * (1.0 + r2);
        f_He0 = 1.0 / dd;
        f_Hep_mol = r1 / dd;
        f_He2p_mol = r1 * r2 / dd;
    } else {
        f_He0 = 1.0;
        f_Hep_mol = 0.0;
        f_He2p_mol = 0.0;
    }

    double n_He_final = n_He * f_He0;
    double n_Hep_final = n_He * f_Hep_mol;
    double n_He2p_final = n_He * f_He2p_mol;
    double n_elec = n_Hp + n_Hep_final + 2.0 * n_He2p_final;

    // Total moles per kg cold gas
    double n_tot = n_H2 + n_H_final + n_Hp + n_He_final
                 + n_Hep_final + n_He2p_final + n_elec;

    // Mole fractions
    mix.X[H2]   = n_H2 / n_tot;
    mix.X[H]    = n_H_final / n_tot;
    mix.X[Hp]   = n_Hp / n_tot;
    mix.X[He]   = n_He_final / n_tot;
    mix.X[Hep]  = n_Hep_final / n_tot;
    mix.X[He2p] = n_He2p_final / n_tot;
    mix.X[Elec] = n_elec / n_tot;

    // Mixture molecular weight
    mix.M_mix = 0.0;
    for (int i = 0; i < NSPECIES; i++)
        mix.M_mix += mix.X[i] * MW[i];

    // Mass fractions
    for (int i = 0; i < NSPECIES; i++)
        mix.Y[i] = mix.X[i] * MW[i] / mix.M_mix;

    // Mixture gas constant
    mix.R_mix = RGAS / mix.M_mix;

    // Density from ideal gas (adequate for these pressures)
    mix.rho = p / (mix.R_mix * T);

    // ── Step 4: Mixture thermodynamic properties ──

    // Frozen specific heat (no reaction contribution)
    double cp_frozen = 0.0;
    for (int i = 0; i < NSPECIES; i++)
        cp_frozen += mix.Y[i] * cp_species(i, T);

    // Equilibrium Cp includes reaction contributions
    // dh/dT|_p = cp_frozen + Σ (dY_i/dT * h_i)
    // Approximate equilibrium Cp by finite difference
    double dT = T * 1e-4;
    if (dT < 0.1) dT = 0.1;
    GasMixture mix_hi;
    // Forward difference: re-solve at T+dT (recursive call limited to property calc)
    {
        // Lightweight re-solve for composition at T+dT
        double Kp_hi = dissociation_kp(T + dT);
        double ad_hi = alpha_d;
        double r = N_He_total / n_H2_0;
        for (int iter = 0; iter < 30; iter++) {
            double denom = (1.0 - ad_hi) * (1.0 + ad_hi + r);
            if (denom < 1e-30) { ad_hi = 1.0 - 1e-10; break; }
            double f = 4.0 * ad_hi * ad_hi * p / denom - Kp_hi;
            double df = 8.0 * ad_hi * p / denom
                      + 4.0 * ad_hi * ad_hi * p *
                        ((1.0 + ad_hi + r) + (1.0 - ad_hi)) / (denom * denom);
            if (fabs(df) < 1e-40) break;
            ad_hi -= f / df;
            ad_hi = std::max(0.0, std::min(1.0 - 1e-12, ad_hi));
        }
        // Recompute moles at T+dT
        double nH2_hi = n_H2_0 * (1.0 - ad_hi);
        double nH_hi = 2.0 * n_H2_0 * ad_hi;
        double ntot_hi = nH2_hi + nH_hi + n_He;  // approximate (ignore ionization change)
        double Mmix_hi = 0.0;
        double Xhi[NSPECIES] = {};
        Xhi[H2] = nH2_hi / ntot_hi;
        Xhi[H] = nH_hi / ntot_hi;
        Xhi[He] = n_He / ntot_hi;
        for (int i = 0; i < NSPECIES; i++) Mmix_hi += Xhi[i] * MW[i];
        double Yhi[NSPECIES] = {};
        for (int i = 0; i < NSPECIES; i++) Yhi[i] = Xhi[i] * MW[i] / Mmix_hi;

        double h_hi = 0.0;
        for (int i = 0; i < NSPECIES; i++)
            h_hi += Yhi[i] * enthalpy_species(i, T + dT);

        double h_lo = 0.0;
        for (int i = 0; i < NSPECIES; i++)
            h_lo += mix.Y[i] * enthalpy_species(i, T);

        mix.Cp = (h_hi - h_lo) / dT;
    }

    // If equilibrium Cp is unreasonable, fall back to frozen
    if (mix.Cp < cp_frozen * 0.5 || mix.Cp > cp_frozen * 100.0)
        mix.Cp = cp_frozen;

    // Mixture enthalpy
    mix.h = 0.0;
    for (int i = 0; i < NSPECIES; i++)
        mix.h += mix.Y[i] * enthalpy_species(i, T);

    // Cv and gamma
    mix.Cv = mix.Cp - mix.R_mix;
    if (mix.Cv < 1.0) mix.Cv = mix.Cp * 0.6;  // safety
    mix.gamma_eff = mix.Cp / mix.Cv;

    // Equilibrium speed of sound
    mix.a_sound = sqrt(mix.gamma_eff * mix.R_mix * T);

    // Store fractions
    mix.alpha_d = alpha_d;
    mix.alpha_iH = alpha_iH;
    mix.alpha_iHe1 = alpha_iHe1;
    mix.alpha_iHe2 = alpha_iHe2;
    mix.n_e = n_e_density;

    return mix;
}


///////////////////////////////////////////////////////////////////////////////
// RANKINE-HUGONIOT WITH REAL GAS
//
// Given freestream conditions and velocity, compute post-shock equilibrium
// state iteratively (since γ depends on post-shock temperature).
///////////////////////////////////////////////////////////////////////////////

struct ShockState {
    GasMixture pre;         // Freestream state
    GasMixture post;        // Post-shock equilibrium state
    double v_pre;           // Freestream velocity (m/s)
    double v_post;          // Post-shock velocity (m/s)
    double Mach_pre;        // Freestream Mach number
    double density_ratio;   // ρ₂/ρ₁
    double pressure_ratio;  // p₂/p₁
    double h_stag;          // Stagnation enthalpy (J/kg)
    int iterations;         // Solver iterations
};

inline ShockState normal_shock(double v_inf, double T_inf, double p_inf,
                                double Y_H2 = 0.862, double Y_He = 0.138)
{
    ShockState shock;

    // Freestream state
    shock.pre = equilibrium_state(T_inf, p_inf, Y_H2, Y_He);
    shock.v_pre = v_inf;
    shock.Mach_pre = v_inf / shock.pre.a_sound;
    shock.h_stag = shock.pre.h + 0.5 * v_inf * v_inf;

    // Iterative solution for post-shock state
    // Start with perfect-gas estimate
    double g = shock.pre.gamma_eff;
    double M = shock.Mach_pre;
    double M2 = M * M;

    // Perfect gas initial guess
    double T2 = T_inf * (2.0*g*M2 - (g-1.0)) * ((g-1.0)*M2 + 2.0)
                / ((g+1.0)*(g+1.0)*M2);
    double p2 = p_inf * (2.0*g*M2 - (g-1.0)) / (g+1.0);
    double rho_ratio = (g+1.0)*M2 / ((g-1.0)*M2 + 2.0);

    // Iterate: use energy conservation with real-gas enthalpy
    // h1 + v1²/2 = h2 + v2²/2
    // ρ1 * v1 = ρ2 * v2  (mass conservation)
    // p1 + ρ1*v1² = p2 + ρ2*v2²  (momentum conservation)

    shock.iterations = 0;
    for (int iter = 0; iter < 50; iter++) {
        shock.iterations = iter + 1;

        GasMixture post = equilibrium_state(T2, p2, Y_H2, Y_He);

        // From mass conservation: v2 = ρ1*v1/ρ2
        double rho1 = shock.pre.rho;
        double rho2 = post.rho;
        double v2 = rho1 * v_inf / rho2;

        // Energy conservation: h1 + v1²/2 = h2 + v2²/2
        double h2_target = shock.h_stag - 0.5 * v2 * v2;

        // Momentum: p2 = p1 + ρ1*v1² - ρ2*v2²
        double p2_momentum = p_inf + rho1*v_inf*v_inf - rho2*v2*v2;

        // Temperature from enthalpy error
        double h_error = h2_target - post.h;
        double T2_new = T2 + h_error / post.Cp;

        // Pressure from momentum
        double p2_new = p2_momentum;

        // Damped update
        double relax = 0.5;
        double dT = relax * (T2_new - T2);
        double dp = relax * (p2_new - p2);

        T2 += dT;
        p2 += dp;

        if (T2 < T_inf) T2 = T_inf * 1.1;
        if (p2 < p_inf) p2 = p_inf * 1.1;

        if (fabs(dT) / T2 < 1e-8 && fabs(dp) / p2 < 1e-8)
            break;
    }

    shock.post = equilibrium_state(T2, p2, Y_H2, Y_He);
    shock.v_post = shock.pre.rho * v_inf / shock.post.rho;
    shock.density_ratio = shock.post.rho / shock.pre.rho;
    shock.pressure_ratio = p2 / p_inf;

    return shock;
}


///////////////////////////////////////////////////////////////////////////////
// CONVENIENCE: Print mixture state
///////////////////////////////////////////////////////////////////////////////

inline void print_mixture(const GasMixture &m, const char *label = "")
{
    if (label[0]) printf("=== %s ===\n", label);
    printf("  T = %.1f K,  p = %.3e Pa,  rho = %.4e kg/m³\n", m.T, m.p, m.rho);
    printf("  M_mix = %.4f g/mol,  R_mix = %.1f J/(kg·K)\n", m.M_mix*1000, m.R_mix);
    printf("  h = %.3e J/kg,  Cp = %.1f J/(kg·K),  γ = %.4f\n", m.h, m.Cp, m.gamma_eff);
    printf("  a = %.1f m/s\n", m.a_sound);
    printf("  Composition (mole fractions):\n");
    const char *names[] = {"H2", "H", "H+", "He", "He+", "He2+", "e-"};
    for (int i = 0; i < NSPECIES; i++) {
        if (m.X[i] > 1e-10)
            printf("    %-6s %.6e\n", names[i], m.X[i]);
    }
    printf("  α_dissoc = %.6f,  α_ion(H) = %.6f\n", m.alpha_d, m.alpha_iH);
    printf("  α_ion(He1) = %.6f,  α_ion(He2) = %.6f\n", m.alpha_iHe1, m.alpha_iHe2);
    printf("  n_e = %.3e m⁻³\n", m.n_e);
}

inline void print_shock(const ShockState &s)
{
    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║  NORMAL SHOCK — REAL GAS (7-species H2/He)      ║\n");
    printf("╚══════════════════════════════════════════════════╝\n\n");
    printf("Freestream: v = %.1f m/s (%.2f km/s), Mach = %.1f\n",
           s.v_pre, s.v_pre/1000, s.Mach_pre);
    print_mixture(s.pre, "PRE-SHOCK");
    printf("\n");
    print_mixture(s.post, "POST-SHOCK (equilibrium)");
    printf("\n  Density ratio:  %.2f  (perfect gas limit: %.2f)\n",
           s.density_ratio,
           (s.pre.gamma_eff+1)*s.Mach_pre*s.Mach_pre /
           ((s.pre.gamma_eff-1)*s.Mach_pre*s.Mach_pre + 2));
    printf("  Pressure ratio: %.1f\n", s.pressure_ratio);
    printf("  h_stagnation:   %.3e J/kg\n", s.h_stag);
    printf("  Post-shock v:   %.1f m/s\n", s.v_post);
    printf("  Converged in %d iterations\n", s.iterations);
}

} // namespace stat_thermo

#endif // STAT_THERMO_HPP
