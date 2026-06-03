// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// TRANSPORT PROPERTIES FOR HIGH-TEMPERATURE GAS MIXTURES
//
// Viscosity, thermal conductivity, and binary diffusion coefficients
// for the 7-species H2/He system (H2, H, H+, He, He+, He2+, e-)
//
// Methods:
//   - Neutral-neutral: Chapman-Enskog with Lennard-Jones collision integrals
//   - Ion-neutral: Langevin polarization model
//   - Charged-charged: Coulomb collision integrals (screened by Debye length)
//   - Mixture averaging: Wilke's rule (viscosity), Eucken (conductivity)
//
// References:
//   Gupta, Yos, Thompson, Lee — NASA RP-1232 (1990)
//   Hirschfelder, Curtiss, Bird — "Molecular Theory of Gases & Liquids"
//   Capitelli et al. — "Transport Properties of High Temperature Air" (2000)
//   Svehla — NASA TR R-132 (1962) collision integral tables
///////////////////////////////////////////////////////////////////////////////

#ifndef TRANSPORT_PROPERTIES_HPP
#define TRANSPORT_PROPERTIES_HPP

#include "stat_thermo.hpp"
#include <cmath>

namespace transport {

using namespace stat_thermo;

///////////////////////////////////////////////////////////////////////////////
// LENNARD-JONES PARAMETERS (neutral species)
///////////////////////////////////////////////////////////////////////////////

// σ (collision diameter, Å) and ε/k (well depth, K)
// Source: Svehla NASA TR R-132, Hirschfelder et al.
struct LJParams { double sigma; double eps_k; };

// Index: H2=0, H=1, He=3 (matching stat_thermo species indices for neutrals)
static const LJParams LJ_H2  = {2.827, 59.7};
static const LJParams LJ_H   = {2.708, 37.0};   // from Capitelli
static const LJParams LJ_He  = {2.551, 10.22};

// Combining rules for unlike pairs (Lorentz-Berthelot)
inline LJParams lj_combine(const LJParams &a, const LJParams &b)
{
    return { 0.5 * (a.sigma + b.sigma), sqrt(a.eps_k * b.eps_k) };
}


///////////////////////////////////////////////////////////////////////////////
// COLLISION INTEGRALS Ω^(l,s)*
//
// Neufeld-Janzen-Aziz fit for Lennard-Jones potential:
//   Ω^(2,2)* = A/T*^B + C/exp(D*T*) + E/exp(F*T*) + G/exp(H*T*)
//   where T* = kT/ε
//
// Coefficients from Neufeld, Janzen, Aziz (1972) J. Chem. Phys. 57:1100
///////////////////////////////////////////////////////////////////////////////

// Ω^(1,1)* — for diffusion
inline double omega_11_star(double T_star)
{
    double A = 1.06036, B = 0.15610, C = 0.19300, D = 0.47635;
    double E = 1.03587, F = 1.52996, G = 1.76474, H = 3.89411;
    return A / pow(T_star, B) + C / exp(D * T_star)
         + E / exp(F * T_star) + G / exp(H * T_star);
}

// Ω^(2,2)* — for viscosity
inline double omega_22_star(double T_star)
{
    double A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320;
    double E = 2.16178, F = 2.43787;
    return A / pow(T_star, B) + C / exp(D * T_star) + E / exp(F * T_star);
}


///////////////////////////////////////////////////////////////////////////////
// COULOMB COLLISION INTEGRALS (charged species)
//
// For ion-ion and ion-electron interactions, screened Coulomb potential.
// Uses Debye length for screening.
//
// Ω_C^(l,s) = (4π/kT)^2 * (Z1*Z2*e²/(4πε0))^2 * ln(Λ) / (16πn)
// Simplified: use tabulated Coulomb logarithm
///////////////////////////////////////////////////////////////////////////////

inline double debye_length(double T, double n_e)
{
    if (n_e < 1.0) return 1e10;  // no plasma
    double eps0 = 8.854187817e-12;
    return sqrt(eps0 * KB * T / (n_e * EV_TO_J * EV_TO_J / (KB * T) * n_e));
    // Simplified: λ_D = sqrt(ε₀ kT / (n_e e²))
}

inline double coulomb_log(double T, double n_e)
{
    if (n_e < 1.0) return 10.0;
    // ln(Λ) = ln(12π n_e λ_D³)
    // Simplified: ln(Λ) ≈ 23 - ln(sqrt(n_e) / T^(3/2))  (for T in eV, n_e in cm⁻³)
    double T_eV = T * KB / EV_TO_J;
    double n_cgs = n_e * 1e-6;  // m⁻³ to cm⁻³
    double lnL = 23.0 - log(sqrt(n_cgs) / pow(T_eV, 1.5));
    if (lnL < 2.0) lnL = 2.0;
    if (lnL > 30.0) lnL = 30.0;
    return lnL;
}


///////////////////////////////////////////////////////////////////////////////
// SINGLE-SPECIES VISCOSITY
///////////////////////////////////////////////////////////////////////////////

// Chapman-Enskog viscosity for a pure gas:
// μ = 5/(16) * sqrt(π m kT) / (π σ² Ω^(2,2)*)
// Returns Pa·s
inline double viscosity_lj(double T, double mass_kg, const LJParams &lj)
{
    double T_star = T / lj.eps_k;
    double omega22 = omega_22_star(T_star);
    double sigma_m = lj.sigma * 1e-10;  // Å to m
    return 5.0 / 16.0 * sqrt(PI * mass_kg * KB * T)
           / (PI * sigma_m * sigma_m * omega22);
}

// Charged species viscosity (Coulomb collisions)
// μ_ion ≈ 0.96 * (kT)^(5/2) * sqrt(m_ion) / (Z^4 * e^4 * ln(Λ))
inline double viscosity_coulomb(double T, double mass_kg, int Z, double n_e)
{
    double lnL = coulomb_log(T, n_e);
    double e = 1.602176634e-19;
    double eps0 = 8.854187817e-12;
    // Spitzer viscosity
    double mu = 0.96 * pow(KB * T, 2.5) * sqrt(mass_kg)
                / (pow(Z, 4) * pow(e, 4) / (16.0 * PI * PI * eps0 * eps0) * lnL);
    return mu;
}


///////////////////////////////////////////////////////////////////////////////
// SINGLE-SPECIES THERMAL CONDUCTIVITY
///////////////////////////////////////////////////////////////////////////////

// Eucken correction for polyatomic molecules:
// k = μ * (f_trans * Cv_trans + f_int * Cv_int) / M
// where f_trans = 5/2, f_int = 1 (approximate)
//
// For monatomic: k = (15/4) * (R/M) * μ  (exact Chapman-Enskog)

inline double conductivity_mono(double mu, double M_kg_mol)
{
    return 15.0 / 4.0 * RGAS / M_kg_mol * mu;
}

inline double conductivity_diatomic(double T, double mu, double M_kg_mol)
{
    // Eucken: k = μ/M * (5/2 * Cv_trans + Cv_rot + Cv_vib)
    // Cv_trans = 3/2 R, Cv_rot ≈ R, Cv_vib = R*(θ/T)²*exp(θ/T)/(exp(θ/T)-1)²
    double Cv_t = 1.5 * RGAS;
    double Cv_r = RGAS;  // fully excited
    double x = THETA_VIB_H2 / T;
    double Cv_v = (x > 50.0) ? 0.0 : RGAS * x * x * exp(x) / pow(exp(x) - 1.0, 2);
    return mu / M_kg_mol * (2.5 * Cv_t + Cv_r + Cv_v);
}

// Electron thermal conductivity (Spitzer)
inline double conductivity_electron(double T, double n_e)
{
    double lnL = coulomb_log(T, n_e);
    // k_e = 3.2 * (kT)^(5/2) / (m_e^(1/2) * e^4 * lnΛ) * (4πε₀)²
    double e = 1.602176634e-19;
    double eps0 = 8.854187817e-12;
    return 3.2 * pow(KB * T, 2.5) / (sqrt(ME) * pow(e, 4)
           / (16 * PI * PI * eps0 * eps0) * lnL);
}


///////////////////////////////////////////////////////////////////////////////
// MIXTURE TRANSPORT PROPERTIES
///////////////////////////////////////////////////////////////////////////////

struct TransportProps {
    double mu;          // mixture viscosity (Pa·s)
    double k_thermal;   // mixture thermal conductivity (W/(m·K))
    double Pr;          // Prandtl number
};

inline TransportProps mixture_transport(const GasMixture &mix)
{
    TransportProps tp;

    // Per-species viscosities
    double mu_sp[NSPECIES];
    double k_sp[NSPECIES];

    for (int i = 0; i < NSPECIES; i++) {
        if (mix.X[i] < 1e-15) {
            mu_sp[i] = 0;
            k_sp[i] = 0;
            continue;
        }

        switch (i) {
            case H2:
                mu_sp[i] = viscosity_lj(mix.T, MASS[H2], LJ_H2);
                k_sp[i] = conductivity_diatomic(mix.T, mu_sp[i], MW[H2]);
                break;
            case H:
                mu_sp[i] = viscosity_lj(mix.T, MASS[H], LJ_H);
                k_sp[i] = conductivity_mono(mu_sp[i], MW[H]);
                break;
            case He:
                mu_sp[i] = viscosity_lj(mix.T, MASS[He], LJ_He);
                k_sp[i] = conductivity_mono(mu_sp[i], MW[He]);
                break;
            case Hp:
                mu_sp[i] = viscosity_coulomb(mix.T, MASS[Hp], 1, mix.n_e);
                k_sp[i] = conductivity_mono(mu_sp[i], MW[Hp]);
                break;
            case Hep:
                mu_sp[i] = viscosity_coulomb(mix.T, MASS[Hep], 1, mix.n_e);
                k_sp[i] = conductivity_mono(mu_sp[i], MW[Hep]);
                break;
            case He2p:
                mu_sp[i] = viscosity_coulomb(mix.T, MASS[He2p], 2, mix.n_e);
                k_sp[i] = conductivity_mono(mu_sp[i], MW[He2p]);
                break;
            case Elec:
                mu_sp[i] = 0;  // electron viscosity negligible
                k_sp[i] = conductivity_electron(mix.T, mix.n_e);
                break;
        }
    }

    // Wilke's mixing rule for viscosity
    // μ_mix = Σ (X_i * μ_i) / Σ (X_j * Φ_ij)
    tp.mu = 0;
    for (int i = 0; i < NSPECIES; i++) {
        if (mix.X[i] < 1e-15 || mu_sp[i] < 1e-30) continue;
        double denom = 0;
        for (int j = 0; j < NSPECIES; j++) {
            if (mix.X[j] < 1e-15) continue;
            if (mu_sp[j] < 1e-30) continue;
            double phi = pow(1.0 + sqrt(mu_sp[i] / mu_sp[j])
                        * pow(MW[j] / MW[i], 0.25), 2)
                        / sqrt(8.0 * (1.0 + MW[i] / MW[j]));
            denom += mix.X[j] * phi;
        }
        if (denom > 0)
            tp.mu += mix.X[i] * mu_sp[i] / denom;
    }

    // Same mixing rule for conductivity (approximate)
    tp.k_thermal = 0;
    for (int i = 0; i < NSPECIES; i++) {
        if (mix.X[i] < 1e-15 || k_sp[i] < 1e-30) continue;
        double denom = 0;
        for (int j = 0; j < NSPECIES; j++) {
            if (mix.X[j] < 1e-15) continue;
            if (mu_sp[j] < 1e-30 && k_sp[j] < 1e-30) continue;
            double mu_j = (mu_sp[j] > 0) ? mu_sp[j] : 1e-10;
            double mu_i = (mu_sp[i] > 0) ? mu_sp[i] : 1e-10;
            double phi = pow(1.0 + sqrt(mu_i / mu_j)
                        * pow(MW[j] / MW[i], 0.25), 2)
                        / sqrt(8.0 * (1.0 + MW[i] / MW[j]));
            denom += mix.X[j] * phi;
        }
        if (denom > 0)
            tp.k_thermal += mix.X[i] * k_sp[i] / denom;
    }

    // Prandtl number
    tp.Pr = (tp.k_thermal > 0) ? mix.Cp * tp.mu / tp.k_thermal : 0.7;

    return tp;
}


///////////////////////////////////////////////////////////////////////////////
// CONVENIENCE: Print transport properties
///////////////////////////////////////////////////////////////////////////////

inline void print_transport(const TransportProps &tp, const char *label = "")
{
    if (label[0]) printf("  %s:\n", label);
    printf("    μ = %.4e Pa·s,  k = %.4e W/(m·K),  Pr = %.3f\n",
           tp.mu, tp.k_thermal, tp.Pr);
}

} // namespace transport

#endif // TRANSPORT_PROPERTIES_HPP
