// Copyright (c) 2026 Darrell Thomas. MIT License.
#include "Sim/EnvironmentComponent.h"
#include <cmath>

UEnvironmentComponent::UEnvironmentComponent() { PrimaryComponentTick.bCanEverTick = false; }

// Planetary constants: GM (m³/s²) and radius (m)
// Order must match ECelestialBody enum
static const double GM_TABLE[] = {
	3.986004418e14,  // Earth
	4.9028e12,       // Moon
	3.24859e14,      // Venus
	4.282837e13,     // Mars
	1.26686534e17,   // Jupiter
	3.7931187e16,    // Saturn
	8.978138e12,     // Titan
	6.835100e15,     // Neptune
	5.793951e15,     // Uranus
	1.32712440018e20 // Sun
};
static const double RADIUS_TABLE[] = {
	6378137.0,       // Earth
	1738100.0,       // Moon
	6051800.0,       // Venus
	3396200.0,       // Mars
	71492000.0,      // Jupiter
	60268000.0,      // Saturn
	2574700.0,       // Titan
	24764000.0,      // Neptune
	25559000.0,      // Uranus
	695700000.0      // Sun
};
static const double G0_TABLE[] = {
	9.80665,         // Earth
	1.62,            // Moon
	8.87,            // Venus
	3.71,            // Mars
	24.79,           // Jupiter
	10.44,           // Saturn
	1.352,           // Titan
	11.15,           // Neptune
	8.87,            // Uranus
	274.0            // Sun
};

// ─────────────────────────────────────────────────────────
// Titan atmosphere: data-driven model from Huygens HASI
// 273-point lookup table from actual descent measurements
// Source: NASA PDS HP-SSA-HASI-2-3-4-MISSION-V1.1
// Reference: Fulchignoni et al. 2005, Nature 438:785
// ─────────────────────────────────────────────────────────

// Include the full HASI lookup table
// This header provides atmosphere_titan_hasi() inline function
// with binary search + log-linear interpolation
#include "../../../../PlanetaryModels/titan_hasi_table.hpp"

// Titan gas properties
static const double TITAN_RGAS = 296.8;   // J/(kg·K) — N2-dominated
static const double TITAN_GAMMA = 1.40;   // ratio of specific heats

double UEnvironmentComponent::GetGM() const { return GM_TABLE[(int)Body]; }
double UEnvironmentComponent::GetSurfaceRadius() const { return RADIUS_TABLE[(int)Body]; }

double UEnvironmentComponent::GetDensityAtAltitude(double AltitudeM) const
{
	if (Body == ECelestialBody::Titan)
	{
		double rho, press, tempk;
		atmosphere_titan_hasi(rho, press, tempk, AltitudeM);
		return rho;
	}

	// Placeholder for other bodies — full implementation will call planetary_environment models
	return 0.0;
}

FVector UEnvironmentComponent::Evaluate(const FVector& Position)
{
	double R = Position.Size();
	double Alt = R - GetSurfaceRadius();

	// Gravity: point mass (harmonics in full implementation)
	double GM = GetGM();
	if (R > 1.0)
	{
		FVector Grav = -Position * (GM / (R * R * R));
		Gravity.Acceleration = Grav;
		Gravity.Magnitude = Grav.Size();
	}

	// Atmosphere
	if (Body == ECelestialBody::Titan)
	{
		double rho, press, tempk;
		atmosphere_titan_hasi(rho, press, tempk, Alt);
		Atmosphere.Density = rho;
		Atmosphere.Pressure = press;
		Atmosphere.Temperature = tempk;
		Atmosphere.SpeedOfSound = std::sqrt(TITAN_GAMMA * TITAN_RGAS * tempk);
		// Wind: simplified from DWE data — nearly calm at surface
		// Full implementation will interpolate ZONALWIND.TAB
		Atmosphere.Wind = FVector(0.2, 0, 0); // ~0.2 m/s eastward at surface
	}
	else
	{
		Atmosphere.Density = GetDensityAtAltitude(Alt);
	}

	return Gravity.Acceleration;
}
