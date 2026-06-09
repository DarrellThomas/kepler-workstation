// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// TitanAtmosphereConfig.h — Atmosphere visual configuration for Titan
//
// Provides parameter structs for configuring UE5's SkyAtmosphere,
// ExponentialHeightFog, DirectionalLight, and post-processing to
// accurately represent Titan's orange tholin haze environment.
//
// Visual reference: Huygens DISR POSTER_F (surface panorama, Jan 14 2005)
//
// Key atmospheric properties:
//   - 98.4% N2 + 1.6% CH4, 1.47 bar surface pressure
//   - Tholin haze from UV photochemistry in upper atmosphere
//   - Orange color from tholin Mie scattering + methane blue absorption
//   - Surface visibility: 1-3 km (Huygens descent imagery)
//   - Sunlight at surface: 0.1% of Earth (perpetual deep twilight)
//   - No stars visible through haze
//   - Saturn visible as ghostly glow (5.2° angular diameter)

#pragma once

#include "CoreMinimal.h"
#include "TitanAtmosphereConfig.generated.h"

/**
 * Configuration for Titan's atmospheric visual appearance in UE5.
 * Apply these values to SkyAtmosphere, HeightFog, and lighting actors.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FTitanAtmosphereParams
{
	GENERATED_BODY()

	// ── SkyAtmosphere ──

	/** Titan's radius (km) — used as SkyAtmosphere BottomRadius */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere")
	float PlanetRadiusKm = 2575.0f;

	/** Height of modeled atmosphere (km) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere")
	float AtmosphereHeightKm = 600.0f;

	// Rayleigh scattering (N2 molecular scattering)
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Rayleigh")
	float RayleighScatteringScale = 1.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Rayleigh")
	FLinearColor RayleighScattering = FLinearColor(0.02f, 0.04f, 0.08f, 1.0f);

	/** Pressure scale height (km) — Titan: R*T/g = 296.8*93.5/1352 ≈ 20.5 km */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Rayleigh")
	float RayleighExponentialDistribution = 21.0f;

	// Mie scattering (tholin haze particles — dominates Titan's appearance)
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	float MieScatteringScale = 4.0f;

	/** Orange-tinted forward scattering from tholin aerosols */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	FLinearColor MieScattering = FLinearColor(0.8f, 0.4f, 0.1f, 1.0f);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	float MieAbsorptionScale = 3.0f;

	/** Tholins absorb blue/UV strongly */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	FLinearColor MieAbsorption = FLinearColor(0.1f, 0.4f, 0.9f, 1.0f);

	/** Strong forward scattering (typical for ~1μm tholin particles) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	float MieAnisotropy = 0.76f;

	/** Haze concentrated in lowest few km */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere|Mie")
	float MieExponentialDistribution = 5.0f;

	/** Multi-scattering factor — critical in thick atmosphere */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere")
	float MultiScatteringFactor = 2.0f;

	/** Ground albedo: dark tholin-covered surface */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "SkyAtmosphere")
	FColor GroundAlbedo = FColor(89, 46, 20);

	// ── Exponential Height Fog ──

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	float FogDensity = 0.015f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	float FogHeightFalloff = 0.001f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	FLinearColor FogInscatteringColor = FLinearColor(0.6f, 0.35f, 0.15f, 1.0f);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	float FogMaxOpacity = 0.95f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	bool bVolumetricFog = true;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	float VolumetricFogScatteringDistribution = 0.7f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Fog")
	FColor VolumetricFogAlbedo = FColor(180, 110, 40);

	// ── Directional Light (Sun) ──

	/** Titan receives ~0.1% of Earth's surface sunlight */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Lighting")
	float SunIntensityLux = 0.5f;

	/** Sun color: orange-filtered through tholin haze */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Lighting")
	FLinearColor SunColor = FLinearColor(1.0f, 0.7f, 0.3f, 1.0f);

	/** Sun disk angular size (degrees) — tiny at 9.54 AU */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Lighting")
	float SunSourceAngle = 0.06f;

	// ── Saturn Sky Element ──

	/** Saturn's angular diameter from Titan (degrees) — 10× the Moon from Earth */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Saturn")
	float SaturnAngularDiameterDeg = 5.2f;

	/** Saturn elevation above horizon (degrees) — depends on observer longitude */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Saturn")
	float SaturnElevationDeg = 30.0f;

	/** Saturn azimuth (degrees from north) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Saturn")
	float SaturnAzimuthDeg = 0.0f;

	// ── Post-Processing ──

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "PostProcess")
	float AutoExposureMinBrightness = 0.1f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "PostProcess")
	float AutoExposureMaxBrightness = 2.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "PostProcess")
	float BloomIntensity = 1.5f;

	/** Orange color grading — global saturation toward warm tones */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "PostProcess")
	FVector4 GlobalSaturation = FVector4(1.0f, 0.8f, 0.5f, 1.0f);

	// ── Audio ──

	/** Speed of sound on Titan (m/s) */
	static constexpr float SpeedOfSoundMS = 197.0f;

	/** Pitch shift factor: Titan / Earth speed of sound */
	static constexpr float AudioPitchShift = 197.0f / 343.0f; // 0.574

	/** Sound attenuation distance multiplier (dense air = carries further) */
	static constexpr float SoundAttenuationScale = 2.5f;
};

/**
 * Titan day/night cycle parameters.
 * Titan's day = 15.95 Earth days (tidally locked to Saturn).
 * In practice: very gradual lighting changes.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FTitanDayNightParams
{
	GENERATED_BODY()

	/** Titan's rotational period (seconds) */
	static constexpr double TitanDaySeconds = 15.95 * 86400.0; // 1,378,080 seconds

	/** In-game time multiplier (1.0 = real-time, 60 = 1 Titan day per ~6.4 hours) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "DayNight")
	float TimeScale = 60.0f;

	/** Current solar elevation angle (radians) — drives light intensity */
	UPROPERTY(BlueprintReadOnly, Category = "DayNight")
	float SolarElevation = 0.3f;
};
