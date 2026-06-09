// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// EnvironmentComponent.h — Planetary atmosphere and gravity
//
// Provides atmospheric state (density, pressure, temperature, speed of sound)
// and gravitational acceleration at the vehicle's current position.
// Supports multiple celestial bodies (Earth, Venus, Mars, Jupiter, Saturn,
// Titan, Neptune, Uranus, Moon).
// Gravity includes zonal harmonics (J2-J10 depending on body).
//
// Data sources:
//   Earth atmo: US Standard Atmosphere 1976
//   Mars atmo:  NASA GRC (Mars Global Surveyor 1996)
//   Venus atmo: Venera/Pioneer Venus probe data
//   Jupiter atmo: Galileo Probe ASI (Seiff et al. 1998)
//   Titan atmo: Huygens HASI (Fulchignoni et al. 2005) — 273-point data-driven lookup
//   Earth gravity: WGS-84/EGM96, J2-J6
//   Mars gravity: GMM-3 (Genova et al. 2016), J2-J6
//   Venus gravity: Magellan MGNP180U (Konopliv et al. 1999), J2-J4
//   Jupiter gravity: Juno (Iess et al. 2018), J2-J10
//   Titan gravity: Cassini (Iess et al. 2010), J2
//   Moon gravity: GRAIL (Zuber/Konopliv 2013), J2-J6

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Sim/SimTypes.h"
#include "EnvironmentComponent.generated.h"

UCLASS(ClassGroup=(Simulation), meta=(BlueprintSpawnableComponent))
class CADACVIS_API UEnvironmentComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UEnvironmentComponent();

	/** Which celestial body's environment to use */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Environment")
	ECelestialBody Body = ECelestialBody::Earth;

	/** Enable non-standard atmosphere (weather deck) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Environment")
	bool bUseCustomAtmosphere = false;

	/** Current atmospheric state (updated each call to Evaluate) */
	UPROPERTY(BlueprintReadOnly, Category = "Environment")
	FAtmoState Atmosphere;

	/** Current gravity (updated each call to Evaluate) */
	UPROPERTY(BlueprintReadOnly, Category = "Environment")
	FGravityState Gravity;

	/**
	 * Evaluate atmosphere and gravity at a given position.
	 * Call this once per simulation step.
	 * @param Position — vehicle position in planet-centered frame (m)
	 * @return Gravity acceleration vector (m/s²)
	 */
	UFUNCTION(BlueprintCallable, Category = "Environment")
	FVector Evaluate(const FVector& Position);

	/** Get atmospheric density at arbitrary altitude (m) */
	UFUNCTION(BlueprintCallable, Category = "Environment")
	double GetDensityAtAltitude(double AltitudeM) const;

	/** Get surface radius for current body (m) */
	UFUNCTION(BlueprintCallable, Category = "Environment")
	double GetSurfaceRadius() const;

	/** Get GM for current body (m³/s²) */
	UFUNCTION(BlueprintCallable, Category = "Environment")
	double GetGM() const;
};
