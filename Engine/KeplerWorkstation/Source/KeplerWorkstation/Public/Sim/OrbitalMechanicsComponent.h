// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// OrbitalMechanicsComponent.h — N-body propagation and orbital tools
//
// Provides:
//   - N-body gravitational propagation (Sun + planets)
//   - JPL DE440 ephemeris access
//   - Lambert solver for transfer orbit design
//   - Orbital element computation (a, e, i, RAAN, AoP, TA)
//   - Sphere of influence detection
//   - Coordinate frame transforms (ECI ↔ heliocentric ↔ planetocentric)
//
// This is the backbone for interplanetary mission simulation.
// Used by the Bucket/Dipper constellation for VEEJ trajectory computation.
//
// Data: JPL DE440 ephemeris (1550-2650 CE), public domain.

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Sim/SimTypes.h"
#include "OrbitalMechanicsComponent.generated.h"

/** Classical orbital elements */
USTRUCT(BlueprintType)
struct CADACVIS_API FOrbitalElements
{
	GENERATED_BODY()

	UPROPERTY(BlueprintReadOnly) double SemiMajorAxis = 0;     // m
	UPROPERTY(BlueprintReadOnly) double Eccentricity = 0;
	UPROPERTY(BlueprintReadOnly) double Inclination = 0;       // rad
	UPROPERTY(BlueprintReadOnly) double RAAN = 0;              // rad (Right Ascension of Ascending Node)
	UPROPERTY(BlueprintReadOnly) double ArgPeriapsis = 0;      // rad
	UPROPERTY(BlueprintReadOnly) double TrueAnomaly = 0;       // rad
	UPROPERTY(BlueprintReadOnly) double Period = 0;            // seconds
	UPROPERTY(BlueprintReadOnly) double Apoapsis = 0;          // m (from center)
	UPROPERTY(BlueprintReadOnly) double Periapsis = 0;         // m (from center)
};

/** Result of a Lambert transfer computation */
USTRUCT(BlueprintType)
struct CADACVIS_API FLambertResult
{
	GENERATED_BODY()

	UPROPERTY(BlueprintReadOnly) FVector DepartureVelocity = FVector::ZeroVector; // km/s heliocentric
	UPROPERTY(BlueprintReadOnly) FVector ArrivalVelocity = FVector::ZeroVector;   // km/s heliocentric
	UPROPERTY(BlueprintReadOnly) double VInfDepart = 0;    // km/s
	UPROPERTY(BlueprintReadOnly) double VInfArrive = 0;    // km/s
	UPROPERTY(BlueprintReadOnly) double C3Depart = 0;      // km²/s²
	UPROPERTY(BlueprintReadOnly) double TOFDays = 0;
	UPROPERTY(BlueprintReadOnly) bool bValid = false;
};

UCLASS(ClassGroup=(Simulation), meta=(BlueprintSpawnableComponent))
class CADACVIS_API UOrbitalMechanicsComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UOrbitalMechanicsComponent();

	/** Load JPL DE440 ephemeris from file */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	bool LoadEphemeris(const FString& FilePath);

	/** Get planet position at Julian Date (km, SSB J2000) */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	FVector GetPlanetPosition(ECelestialBody Body, double JulianDate);

	/** Get planet velocity at Julian Date (km/s, SSB J2000) */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	FVector GetPlanetVelocity(ECelestialBody Body, double JulianDate);

	/** Solve Lambert's problem between two positions */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	FLambertResult SolveLambert(const FVector& R1, const FVector& R2,
	                            double TOFSeconds, double GM);

	/** Compute orbital elements from state vector */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	FOrbitalElements StateToElements(const FVector& Position, const FVector& Velocity, double GM);

	/** N-body propagate a state vector by DeltaTime seconds */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	void PropagateNBody(UPARAM(ref) FVector& Position, UPARAM(ref) FVector& Velocity,
	                    double JulianDate, double DeltaTime);

	/** Convert Julian Date to calendar year (decimal) */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	static double JDToYear(double JD);

	/** Convert calendar year to Julian Date */
	UFUNCTION(BlueprintCallable, Category = "Orbital")
	static double YearToJD(double Year);
};
