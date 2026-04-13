// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.
//
// SimTypes.h — Core types for aerospace vehicle simulation
//
// Cleanroom implementation of 6DOF aerospace simulation framework.
// No third-party simulation code. Physics from published textbooks:
//   - Stevens & Lewis, "Aircraft Control and Simulation" (3rd ed)
//   - Zipfel, "Modeling and Simulation of Aerospace Vehicle Dynamics" (theory only)
//   - Bate, Mueller & White, "Fundamentals of Astrodynamics"
//   - Vallado, "Fundamentals of Astrodynamics and Applications"
//
// All simulation state is expressed in SI units:
//   Position: meters, Velocity: m/s, Angles: radians (internal), Force: Newtons
//   Mass: kg, Time: seconds
//
// Coordinate frames:
//   ECI  — Earth-Centered Inertial (J2000)
//   ECEF — Earth-Centered Earth-Fixed
//   NED  — North-East-Down (local tangent plane)
//   Body — Vehicle body frame (X-forward, Y-right, Z-down)
//   SSB  — Solar System Barycentric (for interplanetary)

#pragma once

#include "CoreMinimal.h"
#include "SimTypes.generated.h"

/**
 * Complete state vector for a rigid body vehicle.
 * Contains translational + rotational state, mass properties,
 * and derived quantities updated each integration step.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FVehicleState
{
	GENERATED_BODY()

	// --- Translational state ---

	/** Position in reference frame (m) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) FVector Position = FVector::ZeroVector;

	/** Velocity in reference frame (m/s) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) FVector Velocity = FVector::ZeroVector;

	/** Total acceleration (m/s²) — computed, not integrated */
	UPROPERTY(BlueprintReadOnly) FVector Acceleration = FVector::ZeroVector;

	// --- Rotational state ---

	/** Orientation quaternion (body wrt reference frame) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) FQuat Orientation = FQuat::Identity;

	/** Angular velocity in body frame (rad/s) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) FVector AngularVelocity = FVector::ZeroVector;

	/** Angular acceleration in body frame (rad/s²) — computed */
	UPROPERTY(BlueprintReadOnly) FVector AngularAcceleration = FVector::ZeroVector;

	// --- Euler angles (derived from quaternion, for display/output) ---

	/** Yaw/Heading (rad) */
	UPROPERTY(BlueprintReadOnly) double Yaw = 0;

	/** Pitch (rad) */
	UPROPERTY(BlueprintReadOnly) double Pitch = 0;

	/** Roll (rad) */
	UPROPERTY(BlueprintReadOnly) double Roll = 0;

	// --- Mass properties ---

	/** Total mass (kg) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) double Mass = 1.0;

	/** Moments of inertia about body axes (kg·m²) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite) FVector MoI = FVector(1, 1, 1);

	// --- Derived / computed ---

	/** Speed (m/s) */
	UPROPERTY(BlueprintReadOnly) double Speed = 0;

	/** Mach number */
	UPROPERTY(BlueprintReadOnly) double Mach = 0;

	/** Altitude above reference surface (m) */
	UPROPERTY(BlueprintReadOnly) double Altitude = 0;

	/** Dynamic pressure (Pa) */
	UPROPERTY(BlueprintReadOnly) double DynamicPressure = 0;

	/** Angle of attack (rad) */
	UPROPERTY(BlueprintReadOnly) double Alpha = 0;

	/** Sideslip angle (rad) */
	UPROPERTY(BlueprintReadOnly) double Beta = 0;

	/** Total load factor (g) */
	UPROPERTY(BlueprintReadOnly) double LoadFactor = 0;

	/** Simulation time (s) */
	UPROPERTY(BlueprintReadOnly) double SimTime = 0;

	/** Update derived quantities from primary state */
	void UpdateDerived(double AtmoDensity, double SpeedOfSound, double SurfaceRadius);
};

/**
 * Atmospheric conditions at the vehicle's current position.
 * Computed by the environment component each step.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FAtmoState
{
	GENERATED_BODY()

	UPROPERTY(BlueprintReadOnly) double Density = 0;        // kg/m³
	UPROPERTY(BlueprintReadOnly) double Pressure = 0;       // Pa
	UPROPERTY(BlueprintReadOnly) double Temperature = 0;    // K
	UPROPERTY(BlueprintReadOnly) double SpeedOfSound = 0;   // m/s
	UPROPERTY(BlueprintReadOnly) FVector Wind = FVector::ZeroVector; // m/s in NED
};

/**
 * Gravitational field at the vehicle's current position.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FGravityState
{
	GENERATED_BODY()

	UPROPERTY(BlueprintReadOnly) FVector Acceleration = FVector(0, 0, 0); // m/s²
	UPROPERTY(BlueprintReadOnly) double Magnitude = 0; // m/s²
};

/**
 * Forces and moments acting on the vehicle (body frame).
 * Each physics component contributes to these, then Newton/Euler integrates.
 */
USTRUCT(BlueprintType)
struct CADACVIS_API FForcesMoments
{
	GENERATED_BODY()

	UPROPERTY(BlueprintReadWrite) FVector Force = FVector::ZeroVector;   // N (body frame)
	UPROPERTY(BlueprintReadWrite) FVector Moment = FVector::ZeroVector;  // N·m (body frame)

	void Reset() { Force = FVector::ZeroVector; Moment = FVector::ZeroVector; }
	void Add(const FForcesMoments& Other) { Force += Other.Force; Moment += Other.Moment; }
};

/**
 * Celestial body identifier for multi-body simulations.
 */
UENUM(BlueprintType)
enum class ECelestialBody : uint8
{
	Earth    UMETA(DisplayName = "Earth"),
	Moon     UMETA(DisplayName = "Moon"),
	Venus    UMETA(DisplayName = "Venus"),
	Mars     UMETA(DisplayName = "Mars"),
	Jupiter  UMETA(DisplayName = "Jupiter"),
	Saturn   UMETA(DisplayName = "Saturn"),
	Sun      UMETA(DisplayName = "Sun"),
};

/**
 * Reference frame for position/velocity.
 */
UENUM(BlueprintType)
enum class ESimFrame : uint8
{
	ECI      UMETA(DisplayName = "Earth-Centered Inertial"),
	ECEF     UMETA(DisplayName = "Earth-Centered Earth-Fixed"),
	NED      UMETA(DisplayName = "North-East-Down"),
	Body     UMETA(DisplayName = "Body"),
	SSB      UMETA(DisplayName = "Solar System Barycentric"),
};
