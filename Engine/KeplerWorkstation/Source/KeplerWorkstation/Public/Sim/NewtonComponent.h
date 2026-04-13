// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// NewtonComponent.h — Translational equations of motion (F = ma)
//
// Integrates translational state (position, velocity) using
// accumulated forces, gravity, and vehicle mass.
//
// Integration: RK4 by default, Euler-midpoint option for speed.
// Supports flat-earth (NED) and round-earth (ECI/ECEF) frames.
//
// Reference: Any introductory dynamics text. F = ma is not copyrightable.

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Sim/SimTypes.h"
#include "NewtonComponent.generated.h"

UENUM(BlueprintType)
enum class EIntegrationMethod : uint8
{
	Euler        UMETA(DisplayName = "Euler (1st order)"),
	EulerMidpoint UMETA(DisplayName = "Euler-Midpoint (2nd order)"),
	RK4          UMETA(DisplayName = "Runge-Kutta 4 (4th order)"),
};

UCLASS(ClassGroup=(Simulation), meta=(BlueprintSpawnableComponent))
class CADACVIS_API UNewtonComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UNewtonComponent();

	/** Integration method */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Integration")
	EIntegrationMethod Method = EIntegrationMethod::RK4;

	/** Integration time step (seconds) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Integration")
	double TimeStep = 0.001;

	/**
	 * Integrate one step: advance position and velocity.
	 * @param State      — current vehicle state (modified in place)
	 * @param Forces     — total forces in body frame (N)
	 * @param Gravity    — gravitational acceleration in reference frame (m/s²)
	 * @param DeltaTime  — time step (s), overrides TimeStep if > 0
	 */
	UFUNCTION(BlueprintCallable, Category = "Newton")
	void Step(UPARAM(ref) FVehicleState& State,
	          const FForcesMoments& Forces,
	          const FVector& GravityAccel,
	          double DeltaTime = 0);
};
