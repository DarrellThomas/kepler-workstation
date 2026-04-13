// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// EulerComponent.h — Rotational equations of motion
//
// Integrates rotational state (quaternion, angular velocity) using
// accumulated moments and moments of inertia.
//
// Euler's equations: I * alpha = M - omega × (I * omega)
// Quaternion kinematics: dq/dt = 0.5 * q * omega_quat
//
// Reference: Goldstein "Classical Mechanics", any rigid body dynamics text.

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Sim/SimTypes.h"
#include "EulerComponent.generated.h"

UCLASS(ClassGroup=(Simulation), meta=(BlueprintSpawnableComponent))
class CADACVIS_API UEulerComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UEulerComponent();

	/**
	 * Integrate one step: advance orientation and angular velocity.
	 * @param State      — current vehicle state (modified in place)
	 * @param Moments    — total moments in body frame (N·m)
	 * @param DeltaTime  — time step (s)
	 */
	UFUNCTION(BlueprintCallable, Category = "Euler")
	void Step(UPARAM(ref) FVehicleState& State,
	          const FVector& Moments,
	          double DeltaTime);
};
