// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// TitanPlayerCharacter.h — First-person character for Titan surface exploration
//
// Physics:
//   Gravity: 1.352 m/s² (0.138 × Earth)
//   Surface air density: 5.34 kg/m³ (4.4 × Earth)
//   Terminal velocity: 7.6 m/s (can't die from falling)
//   Speed of sound: 197 m/s
//
// Data source: Huygens HASI (Fulchignoni et al. 2005)

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Character.h"
#include "Sim/EnvironmentComponent.h"
#include "TitanPlayerCharacter.generated.h"

UCLASS()
class CADACVIS_API ATitanPlayerCharacter : public ACharacter
{
	GENERATED_BODY()

public:
	ATitanPlayerCharacter();

	virtual void BeginPlay() override;
	virtual void Tick(float DeltaTime) override;
	virtual void SetupPlayerInputComponent(class UInputComponent* PlayerInputComponent) override;

	// ── Environment ──

	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Titan")
	UEnvironmentComponent* Environment;

	// ── Titan Physics Constants ──

	/** Surface gravity (m/s²) */
	static constexpr double TitanGravity = 1.352;

	/** Earth gravity (m/s²) */
	static constexpr double EarthGravity = 9.80665;

	/** Gravity scale for UE5 CharacterMovement (Titan/Earth) */
	static constexpr float TitanGravityScale = TitanGravity / EarthGravity; // 0.138

	/** Surface atmospheric density (kg/m³) */
	static constexpr double TitanSurfaceDensity = 5.34;

	/** Player drag area (Cd × A, m²) — person in thermal suit */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Titan|Physics")
	float DragCdA = 0.7f;

	/** Player mass for drag calculation (kg) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Titan|Physics")
	float PlayerMassKg = 80.0f;

	// ── HUD Data (read by UI blueprint) ──

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentTemperatureK = 93.5f;

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentPressureBar = 1.467f;

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentAltitudeM = 0.0f;

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentWindSpeedMS = 0.2f;

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentSpeedMS = 0.0f;

	UPROPERTY(BlueprintReadOnly, Category = "Titan|HUD")
	float CurrentDensityKgM3 = 5.34f;

protected:
	/** Apply atmospheric drag force opposing movement */
	void ApplyAtmosphericDrag(float DeltaTime);

	/** Update HUD telemetry from environment */
	void UpdateTelemetry();

	// ── Input ──
	void MoveForward(float Value);
	void MoveRight(float Value);
	void LookUp(float Value);
	void LookRight(float Value);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Titan|Input")
	float LookSensitivity = 1.0f;
};
