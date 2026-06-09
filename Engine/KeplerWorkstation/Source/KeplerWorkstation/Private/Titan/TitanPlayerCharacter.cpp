// Copyright (c) 2026 Darrell Thomas. MIT License.
#include "Titan/TitanPlayerCharacter.h"
#include "GameFramework/CharacterMovementComponent.h"
#include "Camera/CameraComponent.h"
#include "Components/CapsuleComponent.h"

ATitanPlayerCharacter::ATitanPlayerCharacter()
{
	PrimaryActorTick.bCanEverTick = true;

	// Attach environment component for Titan atmosphere/gravity queries
	Environment = CreateDefaultSubobject<UEnvironmentComponent>(TEXT("TitanEnvironment"));
	Environment->Body = ECelestialBody::Titan;
}

void ATitanPlayerCharacter::BeginPlay()
{
	Super::BeginPlay();

	UCharacterMovementComponent* MoveComp = GetCharacterMovement();
	if (!MoveComp) return;

	// ── Titan gravity: 0.138× Earth ──
	MoveComp->GravityScale = TitanGravityScale;

	// ── Movement tuning for low-g + dense atmosphere ──

	// Walking: slow in thermal suit, ~2 m/s max
	MoveComp->MaxWalkSpeed = 200.0f;         // cm/s
	MoveComp->MaxWalkSpeedCrouched = 100.0f;

	// Jumping: in 0.14g, v² = 2gh → h = v²/(2g)
	// With JumpZ=300 cm/s: h = 300²/(2×135.4) = 332 cm = 3.3m
	MoveComp->JumpZVelocity = 300.0f;        // cm/s

	// Dense atmosphere provides air control while airborne
	MoveComp->AirControl = 0.5f;

	// Atmospheric drag while falling (dense air = terminal velocity limit)
	// Terminal velocity: sqrt(2mg/(rho×Cd×A)) = sqrt(2×80×1.352/(5.34×0.7)) = 7.6 m/s
	MoveComp->FallingLateralFriction = 0.3f;

	// Ground friction: tholin-coated ice is moderately slippery
	MoveComp->GroundFriction = 5.0f;
	MoveComp->BrakingDecelerationWalking = 400.0f;

	// No fall damage — terminal velocity is only 7.6 m/s (27 km/h)
	// Can't die from falling at any height on Titan
}

void ATitanPlayerCharacter::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

	ApplyAtmosphericDrag(DeltaTime);
	UpdateTelemetry();
}

void ATitanPlayerCharacter::ApplyAtmosphericDrag(float DeltaTime)
{
	// Titan surface density: 5.34 kg/m³ — 4.4× Earth
	// Moving through this feels like wading through thick air
	//
	// F_drag = 0.5 × Cd×A × ρ × v²
	// At 2 m/s:  F = 0.5 × 0.7 × 5.34 × 4.0 = 7.5 N (noticeable)
	// At 5 m/s:  F = 0.5 × 0.7 × 5.34 × 25  = 46.7 N (significant)
	// At 7.6 m/s (terminal): F = mg = 108 N (equilibrium)

	UCharacterMovementComponent* MoveComp = GetCharacterMovement();
	if (!MoveComp) return;

	FVector Velocity = MoveComp->Velocity;
	float SpeedCmS = Velocity.Size();

	if (SpeedCmS < 1.0f) return; // negligible

	float SpeedMS = SpeedCmS / 100.0f;

	// Get local atmospheric density from environment component
	float LocalDensity = TitanSurfaceDensity;
	if (Environment)
	{
		LocalDensity = Environment->Atmosphere.Density;
		if (LocalDensity < 0.001) LocalDensity = TitanSurfaceDensity;
	}

	// Drag force (Newtons)
	float DragForceN = 0.5f * DragCdA * LocalDensity * SpeedMS * SpeedMS;

	// Drag deceleration (m/s²)
	float DragDecelMS2 = DragForceN / PlayerMassKg;

	// Convert to UE units (cm/s²) and apply as opposing force
	// Clamp to prevent oscillation: don't decelerate more than current speed
	float MaxDecel = SpeedMS / DeltaTime; // would stop in one frame
	DragDecelMS2 = FMath::Min(DragDecelMS2, MaxDecel * 0.5f);

	FVector DragDir = -Velocity.GetSafeNormal();
	float DragForceCmS2 = DragDecelMS2 * 100.0f; // convert m/s² to cm/s²

	// Apply as velocity change (impulse approach for stability)
	FVector DragImpulse = DragDir * DragForceCmS2 * DeltaTime;
	MoveComp->Velocity += DragImpulse;
}

void ATitanPlayerCharacter::UpdateTelemetry()
{
	if (!Environment) return;

	// Query environment at current position
	FVector Pos = GetActorLocation();
	// Convert UE position to planet-centered (for now, assume flat surface)
	// In full implementation, this would account for curved planetary surface
	double AltM = Pos.Z / 100.0; // UE cm to meters, Z-up
	if (AltM < 0) AltM = 0;

	CurrentAltitudeM = AltM;
	CurrentTemperatureK = Environment->Atmosphere.Temperature;
	CurrentPressureBar = Environment->Atmosphere.Pressure / 1e5;
	CurrentDensityKgM3 = Environment->Atmosphere.Density;
	CurrentWindSpeedMS = Environment->Atmosphere.Wind.Size();
	CurrentSpeedMS = GetCharacterMovement()->Velocity.Size() / 100.0f;
}

void ATitanPlayerCharacter::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
	Super::SetupPlayerInputComponent(PlayerInputComponent);

	PlayerInputComponent->BindAxis("MoveForward", this, &ATitanPlayerCharacter::MoveForward);
	PlayerInputComponent->BindAxis("MoveRight", this, &ATitanPlayerCharacter::MoveRight);
	PlayerInputComponent->BindAxis("LookUp", this, &ATitanPlayerCharacter::LookUp);
	PlayerInputComponent->BindAxis("Turn", this, &ATitanPlayerCharacter::LookRight);
	PlayerInputComponent->BindAction("Jump", IE_Pressed, this, &ACharacter::Jump);
}

void ATitanPlayerCharacter::MoveForward(float Value)
{
	if (FMath::Abs(Value) > 0.01f)
	{
		FVector Forward = GetActorForwardVector();
		AddMovementInput(Forward, Value);
	}
}

void ATitanPlayerCharacter::MoveRight(float Value)
{
	if (FMath::Abs(Value) > 0.01f)
	{
		FVector Right = GetActorRightVector();
		AddMovementInput(Right, Value);
	}
}

void ATitanPlayerCharacter::LookUp(float Value)
{
	AddControllerPitchInput(Value * LookSensitivity);
}

void ATitanPlayerCharacter::LookRight(float Value)
{
	AddControllerYawInput(Value * LookSensitivity);
}
