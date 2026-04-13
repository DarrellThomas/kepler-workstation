// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// MagneticScoopComponent.h — Magnetic funnel atmospheric harvesting
//
// Models a superconducting solenoid that creates a magnetic funnel
// to deflect ionized atmospheric gas into a collection tank.
// No physical contact with atmosphere — the magnetic field IS the scoop.
//
// Physics:
//   - Bow shock ionization: KE_molecule > ionization_energy at v > 40 km/s
//   - Magnetic pressure: P_B = B²/(2μ₀) must exceed ram pressure P_ram = ½ρv²
//   - Dipole field scaling: B(r) = B_coil * (r_coil/r)³
//   - Effective collection area: A_eff = π * R_mag²
//   - Mass flow: ṁ = ρ * v * A_eff * η

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Sim/SimTypes.h"
#include "MagneticScoopComponent.generated.h"

UCLASS(ClassGroup=(Simulation), meta=(BlueprintSpawnableComponent))
class CADACVIS_API UMagneticScoopComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UMagneticScoopComponent();

	/** Coil magnetic field strength at bore (Tesla) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Scoop")
	double CoilField = 20.0;

	/** Coil radius (meters) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Scoop")
	double CoilRadius = 3.0;

	/** Collection efficiency (0-1) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Scoop")
	double CollectionEfficiency = 0.5;

	/** Tank capacity (kg) */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Scoop")
	double TankCapacity = 33000.0;

	/** Current H2 in tank (kg) */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double H2InTank = 0;

	/** Current He in tank (kg) */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double HeInTank = 0;

	/** Effective magnetic funnel radius (m) — computed */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double FunnelRadius = 0;

	/** Effective collection area (m²) — computed */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double EffectiveArea = 0;

	/** Instantaneous mass flow rate (kg/s) — computed */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double MassFlowRate = 0;

	/** Is the gas self-ionizing at current velocity? */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	bool bSelfIonizing = false;

	/** Hull heating from magnetic leakage (kW/cm²) */
	UPROPERTY(BlueprintReadOnly, Category = "Scoop")
	double HullHeating = 0;

	/**
	 * Evaluate scoop performance at current conditions.
	 * @param Density     — atmospheric density (kg/m³)
	 * @param Velocity    — orbital velocity (m/s)
	 * @param Temperature — atmospheric temperature (K)
	 * @param DeltaTime   — time step (s) for tank filling
	 * @return Mass harvested this step (kg)
	 */
	UFUNCTION(BlueprintCallable, Category = "Scoop")
	double Evaluate(double Density, double Velocity, double Temperature, double DeltaTime);

	/** Is the tank full? */
	UFUNCTION(BlueprintCallable, Category = "Scoop")
	bool IsTankFull() const { return (H2InTank + HeInTank) >= TankCapacity; }

	/** Reset tank to empty */
	UFUNCTION(BlueprintCallable, Category = "Scoop")
	void EmptyTank() { H2InTank = 0; HeInTank = 0; }
};
