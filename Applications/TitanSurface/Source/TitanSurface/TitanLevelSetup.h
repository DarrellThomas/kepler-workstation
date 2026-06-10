// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// TitanLevelSetup.h — Configures atmosphere, fog, and lighting for Titan
//
// Drop this actor into any level. On BeginPlay it finds (or spawns)
// SkyAtmosphere, ExponentialHeightFog, and DirectionalLight actors
// and applies FTitanAtmosphereParams values from TitanAtmosphereConfig.h.
//
// All visual parameters derived from Huygens HASI/DISR mission data.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Titan/TitanAtmosphereConfig.h"
#include "TitanLevelSetup.generated.h"

class ASkyAtmosphere;
class AExponentialHeightFog;
class ADirectionalLight;

UCLASS()
class TITANSURFACE_API ATitanLevelSetup : public AActor
{
	GENERATED_BODY()

public:
	ATitanLevelSetup();

	virtual void BeginPlay() override;

	/** Atmosphere visual parameters — all derived from Cassini-Huygens data */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Titan")
	FTitanAtmosphereParams AtmoParams;

	/** Day/night cycle parameters */
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Titan")
	FTitanDayNightParams DayNightParams;

protected:
	void ConfigureSkyAtmosphere();
	void ConfigureFog();
	void ConfigureSunLight();

	UPROPERTY()
	USkyAtmosphereComponent* SkyAtmoComp;

	UPROPERTY()
	UExponentialHeightFogComponent* FogComp;

	UPROPERTY()
	UDirectionalLightComponent* SunComp;
};
