// Copyright (c) 2026 Darrell Thomas. MIT License.
#include "TitanLevelSetup.h"
#include "Components/SkyAtmosphereComponent.h"
#include "Components/ExponentialHeightFogComponent.h"
#include "Components/DirectionalLightComponent.h"
#include "Engine/World.h"
#include "Kismet/GameplayStatics.h"
#include "Atmosphere/AtmosphericFog.h"

ATitanLevelSetup::ATitanLevelSetup()
{
	PrimaryActorTick.bCanEverTick = false;

	// Root scene component
	USceneComponent* Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);

	// SkyAtmosphere — owned by this actor so it always exists
	SkyAtmoComp = CreateDefaultSubobject<USkyAtmosphereComponent>(TEXT("SkyAtmosphere"));
	SkyAtmoComp->SetupAttachment(Root);

	// Exponential Height Fog
	FogComp = CreateDefaultSubobject<UExponentialHeightFogComponent>(TEXT("HeightFog"));
	FogComp->SetupAttachment(Root);

	// Directional Light (Sun)
	SunComp = CreateDefaultSubobject<UDirectionalLightComponent>(TEXT("SunLight"));
	SunComp->SetupAttachment(Root);
}

void ATitanLevelSetup::BeginPlay()
{
	Super::BeginPlay();

	ConfigureSkyAtmosphere();
	ConfigureFog();
	ConfigureSunLight();
}

void ATitanLevelSetup::ConfigureSkyAtmosphere()
{
	if (!SkyAtmoComp) return;

	// Planet geometry
	SkyAtmoComp->BottomRadius = AtmoParams.PlanetRadiusKm;
	SkyAtmoComp->AtmosphereHeight = AtmoParams.AtmosphereHeightKm;

	// Rayleigh scattering (N2 molecular)
	SkyAtmoComp->RayleighScatteringScale = AtmoParams.RayleighScatteringScale;
	SkyAtmoComp->RayleighScattering = AtmoParams.RayleighScattering;
	SkyAtmoComp->RayleighExponentialDistribution = AtmoParams.RayleighExponentialDistribution;

	// Mie scattering (tholin haze — dominates Titan's orange appearance)
	SkyAtmoComp->MieScatteringScale = AtmoParams.MieScatteringScale;
	SkyAtmoComp->MieScattering = AtmoParams.MieScattering;
	SkyAtmoComp->MieAbsorptionScale = AtmoParams.MieAbsorptionScale;
	SkyAtmoComp->MieAbsorption = AtmoParams.MieAbsorption;
	SkyAtmoComp->MieAnisotropy = AtmoParams.MieAnisotropy;
	SkyAtmoComp->MieExponentialDistribution = AtmoParams.MieExponentialDistribution;

	// Multi-scattering: critical for thick haze atmosphere
	SkyAtmoComp->MultiScatteringFactor = AtmoParams.MultiScatteringFactor;

	// Ground albedo: dark tholin surface
	SkyAtmoComp->GroundAlbedo = AtmoParams.GroundAlbedo;
}

void ATitanLevelSetup::ConfigureFog()
{
	if (!FogComp) return;

	// Titan's lower atmosphere: 1-3 km visibility through tholin haze
	FogComp->FogDensity = AtmoParams.FogDensity;
	FogComp->FogHeightFalloff = AtmoParams.FogHeightFalloff;
	FogComp->FogInscatteringColor = AtmoParams.FogInscatteringColor;
	FogComp->FogMaxOpacity = AtmoParams.FogMaxOpacity;

	// Volumetric fog for depth and atmosphere
	FogComp->bEnableVolumetricFog = AtmoParams.bVolumetricFog;
	FogComp->VolumetricFogScatteringDistribution = AtmoParams.VolumetricFogScatteringDistribution;
	FogComp->VolumetricFogAlbedo = AtmoParams.VolumetricFogAlbedo;
}

void ATitanLevelSetup::ConfigureSunLight()
{
	if (!SunComp) return;

	// Titan receives ~0.1% of Earth's surface sunlight at 9.54 AU
	SunComp->Intensity = AtmoParams.SunIntensityLux;
	SunComp->LightColor = AtmoParams.SunColor.ToFColor(true);
	SunComp->LightSourceAngle = AtmoParams.SunSourceAngle;

	// Sun is an atmospheric light source for SkyAtmosphere
	SunComp->bAtmosphereSunLight = true;

	// Low elevation: perpetual twilight feel
	// Solar elevation from day/night params, azimuth arbitrary
	FRotator SunDir(-FMath::RadiansToDegrees(DayNightParams.SolarElevation), 180.0f, 0.0f);
	SunComp->SetWorldRotation(SunDir);
}
