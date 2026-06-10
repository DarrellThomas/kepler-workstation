// Copyright (c) 2026 Darrell Thomas. MIT License.
//
// TitanGameMode.h — Game mode for Titan surface exploration
//
// Spawns ATitanPlayerCharacter as the default pawn so the player
// experiences 0.14g gravity, atmospheric drag, and HASI-driven telemetry.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/GameModeBase.h"
#include "TitanGameMode.generated.h"

UCLASS()
class TITANSURFACE_API ATitanGameMode : public AGameModeBase
{
	GENERATED_BODY()

public:
	ATitanGameMode();
};
