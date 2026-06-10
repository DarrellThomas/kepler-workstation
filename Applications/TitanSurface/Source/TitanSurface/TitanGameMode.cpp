// Copyright (c) 2026 Darrell Thomas. MIT License.
#include "TitanGameMode.h"
#include "Titan/TitanPlayerCharacter.h"
#include "GameFramework/PlayerController.h"

ATitanGameMode::ATitanGameMode()
{
	DefaultPawnClass = ATitanPlayerCharacter::StaticClass();

	// First-person: let the pawn handle its own camera
	PlayerControllerClass = APlayerController::StaticClass();
}
