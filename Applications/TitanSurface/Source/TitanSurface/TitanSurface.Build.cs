using UnrealBuildTool;

public class TitanSurface : ModuleRules
{
	public TitanSurface(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;
		PublicDependencyModuleNames.AddRange(new string[] {
			"Core", "CoreUObject", "Engine", "InputCore",
			"Landscape", "Foliage", "KeplerWorkstation"
		});
	}
}
