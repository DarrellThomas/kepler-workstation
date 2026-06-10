using UnrealBuildTool;
public class TitanSurfaceEditorTarget : TargetRules
{
	public TitanSurfaceEditorTarget(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Editor;
		DefaultBuildSettings = BuildSettingsVersion.V4;
		IncludeOrderVersion = EngineIncludeOrderVersion.Unreal5_4;
		ExtraModuleNames.Add("TitanSurface");
	}
}
