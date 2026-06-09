# Stand on Titan

Walk on the surface of Saturn's moon Titan using real Cassini-Huygens mission data.

## What This Is

A UE5 application that places you at the Huygens landing site (10.6°S, 192°W) on Titan.
Everything you see, feel, and hear is derived from actual NASA mission data:

- **Terrain**: GTDR radar altimetry heightmap (Cassini RADAR, Zebker et al. 2013)
- **Surface texture**: SAR radar mosaic (351 m/px) + ISS optical mosaic (450 m/px)
- **Atmosphere**: HASI temperature/pressure/density profile (2,727 data points from Huygens descent)
- **Wind**: DWE Doppler wind measurements (2,915 points from Huygens descent)
- **Gravity**: 1.352 m/s² (0.14g) from Cassini orbital tracking

## How It Feels

- You weigh 11 kg (80 kg person at 0.14g). Movement is bouncy.
- The air is 4.4× denser than Earth's. Walking feels like wading.
- Terminal velocity is 27 km/h. You can't die from falling.
- Jump 3.3 meters high from a standing start.
- The sky is orange in every direction. No stars visible. Perpetual twilight.
- Saturn hangs in the haze — a ghostly ringed glow 10× the size of our Moon.
- The wind is almost still: 0.2 m/s (0.4 knots). You barely feel it.
- The ground is dark orange-brown tholin organic material over water ice.
- Temperature: -180°C. Pressure: 1.47 bar. You need a thermal suit and O2.

## Data Pipeline

### 1. Heightmap
```bash
cd Tools/titan_heightmap
python3 build_titan_heightmap.py
# Output: output/titan_huygens_heightmap.r16 (2017×2017, uint16)
```

### 2. Textures
```bash
python3 extract_titan_textures.py
# Output: output/T_TitanSAR_D.png  (colorized SAR, 4096×4096)
#         output/T_TitanISS_D.png  (colorized ISS, 4096×4096)
#         output/T_TitanSAR_N.png  (normal map, 4096×4096)
```

### 3. UE5 Import
1. Create a new UE5 Level with World Partition
2. Landscape Mode → Import → select `titan_huygens_heightmap.r16`
3. Set Scale X/Y per `titan_huygens_heightmap.json` metadata
4. Apply `M_TitanSurface` material with SAR/ISS textures
5. Configure SkyAtmosphere using values from `TitanAtmosphereConfig.h`
6. Set PlayerCharacter to `TitanPlayerCharacter`

## Source Code

| File | Purpose |
|------|---------|
| `Engine/.../Public/Titan/TitanPlayerCharacter.h` | 0.14g movement, atmospheric drag |
| `Engine/.../Public/Titan/TitanAtmosphereConfig.h` | Sky, fog, lighting parameters |
| `Engine/.../Public/Sim/SimTypes.h` | Titan added to ECelestialBody enum |
| `Engine/.../Private/Sim/EnvironmentComponent.cpp` | HASI atmosphere lookup wired in |
| `PlanetaryModels/titan_hasi_table.hpp` | 273-point atmosphere from Huygens |
| `Tools/titan_heightmap/build_titan_heightmap.py` | GTDR → heightmap pipeline |
| `Tools/titan_heightmap/extract_titan_textures.py` | SAR/ISS → texture pipeline |

## Data Sources

- Huygens HASI: NASA PDS HP-SSA-HASI-2-3-4-MISSION-V1.1
- Huygens DWE: NASA PDS HP-SSA-DWE-2-3-DESCENT-V1.0
- Cassini RADAR GTDR: NASA PDS CO-SSA-RADAR-5-GTDR-V1.0
- Cassini SAR Mosaic: USGS Astrogeology (351 m/px)
- Cassini ISS Mosaic: USGS Astrogeology (450 m/px)
- Fulchignoni et al. (2005), "In situ measurements of Titan's environment", Nature 438:785

## What's Next

If standing on Titan feels magical, Phase 1 adds:
- Ice cave building (carve the walls, insulate, heat)
- O2 and heat survival meters
- Methane weather events
- Pedal-powered flight

*Kepler Workstation — MIT License. Copyright © 2026 Darrell Thomas.*
