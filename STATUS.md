# Kepler Workstation — Status

Aerospace simulation and mission design platform. Real physics, real data, open source.

## Last Touched

2026-06-10 — Landed TitanSurface UE5 project (Phase 0); gitignore cleanup for generated artifacts

## Quick Start

```bash
# Viewer (Three.js, no build step)
cd Viewer && python3 -m http.server 8090
# Live: https://darrellthomas.github.io/kepler-workstation/

# Standalone C++ (no UE5 needed)
cd PlanetaryModels && g++ -std=c++11 -o test test_models.cpp && ./test
cd Ephemeris && g++ -std=c++11 -O2 -o test test_propagator.cpp && ./test
cd Ephemeris/cr3bp && g++ -std=c++11 -O2 -o test test_cr3bp.cpp && ./test
cd Applications/OrbitalTransfer
g++ -std=c++11 -O2 -o hohmann_demo hohmann_demo.cpp && ./hohmann_demo
cd Applications/HypersonicEntry
g++ -std=c++11 -O2 -o entry_corridor entry_corridor.cpp && ./entry_corridor
cd Applications/SaturnMission
g++ -std=c++11 -O2 -o launch_mission launch_mission.cpp && ./launch_mission

# Python bindings + notebooks (pip install numpy matplotlib)
cd notebooks
jupyter lab  # or open .py files in VS Code
```

## Progress

### Milestone 1: The Demo
- [x] Three.js web viewer — interactive constellation, time scrubber, telemetry HUD
- [x] Deploy to GitHub Pages — auto-deploys on push, live at darrellthomas.github.io/kepler-workstation
- [ ] Screen-record 15s clip for X post
- [ ] Craft the X post

### Milestone 2: The Workstation
- [-] UE5 plugin — headers + component architecture done, implementations are stubs with TODOs
- [x] Jupiter Harvester app — constellation sim, magnetic scoop, aeropass, manifold analysis (4,794 LOC)
- [x] Saturn Harvester app — full sim, Titan base, head-to-head vs Jupiter (1,963 LOC)
- [x] OrbitalTransfer app — Hohmann demo + porkchop plot generator (809 LOC)
- [x] HypersonicEntry app — 7-body entry corridor analysis, 3-DOF RK4 sim (516 LOC)
- [x] SaturnMission app — Earth→Saturn→Earth launch script, ISRU analysis (535 LOC)
- [x] CR3BP dynamics — halo orbits, manifolds, patched 3-body, Python validation
- [-] TitanSurface app — UE5 project + game module landed (TitanGameMode, TitanLevelSetup,
      TitanPlayerCharacter, HASI atmosphere); data pipeline done; needs editor session for
      landscape import + level assembly (see Applications/TitanSurface/README.md)
- [ ] GPU compute (CUDA kernels) — not started

### Milestone 3: The Education
- [x] Tutorial docs — 5-doc progression from "Your first orbit" to "Build a constellation"
- [x] Jupyter notebooks — 5 interactive notebooks with matplotlib viz
- [x] Python bindings — `kepler` package: orbit, transfer, atmosphere, entry (6 modules, 19 KB)
- [x] Standalone mode — all models + apps compile with g++, 13 working binaries

### Milestone 4: The Message
- [x] Git history — 18 commits showing progression from initial commit to notebooks
- [x] BUILD_LOG.md — origin story, helium crisis, magnetic scoop pivot, collaboration model
- [-] Attribution — README has acknowledgments, needs Claude co-author credit

## Next Up

1. **TitanSurface editor session** — build TitanSurfaceEditor, import heightmap as landscape,
   apply SAR/ISS material, drop TitanLevelSetup actor, walk on Titan (Phase 0 finish line)
2. **Claude co-author credit** — README attribution fix
3. **Screen-record 15s clip** — X post visual (camera: behind Jupiter → harvesters arrive → scooping → return)
4. **Craft the X post** — hook: helium crisis + physics insight, visual, link, tag @elonmusk
5. **UE5 plugin** — fill stub implementations in Environment and OrbitalMechanics components
6. **GPU compute** — CUDA kernels for batch Lambert and Monte Carlo

## Blocked / Waiting

- UE5 plugin integration needs a session dedicated to porting standalone libs into UE5 components
- No JPL DE440 binary in repo (.gitignored) — users must download separately (~100 MB)
- X post needs video clip before crafting
- ~~GPU permissions for UE5 rendering~~ — resolved, darrell is in video+render groups

## Key Findings

- **Magnetic scoop at Jupiter works** (42 km/s → H2 self-ionizes at 18 eV > 15.4 eV threshold)
- **Magnetic scoop at Saturn fails** (26 km/s → H2 KE only ~7 eV, below ionization threshold)
- **Titan ISRU is the Saturn-system answer**: pump liquid methane from lakes
- **Enceladus geysers are impractical for ISRU**: plume expands into vacuum, density falls ~1/r²
- **Titan vs Mars for habitation**: Titan safer (no radiation, gentle pressure diff.), Mars better for psychology

## Key Files

| File | What |
|------|------|
| `Viewer/index.html` | Three.js constellation viewer (491 LOC) |
| `PlanetaryModels/planetary_environment.hpp` | 8-body atmosphere + gravity (1,054 LOC) |
| `PlanetaryModels/stat_thermo.hpp` | H2/He real-gas thermodynamics (710 LOC) |
| `Ephemeris/lambert.hpp` | Lambert transfer solver (248 LOC) |
| `Ephemeris/cr3bp/cr3bp.hpp` | Three-body dynamics (675 LOC) |
| `Applications/JupiterHarvester/constellation_sim.cpp` | 170-harvester fleet sim (422 LOC) |
| `Applications/OrbitalTransfer/porkchop_plot.cpp` | Launch window sweep + CSV (461 LOC) |
| `Applications/HypersonicEntry/entry_corridor.cpp` | 7-body entry corridor, 3-DOF RK4 (516 LOC) |
| `Applications/SaturnMission/launch_mission.cpp` | Earth→Saturn→Earth mission script (535 LOC) |
| `bindings/python/kepler/` | Python package — orbit, transfer, atmosphere, entry |
| `notebooks/01-your-first-orbit.py` | Interactive orbit viz + 8-body comparison |
| `notebooks/02-transfer-to-mars.py` | Porkchop contour plot, launch windows |
| `notebooks/03-atmospheric-entry.py` | Corridor sweep, trajectory comparison |
| `notebooks/04-design-a-mission.py` | Saturn round trip, Isp sensitivity |
| `notebooks/05-build-a-constellation.py` | Fleet sim, attrition, economics dashboard |
| `docs/tutorial-01-your-first-orbit.md` | Vis-viva, circular orbits, Kepler |
| `docs/tutorial-05-build-a-constellation.md` | Fleet economics, Jupiter's unique advantage |
| `BUILD_LOG.md` | Level 4 origin story (158 lines) |
| `.claude/plans/master-plan.md` | Full roadmap with 4 levels |

## Stats

- **~14,000 LOC** C++ (standalone + UE5 stubs) + **~700 LOC** Python
- **13 working C++ binaries** + **5 Jupyter notebooks**
- **8 atmospheric bodies** modeled (Earth, Venus, Mars, Jupiter, Saturn, Titan, Neptune, Uranus)
- **23 git commits**, single main branch
- **Remote:** `git@github.com:DarrellThomas/kepler-workstation.git`
- **Live:** `https://darrellthomas.github.io/kepler-workstation/`
