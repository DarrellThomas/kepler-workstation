# Kepler Workstation

An open-source aerospace simulation and mission design platform built on Unreal Engine 5.

Design orbits. Simulate atmospheres. Plan interplanetary missions. From LEO to Jupiter and back.

## What it does

- **6DOF vehicle simulation** — rigid body dynamics with translational and rotational equations of motion
- **Planetary environments** — atmosphere and gravity models for Earth, Venus, Mars, Jupiter, and Moon
- **Orbital mechanics** — N-body propagation, Lambert transfer solver, JPL DE440 ephemeris (1550-2650 CE)
- **3D visualization** — real-time rendering in Unreal Engine 5 with flight data overlay
- **GPU acceleration** — CUDA support for Monte Carlo analysis and batch trajectory computation

## Who it's for

- **Students** learning orbital mechanics, atmospheric flight, or aerospace engineering
- **Engineers** designing missions, analyzing trajectories, or validating simulations
- **Researchers** studying planetary atmospheres, gravity fields, or interplanetary transfers
- **Anyone** curious about how things fly through space

## Planetary Models

Real data from real missions:

| Body | Atmosphere | Gravity | Source |
|------|-----------|---------|--------|
| Earth | US Standard 1976 | WGS-84 J2-J6 | NOAA, NIMA |
| Venus | Venera/Pioneer Venus | Magellan J2-J4 | NASA PDS |
| Mars | Mars Global Surveyor | GMM-3 J2-J6 | NASA GRC, Genova et al. |
| Jupiter | Galileo Probe ASI | Juno J2-J10 | Seiff et al., Iess et al. |
| Moon | — (vacuum) | GRAIL J2-J6 | Zuber/Konopliv et al. |

## Quick Start

```bash
# Clone
git clone https://github.com/[your-org]/kepler-workstation.git

# Planetary models (standalone, no UE5 needed)
cd PlanetaryModels
g++ -std=c++11 -o test_models test_models.cpp && ./test_models

# Solar system propagator (standalone)
cd ../Ephemeris
# Download JPL DE440 (~100 MB):
wget https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440
cd ../Applications/OrbitalTransfer
g++ -std=c++11 -O2 -o lambert_demo lambert_demo.cpp && ./lambert_demo

# Full UE5 workstation (requires UE 5.4+ source build)
# See Engine/README.md for setup
```

## Architecture

```
kepler-workstation/
├── Engine/                      UE5 plugin (simulation + visualization)
│   └── KeplerWorkstation/
│       └── Source/
│           └── KeplerWorkstation/
│               ├── Public/Sim/  Simulation components (headers)
│               └── Private/Sim/ Implementations
├── PlanetaryModels/             Standalone atmosphere + gravity (C++)
├── Ephemeris/                   JPL DE440 reader + N-body propagator
├── Applications/                Example missions
│   ├── JupiterHarvester/        Magnetic scoop helium collection
│   ├── HypersonicEntry/         Atmospheric entry simulation
│   └── OrbitalTransfer/         Lambert transfers, gravity assists
└── Viewer/                      Three.js web viewer for demos
```

## Simulation Components

All cleanroom implementations from published physics. No proprietary code.

| Component | What it does |
|-----------|-------------|
| `EnvironmentComponent` | Atmosphere density/pressure/temperature + gravity with zonal harmonics |
| `NewtonComponent` | Translational EOM (F=ma), RK4 integration |
| `EulerComponent` | Rotational EOM, quaternion kinematics |
| `OrbitalMechanicsComponent` | N-body propagation, Lambert solver, orbital elements |
| `MagneticScoopComponent` | Magnetic funnel collection physics |

## License

MIT License. Copyright (c) 2026 Darrell Thomas.

Use it freely. Learn from it. Build on it. Go somewhere.

## Acknowledgments

- Peter H. Zipfel — whose textbooks and teaching inspired the mathematical framework
- JPL — DE440 planetary ephemeris (public domain)
- NASA — Galileo, Juno, Magellan, GRAIL, Mars Global Surveyor mission data
- Epic Games — Unreal Engine 5
