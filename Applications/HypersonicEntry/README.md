# Hypersonic Entry — Atmospheric Entry Corridor Analysis

Computes viable entry corridors for 7 planetary bodies using real
atmosphere models and 3-DOF ballistic trajectory simulation.

Part of Kepler Workstation.

## What's Here

| File | What it does |
|------|-------------|
| `entry_corridor.cpp` | Entry corridor sweep for Earth, Mars, Jupiter, Venus, Titan, Neptune, Uranus |
| `README.md` | This file |

## Quick Start

```bash
g++ -std=c++11 -O2 -o entry_corridor entry_corridor.cpp
./entry_corridor
```

No dependencies. Standalone. Bundled atmosphere models for all 7 bodies.

## What It Computes

For each body, sweeps entry flight path angles (γ) from -0.5° to -60°
and simulates a 3-DOF ballistic entry trajectory to find:

- **Skip-out boundary** — shallowest angle that captures. Any shallower
  and the vehicle skips out of the atmosphere.
- **G-limit boundary** — steepest angle before exceeding structural
  acceleration limits (default: 10 g).
- **Heat-rate boundary** — steepest angle before exceeding TPS heat flux
  limits (default: 1500 kW/m², carbon phenolic).
- **Corridor width** — total viable entry angle range.

## Bodies Analyzed

| Body | Entry v (km/s) | Interface (km) | Atmosphere source |
|------|---------------|----------------|-------------------|
| Earth | 11.0 | 120 | US Standard Atmosphere 1976 |
| Mars | 6.0 | 125 | Mars Global Surveyor (NASA GRC) |
| Jupiter | 48.0 | 500 | Galileo Probe ASI (Seiff et al. 1998) |
| Venus | 11.5 | 150 | Venera/Pioneer Venus |
| Titan | 6.5 | 500 | Cassini-Huygens (Niemann et al. 2005) |
| Neptune | 23.0 | 500 | Voyager 2 RSS (Lindal et al. 1990) |
| Uranus | 22.0 | 500 | Voyager 2 RSS (Lindal et al. 1987) |

## Vehicle Model

Default: Apollo/Orion-class blunt-body capsule.

| Parameter | Value |
|-----------|-------|
| Mass | 5,500 kg |
| Reference area | 12.57 m² (4m diameter) |
| CD | 1.3 (blunt body) |
| Nose radius | 2.0 m |
| G-limit | 10 g |
| Heat-rate limit | 1,500 kW/m² |

## Physics

- **Trajectory**: 3-DOF (altitude, velocity, flight path angle), RK4 integration
- **Drag**: D = ½ρv²CD A, ballistic entry (no lift)
- **Gravity**: inverse-square, g = GM/(R+h)²
- **Heating**: Sutton-Graves stagnation heat flux correlation
  q = 1.83×10⁻⁴ √(ρ/r_n) v³ kW/m²
- **Atmosphere**: real layered models for each body (piecewise T, hydrostatic p, ideal gas ρ)

## Interpreting Results

- **Skip-out ~ -5° to -7°** on most bodies means you need at least
  that steep an entry to capture
- **G-limit ~ -15° to -30°** varies dramatically with entry velocity
  and atmospheric scale height
- **Jupiter** is the hardest entry — 48 km/s means enormous heating,
  very narrow corridor
- **Mars** is relatively gentle — low velocity, thin atmosphere gives
  wide corridors
- **Titan** is unique — dense N2 atmosphere (1.47 bar surface) but
  low gravity (1.35 m/s²), producing a different entry regime

## References

- Allen & Eggers (1958) — "A Study of the Motion and Aerodynamic
  Heating of Ballistic Missiles Entering the Earth's Atmosphere at
  High Supersonic Speeds", NACA TR 1381
- Sutton & Graves (1971) — "A General Stagnation-Point Convective-
  Heating Equation for Arbitrary Gas Mixtures", NASA TR R-382
- Chapman (1959) — "An Approximate Analytical Method for Studying
  Entry into Planetary Atmospheres", NASA TR R-11
- US Standard Atmosphere 1976 — NOAA-S/T 76-1562
- Galileo Probe — Seiff et al., Science 272:844-845, 1996
- Voyager 2 RSS — Lindal et al., JGR 92:14987 (1987), JGR 95:16949 (1990)
- Cassini-Huygens — Niemann et al., Nature 438:779-784, 2005

## Next Steps

- Lifting entry (add L/D ratio parameter)
- Skip-entry trajectories (Apollo-style)
- Aerocapture analysis (single-pass capture with ΔV savings)
- Heat shield sizing (integrated heat load, not just peak)
