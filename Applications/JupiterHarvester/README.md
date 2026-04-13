# THE BUCKET

## One sentence

An 8-tonne autonomous magnetic scoop flies to Jupiter, fills up on H2/He,
and flies home. Repeat forever.

## The Vehicle

- **Dry mass:** 8 tonnes
- **Moving parts:** 3 valves, 3 reaction wheels
- **Magnet:** 20 T HTS superconducting solenoid, 3m radius
- **Tank:** COPV, ~500 m³ (~10m diameter sphere), PTFE liner (no H2 embrittlement)
- **Thruster:** Hall effect, H2 propellant, Isp = 3000s (Starlink heritage)
- **Power:** H2 fuel cell + solar panels
- **Brain:** Autonomous nav/comms (Starlink-class avionics)
- **Heat shield:** None. The magnetic field IS the shield.
- **Size:** 12m × 12m × 15m. A large satellite, not a spaceship.

## The Loop

```
EARTH DEPOT ──→ Jupiter (5 yr) ──→ SCOOP (2 yr) ──→ EARTH DEPOT (3 yr)
     ↑                                                        │
     └────────────────── refuel, repeat ──────────────────────┘
```

Total loop: ~10 years. Vehicle does this until it breaks.

## How It Works

1. **Launch:** 10 Buckets per Starship. Deploy on VEEJ trajectory.
2. **Cruise to Jupiter:** 5 years. Venus gravity assists reduce departure ΔV to 3.5 km/s.
3. **Jupiter capture:** Small burn (0.5 km/s) + magnetic drag in upper atmosphere.
4. **Scoop:** Elliptical orbit, periapsis at 300-400 km above 1-bar level.
   - Magnetic funnel (30-50m radius) sweeps ionized H2/He into tank.
   - Bow shock self-ionizes the gas at 42+ km/s (18 eV > 15.4 eV ionization).
   - Vehicle hull stays behind magnetic shield — near vacuum, no heating.
   - ~10 kg per pass × 1000+ passes over 2 years = 33 tonnes.
   - Every pass adds fuel. The vehicle refuels itself.
5. **Jupiter departure:** Oberth burn at periapsis (5 km/s), coast to escape.
6. **Return cruise:** 3 years. Hall thruster burns captured H2.
7. **Earth arrival:** Aerocapture. Free.
8. **Depot delivery:** 21.7 tonnes (13.4t H2 + 8.2t He) to orbital depot.
9. **Refuel and repeat.**

## Per-Loop Numbers

| Item | Mass |
|------|------|
| Scooped at Jupiter | 33 tonnes (24.8t H2, 8.2t He) |
| H2 burned as fuel | 11.3 tonnes |
| H2 delivered | 13.4 tonnes |
| He delivered | 8.2 tonnes |
| **Total delivered** | **21.7 tonnes** |
| Delivery ratio | **2.7× dry mass** |

## ΔV Budget

| Maneuver | ΔV |
|----------|-----|
| Earth departure (VEEJ) | 3.5 km/s |
| Jupiter capture | 0.5 km/s |
| Jupiter departure (Oberth) | 5.0 km/s |
| Earth arrival (aerocapture) | 0.0 km/s |
| Margin | 0.5 km/s |
| **Total per loop** | **9.5 km/s** |

Hall thruster Isp = 3000s → mass ratio = 1.38 → closes with 13.4t H2 payload.

## The Fleet

| Parameter | Value |
|-----------|-------|
| Buckets per Starship | 10 |
| Launch window | every 399 days |
| Pipeline depth | 10 windows × 10 buckets = **100 buckets** |
| Attrition | 10% per loop (planned) |
| Annual delivery | **179 tonnes/year** |
| Build rate | 9/year (replaces losses + growth) |

## Economics

| Item | Value |
|------|-------|
| Cost per bucket | $245M (build + launch) |
| Fleet startup (100 buckets) | $24.5B |
| Annual replacement | $2.25B/year |
| Value at $1500/kg avoided | $269M/year |
| Breakeven | When orbital propellant enables Mars/asteroid missions |

The business case isn't $/kg delivered. It's: **what missions become possible
when there are gas stations in space?**

## Key Physics

1. **Self-ionization:** At 42 km/s, H2 kinetic energy (18 eV) exceeds ionization
   energy (15.4 eV). The bow shock creates plasma. No ionizer needed.

2. **Magnetic deflection:** Ram pressure at 350 km is ~2.8 Pa. A 20T coil deflects
   this with a 30m radius funnel. Effective area: 3000+ m².

3. **No thermal problem:** Hull sits behind the magnetic shield in near-vacuum.
   No heat shield. No ablation. No thermal cycling. Unlimited reuse.

4. **H2 embrittlement solved:** COPV tanks with PTFE liner. Zero metal-H2 contact.
   Active cryocooler (fuel cell powered) manages boiloff.

5. **Gravity assists:** VEEJ trajectory (Venus-Earth-Earth-Jupiter) reduces
   departure ΔV from 6.3 to 3.5 km/s. Venus flyby is free.

## What This Is NOT

- Not a Europa lander. No landing on anything. Pure orbital mechanics.
- Not a chemical rocket. H2-only propulsion (Hall thruster).
- Not a one-off mission. A permanent infrastructure loop.
- Not a spaceship. A bucket.

## Simulation Code

All models validated:
- `planetary_models/planetary_environment.hpp` — Jupiter/Venus/Earth/Mars/Moon atmosphere + gravity
- `solar_system/jpl_ephemeris.hpp` — JPL DE440 ephemeris reader
- `solar_system/nbody_propagator.hpp` — N-body integration, Sun + all planets
- `solar_system/lambert.hpp` — Transfer orbit solver
- `solar_system/veej_simulation.cpp` — VEEJ trajectory with real ephemeris
- `solar_system/magnetic_scoop_sim.cpp` — Magnetic scoop harvest campaign

## Open Questions

1. Magnetic field leakage rate — drives residual hull heating
2. HTS coil mass at 20T — current tech vs near-term
3. H2 Hall thruster performance data at 3000s Isp
4. Autonomous navigation over 10-year loops
5. Cryogenic H2 storage over 3-year return cruise (boiloff budget)
6. Jupiter radiation belt dose on electronics
7. Optimal scoop altitude (300 vs 350 vs 400 km) — bigger magnet vs more heat

## The Tagline

*"The solar system has gas stations. We're building the pumps."*
