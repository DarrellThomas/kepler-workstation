# Tutorial 5: Build a Constellation

**Prerequisites:** Tutorials 1–4.

**You'll learn:** How to turn a single mission into a fleet — and what happens
when you run the numbers for 36 years.

**Time:** 25 minutes.

---

## From One Bucket to One Hundred

A single Jupiter Helium Harvester ("Bucket") delivers 8.2 tonnes of helium
per decade-long Earth-Jupiter-Earth loop. That's useful — it proves the
physics. But 8.2 tonnes/year doesn't change the world.

A fleet does.

The concept: launch 10 Buckets per Starship during each Earth-Jupiter synodic
window (every 399 days). Over 10 windows, you have 100 Buckets in the pipeline.
Each one runs the loop independently. Some fail — attrition is planned at 10%
per loop. Replacements keep the fleet at steady state.

---

## The Fleet Physics

Run the constellation simulation:

```bash
cd Applications/JupiterHarvester
g++ -std=c++11 -O2 -o constellation constellation_sim.cpp \
    -I../../Ephemeris -I../../PlanetaryModels
# (requires DE440 ephemeris at ../../Ephemeris/linux_p1550p2650.440)
./constellation
```

The simulation tracks:
- **170 harvesters** (17 windows × 10/window for the expanded fleet)
- **36 years** of operations
- Each harvester: outbound → scoop → return → depot → repeat
- Attrition at each phase
- Cumulative helium delivered

---

## The Fleet Numbers

| Parameter | Value |
|-----------|-------|
| Buckets in fleet | 100 (steady state) |
| Buckets launched per window | 10 |
| Launch cadence | Every 399 days (Earth-Jupiter synodic) |
| Pipeline depth | 10 windows |
| Loop duration | ~10 years per Bucket |
| Attrition rate | 10% per loop (planned) |
| Annual helium delivery | **179 tonnes** |
| Global helium demand | ~6,000 tonnes/year (pre-Ras Laffan) |
| Fleet contribution | ~3% of global demand |

179 tonnes of helium per year, forever. At $1,500/kg (post-crisis spot price),
that's $269 million/year in helium alone. The hydrogen delivered alongside
(13.4 tonnes per Bucket per loop) is orbital propellant — "gas stations in
space."

---

## The Economics

| Item | Cost |
|------|------|
| Bucket build cost | $200M (HTS magnet, Hall thruster, avionics, tank) |
| Launch cost (Starship, 10-pack) | $45M ($4.5M per Bucket) |
| Total per Bucket | $245M |
| Fleet startup (100 Buckets) | $24.5B |
| Annual replacement (9/year at 10% attrition) | $2.25B/year |
| Annual revenue at $1,500/kg He | $269M/year |

The business case does NOT close on helium alone. But it was never about
helium revenue. The question is: **what missions become possible when there are
gas stations in space?**

The hydrogen propellant delivered to Earth orbit enables Mars missions, lunar
base logistics, asteroid mining, and deep-space telescopes with orbital
servicing. The helium is a byproduct. The orbital propellant infrastructure is
the product.

---

## Attrition — Why Buckets Fail

| Failure mode | Rate | Mitigation |
|-------------|------|-----------|
| Jupiter radiation dose | 4% per loop | Shielded avionics bay |
| Micrometeoroid strike | 2% | Redundant systems, Whipple shield |
| HTS quench | 1.5% | Active thermal management |
| Hall thruster wear | 1.5% | Designed for multi-loop life |
| Navigation error | 0.5% | Autonomous star tracker + Earth updates |
| Unknown | 0.5% | Margin in fleet sizing |
| **Total** | **10%** | Planned replacement pipeline |

At 10% attrition per loop and 10-year loops, you need to replace 10 Buckets
per year. The factory builds 9 per year (steady-state), with surge capacity
for bad years. The pipeline has slack — you can absorb a 15% attrition year
by dipping into the replacement queue.

---

## Why It Works at Jupiter

Jupiter is unique in the solar system for the magnetic scoop concept:

1. **Orbital velocity at scoop altitude: 42 km/s.** This is the key. H2
   molecules at 42 km/s have 18 eV of kinetic energy — above the 15.4 eV
   ionization threshold. The bow shock self-ionizes the incoming gas. No
   ionizer needed.

2. **Magnetic field: 4.2 Gauss at cloud tops.** The strongest planetary
   magnetic field in the solar system. The 20T superconducting coil on the
   Bucket generates a 30–50 meter funnel that sweeps ionized H2/He into
   the collection tank.

3. **Oberth effect at departure.** Burning at periapsis (closest approach)
   gives the maximum kinetic energy gain. Jupiter's deep gravity well means
   a 5 km/s burn at periapsis produces enough energy to escape to Earth.

4. **Atmospheric composition: 10% helium by number density.** 25% by mass.
   The highest helium concentration of any planet. Earth's atmosphere is
   0.0005% helium.

This combination — high velocity, strong magnetic field, deep gravity well,
helium-rich atmosphere — doesn't exist anywhere else. Saturn is too slow
(26 km/s, won't ionize). Uranus and Neptune are too far, too cold, too slow.
Jupiter is the only gas station in the solar system that opens itself when
you show up at the right speed.

---

## The Big Picture

A single Bucket is an engineering marvel — 8 tonnes, 3 moving parts, 1 magnet,
autonomous for a decade. A fleet of 100 is infrastructure. Infrastructure
changes what's possible.

Before the transcontinental railroad, crossing North America took months and
killed people. After, it took days and carried freight. The railroad didn't
make money on passenger tickets — it made money on what the passengers built
at the other end.

The Bucket fleet is the railroad. The helium pays for the track. What gets
built at the other end — orbital manufacturing, asteroid mining, deep-space
science — is the business.

---

## What You've Learned

| Tutorial | Concept | Tool |
|----------|---------|------|
| 1. Your First Orbit | Vis-viva, circular orbits | `orbit.cpp` |
| 2. Transfer to Mars | Hohmann, Lambert, porkchop | `hohmann_demo`, `porkchop_plot` |
| 3. Atmospheric Entry | Drag, heating, corridor | `entry_corridor` |
| 4. Design a Mission | ΔV budgets, ISRU, rocket eq | `launch_mission` |
| 5. Build a Constellation | Fleet sim, attrition, economics | `constellation_sim` |

You can now:

- Compute any circular orbit around any body
- Find the optimal transfer window between any two planets
- Determine whether a vehicle survives atmospheric entry
- Design a complete mission and check if it closes
- Simulate a fleet of autonomous vehicles over decades
- Understand why ISRU is the single biggest lever in spaceflight economics

Now go build something.

---

*Kepler Workstation — MIT License. Real physics, real data, open source.*
