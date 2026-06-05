# Tutorial 3: Atmospheric Entry

**Prerequisites:** Tutorials 1 and 2.

**You'll learn:** What happens when an orbit hits an atmosphere, and how to
survive it.

**Time:** 25 minutes.

---

## The Problem

You've spent 259 days coasting from Earth to Mars. You're approaching at
6 km/s. The atmosphere ahead is thin — less than 1% of Earth's surface
pressure — but at 6 km/s, even thin air becomes a wall of plasma.

You have one shot. Too shallow and you skip off the atmosphere into deep space.
Too steep and the g-loads crush your spacecraft. The viable range of entry
angles is called the *entry corridor*.

---

## The Physics

A blunt-body capsule entering an atmosphere experiences:

**Drag deceleration:**
```
a_drag = ½ρv² / β
```
where `β = m / (C_D · A)` is the ballistic coefficient. Lower β = more drag
per unit mass = gentler deceleration.

**Stagnation heating** (Sutton-Graves correlation):
```
q = 1.83×10⁻⁴ · √(ρ / r_n) · v³   kW/m²
```
where `r_n` is the nose radius. Larger nose = lower heat flux.

**Gravity:**
```
g = GM / (R + h)²
```

The atmosphere provides the braking. The vehicle absorbs the heat.

---

## The Corridor

Run the entry corridor analysis:

```bash
cd Applications/HypersonicEntry
g++ -std=c++11 -O2 -o entry_corridor entry_corridor.cpp
./entry_corridor
```

This simulates ballistic entry for 7 planetary bodies. For Earth at 11 km/s
(lunar return speed):

| Boundary | Angle | Meaning |
|----------|-------|---------|
| Skip-out | ~-5.5° | Any shallower and you bounce off |
| G-limit | ~-8° to -10° | Steep enough to exceed 10g |
| Heat limit | ~-7° to -9° | Steep enough to melt TPS |

**Corridor width: ~2–4°.** That's what Apollo had to hit — a 4° window from
384,000 km away, with 1960s navigation. They hit it every time.

---

## How Each Body Compares

| Body | Entry v (km/s) | Atmosphere | Hardest part |
|------|---------------|------------|-------------|
| Earth | 11.0 | US Standard 1976, 0–85 km | Moderate heating, narrow corridor |
| Mars | 6.0 | CO₂, thin, 0–125 km | Low density = need steep entry to capture |
| Jupiter | 48.0 | H₂/He, deep, 0–500 km | Extreme heating — Galileo probe hit 230g |
| Venus | 11.5 | CO₂, 92 bar at surface | Dense atmo = easy capture, brutal surface |
| Titan | 6.5 | N₂, 1.47 bar, 0–500 km | Low gravity, dense atmo — easiest entry |
| Neptune | 23.0 | H₂/He, 0–500 km | High velocity, cold atmo |
| Uranus | 22.0 | H₂/He, coldest atmo | Similar to Neptune, slightly gentler |

---

## Apollo vs Mars vs Jupiter

**Apollo (Earth, 11 km/s):**
- Blunt-body capsule (β ~ 350 kg/m²)
- Ablative TPS (AVCOAT)
- Peak g: ~6.5g
- Peak heating: ~1,000 kW/m²
- Corridor: ~4°

**Mars Science Laboratory (6 km/s):**
- 70° sphere-cone aeroshell
- Tiled PICA TPS
- Peak g: ~13g
- Peak heating: ~200 W/cm² (~2,000 kW/m²)
- Guided entry (lifting body) widened corridor

**Galileo Probe (Jupiter, 47.4 km/s):**
- 45° sphere-cone
- Carbon phenolic TPS (thick)
- Peak g: 230g
- Peak heating: ~40,000 kW/m²
- The most demanding entry ever attempted. It survived.

---

## Aerocapture

What if you don't land? Fly through the atmosphere once, shed enough velocity
to enter orbit, and exit before you hit the ground. That's aerocapture — zero
propellant for orbital insertion.

The same physics applies, but you exit the atmosphere at exactly the right
velocity. Aerocapture saves ~2–5 km/s of ΔV compared to propulsive capture.
For Saturn missions, it's the difference between feasible and impossible.

The corridor is narrower than for landing — you need precise control of
the exit velocity. But the propellant savings are enormous.

---

## Key Insight

Entry corridor width determines everything about mission design. Narrow
corridor = high navigation accuracy required = more risk. Wide corridor =
more margin = safer mission. Jupiter has the narrowest corridor in the
solar system. Titan has one of the widest — dense atmosphere, low gravity,
gentle entry velocity.

The vehicle's ballistic coefficient is your design lever. Make it smaller
(more drag per kg) and the corridor widens. That's why entry capsules are
blunt — maximizing drag is more important than minimizing mass.

---

## Next

[Tutorial 4: Design a Mission →](tutorial-04-design-a-mission.md)
