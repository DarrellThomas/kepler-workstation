# Tutorial 4: Design a Mission

**Prerequisites:** Tutorials 1–3.

**You'll learn:** How to combine orbits, transfers, and entries into a complete
mission — and whether you need to make fuel along the way.

**Time:** 30 minutes.

---

## The Full Stack

A mission from Earth surface to another planet and back has six phases:

| Phase | What happens | ΔV source |
|-------|-------------|-----------|
| 1. Launch | Surface → parking orbit | Launch vehicle |
| 2. Departure | Parking orbit → interplanetary transfer | Upper stage |
| 3. Capture | Arrival → target orbit | Aerocapture or propulsive |
| 4. ISRU (optional) | Collect propellant for return | In-situ resources |
| 5. Return departure | Target orbit → Earth transfer | ISRU propellant |
| 6. Earth return | Arrival → landing | Aerocapture or direct entry |

Each phase has a ΔV cost. The sum determines the propellant mass. The
propellant mass determines whether the mission closes.

---

## Run the Saturn Mission

```bash
cd Applications/SaturnMission
g++ -std=c++11 -O2 -o launch_mission launch_mission.cpp
./launch_mission
```

This simulates Earth → Saturn → Earth with a 25-tonne payload (Orion-class
crew module + habitat), nuclear-thermal propulsion (Isp = 900s), and
aerocapture at both ends.

The ΔV budget:

| Phase | ΔV (km/s) | Propellant (kg) |
|-------|-----------|-----------------|
| Earth departure (LEO → Saturn) | ~7.6 | ~67,000 |
| Saturn capture (aerocapture) | 0 | 0 |
| Saturn departure (→ Earth) | ~7.6 | ~67,000 |
| Earth return (direct entry) | 0 | 0 |
| **Total** | **~15.2** | **~134,000** |

Without ISRU, you'd need 134 tonnes of propellant in LEO for 25 tonnes of
payload — a 6.4:1 ratio. With Starship-class heavy lift (~100t to LEO), you'd
need two launches just for the propellant, plus one for the payload.

---

## ISRU Changes Everything

If you can make propellant at the destination, you only carry the outbound fuel:

| Phase | ΔV (km/s) | Propellant (kg) |
|-------|-----------|-----------------|
| Earth departure | ~7.6 | ~67,000 |
| Saturn capture (aerocapture) | 0 | 0 |
| **ISRU at Titan** | — | Refuel for return |
| Saturn departure | ~7.6 | ~0 (ISRU) |
| Earth return | 0 | 0 |
| **Total from Earth** | **~7.6** | **~67,000** |

Now it fits in a single Starship launch. ISRU cuts the LEO mass nearly in half.

---

## Where to Get Fuel in the Saturn System

**Titan — the gas station:**

| Resource | How | Why |
|----------|-----|-----|
| Liquid methane | Submerged pump in Ligeia Mare | 450 kg/m³, 94 K, 1.47 bar |
| Water ice | Mine crust, melt, electrolyze | Titan IS 50% water ice by mass |
| Nitrogen | Ambient atmosphere at 1.47 bar | Fertilizer, habitat pressurization |

Titan has liquid methane on the surface. Not frozen. Not trapped in rock.
Liquid. Pump it directly into your tanks. Saturn arrival at Titan is 0 km/s
(aerocapture in Titan's atmosphere). Surface to orbit is 2.6 km/s.

**What doesn't work:**
- **Saturn atmospheric scoop.** Unlike Jupiter, Saturn's orbital velocity
  (26 km/s) is too slow for H2 self-ionization. The magnetic scoop design
  requires 42+ km/s to auto-ionize the incoming gas. Saturn is 60% slower.
- **Enceladus geysers.** The plume expands into vacuum. At any safe flyby
  altitude, density is ~10⁻⁸ kg/m³. A 50m-radius collector nets ~20 grams
  per pass. Millions of passes needed. Not practical.

---

## Gravity Assists

Sometimes the direct route costs too much fuel. A gravity assist steals
orbital energy from a planet.

**VEEJ (Venus-Earth-Earth-Jupiter):**
- Used by Galileo, Juno, and proposed for Bucket harvester
- Earth departure ΔV drops from 6.3 to 3.5 km/s
- Adds ~3 years to transit
- Free energy — the planet gives it, Newton collects

**How it works:** Fly behind a planet in its orbit → planet's gravity pulls you
forward → you gain velocity relative to the Sun. Fly ahead → you lose velocity.
The planet's orbit is trillions of times more massive than your spacecraft — the
planet doesn't notice. You get km/s of ΔV for zero propellant.

---

## The Rocket Equation

Everything flows from Tsiolkovsky:

```
ΔV = Isp · g₀ · ln(m₀ / mf)
```

Where:
- `Isp` is specific impulse (seconds) — how efficiently the engine uses propellant
- `g₀` is 9.80665 m/s² (standard gravity)
- `m₀` is initial mass
- `mf` is final mass (after burning propellant)

Rearranged for mass ratio:
```
m₀ / mf = exp(ΔV / (Isp · g₀))
```

At Isp = 900s (nuclear-thermal), ΔV = 7.6 km/s:
```
MR = exp(7600 / (900 · 9.80665)) = exp(0.861) = 2.37
```

You need 2.37 kg at departure for every 1 kg that arrives. The remaining
1.37 kg is propellant.

At Isp = 450s (chemical, methane/oxygen):
```
MR = exp(7600 / (450 · 9.80665)) = exp(1.722) = 5.60
```

You need 5.6 kg for every 1 kg delivered. This is why high-Isp propulsion
(nuclear, electric) matters for outer planet missions.

---

## Does Your Mission Close?

A mission "closes" if the initial mass in LEO is within your launch vehicle's
capability. Check:

1. Sum the ΔV for each phase
2. Compute the mass ratio: exp(ΣΔV / (Isp · g₀))
3. Initial mass = payload · MR · (1 + inert_fraction · (MR - 1))
4. If initial mass < launch vehicle capacity: mission closes

For Saturn with ISRU: 67 tonnes in LEO. Closes on a single Starship launch.
Without ISRU: 134 tonnes. Needs two launches and orbital assembly. ISRU is not
optional for outer planet missions — it's the single largest lever in the
mission design.

---

## Next

[Tutorial 5: Build a Constellation →](tutorial-05-build-a-constellation.md)
