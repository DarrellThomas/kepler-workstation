# Tutorial 2: Transfer to Mars

**Prerequisite:** [Tutorial 1: Your First Orbit](tutorial-01-your-first-orbit.md)

**You'll learn:** How to get from Earth to Mars using the least energy possible.

**Time:** 20 minutes.

---

## The Problem

You're in a circular orbit around the Sun at 1 AU (Earth). You want to reach a
circular orbit at 1.524 AU (Mars). How do you do it with the minimum fuel?

Answer: the Hohmann transfer — an elliptical orbit that touches both circles.

---

## The Hohmann Transfer

Walter Hohmann solved this in 1925. Fire your engine once to enter an elliptical
transfer orbit. Coast for half an orbit. Fire again to circularize at
the target.

```
      ┌─────────────────────────────────────────┐
      │            Transfer Ellipse              │
      │  ┌─────┐                    ┌─────┐     │
      │  │ Sun │═══════════════════│ Mars│     │
      │  └─────┘                    └─────┘     │
      │  Earth                                  │
      └─────────────────────────────────────────┘
```

**The math:**

Transfer orbit semi-major axis:
```
a_transfer = (r_earth + r_mars) / 2
```

Vis-viva at departure (perihelion of transfer orbit):
```
v_depart = √(GM_sun · (2/r_earth - 1/a_transfer))
```

Earth's orbital velocity:
```
v_earth = √(GM_sun / r_earth)
```

The extra velocity needed (v_infinity):
```
v_inf = v_depart - v_earth
C3 = v_inf²
```

Time of flight (half the transfer orbit period):
```
TOF = π · √(a_transfer³ / GM_sun)
```

---

## Run the Numbers

```bash
cd Applications/OrbitalTransfer
g++ -std=c++11 -O2 -o hohmann_demo hohmann_demo.cpp
./hohmann_demo
```

The ideal Hohmann for Earth→Mars:

| Parameter | Value |
|-----------|-------|
| Transfer a | 1.262 AU |
| C3 departure | ~8.7 km²/s² |
| ΔV departure | ~3.6 km/s (from 300 km LEO) |
| ΔV arrival | ~2.1 km/s |
| ΔV total | ~5.7 km/s |
| Time of flight | ~259 days (8.5 months) |

259 days. That's why Mars missions launch in specific windows — the planets have
to be in the right positions when you arrive. Mars needs to be about 44° ahead
of Earth at departure, since Earth moves faster and "catches up" during the
transfer.

---

## Real Transfers Aren't Perfect

The ideal Hohmann assumes circular, coplanar orbits. Reality:

- Earth's orbit eccentricity: 0.017 (nearly circular)
- Mars' orbit eccentricity: 0.093 (noticeably elliptical)
- Mars' inclination to ecliptic: 1.85°
- Both orbits precess slowly

These imperfections mean the real C3 varies by 10–30% depending on the specific
departure date. A porkchop plot shows this.

---

## The Porkchop Plot

```bash
cd Applications/OrbitalTransfer
g++ -std=c++11 -O2 -o porkchop_plot porkchop_plot.cpp
./porkchop_plot
```

This sweeps departure dates (2026–2036) and times of flight (50–450 days),
solving Lambert's problem for each combination. Output:

- `porkchop_earth_mars.csv` — every valid transfer
- Top 5 windows printed to terminal
- Best C3: ~8.7 km²/s² at TOF ~260 days (the Hohmann minimum)

The "porkchop" name comes from the contour shape: low C3 in the center (the
porkchop), high C3 forming the "bone" around it. Launch windows appear every
~780 days (the Earth-Mars synodic period).

Plot it with Python:
```python
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv('porkchop_earth_mars.csv')
plt.tricontourf(df['departure_year'], df['tof_days'], df['c3_depart'],
                levels=20, cmap='viridis_r')
plt.colorbar(label='C3 (km²/s²)')
plt.xlabel('Departure Year'); plt.ylabel('Time of Flight (days)')
plt.show()
```

---

## Lambert's Problem

The Hohmann transfer works for circular orbits. For arbitrary start and end
positions, we need Lambert's problem: given two position vectors and a time of
flight, find the transfer orbit.

The universal variable formulation (Bate, Mueller & White, 1971) solves this
via Newton-Raphson iteration. It works for elliptic, parabolic, and hyperbolic
transfers. The code is in `Ephemeris/lambert.hpp`.

This is what NASA uses for preliminary mission design. The same algorithm that
routed Voyager, Cassini, and Perseverance.

---

## Next

[Tutorial 3: Atmospheric Entry →](tutorial-03-atmospheric-entry.md)
