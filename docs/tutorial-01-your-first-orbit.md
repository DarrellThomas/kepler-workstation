# Tutorial 1: Your First Orbit

**Audience:** Anyone who can run a terminal. No aerospace background needed.

**You'll learn:** What an orbit actually is, why it works, and how to compute one.

**Time:** 15 minutes.

---

## What Is an Orbit?

An orbit is falling and missing the ground.

Isaac Newton drew the picture in 1687: fire a cannonball horizontally from a
mountaintop. At low speed it falls to Earth. At higher speed it goes further
before hitting. At a certain speed — *orbital velocity* — the Earth curves away
beneath the cannonball at the same rate the cannonball falls. It keeps falling
forever. That's an orbit.

The same physics governs the ISS, the Moon, and every satellite. Gravity pulls
down. Horizontal velocity keeps you missing.

---

## The Vis-Viva Equation

For a circular orbit around a planet:

```
v = √(GM / r)
```

Where:
- `v` is orbital velocity (m/s)
- `GM` is the planet's gravitational parameter (m³/s²). Earth: 3.986×10¹⁴
- `r` is the distance from the planet's center (m). Earth radius + altitude

For the ISS at 408 km altitude:
```
r = 6,371,000 + 408,000 = 6,779,000 m
v = √(3.986×10¹⁴ / 6,779,000) = 7,668 m/s = 7.67 km/s
```

That's 27,600 km/h. The ISS circles Earth every 93 minutes.

---

## Compute It Yourself

No spacecraft required. Open a terminal:

```bash
cd PlanetaryModels
g++ -std=c++11 -o orbit orbit.cpp
```

Create `orbit.cpp`:

```cpp
#include <cstdio>
#include <cmath>

int main() {
    double GM   = 3.986004418e14;  // Earth, m³/s²
    double R    = 6378137.0;       // Earth radius, m
    double alt  = 408000.0;        // ISS altitude, m
    double r    = R + alt;

    double v    = sqrt(GM / r);
    double T    = 2.0 * M_PI * r / v;  // period in seconds

    printf("Orbit: %.0f km altitude\n", alt/1000);
    printf("Velocity:   %.2f km/s  (%5.0f km/h)\n", v/1000, v*3.6/1000);
    printf("Period:     %.1f minutes\n", T/60);

    return 0;
}
```

```
$ ./orbit
Orbit: 408 km altitude
Velocity:   7.67 km/s  (27602 km/h)
Period:     92.5 minutes
```

---

## What Changes With Altitude?

Higher orbits are slower. Counterintuitive but true — vis-viva says v decreases
with √(1/r). But higher orbits have more *energy* because potential energy
increases faster than kinetic energy decreases.

| Orbit | Altitude | Velocity | Period |
|-------|----------|----------|--------|
| ISS | 408 km | 7.67 km/s | 93 min |
| Hubble | 540 km | 7.59 km/s | 96 min |
| GPS | 20,200 km | 3.87 km/s | 12 hrs |
| Geostationary | 35,786 km | 3.07 km/s | 24 hrs |
| Moon | 384,400 km | 1.02 km/s | 27 days |

Run the numbers yourself — change `alt` in the code and recompile.

---

## Kepler's Third Law

The orbital period squared is proportional to the semi-major axis cubed:

```
T² ∝ a³
```

Or more precisely:
```
T = 2π √(a³ / GM)
```

Plug in Earth's GM and the Moon's distance (384,400 km) and you get 27.3 days.
This is how we know the mass of every planet with a moon — watch the moon orbit,
measure the period and distance, solve for GM.

---

## Why This Matters

Every mission in Kepler Workstation starts here. The Hohmann transfer (Tutorial
2) is two tangential burns between circular orbits. The entry corridor (Tutorial
3) is what happens when an orbit intersects an atmosphere. The constellation
(Tutorial 5) is thousands of orbits stitched together across decades.

Orbital mechanics is not rocket science. Rocket science is rocket science.
Orbital mechanics is geometry plus one equation.

---

## Next

[Tutorial 2: Transfer to Mars →](tutorial-02-transfer-to-mars.md)
