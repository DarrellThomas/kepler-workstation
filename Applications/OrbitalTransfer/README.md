# Orbital Transfer — Interplanetary Mission Design

Earth-to-Mars transfer analysis using real JPL DE440 ephemeris and a
Lambert transfer solver. Part of Kepler Workstation.

## What's Here

| File | What it does |
|------|-------------|
| `hohmann_demo.cpp` | Earth→Mars Hohmann transfer: ideal (textbook) vs real (ephemeris) |
| `porkchop_plot.cpp` | Full departure/TOF sweep, CSV output, best launch windows |
| `README.md` | This file |

## Quick Start

```bash
# 1. Hohmann transfer demo (works without ephemeris for ideal portion)
g++ -std=c++11 -O2 -o hohmann_demo hohmann_demo.cpp
./hohmann_demo

# With real ephemeris:
g++ -std=c++11 -O2 -DHAS_EPHEMERIS -o hohmann_demo hohmann_demo.cpp
./hohmann_demo

# 2. Porkchop plot (recommended: use ephemeris for real data)
g++ -std=c++11 -O2 -o porkchop_plot porkchop_plot.cpp
./porkchop_plot
# Outputs porkchop_earth_mars.csv + prints top 5 windows

# With real ephemeris:
g++ -std=c++11 -O2 -DHAS_EPHEMERIS -o porkchop_plot porkchop_plot.cpp
./porkchop_plot
```

The DE440 ephemeris file is not in the repo (.gitignored, ~100 MB).
Download it once:

```bash
cd ../../Ephemeris
wget https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440
```

## hohmann_demo — What You'll See

**Ideal Hohmann** (circular, coplanar — the textbook answer):
- Transfer orbit: 1.000–1.524 AU, eccentricity ~0.21
- C3 departure: ~8.6 km²/s²
- ΔV total: ~5.7 km/s (from 300 km parking orbits)
- Time of flight: ~259 days (8.5 months)

**Real transfer** (JPL DE440 + Lambert solver):
- Solar system isn't circular, coplanar, or static
- C3 varies with the Earth-Mars relative geometry
- Best windows cluster near opposition (every ~780 days)
- Real ΔV typically 10–20% above ideal due to inclination and eccentricity

## porkchop_plot — What You'll See

Sweeps departure dates (2026–2036, 20-day steps) and times of flight
(50–450 days, 10-day steps), solving Lambert's problem for each.

**Outputs:**
- `porkchop_earth_mars.csv` — every valid transfer with C3, v_inf, ΔV
- Top 5 windows printed to stdout (minimum C3)
- Summary statistics (C3 range, TOF range, valid window count)

**Plotting the porkchop** (Python):
```python
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv('porkchop_earth_mars.csv')
plt.tricontourf(df['departure_year'], df['tof_days'], df['c3_depart'],
                levels=20, cmap='viridis_r')
plt.colorbar(label='C3 (km²/s²)')
plt.xlabel('Departure Year')
plt.ylabel('Time of Flight (days)')
plt.title('Earth→Mars Porkchop Plot')
plt.show()
```

## How It Works

1. **Lambert solver** (`Ephemeris/lambert.hpp`) — given two positions and a
   time of flight, computes the transfer orbit's velocity vectors. Universal
   variable formulation from Bate, Mueller & White (1971) with Battin's
   improvement for near-180° transfers.

2. **JPL DE440 ephemeris** (`Ephemeris/jpl_ephemeris.hpp`) — Chebyshev
   polynomial interpolation of JPL's planetary positions. Covers 1550–2650 CE.
   The gold standard for interplanetary mission design.

3. **Porkchop plot** — NASA's standard visualization for launch window
   analysis. Sweep departure date (x-axis) vs time of flight (y-axis),
   color by C3 or ΔV. The "porkchop" shape comes from the synodic period
   of the two planets.

## The Math

**Hohmann transfer** (circular, coplanar):
- Transfer semi-major axis: a = (r₁ + r₂) / 2
- Vis-viva: v = √(μ · (2/r − 1/a))
- Time of flight: TOF = π · √(a³/μ)
- C3 = (v_transfer − v_planet)²

**Lambert's problem** (general case):
- Given r₁, r₂, TOF, and μ — find v₁, v₂
- Solved via Newton-Raphson on the universal variable z
- Works for elliptic, parabolic, and hyperbolic transfers

## References

- Bate, Mueller & White — *Fundamentals of Astrodynamics* (1971)
- Battin — *An Introduction to the Mathematics and Methods of Astrodynamics* (1999)
- JPL — DE440 Planetary Ephemeris (Park et al. 2021, AJ 161:105)
- NASA — *Basics of Space Flight* (online tutorial)

## Next Steps

- Venus and Mars gravity assist trajectories
- Earth→Jupiter transfer (already in JupiterHarvester)
- Multi-revolution Lambert solver
- Python bindings for interactive notebook use
