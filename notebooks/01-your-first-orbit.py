# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Tutorial 1: Your First Orbit
# 
# **Kepler Workstation** — Interactive Orbital Mechanics
# 
# This notebook computes circular orbits around any planet. Change the
# altitude and see what happens to velocity and period.

# %%
import sys
sys.path.insert(0, '../../bindings/python')
import numpy as np
import matplotlib.pyplot as plt
from kepler.orbit import circular_velocity, circular_period
from kepler.constants import *

# %% [markdown]
# ## ISS Orbit
# 
# The International Space Station orbits at ~408 km. Compute its velocity
# and period using the vis-viva equation.

# %%
alt = 408e3  # meters
r = R_EARTH + alt
v = circular_velocity(GM_EARTH, r)
T = circular_period(GM_EARTH, r)

print(f"ISS at {alt/1000:.0f} km altitude:")
print(f"  Velocity: {v/1000:.2f} km/s  ({v*3.6/1000:.0f} km/h)")
print(f"  Period:   {T/60:.1f} minutes  ({T/3600:.2f} hours)")

# %% [markdown]
# ## Orbit at Any Altitude
# 
# Sweep altitudes from 200 km (LEO minimum) to geostationary (35,786 km)
# and see how velocity changes. Higher orbits are slower — counterintuitive
# but true.

# %%
alts_km = np.linspace(200, 36000, 200)
alts_m = alts_km * 1000
velocities = np.array([circular_velocity(GM_EARTH, R_EARTH + a) / 1000 for a in alts_m])
periods = np.array([circular_period(GM_EARTH, R_EARTH + a) / 3600 for a in alts_m])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(alts_km, velocities, 'b-', linewidth=2)
ax1.axvline(408, color='red', linestyle='--', alpha=0.5, label='ISS (408 km)')
ax1.axvline(35786, color='green', linestyle='--', alpha=0.5, label='GEO (35,786 km)')
ax1.set_xlabel('Altitude (km)'); ax1.set_ylabel('Velocity (km/s)')
ax1.set_title('Orbital Velocity vs Altitude'); ax1.grid(True, alpha=0.3); ax1.legend()

ax2.plot(alts_km, periods, 'r-', linewidth=2)
ax2.axhline(24, color='green', linestyle='--', alpha=0.5, label='24 hours (GEO)')
ax2.set_xlabel('Altitude (km)'); ax2.set_ylabel('Period (hours)')
ax2.set_title('Orbital Period vs Altitude'); ax2.grid(True, alpha=0.3); ax2.legend()

plt.tight_layout(); plt.show()

# %% [markdown]
# ## Compare Every Planet
# 
# What's the circular orbit velocity at 500 km altitude around each body?

# %%
bodies = {
    'Earth':   (GM_EARTH,   R_EARTH,   G0_EARTH),
    'Mars':    (GM_MARS,    R_MARS,    G0_MARS),
    'Jupiter': (GM_JUPITER, R_JUPITER, G0_JUPITER),
    'Venus':   (GM_VENUS,   R_VENUS,   G0_VENUS),
    'Titan':   (GM_TITAN,   R_TITAN,   G0_TITAN),
    'Neptune': (GM_NEPTUNE, R_NEPTUNE, G0_NEPTUNE),
    'Uranus':  (GM_URANUS,  R_URANUS,  G0_URANUS),
    'Moon':    (GM_MOON,    R_MOON,    G0_MOON),
}

alt_km = 500
print(f"{'Body':<12} {'v (km/s)':>10} {'Period':>12} {'g (m/s²)':>10}")
print("-" * 48)
for name, (gm, r, g0) in bodies.items():
    v = circular_velocity(gm, r + alt_km*1000) / 1000
    T = circular_period(gm, r + alt_km*1000) / 60
    print(f"{name:<12} {v:>10.3f} {T:>9.1f} min {g0:>10.2f}")

# %% [markdown]
# ## Next
# 
# [Notebook 2: Transfer to Mars →](02-transfer-to-mars.py)
