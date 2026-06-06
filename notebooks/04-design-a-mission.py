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
# # Tutorial 4: Design a Mission
#
# **Kepler Workstation** — Earth → Saturn → Earth
#
# Full mission design: launch, interplanetary transfer, ISRU analysis,
# and the rocket equation. Does the mission close?

# %%
import sys
sys.path.insert(0, '../../bindings/python')
import numpy as np
import matplotlib.pyplot as plt
from kepler.constants import (
    GM_SUN, GM_EARTH, GM_SATURN, G0_TITAN,
    R_EARTH, R_SATURN, SMA_EARTH, SMA_SATURN,
)
from kepler.orbit import hohmann_transfer, rocket_equation

# %% [markdown]
# ## Mission Parameters
#
# Configure your mission. Change these and re-run to see how the
# feasibility changes.

# %%
# Vehicle
payload_kg   = 25000   # Orion-class crew module + habitat
isp_s        = 900     # Nuclear-thermal (H2 propellant)

# Parking orbits
alt_leo_km   = 300
alt_sat_km   = 1000

print(f"Mission: {payload_kg/1000:.0f} t payload, Isp={isp_s}s NTR")
print(f"LEO parking: {alt_leo_km} km, Saturn parking: {alt_sat_km} km")

# %% [markdown]
# ## Phase 1: Earth → Saturn Transfer
#
# Hohmann transfer from Earth to Saturn (circular-coplanar approximation).

# %%
h_out = hohmann_transfer(SMA_EARTH, SMA_SATURN, GM_SUN, GM_EARTH, GM_SATURN,
                          R_EARTH, R_SATURN, alt_park=alt_leo_km*1000)

print("Outbound Transfer — Earth → Saturn")
print(f"  C3 departure:      {h_out['c3_depart']/1e6:.1f} km²/s²")
print(f"  v_inf departure:   {h_out['v_inf_depart']/1000:.2f} km/s")
print(f"  v_inf arrival:     {h_out['v_inf_arrive']/1000:.2f} km/s")
print(f"  ΔV departure:      {h_out['dv_depart']/1000:.2f} km/s")
print(f"  ΔV arrival:        {h_out['dv_arrive']/1000:.2f} km/s")
print(f"  Time of flight:    {h_out['tof_days']:.0f} days ({h_out['tof_days']/365.25:.2f} years)")

# Propellant for departure
depart = rocket_equation(h_out['dv_depart'], isp_s, payload_kg)
print(f"\n  Propellant (departure): {depart['propellant_kg']/1000:.0f} t")
print(f"  Mass in LEO:            {depart['initial_kg']/1000:.0f} t")

# %% [markdown]
# ## Phase 2: Saturn Arrival
#
# Saturn aerocapture — free ΔV. The atmosphere does the work.

# %%
# Propulsive capture for comparison
dv_capture_prop = h_out['dv_arrive']
capture_prop = rocket_equation(dv_capture_prop, isp_s,
                                payload_kg + depart['propellant_kg'])

print("Capture options at Saturn:")
print(f"  Propulsive capture ΔV:  {dv_capture_prop/1000:.2f} km/s")
print(f"  Propellant (propulsive): {capture_prop['propellant_kg']/1000:.0f} t")
print("  Aerocapture:            0 km/s — FREE")
print(f"\nSavings from aerocapture: {capture_prop['propellant_kg']/1000:.0f} t propellant")

# %% [markdown]
# ## Phase 3: ISRU at Saturn
#
# Can we refuel for the return trip?

# %%
# Saturn magnetic scoop analysis
v_peri = np.sqrt(GM_SATURN / (R_SATURN + 500e3))
h2_mass = 3.34e-27  # kg
ke_j = 0.5 * h2_mass * v_peri**2
ke_ev = ke_j / 1.602e-19

print("Saturn Magnetic Scoop Feasibility:")
print(f"  Orbital velocity at 500 km: {v_peri/1000:.1f} km/s")
print(f"  H2 kinetic energy:          {ke_ev:.1f} eV")
print("  Ionization threshold:       15.4 eV")
print(f"  SCOOP WORKS:                {'YES' if ke_ev >= 15.4 else 'NO — too slow'}")
print()
print("Titan ISRU (always available):")
print(f"  Surface gravity:  {G0_TITAN:.2f} m/s²")
print("  Atmosphere:       1.47 bar N2/CH4 at 94 K")
print("  Methane lakes:    pump liquid CH4 directly")
print("  Surface → orbit:  ~2.6 km/s")

# %% [markdown]
# ## Phase 4: Return Transfer
#
# Saturn → Earth Hohmann (same ΔV, reversed).

# %%
h_ret = hohmann_transfer(SMA_SATURN, SMA_EARTH, GM_SUN, GM_SATURN, GM_EARTH,
                          R_SATURN, R_EARTH, alt_park=alt_sat_km*1000)

print("Return Transfer — Saturn → Earth")
print(f"  ΔV departure:      {h_ret['dv_depart']/1000:.2f} km/s")
print(f"  v_inf arrival:     {h_ret['v_inf_arrive']/1000:.2f} km/s")
print("  ΔV arrival:        0 km/s (direct entry)")
print(f"  Time of flight:    {h_ret['tof_days']:.0f} days ({h_ret['tof_days']/365.25:.2f} years)")

# %% [markdown]
# ## Mission Summary
#
# Does it close? Compare with and without ISRU.

# %%
# Without ISRU: carry all propellant from Earth
dv_no_isru = h_out['dv_depart'] + h_ret['dv_depart']  # aerocapture at both ends = 0
mass_no_isru = rocket_equation(dv_no_isru, isp_s, payload_kg)

# With ISRU: only carry outbound propellant
dv_isru = h_out['dv_depart']
mass_isru = rocket_equation(dv_isru, isp_s, payload_kg)

print("=" * 60)
print("  MISSION SUMMARY — Earth → Saturn → Earth")
print("=" * 60)
print()
print(f"  {'':<30} {'Without ISRU':>15} {'With ISRU':>15}")
print(f"  {'Outbound ΔV':<30} {h_out['dv_depart']/1000:>13.1f} km/s {h_out['dv_depart']/1000:>13.1f} km/s")
print(f"  {'Saturn capture':<30} {'aerocapture (free)':>15} {'aerocapture (free)':>15}")
print(f"  {'Return ΔV':<30} {h_ret['dv_depart']/1000:>13.1f} km/s {'0 (ISRU refuel)':>15}")
print(f"  {'Earth return':<30} {'direct entry (free)':>15} {'direct entry (free)':>15}")
print(f"  {'Total ΔV':<30} {dv_no_isru/1000:>13.1f} km/s {dv_isru/1000:>13.1f} km/s")
print()
print(f"  {'Propellant':<30} {mass_no_isru['propellant_kg']/1000:>13.0f} t {mass_isru['propellant_kg']/1000:>13.0f} t")
print(f"  {'Inert structure':<30} {mass_no_isru['inert_kg']/1000:>13.0f} t {mass_isru['inert_kg']/1000:>13.0f} t")
print(f"  {'Payload':<30} {payload_kg/1000:>13.0f} t {payload_kg/1000:>13.0f} t")
print(f"  {'INITIAL MASS IN LEO':<30} {mass_no_isru['initial_kg']/1000:>13.0f} t {mass_isru['initial_kg']/1000:>13.0f} t")
print()
starship_capacity = 100e3
print(f"  Starship capacity: {starship_capacity/1000:.0f} t")
print(f"  Without ISRU: {'CLOSES' if mass_no_isru['initial_kg'] < starship_capacity else 'NEEDS MULTIPLE LAUNCHES'}")
print(f"  With ISRU:    {'CLOSES' if mass_isru['initial_kg'] < starship_capacity else 'NEEDS MULTIPLE LAUNCHES'}")

# %% [markdown]
# ## Isp Sensitivity
#
# How does the propellant mass change with Isp? Chemical vs nuclear vs electric.

# %%
isps = np.linspace(300, 3000, 50)
prop_masses = np.array([rocket_equation(dv_no_isru, isp, payload_kg)['propellant_kg']/1000 for isp in isps])

fig, ax = plt.subplots(figsize=(12, 5))
ax.semilogy(isps, prop_masses, 'b-', linewidth=2)
ax.axhline(payload_kg/1000, color='green', linestyle='--', alpha=0.5, label=f'Payload ({payload_kg/1000:.0f} t)')
ax.axhline(100, color='red', linestyle='--', alpha=0.5, label='Starship capacity (100 t)')
ax.axvline(450, color='orange', linestyle=':', alpha=0.7, label='Chemical (450s)')
ax.axvline(900, color='purple', linestyle=':', alpha=0.7, label='NTR (900s)')
ax.axvline(3000, color='brown', linestyle=':', alpha=0.7, label='Hall thruster (3000s)')
ax.set_xlabel('Specific Impulse (seconds)')
ax.set_ylabel('Propellant Mass (tonnes)')
ax.set_title('Propellant Mass vs Isp — Saturn Round Trip')
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Next
#
# [Notebook 5: Build a Constellation →](05-build-a-constellation.py)
