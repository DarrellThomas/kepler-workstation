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
# # Tutorial 2: Transfer to Mars
# 
# **Kepler Workstation** — Earth → Mars Porkchop Plot
# 
# Computes the optimal transfer windows between Earth and Mars using
# a Lambert solver and displays the classic "porkchop" contour plot.

# %%
import sys
sys.path.insert(0, '../../bindings/python')
import numpy as np
import matplotlib.pyplot as plt
from kepler.constants import GM_SUN, GM_EARTH, GM_MARS, SMA_EARTH, SMA_MARS, R_EARTH, R_MARS
from kepler.orbit import hohmann_transfer
from kepler.transfer import porkchop_sweep

# %% [markdown]
# ## Ideal Hohmann Transfer
# 
# First, compute the textbook Hohmann (circular, coplanar) as a reference.

# %%
h = hohmann_transfer(SMA_EARTH, SMA_MARS, GM_SUN, GM_EARTH, GM_MARS,
                      R_EARTH, R_MARS, alt_park=300e3)

print("Ideal Hohmann — Earth → Mars")
print(f"  C3 departure:   {h['c3_depart']/1e6:.2f} km²/s²")
print(f"  v_inf departure: {h['v_inf_depart']/1000:.2f} km/s")
print(f"  v_inf arrival:   {h['v_inf_arrive']/1000:.2f} km/s")
print(f"  ΔV departure:    {h['dv_depart']/1000:.2f} km/s")
print(f"  ΔV arrival:      {h['dv_arrive']/1000:.2f} km/s")
print(f"  ΔV total:        {h['dv_total']/1000:.2f} km/s")
print(f"  Time of flight:  {h['tof_days']:.0f} days")

# %% [markdown]
# ## Porkchop Plot — Earth → Mars 2026–2036
# 
# Sweep departure dates and times of flight, solving Lambert's problem
# for each combination. The contour shape gives the plot its name.

# %%
departure_years = np.linspace(2026, 2036, 80)
tof_days = np.linspace(50, 450, 60)

c3_grid, vinf_grid = porkchop_sweep(None, None, None, None, tof_days, departure_years)

# Mask invalid entries
c3_masked = np.ma.masked_invalid(c3_grid)

fig, ax = plt.subplots(figsize=(12, 7))
levels = np.arange(0, 100, 2)
cs = ax.contourf(departure_years, tof_days, c3_masked, levels=levels,
                  cmap='viridis_r', extend='max')
c = ax.contour(departure_years, tof_days, c3_masked, levels=[10, 20, 30, 50, 80],
                colors='white', linewidths=0.5, alpha=0.5)
ax.clabel(c, inline=True, fontsize=8, fmt='%.0f')

# Mark the ideal Hohmann TOF
ax.axhline(259, color='red', linestyle='--', alpha=0.7, label=f'Ideal Hohmann TOF ({h["tof_days"]:.0f} days)')

cbar = plt.colorbar(cs, ax=ax, label='C3 (km²/s²)')
ax.set_xlabel('Departure Year')
ax.set_ylabel('Time of Flight (days)')
ax.set_title('Earth → Mars Porkchop Plot   (Kepler Workstation)')
ax.legend(loc='upper right')
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Best Windows
# 
# Find the global C3 minimum across the entire sweep.

# %%
min_idx = np.unravel_index(np.nanargmin(c3_grid), c3_grid.shape)
best_tof = tof_days[min_idx[0]]
best_year = departure_years[min_idx[1]]
best_c3 = c3_grid[min_idx]

print("Best transfer window:")
print(f"  Departure: {best_year:.2f}")
print(f"  TOF:       {best_tof:.0f} days")
print(f"  C3:        {best_c3:.2f} km²/s²")
print(f"  Arrival:   {best_year + best_tof/365.25:.2f}")

# %% [markdown]
# ## Launch Windows Over Time
# 
# Every ~780 days (Earth-Mars synodic period), a new launch window opens.
# Plot the minimum C3 for each departure year.

# %%
min_c3_per_year = np.nanmin(c3_grid, axis=0)

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(departure_years, min_c3_per_year, 'b-', linewidth=1.5)
ax.set_xlabel('Departure Year')
ax.set_ylabel('Minimum C3 (km²/s²)')
ax.set_title('Best Earth→Mars C3 by Departure Year')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Next
# 
# [Notebook 3: Atmospheric Entry →](03-atmospheric-entry.py)
