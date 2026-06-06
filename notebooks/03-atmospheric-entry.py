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
# # Tutorial 3: Atmospheric Entry
#
# **Kepler Workstation** — Entry Corridor Analysis
#
# Simulates ballistic atmospheric entry for 7 planetary bodies. Sweeps
# entry flight path angles to find the viable corridor — not too shallow
# (skip out), not too steep (g-limit or heat-limit exceeded).

# %%
import sys
sys.path.insert(0, '../../bindings/python')
import matplotlib.pyplot as plt
from kepler.entry import simulate_entry, find_corridor

# %% [markdown]
# ## Entry Corridor for Every Body
#
# Run the corridor analysis for all 7 bodies. This takes ~10–20 seconds.

# %%
missions = [
    ('Earth',   120e3, 11.0e3),
    ('Mars',    125e3, 6.0e3),
    ('Jupiter', 500e3, 48.0e3),
    ('Venus',   150e3, 11.5e3),
    ('Titan',   500e3, 6.5e3),
    ('Neptune', 500e3, 23.0e3),
    ('Uranus',  500e3, 22.0e3),
]

print(f"{'Body':<12} {'v_entry':>8} {'Skip-out':>10} {'G-limit':>10} {'Q-limit':>10} {'Width':>8}")
print("-" * 62)

results = {}
for body, alt, vel in missions:
    cr = find_corridor(body, alt, vel)
    results[body] = cr
    print(f"{body:<12} {vel/1000:>6.1f}  {cr['skip_out_deg']:>8.1f}° {cr['g_limit_deg']:>8.1f}° {cr['q_limit_deg']:>8.1f}° {cr['corridor_width']:>6.1f}°")

# %% [markdown]
# ## Corridor Visualization — Earth
#
# Plot the complete corridor sweep: skip-out boundary at the shallow end,
# g-limit and heat-rate boundaries at the steep end.

# %%
cr = results['Earth']
angles = [s[0] for s in cr['sweep']]
outcomes = [s[1] for s in cr['sweep']]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Peak g-load vs entry angle
g_vals = [s[1]['g_max'] for s in cr['sweep']]
ax1.plot(angles, g_vals, 'b-', linewidth=2)
ax1.axhline(10, color='red', linestyle='--', alpha=0.7, label='10g limit')
ax1.axvline(cr['skip_out_deg'], color='green', linestyle='--', alpha=0.7, label=f'Skip-out: {cr["skip_out_deg"]:.1f}°')
ax1.set_xlabel('Entry Flight Path Angle (°)')
ax1.set_ylabel('Peak g-load')
ax1.set_title('Peak Deceleration vs Entry Angle — Earth')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Heat rate vs entry angle
q_vals = [s[1]['q_max_kw'] for s in cr['sweep']]
ax2.plot(angles, q_vals, 'r-', linewidth=2)
ax2.axhline(1500, color='red', linestyle='--', alpha=0.7, label='1500 kW/m² limit')
ax2.axvline(cr['skip_out_deg'], color='green', linestyle='--', alpha=0.7, label=f'Skip-out: {cr["skip_out_deg"]:.1f}°')
ax2.set_xlabel('Entry Flight Path Angle (°)')
ax2.set_ylabel('Peak Heat Rate (kW/m²)')
ax2.set_title('Peak Heating vs Entry Angle — Earth')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()

print(f"Earth entry corridor: {cr['skip_out_deg']:.1f}° to {cr['g_limit_deg']:.1f}°  (width: {cr['corridor_width']:.1f}°)")

# %% [markdown]
# ## Trajectory Comparison — Good vs Bad Entry
#
# Compare a successful entry (-6°) vs one that skips out (-3.5°).

# %%
good = simulate_entry('earth', 120e3, 11.0e3, -6.0)
bad  = simulate_entry('earth', 120e3, 11.0e3, -3.5)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(good['trajectory']['vel_kms'], good['trajectory']['alt_km'], 'g-', linewidth=2, label='-6.0° (capture)')
ax1.plot(bad['trajectory']['vel_kms'], bad['trajectory']['alt_km'], 'r-', linewidth=2, label='-3.5° (skip-out)')
ax1.set_xlabel('Velocity (km/s)')
ax1.set_ylabel('Altitude (km)')
ax1.set_title('Entry Trajectories — Earth')
ax1.invert_xaxis()
ax1.grid(True, alpha=0.3)
ax1.legend()

ax2.plot(good['trajectory']['vel_kms'], good['trajectory']['g_load'], 'g-', linewidth=2, label='-6.0°')
ax2.plot(bad['trajectory']['vel_kms'], bad['trajectory']['g_load'], 'r-', linewidth=2, label='-3.5°')
ax2.set_xlabel('Velocity (km/s)')
ax2.set_ylabel('g-load')
ax2.set_title('Deceleration Profile — Earth')
ax2.invert_xaxis()
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Next
#
# [Notebook 4: Design a Mission →](04-design-a-mission.py)
