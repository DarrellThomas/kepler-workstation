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
# # Tutorial 5: Build a Constellation
# 
# **Kepler Workstation** — Jupiter Helium Harvester Fleet
# 
# Simulates a fleet of 100 autonomous "Bucket" harvesters operating between
# Earth and Jupiter. Tracks launches, attrition, deliveries, and economics
# over 36 years.

# %%
import numpy as np
import matplotlib.pyplot as plt

# %% [markdown]
# ## Fleet Parameters
# 
# One Bucket: 8 tonnes dry mass, 20T superconducting magnet, Hall thruster.
# Fleet: 10 buckets per Starship launch, every 399 days (synodic period).

# %%
# Single Bucket
dry_mass_kg      = 8000
tank_cap_kg      = 33000
isp_s            = 3000  # Hall thruster
loop_years       = 10    # Earth→Jupiter→Earth
he_per_loop_kg   = 8200  # helium delivered per loop
h2_per_loop_kg   = 13400 # hydrogen delivered per loop

# Fleet
buckets_per_launch  = 10
synodic_period_days = 399
launch_years        = 36
attrition_per_loop  = 0.10  # 10% per loop

print("The Bucket — Jupiter Helium Harvester")
print(f"  Dry mass:       {dry_mass_kg/1000:.0f} t")
print("  Magnet:         20T HTS superconducting")
print(f"  Thruster:       Hall effect, Isp={isp_s}s")
print(f"  Loop duration:  {loop_years} years")
print(f"  He per loop:    {he_per_loop_kg/1000:.1f} t")
print(f"  H2 per loop:    {h2_per_loop_kg/1000:.1f} t")

# %% [markdown]
# ## Launch Schedule
# 
# Launch 10 Buckets per synodic window. After 10 windows (~11 years),
# the first Buckets start returning.

# %%
window_days = synodic_period_days
n_windows = int(launch_years * 365.25 / window_days)
launch_jd = np.arange(n_windows) * window_days
launch_year = 2026 + launch_jd / 365.25

buckets_launched = np.zeros(n_windows)
buckets_launched[:10] = buckets_per_launch  # first 10 windows at full rate
# After that, replacement rate
for i in range(10, n_windows):
    buckets_launched[i] = buckets_per_launch

total_launched = np.sum(buckets_launched)

print(f"Launch schedule: {total_launched:.0f} Buckets over {n_windows} windows ({launch_years} years)")
print(f"Launch cadence:  {buckets_per_launch} per Starship, every {synodic_period_days} days")

# %% [markdown]
# ## Fleet Simulation
# 
# Track every Bucket through its lifecycle. Each loop takes 10 years.
# Attrition removes 10% per loop. Replacements keep the fleet at steady state.

# %%
class Bucket:
    def __init__(self, launch_year):
        self.launch_year = launch_year
        self.loops_completed = 0
        self.active = True
        self.he_delivered = 0
        self.h2_delivered = 0
        self.return_year = launch_year + loop_years  # first return
    
    def complete_loop(self):
        if not self.active:
            return
        self.loops_completed += 1
        self.he_delivered += he_per_loop_kg
        self.h2_delivered += h2_per_loop_kg
        self.return_year += loop_years
        # Attrition check
        if np.random.random() < attrition_per_loop:
            self.active = False

# Run simulation over 36 years
np.random.seed(42)
buckets = []
for i, ly in enumerate(launch_year):
    for _ in range(int(buckets_launched[i])):
        buckets.append(Bucket(ly))

# Simulate year by year
sim_years = np.arange(2026, 2026 + launch_years + 1)
annual_he = np.zeros(len(sim_years))
active_fleet = np.zeros(len(sim_years))
cumulative_he = np.zeros(len(sim_years))

for y_idx, year in enumerate(sim_years):
    for b in buckets:
        if not b.active:
            continue
        # Check if this bucket returns this year
        if abs(b.return_year - year) < 0.5 and b.loops_completed == 0:
            b.complete_loop()
        elif abs(b.return_year - year) < 0.5 and b.active:
            b.complete_loop()
    
    # Count active and deliveries
    active_fleet[y_idx] = sum(1 for b in buckets if b.active)
    if y_idx > 0:
        annual_he[y_idx] = sum(b.he_delivered for b in buckets 
                               if abs(b.return_year - year) < 0.5 and b.loops_completed > 0)
    cumulative_he[y_idx] = cumulative_he[y_idx-1] + annual_he[y_idx]

# %% [markdown]
# ## Fleet Performance

# %%
# Find steady-state years (after initial ~15 year ramp)
steady_mask = sim_years >= 2041
steady_he = annual_he[steady_mask]
avg_annual_he = np.mean(steady_he) if len(steady_he) > 0 else 0

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Active fleet
ax = axes[0, 0]
ax.plot(sim_years, active_fleet, 'b-', linewidth=2)
ax.axvline(2026 + 10, color='red', linestyle='--', alpha=0.5, label='First returns (~2036)')
ax.set_xlabel('Year')
ax.set_ylabel('Active Buckets')
ax.set_title('Active Fleet Size')
ax.grid(True, alpha=0.3)
ax.legend()

# Annual helium
ax = axes[0, 1]
ax.bar(sim_years, annual_he/1000, width=0.8, color='gold', edgecolor='orange', alpha=0.8)
ax.axhline(179, color='red', linestyle='--', alpha=0.7, label='Steady-state: 179 t/yr')
ax.set_xlabel('Year')
ax.set_ylabel('Helium Delivered (t/yr)')
ax.set_title('Annual Helium Delivery')
ax.grid(True, alpha=0.3)
ax.legend()

# Cumulative helium
ax = axes[1, 0]
ax.plot(sim_years, cumulative_he/1000, 'g-', linewidth=2)
ax.set_xlabel('Year')
ax.set_ylabel('Cumulative Helium (t)')
ax.set_title('Cumulative Helium Delivered')
ax.grid(True, alpha=0.3)

# Fleet composition
ax = axes[1, 1]
loops_data = [b.loops_completed for b in buckets if b.active]
ax.hist(loops_data, bins=range(0, max(loops_data)+2), color='steelblue', edgecolor='white')
ax.set_xlabel('Loops Completed')
ax.set_ylabel('Active Buckets')
ax.set_title('Fleet Experience Distribution')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Fleet Performance (steady state, {sim_years[steady_mask][0]:.0f}+):")
print(f"  Active Buckets:     {np.mean(active_fleet[steady_mask]):.0f}")
print(f"  Annual He:          {avg_annual_he/1000:.0f} t/yr")
print(f"  Cumulative He:      {cumulative_he[-1]/1000:.0f} t over {launch_years} years")
print("  Global He demand:   ~6,000 t/yr (pre-crisis)")
print(f"  Fleet share:        {avg_annual_he/1000/60:.1f}%")

# %% [markdown]
# ## Economics

# %%
bucket_cost    = 200e6  # $200M per Bucket
launch_cost_10 = 45e6   # $45M per Starship (10 Buckets)
total_cost_per_bucket = bucket_cost + launch_cost_10/buckets_per_launch
fleet_cost     = total_launched * total_cost_per_bucket
annual_replace = buckets_per_launch * total_cost_per_bucket  # per window
annual_revenue = avg_annual_he * 1500  # $1,500/kg He

print("Economics:")
print(f"  Cost per Bucket:      ${total_cost_per_bucket/1e6:.0f}M")
print(f"  Fleet startup:        ${fleet_cost/1e9:.1f}B")
print(f"  Annual replacement:   ${annual_replace/1e9:.2f}B/yr")
print(f"  Annual He revenue:    ${annual_revenue/1e6:.0f}M/yr  (@ $1,500/kg)")
print()
print("  The business case does NOT close on helium alone.")
print("  The hydrogen propellant delivered to orbit enables everything else.")
print("  Gas stations in space. The helium is the byproduct.")

# %% [markdown]
# ## Key Insight
# 
# Jupiter is the only planet where the magnetic scoop works:
# 
# - **42 km/s** orbital velocity → H2 self-ionizes (18 eV > 15.4 eV threshold)
# - **4.2 Gauss** magnetic field — strongest in the solar system
# - **10% helium** by number density (25% by mass)
# - **Deep gravity well** — Oberth effect amplifies departure burn
# 
# Saturn is too slow (26 km/s, won't ionize). Uranus and Neptune are too
# far and too cold. Jupiter is the only gas station in the solar system
# that opens itself when you show up at the right speed.

# %% [markdown]
# ---
# 
# ## What You've Learned
# 
# | Notebook | Concept | Python Module |
# |----------|---------|---------------|
# | 1. Your First Orbit | Vis-viva, circular orbits | `kepler.orbit` |
# | 2. Transfer to Mars | Hohmann, Lambert, porkchop | `kepler.transfer` |
# | 3. Atmospheric Entry | Drag, heating, corridor | `kepler.entry` |
# | 4. Design a Mission | ΔV budgets, ISRU, rocket eq | `kepler.orbit` |
# | 5. Build a Constellation | Fleet sim, attrition, economics | (numpy) |
# 
# You can now compute any orbit, find optimal transfer windows, analyze
# atmospheric entry, design complete missions, and simulate autonomous
# fleets over decades — all from Python.
# 
# *Kepler Workstation — MIT License. Real physics, real data, open source.*
