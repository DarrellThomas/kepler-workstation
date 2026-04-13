# Kepler Workstation — Three.js Constellation Viewer

## What it does

Interactive 3D visualization of the Jupiter Helium Harvester constellation.
170 autonomous vehicles streaming between Earth and Jupiter over 26 years,
using real JPL DE440 ephemeris data and validated orbital mechanics.

Single HTML file. No build step. Runs in any modern browser.

## How to run

```bash
cd Viewer
python3 -m http.server 8081
# Open http://localhost:8081
```

## Current features

- **Solar system**: Sun, Venus, Earth, Mars, Jupiter with orbital paths
- **Jupiter**: Real NASA texture, slow rotation
- **Harvester fleet**: 170 dots colored by state (red=outbound, yellow=filling, green=full)
- **Camera**: Starts behind Jupiter, tracks selected object, user can break free with drag
- **Click interaction**: Click any harvester or planet to lock camera + show telemetry
- **Engineering dashboard**: Green-on-black telemetry panel with real physics data
  - Navigation: state, altitude, velocity
  - Atmosphere: density, temperature
  - Magnetic scoop: field status, funnel radius, flow rate, hull temperature
  - Tanks: H2, He, fill percentage
  - Propulsion: thruster status, remaining delta-V
  - Planet data when clicking planets
- **Time controls**: Play/pause, speed (1/8× to 10×), timeline scrubber
- **HUD**: Year, launched, active, delivered, lost, cumulative He tonnes

## Data source

`constellation_data.json` — exported from `export_constellation.cpp`
Contains planet positions (10-day intervals) and harvester waypoints
computed from real DE440 ephemeris + Lambert transfer solutions.

## Known issues / TODO

### Visual
- [ ] Harvester dots near Jupiter cluster together — need slight orbital spread
- [ ] Jupiter texture loads from Wikipedia CDN — should bundle locally or use NASA source
- [ ] No trail lines showing harvester paths (removed for clutter, could be optional toggle)
- [ ] No visual distinction between individual harvesters beyond color
- [ ] Planet labels missing (removed for clutter — add as optional toggle)
- [ ] Sun glow is basic — could use bloom post-processing
- [ ] No grid or scale reference — hard to judge distances

### Interaction
- [ ] Double-click to reset camera to default (behind Jupiter)
- [ ] Keyboard shortcuts: space=play/pause, +/-=speed, tab=cycle through harvesters
- [ ] Hover tooltip showing harvester ID before clicking
- [ ] Click empty space should deselect cleanly (sometimes sticky)
- [ ] Pinch-to-zoom on mobile not tested

### Dashboard
- [ ] Dashboard data is simulated from phase — should interpolate from actual sim output
- [ ] Random jitter on hull temperature (placeholder) — should be deterministic
- [ ] No graph/chart for historical data (tank fill over time, etc.)
- [ ] No mission timeline showing upcoming events for selected harvester

### Data
- [ ] Only 170 harvesters (17 launch windows × 10) — expand to full 330 for 36-year sim
- [ ] Waypoint interpolation is linear — should be curved (conic arcs)
- [ ] Scooping harvesters all sit on top of Jupiter — need orbital offset per pass
- [ ] Return trajectories are straight-line interpolated — should follow Lambert solutions
- [ ] No Ganymede gravity assist modeled in trajectory export

### Performance
- [ ] 170 spheres is fine — may need instanced meshes at 1000+ harvesters
- [ ] JSON is 421 KB — fine for now, would need binary format at scale

### Demo video
- [ ] Screen-record 15-second clip for X post
- [ ] Camera path: start behind Jupiter → harvesters arrive → scooping → return stream
- [ ] Kling AI prompt ready for close-up collection sequence (see kling_prompt.md)
- [ ] Dashboard overlay to composite on Kling video output
- [ ] Need to sync dashboard numbers with video timing

### Deployment
- [ ] Push to GitHub repo
- [ ] Enable GitHub Pages on Viewer/ directory
- [ ] Clean URL for sharing
- [ ] Test on mobile browsers
- [ ] Add Open Graph meta tags for X post embed preview

## Architecture

```
Viewer/
├── index.html                 ← Single-file Three.js app (everything inline)
├── constellation_data.json    ← 421 KB trajectory data from sim
├── kling_prompt.md            ← AI video generation prompt
└── README.md                  ← This file
```

No dependencies except Three.js loaded from CDN.
No build step. No node_modules. No webpack. Just open and go.
