# Kepler Workstation — Master Plan
## Working directory: /data/src/kepler-workstation/

## Context

We've built the physics, the simulations, and the cleanroom framework in a single session.
Now we need to turn it into something the world can see, use, learn from, and be inspired by.

### The KSP Opportunity
Kerbal Space Program proved that real orbital mechanics + accessible UI = millions of users.
KSP2 is dead (studio shut down May 2024, IP sold twice, now owned by "Fictions" who have
no announced plans). The community — millions of players who intuitively understand orbital
mechanics — has no home. Many of them are now aerospace engineering students and professionals.

Kepler Workstation is NOT a game. But it fills the void KSP left: a real-physics sandbox
where anyone can design missions, fly trajectories, and explore the solar system.
The difference: we use REAL data (JPL ephemeris, NASA atmosphere models, Juno gravity),
REAL engineering (6DOF, multi-body propagation), and REAL visualization (UE5).

KSP players are our Level 2 users. They already get it. They'll find us.

The project operates on four nested levels simultaneously:

**Level 1 — "Neat demo"**
The Jupiter Helium Harvester constellation. 170 autonomous vehicles streaming between
Earth and Jupiter. Real ephemeris. Real physics. A 15-second animation that stops
Elon Musk's thumb mid-scroll.

**Level 2 — "I can make neat demos with this"**
Kepler Workstation as a tool. Anyone downloads it, plugs in their own mission,
and gets the same quality simulation + visualization. The workstation is the product.
The Jupiter demo is just one application.

**Level 3 — "I could build something like this"**
Open source, MIT licensed, exhaustively documented. A high school student reads the
code and understands orbital mechanics. A university professor assigns it as coursework.
An indie developer forks it and builds their own space game on real physics.
The code teaches. The architecture inspires.

**Level 4 — "The entire landscape has changed"**
This entire codebase — planetary models, N-body propagator, Lambert solver, magnetic
scoop physics, UE5 framework, constellation fleet sim, business analysis — was built
in a single human-AI session. The README says so. The git history proves it.
The message: high-quality engineering software is now accessible to everyone.
Not because the tools got easier. Because the collaboration got smarter.

---

## What We Have (built today, April 13 2026)

### Repo: /data/src/kepler-workstation/ (git initialized, GPG signed)
- `Engine/KeplerWorkstation/` — UE5 cleanroom simulation components (6 files, compiles)
- `PlanetaryModels/` — 5-body atmosphere + gravity models (validated)
- `Ephemeris/` — JPL DE440 reader, N-body propagator, Lambert solver (validated)
- `Applications/JupiterHarvester/` — Full constellation sim (170 buckets, 36 years)
- `Viewer/` — constellation_data.json (exported, ready for Three.js)
- `LICENSE` — MIT, Copyright 2026 Darrell Thomas

### Working space: /data/src/zipfel/ (private, not published)
- CADAC++ v4 (Peter Zipfel's code, paid license — reference only)
- UE5 CadacVis plugin (compiled, headless-validated)
- Working notes, analysis, business cases

---

## Deliverables (prioritized)

### Milestone 1: The Demo (Level 1) — ship ASAP
**Goal:** One X post that Elon can't ignore.

1. **Three.js web viewer** (`Viewer/index.html`)
   - Loads constellation_data.json
   - Solar system: Sun, Venus, Earth, Mars, Jupiter orbits
   - Bucket swarm: colored dots streaming between planets
   - Time scrubber + play/pause + speed control
   - Stats overlay: year, buckets active, He delivered counter
   - Jupiter zoom: see scoop orbits
   - Must look stunning at 15 seconds of scrolling
   - Single HTML file, no build step, CDN Three.js

2. **Screen-record 15-second clip** for X post embed
   - Camera path: Earth → pull back → swarm in transit → Jupiter → counter ticking
   - MP4, no audio, optimized for X mobile autoplay

3. **Deploy to public URL** (GitHub Pages from the repo)
   - Clean URL: `[username].github.io/kepler-workstation`
   - Loads instantly, works on mobile

4. **Craft the X post**
   - Hook: helium crisis + physics insight
   - Visual: 15-second video
   - Link: interactive viewer
   - Tag: @elonmusk
   - Anchor: "Starlink for helium"

### Milestone 2: The Workstation (Level 2)
**Goal:** Other people can use Kepler Workstation for their own missions.

5. **Complete UE5 plugin** (`Engine/KeplerWorkstation/`)
   - Fill in stub implementations (Environment, OrbitalMechanics components)
   - Port planetary_environment.hpp into EnvironmentComponent properly
   - Port jpl_ephemeris/nbody/lambert into OrbitalMechanicsComponent
   - UE5 project template with example scene
   - Build.cs, .uplugin, properly structured for distribution

6. **Example applications** (`Applications/`)
   - JupiterHarvester: already done (constellation sim)
   - OrbitalTransfer: Earth-Mars Hohmann demo, porkchop plot generator
   - HypersonicEntry: Earth/Mars/Jupiter entry corridor analysis
   - Each with its own README, input files, expected output

7. **GPU compute** (`Engine/KeplerWorkstation/Source/.../GPU/`)
   - CUDA kernels: batch Lambert solver, Monte Carlo propagation
   - UE5 compute shader fallback for non-CUDA
   - 1000× speedup on constellation fleet sim

### Milestone 3: The Education (Level 3)
**Goal:** The code teaches aerospace engineering.

8. **Documentation as curriculum**
   - Each component header: WHAT it does, WHY it works, WHERE the math comes from
   - `docs/` folder: tutorial progression
     - Tutorial 1: "Your first orbit" (circular orbit, vis-viva, period)
     - Tutorial 2: "Transfer to Mars" (Hohmann, Lambert, porkchop)
     - Tutorial 3: "Atmospheric entry" (drag, heating, TPS)
     - Tutorial 4: "Design a mission" (full VEEJ trajectory)
     - Tutorial 5: "Build a constellation" (fleet, attrition, economics)
   - Each tutorial: theory → code → run it → see the result in 3D

9. **Interactive notebooks** (Jupyter or Observable)
   - Python bindings for the C++ models
   - Plot orbits, atmospheres, gravity fields inline
   - Students modify parameters and see results immediately

10. **Standalone mode** (no UE5 required)
    - PlanetaryModels + Ephemeris compile standalone with g++
    - Terminal-based demos that output CSV
    - Lower barrier to entry than full UE5 install

### Milestone 4: The Message (Level 4)
**Goal:** Demonstrate what AI + human collaboration produces.

11. **The git history tells the story**
    - Atomic commits showing the progression
    - From "what code was ported from Storage_3" to "constellation of 170 vehicles"
    - In one session. With one human and one AI.

12. **BUILD_LOG.md** — honest narrative
    - Origin: Peter Zipfel's CADAC class, University of Florida, early 1990s.
      First exposure to tensor-based 6DOF aerospace simulation.
      (link: https://www.zipfel.us)
    - Trigger: March 2026. Iran bombs Qatar's Ras Laffan helium facility.
      30% of global helium supply gone overnight. MRI machines going dark.
      Chip fabs rationing. Helium prices doubled in weeks.
    - The question: Jupiter's atmosphere is 10% helium. Can we go get it?
    - The build: One session. Human + AI. Real JPL ephemeris, real atmospheric
      models from Galileo/Juno/Magellan probes, real orbital mechanics.
    - The pivots: ram-fill scoop → thermal kills it at every altitude →
      magnetic scoop (bow shock self-ionizes at 42 km/s) → it works
    - The answer: yes. 8 tonnes, 3 valves, 1 magnet. Fleet of 100.
      179 tonnes of helium per year, forever.

13. **Attribution**
    - Peter H. Zipfel (https://www.zipfel.us): whose textbooks and university
      teaching inspired the mathematical framework
    - NASA/JPL: mission data (Galileo, Juno, GRAIL, Magellan, MGS)
    - Claude: AI engineering partner (co-author on commits)
    - The collaboration model is the point

---

## Architecture Decisions

### Why UE5 (not Godot, Unity, or custom renderer)
- Full source build available (we have it at /data/src/UE/)
- C++ native (matches simulation code, no language boundary)
- Large World Coordinates (aerospace scales without precision loss)
- Proven plugin architecture
- GPU compute pipeline (Niagara, compute shaders)
- Industry credibility

### Why cleanroom (not fork CADAC++)
- Licensing: Peter's code is paid AIAA product
- Simplicity: modern C++ is cleaner than 2003-era code
- UE5 native: use FVector/FQuat/FMatrix, not custom Matrix class
- Teachability: students learn UE5 AND physics simultaneously
- IP clarity: MIT license, unambiguous ownership

### Why standalone + UE5 (not UE5-only)
- PlanetaryModels and Ephemeris work without UE5
- Lower barrier: g++ is free, UE5 is a large download
- CI/CD: standalone tests run in any pipeline
- Portability: models reusable in Python, Rust, web (via WASM)

### File structure for the repo
```
kepler-workstation/
├── LICENSE (MIT, Darrell Thomas)
├── README.md
├── BUILD_LOG.md                      ← Level 4: the human story
├── .gitignore (Private/, MissileSim/, *.440)
│
├── Engine/KeplerWorkstation/         ← Level 2: the workstation
│   ├── KeplerWorkstation.uplugin
│   └── Source/KeplerWorkstation/
│       ├── Public/Sim/               ← headers (public API)
│       └── Private/Sim/              ← implementations
│
├── PlanetaryModels/                  ← Level 3: standalone, educational
│   ├── planetary_environment.hpp
│   └── test_models.cpp
│
├── Ephemeris/                        ← Level 3: standalone, educational
│   ├── jpl_ephemeris.hpp
│   ├── nbody_propagator.hpp
│   ├── lambert.hpp
│   └── test_propagator.cpp
│
├── Applications/                     ← Level 1: the demos
│   ├── JupiterHarvester/
│   │   ├── README.md (THE BUCKET)
│   │   ├── bucket_sim.cpp
│   │   ├── constellation_sim.cpp
│   │   ├── magnetic_scoop_sim.cpp
│   │   └── export_constellation.cpp
│   ├── OrbitalTransfer/              ← future
│   └── HypersonicEntry/              ← future
│
├── Viewer/                           ← Level 1: the demo visual
│   ├── index.html (Three.js)
│   └── constellation_data.json
│
└── docs/                             ← Level 3: tutorials
    ├── tutorial-01-first-orbit.md
    ├── tutorial-02-transfer-to-mars.md
    └── ...
```

---

## Immediate Next Actions (priority order)

1. **Build Three.js viewer** (Viewer/index.html) — this is the demo
2. **Deploy to GitHub Pages** — make it clickable
3. **Screen-record 15s video** — the X post visual
4. **Craft the X post** — the hook
5. **Fill UE5 component stubs** — make the workstation actually work
6. **Write Tutorial 1** — "Your first orbit" — the educational entry point
7. **Write BUILD_LOG.md** — the Level 4 narrative

---

## Verification

- Three.js viewer loads constellation_data.json and renders animated solar system
- Standalone models compile and validate: `cd PlanetaryModels && g++ -std=c++11 -o test test_models.cpp && ./test`
- Standalone ephemeris works: `cd Ephemeris && g++ -std=c++11 -O2 -o test test_propagator.cpp && ./test`
- UE5 plugin compiles: `Build.sh KeplerWorkstationEditor Linux Development`
- GitHub Pages serves the viewer at public URL
- The entire repo can be cloned and used by someone who has never seen it before
