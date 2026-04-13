# Kepler Workstation — Task List

## Carry-over from April 13, 2026 session

### Milestone 1: The Demo (ship ASAP)

- [x] Export constellation trajectory data to JSON
- [x] Build cleanroom UE5 simulation framework
- [x] Build Three.js web viewer (Viewer/index.html)
- [ ] Test Three.js viewer locally (need desktop — `python3 -m http.server` in Viewer/)
- [ ] Deploy viewer to GitHub Pages (need GitHub repo created first)
- [ ] Screen-record 15-second demo video for X post
- [ ] Craft the X post (text + video + link + @elonmusk)

### Milestone 2: The Workstation

- [ ] Integrate N-body propagator into UE5 OrbitalMechanicsComponent (currently stub)
- [ ] Integrate planetary models into UE5 EnvironmentComponent (currently stub)
- [ ] Integrate Lambert solver into UE5 OrbitalMechanicsComponent (currently stub)
- [ ] Integrate JPL DE440 reader into UE5 OrbitalMechanicsComponent (currently stub)
- [ ] Integrate magnetic scoop into UE5 MagneticScoopComponent (working but needs UE5 wiring)
- [ ] Create UE5 project template (KeplerWorkstation.uproject) with example scene
- [ ] Build OrbitalTransfer example application (Earth-Mars Hohmann, porkchop plot)
- [ ] Build HypersonicEntry example application
- [ ] Add GPU compute support (CUDA kernels for Monte Carlo and batch Lambert)

### Milestone 3: Education

- [ ] Write Tutorial 1: "Your first orbit" (vis-viva, circular orbit, period)
- [ ] Write Tutorial 2: "Transfer to Mars" (Hohmann, Lambert, porkchop)
- [ ] Write Tutorial 3: "Atmospheric entry" (drag, heating, TPS)
- [ ] Write Tutorial 4: "Design a mission" (full VEEJ trajectory)
- [ ] Write Tutorial 5: "Build a constellation" (fleet, attrition, economics)
- [ ] Add inline documentation to all component headers (WHAT/WHY/WHERE)

### Milestone 4: The Message

- [ ] Write BUILD_LOG.md (origin: Zipfel class → Qatar bombing → Jupiter → magnetic scoop)
- [ ] Structure git history into meaningful atomic commits for the public narrative

### Open Engineering Questions

- [ ] Validate magnetic scoop efficiency with MHD simulation (50% assumed, may be optimistic)
- [ ] 20T HTS coil at 3m radius — mass/power/quench analysis
- [ ] Magnetic leakage rate → hull heating (1% assumed, needs field topology analysis)
- [ ] Jupiter radiation belt dose on electronics per scoop pass
- [ ] H2 Hall thruster at Isp 3000s — lab data exists, not flight-proven
- [ ] Aerocapture at Earth at 14.7 km/s — TPS design for reusable vehicle
- [ ] LH2 boiloff over 3-year return cruise — thermal analysis needed

### Naming

- [ ] Decide vehicle name: Dipper? Prometheus? Shearwater? (owner's choice)
- [ ] Decide game engine: UE5 vs Godot vs other (see analysis below)

### Game Engine Decision (pending)

If we want KSP-genre accessibility (not just engineering workstation), evaluate:
- UE5: current choice, compiled, powerful but heavy
- Godot: open source, lighter, growing community, MIT licensed like us
- Bevy (Rust): ECS architecture, modern, but immature
- Custom WebGPU: runs in browser natively, no download needed
- See .claude/ENGINE_ANALYSIS.md when created
