# Titan Survival — Technical Scope

## Phase 0: Proof of Concept (Vertical Slice)

**Goal**: One player, one ice cave, walk around on Titan. Prove it feels right.

### Deliverables
- [ ] Terrain chunk from GTDR/SAR data (10km × 10km area near Huygens landing site)
- [ ] Titan atmosphere shader (orange haze, twilight lighting, distance fog)
- [ ] Player controller at 0.14g (bouncy movement, low terminal velocity)
- [ ] Basic ice cave: carved interior room, thermal door, outside/inside transition
- [ ] Temperature system: cold ticks down outside, warm inside
- [ ] O2 system: meter drains outside, replenishes inside
- [ ] Day/night cycle (16 Earth days, orange twilight → pitch black)
- [ ] Saturn in the sky (stationary, visible through haze)
- [ ] Ambient sound: wind at 1 m/s, dense atmosphere acoustics (lower pitch, speed of sound 197 m/s)

### Technical Decisions Needed
- **Engine**: UE5 vs Godot 4
  - UE5: better lighting (Lumen), better large terrain (Nanite), harder to mod
  - Godot 4: open source, faster iteration, aligns with MIT license, needs custom shaders
- **Terrain pipeline**: GTDR heightmap → game engine terrain
  - SAR mosaic as surface texture (dark tholins, ice outcrops)
  - ISS mosaic for distant terrain color
- **Scale**: 1:1 Titan (16,177 km circumference) vs scaled region
  - Start with 10×10 km tile. Expand later with streaming.

### Estimated Effort
- Terrain import + shader: 2 weeks
- Player controller + physics: 1 week
- Ice cave prototype: 2 weeks
- Survival meters (O2, heat): 1 week
- Atmosphere visuals + audio: 2 weeks
- **Total Phase 0: ~8 weeks**

---

## Phase 1: Core Survival Loop

**Goal**: Playable survival from landing to stable base.

### Deliverables
- [ ] Drop pod landing sequence (entry from orbit, atmosphere braking)
- [ ] RTG starter item (initial power source)
- [ ] Ice mining tool (carve the walls)
- [ ] Electrolyzer (ice → O2 + H2)
- [ ] Crafting system: thermal suit, O2 bottles, insulation panels, thermal doors
- [ ] Base building: carve ice, place insulation, install O2 vents, seal rooms
- [ ] Heat system: RTG radiates warmth, CH4 burner for larger areas
- [ ] Food: basic algae farm with LED grow lights
- [ ] Base decay: unheated rooms slowly freeze over (ice encroachment)
- [ ] Death/respawn mechanics

### Technical Challenges
- **Voxel or mesh carving** for ice caves: need destructible/constructible terrain
  - UE5: Voxel Plugin or custom SDF approach
  - Godot: VoxelTools addon or custom
- **Thermal simulation**: heat propagation through ice walls, insulation layers
  - Simplified: per-room temperature, insulation rating, heat sources
  - Advanced: 3D heat diffusion grid (GPU compute)
- **Gas simulation**: O2 concentration per room, leaks through unsealed doors
  - Simplified: per-room O2 level, door open = drain rate

### Estimated Effort
- Ice cave building system: 4 weeks
- Crafting + items: 2 weeks
- Survival systems (heat, O2, food): 3 weeks
- Base decay: 1 week
- Drop pod / landing: 2 weeks
- Polish + balancing: 2 weeks
- **Total Phase 1: ~14 weeks**

---

## Phase 2: Transportation & Exploration

**Goal**: Leave the base, explore Titan.

### Deliverables
- [ ] Pedal-powered aircraft (0.14g + 4.4× density aerodynamics)
- [ ] Flight physics from HASI atmosphere model + DWE winds
- [ ] Balloon system (hot N2, variable altitude, cargo hauling)
- [ ] Expanded terrain: multiple 10×10 km tiles, streaming
- [ ] Methane lakes (liquid surface, Kraken Mare region)
- [ ] Methane rain weather events
- [ ] Surface features: tholin dunes, ice pebble riverbeds, cryovolcano terrain
- [ ] Points of interest: Huygens landing site, Sotra Patera, Mithrim Montes
- [ ] Lake boat (methane surface vessel)
- [ ] Map/navigation system

### Technical Challenges
- **Flight model**: Real aerodynamics at Re ~2M, low speed, dense atmosphere
  - Lift/drag from HASI density profile
  - Wind from DWE data (altitude-dependent)
- **Terrain streaming**: Titan is 16,177 km around. Need LOD + streaming.
  - Heightmap from GTDR, texture from SAR/ISS
  - Only render ~50 km radius at detail
- **Liquid methane lakes**: Reflective surface, low waves, shore interaction
  - Mirror-smooth most of the time (< 3mm waves measured by Cassini)

### Estimated Effort
- Pedal aircraft + flight model: 4 weeks
- Balloon system: 2 weeks
- Terrain streaming: 4 weeks
- Methane lakes: 3 weeks
- Weather system: 2 weeks
- Points of interest: 2 weeks
- **Total Phase 2: ~17 weeks**

---

## Phase 3: Colony Systems

**Goal**: Build a real colony. Sabatier loop. Centrifuge. Agriculture.

### Deliverables
- [ ] Sabatier reactor (CO2 + H2 → CH4 + H2O closed loop)
- [ ] Methane pipeline system (lake → base)
- [ ] Nuclear reactor (1 MW, late-game power upgrade)
- [ ] Centrifuge gym (50-100m radius, 1g exercise)
- [ ] Advanced agriculture: multiple crop types, insect farms, tholin processing
- [ ] Methane-powered aircraft (late game, long range)
- [ ] Upper atmosphere floating habitat (balloon city concept)
- [ ] Research system: study tholins, model chemistry, unlock tech tree

### Estimated Effort
- Colony systems: 6 weeks
- Tech tree / research: 3 weeks
- Advanced vehicles: 4 weeks
- Floating habitat: 3 weeks
- **Total Phase 3: ~16 weeks**

---

## Phase 4: Multiplayer

**Goal**: Private server, 20-50 players.

### Deliverables
- [ ] Dedicated server architecture
- [ ] Player authentication + server browser
- [ ] Shared world persistence (bases, terrain modifications)
- [ ] Territory/claim system
- [ ] Cooperative building
- [ ] Chat / voice (with 75-min delay to "Earth" chat channel for flavor)
- [ ] Admin tools (server config, player management)

### Estimated Effort
- Netcode + server: 8 weeks
- Persistence: 4 weeks
- Social systems: 3 weeks
- Admin tools: 2 weeks
- **Total Phase 4: ~17 weeks**

---

## Phase 5: Polish & Release

- [ ] Tutorial / onboarding (learn survival mechanics)
- [ ] Sound design (dense atmosphere acoustics, wind, machinery)
- [ ] Music (ambient, isolation, wonder)
- [ ] UI/UX polish
- [ ] Performance optimization
- [ ] Steam Early Access launch
- [ ] Modding support

---

## Total Estimated Scope

| Phase | Weeks | Description |
|-------|-------|-------------|
| 0 - Vertical Slice | 8 | Walk on Titan, ice cave, survive |
| 1 - Core Survival | 14 | Landing → stable base |
| 2 - Transportation | 17 | Flight, balloons, lakes, exploration |
| 3 - Colony Systems | 16 | Sabatier, nuclear, centrifuge, agriculture |
| 4 - Multiplayer | 17 | Private servers, 20-50 players |
| 5 - Polish | 8 | Audio, tutorial, optimization, launch |
| **Total** | **~80 weeks** | |

Phase 0 is the gate. If walking on Titan feels magical — the orange light, the bouncy gravity, the dense air pushing back on you, Saturn hanging in the haze — then the rest follows.

---

## Competitive Advantage

1. **Real data**: 2.9 GB of Cassini-Huygens data. No other game has this.
2. **Real physics**: HASI atmosphere, DWE winds, actual aerodynamics.
3. **Novel setting**: Nobody has done walkable Titan. First mover.
4. **KSP + Rust crossover audience**: Millions of potential players.
5. **Educational**: Teaches real planetary science through gameplay.
6. **Open source foundation**: Kepler Workstation (MIT) provides the simulation core.
