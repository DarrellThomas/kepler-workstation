# Engine Analysis: What should Kepler Workstation run on?

## The Question

If the goal is KSP-genre accessibility — "curious high school student to NASA engineer" —
is UE5 the right framework? Or is there something better?

## Requirements

1. **Real-time 3D visualization** of orbits, planets, vehicles at solar system scale
2. **Precision math** — double precision throughout (no float32 position jitter at large distances)
3. **Cross-platform** — Linux, Windows, Mac minimum. Web is a huge plus.
4. **Accessible** — low barrier to download, install, and run
5. **Moddable** — users should be able to add their own vehicles, missions, planets
6. **Open source compatible** — MIT license, no proprietary dependencies in core sim
7. **GPU compute** — Monte Carlo, batch trajectories, fleet simulation
8. **Good-looking** — needs to be visually compelling, not just functional

## Candidates

### UE5 (current choice)

| Pros | Cons |
|------|------|
| Already compiled and working | 200+ GB install (source build) |
| C++ native (matches sim code) | Proprietary (Epic license, not MIT) |
| Large World Coordinates (double precision) | Massive complexity for what we need |
| Industry-standard rendering | Overkill for orbit lines and dots |
| Plugin architecture proven | Steep learning curve for contributors |
| GPU compute via Niagara/compute shaders | Build times measured in hours |
| Already have source build at /data/src/UE/ | Can't ship UE5 as part of our MIT repo |

**Verdict:** Powerful but heavy. Good for the engineering workstation (Level 2).
Bad for accessibility (Level 3). Contributors need UE5 source access.

### Godot 4

| Pros | Cons |
|------|------|
| **MIT licensed** — matches our project | No built-in double precision (float32 world) |
| ~100 MB download (vs 200 GB for UE5) | Smaller ecosystem than UE5 |
| GDScript + C++ + C# support | Less mature rendering pipeline |
| Huge indie community (KSP audience!) | No Large World Coordinates natively |
| Runs on Linux/Win/Mac/Web (HTML5 export!) | GPU compute less mature |
| Easy to mod and extend | Would need custom double-precision layer |
| Active development, growing fast | |

**Verdict:** Perfect fit for accessibility. MIT licensed like us. Web export is killer.
Double precision needs a custom solution (store positions as doubles, render with
camera-relative float32 — same trick KSP used). The KSP community already knows Godot
because several KSP-inspired projects use it (Spaceflight Simulator, Juno: New Origins).

### Bevy (Rust)

| Pros | Cons |
|------|------|
| ECS architecture (perfect for simulation) | Very immature (pre-1.0) |
| Rust safety + performance | Small community |
| MIT/Apache dual licensed | No visual editor |
| WebGPU support (browser-native) | Steep Rust learning curve |
| Double precision possible natively | Limited asset pipeline |
| Hot reloading | Not ready for production |

**Verdict:** Architecturally ideal but too immature. Revisit in 2-3 years.

### Custom WebGPU (Three.js / raw WebGPU)

| Pros | Cons |
|------|------|
| **Runs in any browser** — zero install | Limited to browser capabilities |
| Three.js already working (we built the viewer) | No real "engine" — build everything |
| WebGPU for compute (GPU in browser) | No editor, no scene graph |
| Instant sharing — just a URL | Performance ceiling lower than native |
| JavaScript/TypeScript — most developers know it | Double precision in JS is native (float64!) |
| Perfect for Level 1 demo | Hard to scale to full workstation |

**Verdict:** Already proven for the demo viewer. Could be the primary platform for
the educational tier (Level 3). JS has native float64 — no precision issues!
WebGPU for compute is shipping in Chrome/Edge now.

### Hybrid Approach (recommended investigation)

```
Level 1 (Demo):      Three.js viewer (DONE — Viewer/index.html)
Level 2 (Tool):      Godot 4 with C++ GDExtension for sim engine
Level 3 (Education): Three.js + Observable/Jupyter notebooks
Level 4 (Message):   The code itself (any engine, MIT licensed)
```

OR:

```
All levels:          Godot 4
                     - MIT licensed (can ship as part of our repo)
                     - Web export (Level 1 demo runs in browser)
                     - Desktop app (Level 2 workstation)
                     - GDScript for accessibility (Level 3 education)
                     - C++ GDExtension for performance (sim engine)
```

## Key Insight

UE5 is the **best renderer** but the **worst choice for accessibility**.
The person we most want to reach — the high school student — is never going to
download 200 GB of UE5 source and compile it.

Godot + web export means that student clicks a link and is designing orbits in 5 seconds.

## Decision Needed

The sim engine (physics, propagation, atmospheric models) is engine-agnostic.
It's C++ headers. It can plug into UE5, Godot, Bevy, or compile standalone.

The question is: what RENDERS it?

Options:
A. **UE5 only** — powerful but exclusive
B. **Godot only** — accessible but less pretty
C. **Web only** — most accessible but limited
D. **Multi-target** — sim engine in C++, render backends for UE5 + Godot + Web
E. **Godot primary, UE5 optional** — Godot for distribution, UE5 for premium vis

Recommendation: **E** — Godot primary, UE5 optional. Or possibly **C** pushed to its
limits with WebGPU. The web viewer we already built is surprisingly capable.

## Action Item

Before deciding: try the Three.js viewer on desktop. If it looks good enough,
maybe web IS the answer and we don't need a game engine at all.
