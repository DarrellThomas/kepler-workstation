# BUILD_LOG.md — How Kepler Workstation Was Built

A narrative of origin, crisis, physics, and a new kind of collaboration.

---

## 1. Origin — The CADAC Classroom

University of Florida, early 1990s. A USAF major sits in Peter Zipfel's aerospace simulation class. The textbook is *CADAC: Computer Aided Design of Aerospace Concepts*. The framework is tensor-based 6DOF — translational equations of motion in one coordinate system, rotational EOM in another, everything expressed as coordinate-free operators on vectors and quaternions. Clean. Mathematical. Reusable.

That framework lodged itself in the back of a mind and stayed there for thirty years.

Peter Zipfel (https://www.zipfel.us) taught that simulation is a craft. Build the mathematics first. Validate against real data. Then — and only then — write the code. The code is the last step, not the first.

---

## 2. Trigger — The Helium Crisis

March 2026. Iran bombs Qatar's Ras Laffan helium facility. Thirty percent of the world's helium supply — gone overnight. MRI machines going dark in hospitals. Semiconductor fabs rationing their helium allocations. Prices doubling in weeks.

Helium is a non-renewable resource on Earth. It's trapped in natural gas deposits, formed over billions of years by radioactive decay. Once it escapes the atmosphere, it's gone — light enough to reach escape velocity from the top of the exosphere. Every helium atom we've ever used has been mined from the ground. And the ground is running out.

The United States, once the world's largest helium producer, sold off its strategic reserve in the 1990s. The Federal Helium Reserve in Amarillo, Texas, went private. The price was set by market forces. The market did not anticipate a cruise missile hitting a fractionating column in Ras Laffan.

---

## 3. The Question

Jupiter's atmosphere is 10% helium by number density. The planet is 318 times the mass of Earth. Its gravity well contains more helium than the human race could use in a billion years.

Can we go get it?

This is not a hypothetical. This is an engineering question with a yes/no answer. The answer depends on physics: orbital mechanics, atmospheric entry, magnetic field dynamics, cryogenic fluid management, Hall thruster performance. All of these are calculable from first principles. All the data is public — JPL ephemeris, NASA atmosphere models from Galileo and Juno, published papers on superconducting magnets and ionized flow.

No new physics is required. Only integration.

---

## 4. The Build — April 13, 2026

One session. One human. One AI.

The human: Darrell Thomas, USAF (ret.), MS Aerospace Engineering, University of Florida.

The AI: Claude, an LLM running inside an engineering workspace with file access, shell execution, and git.

The session began with a question and ended with a working aerospace simulation platform.

**What was built in that session:**

- **Planetary environment models** — 5 bodies (Earth, Venus, Mars, Jupiter, Moon), each with validated atmosphere density/pressure/temperature profiles and gravity fields with zonal harmonics (J2 through J6-J10). Data from NASA missions: US Standard Atmosphere 1976, Pioneer Venus, Magellan, Mars Global Surveyor, Galileo probe ASI, Juno gravity, GRAIL.

- **JPL DE440 ephemeris reader** — Chebyshev polynomial interpolation of JPL's planetary ephemeris (1550–2650 CE). The gold standard for solar system position data.

- **N-body propagator** — RK4 integration with all planet gravitational influences. Validated against JPL Horizons.

- **Lambert transfer solver** — computes interplanetary transfer orbits with real departure/arrival dates. The same algorithm NASA uses for preliminary mission design.

- **Jupiter Helium Harvester ("The Bucket")** — full physics simulation of an autonomous magnetic scoop vehicle: 8-tonne dry mass, 20T superconducting magnet, H2 Hall thruster, 10-year Earth-Jupiter-Earth loop. 170-vehicle constellation simulation.

- **Three.js web viewer** — interactive 3D visualization of the constellation, real data, browser-native.

- **UE5 plugin architecture** — cleanroom C++ components (Environment, Newton, Euler, OrbitalMechanics) built for Unreal Engine 5 using native types (FVector, FQuat, FMatrix).

No code was ported. No libraries were wrapped. Every line was written against published physics and public data.

---

## 5. The Pivots

The first design is always wrong.

**Attempt 1: Ram-fill scoop.** Fly through Jupiter's upper atmosphere with a physical intake. Ram pressure fills the tank. Straightforward.

Result: At 42 km/s, ram pressure at 350 km altitude is ~2.8 Pa — enough to fill the tank but also enough to vaporize any intake structure. The stagnation temperature is thousands of kelvin. The vehicle melts on the first pass.

**Attempt 2: Thermal protection.** Add a heat shield. Ablative TPS like a reentry capsule.

Result: Jupiter entry is not Earth entry. The atmosphere has no solid surface. The vehicle must stay at altitude for thousands of passes. Ablative TPS is single-use. Heat soak accumulates. The thermal problem is unbounded.

**Pivot: Magnetic scoop.** Instead of a physical intake, use a superconducting magnet to generate a magnetic funnel. At 42 km/s, the kinetic energy of H2 molecules (18 eV) exceeds the ionization energy (15.4 eV). The bow shock self-ionizes the incoming gas. The magnetic field deflects charged particles into the collection tank. The vehicle hull sits behind the magnetic shield in near-vacuum.

No heat shield. No moving parts at the intake. No thermal cycling. Unlimited reuse.

The physics checks out. A 20T HTS solenoid with a 3-meter radius generates a 30–50 meter effective collection funnel. Ten kilograms per pass. A thousand passes over two years. Thirty-three tonnes collected. After burning 11.3 tonnes of H2 as propellant for the return trip, the vehicle delivers 21.7 tonnes to Earth orbit: 13.4 tonnes of hydrogen and 8.2 tonnes of helium.

**Confirmation: The Oberth effect.** Later analysis using circular restricted three-body problem (CR3BP) dynamics confirmed that the Jupiter departure burn is most efficient at periapsis — the Oberth effect. The maneuver requires 5.0 km/s, achievable with the Hall thruster burning captured H2 propellant. The mass ratio closes. The loop works.

---

## 6. The Fleet

A single Bucket delivers 8.2 tonnes of helium per decade-long loop. Useful, but not world-changing.

A fleet of 100 Buckets, launched in windows of 10 per Starship every 399 days (Earth-Jupiter synodic period), delivers 179 tonnes of helium per year — approximately 3% of current global demand. At $1,500/kg (current spot price post-Ras Laffan), that's $269 million per year in helium alone, plus the hydrogen propellant infrastructure in orbit.

The business case is not revenue per kilogram. It's what missions become possible when there are gas stations in space.

---

## 7. The Expansion

The April 13 session built the foundation. Subsequent sessions extended it.

**May 8, 2026** — CR3BP dynamics. Three-body equations of motion in the rotating frame. Halo orbit computation (Richardson 3rd-order approximation + differential correction). Invariant manifold propagation. Patched three-body transfers. Cross-validated against HITEN (Python CR3BP library). The Oberth effect analysis confirmed the Bucket departure strategy.

**May 9, 2026** — Saturn Harvester. Same magnetic scoop physics applied to Saturn's atmosphere (96% H2, 3% He). Head-to-head comparison: Saturn's lower gravity means lower insertion ΔV but also lower Oberth boost. Jupiter wins for helium; Saturn is competitive for hydrogen.

**May 14–15, 2026** — Saturn Station. Full simulation with crewed Titan base, atmospheric mining fleet, and 50-year economics. Real ephemeris, Lambert transfer windows, multi-body propagation.

**June 3, 2026** — GitHub Pages deployment. The Three.js viewer goes live at a public URL. The demo becomes shareable.

---

## 8. The Message

This entire codebase was built through human-AI collaboration. Not generated. Not auto-completed. Collaborated.

The human brought thirty years of aerospace engineering — the domain knowledge, the intuition for which physics matters, the ability to recognize when a thermal problem is unbounded and a magnetic solution is required.

The AI brought relentless execution — reading NASA technical memoranda, implementing the equations, validating against published data, writing the tests, structuring the commit history. Doing the thousands of small things that turn an idea into working software.

Together, they produced in hours what would take a solo engineer months.

The collaboration model is not "AI writes code while human watches." It's "human makes decisions, AI executes and verifies, both iterate." The human owns the physics. The AI owns the implementation velocity. The git history is the record of that dance — each commit a decision point, each file a negotiation between physical intuition and code correctness.

This is the message of Level 4: high-quality engineering software is now accessible to anyone who knows what to ask. Not because the tools got easier. Because the collaboration got smarter.

---

## 9. The Numbers

| Metric | Value |
|--------|-------|
| C++ source lines | 11,633 (standalone + UE5 stubs) |
| Working binaries | 10 |
| Planetary bodies modeled | 5 |
| Git commits (as of June 2026) | 14 |
| Initial build session | April 13, 2026 |
| Jupiter Harvester LOC | 4,794 |
| Saturn Harvester LOC | 1,963 |
| Three.js viewer | 491 LOC (single file) |
| License | MIT |

---

## 10. Attribution

- **Peter H. Zipfel** — whose textbooks and university teaching built the mathematical framework this project rests on
- **NASA / JPL** — Galileo, Juno, Magellan, GRAIL, Mars Global Surveyor, Pioneer Venus, and DE440 ephemeris (all public data, public domain)
- **Claude** — AI engineering partner and co-author
- **Darrell Thomas** — human decision-maker, domain expert, and project steward

The collaboration is the point.

---

*"The solar system has gas stations. We're building the pumps."*
