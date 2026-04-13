# Kepler Workstation — Claude Code Instructions

## Project

Open-source aerospace simulation and mission design platform built on Unreal Engine 5.
Copyright (c) 2026 Darrell Thomas. MIT License.

Working directory: `/data/src/kepler-workstation/`
Private working space (not published): `/data/src/zipfel/`

## Architecture

- `Engine/KeplerWorkstation/` — UE5 plugin with simulation components
- `PlanetaryModels/` — Standalone atmosphere + gravity models (5 bodies)
- `Ephemeris/` — JPL DE440 reader, N-body propagator, Lambert solver
- `Applications/` — Mission demos (JupiterHarvester, OrbitalTransfer, etc.)
- `Viewer/` — Three.js web viewer for demos
- `docs/` — Tutorials and documentation

## Conventions

- All simulation code: Copyright Darrell Thomas, MIT License
- Cleanroom implementation — zero third-party simulation code
- Physics from published textbooks and NASA public data
- UE5 components use FVector/FQuat/FMatrix (not custom math classes)
- Standalone code compiles with `g++ -std=c++11`
- Commit messages: descriptive, co-authored with Claude

## Key Context

- Owner is a USAF Major (ret), MS Aerospace from University of Florida
- Took Peter Zipfel's CADAC simulation class at UF in early 1990s
- Peter's work inspires the mathematical framework (tensors, 6DOF) but no code is used
- The Jupiter Harvester ("The Bucket") is the flagship demo application
- Target audience: high school students through NASA engineers
- KSP community is a key audience (millions of orbital mechanics enthusiasts with no game)

## Current Plan

Read `.claude/plans/master-plan.md` for the full roadmap.

Priority: Three.js viewer → GitHub Pages → X post → UE5 stubs → tutorials

## .gitignore Rules

- `Private/` — working notes, not published
- `Applications/MissileSim/` — export control sensitivity
- `Ephemeris/*.440` — large binary (users download separately)

## What NOT to do

- Do not include any Peter Zipfel source code in this repo
- Do not reference private working notes or AIAA purchase history
- Do not use CADAC++ variable names, data structures, or code patterns
- Do not commit large binary files (ephemeris data, UE5 binaries)
