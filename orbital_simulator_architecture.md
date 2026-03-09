# Patch‑Conic & N‑Body Orbital Simulator

### Architecture, Physics Models, and Development Plan

Author: Student Project\
Language: Python\
Purpose: Class project orbital simulator with both patched‑conic
trajectory modeling and full N‑body propagation, plus a spatial grid for
visualization of gravitational potential and other environmental fields.

------------------------------------------------------------------------

# 1. Project Overview

This simulator models spacecraft trajectories using two different
approaches:

1.  **Patched‑Conic Approximation**
2.  **Full N‑Body Numerical Simulation**

The project also includes a **3D spatial grid system** for visualization
of fields such as:

-   gravitational potential
-   gravitational acceleration
-   solar radiation flux
-   other environmental quantities

The grid is primarily for visualization and analysis rather than physics
propagation.

------------------------------------------------------------------------

# 2. Physics Background

## 2.1 Newtonian Gravity

Gravitational acceleration between two bodies:

a = −GM r / \|r\|³

Where:

-   G = gravitational constant
-   M = mass of attracting body
-   r = vector from body to spacecraft

------------------------------------------------------------------------

## 2.2 N‑Body Gravity

For multiple bodies:

a = Σ −GMᵢ (r − rᵢ) / \|r − rᵢ\|³

This allows simulation of full gravitational interactions.

------------------------------------------------------------------------

## 2.3 Gravitational Potential

Scalar potential:

Φ = −GM / r

Multiple bodies:

Φ = Σ −GMᵢ / rᵢ

Acceleration is derived from the gradient:

a = −∇Φ

------------------------------------------------------------------------

# 3. Patched‑Conic Method

Patched conics approximate multi‑body motion by switching between
two‑body systems.

Regions:

Planet SOI → planet‑spacecraft two‑body problem\
Outside SOI → sun‑spacecraft two‑body problem

------------------------------------------------------------------------

## 3.1 Sphere of Influence

SOI radius:

r_SOI = a (m/M)\^(2/5)

Where:

-   a = planet orbital semimajor axis
-   m = planet mass
-   M = solar mass

------------------------------------------------------------------------

## 3.2 Patch‑Conic Algorithm

1.  Determine current sphere of influence
2.  Propagate Keplerian orbit around central body
3.  Detect SOI crossing
4.  Switch reference frame to new central body

------------------------------------------------------------------------

# 4. Gravity Assist Physics

Gravity assists are modeled as hyperbolic flybys.

Turning angle:

δ = 2 asin(1 / (1 + rp v∞² / μ))

Where:

-   rp = periapsis distance
-   v∞ = hyperbolic excess velocity
-   μ = gravitational parameter

------------------------------------------------------------------------

# 5. Simulation Architecture

The simulator should follow a layered architecture.

Visualization Layer\
Simulation Engine\
Physics Models\
Data Objects

This structure mirrors professional astrodynamics software.

------------------------------------------------------------------------

# 6. Core Data Objects

## 6.1 Celestial Body

Represents planets, moons, or stars.

Fields:

-   name
-   mass
-   radius
-   position
-   velocity
-   gravitational parameter μ
-   optional J2 coefficient
-   optional luminosity

Example structure:

class Body: def **init**(self,name,mass,radius,position,velocity):
self.name=name self.mass=mass self.radius=radius self.position=position
self.velocity=velocity self.mu=G\*mass

------------------------------------------------------------------------

## 6.2 Spacecraft

Fields:

-   mass
-   position
-   velocity
-   burn schedule

------------------------------------------------------------------------

## 6.3 Simulation State

Tracks the entire system state.

class State: def **init**(self): self.time=0 self.bodies=\[\]
self.spacecraft=\[\]

------------------------------------------------------------------------

# 7. Physics Modules

Suggested modules:

physics/ gravity.py potential.py radiation.py

These compute accelerations and environmental quantities.

------------------------------------------------------------------------

# 8. Propagation Engines

propagators/

-   kepler.py
-   patched_conic.py
-   nbody.py

------------------------------------------------------------------------

## 8.1 RK4 Integrator

Runge‑Kutta 4th order integration.

Steps:

k1 = f(x)\
k2 = f(x + dt/2)\
k3 = f(x + dt/2)\
k4 = f(x + dt)

x_next = x + dt/6\*(k1 + 2k2 + 2k3 + k4)

------------------------------------------------------------------------

# 9. Spatial Grid for Visualization

A 3D lattice representing space.

Nx × Ny × Nz cells

Example:

100 × 100 × 100

Grid spacing example:

dx = 10000 km

Coordinate mapping:

x = i \* dx\
y = j \* dy\
z = k \* dz

------------------------------------------------------------------------

## 9.1 Grid Data Fields

Possible stored fields:

-   gravitational potential
-   solar luminosity
-   radiation pressure
-   density fields

------------------------------------------------------------------------

## 9.2 Grid Class Example

class SpaceGrid: def **init**(self,nx,ny,nz,spacing): self.nx=nx
self.ny=ny self.nz=nz self.spacing=spacing
self.potential=np.zeros((nx,ny,nz))

------------------------------------------------------------------------

# 10. Computing Potential on the Grid

Algorithm:

For each cell:

1.  compute position
2.  compute distance to each body
3.  sum potentials

Φ_total = Σ −GM/r

------------------------------------------------------------------------

# 11. Simulation Engine

Main loop:

while simulation_running:

    propagate_planets(dt)
    propagate_spacecraft(dt)
    update_grid_fields()
    detect_events()
    render()

------------------------------------------------------------------------

# 12. Event Detection

Events include:

-   SOI crossing
-   periapsis passage
-   scheduled burns
-   gravity assists

Example:

if distance \< body.SOI: switch_reference_frame()

------------------------------------------------------------------------

# 13. Visualization

Recommended Python tools:

-   matplotlib
-   plotly
-   pyvista
-   vispy

Possible visualizations:

-   trajectory paths
-   gravity wells
-   potential field maps
-   spheres of influence

------------------------------------------------------------------------

# 14. Suggested Project Directory

orbital_sim/

bodies/ body.py spacecraft.py

physics/ gravity.py potential.py

propagators/ kepler.py nbody.py patched_conic.py

environment/ grid.py fields.py

mission/ burns.py flyby.py

visualization/ plot_orbits.py grid_renderer.py

simulation.py

------------------------------------------------------------------------

# 15. Recommended Development Order

Phase 1 Two‑body orbit propagation

Phase 2 Burn modeling

Phase 3 Patch‑conic switching

Phase 4 N‑body propagation

Phase 5 Spatial grid generation

Phase 6 Gravitational potential maps

Phase 7 Visualization tools

------------------------------------------------------------------------

# 16. Performance Considerations

Potential bottlenecks:

-   grid resolution
-   N‑body force calculations

Possible solutions:

-   numpy vectorization
-   numba acceleration
-   reduced grid update frequency

------------------------------------------------------------------------

# 17. Future Advanced Features

Possible improvements:

-   J2 planetary oblateness
-   solar radiation pressure
-   atmospheric drag
-   SPICE ephemeris data
-   adaptive time‑step integrators
-   GPU acceleration
-   spherical harmonic gravity models

------------------------------------------------------------------------

# 18. References

Recommended textbooks:

Fundamentals of Astrodynamics --- Bate, Mueller, White

Orbital Mechanics for Engineering Students --- Curtis

Fundamentals of Astrodynamics and Applications --- Vallado

Recommended software for inspiration:

GMAT\
Orekit\
Poliastro\
Open Space Toolkit
