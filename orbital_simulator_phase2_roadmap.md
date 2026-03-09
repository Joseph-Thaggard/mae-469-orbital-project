# Orbital Simulator -- Phase 2 Roadmap

### Advanced Mission Design & Trajectory Optimization

This document extends the Phase‑1 simulator architecture with advanced
mission‑design capabilities such as Lambert solvers, gravity assists,
trajectory optimization, and ephemeris integration.

Language: Python Purpose: Expand the simulator into a mission analysis
and trajectory design tool.

------------------------------------------------------------------------

# 1. Phase‑2 Objectives

Phase‑2 transforms the simulator from a trajectory propagator into a
**mission design system**.

New capabilities:

-   Lambert transfer solver
-   Porkchop launch window plots
-   Gravity assist targeting
-   Trajectory optimization
-   SPICE ephemeris integration
-   Lagrange point / CR3BP modeling
-   Advanced visualization tools

These features mirror those used in professional astrodynamics tools.

------------------------------------------------------------------------

# 2. High‑Level Architecture

The architecture now includes a new **Mission Design Layer**.

Visualization Layer\
Mission Design Layer\
Simulation Engine\
Physics Models\
Data Objects

------------------------------------------------------------------------

# 3. New Directory Structure

orbital_sim/

bodies/\
environment/\
physics/\
propagators/

mission/\
burns.py\
flyby.py

mission_design/\
lambert.py\
porkchop.py\
flyby_targeting.py

optimization/\
trajectory_optimizer.py

visualization/\
porkchop_plot.py\
trajectory_plots.py

simulation.py

------------------------------------------------------------------------

# 4. Lambert Transfer Solver

## Purpose

Solve for the orbit connecting two points in space in a given time.

Inputs:

-   r1: departure position
-   r2: arrival position
-   tof: time of flight
-   μ: gravitational parameter

Outputs:

-   v1: departure velocity
-   v2: arrival velocity

------------------------------------------------------------------------

## Module

mission_design/lambert.py

Example interface:

def lambert_solver(r1, r2, tof, mu): return v1, v2

------------------------------------------------------------------------

# 5. Porkchop Plot Generator

Porkchop plots visualize optimal launch windows.

Axes:

x‑axis → launch date\
y‑axis → arrival date

Color map → required Δv

------------------------------------------------------------------------

## Workflow

1.  Generate departure dates
2.  Generate arrival dates
3.  Solve Lambert for each pair
4.  Compute Δv
5.  Plot contour map

------------------------------------------------------------------------

## Module

mission_design/porkchop.py

Example:

def generate_porkchop(departure_range, arrival_range):

------------------------------------------------------------------------

# 6. Gravity Assist Targeting

Gravity assists allow spacecraft to change velocity using planetary
flybys.

Turning angle:

δ = 2 asin(1 / (1 + rp v∞² / μ))

Variables:

-   incoming velocity v∞
-   periapsis distance
-   planetary μ

Outputs:

-   outgoing velocity vector

------------------------------------------------------------------------

## Module

mission_design/flyby_targeting.py

Example:

def compute_flyby(v_inf_in, rp, mu):

------------------------------------------------------------------------

# 7. Trajectory Optimization

Optimization determines the best mission parameters.

Optimization variables may include:

-   launch time
-   burn magnitude
-   burn direction
-   flyby altitude
-   transfer time

Objectives:

-   minimize total Δv
-   minimize time of flight

------------------------------------------------------------------------

## Optimization Algorithms

Recommended algorithms:

### Differential Evolution

From SciPy:

scipy.optimize.differential_evolution

Advantages:

-   handles nonlinear problems
-   works with flyby trajectories
-   simple to implement

------------------------------------------------------------------------

## Cost Function Example

def trajectory_cost(parameters):

    trajectory = simulate(parameters)

    return total_delta_v

------------------------------------------------------------------------

# 8. Example Optimization Problem

Example mission:

Earth → Venus → Earth → Jupiter

Variables:

-   launch date
-   transfer time
-   flyby altitude
-   arrival time

Objective:

minimize Δv.

------------------------------------------------------------------------

# 9. SPICE Ephemeris Integration

To improve realism, load real planetary data.

Library:

spiceypy

Capabilities:

-   planetary positions
-   moon positions
-   spacecraft ephemerides

------------------------------------------------------------------------

## Example Usage

spiceypy.spkezr()

Returns position and velocity vectors.

------------------------------------------------------------------------

# 10. Lagrange Points and CR3BP

Add restricted three‑body modeling.

Useful for missions involving:

-   Earth‑Sun L1/L2
-   Earth‑Moon L1/L2
-   halo orbits

Important points:

L1\
L2\
L3\
L4\
L5

------------------------------------------------------------------------

## CR3BP Equations

x'' − 2y' = ∂Ω/∂x\
y'' + 2x' = ∂Ω/∂y\
z'' = ∂Ω/∂z

Where Ω is the pseudo‑potential.

------------------------------------------------------------------------

# 11. Advanced Visualization

Add specialized plots:

Trajectory plots\
Porkchop plots\
Gravity well visualizations\
SOI boundary displays

Libraries:

-   matplotlib
-   plotly
-   pyvista

------------------------------------------------------------------------

# 12. Phase‑2 Implementation Milestones

## Milestone 1

Implement Lambert solver.

------------------------------------------------------------------------

## Milestone 2

Create porkchop plot generator.

------------------------------------------------------------------------

## Milestone 3

Implement gravity assist targeting.

------------------------------------------------------------------------

## Milestone 4

Add trajectory optimization.

------------------------------------------------------------------------

## Milestone 5

Integrate SPICE planetary ephemeris.

------------------------------------------------------------------------

## Milestone 6

Add Lagrange point / CR3BP dynamics.

------------------------------------------------------------------------

## Milestone 7

Implement advanced visualizations.

------------------------------------------------------------------------

# 13. Final Phase‑2 Capabilities

The simulator will be able to:

-   simulate spacecraft trajectories
-   compute interplanetary transfers
-   generate launch window porkchop plots
-   model gravity assists
-   optimize complex trajectories
-   visualize gravitational environments
-   load real planetary ephemeris data

------------------------------------------------------------------------

# 14. Recommended References

Fundamentals of Astrodynamics --- Bate, Mueller, White

Orbital Mechanics for Engineering Students --- Curtis

Fundamentals of Astrodynamics and Applications --- Vallado
