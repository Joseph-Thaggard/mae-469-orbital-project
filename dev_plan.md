**MAE 469/569: Intro to Astrodynamics**  
**Project Team: [insert name here]**

Goal is to create a Europa Clipper style orbital simulation that can do calculations for an Earth-Mars-Earth-Jupiter flight plan. Two propulsive burns are needed; one from Earth to Mars, one from Mars to the Earth-Jupiter gravity assist. The Jupiter approach can be blunt, seen as an impact rather as an orbital insertion. 

**Development Plan: Required**
Create a porkchop plot for the Earth-Mars orbit for the 2026-2046 period, though the function to generate this should be specific to target position (relative to Mars) and velocity for flyby. Can be two functions: transfer-plot and flyby-plot.
Could be second use of porkchop plots, which would be something like maneuver-plot. Difference is a loosened criteria for starting state vectors, though it’s intrinsically similar to a standard transfer-plot, as both are on orbits, just the spacecraft is non-standard. Maneuver-plot should also consider flybys without propulsive changes.
Consider backwards-planning the DSM burn and forward-planning the Earth departure. Should have a way of using porkchop plots (or dozens of them) to find an ideal solution.
Inital and final altitudes are provided. 
Need to (want to) create a global gravitational potential map and vector field for gradient gravitational potentail
Force can be: -m*del(G)
Planetary bodies are considered 2-body entities so they can use fixed orbits 
Kepler's laws and other equations can be used to propogate orbits or generate time-location fetch data
For 2026 to 2046, that's considered 630800000sec (630.8M), so dt terms should be flexible or no longer than 1sec
1sec minimum dt for fixed cases
Flexible cases can use 0.1s for propulsive burns and variable dt outside (think KSP)
Presume no atmosphere
Calculation checks should include radius of each planet to prevent collisions 
Spacecraft should use F=-m*del(G) per each step
Spacecraft actions should be scheduled, manually entered at checkpoints, or paused and entered
Deep-space maneuvers (DSMs) can use a standard orbit transfer 
Look into flyby mechanics
Position tracking can be done with a "states.csv" or other array file, such as "positions.csv", "velocities.csv", etc. Breaking each element into a different file is an option to keep it human-readable, though .pvd or other ParaView-linking files will be useful for visualization in bulk. JSON files can work in an "outputs" folder. 


**Optional/Fun Things**
Consider burns in existing gravity wells
Determine local sphere of infuence definitions and implications for orbial sims