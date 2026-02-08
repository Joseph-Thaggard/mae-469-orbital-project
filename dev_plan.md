**MAE 469/569: Intro to Astrodynamics**  
**Project Team: [insert name here]**

Goal is to create a Europa Clipper style orbital simulation that can do calculations for an Earth-Mars-Earth-Jupiter flight plan. Two propulsive burns are needed; one from Earth to Mars, one from Mars to the Earth-Jupiter gravity assist. The Jupiter approach can be blunt, seen as an impact rather as an orbital insertion. 

**Development Plan**
Look at a forward-step (accel-vel-dist) code (assess, may not be required)  
Place all major planets and verify state vectors for small and large durations (i+1s, i+1e9s)  
Create functions to bounce between orbital elements and r and v vectors.   
Force handling may not be required 
Create a porkchop plot for the Earth-Mars orbit for the 2026-2046 period, though the function to generate this should be specific to target position (relative to Mars) and velocity for flyby. Can be two functions: transfer-plot and flyby-plot.
Could be second use of porkchop plots, which would be something like maneuver-plot. Difference is a loosened criteria for starting state vectors, though it’s intrinsically similar to a standard transfer-plot, as both are on orbits, just the spacecraft is non-standard. Maneuver-plot should also consider flybys without propulsive changes.
Consider backwards-planning the DSM burn and forward-planning the Earth departure. Should have a way of using porkchop plots (or dozens of them) to find an ideal solution.
This sets the initial and final altitudes. 