### 4/7/2026

**Progress audit against dev_plan.md and phase2_roadmap.md**

---

**COMPLETED from dev_plan.md (Required):**

- Body classes (`Body`, `Spacecraft`) with name, mass, velocity, position, radius, mu — `bodies/`
- Spacecraft orbital placement factory `Spacecraft.in_orbit` using parent body radius to avoid placement errors
- Gravitational potential field computation on grid (`Grid.compute_potential`, `physics/potential.py`)
- Potential field heatmap visualization (pcolormesh with LogNorm, coded in simulation.py, commented out)
- N-body propagation for all bodies including spacecraft (`physics/gravity.py → gravity_nbody`)
- 2-body propagation (Sun-fixed) (`gravity_2body`)
- Symplectic Euler integrator (`propagators/classic.py`)
- SOI detection and automatic frame switching (`physics/soi.py → check_soi`)
- SOI-aware sub-stepping (`propagate_soi`) — coarse dt for interplanetary, fine dt_sub inside SOI
- Instant burns (`mission/burns.py → instant_burn`)
- Scheduled burns (`add_burn`, `burn_schedule`) — implemented, not yet wired into simulation loop
- Time-resolved burn class stub (`burn.apply`) — placeholder, not yet implemented
- Collision/SOI event detector stub (`solvers/detect_events.py → check_event`) — exists, not wired into loop
- Interactive 3D solar-grid visualization with body trails (`visualization/grid_render.py → render_grid`)
- Body-frame SOI visualization with trail history (`render_body_frame`)
- Step-by-step orbit log with position, velocity, SOI membership (`diags/log.py → log_orbits`)
- Orbit completion detection (heliocentric) with apoapsis/periapsis/eccentricity/semi-major axis
- Spacecraft thrust attribute added (`spacecraft.thrust = 1500 N` placeholder)
- Simulation time variable (`t_sim = i * dt`) tracked in loop

**PARTIALLY DONE:**

- Burn scheduling: `add_burn`/`burn_schedule` exist in `mission/burns.py` but `burn_schedule` is not called each step in `simulation.py` loop — burns will not execute automatically unless wired in
- Event detection: `solvers/detect_events.py` has `check_event` but it has a bug (uses `body.radius * (mass ratio)^(2/5)` instead of orbital SOI formula from `check_soi`), and is not called in the main loop
- State saving: no `states.csv` / `positions.csv` / JSON output yet — only `output_log.txt` and `output_orbits.txt`
- Kepler propagator: planned in dev_plan, not yet implemented — all propagation is currently numerical n-body
- RK4 integrator: planned in dev_plan, not yet implemented — using symplectic Euler only

**NOT STARTED from dev_plan.md (Required):**

- Porkchop plot (Earth-Mars 2026-2046) — no `mission_design/` directory exists yet
- Lambert solver — not implemented
- Gravity assist / flyby mechanics (`mission/flyby.py` is empty)
- DSM burn planning
- Variable dt (KSP-style adaptive stepping outside SOI) — currently fixed dt
- Atmospheric exclusion radius check for collisions — `check_event` checks radius but is not called
- State file output (`states.csv`, `positions.csv`, `.pvd`, JSON) — `output/save.py` is a stub
- Kepler's equation solver (`KeplerNR`, `M2TA`) — planned in dev_plan, not started
- Trajectory plot (orbital ellipse preview from Kepler elements)

**NOT STARTED from phase2_roadmap.md:**

- Lambert solver (`mission_design/lambert.py`)
- Porkchop plot generator (`mission_design/porkchop.py`)
- Gravity assist targeting (`mission_design/flyby_targeting.py`)
- Trajectory optimization (differential evolution)
- SPICE ephemeris integration (`spiceypy`)
- Lagrange point / CR3BP modeling
- Advanced visualization (plotly, pyvista)

---

**Immediate priorities suggested by current state:**
1. Wire `burn_schedule` into the simulation loop so scheduled burns execute automatically
2. Wire `check_event` (or replace with `check_soi` logic) into the loop for collision detection
3. Implement state file output so simulation data persists beyond terminal/plot
4. Implement Lambert solver — blocks porkchop plots and DSM planning

---

### 4/1/2026

Added two .md files based on a ChatGPT conversation on potential methods we can use for this project. Some of the "advanced" techniques will be required as the conversation didn't include context of the project. 

### 4/1/2026

Making a general plan as to how we can build this project. Here are the thoughts:

Developer notes:
- When calling an attribute, check what function set it's from. For example, Earth's velocity will be stored as `earth[4]` if given a variable `earth` that is declared with `earth = bodies.Body('Earth', 5.972e24, 6371e3, np.array([1.496e11, 0, 0]), np.array([0, 29780, 0]))`. 

Initialization:
- Add bodies
    - Includes name, mass, velocity, position, radius, SOI
    - Uses defined class structure in /bodies/body.py
    - Spacecraft version is /bodies/spacecraft.py
    - These can be called by using "bodies.Body" where "bodies" is the folder and "Body" is the class 
- Add fields / diagnostic maps
    - Calculate potential 
    - Plot positions and trace vectors for current orbit for Sun-centered (kepler prediction)
Set t=0 (adjust for epoch if needed)
Propogate bodies:
- Kepler:
    - Use timestep (t+dt) to determine position of each body
    - Use solvers (mean to true anomaly) as functions to reduce error complexity (function can be called M2TA)
    - Make a dedicated Newton-Raphson for this problem (KeplerNR)
- Spacecraft-Nbody:
    - Use RK4 to determine positions based on current velocity. Use t+dt.
- True Nbody
    - Do spacecraft-Nbody for all bodies
Update state:
- Update positions
- Update velocities
Detect events:
- Check for SOI:
    - If within SOI, note within events
    - Swap from NBody to Kepler for spacecraft propogation, with body as fixed
- Check for burn time
Update burn:
- Instant:
    - Add all velocity immediately 
- Time-resolved:
    - Add velocity as (F/m_spacecraft)*dt
    - Splice burn between timesteps if there's no easy division
    - Consider writing dt as smaller during high-interest periods (burns) to resolve accuracy
    - Continue each step within the burn duration or dV delivery
- Adjust spacecraft mass (if known) using Isp equation
Add fields:
- Potential field
- Trajectory plot (with SOI interruptions)
Save bodies:
- Write current state to a "state.json" file so it's not just saved internally
- Write lines such as position/vel of all bodies to a .csv 
- Can write to body-specific files like sun_save.csv
Advance time:
- t=t+dt

In a bit more of a literal example:
```python 
state = sim.init(initial_state)
sim.start()     #indicating some start/end conditions
    sim.propagate(state)
    sim.update(state)
    sim.detect(state)
    sim.action(state)
    sim.fields(state)
    sim.save(state)
    sim.advance(state)
sim.stop()
```




 