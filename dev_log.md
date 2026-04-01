### 3/9/2026

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
- Check for SOI interruptions 
- Check for burn time
Update burn:
- Instant:
    - Add all velocity immediately 
- Time-resolved:
    - Add velocity as (F/m_spacecraft)*dt
    - Splice burn between timesteps if there's no easy division
    - Consider writing dt as smaller during high-interest periods (burns) to resolve accuracy
    - Continue each step within the burn duration or dV delivery
- Adjust mass (if known) using Isp equation
Save bodies:
- Write current state to a "state.json" file so it's not just saved internally
- Write lines such as position/vel of all bodies to a .csv 
- Can write to body-specific files like sun_save.csv



 