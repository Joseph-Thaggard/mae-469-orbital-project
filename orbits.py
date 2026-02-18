# Package imports

#import matplotlib 
 
# HIGHLY simplistic to not account for input files yet
# Minimum descriptors of orbits using predefined characteristics
# Eventually use "add_object" or some other method to add an object 
#   class instead of manually defining them. 
# Use lists instead to give characteristics (or lists of lists)

## Simulation parameters

dt = 1
csv_path = "~/orbit_characteristics.csv"
G = 6.6743*10^-11

# Given data: Planet,  AU, eccentricity, i to eccliptic, Ω deg, ω deg,  θdeg 

# Initialize objects

# Fun equations I found: https://kyleniemeyer.github.io/space-systems-notes/orbital-mechanics/two-body-problems.html 
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def rhs(t, z):
    # 3D orbital motion ODE
    mu = 398.6e3 # km^3/s^2
    r = np.sqrt(z[0]**2 + z[1]**2 + z[2]**2)
    dzdt = np.zeros(6)
    dzdt[0] = z[3]
    dzdt[1] = z[4]
    dzdt[2] = z[5]
    dzdt[3] = (-mu/r**3) * z[0]
    dzdt[4] = (-mu/r**3) * z[1]
    dzdt[5] = (-mu/r**3) * z[2]
    return dzdt
    
r0 = [, , ,]
v0 = [, , ,]
T = .

sol = solve_ivp(rhs, [0, 2*T], np.array(r0 + v0))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot3D(sol.y[0,:], sol.y[1,:], sol.y[2.:])

# Need to add if nbody = true or nbody = false at the end
# Consider reformating .csv to JSON
# Presume we store values like:
# [mass],[position],[velocty],[acceleration],nbody
# I'm a little dumb and forgot how to do python nested vectors

# Force calculator for nbody-enable objects (spacecraft), option 1
for object in sim:
    # Checking for nbody simulation check
    if sim[object][end] == True:
        for object2 in sim:
            # Get position
            target = object
            pos1 = sim[target][1]
            pos2 = sim[object2][1]
            mass1 = sim[target][0]
            mass2 = sim[object2][0]
            dpos = pos2-pos1
            force = G*(mass2*mass1)/(dpos^2)

# Gravitational potential calculator for all objects: 
# Will be calculated at each timestep for all objects
# Calculated at each cell given the current position
# Equation of gravitation potential: -GM/r
# Each object adds to the potential map
