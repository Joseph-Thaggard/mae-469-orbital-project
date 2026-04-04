import numpy.linalg as la

def gravity_2body(bodies):
    """Compute gravitational acceleration on each body due to bodies[0] (the central body, e.g. Sun).
    Resets acceleration on all bodies before computing."""
    import numpy as np
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    central = bodies[0]
    for body in bodies:
        body.acceleration = np.zeros(3)  # reset each timestep
    for body in bodies[1:]:
        r = body.position - central.position  # vector from central to body
        r_mag = la.norm(r)
        body.acceleration -= G * central.mass * r / r_mag**3  # points toward central body

def gravity_nbody(bodies):
    """Compute gravitational acceleration on every body due to all other bodies.
    Resets acceleration on all bodies before computing."""
    import numpy as np
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    for body in bodies:
        body.acceleration = np.zeros(3)  # reset each timestep
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            rij = bodies[j].position - bodies[i].position  # vector from i to j
            r_mag = la.norm(rij)
            bodies[i].acceleration += G * bodies[j].mass * rij / r_mag**3  # i pulled toward j
            bodies[j].acceleration -= G * bodies[i].mass * rij / r_mag**3  # j pulled toward i
