import numpy as np
import numpy.linalg as la

def gravity_2body(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    a = np.zeros_like(bodies.positions)  # Initialize acceleration array
    for i in range(1, len(bodies)):
        ri = bodies.positions[i] - bodies.positions[0]  # Vector from the central body to the current body
        M0 = bodies.masses[0]  # Mass of the central body
        a[i] += G*M0*ri/(la.norm(ri)**2)
    
    return a

def gravity_nbody(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    a = np.zeros_like(bodies.positions)  # Initialize acceleration array
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            rij = bodies.positions[j] - bodies.positions[i]  # Vector from body i to body j
            Mj = bodies.masses[j]  # Mass of body j
            Mi = bodies.masses[i]  # Mass of body i
            a[i] += G*Mj*rij/(la.norm(rij)**2)  # Acceleration on body i due to body j
            a[j] -= G*Mi*rij/(la.norm(rij)**2)  # Acceleration on body j due to body i (equal and opposite)
    
    return a
