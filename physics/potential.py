import numpy as np

def potential_2body(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    V = np.zeros(len(bodies))  # Initialize potential array
    for i in range(1, len(bodies)):
        ri = bodies.positions[i] - bodies.positions[0]  # Vector from the central body to the current body
        M0 = bodies.masses[0]  # Mass of the central body
        V[i] -= G*M0/ri  # Potential on body i due to the central body
    
    return V


def potential_nbody(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    V = np.zeros(len(bodies))  # Initialize potential array
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            rij = bodies.positions[j] - bodies.positions[i]  # Vector from body i to body j
            Mj = bodies.masses[j]  # Mass of body j
            Mi = bodies.masses[i]  # Mass of body i
            V[i] -= G*Mj/rij  # Potential on body i due to body j
            V[j] -= G*Mi/rij  # Potential on body j due to body i (equal and opposite)
    
    return V

