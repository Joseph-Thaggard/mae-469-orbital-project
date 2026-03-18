import numpy as np
import numpy.linalg as la

def potential_2body(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    V = np.zeros(len(bodies))  # Initialize potential array
    for i in range(1, len(bodies)):
        ri = bodies[i].position - bodies[0].position  # Vector from the central body to the current body
        M0 = bodies[0].mass  # Mass of the central body
        V[i] -= G*M0/la.norm(ri)  # Potential on body i due to the central body
    
    return V


def potential_nbody(bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    V = np.zeros(len(bodies))  # Initialize potential array
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            print(f"Computing potential between {bodies[i].name} and {bodies[j].name}")
            print(f"Positions: {bodies[i].position} and {bodies[j].position}")
            rij = bodies[j].position - bodies[i].position  # Vector from body i to body j
            Mj = bodies[j].mass  # Mass of body j
            Mi = bodies[i].mass  # Mass of body i
            V[i] -= G*Mj/la.norm(rij)  # Potential on body i due to body j
            V[j] -= G*Mi/la.norm(rij)  # Potential on body j due to body i (equal and opposite)
    
    return V

