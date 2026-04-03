import numpy as np
import numpy.linalg as la
import bodies
import diags as d

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
            #print(f"Computing potential between {bodies[i].name} and {bodies[j].name}")
            #print(f"Positions: {bodies[i].position} and {bodies[j].position}")
            rij = bodies[j].position - bodies[i].position  # Vector from body i to body j
            d.log(f"Vector rij: {rij}, Norm: {la.norm(rij)}")
            Mj = bodies[j].mass  # Mass of body j
            Mi = bodies[i].mass  # Mass of body i
            V[i] -= G*Mj/la.norm(rij)  # Potential on body i due to body j
            V[j] -= G*Mi/la.norm(rij)  # Potential on body j due to body i (equal and opposite)
    #print(f"Computed potentials: {V}")
    return V

def potential_grid(grid,pos, bodies):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    point = np.array([
        grid.geometry[0][0] + pos[0]*grid.spacing,
        grid.geometry[1][0] + pos[1]*grid.spacing,
        grid.geometry[2][0] + pos[2]*grid.spacing
    ])  # Physical position of the grid point
    V = 0  # Initialize potential at this grid point
    for body in bodies:
        if np.all(point-body.position <= body.radius):
            d.log(f"Warning: Grid point coincides within the size of {body.name}. Assigning to zero.")
            V -= 0  # Assign zero potential if the grid point coincides with the body's position
            continue
        r = point - body.position  # Vector from the body to the grid point
        V -= G*body.mass/la.norm(r)  # Potential contribution from this body
    return V

