import numpy as np
class Spacecraft: 
    def __init__(self,name,mass,radius,position,velocity):
        self.name=name
        self.mass=mass
        self.radius=radius
        self.position=np.array(position, dtype=float)
        self.velocity=np.array(velocity, dtype=float)
        self.mu=6.67430e-11*mass
        self.acceleration=np.zeros(3)

    def assign_grid_properties(self,grid):
        origin = np.array([grid.geometry[0][0], grid.geometry[1][0], grid.geometry[2][0]])
        self.grid_position = self.position - origin  # position relative to grid origin (0,0,0) corner
        self.grid_index = np.round(self.grid_position / grid.spacing).astype(int)
        


