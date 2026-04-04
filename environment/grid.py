import numpy as np
import numpy.linalg as la
import diags as d
import bodies as b
class Grid: 
    def __init__(self,ncell,spacing,geometry): 
        self.nx=ncell[0]
        self.ny=ncell[1]
        self.nz=ncell[2]
        self.spacing=spacing
        self.geometry=geometry
        self.potential=np.zeros((self.nx,self.ny,self.nz))
        
    def place_bodies(self,bodies):
        for body in bodies:
            body.assign_grid_properties(self)
            #print(f"Placing {body.name} at grid index {body.grid_index} corresponding to position {body.grid_position}")
            # Check if the body index is within the grid bounds
            if np.any(body.grid_index < 0) or np.any(body.grid_index >= [self.nx, self.ny, self.nz]):
                print(f"Warning: {body.name} is out of grid bounds at index {body.grid_index}. Skipping placement.")
                continue
            
        
    def compute_potential(self,bodies):
        import physics.potential as potential
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    # Presumes scalar at each grid point, remove np.norm for vector use later
                    self.potential[i,j,k]=la.norm(potential.potential_grid(self,np.array([i,j,k]), bodies))
                    #d.log(f"Computed potential at grid point ({i}, {j}, {k}): {self.potential[i,j,k]}")
        return self.potential