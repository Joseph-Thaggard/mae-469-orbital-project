import numpy as np
import numpy.linalg as la
import diags as d
class SpaceGrid: 
    def __init__(self,nx,ny,nz,spacing): 
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.spacing=spacing
        self.potential=np.zeros((nx,ny,nz))
        
        
    def compute_potential(self,bodies):
        import physics.potential as potential
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    # Presumes scalar at each grid point, remove np.norm for vector use later
                    self.potential[i,j,k]=la.norm(potential.potential_grid(self,np.array([i,j,k]), bodies))
                    d.log(f"Computed potential at grid point ({i}, {j}, {k}): {self.potential[i,j,k]}")
        return self.potential