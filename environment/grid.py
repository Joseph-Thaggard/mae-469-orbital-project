import numpy as np
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
                    self.potential[i,j,k]=potential.potential_nbody(bodies)