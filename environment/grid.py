import numpy as np
class SpaceGrid: 
    def __init__(self,nx,ny,nz,spacing): 
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.spacing=spacing
        self.potential=np.zeros((nx,ny,nz))