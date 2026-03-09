import numpy as np
#import matplotlib
## Forming a grid of 10 points between 0 and 100
# Number of points in grid per dimension
dx = 10
dy = 10
dz = 10

# Setting bounds of grid (arbitrary 0 to 100)
xmin = 0
ymin = 0
zmin = 0
xmax = 100
ymax = 100
zmax = 100
# Establishing single dimension size
xvec = np.linspace(xmin,xmax,dx)
print(len(xvec))
yvec = np.linspace(ymin,ymax,dy)
zvec = np.linspace(zmin,zmax,dz)

# Building mesh directly from methods shown in xyz coords, prefilled with zeros
x, y, z = np.meshgrid(np.zeros(len(xvec)),np.zeros(len(xvec)),np.zeros(len(xvec)), indexing='ij')
# Stacking grids for a single object
idxmesh = np.stack([x,y,z],axis=-1)

# Add value of 10 at 5,5,5
print(idxmesh[5,5,5])
idxmesh[5,5,5] = 10
print(idxmesh[5,5,5])
slice = idxmesh[:,:,5]
print(slice)
print(np.size(idxmesh))

