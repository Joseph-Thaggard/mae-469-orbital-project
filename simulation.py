import numpy as np
import inspect
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Import code files
import bodies
import environment
from environment.grid import Grid
import mission
import physics as p
import propagators
import visualization as v
import diags as d

dt = 1000

sun = bodies.Body('Sun', 1.989e30, 696340e3, np.array([0,0,0]), np.array([0,0,0]))
earth = bodies.Body('Earth', 5.972e24, 6371e3, np.array([1.496e11, 2e11, 0]), np.array([0, 29780, 0]))
mars = bodies.Body('Mars', 0.64171e24, 3389.5e3, np.array([2.279e11, 4e11, 0]), np.array([0, 24077, 0]))
d.log_start()
bodies_list = [sun, earth, mars]

# Input list, currently not used
xmin = -5e11  # minimum x coordinate in meters
xmax = 5e11   # maximum x coordinate in meters
ymin = -5e11  # minimum y coordinate in meters
ymax = 5e11   # maximum y coordinate in meters  
zmin = -5e11  # minimum z coordinate in meters
zmax = 5e11   # maximum z coordinate in meters
xbound = [xmin,xmax]
ybound = [ymin,ymax]
zbound = [zmin,zmax]
geometry = np.array([xbound, ybound, zbound])  # Physical size of the grid in meters

# Input list, currently used
spacing = 1e10  # spacing between grid points in meters (5 million km)
ncell_x = 100  # number of grid cells in x direction
ncell_y = 100  # number of grid cells in y direction
ncell_z = 100  # number of grid cells in z direction
ncell = np.array([ncell_x, ncell_y, ncell_z])   # number of grid cells in each direction (sim resolution)
grid = environment.Grid(ncell, spacing, geometry)
grid.place_bodies(bodies_list)
#p.gravity_2body(bodies_list)
#dG = Grid.compute_potential(grid,bodies_list)
print(f"Potential grid shape: {grid.potential.shape}")
v.render_grid(grid, bodies_list)

# Heatmap of gravitational potential at z=0 slice
#z_slice = dG[:, :, 0]  # shape (nx, ny)
#fig, ax = plt.subplots()
#c = ax.pcolormesh(z_slice.T, cmap='coolwarm', norm=matplotlib.colors.LogNorm())
#fig.colorbar(c, ax=ax, label='Gravitational Potential (J/kg)')
#ax.set_xlabel('X grid index')
#ax.set_ylabel('Y grid index')
#ax.set_title('Gravitational Potential at z=0')
#plt.tight_layout()
#plt.show()

# Print all attributes of the sun object
attrs = vars(sun)
print(', '.join("%s: %s" % item for item in attrs.items()))

