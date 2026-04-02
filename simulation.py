import numpy as np
import inspect
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Import code files
import bodies
import environment
from environment.grid import SpaceGrid
import mission
import physics 
import propagators
import visualization
import diags as d



sun = bodies.Body('Sun', 1.989e24, 696340e3, np.array([0,0,0]), np.array([0,0,0]))
earth = bodies.Body('Earth', 5.972e24, 6371e3, np.array([1.496e11, 2e11, 0]), np.array([0, 29780, 0]))
mars = bodies.Body('Mars', 0.64171e24, 3389.5e3, np.array([2.279e11, 4e11, 0]), np.array([0, 24077, 0]))
d.log_start()
bodies_list = [sun, earth, mars]
grid = environment.SpaceGrid(100, 100, 100, 5e9)
dG = SpaceGrid.compute_potential(grid,bodies_list)
d.log(f"Potential grid shape: {grid.potential.shape}")

# Heatmap of gravitational potential at z=0 slice
z_slice = dG[:, :, 0]  # shape (nx, ny)
fig, ax = plt.subplots()
c = ax.pcolormesh(z_slice.T, cmap='coolwarm', norm=matplotlib.colors.LogNorm())
fig.colorbar(c, ax=ax, label='Gravitational Potential (J/kg)')
ax.set_xlabel('X grid index')
ax.set_ylabel('Y grid index')
ax.set_title('Gravitational Potential at z=0')
plt.tight_layout()
plt.show()

# Print all attributes of the sun object
attrs = vars(sun)
d.log(', '.join("%s: %s" % item for item in attrs.items()))

