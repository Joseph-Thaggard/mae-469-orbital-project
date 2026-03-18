import numpy as np
import inspect
import matplotlib
import matplotlib.pyplot as plt

# Import code files
import bodies
import environment
from environment.grid import SpaceGrid
import mission
import physics 
import propagators
import visualization


sun = bodies.Body('Sun', 1.989e30, 696340e3, np.array([0,0,0]), np.array([0,0,0]))
earth = bodies.Body('Earth', 5.972e24, 6371e3, np.array([1.496e11, 0, 0]), np.array([0, 29780, 0]))
mars = bodies.Body('Mars', 0.64171e24, 3389.5e3, np.array([2.279e11, 0, 0]), np.array([0, 24077, 0]))

bodies_list = [sun, earth, mars]
grid = environment.SpaceGrid(100, 100, 100, 1e9)
dG = SpaceGrid.compute_potential(grid,bodies_list)
plt.plot(dG)
plt.show()

# Print all attributes of the sun object
attrs = vars(sun)
print(', '.join("%s: %s" % item for item in attrs.items()))

