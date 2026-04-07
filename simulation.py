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
import visualization as v
import diags as d
import mission as m


dt = 1e6
max_steps = int(1e8)

sun = bodies.Body('Sun', 1.989e30, 696340e3, np.array([0,0,0]), np.array([0,0,0]))
earth = bodies.Body('Earth', 5.972e24, 6371e3, np.array([1.496e11, 0, 0]), np.array([0, 29780, 0]))
mars = bodies.Body('Mars', 0.64171e24, 3389.5e3, np.array([2.279e11, 0, 0]), np.array([0, 24077, 0]))
ship = bodies.Spacecraft.in_orbit('Spaceship', 1000, 10, earth, 400e3, extra_v=0)  # circular LEO; set extra_v > 3176 to escape Earth SOI
bodies_list = [sun, earth, mars, ship]

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
spacecraft_list = [ship]
m.instant_burn(ship, np.array([0, 0, 2000]))  # example burn: +500 m/s in tangential direction for 1e6 s (should be ~3333 m/s actual delta-v)
# dt_sub = dt / N_SOI_SUBSTEPS when inside an SOI (100 s for dt=1e6, n=10000)
# SOI_RENDER_EVERY: render SOI plot every N sub-steps (~once per LEO orbit at 59)
N_SOI_SUBSTEPS  = 10000
SOI_RENDER_EVERY = 10

for i in range(1,max_steps):
    t_sim = i * dt
    soi_map = p.propagate_soi(bodies_list, spacecraft_list, grid, dt,
                               n_substeps=N_SOI_SUBSTEPS,
                               render_every=SOI_RENDER_EVERY,
                               soi_render_fn=v.render_body_frame)
    v.render_grid(grid, bodies_list)
    d.log_orbits(bodies_list, soi_map)


