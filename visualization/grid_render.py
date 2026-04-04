import atexit
from collections import deque
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  # noqa: F401 — registers 3D projection
from environment.grid import Grid

TRAIL_MAX = 200

def _keep_open():
    if _fig is not None:
        plt.ioff()
        plt.show()  # blocking — holds the window open until user closes it

atexit.register(_keep_open)

# Module-level persistent state
_fig = None
_ax = None
_trail = {}   # {body.name: {'indices': [[ix,iy,iz], ...], 'color': ...}}

def render_grid(grid, bodies):
    global _fig, _ax, _trail

    # First call: enable interactive mode and create the window
    if _fig is None:
        plt.ion()
        _fig = plt.figure()
        _ax = _fig.add_subplot(111, projection='3d')

    _ax.cla()  # clear previous frame without closing the window

    for body in bodies:
        body.assign_grid_properties(grid)  # refresh index from current physical position
        ix, iy, iz = body.grid_index

        # Plot current position and capture its color for trail matching
        sc = _ax.scatter(ix, iy, iz, s=80, label=body.name)
        color = sc.get_facecolor()[0]
        _ax.text(ix, iy, iz, f"  {body.name}", fontsize=8)

        # Store trail entry (deque auto-drops oldest when full)
        if body.name not in _trail:
            _trail[body.name] = {'indices': deque(maxlen=TRAIL_MAX), 'color': color}
        _trail[body.name]['indices'].append([ix, iy, iz])
        _trail[body.name]['color'] = color  # update in case color assignment shifts

        # Draw trail dots for all previous positions
        trail_pts = np.array(_trail[body.name]['indices'])
        if len(trail_pts) > 1:
            _ax.scatter(trail_pts[:-1, 0], trail_pts[:-1, 1], trail_pts[:-1, 2],
                        s=10, color=color, alpha=0.4, zorder=1)

        speed = la.norm(body.velocity)
        if speed > 0:
            v_hat = body.velocity / speed
            arrow_len = max(grid.nx, grid.ny, grid.nz) * 0.1
            vx, vy, vz = v_hat * arrow_len
            _ax.quiver(ix, iy, iz, vx, vy, vz, color='gray', arrow_length_ratio=0.3)

    _ax.set_xlim(0, grid.nx)
    _ax.set_ylim(0, grid.ny)
    _ax.set_zlim(0, grid.nz)
    _ax.set_xlabel('X grid index')
    _ax.set_ylabel('Y grid index')
    _ax.set_zlabel('Z grid index')
    _ax.set_title('Body Positions and Velocity Vectors')
    _ax.legend()
    plt.tight_layout()
    plt.pause(0.001)  # non-blocking draw update
