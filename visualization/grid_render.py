import atexit
from collections import deque
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  # noqa: F401 — registers 3D projection
from environment.grid import Grid

TRAIL_MAX = 100

# --- global grid view state ---
_fig = None
_ax  = None
_trail = {}   # {body.name: {'indices': deque, 'color': ...}}

# --- per-body frame view state ---
_body_figs = {}   # {focus_body.name: {'fig': fig, 'ax': ax, 'trail': {body.name: {'pts': deque, 'color': ...}}}}

def _keep_open():
    any_open = _fig is not None or bool(_body_figs)
    if any_open:
        plt.ioff()
        plt.show()  # blocking — holds all windows open until user closes them

atexit.register(_keep_open)


def render_grid(grid, bodies):
    global _fig, _ax

    if _fig is None:
        plt.ion()
        _fig = plt.figure()
        _ax = _fig.add_subplot(111, projection='3d')

    _ax.cla()

    for body in bodies:
        # Spacecraft inside an SOI are managed by check_soi — skip them here
        if not getattr(body, 'track_on_solar_grid', True):
            continue
        body.assign_grid_properties(grid)  # refresh from current physical position
        ix, iy, iz = body.grid_index
        if ix < 0 or ix >= grid.nx or iy < 0 or iy >= grid.ny or iz < 0 or iz >= grid.nz:
            continue  # body is outside the grid bounds — don't plot

        sc = _ax.scatter(ix, iy, iz, s=80, label=body.name)
        color = sc.get_facecolor()[0]
        _ax.text(ix, iy, iz, f"  {body.name}", fontsize=8)

        if body.name not in _trail:
            _trail[body.name] = {'indices': deque(maxlen=TRAIL_MAX), 'color': color}
        _trail[body.name]['indices'].append([ix, iy, iz])
        _trail[body.name]['color'] = color

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
    plt.pause(0.001)


def render_body_frame(focus_body, bodies, radius):
    """Render a separate 3D plot in physical coordinates centered on focus_body.
    Only bodies within `radius` meters of focus_body are shown.
    Intended to be called by the SOI detector or manually with a known radius.

    Args:
        focus_body: the body to center the view on
        bodies:     full list of simulation bodies
        radius:     half-width of the view in meters (e.g. SOI radius)
    """
    global _body_figs

    if focus_body.name not in _body_figs:
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d')
        _body_figs[focus_body.name] = {'fig': fig, 'ax': ax, 'trail': {}}

    state = _body_figs[focus_body.name]
    ax    = state['ax']
    ax.cla()

    cx, cy, cz = focus_body.position

    for body in bodies:
        rx = body.position[0] - cx
        ry = body.position[1] - cy
        rz = body.position[2] - cz

        # Skip bodies outside the view radius
        if la.norm([rx, ry, rz]) > radius:
            continue

        sc    = ax.scatter(rx, ry, rz, s=80, label=body.name)
        color = sc.get_facecolor()[0]
        ax.text(rx, ry, rz, f"  {body.name}", fontsize=8)

        # Trail in focus body's frame
        if body.name not in state['trail']:
            state['trail'][body.name] = {'pts': deque(maxlen=TRAIL_MAX), 'color': color}
        state['trail'][body.name]['pts'].append([rx, ry, rz])
        state['trail'][body.name]['color'] = color

        trail_pts = np.array(state['trail'][body.name]['pts'])
        if len(trail_pts) > 1:
            ax.scatter(trail_pts[:-1, 0], trail_pts[:-1, 1], trail_pts[:-1, 2],
                       s=10, color=color, alpha=0.4, zorder=1)

        speed = la.norm(body.velocity)
        if speed > 0:
            v_hat     = body.velocity / speed
            arrow_len = radius * 0.1
            vx, vy, vz = v_hat * arrow_len
            ax.quiver(rx, ry, rz, vx, vy, vz, color='gray', arrow_length_ratio=0.3)

    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(-radius, radius)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title(f'{focus_body.name} frame  radius={radius:.2e} m')
    ax.legend()
    plt.tight_layout()
    plt.pause(0.001)
