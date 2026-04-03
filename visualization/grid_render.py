import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d  # noqa: F401 — registers 3D projection
from environment.grid import Grid

def render_grid(grid, bodies):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for body in bodies:
        if body.grid_index is None:
            print(f"Warning: {body.name} has no grid index. Call grid.place_bodies() first.")
            continue

        ix, iy, iz = body.grid_index

        # Plot body position
        ax.scatter(ix, iy, iz, s=80, label=body.name)
        ax.text(ix, iy, iz, f"  {body.name}", fontsize=8)

        # Scale velocity vector to grid units so it's visible relative to grid size
        speed = np.linalg.norm(body.velocity)
        if speed > 0:
            v_hat = body.velocity / speed
            arrow_len = max(grid.nx, grid.ny, grid.nz) * 0.1  # 10% of grid span
            vx, vy, vz = v_hat * arrow_len
            ax.quiver(ix, iy, iz, vx, vy, vz, color='gray', arrow_length_ratio=0.3)

    ax.set_xlim(0, grid.nx)
    ax.set_ylim(0, grid.ny)
    ax.set_zlim(0, grid.nz)
    ax.set_xlabel('X grid index')
    ax.set_ylabel('Y grid index')
    ax.set_zlabel('Z grid index')
    ax.set_title('Body Positions and Velocity Vectors')
    ax.legend()
    plt.tight_layout()
    plt.show()
