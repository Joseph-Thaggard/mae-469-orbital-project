import numpy as np

class Spacecraft:
    def __init__(self, name, mass, radius, position, velocity):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.mu = 6.67430e-11 * mass
        self.acceleration = np.zeros(3)
        self.thrust = 1500 # Placeholder 1.5kN thrust for burn calculations; set as class property later

    @classmethod
    def in_orbit(cls, name, mass, sc_radius, parent_body, altitude, extra_v=0.0):
        """Place spacecraft in a circular orbit at a given altitude above parent_body's surface.

        Positions the spacecraft along +x from the parent and gives it +y tangential velocity
        (counterclockwise orbit). Use extra_v to add velocity above circular (e.g. for escape).

        Args:
            name:        spacecraft name
            mass:        spacecraft mass (kg)
            sc_radius:   spacecraft physical radius (m)
            parent_body: Body object to orbit — uses parent_body.radius and parent_body.mu
            altitude:    altitude above parent_body's surface (m)
            extra_v:     additional tangential velocity above circular orbit speed (m/s)
                         Earth escape from 400 km: extra_v > ~3176 m/s (v_esc - v_circ)

        Returns:
            Spacecraft instance
        """
        r = parent_body.radius + altitude               # orbital radius from parent center
        v_circ = np.sqrt(parent_body.mu / r)            # circular orbital velocity

        position = parent_body.position + np.array([r, 0.0, 0.0])
        velocity = parent_body.velocity + np.array([0.0, v_circ + extra_v, 0.0])

        return cls(name, mass, sc_radius, position, velocity)

    def assign_grid_properties(self, grid):
        origin = np.array([grid.geometry[0][0], grid.geometry[1][0], grid.geometry[2][0]])
        self.grid_position = self.position - origin
        self.grid_index = np.round(self.grid_position / grid.spacing).astype(int)
