import numpy as np
class burn:
    def __init__(self, spaceship):
        self.spaceship = spaceship
        self.burn_schedule = np.array([])
        self.burn_vectors = np.array([])

    def apply(self):
        """Apply the burn to the spaceship."""
        # Finish later
        pass

def instant_burn(spaceship,dv_vector):
    """Apply an instantaneous burn to the spaceship, changing its velocity by a specified delta-v."""
    spaceship.velocity += dv_vector
    print(f"Applied instant burn: Δv = {dv_vector} m/s. New velocity: {spaceship.velocity} m/s.")
    pass

def add_burn(spaceship,dv_vector,burn_start):
    """Schedule a burn to be applied at a specific simulation time."""
    spaceship.burn_schedule = np.append(spaceship.burn_schedule, burn_start)
    spaceship.burn_vectors = np.vstack([spaceship.burn_vectors, dv_vector])
    pass

def burn_schedule(spaceship,t_sim):
    """Schedule of burns to apply at specific simulation times."""
    for burn in spaceship.burn_schedule:
        if t_sim == burn:
            # Apply corresponding burn vector from burn.burn_vectors
            instant_burn(spaceship, spaceship.burn_vectors[spaceship.burn_schedule == burn][0])
    pass