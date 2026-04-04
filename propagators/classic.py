import numpy as np

def propogate_velocity(bodies, dt):
    for body in bodies:
        body.velocity += body.acceleration * dt
        #print(f"Updated velocity of {body.name}: {body.velocity}")


def propogate_position(bodies, dt):
    for body in bodies:
        body.position += body.velocity * dt
        #print(f"Updated position of {body.name}: {body.position}")