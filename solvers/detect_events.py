import numpy as np

def check_event(bodies, spacecraft, t_sim):
    """Check for events such as entering/exiting SOI, collisions, or burn execution."""
    events = []
    
    # Check for SOI entry/exit
    for body in bodies:
        r_vec = spacecraft.position - body.position
        r_mag = np.linalg.norm(r_vec)
        soi_radius = body.radius * (body.mass / bodies[0].mass)**(2/5)  # Simplified SOI radius
        
        if r_mag < soi_radius:
            events.append(f"Entered SOI of {body.name} at t={t_sim:.2f} s")
        else:
            events.append(f"Exited SOI of {body.name} at t={t_sim:.2f} s")
    
    # Check for collisions
    for body in bodies:
        r_vec = spacecraft.position - body.position
        r_mag = np.linalg.norm(r_vec)
        
        if r_mag < (spacecraft.radius + body.radius):
            events.append(f"Collision with {body.name} at t={t_sim:.2f} s")
    
    return events
    
    