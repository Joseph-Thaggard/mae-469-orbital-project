import numpy as np
import numpy.linalg as la

def check_soi(bodies, spacecraft_list, grid):
    """Determine which body's sphere of influence each spacecraft is in.

    SOI radius formula: r_soi = a * (m_planet / m_central)^(2/5)

    Sets spacecraft.track_on_solar_grid:
        True  — spacecraft is outside all SOIs; grid_index is updated for solar plot
        False — spacecraft is inside a body's SOI; grid_index is set to None (hidden from solar grid)

    Args:
        bodies:          full bodies list, bodies[0] is the central body (Sun)
        spacecraft_list: list of Spacecraft objects to check
        grid:            the SpaceGrid, used to compute grid_index when outside all SOIs

    Returns:
        dict {spacecraft.name: (body, r_soi) if inside an SOI, else (None, None)}
    """
    central    = bodies[0]
    mu_central = central.mu
    planets    = [b for b in bodies[1:] if b not in spacecraft_list]

    # Compute each planet's SOI radius from its current orbital elements
    soi_radii = {}
    for planet in planets:
        r     = planet.position - central.position
        r_mag = la.norm(r)
        v_mag = la.norm(planet.velocity)
        a     = 1.0 / (2.0/r_mag - v_mag**2/mu_central)   # vis-viva semi-major axis
        if a <= 0:
            continue  # hyperbolic — no bounded SOI
        r_soi = a * (planet.mass / central.mass) ** (2.0/5.0)
        soi_radii[planet.name] = (planet, r_soi)

    result = {}
    for sc in spacecraft_list:
        # Find the body whose SOI the spacecraft is deepest inside (smallest d/r_soi ratio)
        best_body  = None
        best_ratio = 1.0  # threshold: ratio < 1 means inside SOI

        for name, (planet, r_soi) in soi_radii.items():
            d     = la.norm(sc.position - planet.position)
            ratio = d / r_soi
            if ratio < best_ratio:
                best_ratio = ratio
                best_body  = planet

        if best_body is not None:
            sc.track_on_solar_grid = False
            sc.grid_index          = None
            result[sc.name]        = (best_body, soi_radii[best_body.name][1])
        else:
            sc.track_on_solar_grid = True
            sc.assign_grid_properties(grid)
            result[sc.name]        = (None, None)

    return result


def propagate_soi(bodies, spacecraft_list, grid, dt,
                  n_substeps=10000, render_every=59, soi_render_fn=None):
    """Check SOI membership and propagate all bodies for one timestep.

    When any spacecraft is inside a body's SOI, dt is subdivided into n_substeps
    smaller steps to resolve near-body orbits. Otherwise a single full-dt step is used.

    Args:
        bodies:          full bodies list, bodies[0] is the central body
        spacecraft_list: list of Spacecraft objects
        grid:            SpaceGrid for SOI detection and grid index assignment
        dt:              outer timestep (s)
        n_substeps:      sub-step count when inside an SOI (dt_sub = dt / n_substeps)
        render_every:    call soi_render_fn every N sub-steps (default ~once per LEO orbit)
        soi_render_fn:   optional callable(focus_body, bodies, view_radius) for SOI plots

    Returns:
        soi_map: dict {spacecraft.name: (body_or_None, r_soi_or_None)}
    """
    from .gravity import gravity_nbody
    from propagators.classic import propogate_velocity, propogate_position

    soi_map = check_soi(bodies, spacecraft_list, grid)
    in_soi  = any(body is not None for body, _ in soi_map.values())

    if in_soi:
        dt_sub = dt / n_substeps
        for sub in range(n_substeps):
            gravity_nbody(bodies)
            propogate_velocity(bodies, dt_sub)
            propogate_position(bodies, dt_sub)
            if soi_render_fn is not None and sub % render_every == 0:
                for sc in spacecraft_list:
                    focus_body, _ = soi_map[sc.name]
                    if focus_body is not None:
                        sc_dist     = la.norm(sc.position - focus_body.position)
                        view_radius = max(sc_dist * 5, focus_body.radius * 2)
                        soi_render_fn(focus_body, bodies, view_radius)
    else:
        gravity_nbody(bodies)
        propogate_velocity(bodies, dt)
        propogate_position(bodies, dt)

    return soi_map
