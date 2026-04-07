import os
import numpy as np

# --- simple log ---

def log_start():
    if os.path.exists("output_log.txt"):
        os.remove("output_log.txt")
    with open("output_log.txt", "w") as f:
        f.write("Log started\n")

def log(string):
    with open("output_log.txt", "a") as f:
        f.write(string + "\n")

# --- orbit log ---

_orbits_initialized = False
_step = 0
_prev_y  = {}       # body.name -> previous y relative to central body
_orbit_count = {}   # body.name -> completed orbit count

def log_orbits(bodies, soi_map=None):
    global _orbits_initialized, _step, _prev_y, _orbit_count

    mode = "w" if not _orbits_initialized else "a"
    _orbits_initialized = True
    _step += 1

    mu = bodies[0].mu  # gravitational parameter of central body (Sun)

    with open("output_orbits.txt", mode) as f:
        if mode == "w":
            f.write("Orbit Log\n")
            f.write("Bodies: " + ", ".join(b.name for b in bodies) + "\n\n")

        f.write(f"Step {_step}\n")

        for body in bodies[1:]:
            r = body.position - bodies[0].position
            r_mag = np.linalg.norm(r)
            v_mag = np.linalg.norm(body.velocity)

            # SOI annotation for this body
            soi_label = ""
            if soi_map is not None and body.name in soi_map:
                soi_body, _ = soi_map[body.name]
                soi_label = f" [SOI: {soi_body.name}]" if soi_body is not None else " [SOI: interplanetary]"

            # Per-step: position and velocity
            f.write(
                f"  {body.name}{soi_label}: "
                f"pos=[{body.position[0]:.4e}, {body.position[1]:.4e}, {body.position[2]:.4e}] m  "
                f"vel=[{body.velocity[0]:.4e}, {body.velocity[1]:.4e}, {body.velocity[2]:.4e}] m/s  "
                f"|r|={r_mag:.4e} m  |v|={v_mag:.4e} m/s\n"
            )

            # Orbit completion: y relative to Sun crosses 0 upward on the positive-x side
            prev_y = _prev_y.get(body.name)
            if prev_y is not None and prev_y < 0 and r[1] >= 0 and r[0] > 0:
                _orbit_count[body.name] = _orbit_count.get(body.name, 0) + 1
                a = 1.0 / (2.0/r_mag - v_mag**2/mu)
                h = np.linalg.norm(np.cross(r, body.velocity))
                if a > 0:
                    e = np.sqrt(max(0.0, 1.0 - h**2 / (mu * a)))
                    periapsis = a * (1.0 - e)
                    apoapsis  = a * (1.0 + e)
                    orbit_type = "elliptic"
                    apoapsis_str = f"{apoapsis:.4e} m"
                else:
                    e = np.sqrt(1.0 + h**2 / (mu * abs(a)))
                    periapsis = a * (1.0 - e)   # a<0, e>1 → positive periapsis
                    orbit_type = "hyperbolic"
                    apoapsis_str = "inf"
                f.write(
                    f"  *** {body.name} orbit #{_orbit_count[body.name]} complete ({orbit_type}): "
                    f"a={a:.4e} m  e={e:.4f}  "
                    f"periapsis={periapsis:.4e} m  apoapsis={apoapsis_str}\n"
                )
            _prev_y[body.name] = r[1]
