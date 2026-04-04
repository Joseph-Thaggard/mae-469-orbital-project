import os
import numpy as np

def log_start():
    if os.path.exists("output_log.txt"):
        os.remove("output_log.txt")
    with open("output_log.txt", "w") as log_file:
        log_file.write("Log started\n")

def log(string):
    with open("output_log.txt", "a") as log_file:
        log_file.write(string + "\n")
        
def log_orbits(bodies):
    if os.path.exists("output_orbits.txt"):
        os.remove("output_orbits.txt")
    with open("output_orbits.txt", "w") as log_file:
        log_file.write("Orbit log started\n")
        log_file.write("Bodies List: " + ", ".join([body.name for body in bodies]) + "\n")
        log_file.write("Body Name, Apoapsis (m), Periapsis (m), Semi-Major Axis (m), Eccentricity\n")
        mu = bodies[0].mu  # gravitational parameter of central body
        for body in bodies[1:]:  # skip central body — it has no orbit around itself
            r = body.position - bodies[0].position  # vector from central body to current body
            r_mag = np.linalg.norm(r)
            v_mag = np.linalg.norm(body.velocity)
            a = 1 / (2/r_mag - v_mag**2/mu)  # semi-major axis via vis-viva
            h = np.linalg.norm(np.cross(r, body.velocity))  # specific angular momentum
            e = np.sqrt(1 - h**2 / (mu * a))  # eccentricity
            periapsis = a * (1 - e)
            apoapsis = a * (1 + e)
            log_file.write(f"{body.name}, {apoapsis:.2e}, {periapsis:.2e}, {a:.2e}, {e:.4f}\n")