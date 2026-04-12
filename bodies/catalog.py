"""
catalog.py — Solar system body definitions at J2000 epoch (2000-Jan-01.5 TT)

Initial states use the mean ecliptic orbital elements from the project spec
(Standish 1992 / JPL Planetary Fact Sheet), converted to heliocentric Cartesian
via elements_to_rv().  All values in SI (meters, m/s).

J2000 epoch definition:
    January 1, 2000 at 12:00:00 TT ≈ January 1, 2000 11:58:56 UTC
    (the project spec rounds this to "11:58 am UTC")

Usage:
    from bodies.catalog import solar_system, make_earth, make_mars, J2000_EPOCH_UTC
    from solvers.solver import kepler_propagate

    earth = make_earth()
    mu_sun = make_sun().mu

    # Propagate to a target time (seconds from J2000):
    dt = (target_datetime - J2000_EPOCH_UTC).total_seconds()
    r_t, v_t = kepler_propagate(earth.position, earth.velocity, dt, mu_sun)

Orbital elements source (project document, Table 1):
    Planet  a (AU)    e         i (deg)   Ω (deg)     ω (deg)      θ (deg)
    Mercury 0.387099  0.205631  7.00487   48.33167    29.12478    174.7944
    Venus   0.723332  0.006773  3.39471   76.68069    54.85229     50.44675
    Earth   1.000000  0.01671   0.00005  -11.26064   114.20783    -2.48284
    Mars    1.523662  0.093412  1.85061   49.57854   286.4623      19.41248
    Jupiter 5.203363  0.048393  1.3053   100.55615   -85.8023      19.55053
    Saturn  9.537070  0.054151  2.48446  113.71504   -21.2831     -42.4876
    Uranus  19.19126  0.047168  0.76986   74.22988    96.73436    142.2679
    Neptune 30.06896  0.008586  1.76917  131.72169   -86.75034    259.9087
    Pluto   39.48169  0.248808  17.14175 110.30347   113.76329     14.86205
"""

import numpy as np
from bodies.body import Body
from solvers.solver import elements_to_rv

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
G   = 6.67430e-11          # m³ kg⁻¹ s⁻²
AU  = 1.495978707e11       # m per AU (IAU 2012 definition)
DEG = np.pi / 180.0        # degrees → radians

# J2000 epoch as a datetime (UTC) for external reference
from datetime import datetime, timezone
J2000_EPOCH_UTC = datetime(2000, 1, 1, 11, 58, 56, tzinfo=timezone.utc)

# Solar gravitational parameter
MU_SUN = 1.32712440018e20  # m³/s²  (IAU 2015)

# ---------------------------------------------------------------------------
# Helper — build Body from J2000 Keplerian elements
# ---------------------------------------------------------------------------

def _body_from_elements(name, mass, radius_m, a_au, e, i_deg, raan_deg, argp_deg, nu_deg):
    """Create a Body at its J2000 heliocentric state vector.

    Args:
        name:       body name string
        mass:       body mass (kg)
        radius_m:   mean radius (m)
        a_au:       semi-major axis (AU)
        e:          eccentricity
        i_deg:      inclination (degrees)
        raan_deg:   RAAN / longitude of ascending node Ω (degrees)
        argp_deg:   argument of periapsis ω (degrees)
        nu_deg:     true anomaly θ at J2000 (degrees)

    Returns:
        Body with position and velocity set to J2000 heliocentric ecliptic values.
    """
    a    = a_au * AU
    i    = i_deg    * DEG
    raan = raan_deg * DEG
    argp = argp_deg * DEG
    nu   = nu_deg   * DEG

    r, v = elements_to_rv(a, e, i, raan, argp, nu, MU_SUN)
    return Body(name, mass, radius_m, r, v)


# ---------------------------------------------------------------------------
# Individual body constructors — J2000 epoch
# ---------------------------------------------------------------------------

def make_sun():
    b = Body('Sun', 1.989e30, 696_340e3,
             np.array([0.0, 0.0, 0.0]),
             np.array([0.0, 0.0, 0.0]))
    b.mu = MU_SUN   # use precise value
    return b

def make_mercury():
    return _body_from_elements(
        'Mercury', 3.301e23, 2_439.7e3,
        a_au=0.387099, e=0.205631, i_deg=7.00487,
        raan_deg=48.33167, argp_deg=29.12478, nu_deg=174.7944)

def make_venus():
    return _body_from_elements(
        'Venus', 4.867e24, 6_051.8e3,
        a_au=0.723332, e=0.006773, i_deg=3.39471,
        raan_deg=76.68069, argp_deg=54.85229, nu_deg=50.44675)

def make_earth():
    return _body_from_elements(
        'Earth', 5.972e24, 6_371e3,
        a_au=1.000000, e=0.01671, i_deg=0.00005,
        raan_deg=-11.26064, argp_deg=114.20783, nu_deg=-2.48284)

def make_mars():
    return _body_from_elements(
        'Mars', 6.4171e23, 3_389.5e3,
        a_au=1.523662, e=0.093412, i_deg=1.85061,
        raan_deg=49.57854, argp_deg=286.4623, nu_deg=19.41248)

def make_jupiter():
    return _body_from_elements(
        'Jupiter', 1.898e27, 71_492e3,
        a_au=5.203363, e=0.048393, i_deg=1.3053,
        raan_deg=100.55615, argp_deg=-85.8023, nu_deg=19.55053)

def make_saturn():
    return _body_from_elements(
        'Saturn', 5.683e26, 58_232e3,
        a_au=9.537070, e=0.054151, i_deg=2.48446,
        raan_deg=113.71504, argp_deg=-21.2831, nu_deg=-42.4876)

def make_uranus():
    return _body_from_elements(
        'Uranus', 8.681e25, 25_362e3,
        a_au=19.19126, e=0.047168, i_deg=0.76986,
        raan_deg=74.22988, argp_deg=96.73436, nu_deg=142.2679)

def make_neptune():
    return _body_from_elements(
        'Neptune', 1.024e26, 24_622e3,
        a_au=30.06896, e=0.008586, i_deg=1.76917,
        raan_deg=131.72169, argp_deg=-86.75034, nu_deg=259.9087)

def make_pluto():
    return _body_from_elements(
        'Pluto', 1.309e22, 1_188.3e3,
        a_au=39.48169, e=0.248808, i_deg=17.14175,
        raan_deg=110.30347, argp_deg=113.76329, nu_deg=14.86205)


# ---------------------------------------------------------------------------
# Convenience bundles
# ---------------------------------------------------------------------------

def solar_system():
    """Return (sun, mercury, venus, earth, mars, jupiter, saturn) at J2000."""
    return (make_sun(), make_mercury(), make_venus(),
            make_earth(), make_mars(), make_jupiter(), make_saturn())

def all_planets():
    """Return all 8 planets + Pluto at J2000, as a dict keyed by name."""
    return {
        'Mercury': make_mercury(),
        'Venus':   make_venus(),
        'Earth':   make_earth(),
        'Mars':    make_mars(),
        'Jupiter': make_jupiter(),
        'Saturn':  make_saturn(),
        'Uranus':  make_uranus(),
        'Neptune': make_neptune(),
        'Pluto':   make_pluto(),
    }


# ---------------------------------------------------------------------------
# Minimum safe periapsis distances for flyby/capture
# ---------------------------------------------------------------------------

FLYBY_PERIAPSIS_MIN = {
    'Mercury': make_mercury().radius + 200e3,
    'Venus':   make_venus().radius   + 300e3,
    'Earth':   make_earth().radius   + 300e3,
    'Mars':    make_mars().radius    + 200e3,
    'Jupiter': make_jupiter().radius + 500e3,
    'Saturn':  make_saturn().radius  + 500e3,
}
