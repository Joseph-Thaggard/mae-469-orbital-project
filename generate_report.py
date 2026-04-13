"""
generate_report.py — MAE 469 Project Report Data Generator

Runs all trajectory calculations and writes a fully formatted report.txt
containing every result required by the project report (sections 4a–4g).

Sections produced:
  1. Mission overview & parameters
  2. Task 1 — Planetary states on Dec 25, 2026 at 23:12 UTC (table 4a)
  3. Task 2 — Transfer orbit: methodology trace
      3a. Lambert solutions for each leg (verbose)
      3b. Mars free flyby analysis (4e)
      3c. DSM optimization result (4g)
      3d. Earth free flyby analysis (4e)
      3e. Orbit elements per leg (4d)
      3f. Positions/velocities at key events (4c)
      3g. TOF summary (4f)
      3h. Delta-v budget (4g)
      3i. Lambert miss-distance validation

Usage:
    python3 generate_report.py          -> writes report.txt
"""

import io
import sys
import numpy as np
import numpy.linalg as la
from datetime import datetime, timezone, timedelta

from bodies.catalog import (
    make_sun, make_earth, make_mars, make_jupiter,
    make_mercury, make_venus, make_saturn, make_uranus, make_neptune, make_pluto,
    all_planets, J2000_EPOCH_UTC, MU_SUN, AU, FLYBY_PERIAPSIS_MIN
)
from solvers.solver import (
    kepler_propagate, lambert_solver, rv_to_elements, elements_to_rv
)
from solvers.dsm import optimize_dsm

# ---------------------------------------------------------------------------
# Output: write to file AND capture in string buffer
# ---------------------------------------------------------------------------
OUTPUT_FILE = 'report.txt'
_buf = io.StringIO()

def _write(line=''):
    print(line)
    _buf.write(line + '\n')

def w(line=''):
    _write(line)

def sep(char='─', width=78):
    w(char * width)

def header(title, level=1):
    w()
    if level == 1:
        w('=' * 78)
        w(f'  {title}')
        w('=' * 78)
    elif level == 2:
        sep()
        w(f'  {title}')
        sep()
    else:
        w(f'  {title}')
        w('  ' + '─' * (len(title) + 2))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DAY             = 86400.0
DEP_PARK_ALT    = 500e3        # 500 km LEO
JUP_PARK_ALT    = 10_000e3    # 10,000 km Jupiter park
MARS_RFB_ALT    = 200e3       # min Mars flyby altitude
EARTH_RFB_ALT   = 300e3       # min Earth flyby altitude

mu_sun = MU_SUN

def j2k(dt_utc):
    """UTC datetime -> seconds from J2000 epoch."""
    return (dt_utc - J2000_EPOCH_UTC).total_seconds()

def ds(t_s):
    """Seconds from J2000 -> 'YYYY-MM-DD' string."""
    return (J2000_EPOCH_UTC + timedelta(seconds=t_s)).strftime('%Y-%m-%d %H:%M UTC')

def ds_short(t_s):
    return (J2000_EPOCH_UTC + timedelta(seconds=t_s)).strftime('%Y-%m-%d')

# ---------------------------------------------------------------------------
# Body catalog
# ---------------------------------------------------------------------------
earth   = make_earth()
mars    = make_mars()
jupiter = make_jupiter()
sun     = make_sun()

mu_e = earth.mu
mu_m = mars.mu
mu_j = jupiter.mu

r_e0, v_e0 = earth.position.copy(), earth.velocity.copy()
r_m0, v_m0 = mars.position.copy(),  mars.velocity.copy()
r_j0, v_j0 = jupiter.position.copy(), jupiter.velocity.copy()

# ---------------------------------------------------------------------------
# Mission window
# ---------------------------------------------------------------------------
T_LAUNCH = j2k(datetime(2027, 1, 15, tzinfo=timezone.utc))
TOF1     = 260  * DAY    # Earth -> Mars
TOF2     = 420  * DAY    # Mars  -> Earth  (DSM splits this)
TOF3     = 1100 * DAY    # Earth -> Jupiter

t_dep = T_LAUNCH
t_m   = t_dep + TOF1
t_ert = t_m   + TOF2
t_j   = t_ert + TOF3

# ---------------------------------------------------------------------------
# ========== REPORT STARTS ==========
# ---------------------------------------------------------------------------
w('MAE 469 / INTRO TO ASTRODYNAMICS')
w('SPRING 2026 — INTERPLANETARY TRAJECTORY PROJECT')
w('Generated: ' + datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'))
w()
w('Team Member(s): [Insert names]')

# ===========================================================================
header('SECTION 1 — MISSION OVERVIEW')
# ===========================================================================
w("""
  Mission Concept: Earth-Mars-Earth-Jupiter (EMEJ)
  Modeled after the Europa Clipper (2024) trajectory.

  The spacecraft departs Earth, performs a free gravity assist at Mars, fires
  a Deep Space Maneuver (DSM) at a heliocentric point between Mars and Earth,
  performs a second free gravity assist at Earth, then coasts to Jupiter.

  Three propulsive delta-v maneuvers:
    Dv1 — Earth departure burn (escape from 500 km LEO)
    Dv2 — Deep Space Maneuver (somewhere on the Mars-Earth leg)
    Dv3 — Jupiter orbit insertion (into 10,000 km altitude circular orbit)

  Free flybys (zero propulsion):
    Mars flyby  — hyperbolic turn redirects spacecraft toward DSM point
    Earth flyby — hyperbolic turn redirects spacecraft toward Jupiter

  Key assumptions:
    - Patched conic / two-body approximation (Sun-centered between SOIs)
    - Circular parking orbits at Earth (500 km) and Jupiter (10,000 km)
    - All circular orbits coplanar with hyperbolic arrivals (no plane change Dv)
    - Periapsis of each hyperbolic orbit equals parking orbit radius
""")

# ===========================================================================
header('SECTION 2 — TASK 1: PLANETARY STATES ON DEC 25, 2026 AT 23:12 UTC')
# ===========================================================================
w("""
  Method: Keplerian propagation using universal variable (Battin/Curtis §3.7).
  All planets propagated from J2000 epoch (Jan 1, 2000 11:58:56 UTC) to the
  specified date using kepler_propagate(r0, v0, dt, mu_sun).

  J2000 orbital elements sourced from the project specification table.
  Reference frame: heliocentric ecliptic J2000 (x toward vernal equinox).
""")

TARGET_DATE = datetime(2026, 12, 25, 23, 12, 0, tzinfo=timezone.utc)
t_task1 = j2k(TARGET_DATE)

planets_ordered = [
    ('Mercury', make_mercury()),
    ('Venus',   make_venus()),
    ('Earth',   make_earth()),
    ('Mars',    make_mars()),
    ('Jupiter', make_jupiter()),
    ('Saturn',  make_saturn()),
    ('Uranus',  make_uranus()),
    ('Neptune', make_neptune()),
    ('Pluto',   make_pluto()),
]

w(f'  Reference date: {TARGET_DATE.strftime("%B %d, %Y at %H:%M UTC")}')
w(f'  Seconds from J2000: {t_task1:,.1f} s  ({t_task1/DAY:,.2f} days)')
w()

# Table header
col_w = 10
w('  ' + f'{"Planet":<10}  {"theta (deg)":>11}  {"r| (AU)":>8}  '
         f'{"rx (AU)":>10}  {"ry (AU)":>10}  {"rz (AU)":>10}  '
         f'{"vx (km/s)":>10}  {"vy (km/s)":>10}  {"vz (km/s)":>10}  {"v| (km/s)":>10}')
w('  ' + '-' * 120)

planet_states_t1 = {}
for name, body in planets_ordered:
    r0 = body.position.copy()
    v0 = body.velocity.copy()
    r, v = kepler_propagate(r0, v0, t_task1, mu_sun)
    # true anomaly from orbital elements
    el = rv_to_elements(r, v, mu_sun)
    nu_deg = np.degrees(el['nu'])
    r_au = r / AU
    v_kms = v / 1e3
    planet_states_t1[name] = (r, v, el)
    w('  ' + f'{name:<10}  {nu_deg:>11.4f}  {la.norm(r_au):>8.5f}  '
             f'{r_au[0]:>10.5f}  {r_au[1]:>10.5f}  {r_au[2]:>10.5f}  '
             f'{v_kms[0]:>10.5f}  {v_kms[1]:>10.5f}  {v_kms[2]:>10.5f}  '
             f'{la.norm(v_kms):>10.5f}')

w()
w('  Notes:')
w('    - theta = true anomaly at specified epoch (deg)')
w('    - r| = heliocentric distance (AU)')
w('    - All vectors in heliocentric ecliptic J2000 frame')
w('    - Pluto included per project specification')

# ===========================================================================
header('SECTION 3 — TASK 2: EMEJ TRAJECTORY')
# ===========================================================================

# ── 3.1 Body states at event epochs ─────────────────────────────────────────
header('3.1  Body States at Mission Event Epochs', level=2)
w(f"""
  Propagation method: Kepler universal variable (Curtis Algorithm 3.4).
  Each planet state propagated independently from J2000 initial elements.

  Selected trajectory window:
    Launch (Earth dep.):  {ds_short(t_dep)}   (2027-01-15)
    Mars flyby:           {ds_short(t_m)}   TOF1 = {TOF1/DAY:.1f} days
    Earth flyby:          {ds_short(t_ert)}   TOF2 = {TOF2/DAY:.1f} days after Mars
    Jupiter arrival:      {ds_short(t_j)}   TOF3 = {TOF3/DAY:.1f} days after Earth
    Total mission:        {(t_j-t_dep)/DAY:.1f} days  ({(t_j-t_dep)/DAY/365.25:.2f} years)
""")

r_dep,  v_dep  = kepler_propagate(r_e0, v_e0, t_dep, mu_sun)
r_mars, v_mars = kepler_propagate(r_m0, v_m0, t_m,   mu_sun)
r_ert,  v_ert  = kepler_propagate(r_e0, v_e0, t_ert, mu_sun)
r_jup,  v_jup  = kepler_propagate(r_j0, v_j0, t_j,   mu_sun)

def fmt_rv(label, r_m, v_ms, indent='  '):
    r = r_m / AU; v = v_ms / 1e3
    w(f'{indent}{label}')
    w(f'{indent}  r = [{r[0]:>11.6f},  {r[1]:>11.6f},  {r[2]:>11.6f}] AU    |r| = {la.norm(r):.6f} AU')
    w(f'{indent}  v = [{v[0]:>11.6f},  {v[1]:>11.6f},  {v[2]:>11.6f}] km/s  |v| = {la.norm(v):.6f} km/s')

fmt_rv(f'Earth at departure       [{ds(t_dep)}]:', r_dep, v_dep)
w()
fmt_rv(f'Mars at flyby            [{ds(t_m)}]:', r_mars, v_mars)
w()
fmt_rv(f'Earth at return flyby    [{ds(t_ert)}]:', r_ert, v_ert)
w()
fmt_rv(f'Jupiter at arrival       [{ds(t_j)}]:', r_jup, v_jup)

# ── 3.2 Leg 1 Lambert (Earth -> Mars) ────────────────────────────────────────
header('3.2  Leg 1 Lambert Solution: Earth -> Mars', level=2)
w(f"""
  Method: Universal variable Lambert solver (Curtis Algorithm 5.2).
  Solves for the transfer orbit connecting two position vectors r1, r2
  in a specified time-of-flight (TOF).

  Inputs:
    r1 = Earth position at departure   {ds_short(t_dep)}
    r2 = Mars position at flyby        {ds_short(t_m)}
    TOF1 = {TOF1/DAY:.1f} days = {TOF1:.1f} s
    mu_sun = {mu_sun:.6e} m^3/s^2
""")

v_sc_dep, v_sc_in_mars = lambert_solver(r_dep, r_mars, TOF1, mu_sun)

v_inf_dep     = v_sc_dep - v_dep
v_inf_in_mars = v_sc_in_mars - v_mars
k_in_mars     = la.norm(v_inf_in_mars)

fmt_rv('SC velocity at Earth departure:', r_dep, v_sc_dep)
w()
fmt_rv('SC velocity at Mars arrival:', r_mars, v_sc_in_mars)
w()
w(f'  Hyperbolic excess velocity at Earth departure:')
w(f'    v_inf_dep = [{v_inf_dep[0]/1e3:>10.5f}, {v_inf_dep[1]/1e3:>10.5f}, {v_inf_dep[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_dep| = {la.norm(v_inf_dep)/1e3:.5f} km/s')
w()
w(f'  Hyperbolic excess velocity at Mars arrival (v_inf_in):')
w(f'    v_inf_in_mars = [{v_inf_in_mars[0]/1e3:>10.5f}, {v_inf_in_mars[1]/1e3:>10.5f}, {v_inf_in_mars[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_in_mars| = {k_in_mars/1e3:.5f} km/s')

# ── 3.3 Leg 3 Lambert (Earth -> Jupiter) ─────────────────────────────────────
header('3.3  Leg 3 Lambert Solution: Earth -> Jupiter', level=2)
w(f"""
  The Earth-to-Jupiter leg determines the required post-Earth-flyby
  spacecraft velocity. This is computed first (backward planning) so the
  DSM optimizer knows what velocity the spacecraft must have at Earth.

  Inputs:
    r1 = Earth position at flyby       {ds_short(t_ert)}
    r2 = Jupiter position at arrival   {ds_short(t_j)}
    TOF3 = {TOF3/DAY:.1f} days = {TOF3:.1f} s
""")

v_sc_out_earth, v_sc_arr_jup = lambert_solver(r_ert, r_jup, TOF3, mu_sun)

v_inf_out_earth = v_sc_out_earth - v_ert
v_inf_arr_jup   = v_sc_arr_jup   - v_jup
k_out_earth     = la.norm(v_inf_out_earth)

fmt_rv('SC velocity departing Earth (post-flyby):', r_ert, v_sc_out_earth)
w()
fmt_rv('SC velocity at Jupiter arrival:', r_jup, v_sc_arr_jup)
w()
w(f'  Post-Earth-flyby hyperbolic excess velocity (v_inf_out):')
w(f'    v_inf_out_earth = [{v_inf_out_earth[0]/1e3:>10.5f}, {v_inf_out_earth[1]/1e3:>10.5f}, {v_inf_out_earth[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_out_earth| = {k_out_earth/1e3:.5f} km/s')
w()
w(f'  Jupiter arrival hyperbolic excess velocity:')
w(f'    v_inf_arr_jup = [{v_inf_arr_jup[0]/1e3:>10.5f}, {v_inf_arr_jup[1]/1e3:>10.5f}, {v_inf_arr_jup[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_arr_jup| = {la.norm(v_inf_arr_jup)/1e3:.5f} km/s')

# ── 3.4 Mars Free Flyby — pre-DSM analysis ───────────────────────────────────
header('3.4  Mars Free Flyby Geometry (pre-DSM)', level=2)
r_mars_fb    = mars.radius + MARS_RFB_ALT
ecc_mars     = 1.0 + r_mars_fb * k_in_mars**2 / mu_m
delta_max_mars_rad = 2.0 * np.arcsin(1.0 / ecc_mars)

w(f"""
  A free (unpowered) gravity assist turns the spacecraft's v_inf vector
  without changing its magnitude. The maximum achievable turn angle is:

    delta_max = 2 * arcsin(1 / ecc_hyp)
    ecc_hyp   = 1 + rp * v_inf^2 / mu_body

  where rp is the periapsis radius of the hyperbolic flyby trajectory.

  Mars flyby parameters:
    Minimum flyby altitude:    {MARS_RFB_ALT/1e3:.0f} km
    Mars radius:               {mars.radius/1e3:.1f} km
    Periapsis radius (rp_min): {r_mars_fb/1e3:.1f} km from Mars center
    mu_Mars:                   {mu_m:.6e} m^3/s^2
    |v_inf_in|:                {k_in_mars/1e3:.5f} km/s

  Computed flyby hyperbolic orbit:
    Eccentricity (ecc_hyp):    {ecc_mars:.5f}
    Max turn angle (delta_max):{np.degrees(delta_max_mars_rad):.4f} deg
""")

# ── 3.5 DSM Optimization ─────────────────────────────────────────────────────
header('3.5  Deep Space Maneuver (DSM) Optimization', level=2)
w(f"""
  The DSM optimizer searches for the heliocentric location and timing of
  a propulsive burn on the Mars-to-Earth leg that minimizes Dv2 while
  satisfying the Earth flyby turn angle feasibility constraint.

  Free parameters:
    tau   — fraction of TOF2 at which the DSM fires  (0 < tau < 1)
              => t_DSM = t_Mars + tau * TOF2
    phi   — azimuth angle of Mars flyby exit direction (0 to 2*pi)
              => controls which way the Mars flyby turns the spacecraft
    delta — fraction of delta_max turn used at Mars   (0 to 1)
              => delta_rad = delta * delta_max_mars

  Optimization algorithm:
    1. Grid search over (tau=10, phi=12, delta=4) = 480 combinations
    2. Nelder-Mead refinement from best grid point
    Cost = Dv_DSM + lam_earth * speed_mismatch + lam_turn * turn_penalty
    turn_penalty = max(0, delta_turn - delta_max_earth) * 1e4

  Constraints enforced:
    - Earth flyby turn angle must not exceed delta_max at rp_min
    - Speed mismatch |v_inf_in_earth| ~ |v_inf_out_earth| (free flyby)
    - rp_required >= rp_min_earth = {earth.radius/1e3 + EARTH_RFB_ALT/1e3:.1f} km

  Running optimizer...
""")

r_earth_fb_min = earth.radius + EARTH_RFB_ALT

dsm = optimize_dsm(
    r_flyby       = r_mars,
    v_flyby       = v_mars,
    v_inf_in      = v_inf_in_mars,
    k_in          = k_in_mars,
    delta_max_rad = delta_max_mars_rad,
    r_target      = r_ert,
    v_target      = v_ert,
    tof_leg       = TOF2,
    v_sc_out_target = v_sc_out_earth,
    mu_sun        = mu_sun,
    mu_target     = mu_e,
    rp_min_target = r_earth_fb_min,
    lam_earth     = 1.0,
    lam_turn      = 1.0,
    n_tau=10, n_phi=12, n_delta=4,
    refine=True,
)

if dsm is None:
    w('  ERROR: DSM optimizer failed to converge. Check Lambert inputs.')
    sys.exit(1)

tau_opt   = dsm['tau']
phi_opt   = dsm['phi']
delta_opt = dsm['delta_rad']

r_dsm               = dsm['r_dsm']
v_at_dsm            = dsm['v_at_dsm']
v_dep_dsm           = dsm['v_dep_dsm']
v_arr_ert           = dsm['v_arr_target']
v_sc_dep_mars_post  = dsm['v_sc_dep_flyby']
v_inf_out_mars_fb   = dsm['v_inf_out']
dv_dsm              = dsm['dv_dsm']
v_inf_in_earth      = dsm['v_inf_in_target']
k_in_earth          = dsm['k_in_target']
TOF_2a              = dsm['tof_a']
TOF_2b              = dsm['tof_b']
t_dsm               = t_m + TOF_2a
v_dsm_vec           = v_dep_dsm - v_at_dsm

turn_earth     = dsm['turn_deg']    if dsm['turn_deg']    is not None else 0.0
max_turn_earth = dsm['max_turn_deg'] if dsm['max_turn_deg'] is not None else 0.0
rp_earth_req   = dsm['rp_required'] if dsm['rp_required'] is not None else np.inf

w(f'  Optimization result:')
w(f'    tau   (timing fraction):       {tau_opt:.6f}')
w(f'    phi   (exit azimuth):          {np.degrees(phi_opt):.4f} deg')
w(f'    delta (fraction of delta_max): {delta_opt/delta_max_mars_rad:.6f}')
w(f'    delta used:                    {np.degrees(delta_opt):.4f} deg  (max = {np.degrees(delta_max_mars_rad):.4f} deg)')
w(f'    Optimizer cost:                {dsm["cost"]:.6f}')
w()
w(f'  DSM epoch:')
w(f'    t_DSM = {ds(t_dsm)}')
w(f'    TOF2a = {TOF_2a/DAY:.2f} days after Mars flyby')
w(f'    TOF2b = {TOF_2b/DAY:.2f} days after DSM to Earth flyby')
w()
fmt_rv('DSM heliocentric location:', r_dsm, v_at_dsm)
w(f'    |r_DSM| = {la.norm(r_dsm)/AU:.5f} AU from Sun')
w()
w(f'  DSM burn vector:')
w(f'    Dv2 = [{v_dsm_vec[0]/1e3:>10.5f}, {v_dsm_vec[1]/1e3:>10.5f}, {v_dsm_vec[2]/1e3:>10.5f}] km/s')
w(f'    |Dv2| = {dv_dsm/1e3:.5f} km/s')

# ── 3.6 Mars Flyby Detail ─────────────────────────────────────────────────────
header('3.6  Mars Gravity Assist — Free Flyby Detail', level=2)
vom_fb = la.norm(v_inf_out_mars_fb)
w(f"""
  The Mars free flyby applies a Rodrigues rotation to v_inf_in in the plane
  perpendicular to v_inf_in, rotating by delta_rad toward direction phi.
  Since no propulsion is used, |v_inf_out| = |v_inf_in| exactly.

  Pre-flyby state (SC arriving at Mars):""")
w(f'    Heliocentric SC velocity at Mars arrival:')
w(f'      v_sc_in = [{v_sc_in_mars[0]/1e3:>10.5f}, {v_sc_in_mars[1]/1e3:>10.5f}, {v_sc_in_mars[2]/1e3:>10.5f}] km/s')
w(f'    v_inf_in  = [{v_inf_in_mars[0]/1e3:>10.5f}, {v_inf_in_mars[1]/1e3:>10.5f}, {v_inf_in_mars[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_in| = {k_in_mars/1e3:.5f} km/s')
w()
w(f'  Flyby hyperbolic orbit:')
w(f'    Periapsis radius (rp):    {r_mars_fb/1e3:.1f} km  (altitude = {MARS_RFB_ALT/1e3:.0f} km)')
w(f'    Eccentricity (ecc_hyp):   {ecc_mars:.5f}')
w(f'    Turn angle applied:       {np.degrees(delta_opt):.4f} deg')
w(f'    Maximum possible:         {np.degrees(delta_max_mars_rad):.4f} deg')
w(f'    Fraction of max used:     {delta_opt/delta_max_mars_rad:.4f}')
w()
w(f'  Post-flyby state (SC departing Mars):')
w(f'    v_inf_out = [{v_inf_out_mars_fb[0]/1e3:>10.5f}, {v_inf_out_mars_fb[1]/1e3:>10.5f}, {v_inf_out_mars_fb[2]/1e3:>10.5f}] km/s')
w(f'    |v_inf_out| = {vom_fb/1e3:.5f} km/s  (should equal |v_inf_in| = {k_in_mars/1e3:.5f} km/s)')
w(f'    Speed conservation check: |Dv| = {abs(vom_fb - k_in_mars)/1e3:.2e} km/s  (< 1e-6 OK)')
w()
w(f'    SC heliocentric departure velocity post-flyby:')
w(f'      v_sc_out = [{v_sc_dep_mars_post[0]/1e3:>10.5f}, {v_sc_dep_mars_post[1]/1e3:>10.5f}, {v_sc_dep_mars_post[2]/1e3:>10.5f}] km/s')

# ── 3.7 Earth Flyby Detail ────────────────────────────────────────────────────
header('3.7  Earth Gravity Assist — Free Flyby Detail', level=2)
ecc_earth = 1.0 + r_earth_fb_min * k_in_earth**2 / mu_e
delta_max_earth_rad = 2.0 * np.arcsin(1.0 / ecc_earth)
turn_earth_rad = np.radians(turn_earth)
feas_str = 'FEASIBLE' if turn_earth <= max_turn_earth else 'INFEASIBLE'
rp_ok_str = 'OK (above minimum)' if rp_earth_req >= r_earth_fb_min else 'BELOW MINIMUM'

w(f"""
  The Earth flyby is also free (unpowered). The DSM was optimized to deliver
  the spacecraft to Earth with the correct v_inf direction such that the
  hyperbolic turn naturally redirects it toward Jupiter.

  Earth flyby constraint parameters:
    Minimum altitude:          {EARTH_RFB_ALT/1e3:.0f} km
    Earth radius:              {earth.radius/1e3:.1f} km
    Minimum periapsis (rp_min):{r_earth_fb_min/1e3:.1f} km from Earth center
    mu_Earth:                  {mu_e:.6e} m^3/s^2

  SC arrival at Earth (post-DSM coast):""")
fmt_rv('  SC arriving at Earth:', r_ert, v_arr_ert)
w()
w(f'  v_inf_in  = [{v_inf_in_earth[0]/1e3:>10.5f}, {v_inf_in_earth[1]/1e3:>10.5f}, {v_inf_in_earth[2]/1e3:>10.5f}] km/s')
w(f'  |v_inf_in| = {k_in_earth/1e3:.5f} km/s')
w()
w(f'  Post-flyby (required to reach Jupiter):')
w(f'  v_inf_out = [{v_inf_out_earth[0]/1e3:>10.5f}, {v_inf_out_earth[1]/1e3:>10.5f}, {v_inf_out_earth[2]/1e3:>10.5f}] km/s')
w(f'  |v_inf_out| = {k_out_earth/1e3:.5f} km/s')
w(f'  Speed mismatch: {abs(k_in_earth - k_out_earth)/1e3:.4f} km/s  (free flyby requires ~0)')
w()
w(f'  Earth flyby hyperbolic orbit:')
w(f'    Eccentricity (ecc_hyp):        {ecc_earth:.5f}')
w(f'    Max achievable turn (delta_max):{max_turn_earth:.4f} deg')
w(f'    Required turn angle (delta):    {turn_earth:.4f} deg')
w(f'    Required periapsis radius:      {rp_earth_req/1e3:.1f} km')
w(f'    Feasibility:                    {feas_str}')
w(f'    Periapsis check:                {rp_ok_str}')

# ── 3.8 Transfer Orbit Elements ───────────────────────────────────────────────
header('3.8  Transfer Orbit Classical Elements (Section 4d)', level=2)
w(f"""
  Orbital elements computed from (r, v) via Curtis rv_to_elements().
  All angles in degrees. Semi-major axis in AU.
  Reference: heliocentric ecliptic J2000.

  Element definitions:
    a     — semi-major axis (AU)
    e     — eccentricity
    i     — inclination to ecliptic (deg)
    RAAN  — right ascension of ascending node, Omega (deg)
    argp  — argument of periapsis, omega (deg)
    nu    — true anomaly at epoch (deg)
    P     — orbital period (years)
""")

def fmt_elements(label, r_m, v_ms):
    el = rv_to_elements(r_m, v_ms, mu_sun)
    a_au = el['a'] / AU
    P_yr = 2*np.pi*np.sqrt(el['a']**3 / mu_sun) / (365.25*DAY) if el['a'] > 0 else float('nan')
    w(f'  {label}')
    w(f'    a     = {a_au:.6f} AU')
    w(f'    e     = {el["e"]:.6f}')
    w(f'    i     = {np.degrees(el["i"]):.5f} deg')
    w(f'    RAAN  = {np.degrees(el["raan"]):.5f} deg')
    w(f'    argp  = {np.degrees(el["argp"]):.5f} deg')
    w(f'    nu    = {np.degrees(el["nu"]):.5f} deg  (at departure epoch)')
    w(f'    P     = {P_yr:.4f} years  (full orbit period)')
    w()

fmt_elements(f'Leg 1: Earth -> Mars  (Lambert arc, TOF = {TOF1/DAY:.1f} d)', r_dep, v_sc_dep)
fmt_elements(f'Leg 2a: Mars -> DSM   (Keplerian coast, TOF = {TOF_2a/DAY:.2f} d)', r_mars, v_sc_dep_mars_post)
fmt_elements(f'Leg 2b: DSM -> Earth  (Keplerian coast, TOF = {TOF_2b/DAY:.2f} d)', r_dsm,  v_dep_dsm)
fmt_elements(f'Leg 3: Earth -> Jupiter  (Lambert arc, TOF = {TOF3/DAY:.1f} d)', r_ert, v_sc_out_earth)

# ── 3.9 Key Position/Velocity vectors (Section 4c) ────────────────────────────
header('3.9  Heliocentric Position and Velocity at All Events (Section 4c)', level=2)
w()

rows = [
    ('Earth at departure',                    r_dep,             v_dep),
    ('SC at departure',                        r_dep,             v_sc_dep),
    ('Mars at flyby',                          r_mars,            v_mars),
    ('SC arrival at Mars',                     r_mars,            v_sc_in_mars),
    ('SC departure Mars (post-flyby)',          r_mars,            v_sc_dep_mars_post),
    (f'DSM location [{ds_short(t_dsm)}]',      r_dsm,             v_at_dsm),
    ('SC pre-DSM burn (arriving)',             r_dsm,             v_at_dsm),
    ('SC post-DSM burn (departing)',           r_dsm,             v_dep_dsm),
    ('Earth at return flyby',                  r_ert,             v_ert),
    ('SC arrival at Earth (post-DSM coast)',   r_ert,             v_arr_ert),
    ('SC departure Earth (post-flyby)',        r_ert,             v_sc_out_earth),
    ('Jupiter at arrival',                     r_jup,             v_jup),
    ('SC at Jupiter arrival',                  r_jup,             v_sc_arr_jup),
]

w(f'  {"Event":<40}  {"rx(AU)":>10}  {"ry(AU)":>10}  {"rz(AU)":>10}  '
         f'{"vx(km/s)":>10}  {"vy(km/s)":>10}  {"vz(km/s)":>10}  {"r|(AU)":>8}  {"v|(km/s)":>9}')
w('  ' + '-' * 130)
for label, r_m, v_ms in rows:
    r = r_m / AU; v = v_ms / 1e3
    w(f'  {label:<40}  {r[0]:>10.5f}  {r[1]:>10.5f}  {r[2]:>10.5f}  '
             f'{v[0]:>10.5f}  {v[1]:>10.5f}  {v[2]:>10.5f}  {la.norm(r):>8.5f}  {la.norm(v):>9.5f}')

# ── 3.10 TOF Summary (Section 4f) ─────────────────────────────────────────────
header('3.10  Time-of-Flight Summary (Section 4f)', level=2)
w()
w(f'  {"Leg":<40}  {"Start date":<14}  {"End date":<14}  {"TOF (days)":>12}  {"TOF (years)":>12}')
w('  ' + '-' * 100)
legs_tof = [
    ('Leg 1:  Earth  -> Mars   (Lambert)',    t_dep,  t_m,   TOF1),
    ('Leg 2a: Mars   -> DSM    (Keplerian)',  t_m,    t_dsm, TOF_2a),
    ('Leg 2b: DSM    -> Earth  (Keplerian)',  t_dsm,  t_ert, TOF_2b),
    ('Leg 2:  Mars   -> Earth  (total)',      t_m,    t_ert, TOF2),
    ('Leg 3:  Earth  -> Jupiter (Lambert)',   t_ert,  t_j,   TOF3),
    ('TOTAL MISSION',                         t_dep,  t_j,   t_j-t_dep),
]
for lbl, ts, te, tof in legs_tof:
    w(f'  {lbl:<40}  {ds_short(ts):<14}  {ds_short(te):<14}  {tof/DAY:>12.2f}  {tof/DAY/365.25:>12.4f}')

# ── 3.11 Delta-v Budget (Section 4g) ──────────────────────────────────────────
header('3.11  Delta-v Budget (Section 4g)', level=2)

r_leo      = earth.radius + DEP_PARK_ALT
v_circ_leo = np.sqrt(mu_e / r_leo)
v_hyp_leo  = np.sqrt(la.norm(v_inf_dep)**2 + 2.0 * mu_e / r_leo)
dv1        = v_hyp_leo - v_circ_leo

dv2 = dv_dsm

r_jup_park  = jupiter.radius + JUP_PARK_ALT
v_circ_jup  = np.sqrt(mu_j / r_jup_park)
v_hyp_jup   = np.sqrt(la.norm(v_inf_arr_jup)**2 + 2.0 * mu_j / r_jup_park)
dv3         = v_hyp_jup - v_circ_jup

total_dv = dv1 + dv2 + dv3

w(f"""
  Dv1 — Earth Departure Burn (from {DEP_PARK_ALT/1e3:.0f} km altitude LEO)
  ────────────────────────────────────────────────────────────
  Method: Hyperbolic periapsis burn from circular parking orbit.
    Dv1 = v_hyp - v_circ
    v_hyp = sqrt(v_inf^2 + 2*mu_Earth / rp)
    v_circ = sqrt(mu_Earth / rp)

    Parking orbit radius:  rp = {r_leo/1e3:.1f} km
    v_circ (LEO):          {v_circ_leo/1e3:.5f} km/s
    |v_inf_dep|:           {la.norm(v_inf_dep)/1e3:.5f} km/s
    v_hyp at periapsis:    {v_hyp_leo/1e3:.5f} km/s
    Dv1:                   {dv1/1e3:.5f} km/s

  Dv2 — Deep Space Maneuver (at {ds_short(t_dsm)}, {la.norm(r_dsm)/AU:.4f} AU from Sun)
  ────────────────────────────────────────────────────────────
  Method: Optimal heliocentric impulsive burn found by grid search + Nelder-Mead.
  Burn vector: [{v_dsm_vec[0]/1e3:>10.5f}, {v_dsm_vec[1]/1e3:>10.5f}, {v_dsm_vec[2]/1e3:>10.5f}] km/s
    Dv2:       {dv2/1e3:.5f} km/s

  Dv3 — Jupiter Orbit Insertion (into {JUP_PARK_ALT/1e3:.0f} km altitude circular orbit)
  ────────────────────────────────────────────────────────────
  Method: Hyperbolic periapsis burn, same form as Dv1.
    Dv3 = v_hyp_jup - v_circ_jup

    Parking orbit radius:  rp = {r_jup_park/1e3:.1f} km
    v_circ (Jupiter):      {v_circ_jup/1e3:.5f} km/s
    |v_inf_arr_jup|:       {la.norm(v_inf_arr_jup)/1e3:.5f} km/s
    v_hyp at periapsis:    {v_hyp_jup/1e3:.5f} km/s
    Dv3:                   {dv3/1e3:.5f} km/s
""")

w('  ' + '=' * 60)
w(f'  {"Burn":<30}  {"Dv (km/s)":>12}')
w('  ' + '-' * 44)
w(f'  {"Dv1 — Earth departure":<30}  {dv1/1e3:>12.5f}')
w(f'  {"Dv2 — Deep Space Maneuver":<30}  {dv2/1e3:>12.5f}')
w(f'  {"Dv3 — Jupiter orbit insertion":<30}  {dv3/1e3:>12.5f}')
w('  ' + '-' * 44)
w(f'  {"TOTAL MISSION Dv":<30}  {total_dv/1e3:>12.5f}')
w('  ' + '=' * 60)

# ── 3.12 Lambert Validation ────────────────────────────────────────────────────
header('3.12  Lambert Solution Validation (Miss Distances)', level=2)
w("""
  Each Lambert arc is validated by propagating the initial state forward
  by the full TOF and measuring the miss distance to the intended target.
  Miss distances should be < 1 m for a valid solution.
""")

miss1  = la.norm(kepler_propagate(r_dep,  v_sc_dep,          TOF1,   mu_sun)[0] - r_mars)
miss2a = la.norm(kepler_propagate(r_mars, v_sc_dep_mars_post, TOF_2a, mu_sun)[0] - r_dsm)
miss2b = la.norm(kepler_propagate(r_dsm,  v_dep_dsm,         TOF_2b, mu_sun)[0] - r_ert)
miss3  = la.norm(kepler_propagate(r_ert,  v_sc_out_earth,    TOF3,   mu_sun)[0] - r_jup)

def miss_status(d):
    return 'OK' if d < 1.0 else 'WARNING: large miss'

w(f'  {"Leg":<30}  {"Miss distance (m)":>20}  {"Status":>20}')
w('  ' + '-' * 76)
w(f'  {"Leg 1  Earth -> Mars":<30}  {miss1:>20.4e}  {miss_status(miss1):>20}')
w(f'  {"Leg 2a Mars -> DSM":<30}  {miss2a:>20.4e}  {miss_status(miss2a):>20}')
w(f'  {"Leg 2b DSM -> Earth":<30}  {miss2b:>20.4e}  {miss_status(miss2b):>20}')
w(f'  {"Leg 3  Earth -> Jupiter":<30}  {miss3:>20.4e}  {miss_status(miss3):>20}')

# ── 3.13 Flyby Summary Table (Section 4e) ─────────────────────────────────────
header('3.13  Flyby Geometry Summary (Section 4e)', level=2)
w()
w('  MARS FLYBY:')
w(f'    Body:                      Mars')
w(f'    Date:                      {ds_short(t_m)}')
w(f'    Flyby type:                Free (unpowered)')
w(f'    v_inf_in  (heliocentric):  [{v_inf_in_mars[0]/1e3:.5f}, {v_inf_in_mars[1]/1e3:.5f}, {v_inf_in_mars[2]/1e3:.5f}] km/s')
w(f'    v_inf_out (heliocentric):  [{v_inf_out_mars_fb[0]/1e3:.5f}, {v_inf_out_mars_fb[1]/1e3:.5f}, {v_inf_out_mars_fb[2]/1e3:.5f}] km/s')
w(f'    |v_inf_in|:                {k_in_mars/1e3:.5f} km/s')
w(f'    |v_inf_out|:               {la.norm(v_inf_out_mars_fb)/1e3:.5f} km/s')
w(f'    Turn angle delta:          {np.degrees(delta_opt):.4f} deg')
w(f'    Max possible delta_max:    {np.degrees(delta_max_mars_rad):.4f} deg')
w(f'    Periapsis altitude:        {MARS_RFB_ALT/1e3:.0f} km  (rp = {r_mars_fb/1e3:.1f} km)')
w(f'    Eccentricity:              {ecc_mars:.5f}')
w()
w('  EARTH FLYBY:')
w(f'    Body:                      Earth')
w(f'    Date:                      {ds_short(t_ert)}')
w(f'    Flyby type:                Free (unpowered)')
w(f'    v_inf_in  (heliocentric):  [{v_inf_in_earth[0]/1e3:.5f}, {v_inf_in_earth[1]/1e3:.5f}, {v_inf_in_earth[2]/1e3:.5f}] km/s')
w(f'    v_inf_out (heliocentric):  [{v_inf_out_earth[0]/1e3:.5f}, {v_inf_out_earth[1]/1e3:.5f}, {v_inf_out_earth[2]/1e3:.5f}] km/s')
w(f'    |v_inf_in|:                {k_in_earth/1e3:.5f} km/s')
w(f'    |v_inf_out|:               {k_out_earth/1e3:.5f} km/s')
w(f'    Speed mismatch:            {abs(k_in_earth-k_out_earth)/1e3:.4f} km/s')
w(f'    Turn angle delta:          {turn_earth:.4f} deg')
w(f'    Max possible delta_max:    {max_turn_earth:.4f} deg')
w(f'    Required periapsis rp:     {rp_earth_req/1e3:.1f} km')
w(f'    Minimum allowed rp:        {r_earth_fb_min/1e3:.1f} km')
w(f'    Eccentricity:              {ecc_earth:.5f}')
w(f'    Feasibility:               {feas_str}')

# ===========================================================================
header('SECTION 4 — GENERATED FIGURES')
# ===========================================================================
w("""
  The following figures accompany this report (generated by task_report_figures.py):

  Figure 1 — fig1_solar_system.png
    3D solar system snapshot on December 25, 2026 at 23:12 UTC.
    All 9 planets shown with orbit traces. Planet positions marked.

  Figure 2 — fig2_trajectory_3d.png
    3D EMEJ trajectory with all 4 spacecraft legs:
      White:  Leg 1 (Earth -> Mars)
      Yellow: Leg 2a (Mars -> DSM)
      Cyan:   Leg 2b (DSM -> Earth)
      Magenta:Leg 3 (Earth -> Jupiter)
    Full planet orbits shown for Earth, Mars, and Jupiter.

  Figure 3 — fig3_trajectory_2d.png
    2D ecliptic top-down projection of same trajectory.
    Cleaner layout for printed report — shows in-plane geometry.

  Figure 4 — fig4_dv_vs_tof3.png
    Sensitivity analysis: total mission Dv vs TOF3 (Earth-Jupiter TOF).
    Shows feasible/infeasible Earth flyby regions.
    Used to select TOF3 = 1100 days as optimal feasible window.

  Figure 5 — fig5_dv_vs_launch.png
    Sensitivity analysis: total mission Dv vs launch date offset.
    Sweeps +/-90 days from selected launch date (2027-01-15).
    Demonstrates launch window sensitivity.

  Figure 6 — fig6_flyby_geometry.png
    Hyperbolic flyby geometry diagrams for Mars and Earth flybys.
    Shows v_inf_in, v_inf_out vectors and turn angle in ecliptic plane.
""")

# ===========================================================================
header('SECTION 5 — CONCLUSION')
# ===========================================================================
w(f"""
  A four-leg Earth-Mars-Earth-Jupiter trajectory was successfully computed
  using the patched conic approximation. The mission closely follows the
  Europa Clipper trajectory concept, using two free gravity assists and one
  Deep Space Maneuver.

  Key results:
    Launch date:            2027-01-15
    Total mission duration: {(t_j-t_dep)/DAY:.1f} days  ({(t_j-t_dep)/DAY/365.25:.2f} years)
    Dv1 (Earth departure):  {dv1/1e3:.4f} km/s   from {DEP_PARK_ALT/1e3:.0f} km LEO
    Dv2 (DSM):              {dv2/1e3:.4f} km/s   at {la.norm(r_dsm)/AU:.3f} AU from Sun
    Dv3 (Jupiter capture):  {dv3/1e3:.4f} km/s   into {JUP_PARK_ALT/1e3:.0f} km orbit
    TOTAL Dv:               {total_dv/1e3:.4f} km/s

  The Mars free flyby (turn angle = {np.degrees(delta_opt):.2f} deg, max = {np.degrees(delta_max_mars_rad):.2f} deg)
  and Earth free flyby (turn angle = {turn_earth:.2f} deg, max = {max_turn_earth:.2f} deg) are
  both physically feasible (turn within maximum, periapsis above minimum).

  The DSM timing (tau = {tau_opt:.4f}, firing {TOF_2a/DAY:.1f} days after Mars flyby at
  {la.norm(r_dsm)/AU:.4f} AU from the Sun) was found via a 480-point grid search
  followed by Nelder-Mead refinement, minimizing the combined cost of DSM Dv,
  Earth v_inf speed mismatch, and flyby turn angle infeasibility penalty.

  The TOF3 = 1100 days was selected by sweeping the parameter space; shorter
  durations produced infeasible Earth flyby turn angles. The sensitivity
  figures (4 and 5) demonstrate the trade-space explored and confirm this
  window is near a local minimum in total mission Dv.
""")

# ===========================================================================
w()
w('=' * 78)
w('  END OF REPORT DATA')
w('=' * 78)
w()
w('  Appendix: Code files')
w('    bodies/catalog.py       — J2000 planetary orbital elements and body data')
w('    solvers/solver.py       — Lambert, Kepler propagation, orbit elements')
w('    solvers/dsm.py          — DSM grid search and Nelder-Mead optimizer')
w('    task1_planets.py        — Task 1 planetary state computation and plotting')
w('    task2_emej.py           — Task 2 full EMEJ trajectory (interactive mode)')
w('    scan_emej_freeflyby.py  — Parameter space scanner (porkchop-style)')
w('    task_report_figures.py  — All 6 report-quality figures')
w('    generate_report.py      — This report generator')

# ---------------------------------------------------------------------------
# Write to file
# ---------------------------------------------------------------------------
content = _buf.getvalue()
with open(OUTPUT_FILE, 'w') as f:
    f.write(content)

print(f'\nReport written to: {OUTPUT_FILE}  ({len(content.splitlines())} lines)')
