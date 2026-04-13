"""
task2_emej.py — Task 2: Earth-Mars-Earth-Jupiter (EMEJ) Trajectory

Mission profile (Europa Clipper style):
  Leg 1:  Earth departure → Mars free gravity assist
  Leg 2a: Mars (post-flyby) → DSM point  (Keplerian coast)
  DSM:    Deep space maneuver at a real heliocentric location
  Leg 2b: DSM point → Earth free gravity assist
  Leg 3:  Earth → Jupiter arrival & capture

3 total propulsive Δv:
  Δv1 — Earth departure (escape from 500 km LEO)
  Δv2 — Deep space maneuver (fired at heliocentric point between Mars and Earth)
  Δv3 — Jupiter arrival (captures into 10,000 km altitude circular orbit)

Selected base window (from scan_emej_freeflyby):
  Launch:        2027-01-15
  Mars flyby:    2027-10-02   TOF1 = 260 days
  Earth flyby:   2028-11-25   TOF2 = 420 days
  Jupiter arr:   2031-12-30   TOF3 = 1100 days
  Total mission: 1780 days (~4.87 years)

DSM placement is optimized over:
  tau   — fraction of TOF2 at which DSM fires (0 < tau < 1)
  phi   — exit direction of Mars free flyby (angle in perpendicular plane)
  delta — fraction of max turn used at Mars (0 = no turn, 1 = max δ_max)

Outputs:
  task2_emej_trajectory.png — 3D trajectory plot
  Console — full 4-leg trajectory report
"""

import numpy as np
import numpy.linalg as la
from datetime import datetime, timezone, timedelta

from bodies.catalog import (
    make_sun, make_earth, make_mars, make_jupiter,
    J2000_EPOCH_UTC, MU_SUN, AU, FLYBY_PERIAPSIS_MIN
)
from solvers.solver import (
    kepler_propagate, lambert_solver, rv_to_elements, elements_to_rv
)
from solvers.dsm import optimize_dsm

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DAY = 86400.0

DEP_PARK_ALT = 500e3        # 500 km — LEO departure altitude
JUP_PARK_ALT = 10_000e3     # 10,000 km — Jupiter capture altitude
MARS_RFB_ALT  = 200e3       # minimum Mars flyby altitude
EARTH_RFB_ALT = 300e3       # minimum Earth flyby altitude

# ---------------------------------------------------------------------------
# Body catalog (J2000 initial states)
# ---------------------------------------------------------------------------
earth   = make_earth()
mars    = make_mars()
jupiter = make_jupiter()

mu_sun = MU_SUN
mu_e   = earth.mu
mu_m   = mars.mu
mu_j   = jupiter.mu

r_e0, v_e0 = earth.position.copy(), earth.velocity.copy()
r_m0, v_m0 = mars.position.copy(), mars.velocity.copy()
r_j0, v_j0 = jupiter.position.copy(), jupiter.velocity.copy()

def j2k(dt_utc):
    return (dt_utc - J2000_EPOCH_UTC).total_seconds()

def ds(t_s):
    return (J2000_EPOCH_UTC + timedelta(seconds=t_s)).strftime('%Y-%m-%d')

# ---------------------------------------------------------------------------
# Selected base window
# ---------------------------------------------------------------------------
T_LAUNCH = j2k(datetime(2027, 1, 15, tzinfo=timezone.utc))
TOF1     = 260 * DAY    # Earth → Mars
TOF2     = 420 * DAY    # Mars  → Earth  (DSM splits this)
TOF3     = 1100 * DAY   # Earth → Jupiter  (adjusted for feasible Earth free flyby)

t_dep = T_LAUNCH
t_m   = t_dep + TOF1
t_ert = t_m   + TOF2
t_j   = t_ert + TOF3

# ---------------------------------------------------------------------------
# Body states at each event
# ---------------------------------------------------------------------------
r_dep,  v_dep  = kepler_propagate(r_e0, v_e0, t_dep, mu_sun)
r_mars, v_mars = kepler_propagate(r_m0, v_m0, t_m,   mu_sun)
r_ert,  v_ert  = kepler_propagate(r_e0, v_e0, t_ert, mu_sun)
r_jup,  v_jup  = kepler_propagate(r_j0, v_j0, t_j,   mu_sun)

# ---------------------------------------------------------------------------
# Leg 1: Earth → Mars  (Lambert)
# ---------------------------------------------------------------------------
v_sc_dep, v_sc_in_mars = lambert_solver(r_dep, r_mars, TOF1, mu_sun)

v_inf_dep      = v_sc_dep - v_dep
v_inf_in_mars  = v_sc_in_mars - v_mars
k_in_mars      = la.norm(v_inf_in_mars)

# Mars free flyby constraint
r_mars_fb    = mars.radius + MARS_RFB_ALT
ecc_mars     = 1.0 + r_mars_fb * k_in_mars**2 / mu_m
delta_max_mars_rad = 2.0 * np.arcsin(1.0 / ecc_mars)

# ---------------------------------------------------------------------------
# Leg 3: Earth → Jupiter  (Lambert — fixed by TOF3, determines Earth v∞_out)
# ---------------------------------------------------------------------------
v_sc_out_earth, v_sc_arr_jup = lambert_solver(r_ert, r_jup, TOF3, mu_sun)

v_inf_arr_jup  = v_sc_arr_jup - v_jup
v_inf_out_earth = v_sc_out_earth - v_ert
k_out_earth = la.norm(v_inf_out_earth)

# ---------------------------------------------------------------------------
# DSM optimization: find optimal (tau, phi, delta) over Mars→Earth leg
# ---------------------------------------------------------------------------
print("  Optimizing DSM placement (tau × phi × delta grid + Nelder-Mead) ...")

r_earth_fb_min = earth.radius + EARTH_RFB_ALT   # minimum rp at Earth flyby

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
    mu_target     = mu_e,                # Earth GM — for turn angle feasibility
    rp_min_target = r_earth_fb_min,      # 6671 km — enforce Earth flyby constraint
    lam_earth     = 1.0,                 # penalise speed mismatch
    lam_turn      = 1.0,                 # penalise infeasible turn angle
    n_tau=10, n_phi=12, n_delta=4,
    refine=True,
)

if dsm is None:
    raise RuntimeError("DSM optimizer failed — check Lambert convergence for this window")

tau_opt      = dsm['tau']
phi_opt      = dsm['phi']
delta_opt    = dsm['delta_rad']

r_dsm        = dsm['r_dsm']
v_at_dsm     = dsm['v_at_dsm']         # spacecraft velocity arriving at DSM (pre-burn)
v_dep_dsm    = dsm['v_dep_dsm']        # spacecraft velocity departing DSM (post-burn)
v_arr_ert    = dsm['v_arr_target']     # spacecraft velocity arriving at Earth

v_sc_dep_mars_post = dsm['v_sc_dep_flyby']  # post-flyby Mars departure velocity
v_inf_out_mars_fb  = dsm['v_inf_out']       # post-flyby v∞ at Mars

dv_dsm       = dsm['dv_dsm']
v_inf_in_earth = dsm['v_inf_in_target']
k_in_earth   = dsm['k_in_target']

# DSM split TOFs
TOF_2a = dsm['tof_a']    # Mars → DSM
TOF_2b = dsm['tof_b']    # DSM  → Earth

t_dsm = t_m + TOF_2a     # absolute DSM epoch (seconds from J2000)

# DSM Δv vector
v_dsm_vec = v_dep_dsm - v_at_dsm

# ---------------------------------------------------------------------------
# Earth free flyby analysis — use values pre-computed in DSM optimizer
# ---------------------------------------------------------------------------
r_earth_fb    = earth.radius + EARTH_RFB_ALT
turn_earth    = dsm['turn_deg']    if dsm['turn_deg']    is not None else 0.0
max_turn_earth= dsm['max_turn_deg'] if dsm['max_turn_deg'] is not None else 0.0
rp_earth      = dsm['rp_required'] if dsm['rp_required'] is not None else np.inf

# ---------------------------------------------------------------------------
# Lambert validation (miss distance checks)
# ---------------------------------------------------------------------------
miss1 = la.norm(kepler_propagate(r_dep,  v_sc_dep,         TOF1,  mu_sun)[0] - r_mars)
miss2a= la.norm(kepler_propagate(r_mars, v_sc_dep_mars_post, TOF_2a, mu_sun)[0] - r_dsm)
miss2b= la.norm(kepler_propagate(r_dsm,  v_dep_dsm,         TOF_2b, mu_sun)[0] - r_ert)
miss3 = la.norm(kepler_propagate(r_ert,  v_sc_out_earth,    TOF3,  mu_sun)[0] - r_jup)

# ---------------------------------------------------------------------------
# Δv budget
# ---------------------------------------------------------------------------
r_leo       = earth.radius + DEP_PARK_ALT
v_circ_leo  = np.sqrt(mu_e / r_leo)
v_hyp_leo   = np.sqrt(la.norm(v_inf_dep)**2 + 2.0 * mu_e / r_leo)
dv1         = v_hyp_leo - v_circ_leo

dv2         = dv_dsm

r_jup_park  = jupiter.radius + JUP_PARK_ALT
v_circ_jup  = np.sqrt(mu_j / r_jup_park)
v_hyp_jup   = np.sqrt(la.norm(v_inf_arr_jup)**2 + 2.0 * mu_j / r_jup_park)
dv3         = v_hyp_jup - v_circ_jup

total_dv    = dv1 + dv2 + dv3

# ---------------------------------------------------------------------------
# Print full report
# ---------------------------------------------------------------------------
print("=" * 75)
print("  Task 2: EMEJ Trajectory — Earth-Mars-Earth-Jupiter")
print("  (Europa Clipper Style, MAE 469)")
print("=" * 75)

print(f"""
  Mission Summary
  ───────────────────────────────────────────────────────────────
  Launch date:      {ds(t_dep)}
  Mars flyby:       {ds(t_m)}    TOF1  = {TOF1/DAY:.1f} days
  DSM:              {ds(t_dsm)}    TOF2a = {TOF_2a/DAY:.1f} days after Mars
  Earth flyby:      {ds(t_ert)}    TOF2b = {TOF_2b/DAY:.1f} days after DSM
  Jupiter arrival:  {ds(t_j)}    TOF3  = {TOF3/DAY:.1f} days
  Total mission:    {(t_j-t_dep)/DAY:.1f} days  ({(t_j-t_dep)/DAY/365.25:.2f} years)

  DSM parameters:
    tau   (fraction of Mars→Earth leg) = {tau_opt:.4f}
    phi   (flyby exit direction)        = {np.degrees(phi_opt):.2f}°
    delta (fraction of δ_max used)      = {delta_opt/delta_max_mars_rad:.4f}  ({np.degrees(delta_opt):.2f}° turn)

  Lambert validation:
    Leg 1  miss: {miss1:.2e} m
    Leg 2a miss: {miss2a:.2e} m
    Leg 2b miss: {miss2b:.2e} m
    Leg 3  miss: {miss3:.2e} m
""")

# --- Position/Velocity table ---
print("  Heliocentric Position and Velocity Vectors (AU and km/s):")
print(f"  {'Event':<30}  {'r_x':>10}  {'r_y':>10}  {'r_z':>10}"
      f"  {'v_x':>10}  {'v_y':>10}  {'v_z':>10}  {'|r|AU':>8}  {'|v|km/s':>9}")
print("  " + "-" * 115)

def prv(label, r_m, v_ms):
    r = r_m / AU; v = v_ms / 1e3
    print(f"  {label:<30}  {r[0]:>10.5f}  {r[1]:>10.5f}  {r[2]:>10.5f}"
          f"  {v[0]:>10.5f}  {v[1]:>10.5f}  {v[2]:>10.5f}"
          f"  {la.norm(r):>8.5f}  {la.norm(v):>9.5f}")

prv("Earth at departure",          r_dep,             v_dep)
prv("SC at departure",             r_dep,             v_sc_dep)
prv("Mars at flyby",               r_mars,            v_mars)
prv("SC arrival at Mars",          r_mars,            v_sc_in_mars)
prv("SC departure Mars (post-FB)", r_mars,            v_sc_dep_mars_post)
prv("SC arrival at DSM",           r_dsm,             v_at_dsm)
prv("SC departure DSM (post-burn)",r_dsm,             v_dep_dsm)
prv("Earth at return flyby",       r_ert,             v_ert)
prv("SC arrival at Earth",         r_ert,             v_arr_ert)
prv("SC departure Earth (post-FB)",r_ert,             v_sc_out_earth)
prv("Jupiter at arrival",          r_jup,             v_jup)
prv("SC at Jupiter arrival",       r_jup,             v_sc_arr_jup)

# --- Transfer orbit elements ---
print()
print("  Transfer Orbit Elements:")
print(f"  {'Leg':<30}  {'a (AU)':>8}  {'e':>8}  {'i (°)':>7}  "
      f"{'Ω (°)':>9}  {'ω (°)':>9}  {'θ₀ (°)':>9}")
print("  " + "-" * 88)

def pelem(label, r_m, v_ms):
    el = rv_to_elements(r_m, v_ms, mu_sun)
    print(f"  {label:<30}  {el['a']/AU:>8.4f}  {el['e']:>8.5f}  "
          f"{np.degrees(el['i']):>7.4f}  {np.degrees(el['raan']):>9.4f}  "
          f"{np.degrees(el['argp']):>9.4f}  {np.degrees(el['nu']):>9.4f}")

pelem("Leg 1  (Earth → Mars)",       r_dep,  v_sc_dep)
pelem("Leg 2a (Mars → DSM)",         r_mars, v_sc_dep_mars_post)
pelem("Leg 2b (DSM → Earth)",        r_dsm,  v_dep_dsm)
pelem("Leg 3  (Earth → Jupiter)",    r_ert,  v_sc_out_earth)

# --- Mars flyby ---
vom_fb = la.norm(v_inf_out_mars_fb)
print(f"""
  Mars Gravity Assist — Free Flyby (no propulsion at Mars):
  ─────────────────────────────────────────────────────────
    v∞_in   vector (km/s):  [{v_inf_in_mars[0]/1e3:>9.4f}, {v_inf_in_mars[1]/1e3:>9.4f}, {v_inf_in_mars[2]/1e3:>9.4f}]
    v∞_out  vector (km/s):  [{v_inf_out_mars_fb[0]/1e3:>9.4f}, {v_inf_out_mars_fb[1]/1e3:>9.4f}, {v_inf_out_mars_fb[2]/1e3:>9.4f}]
    |v∞_in|:                {k_in_mars/1e3:.4f} km/s
    |v∞_out| (free flyby):  {vom_fb/1e3:.4f} km/s  (= |v∞_in| ✓)
    Turn angle used δ:      {np.degrees(delta_opt):.2f}°   (max achievable: {np.degrees(delta_max_mars_rad):.2f}°)
    Periapsis altitude:     {MARS_RFB_ALT/1e3:.0f} km  ({r_mars_fb/1e3:.1f} km from Mars center)
""")

# --- DSM ---
print(f"  Deep Space Maneuver (Δv₂):")
print(f"  ──────────────────────────────────────────────────────")
print(f"    Location:          {ds(t_dsm)}, {TOF_2a/DAY:.1f} days after Mars flyby")
print(f"    r_DSM (AU):        [{r_dsm[0]/AU:>9.5f}, {r_dsm[1]/AU:>9.5f}, {r_dsm[2]/AU:>9.5f}]")
print(f"    |r_DSM|:           {la.norm(r_dsm)/AU:.4f} AU from Sun")
print(f"    Δv₂ vector (km/s): [{v_dsm_vec[0]/1e3:>9.4f}, {v_dsm_vec[1]/1e3:>9.4f}, {v_dsm_vec[2]/1e3:>9.4f}]")
print(f"    |Δv₂|:             {dv_dsm/1e3:.4f} km/s")

# --- Earth flyby ---
print(f"""
  Earth Gravity Assist — Free Flyby (no propulsion at Earth):
  ────────────────────────────────────────────────────────────
    v∞_in  vector (km/s):  [{v_inf_in_earth[0]/1e3:>9.4f}, {v_inf_in_earth[1]/1e3:>9.4f}, {v_inf_in_earth[2]/1e3:>9.4f}]
    v∞_out vector (km/s):  [{v_inf_out_earth[0]/1e3:>9.4f}, {v_inf_out_earth[1]/1e3:>9.4f}, {v_inf_out_earth[2]/1e3:>9.4f}]
    |v∞_in|:               {k_in_earth/1e3:.4f} km/s
    |v∞_out|:              {k_out_earth/1e3:.4f} km/s  (|Δ|={abs(k_in_earth-k_out_earth)/1e3:.4f} km/s)
    Turn angle δ:           {turn_earth:.2f}°
    Max achievable:         {max_turn_earth:.2f}°  ({'✓' if turn_earth <= max_turn_earth else '✗ INFEASIBLE'})
    Required periapsis rp:  {rp_earth/1e3:.1f} km  ({'✓ achievable' if rp_earth >= r_earth_fb else '✗ below minimum'})
    Min periapsis:          {r_earth_fb/1e3:.1f} km from Earth center
""")

# --- Δv budget ---
print(f"  Δv Budget:")
print(f"  ═══════════════════════════════════════════════════════")
print(f"  Δv₁ — Earth departure from {DEP_PARK_ALT/1e3:.0f} km altitude LEO:")
print(f"    v∞_dep:          {la.norm(v_inf_dep)/1e3:.4f} km/s")
print(f"    v_circ (LEO):    {v_circ_leo/1e3:.4f} km/s")
print(f"    v_hyp  (LEO):    {v_hyp_leo/1e3:.4f} km/s")
print(f"    Δv₁:             {dv1/1e3:.4f} km/s")
print()
print(f"  Δv₂ — Deep Space Maneuver (at {ds(t_dsm)}, {la.norm(r_dsm)/AU:.3f} AU from Sun):")
print(f"    Δv₂:             {dv2/1e3:.4f} km/s")
print()
print(f"  Δv₃ — Jupiter orbit insertion ({JUP_PARK_ALT/1e3:.0f} km altitude circular):")
print(f"    v∞_arr:          {la.norm(v_inf_arr_jup)/1e3:.4f} km/s")
print(f"    v_circ (Jup):    {v_circ_jup/1e3:.4f} km/s")
print(f"    v_hyp  (Jup):    {v_hyp_jup/1e3:.4f} km/s")
print(f"    Δv₃:             {dv3/1e3:.4f} km/s")
print()
print(f"  ───────────────────────────────────────────────────────")
print(f"  TOTAL Δv:          {total_dv/1e3:.4f} km/s")

# ---------------------------------------------------------------------------
# 3D Trajectory Plot
# ---------------------------------------------------------------------------
print()
print("  Generating 3D trajectory plot ...")

fig = plt.figure(figsize=(16, 12))
ax  = fig.add_subplot(111, projection='3d')

# Full planet orbits
for name, col, r0, v0 in [
    ('Earth',   '#4488ff', r_e0, v_e0),
    ('Mars',    '#cc4400', r_m0, v_m0),
    ('Jupiter', '#c8a064', r_j0, v_j0),
]:
    el = rv_to_elements(r0, v0, mu_sun)
    a, e, i, raan, argp = el['a'], el['e'], el['i'], el['raan'], el['argp']
    nu_vals = np.linspace(0, 2 * np.pi, 500)
    pts = np.array([elements_to_rv(a, e, i, raan, argp, nu, mu_sun)[0]
                    for nu in nu_vals]) / AU
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col, lw=0.8, alpha=0.45, label=f'{name} orbit')

# Spacecraft trajectory arcs
def arc(r0, v0, tof, mu, n=300):
    ts = np.linspace(0, tof, n)
    return np.array([kepler_propagate(r0, v0, t, mu)[0] / AU for t in ts])

leg1  = arc(r_dep,             v_sc_dep,          TOF1,   mu_sun)
leg2a = arc(r_mars,            v_sc_dep_mars_post, TOF_2a, mu_sun)
leg2b = arc(r_dsm,             v_dep_dsm,         TOF_2b, mu_sun)
leg3  = arc(r_ert,             v_sc_out_earth,    TOF3,   mu_sun)

ax.plot(leg1[:,0],  leg1[:,1],  leg1[:,2],  'w-',  lw=2.0, label='Leg 1 (E→Mars)',   alpha=0.95)
ax.plot(leg2a[:,0], leg2a[:,1], leg2a[:,2], 'y-',  lw=2.0, label='Leg 2a (Mars→DSM)', alpha=0.95)
ax.plot(leg2b[:,0], leg2b[:,1], leg2b[:,2], 'c-',  lw=2.0, label='Leg 2b (DSM→E)',   alpha=0.95)
ax.plot(leg3[:,0],  leg3[:,1],  leg3[:,2],  'm-',  lw=2.0, label='Leg 3 (E→Jup)',    alpha=0.95)

# Event markers
events = [
    (f'Earth dep\n{ds(t_dep)}',   r_dep,  '#4488ff'),
    (f'Mars flyby\n{ds(t_m)}',    r_mars, '#ff8844'),
    (f'DSM\n{ds(t_dsm)}',         r_dsm,  '#ffff44'),
    (f'Earth flyby\n{ds(t_ert)}', r_ert,  '#44ff88'),
    (f'Jup arr\n{ds(t_j)}',       r_jup,  '#c8a064'),
]
for label, r_m, col in events:
    r = r_m / AU
    ax.scatter(r[0], r[1], r[2], color=col, s=90, zorder=5, edgecolors='white', lw=0.5)
    ax.text(r[0], r[1], r[2], f'  {label.split(chr(10))[0]}', color=col, fontsize=7)

# Sun
ax.scatter(0, 0, 0, color='yellow', s=250, marker='*', edgecolors='orange', lw=0.5, label='Sun')

ax.set_xlabel('X (AU)', fontsize=9, color='white')
ax.set_ylabel('Y (AU)', fontsize=9, color='white')
ax.set_zlabel('Z (AU)', fontsize=9, color='white')
ax.set_title(
    f'EMEJ Trajectory — Launch {ds(t_dep)},  DSM {ds(t_dsm)}\n'
    f'Total Δv = {total_dv/1e3:.2f} km/s   '
    f'TOFs: {TOF1/DAY:.0f}d + {TOF_2a/DAY:.0f}d + {TOF_2b/DAY:.0f}d + {TOF3/DAY:.0f}d',
    fontsize=11, color='white'
)
ax.legend(loc='upper left', fontsize=7, facecolor='#222233', labelcolor='white')
ax.set_facecolor('#111122')
fig.patch.set_facecolor('#111122')
ax.tick_params(colors='white', labelsize=7)
ax.title.set_color('white')

lim = 7
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-2, 2)

plt.tight_layout()
plt.savefig('task2_emej_trajectory.png', dpi=150, facecolor='#111122', bbox_inches='tight')
print("  Saved: task2_emej_trajectory.png")
plt.show()

print()
print("  Done.")
