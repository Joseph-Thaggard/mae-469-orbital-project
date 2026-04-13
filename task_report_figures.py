"""
task_report_figures.py — Generate all report-quality figures for MAE 469 Project

Figures produced:
  fig1_solar_system.png     — 3D solar system on Dec 25, 2026 (white bg, report quality)
  fig2_trajectory_3d.png    — 3D EMEJ trajectory (white bg, report quality)
  fig3_trajectory_2d.png    — 2D ecliptic top-down view of trajectory
  fig4_dv_vs_tof3.png       — Δv sensitivity: total Δv vs TOF3 (Earth→Jupiter leg)
  fig5_dv_vs_launch.png     — Δv sensitivity: total Δv vs launch date
  fig6_flyby_geometry.png   — Hyperbolic turn angle geometry for Mars and Earth flybys

All figures use white backgrounds suitable for printed reports.
Font size ≥ 10 pt per spec. Captions and numbering added in report separately.
"""

import numpy as np
import numpy.linalg as la
from datetime import datetime, timezone, timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from bodies.catalog import (
    make_earth, make_mars, make_jupiter,
    make_mercury, make_venus, make_saturn, make_uranus, make_neptune, make_pluto,
    J2000_EPOCH_UTC, MU_SUN, AU
)
from solvers.solver import kepler_propagate, lambert_solver, rv_to_elements, elements_to_rv
from solvers.dsm import optimize_dsm

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
DAY = 86400.0
plt.rcParams.update({
    'font.family':  'Arial',
    'font.size':    11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi':   150,
})

# Body instances
earth   = make_earth()
mars    = make_mars()
jupiter = make_jupiter()
mu_sun  = MU_SUN
r_e0, v_e0 = earth.position.copy(), earth.velocity.copy()
r_m0, v_m0 = mars.position.copy(),  mars.velocity.copy()
r_j0, v_j0 = jupiter.position.copy(), jupiter.velocity.copy()

def j2k(dt): return (dt - J2000_EPOCH_UTC).total_seconds()
def ds(t_s): return (J2000_EPOCH_UTC + timedelta(seconds=t_s)).strftime('%Y-%m-%d')

# ========================================================================================
# Shared trajectory computation (same as task2_emej.py, reproduced for self-containment)
# ========================================================================================
T_LAUNCH = j2k(datetime(2027, 1, 15, tzinfo=timezone.utc))
TOF1 = 260*DAY; TOF2 = 420*DAY; TOF3 = 1100*DAY
t_dep = T_LAUNCH; t_m = t_dep+TOF1; t_ert = t_m+TOF2; t_j = t_ert+TOF3

r_dep,  v_dep  = kepler_propagate(r_e0, v_e0, t_dep, mu_sun)
r_mars, v_mars = kepler_propagate(r_m0, v_m0, t_m,   mu_sun)
r_ert,  v_ert  = kepler_propagate(r_e0, v_e0, t_ert, mu_sun)
r_jup,  v_jup  = kepler_propagate(r_j0, v_j0, t_j,   mu_sun)

v_sc_dep, v_sc_in_mars   = lambert_solver(r_dep,  r_mars, TOF1, mu_sun)
v_sc_out_earth, v_sc_arr_jup = lambert_solver(r_ert, r_jup,  TOF3, mu_sun)

v_inf_in_mars = v_sc_in_mars - v_mars
k_in_mars = la.norm(v_inf_in_mars)
r_mars_fb = mars.radius + 200e3
ecc_m = 1.0 + r_mars_fb * k_in_mars**2 / mars.mu
delta_max_mars = 2.0 * np.arcsin(1.0 / ecc_m)
r_earth_fb_min = earth.radius + 300e3

dsm = optimize_dsm(
    r_mars, v_mars, v_inf_in_mars, k_in_mars, delta_max_mars,
    r_ert, v_ert, TOF2, v_sc_out_earth, mu_sun,
    mu_target=earth.mu, rp_min_target=r_earth_fb_min,
    lam_earth=1.0, lam_turn=1.0,
    n_tau=10, n_phi=12, n_delta=4, refine=True,
)

r_dsm              = dsm['r_dsm']
v_at_dsm           = dsm['v_at_dsm']
v_dep_dsm          = dsm['v_dep_dsm']
v_arr_ert          = dsm['v_arr_target']
v_sc_dep_mars_post = dsm['v_sc_dep_flyby']
v_inf_out_mars_fb  = dsm['v_inf_out']
v_inf_in_earth     = dsm['v_inf_in_target']
v_inf_out_earth    = dsm['v_inf_out_target']
v_inf_arr_jup      = v_sc_arr_jup - v_jup
TOF_2a = dsm['tof_a']; TOF_2b = dsm['tof_b']
t_dsm  = t_m + TOF_2a
dv_dsm = dsm['dv_dsm']
v_inf_dep = v_sc_dep - v_dep

r_leo = earth.radius + 500e3
dv1 = np.sqrt(la.norm(v_inf_dep)**2 + 2*earth.mu/r_leo) - np.sqrt(earth.mu/r_leo)
r_jp = jupiter.radius + 10_000e3
dv3 = np.sqrt(la.norm(v_inf_arr_jup)**2 + 2*jupiter.mu/r_jp) - np.sqrt(jupiter.mu/r_jp)
total_dv = dv1 + dv_dsm + dv3

def arc(r0, v0, tof, mu=mu_sun, n=300):
    ts = np.linspace(0, tof, n)
    return np.array([kepler_propagate(r0, v0, t, mu)[0] / AU for t in ts])

def full_orbit(r0, v0, n=500):
    el = rv_to_elements(r0, v0, mu_sun)
    a, e, i, raan, argp = el['a'], el['e'], el['i'], el['raan'], el['argp']
    return np.array([elements_to_rv(a, e, i, raan, argp, nu, mu_sun)[0] / AU
                     for nu in np.linspace(0, 2*np.pi, n)])

print("Computing trajectory arcs...")
leg1  = arc(r_dep,             v_sc_dep,          TOF1)
leg2a = arc(r_mars,            v_sc_dep_mars_post, TOF_2a)
leg2b = arc(r_dsm,             v_dep_dsm,         TOF_2b)
leg3  = arc(r_ert,             v_sc_out_earth,    TOF3)

# Planet orbit traces
orb_earth   = full_orbit(r_e0, v_e0)
orb_mars    = full_orbit(r_m0, v_m0)
orb_jupiter = full_orbit(r_j0, v_j0)

# Planet colors
C = {
    'Earth':   '#1f77b4',
    'Mars':    '#d62728',
    'Jupiter': '#ff7f0e',
    'DSM':     '#8c564b',
    'SC':      '#2ca02c',
}

# ========================================================================================
# FIG 1 — 3D Solar System on Dec 25, 2026 (Task 1, report quality)
# ========================================================================================
print("Generating Fig 1: 3D solar system (Dec 25, 2026)...")

TARGET_UTC = datetime(2026, 12, 25, 23, 12, 0, tzinfo=timezone.utc)
dt_s = (TARGET_UTC - J2000_EPOCH_UTC).total_seconds()

all_bodies = {
    'Mercury': make_mercury(), 'Venus': make_venus(), 'Earth': make_earth(),
    'Mars': make_mars(), 'Jupiter': make_jupiter(), 'Saturn': make_saturn(),
    'Uranus': make_uranus(), 'Neptune': make_neptune(), 'Pluto': make_pluto(),
}
body_colors = {
    'Mercury': '#888888', 'Venus': '#c8a028', 'Earth': '#1f77b4',
    'Mars': '#d62728', 'Jupiter': '#ff7f0e', 'Saturn': '#bcbd22',
    'Uranus': '#17becf', 'Neptune': '#7f7f7f', 'Pluto': '#9467bd',
}

fig1 = plt.figure(figsize=(12, 9))
ax1  = fig1.add_subplot(111, projection='3d')
ax1.set_facecolor('white')
fig1.patch.set_facecolor('white')

for name, body in all_bodies.items():
    r_t, v_t = kepler_propagate(body.position, body.velocity, dt_s, mu_sun)
    el = rv_to_elements(r_t, v_t, mu_sun)
    a, e, i, raan, argp = el['a'], el['e'], el['i'], el['raan'], el['argp']
    orb = np.array([elements_to_rv(a, e, i, raan, argp, nu, mu_sun)[0] / AU
                    for nu in np.linspace(0, 2*np.pi, 400)])
    col = body_colors[name]
    ax1.plot(orb[:, 0], orb[:, 1], orb[:, 2], color=col, lw=0.9, alpha=0.6)
    r = r_t / AU
    ax1.scatter(r[0], r[1], r[2], color=col, s=55, zorder=5,
                edgecolors='k', lw=0.4)
    ax1.text(r[0]+0.3, r[1]+0.3, r[2], name, color=col, fontsize=8, fontweight='bold')

ax1.scatter(0, 0, 0, color='gold', s=250, marker='*', zorder=6, edgecolors='darkorange', lw=0.5)
ax1.text(0.5, 0.5, 0, 'Sun', fontsize=9, fontweight='bold', color='darkorange')

ax1.set_xlabel('X (AU)', labelpad=6)
ax1.set_ylabel('Y (AU)', labelpad=6)
ax1.set_zlabel('Z (AU)', labelpad=6)
ax1.set_title('Solar System — December 25, 2026 at 23:12 UTC', pad=14)
lim = 45
ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim); ax1.set_zlim(-lim/3, lim/3)
ax1.tick_params(labelsize=9)
plt.tight_layout()
plt.savefig('fig1_solar_system.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig1_solar_system.png")
plt.close()

# ========================================================================================
# FIG 2 — 3D EMEJ Trajectory (report quality, white background)
# ========================================================================================
print("Generating Fig 2: 3D EMEJ trajectory...")

fig2 = plt.figure(figsize=(13, 9))
ax2  = fig2.add_subplot(111, projection='3d')
ax2.set_facecolor('white')
fig2.patch.set_facecolor('white')

# Planet orbits
ax2.plot(orb_earth[:,0],   orb_earth[:,1],   orb_earth[:,2],   color=C['Earth'],   lw=1.0, alpha=0.5, label="Earth's orbit")
ax2.plot(orb_mars[:,0],    orb_mars[:,1],     orb_mars[:,2],    color=C['Mars'],    lw=1.0, alpha=0.5, label="Mars' orbit")
ax2.plot(orb_jupiter[:,0], orb_jupiter[:,1],  orb_jupiter[:,2], color=C['Jupiter'], lw=1.0, alpha=0.5, label="Jupiter's orbit")

# Spacecraft arcs
ax2.plot(leg1[:,0],  leg1[:,1],  leg1[:,2],  color='black',      lw=2.2, label='Leg 1: Earth→Mars',   zorder=4)
ax2.plot(leg2a[:,0], leg2a[:,1], leg2a[:,2], color='#e377c2',    lw=2.2, label='Leg 2a: Mars→DSM',    zorder=4)
ax2.plot(leg2b[:,0], leg2b[:,1], leg2b[:,2], color='#17becf',    lw=2.2, label='Leg 2b: DSM→Earth',   zorder=4)
ax2.plot(leg3[:,0],  leg3[:,1],  leg3[:,2],  color='#9467bd',    lw=2.2, label='Leg 3: Earth→Jupiter', zorder=4)

# Event markers
events_3d = [
    ('Earth dep.\n2027-01-15',    r_dep,  C['Earth'],   '^', 90),
    ('Mars flyby\n2027-10-02',    r_mars, C['Mars'],    'o', 90),
    ('DSM\n2027-10-10',           r_dsm,  C['DSM'],     's', 60),
    ('Earth flyby\n2028-11-25',   r_ert,  C['Earth'],   'v', 90),
    ('Jupiter arr.\n2031-11-30',  r_jup,  C['Jupiter'], 'D', 110),
]
for label, r_m, col, mk, sz in events_3d:
    r = r_m / AU
    ax2.scatter(r[0], r[1], r[2], color=col, s=sz, marker=mk, zorder=6,
                edgecolors='k', lw=0.6)
    short = label.split('\n')[0]
    ax2.text(r[0]+0.1, r[1]+0.1, r[2]+0.05, short, fontsize=8, color='black',
             fontweight='bold')

ax2.scatter(0, 0, 0, color='gold', s=250, marker='*', zorder=7,
            edgecolors='darkorange', lw=0.5, label='Sun')

ax2.set_xlabel('X (AU)', labelpad=6)
ax2.set_ylabel('Y (AU)', labelpad=6)
ax2.set_zlabel('Z (AU)', labelpad=6)
ax2.set_title(
    f'EMEJ Trajectory — Earth–Mars–Earth–Jupiter\n'
    f'Launch: 2027-01-15  ·  Total Δv = {total_dv/1e3:.2f} km/s  ·  Mission: {(t_j-t_dep)/DAY:.0f} days',
    pad=14
)
ax2.legend(loc='upper left', framealpha=0.85, fontsize=8)
lim2 = 6.5
ax2.set_xlim(-lim2, lim2); ax2.set_ylim(-lim2, lim2); ax2.set_zlim(-1.5, 1.5)
ax2.tick_params(labelsize=9)
plt.tight_layout()
plt.savefig('fig2_trajectory_3d.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig2_trajectory_3d.png")
plt.close()

# ========================================================================================
# FIG 3 — 2D Ecliptic Top-Down View
# ========================================================================================
print("Generating Fig 3: 2D ecliptic projection...")

fig3, ax3 = plt.subplots(figsize=(10, 9))
ax3.set_facecolor('white')
fig3.patch.set_facecolor('white')
ax3.set_aspect('equal')

# Planet orbits (X-Y only)
ax3.plot(orb_earth[:,0],   orb_earth[:,1],   color=C['Earth'],   lw=1.0, alpha=0.5)
ax3.plot(orb_mars[:,0],    orb_mars[:,1],     color=C['Mars'],    lw=1.0, alpha=0.5)
ax3.plot(orb_jupiter[:,0], orb_jupiter[:,1],  color=C['Jupiter'], lw=1.0, alpha=0.5)

# Spacecraft arcs
ax3.plot(leg1[:,0],  leg1[:,1],  color='black',   lw=2.0, label='Leg 1: Earth→Mars',   zorder=4)
ax3.plot(leg2a[:,0], leg2a[:,1], color='#e377c2', lw=2.0, label='Leg 2a: Mars→DSM',    zorder=4)
ax3.plot(leg2b[:,0], leg2b[:,1], color='#17becf', lw=2.0, label='Leg 2b: DSM→Earth',   zorder=4)
ax3.plot(leg3[:,0],  leg3[:,1],  color='#9467bd', lw=2.0, label='Leg 3: Earth→Jupiter', zorder=4)

# Event markers
offsets = {
    'Earth dep.':   (+0.15, +0.15),
    'Mars flyby':   (+0.10, -0.20),
    'DSM':          (+0.10, -0.20),
    'Earth flyby':  (+0.15, +0.15),
    'Jupiter arr.': (+0.15, +0.10),
}
for label, r_m, col, mk, sz in events_3d:
    r = r_m / AU
    ax3.scatter(r[0], r[1], color=col, s=sz, marker=mk, zorder=6,
                edgecolors='k', lw=0.6)
    short = label.split('\n')[0]
    ox, oy = offsets.get(short, (0.15, 0.15))
    ax3.annotate(f"{short}\n{label.split(chr(10))[1]}", xy=(r[0], r[1]),
                 xytext=(r[0]+ox, r[1]+oy), fontsize=8, color='black',
                 fontweight='bold',
                 arrowprops=dict(arrowstyle='-', color='gray', lw=0.8))

# Planet labels at orbit midpoints
for name, orb, col in [('Earth', orb_earth, C['Earth']),
                        ('Mars', orb_mars, C['Mars']),
                        ('Jupiter', orb_jupiter, C['Jupiter'])]:
    idx = len(orb)//4
    ax3.text(orb[idx,0]+0.1, orb[idx,1]+0.1, name, color=col, fontsize=9,
             fontstyle='italic', alpha=0.8)

# Sun
ax3.scatter(0, 0, color='gold', s=200, marker='*', zorder=7,
            edgecolors='darkorange', lw=0.5)
ax3.text(0.1, 0.1, 'Sun', fontsize=9, fontweight='bold', color='darkorange')

# Add tick marks for scale
ax3.axhline(0, color='lightgray', lw=0.5, zorder=0)
ax3.axvline(0, color='lightgray', lw=0.5, zorder=0)

ax3.set_xlabel('X (AU)', fontsize=11)
ax3.set_ylabel('Y (AU)', fontsize=11)
ax3.set_title(
    'EMEJ Trajectory — Ecliptic Plane Projection\n'
    'Earth–Mars–Earth–Jupiter (2027–2031)',
    fontsize=13
)
ax3.legend(loc='lower right', framealpha=0.9, fontsize=9)
ax3.set_xlim(-7.5, 7.5); ax3.set_ylim(-7.5, 7.5)
ax3.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig('fig3_trajectory_2d.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig3_trajectory_2d.png")
plt.close()

# ========================================================================================
# FIG 4 — Δv Sensitivity: Total Δv vs TOF3 (Earth→Jupiter)
# (Bonus optimization analysis — holds TOF1=260d, TOF2=420d, launch=2027-01-15)
# ========================================================================================
print("Generating Fig 4: Δv sensitivity vs TOF3...")

tof3_range = np.arange(900, 1251, 25)
dv_totals, dv1_arr, dv3_arr, feasible = [], [], [], []

for tof3_d in tof3_range:
    tof3_test = tof3_d * DAY
    t_j_test  = t_ert + tof3_test
    try:
        r_jt, v_jt = kepler_propagate(r_j0, v_j0, t_j_test, mu_sun)
        v_oe, v_aj = lambert_solver(r_ert, r_jt, tof3_test, mu_sun)
    except Exception:
        dv_totals.append(np.nan); dv1_arr.append(np.nan)
        dv3_arr.append(np.nan);   feasible.append(False)
        continue

    v_inf_j = v_aj - v_jt
    dv3_t = np.sqrt(la.norm(v_inf_j)**2 + 2*jupiter.mu/r_jp) - np.sqrt(jupiter.mu/r_jp)

    # Earth flyby check: use v_oe as post-flyby departure direction
    v_inf_in_e  = v_arr_ert - v_ert     # from fixed DSM result
    v_inf_out_e = v_oe - v_ert
    k_ie = la.norm(v_inf_in_e); k_oe = la.norm(v_inf_out_e)
    ecc_e = 1.0 + r_earth_fb_min * k_ie**2 / earth.mu
    max_turn_e = 2.0 * np.arcsin(1.0 / ecc_e)
    cos_d = np.clip(np.dot(v_inf_in_e, v_inf_out_e)/(k_ie*k_oe), -1, 1)
    turn_e = np.arccos(cos_d)
    sin_h = np.sin(turn_e/2)
    rp_req = earth.mu/k_ie**2*(1/sin_h-1) if sin_h > 1e-9 else np.inf
    ok = (turn_e <= max_turn_e) and (rp_req >= r_earth_fb_min) and (abs(k_ie-k_oe) < 500)

    dv_totals.append((dv1 + dv_dsm + dv3_t) / 1e3)
    dv1_arr.append(dv1 / 1e3)
    dv3_arr.append(dv3_t / 1e3)
    feasible.append(ok)

tof3_range_d = tof3_range
dv_totals = np.array(dv_totals)
dv3_arr   = np.array(dv3_arr)
feasible  = np.array(feasible)

fig4, (ax4a, ax4b) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
fig4.patch.set_facecolor('white')

# Top: component breakdown
ax4a.plot(tof3_range_d, dv3_arr, 'o-', color=C['Jupiter'], lw=2, ms=5, label=r'$\Delta v_3$ (Jupiter arrival)')
dv_dsm_line = np.full_like(dv_totals, dv_dsm/1e3)
ax4a.axhline(dv1/1e3, color=C['Earth'], lw=1.5, ls='--', label=rf'$\Delta v_1$ (Earth dep.) = {dv1/1e3:.2f} km/s')
ax4a.axhline(dv_dsm/1e3, color=C['DSM'], lw=1.5, ls=':', label=rf'$\Delta v_2$ (DSM) = {dv_dsm/1e3:.2f} km/s')
ax4a.set_ylabel('Δv Component (km/s)', fontsize=11)
ax4a.set_title('Δv Sensitivity to Earth→Jupiter Transfer Duration\n(Launch: 2027-01-15, TOF1=260d, TOF2=420d)', fontsize=12)
ax4a.legend(fontsize=9)
ax4a.grid(True, alpha=0.3)
ax4a.set_ylim(0, 25)

# Bottom: total Δv with feasibility shading
ax4b.plot(tof3_range_d, dv_totals, 'k-o', lw=2, ms=5, label='Total Δv')
feas_idx = feasible & ~np.isnan(dv_totals)
infeas_idx = ~feasible & ~np.isnan(dv_totals)
ax4b.scatter(tof3_range_d[feas_idx],   dv_totals[feas_idx],   color='green', s=60, zorder=5,
             label='Feasible Earth flyby (OK)')
ax4b.scatter(tof3_range_d[infeas_idx], dv_totals[infeas_idx], color='red',   s=60, marker='x',
             zorder=5, label='Earth flyby infeasible')
# Mark selected
ax4b.axvline(TOF3/DAY, color='navy', lw=1.5, ls='--', label=f'Selected: {TOF3/DAY:.0f}d')
ax4b.scatter([TOF3/DAY], [(dv1+dv_dsm+dv3)/1e3], color='navy', s=120, marker='*', zorder=7)
ax4b.set_xlabel(r'$TOF_3$ — Earth to Jupiter (days)', fontsize=11)
ax4b.set_ylabel('Total Mission Δv (km/s)', fontsize=11)
ax4b.legend(fontsize=9)
ax4b.grid(True, alpha=0.3)
ax4b.set_ylim(30, 42)

plt.tight_layout()
plt.savefig('fig4_dv_vs_tof3.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig4_dv_vs_tof3.png")
plt.close()

# ========================================================================================
# FIG 5 — Δv Sensitivity: Total Δv vs Launch Date
# (Holds TOF1=260d, TOF2=420d, TOF3=1100d — sweeps launch date ±90 days)
# ========================================================================================
print("Generating Fig 5: Δv sensitivity vs launch date...")

launch_offsets = np.arange(-90, 91, 10)   # days from 2027-01-15
dv1_scan, dv3_scan, dv_tot_scan = [], [], []
launch_labels = []

for off in launch_offsets:
    t_dep_s = T_LAUNCH + off*DAY
    t_m_s   = t_dep_s + TOF1
    t_ert_s = t_m_s + TOF2
    t_j_s   = t_ert_s + TOF3
    try:
        r_d, v_d = kepler_propagate(r_e0, v_e0, t_dep_s, mu_sun)
        r_ms, v_ms = kepler_propagate(r_m0, v_m0, t_m_s, mu_sun)
        r_es, v_es = kepler_propagate(r_e0, v_e0, t_ert_s, mu_sun)
        r_js, v_js = kepler_propagate(r_j0, v_j0, t_j_s, mu_sun)
        v_dep_s2, _ = lambert_solver(r_d, r_ms, TOF1, mu_sun)
        v_oe_s, v_aj_s = lambert_solver(r_es, r_js, TOF3, mu_sun)
    except Exception:
        dv1_scan.append(np.nan); dv3_scan.append(np.nan); dv_tot_scan.append(np.nan)
        launch_labels.append(off)
        continue

    vi_dep = v_dep_s2 - v_d
    vi_jup = v_aj_s   - v_js
    dv1_s = np.sqrt(la.norm(vi_dep)**2 + 2*earth.mu/r_leo) - np.sqrt(earth.mu/r_leo)
    dv3_s = np.sqrt(la.norm(vi_jup)**2 + 2*jupiter.mu/r_jp) - np.sqrt(jupiter.mu/r_jp)
    dv1_scan.append(dv1_s/1e3); dv3_scan.append(dv3_s/1e3)
    dv_tot_scan.append((dv1_s + dv_dsm + dv3_s)/1e3)
    launch_labels.append(off)

dv1_scan = np.array(dv1_scan); dv3_scan = np.array(dv3_scan); dv_tot_scan = np.array(dv_tot_scan)

fig5, ax5 = plt.subplots(figsize=(10, 6))
fig5.patch.set_facecolor('white')
ax5.plot(launch_offsets, dv1_scan, 's--', color=C['Earth'], lw=1.5, ms=5, label=r'$\Delta v_1$ (Earth departure)')
ax5.plot(launch_offsets, dv3_scan, 'D--', color=C['Jupiter'], lw=1.5, ms=5, label=r'$\Delta v_3$ (Jupiter arrival)')
ax5.plot(launch_offsets, dv_tot_scan, 'k-o', lw=2.2, ms=5, label=r'Total $\Delta v$ ($\Delta v_1$ + DSM + $\Delta v_3$)')
ax5.axvline(0, color='navy', lw=1.5, ls='--', label='Selected launch (2027-01-15)')
ax5.scatter([0], [(dv1+dv_dsm+dv3)/1e3], color='navy', s=120, marker='*', zorder=7)
ax5.set_xlabel('Days from Selected Launch Date (2027-01-15)', fontsize=11)
ax5.set_ylabel('Δv (km/s)', fontsize=11)
ax5.set_title('Δv Sensitivity to Launch Date\n(TOF1=260d, TOF2=420d, TOF3=1100d fixed)', fontsize=12)
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)
ax5_x2 = ax5.twiny()
ax5_x2.set_xlim(ax5.get_xlim())
tick_offs = [-90, -60, -30, 0, 30, 60, 90]
tick_labels = [(datetime(2027,1,15,tzinfo=timezone.utc) + timedelta(days=o)).strftime('%b %d\n%Y')
               for o in tick_offs]
ax5_x2.set_xticks(tick_offs); ax5_x2.set_xticklabels(tick_labels, fontsize=8)
plt.tight_layout()
plt.savefig('fig5_dv_vs_launch.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig5_dv_vs_launch.png")
plt.close()

# ========================================================================================
# FIG 6 — Flyby Geometry Diagrams (Mars and Earth)
# ========================================================================================
print("Generating Fig 6: Flyby geometry diagrams...")

fig6, (ax6a, ax6b) = plt.subplots(1, 2, figsize=(13, 6))
fig6.patch.set_facecolor('white')

def draw_flyby(ax, v_inf_in, v_inf_out, body_name, body_color, rp_m, radius_m, mu_body, title):
    k_in  = la.norm(v_inf_in); k_out = la.norm(v_inf_out)
    cos_d = np.clip(np.dot(v_inf_in, v_inf_out)/(k_in*k_out), -1, 1)
    turn  = np.degrees(np.arccos(cos_d))
    ecc   = 1.0 + rp_m * k_in**2 / mu_body
    max_turn = np.degrees(2.0 * np.arcsin(1.0/ecc))

    # Normalize to 2D (use x-y components)
    vin  = np.array([v_inf_in[0],  v_inf_in[1]])  / k_in
    vout = np.array([v_inf_out[0], v_inf_out[1]]) / k_out

    scale = 2.2
    ax.annotate('', xy=vin*scale*0.3, xytext=-vin*scale,
                arrowprops=dict(arrowstyle='->', color='steelblue', lw=2))
    ax.annotate('', xy=vout*scale, xytext=vout*scale*0.3,
                arrowprops=dict(arrowstyle='->', color='tomato', lw=2))
    ax.text(-vin[0]*scale*1.05, -vin[1]*scale*1.05,
            f'v∞_in\n{k_in/1e3:.2f} km/s', color='steelblue', fontsize=9,
            ha='center', va='center')
    ax.text(vout[0]*scale*1.15, vout[1]*scale*1.15,
            f'v∞_out\n{k_out/1e3:.2f} km/s', color='tomato', fontsize=9,
            ha='center', va='center')

    # Draw arc indicating turn angle
    ang1 = np.degrees(np.arctan2(-vin[1], -vin[0]))
    ang2 = np.degrees(np.arctan2(vout[1], vout[0]))
    theta_arc = np.linspace(np.radians(ang1), np.radians(ang2), 80)
    r_arc = 0.7
    ax.plot(r_arc*np.cos(theta_arc), r_arc*np.sin(theta_arc), 'k-', lw=1.5)
    mid_arc = (ang1 + ang2) / 2
    ax.text(r_arc*1.2*np.cos(np.radians(mid_arc)), r_arc*1.2*np.sin(np.radians(mid_arc)),
            f'δ = {turn:.1f}°\n(max {max_turn:.1f}°)', fontsize=9, ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', edgecolor='gray', alpha=0.9))

    # Body circle (scaled)
    circ = plt.Circle((0, 0), 0.15, color=body_color, zorder=5, alpha=0.9)
    ax.add_patch(circ)
    ax.text(0, 0, body_name, color='white', fontsize=9, ha='center', va='center',
            fontweight='bold')

    # rp annotation
    ax.annotate('', xy=(0.15, 0), xytext=(0.15+rp_m/radius_m*0.05, 0),
                arrowprops=dict(arrowstyle='<->', color='gray', lw=1))
    ax.text(0.22, -0.15, f'$r_p$ = {rp_m/1e3:.0f} km', fontsize=8, color='gray', ha='left')

    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=11, pad=10)

# Mars flyby
v_inf_in_mars_2d  = np.array([*v_inf_in_mars[:2],  0])
v_inf_out_mars_2d = np.array([*v_inf_out_mars_fb[:2], 0])
draw_flyby(ax6a, v_inf_in_mars_2d, v_inf_out_mars_2d,
           'Mars', C['Mars'], mars.radius+200e3, mars.radius, mars.mu,
           'Mars Gravity Assist (Free Flyby)\n2027-10-02')

# Earth flyby
v_inf_in_earth_2d  = np.array([*v_inf_in_earth[:2],  0])
v_inf_out_earth_2d = np.array([*v_inf_out_earth[:2], 0])
rp_e = dsm['rp_required'] if dsm['rp_required'] else earth.radius+300e3
draw_flyby(ax6b, v_inf_in_earth_2d, v_inf_out_earth_2d,
           'Earth', C['Earth'], rp_e, earth.radius, earth.mu,
           'Earth Gravity Assist (Free Flyby)\n2028-11-25')

fig6.suptitle('Hyperbolic Flyby Geometry — Ecliptic Plane Projection',
              fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('fig6_flyby_geometry.png', dpi=150, bbox_inches='tight', facecolor='white')
print("  Saved: fig6_flyby_geometry.png")
plt.close()

# ========================================================================================
print()
print("All report figures generated:")
print("  fig1_solar_system.png    — 3D solar system, Dec 25 2026")
print("  fig2_trajectory_3d.png   — 3D EMEJ trajectory")
print("  fig3_trajectory_2d.png   — 2D ecliptic projection")
print("  fig4_dv_vs_tof3.png      — Δv sensitivity vs TOF3")
print("  fig5_dv_vs_launch.png    — Δv sensitivity vs launch date")
print("  fig6_flyby_geometry.png  — Mars & Earth flyby geometry")
print()
print(f"Final mission Δv:  {total_dv/1e3:.3f} km/s")
print(f"  Δv1 = {dv1/1e3:.3f} km/s  |  Δv2 = {dv_dsm/1e3:.3f} km/s  |  Δv3 = {dv3/1e3:.3f} km/s")
