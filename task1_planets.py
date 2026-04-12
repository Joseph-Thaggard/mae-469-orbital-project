"""
task1_planets.py — Task 1: Planetary States on December 25, 2026 at 11:12 PM UTC

Outputs:
  - Table of heliocentric position (r), velocity (v), and true anomaly (θ)
    for all 8 planets + Pluto on Dec 25, 2026 at 23:12:00 UTC
  - 3D plot of all planet orbits with markers at their Dec 25, 2026 positions

Method:
  1. Start from J2000 orbital elements (Jan 1, 2000 12:00 TT ≈ 11:58:56 UTC)
  2. Compute elapsed seconds from J2000 to target date
  3. Use kepler_propagate (universal variables) to advance each body
  4. Convert back to orbital elements to get true anomaly at target date
"""

import numpy as np
import numpy.linalg as la
from datetime import datetime, timezone

from bodies.catalog import (
    make_mercury, make_venus, make_earth, make_mars,
    make_jupiter, make_saturn, make_uranus, make_neptune, make_pluto,
    make_sun, J2000_EPOCH_UTC, MU_SUN, AU
)
from solvers.solver import kepler_propagate, rv_to_elements, elements_to_rv

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------------------
# Target date
# ---------------------------------------------------------------------------
TARGET_UTC = datetime(2026, 12, 25, 23, 12, 0, tzinfo=timezone.utc)

# Elapsed seconds from J2000 epoch to target date
dt_seconds = (TARGET_UTC - J2000_EPOCH_UTC).total_seconds()
dt_days    = dt_seconds / 86400.0

print("=" * 72)
print("  Task 1: Planetary States on December 25, 2026 at 11:12 PM UTC")
print("=" * 72)
print(f"  J2000 epoch:    {J2000_EPOCH_UTC.strftime('%Y-%m-%d %H:%M:%S UTC')}")
print(f"  Target date:    {TARGET_UTC.strftime('%Y-%m-%d %H:%M:%S UTC')}")
print(f"  Elapsed time:   {dt_days:.4f} days  ({dt_days/365.25:.4f} Julian years)")

# ---------------------------------------------------------------------------
# Planetary catalog at J2000
# ---------------------------------------------------------------------------
bodies_j2000 = {
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
# Propagate to target date
# ---------------------------------------------------------------------------
states = {}   # name → dict(r, v, r_AU, v_kms, nu_deg)

for name, body in bodies_j2000.items():
    r_t, v_t = kepler_propagate(body.position, body.velocity, dt_seconds, MU_SUN)
    elems = rv_to_elements(r_t, v_t, MU_SUN)
    nu_deg = np.degrees(elems['nu'])
    states[name] = dict(
        r=r_t,
        v=v_t,
        r_AU=r_t / AU,
        v_kms=v_t / 1e3,
        nu_deg=nu_deg,
        r_mag_AU=la.norm(r_t) / AU,
        v_mag_kms=la.norm(v_t) / 1e3,
        a_AU=elems['a'] / AU,
        e=elems['e'],
    )

# ---------------------------------------------------------------------------
# Print table
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print(f"  Planet States at {TARGET_UTC.strftime('%B %d, %Y  %H:%M UTC')}")
print("=" * 72)
print()

# Header
print(f"  {'Planet':<10}  {'θ (deg)':>9}  {'|r| (AU)':>10}  {'|v| (km/s)':>11}")
print(f"  {'-'*10}  {'-'*9}  {'-'*10}  {'-'*11}")
for name, s in states.items():
    print(f"  {name:<10}  {s['nu_deg']:>9.4f}  {s['r_mag_AU']:>10.6f}  {s['v_mag_kms']:>11.5f}")

print()
print("  Full position and velocity vectors (AU and km/s):")
print()
print(f"  {'Planet':<10}  {'r_x (AU)':>12}  {'r_y (AU)':>12}  {'r_z (AU)':>12}"
      f"  {'v_x (km/s)':>12}  {'v_y (km/s)':>12}  {'v_z (km/s)':>12}  {'θ (deg)':>9}")
print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}"
      f"  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*9}")
for name, s in states.items():
    r = s['r_AU']
    v = s['v_kms']
    print(f"  {name:<10}  {r[0]:>12.6f}  {r[1]:>12.6f}  {r[2]:>12.6f}"
          f"  {v[0]:>12.6f}  {v[1]:>12.6f}  {v[2]:>12.6f}  {s['nu_deg']:>9.4f}")

# ---------------------------------------------------------------------------
# 3D Orbit + Position Plot
# ---------------------------------------------------------------------------
print()
print("  Generating 3D orbit plot ...")

fig = plt.figure(figsize=(14, 11))
ax  = fig.add_subplot(111, projection='3d')

# Color map for planets
COLORS = {
    'Mercury': '#b5b5b5',
    'Venus':   '#e8c86e',
    'Earth':   '#4488ff',
    'Mars':    '#cc4400',
    'Jupiter': '#c8a064',
    'Saturn':  '#d4c080',
    'Uranus':  '#88ddff',
    'Neptune': '#4466dd',
    'Pluto':   '#aaaaaa',
}

# Number of points for a full orbit
N_ORB = 500

for name, s in states.items():
    body = bodies_j2000[name]
    elems = rv_to_elements(s['r'], s['v'], MU_SUN)
    a, e, i, raan, argp = elems['a'], elems['e'], elems['i'], elems['raan'], elems['argp']

    # Sample full orbit: nu from 0 to 2π
    nu_vals = np.linspace(0, 2*np.pi, N_ORB)
    orb_pts = np.array([elements_to_rv(a, e, i, raan, argp, nu, MU_SUN)[0]
                        for nu in nu_vals]) / AU

    col = COLORS.get(name, 'white')
    ax.plot(orb_pts[:, 0], orb_pts[:, 1], orb_pts[:, 2],
            color=col, lw=0.8, alpha=0.7)

    # Planet position marker
    r = s['r_AU']
    ax.scatter(r[0], r[1], r[2], color=col, s=60, zorder=5, edgecolors='white', lw=0.5)
    ax.text(r[0], r[1], r[2], f'  {name}', color=col, fontsize=7)

# Sun marker
ax.scatter(0, 0, 0, color='yellow', s=200, zorder=6, marker='*', edgecolors='orange', lw=0.5)
ax.text(0, 0, 0, '  Sun', color='yellow', fontsize=8)

ax.set_xlabel('X (AU)', fontsize=9)
ax.set_ylabel('Y (AU)', fontsize=9)
ax.set_zlabel('Z (AU)', fontsize=9)
ax.set_title(f'Solar System — {TARGET_UTC.strftime("%B %d, %Y  %H:%M UTC")}', fontsize=12)
ax.set_facecolor('#111122')
fig.patch.set_facecolor('#111122')
ax.tick_params(colors='white', labelsize=7)
ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.zaxis.label.set_color('white')
ax.title.set_color('white')

# Make axis limits symmetric and include Pluto
lim = 45
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-lim/3, lim/3)

plt.tight_layout()
plt.savefig('task1_planet_orbits.png', dpi=150, facecolor='#111122', bbox_inches='tight')
print("  Saved: task1_planet_orbits.png")
plt.show()

print()
print("  Done.")
