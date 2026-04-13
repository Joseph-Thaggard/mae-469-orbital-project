"""
scan_emej_freeflyby.py — Find EMEJ windows with a near-free Mars gravity assist
                         and optimized DSM placement.

Two-stage search:

  Stage 1 — Sweep (t0, TOF1) over the full Earth-Mars synodic period.
             For each combination, find the TOF2 (Mars→Earth) that minimises
             the v∞ mismatch at Mars (|v∞_in| - |v∞_out|).

  Stage 2 — For every Stage-1 survivor, sweep TOF3 (Earth→Jupiter).
             For each (survivor, TOF3): optimise DSM placement (tau, phi, delta)
             via dsm_grid_search() to find minimum total Δv with a near-free
             Earth flyby.  Rank by Δv1 + dv_DSM + Δv3.

DSM model: Mars free flyby provides a turn of up to δ_max. Spacecraft coasts
for tau*TOF2, then fires DSM to redirect onto a Lambert arc to Earth. The Earth
flyby v∞ mismatch is penalised in the DSM cost function.

Usage:
    python3 scan_emej_freeflyby.py
"""

import numpy as np
import numpy.linalg as la

from bodies.catalog import make_sun, make_earth, make_mars, make_jupiter, FLYBY_PERIAPSIS_MIN
from solvers.solver import kepler_propagate, lambert_solver
from solvers.dsm import dsm_grid_search

DAY = 86400.0
AU  = 1.496e11

# ============================================================
# CONFIGURATION
# ============================================================

EPOCH_OFFSET_DAYS = 300.0        # days — same as workflow_emej.py

# Stage 1: Earth→Mars leg sweep
T0_MIN, T0_MAX, T0_STEP    = 0.0,   800.0,  5.0    # days from epoch
TOF1_MIN, TOF1_MAX, TOF1_STEP = 150., 500., 10.0   # Earth→Mars TOF (days)

# Stage 1: Mars→Earth sweep for each (t0, TOF1)
TOF2_MIN, TOF2_MAX, TOF2_STEP = 200., 900., 10.0   # Mars→Earth TOF (days)

# Stage 1 filter: keep windows where |Δv∞_mars| < threshold (m/s)
MARS_FREE_THRESHOLD = 500.0   # m/s  (0 = perfect free flyby, ~500 practical)

# Stage 2: Earth→Jupiter sweep (applied to Stage-1 survivors only)
TOF3_MIN, TOF3_MAX, TOF3_STEP = 600., 1600., 50.0  # Earth→Jupiter TOF (days)

# Periapsis altitudes
MARS_ALT   = 400e3     # 400 km above Mars surface
EARTH_ALT  = 300e3     # 300 km — LEO clearance

N_RESULTS = 15
VERBOSE   = True

# ============================================================
# Body setup (identical to workflow_emej.py)
# ============================================================

sun     = make_sun()
earth   = make_earth()
mars    = make_mars()
jupiter = make_jupiter()
mu_sun  = sun.mu

t_epoch = EPOCH_OFFSET_DAYS * DAY
r_e0, v_e0 = kepler_propagate(earth.position,   earth.velocity,   t_epoch, mu_sun)
r_m0, v_m0 = kepler_propagate(mars.position,    mars.velocity,    t_epoch, mu_sun)
r_j0, v_j0 = kepler_propagate(jupiter.position, jupiter.velocity, t_epoch, mu_sun)

class _Body:
    def __init__(self, name, mu, position, velocity):
        self.name = name; self.mu = mu
        self.position = position; self.velocity = velocity

rp_min_mars  = FLYBY_PERIAPSIS_MIN['Mars']
rp_min_earth = earth.radius + EARTH_ALT

def body_state(r0, v0, t):
    return kepler_propagate(r0, v0, t, mu_sun)

def fn_earth1(t): return body_state(r_e0, v_e0, t)
def fn_mars(t):   return body_state(r_m0, v_m0, t)
def fn_earth2(t): return body_state(r_e0, v_e0, t)
def fn_jupiter(t): return body_state(r_j0, v_j0, t)

print("=" * 68)
print("  EMEJ Free-Flyby Window Scanner")
print("=" * 68)
print(f"  Epoch offset:      {EPOCH_OFFSET_DAYS:.0f} days")
print(f"  Mars free-flyby threshold: {MARS_FREE_THRESHOLD:.0f} m/s")
print(f"  Stage-1 grid: {len(np.arange(T0_MIN,T0_MAX,T0_STEP)):.0f} t0 × "
      f"{len(np.arange(TOF1_MIN,TOF1_MAX,TOF1_STEP)):.0f} tof1 × "
      f"{len(np.arange(TOF2_MIN,TOF2_MAX,TOF2_STEP)):.0f} tof2 = "
      f"{len(np.arange(T0_MIN,T0_MAX,T0_STEP))*len(np.arange(TOF1_MIN,TOF1_MAX,TOF1_STEP))*len(np.arange(TOF2_MIN,TOF2_MAX,TOF2_STEP)):.0f} combos")

# ============================================================
# Stage 1 — Find near-free Mars flyby windows
# ============================================================

t0_arr   = np.arange(T0_MIN,   T0_MAX,   T0_STEP)   * DAY
tof1_arr = np.arange(TOF1_MIN, TOF1_MAX, TOF1_STEP) * DAY
tof2_arr = np.arange(TOF2_MIN, TOF2_MAX, TOF2_STEP) * DAY

survivors = []   # (dv_mars_mismatch, t0, tof1, tof2, v_inf_in_vec, v_inf_out_vec, r_mars, r_earth_ret)

n_total   = len(t0_arr)
for i_t0, t0 in enumerate(t0_arr):
    if VERBOSE and i_t0 % max(1, n_total // 10) == 0:
        pct = 100 * i_t0 // n_total
        print(f"  Stage 1 — t0={t0/DAY:.0f} d ({pct}%)", flush=True)

    t1_arr_local = t0 + tof1_arr

    # Departure body state (Earth at t0)
    r_dep, v_dep = fn_earth1(t0)

    for tof1 in tof1_arr:
        t1 = t0 + tof1
        r_mars_t1, v_mars_t1 = fn_mars(t1)

        # Leg 1: Earth → Mars
        try:
            v_sc_out1, v_sc_in1 = lambert_solver(r_dep, r_mars_t1, tof1, mu_sun)
        except Exception:
            continue

        v_inf_in_mars = v_sc_in1 - v_mars_t1
        v_inf_in_mag  = la.norm(v_inf_in_mars)

        # Leg 2: Mars → Earth — sweep TOF2, minimise |v∞_in - v∞_out| at Mars
        best_mismatch = np.inf
        best_tof2     = None
        best_vout     = None
        best_r_earth  = None

        for tof2 in tof2_arr:
            t2 = t1 + tof2
            r_earth_t2, v_earth_t2 = fn_earth2(t2)

            try:
                v_sc_out2, v_sc_in2 = lambert_solver(r_mars_t1, r_earth_t2, tof2, mu_sun)
            except Exception:
                continue

            v_inf_out_mars = v_sc_out2 - v_mars_t1
            mismatch = abs(la.norm(v_inf_out_mars) - v_inf_in_mag)

            if mismatch < best_mismatch:
                best_mismatch = mismatch
                best_tof2     = tof2
                best_vout     = v_inf_out_mars
                best_r_earth  = r_earth_t2

        if best_mismatch < MARS_FREE_THRESHOLD and best_tof2 is not None:
            # Compute departure v∞ to assess feasibility
            v_inf_dep = v_sc_out1 - v_dep
            survivors.append(dict(
                dv_mars_mismatch = best_mismatch,
                t0               = t0,
                tof1             = tof1,
                tof2             = best_tof2,
                v_inf_dep        = v_inf_dep,
                v_inf_in_mars    = v_inf_in_mars,
                v_inf_out_mars   = best_vout,
                r_dep            = r_dep,
                r_mars           = r_mars_t1,
                r_earth_ret      = best_r_earth,
            ))

print(f"\n  Stage 1 done — {len(survivors)} near-free-flyby windows found "
      f"(|Δv∞_mars| < {MARS_FREE_THRESHOLD:.0f} m/s)\n")

if not survivors:
    print("  No windows found — try increasing MARS_FREE_THRESHOLD or widening TOF ranges.")
    raise SystemExit

# ============================================================
# Stage 2 — Sweep TOF3 for each survivor, compute full Δv budget
# ============================================================

tof3_arr = np.arange(TOF3_MIN, TOF3_MAX, TOF3_STEP) * DAY

n_s2_combos = len(survivors) * len(tof3_arr)
print(f"  Stage 2 — Sweeping TOF3 for {len(survivors)} windows × "
      f"{len(tof3_arr)} TOF3 values = {n_s2_combos} combos")
print(f"            DSM grid search (6 tau × 8 phi × 3 delta = 144 evals per combo)\n", flush=True)

full_results = []   # (total_dv, details_dict)

for i_sv, sv in enumerate(survivors):
    t0   = sv['t0']
    tof1 = sv['tof1']
    tof2 = sv['tof2']
    t1   = t0 + tof1
    t2   = t1 + tof2

    r_dep       = sv['r_dep']
    r_mars      = sv['r_mars']
    r_earth_ret = sv['r_earth_ret']
    _, v_earth_t2 = fn_earth2(t2)

    # Mars flyby constraint
    k_in_mars = la.norm(sv['v_inf_in_mars'])
    ecc_m = 1.0 + rp_min_mars * k_in_mars**2 / mars.mu
    delta_max_mars = 2.0 * np.arcsin(1.0 / ecc_m)

    v_inf_dep_mag = la.norm(sv['v_inf_dep'])

    # Departure Δv (Oberth from LEO)
    r_leo = earth.radius + 500e3
    v_hyp_leo = np.sqrt(v_inf_dep_mag**2 + 2.0 * earth.mu / r_leo)
    v_circ_leo = np.sqrt(earth.mu / r_leo)
    dv1 = v_hyp_leo - v_circ_leo

    for tof3 in tof3_arr:
        t3 = t2 + tof3
        r_jup, v_jup = fn_jupiter(t3)

        try:
            v_sc_out_earth, v_sc_in3 = lambert_solver(r_earth_ret, r_jup, tof3, mu_sun)
        except Exception:
            continue

        v_inf_arr_jup = v_sc_in3 - v_jup

        # Jupiter capture Δv
        r_jup_park = jupiter.radius + 10_000e3
        v_hyp_jup = np.sqrt(la.norm(v_inf_arr_jup)**2 + 2.0 * jupiter.mu / r_jup_park)
        v_circ_jup = np.sqrt(jupiter.mu / r_jup_park)
        dv3 = v_hyp_jup - v_circ_jup

        # Mars velocity at flyby (already have r_mars from survivor; recompute v)
        _, v_mars_fly = fn_mars(t1)

        # DSM grid search: optimize (tau, phi, delta) over Mars→Earth leg
        dv_dsm, k_mis, tau_b, phi_b, df_b = dsm_grid_search(
            r_flyby       = r_mars,
            v_flyby       = v_mars_fly,
            v_inf_in      = sv['v_inf_in_mars'],
            k_in          = k_in_mars,
            delta_max_rad = delta_max_mars,
            r_target      = r_earth_ret,
            v_target      = v_earth_t2,
            tof_leg       = tof2,
            v_sc_out_target = v_sc_out_earth,
            mu_sun        = mu_sun,
            lam_earth     = 1.0,    # penalise Earth flyby mismatch equally
            n_tau=6, n_phi=8, n_delta=3,
        )

        total_dv = dv1 + dv_dsm + dv3

        full_results.append((total_dv, dict(
            total_dv       = total_dv,
            dv1            = dv1,
            dv_dsm         = dv_dsm,
            k_mis_earth    = k_mis,
            dv3            = dv3,
            dv_arr_vinf    = la.norm(v_inf_arr_jup),
            tau_dsm        = tau_b,
            phi_dsm        = phi_b,
            delta_frac_dsm = df_b,
            t0             = t0,
            tof1           = tof1,
            tof2           = tof2,
            tof3           = tof3,
            t_dep_abs      = (t_epoch + t0) / DAY,
            t_mars_abs     = (t_epoch + t1) / DAY,
            t_earth_abs    = (t_epoch + t2) / DAY,
            t_jup_abs      = (t_epoch + t3) / DAY,
            k_in_mars      = k_in_mars,
            dv_mars_mismatch = sv['dv_mars_mismatch'],
        )))

full_results.sort(key=lambda x: x[0])
full_results = full_results[:N_RESULTS]

# ============================================================
# Print ranked results
# ============================================================

print("=" * 82)
print(f"  Stage 2 — Top {len(full_results)} EMEJ windows (DSM-optimized, free Mars flyby)")
print("=" * 82)
print(f"  {'Rank':>4}  {'t_dep':>8}  {'TotalΔv':>9}  {'Δv1':>8}  "
      f"{'DSM':>8}  {'Δv3':>8}  {'EarthMis':>9}  {'v∞Jup':>8}")
print(f"  {'':>4}  {'(day)':>8}  {'(km/s)':>9}  {'(km/s)':>8}  "
      f"{'(km/s)':>8}  {'(km/s)':>8}  {'(m/s)':>9}  {'(km/s)':>8}")
print("  " + "-" * 80)

for rank, (total_dv, r) in enumerate(full_results, 1):
    print(f"  {rank:>4}  {r['t_dep_abs']:>8.1f}  {total_dv/1e3:>9.3f}  "
          f"{r['dv1']/1e3:>8.3f}  {r['dv_dsm']/1e3:>8.3f}  {r['dv3']/1e3:>8.3f}  "
          f"{r['k_mis_earth']:>9.1f}  {r['dv_arr_vinf']/1e3:>8.3f}")
    print(f"          TOF1={r['tof1']/DAY:.0f}d  TOF2={r['tof2']/DAY:.0f}d  TOF3={r['tof3']/DAY:.0f}d  "
          f"tau={r['tau_dsm']:.2f}  |Δv∞_Mars|={r['dv_mars_mismatch']:.0f} m/s")
    print(f"          t_mars={r['t_mars_abs']:.0f}d  t_earth={r['t_earth_abs']:.0f}d  "
          f"t_jup={r['t_jup_abs']:.0f}d")
    print()

# ============================================================
# Detailed breakdown of best result
# ============================================================

if full_results:
    best_dv, best = full_results[0]
    print("=" * 70)
    print("  Best Window — Detailed Breakdown (with DSM placement)")
    print("=" * 70)
    print(f"  Earth departure:  day {best['t_dep_abs']:.1f}")
    print(f"  Mars flyby:       day {best['t_mars_abs']:.1f}  (TOF1 = {best['tof1']/DAY:.0f} d)")
    tof2a = best['tau_dsm'] * best['tof2']
    tof2b = (1 - best['tau_dsm']) * best['tof2']
    print(f"  DSM:              day {best['t_mars_abs'] + tof2a/DAY:.1f}  (tau={best['tau_dsm']:.2f}, "
          f"TOF2a={tof2a/DAY:.0f}d, TOF2b={tof2b/DAY:.0f}d)")
    print(f"  Earth flyby:      day {best['t_earth_abs']:.1f}  (TOF2 = {best['tof2']/DAY:.0f} d)")
    print(f"  Jupiter arrival:  day {best['t_jup_abs']:.1f}  (TOF3 = {best['tof3']/DAY:.0f} d)")
    total_days = best['t_jup_abs'] - best['t_dep_abs']
    print(f"  Total mission:    {total_days:.0f} days ({total_days/365.25:.1f} yr)")
    print()
    print(f"  Δv1 (Earth dep):  {best['dv1']/1e3:.3f} km/s")
    print(f"  DSM Δv:           {best['dv_dsm']/1e3:.3f} km/s  (tau={best['tau_dsm']:.2f})")
    print(f"  Δv3 (Jup arr):    {best['dv3']/1e3:.3f} km/s")
    print(f"  TOTAL:            {best_dv/1e3:.3f} km/s")
    print()
    print(f"  Earth flyby v∞ mismatch: {best['k_mis_earth']:.1f} m/s  "
          f"({'free ✓' if best['k_mis_earth'] < 500 else 'powered'})")
    print(f"  Mars v∞_in:  {best['k_in_mars']/1e3:.3f} km/s")
