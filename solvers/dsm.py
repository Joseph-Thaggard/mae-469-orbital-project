"""
solvers/dsm.py — Deep Space Maneuver placement optimizer

Physics model:
  After a free gravity assist at a flyby body (|v∞_out| = |v∞_in|, turn ≤ δ_max),
  the spacecraft coasts on a Keplerian orbit for time τ·TOF_leg, then fires a DSM
  to redirect onto a Lambert arc connecting the DSM point to the next target body.

  Three free parameters per flyby–DSM pair:
    tau   ∈ (0, 1)  — fraction of total leg time at which DSM fires
    phi   ∈ [0, 2π) — direction angle for the flyby exit turn
    delta ∈ [0, 1]  — fraction of max turn used at flyby (1 = full δ_max)

  Objective: minimize DSM Δv + λ × |Earth flyby v∞ speed mismatch|
  Constraint: delta ≤ 1  (turn ≤ δ_max),  0 < tau < 1

Usage (from task2_emej.py):
    from solvers.dsm import optimize_dsm, compute_dsm_states

Usage (from scan):
    from solvers.dsm import dsm_grid_search
"""

import numpy as np
import numpy.linalg as la
from scipy.optimize import minimize

from solvers.solver import kepler_propagate, lambert_solver


# ---------------------------------------------------------------------------
# Build post-flyby v∞ vector
# ---------------------------------------------------------------------------

def build_v_inf_out(v_inf_in, k_in, delta_rad, phi):
    """
    Rotate v_inf_in by delta_rad toward direction phi in the plane
    perpendicular to v_inf_in.  Returns v∞_out with |v∞_out| = k_in.

    Args:
        v_inf_in:  incoming v∞ vector (m/s)
        k_in:      |v_inf_in| (m/s)
        delta_rad: turn angle (rad); must be ≤ δ_max
        phi:       rotation angle in perpendicular plane (rad)

    Returns:
        v_inf_out: outgoing v∞ vector (m/s), |v_inf_out| = k_in
    """
    v_in_unit = v_inf_in / k_in
    arb = np.array([1.0, 0.0, 0.0]) if abs(v_in_unit[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp1 = np.cross(v_in_unit, arb)
    perp1 /= la.norm(perp1)
    perp2 = np.cross(v_in_unit, perp1)
    turn_dir = np.cos(phi) * perp1 + np.sin(phi) * perp2
    return k_in * (np.cos(delta_rad) * v_in_unit + np.sin(delta_rad) * turn_dir)


# ---------------------------------------------------------------------------
# Single DSM evaluation
# ---------------------------------------------------------------------------

def eval_dsm(tau, phi, delta_frac,
             r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
             r_target, v_target, tof_leg,
             v_sc_out_target,
             mu_sun,
             mu_target=None, rp_min_target=None,
             return_states=False):
    """
    Evaluate DSM Δv and target-body flyby feasibility for given parameters.

    Args:
        tau:            DSM timing fraction (0 < tau < 1)
        phi:            flyby exit direction angle (rad)
        delta_frac:     fraction of δ_max used (0 ≤ delta_frac ≤ 1)
        r_flyby:        flyby body position at flyby epoch (m)
        v_flyby:        flyby body velocity at flyby epoch (m/s)
        v_inf_in:       incoming v∞ at flyby body (m/s)
        k_in:           |v_inf_in| (m/s)
        delta_max_rad:  max achievable turn angle at flyby body (rad)
        r_target:       next target body position at arrival (m)
        v_target:       next target body velocity at arrival (m/s)
        tof_leg:        total time from flyby to target body (s)
        v_sc_out_target: spacecraft velocity DEPARTING target body (m/s)
        mu_sun:         solar GM (m³/s²)
        mu_target:      GM of target flyby body (m³/s²); used for turn feasibility
        rp_min_target:  minimum periapsis distance at target body (m)
        return_states:  if True, return full state dict; else (dv_dsm, k_mis, turn_penalty)

    Returns:
        if return_states=False: (dv_dsm, k_mismatch, turn_penalty)
        if return_states=True:  dict with all states, or None on failure
    """
    if not (0.02 < tau < 0.98) or not (0.0 <= delta_frac <= 1.0):
        if return_states:
            return None
        return (1e10, 1e10, 1e10)

    delta_rad = delta_frac * delta_max_rad
    tof_a = tau * tof_leg
    tof_b = (1.0 - tau) * tof_leg

    # Post-flyby spacecraft velocity (free gravity assist)
    v_inf_out = build_v_inf_out(v_inf_in, k_in, delta_rad, phi)
    v_sc_dep_flyby = v_flyby + v_inf_out

    # Propagate to DSM point
    r_dsm, v_at_dsm = kepler_propagate(r_flyby, v_sc_dep_flyby, tof_a, mu_sun)

    # Lambert: DSM → target body
    try:
        v_dep_dsm, v_arr_target = lambert_solver(r_dsm, r_target, tof_b, mu_sun)
    except Exception:
        if return_states:
            return None
        return (1e10, 1e10, 1e10)

    dv_dsm = la.norm(v_dep_dsm - v_at_dsm)

    # Check free flyby at target body (speed mismatch)
    v_inf_in_target  = v_arr_target - v_target
    v_inf_out_target = v_sc_out_target - v_target
    k_in_target  = la.norm(v_inf_in_target)
    k_out_target = la.norm(v_inf_out_target)
    k_mismatch   = abs(k_in_target - k_out_target)

    # Turn feasibility penalty at target body
    turn_penalty = 0.0
    turn_deg = None
    max_turn_deg = None
    rp_required = None
    if mu_target is not None and rp_min_target is not None and k_in_target > 1e-3:
        ecc_tgt = 1.0 + rp_min_target * k_in_target**2 / mu_target
        max_turn_rad = 2.0 * np.arcsin(1.0 / ecc_tgt)
        cos_turn = np.clip(
            np.dot(v_inf_in_target, v_inf_out_target) / (k_in_target * k_out_target),
            -1.0, 1.0
        )
        turn_rad = np.arccos(cos_turn)
        turn_deg = np.degrees(turn_rad)
        max_turn_deg = np.degrees(max_turn_rad)
        sin_half = np.sin(turn_rad / 2.0)
        rp_required = mu_target / k_in_target**2 * (1.0 / sin_half - 1.0) if sin_half > 1e-9 else np.inf
        # Penalty: m/s per degree over the limit
        turn_penalty = max(0.0, turn_rad - max_turn_rad) * 1e4  # 10 km/s per radian over

    if not return_states:
        return dv_dsm, k_mismatch, turn_penalty

    return {
        'dv_dsm':           dv_dsm,
        'k_mismatch':       k_mismatch,
        'turn_penalty':     turn_penalty,
        'tau':              tau,
        'phi':              phi,
        'delta_frac':       delta_frac,
        'delta_rad':        delta_rad,
        'r_dsm':            r_dsm,
        'v_at_dsm':         v_at_dsm,
        'v_dep_dsm':        v_dep_dsm,
        'v_arr_target':     v_arr_target,
        'v_sc_dep_flyby':   v_sc_dep_flyby,
        'v_inf_out':        v_inf_out,
        'v_inf_in_target':  v_inf_in_target,
        'v_inf_out_target': v_inf_out_target,
        'k_in_target':      k_in_target,
        'k_out_target':     k_out_target,
        'tof_a':            tof_a,
        'tof_b':            tof_b,
        'turn_deg':         turn_deg,
        'max_turn_deg':     max_turn_deg,
        'rp_required':      rp_required,
    }


# ---------------------------------------------------------------------------
# Grid search (fast, used by scanner)
# ---------------------------------------------------------------------------

def dsm_grid_search(r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                    r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                    mu_target=None, rp_min_target=None,
                    lam_earth=1.0, lam_turn=1.0, n_tau=6, n_phi=8, n_delta=3):
    """
    Coarse grid search over (tau, phi, delta_frac).  Returns best (dv_dsm, k_mismatch)
    and the winning parameters.  Used in scan for speed.

    Cost = dv_dsm + lam_earth * k_mismatch + lam_turn * turn_penalty
    """
    best_cost   = np.inf
    best_params = (0.3, 0.0, 1.0)

    tau_vals   = np.linspace(0.05, 0.95, n_tau)
    phi_vals   = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    delta_vals = np.linspace(0.5, 1.0, n_delta)

    for tau in tau_vals:
        for phi in phi_vals:
            for df in delta_vals:
                dv, km, tp = eval_dsm(tau, phi, df,
                                      r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                                      r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                                      mu_target, rp_min_target)
                cost = dv + lam_earth * km + lam_turn * tp
                if cost < best_cost:
                    best_cost   = cost
                    best_params = (tau, phi, df)

    tau_opt, phi_opt, df_opt = best_params
    dv_opt, km_opt, _ = eval_dsm(tau_opt, phi_opt, df_opt,
                                  r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                                  r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                                  mu_target, rp_min_target)
    return dv_opt, km_opt, tau_opt, phi_opt, df_opt


# ---------------------------------------------------------------------------
# Full optimizer (grid + Nelder-Mead refinement, used by task2)
# ---------------------------------------------------------------------------

def optimize_dsm(r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                 r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                 mu_target=None, rp_min_target=None,
                 lam_earth=1.0, lam_turn=1.0,
                 n_tau=10, n_phi=12, n_delta=4,
                 refine=True):
    """
    Find (tau, phi, delta_frac) minimizing total cost including DSM Δv,
    Earth v∞ mismatch, and target flyby turn-angle feasibility.

    Args:
        mu_target:      GM of the target flyby body (m³/s²)
        rp_min_target:  minimum periapsis at target body (m)
        lam_earth:      weight for Earth v∞ speed mismatch penalty
        lam_turn:       weight for turn-angle infeasibility penalty
                        (cost += lam_turn * max(0, turn_rad - max_turn_rad) * 1e4)
        n_tau, n_phi, n_delta: grid resolution
        refine: polish with Nelder-Mead

    Returns:
        dict from eval_dsm with return_states=True, plus 'cost' key.
        None if no feasible solution found.
    """
    def _cost(params):
        tau, phi, df = params
        dv, km, tp = eval_dsm(tau, phi, df,
                               r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                               r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                               mu_target, rp_min_target)
        return dv + lam_earth * km + lam_turn * tp

    # --- Grid search ---
    best_cost   = np.inf
    best_params = (0.3, 0.0, 1.0)

    for tau in np.linspace(0.05, 0.95, n_tau):
        for phi in np.linspace(0, 2 * np.pi, n_phi, endpoint=False):
            for df in np.linspace(0.3, 1.0, n_delta):
                c = _cost([tau, phi, df])
                if c < best_cost:
                    best_cost   = c
                    best_params = (tau, phi, df)

    # --- Nelder-Mead refinement ---
    if refine:
        res = minimize(_cost, list(best_params),
                       method='Nelder-Mead',
                       options={'xatol': 0.003, 'fatol': 5.0, 'maxiter': 500})
        if res.fun < best_cost:
            best_cost   = res.fun
            best_params = tuple(res.x)

    tau_opt, phi_opt, df_opt = best_params
    states = eval_dsm(tau_opt, phi_opt, df_opt,
                      r_flyby, v_flyby, v_inf_in, k_in, delta_max_rad,
                      r_target, v_target, tof_leg, v_sc_out_target, mu_sun,
                      mu_target, rp_min_target, return_states=True)
    if states is not None:
        states['cost'] = best_cost
    return states
