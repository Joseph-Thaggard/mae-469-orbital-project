"""
events.py — Simulation event system

Events represent conditions that trigger actions (burns, logging, SOI transitions)
at precise simulation times or when state-based conditions are met.

Event types:
    TimedBurnEvent      — fires at a specific t_sim, splits dt for precision
    DistanceBurnEvent   — fires when spacecraft enters a distance threshold from a body
    SOIEvent            — fires on SOI entry or exit (wraps physics.soi.check_soi)

EventSchedule holds all registered events and drives the step execution.

Usage:
    schedule = EventSchedule()
    schedule.add(TimedBurnEvent(ship, np.array([0, 1000, 0]), t_trigger=3e7))
    schedule.add(DistanceBurnEvent(ship, np.array([0, 500, 0]), earth, distance=1e8))

    # In simulation loop:
    soi_map = schedule.execute_step(bodies_list, spacecraft_list, grid, t_sim, dt,
                                    n_substeps=10000, render_every=59,
                                    soi_render_fn=v.render_body_frame)
"""

import numpy as np
import numpy.linalg as la
from abc import ABC, abstractmethod
from mission.burns import instant_burn


# ---------------------------------------------------------------------------
# Event base class
# ---------------------------------------------------------------------------

class Event(ABC):
    """Base class for all simulation events. Subclass and override is_triggered and execute."""

    def __init__(self, spacecraft):
        self.spacecraft = spacecraft
        self.fired = False          # prevents repeated firing after trigger

    @abstractmethod
    def is_triggered(self, bodies, t_sim, dt):
        """Return True if this event should fire during [t_sim, t_sim+dt)."""

    def t_exact(self, t_sim):
        """Return the exact time within the current step that the event fires.
        For time-based events this is precise. For condition-based events, returns t_sim
        (fire immediately at step start since condition is already met)."""
        return t_sim

    @abstractmethod
    def execute(self, bodies, t_sim):
        """Apply the event's action. Called at t_exact within the step."""

    def __repr__(self):
        return f"{self.__class__.__name__}(sc={self.spacecraft.name}, fired={self.fired})"


# ---------------------------------------------------------------------------
# Timed burn — fires at a specific simulation time
# ---------------------------------------------------------------------------

class TimedBurnEvent(Event):
    """Apply an instantaneous dv_vector to spacecraft at a specific simulation time.

    The step executor splits dt so that propagation stops exactly at t_trigger,
    applies the burn, then resumes for the remainder of dt.

    Args:
        spacecraft:  Spacecraft object to burn
        dv_vector:   np.array([dvx, dvy, dvz]) in m/s (inertial frame)
        t_trigger:   simulation time in seconds at which the burn fires
    """

    def __init__(self, spacecraft, dv_vector, t_trigger):
        super().__init__(spacecraft)
        self.dv_vector = np.array(dv_vector, dtype=float)
        self.t_trigger = float(t_trigger)

    def is_triggered(self, _bodies, t_sim, dt):
        return not self.fired and t_sim <= self.t_trigger < t_sim + dt

    def t_exact(self, t_sim):
        return self.t_trigger

    def execute(self, _bodies, t_sim):
        instant_burn(self.spacecraft, self.dv_vector)
        self.fired = True
        print(f"[Event] TimedBurn on {self.spacecraft.name} at t={t_sim:.2e} s: Δv={self.dv_vector} m/s")


# ---------------------------------------------------------------------------
# Distance burn — fires when spacecraft is within a distance of a body
# ---------------------------------------------------------------------------

class DistanceBurnEvent(Event):
    """Apply an instantaneous dv_vector when spacecraft comes within `distance` meters of `body`.

    Triggers once (self.fired prevents repeat). Fires at the start of the step
    when the condition is first detected — no sub-dt precision (distance crossing
    is not predictable without root-finding).

    Args:
        spacecraft:  Spacecraft object to burn
        dv_vector:   np.array([dvx, dvy, dvz]) in m/s
        body:        Body object to measure distance from
        distance:    trigger distance in meters (fires when |sc - body| <= distance)
    """

    def __init__(self, spacecraft, dv_vector, body, distance):
        super().__init__(spacecraft)
        self.dv_vector = np.array(dv_vector, dtype=float)
        self.body = body
        self.distance = float(distance)

    def is_triggered(self, bodies, t_sim, dt):
        if self.fired:
            return False
        d = la.norm(self.spacecraft.position - self.body.position)
        return d <= self.distance

    def execute(self, bodies, t_sim):
        d = la.norm(self.spacecraft.position - self.body.position)
        instant_burn(self.spacecraft, self.dv_vector)
        self.fired = True
        print(f"[Event] DistanceBurn on {self.spacecraft.name} at t={t_sim:.2e} s "
              f"(d={d:.3e} m from {self.body.name}): Δv={self.dv_vector} m/s")


# ---------------------------------------------------------------------------
# SOI event — fires when spacecraft enters or exits a body's SOI
# ---------------------------------------------------------------------------

class SOIEvent(Event):
    """Fire a callback when a spacecraft enters or exits a body's sphere of influence.

    Does not apply a burn — calls `on_enter` or `on_exit` callables instead,
    so the user can trigger logging, burns, or any other action.

    Args:
        spacecraft:  Spacecraft object to monitor
        body:        Body whose SOI is being monitored
        on_enter:    callable(spacecraft, body, t_sim) called on SOI entry  (optional)
        on_exit:     callable(spacecraft, body, t_sim) called on SOI exit   (optional)
    """

    def __init__(self, spacecraft, body, on_enter=None, on_exit=None):
        super().__init__(spacecraft)
        self.body     = body
        self.on_enter = on_enter
        self.on_exit  = on_exit
        self._was_inside = None     # None = not yet evaluated
        self.fired = False          # SOI events re-fire on each transition, not once-only

    def is_triggered(self, bodies, t_sim, dt):
        from physics.soi import check_soi
        # SOI events use check_soi logic but are evaluated every step
        return True     # always check — execute() decides whether a transition occurred

    def execute(self, bodies, t_sim):
        from physics.soi import check_soi
        # Recompute SOI membership using the full soi machinery
        central   = bodies[0]
        mu_c      = central.mu
        r         = self.body.position - central.position
        r_mag     = la.norm(r)
        v_mag     = la.norm(self.body.velocity)
        a         = 1.0 / (2.0/r_mag - v_mag**2/mu_c)
        if a <= 0:
            return
        r_soi    = a * (self.body.mass / central.mass) ** (2.0/5.0)
        d        = la.norm(self.spacecraft.position - self.body.position)
        inside   = d <= r_soi

        if self._was_inside is None:
            self._was_inside = inside
            return  # no transition on first evaluation

        if inside and not self._was_inside:
            print(f"[Event] {self.spacecraft.name} entered SOI of {self.body.name} at t={t_sim:.2e} s")
            if self.on_enter:
                self.on_enter(self.spacecraft, self.body, t_sim)
        elif not inside and self._was_inside:
            print(f"[Event] {self.spacecraft.name} exited SOI of {self.body.name} at t={t_sim:.2e} s")
            if self.on_exit:
                self.on_exit(self.spacecraft, self.body, t_sim)

        self._was_inside = inside


# ---------------------------------------------------------------------------
# EventSchedule — holds all events and drives step execution
# ---------------------------------------------------------------------------

class EventSchedule:
    """Container for simulation events. Drives precision sub-stepping around timed events.

    Usage:
        schedule = EventSchedule()
        schedule.add(TimedBurnEvent(ship, dv, t_trigger=3e7))
        schedule.add(DistanceBurnEvent(ship, dv, earth, distance=9e8))

        # in loop:
        soi_map = schedule.execute_step(bodies, spacecraft_list, grid, t_sim, dt, ...)
    """

    def __init__(self):
        self.events = []

    def add(self, event):
        """Register an event with the schedule."""
        self.events.append(event)

    def _triggered_this_step(self, bodies, t_sim, dt):
        """Return list of events that trigger during [t_sim, t_sim+dt), sorted by t_exact."""
        triggered = [e for e in self.events if e.is_triggered(bodies, t_sim, dt)]
        return sorted(triggered, key=lambda e: e.t_exact(t_sim))

    def execute_step(self, bodies, spacecraft_list, grid, t_sim, dt,
                     n_substeps=10000, render_every=59, soi_render_fn=None):
        """Execute one full dt step, splitting at timed event boundaries for precision.

        For each TimedBurnEvent in [t_sim, t_sim+dt):
            1. Propagate from t_current to t_burn
            2. Fire the burn
            3. Continue to next event or end of dt

        Condition-based events (distance, SOI) are checked and fired at step start,
        before any propagation, since their exact crossing time is not predicted.

        Args:
            bodies:           full bodies list
            spacecraft_list:  list of Spacecraft objects for SOI tracking
            grid:             SpaceGrid
            t_sim:            current simulation time (start of this step)
            dt:               outer timestep in seconds
            n_substeps:       sub-steps inside SOI (passed to propagate_soi)
            render_every:     SOI render interval (passed to propagate_soi)
            soi_render_fn:    callable(focus_body, bodies, radius) for SOI plots

        Returns:
            soi_map from the final propagate_soi call
        """
        from physics.soi import propagate_soi

        triggered = self._triggered_this_step(bodies, t_sim, dt)

        # Fire condition-based events first (no sub-dt precision needed)
        condition_events = [e for e in triggered if not isinstance(e, TimedBurnEvent)]
        for event in condition_events:
            event.execute(bodies, t_sim)

        # Split dt at each timed event boundary
        timed_events = [e for e in triggered if isinstance(e, TimedBurnEvent)]

        t_current = t_sim
        soi_map   = {}

        for event in timed_events:
            dt_to_event = event.t_exact(t_sim) - t_current
            if dt_to_event > 0:
                soi_map = propagate_soi(bodies, spacecraft_list, grid, dt_to_event,
                                        n_substeps=max(1, int(n_substeps * dt_to_event / dt)),
                                        render_every=render_every,
                                        soi_render_fn=soi_render_fn)
            event.execute(bodies, event.t_exact(t_sim))
            t_current = event.t_exact(t_sim)

        # Propagate remaining time after last event (or full dt if no timed events)
        dt_remaining = (t_sim + dt) - t_current
        if dt_remaining > 0:
            soi_map = propagate_soi(bodies, spacecraft_list, grid, dt_remaining,
                                    n_substeps=max(1, int(n_substeps * dt_remaining / dt)),
                                    render_every=render_every,
                                    soi_render_fn=soi_render_fn)

        return soi_map
