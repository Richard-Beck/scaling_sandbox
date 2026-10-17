"""Exploratory mean-field slip clutches anchored to a fixed laboratory grid.

Traction is force per reference length in the cortex equation. Site density
is per physical substrate length; quadrature accounts for this conversion.
"""
from dataclasses import dataclass, asdict

import numpy as np


@dataclass(frozen=True)
class FixedSubstrateParameters:
    density: float = 10.0  # sites / micrometre of substrate
    spring: float = 1.0  # pN / micrometre per bound clutch
    on_rate: float = 0.1  # / s
    off_rate: float = 0.05  # / s at zero force
    rupture_force: float = 0.5  # pN
    spacing: float = 0.025  # micrometres; numerical substrate lattice

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in asdict(self).values()):
            raise ValueError("Clutch parameters must be finite and positive")


class FixedSubstrateClutch:
    """Site states never move with the cell; only the contact footprint moves.

The off-rate is evaluated from the previous bound-mean extension. Occupancy
and the extension moment use backward Euler, with the new slip velocity
coupled implicitly to cortical velocity and rigid translation.
"""
    def __init__(self, parameters=FixedSubstrateParameters(), keep_history=False):
        self.parameters = parameters
        self.keep_history = keep_history
        self.ids = np.empty(0, dtype=int)
        self.n = np.empty(0)
        self.q = np.empty(0)
        self.translation = 0.0
        self.speed = 0.0
        self.max_force_residual = 0.0
        self.fields = {}
        self.history = []

    def _contact(self, edges):
        h = self.parameters.spacing
        first = int(np.ceil(edges[0] / h - 0.5))
        last = int(np.floor(edges[-1] / h - 0.5))
        ids = np.arange(first, last + 1)
        positions = (ids + 0.5) * h
        if len(ids) < 3:
            raise ValueError("Substrate lattice is too coarse for this cell")
        # Control intervals tile the footprint exactly; interiors have width h.
        weights = np.diff(np.r_[edges[0], 0.5 * (positions[1:] + positions[:-1]), edges[-1]])
        bins = np.clip(np.searchsorted(edges, positions, side="right") - 1, 0, len(edges)-2)
        if np.unique(bins).size != len(edges)-1:
            raise ValueError("Each cell bin must contain substrate sites; reduce spacing")
        n, q = np.zeros(len(ids)), np.zeros(len(ids))
        if self.ids.size:
            old = ids - self.ids[0]
            valid = (old >= 0) & (old < len(self.ids))
            n[valid], q[valid] = self.n[old[valid]], self.q[old[valid]]
        return ids, weights, bins, n, q

    def initialize(self, edges):
        self.ids, _, _, self.n, self.q = self._contact(edges)
        p = self.parameters
        self.n[:] = p.on_rate / (p.on_rate + p.off_rate)
        cells = len(edges)-1
        self.fields = {name: np.zeros(cells) for name in (
            "clutch_bound", "clutch_extension", "clutch_traction", "clutch_traction_lab_density")}
        self.fields["clutch_bound"][:] = self.n[0]

    def _prepare(self, edges, dx, dt):
        p = self.parameters
        ids, weights, bins, n_old, q_old = self._contact(edges)
        extension = np.divide(q_old, n_old, out=np.zeros_like(q_old), where=n_old > 1e-14)
        # Fail explicitly rather than silently capping a physically divergent rate.
        exponent = np.abs(p.spring * extension) / p.rupture_force
        if np.max(exponent) > 600:
            raise FloatingPointError("Clutch rupture rate overflow; reduce timestep")
        off = p.off_rate * np.exp(exponent)
        n_new = (n_old + dt * p.on_rate) / (1 + dt * (p.on_rate + off))
        q_memory = q_old / (1 + dt * off)
        q_slip = dt * n_new / (1 + dt * off)
        cells = len(edges)-1
        def deposit(value):
            return np.bincount(bins, weights=weights*value, minlength=cells) / dx
        memory = deposit(p.density * p.spring * q_memory)
        drag = deposit(p.density * p.spring * q_slip)
        return ids, weights, bins, n_new, q_memory, q_slip, memory, drag

    def _commit(self, prepared, velocity, speed, dx, dt):
        ids, weights, bins, n_new, q_memory, q_slip, memory, drag = prepared
        p = self.parameters
        cells = len(velocity)
        def deposit(value):
            return np.bincount(bins, weights=weights*value, minlength=cells) / dx
        self.speed = float(speed)
        self.ids, self.n = ids, n_new
        self.q = q_memory + q_slip * (velocity[bins] + self.speed)
        traction = deposit(p.density * p.spring * self.q)
        residual = abs(float(np.sum(traction) * dx))
        self.max_force_residual = max(self.max_force_residual, residual)
        if residual > 1e-8 * max(1.0, float(np.sum(np.abs(traction)) * dx)):
            raise FloatingPointError("Net substrate force did not balance")
        coverage = deposit(np.ones(len(ids)))
        occupancy = deposit(self.n) / coverage
        moment = deposit(self.q) / coverage
        self.fields = {
            "clutch_bound": occupancy,
            "clutch_extension": np.divide(moment, occupancy, out=np.zeros(cells), where=occupancy > 1e-14),
            "clutch_traction": traction,
            "clutch_traction_lab_density": traction / coverage,
        }
        self.translation += dt * self.speed
        return velocity

    def advance(self, *, edges, dx, dt, diagonal, viscosity, rhs, boundary, solve):
        prepared = self._prepare(edges, dx, dt)
        memory, drag = prepared[-2:]
        base = solve(diagonal + drag, viscosity, rhs - memory, dx, *boundary)
        response = solve(diagonal + drag, viscosity, -drag, dx, 0.0, 0.0)
        denominator = np.sum(drag * (1 + response))
        if denominator <= 0 or not np.isfinite(denominator):
            raise FloatingPointError("Singular free-translation clutch balance")
        self.speed = -float(np.sum(memory + drag * base) / denominator)
        velocity = base + response * self.speed
        return self._commit(prepared, velocity, self.speed, dx, dt)

    def advance_passive_edge(self, *, edges, dx, dt, **mechanics):
        from .passive_edge import solve_passive_edge
        prepared = self._prepare(edges, dx, dt)
        memory, drag = prepared[-2:]
        result = solve_passive_edge(dx=dx, dt=dt, memory=memory, drag=drag, **mechanics)
        u, velocity, speed, boundary_u, boundary_v, residuals = result
        self._commit(prepared, velocity, speed, dx, dt)
        return u, velocity, boundary_u, boundary_v, residuals

    def advance_local_rac_edges(self, *, edges, dx, dt, **mechanics):
        from .passive_edge import solve_local_rac_edges
        prepared = self._prepare(edges, dx, dt)
        memory, drag = prepared[-2:]
        u, velocity, speed, boundary_u, boundary_v, residuals = solve_local_rac_edges(
            dx=dx, dt=dt, memory=memory, drag=drag, **mechanics)
        self._commit(prepared, velocity, speed, dx, dt)
        return u, velocity, boundary_u, boundary_v, residuals

    def record(self, time):
        if self.keep_history:
            self.history.append((float(time), self.ids.copy(), self.n.copy(), self.q.copy()))

    def save_history(self, path):
        first = min(ids[0] for _, ids, _, _ in self.history)
        last = max(ids[-1] for _, ids, _, _ in self.history)
        n = np.full((len(self.history), last-first+1), np.nan)
        q = np.full_like(n, np.nan)
        for row, (_, ids, bound, moment) in enumerate(self.history):
            n[row, ids-first], q[row, ids-first] = bound, moment
        np.savez_compressed(path, time=np.array([r[0] for r in self.history]),
                            substrate_x=(np.arange(first,last+1)+0.5)*self.parameters.spacing,
                            bound_fraction=n, extension_moment=q)
