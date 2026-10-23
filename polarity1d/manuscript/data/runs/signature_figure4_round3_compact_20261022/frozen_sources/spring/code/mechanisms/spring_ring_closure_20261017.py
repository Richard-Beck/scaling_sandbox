"""Proposed Spring ring closure: algebraic checks only, no trajectory runs.

The local kinetic/force laws are imported from the active native adapter.
The rotationally invariant mean-force circumference closure is NEW; it is
specified for review in SPRING_RING_CLOSURE_20261017.md, not a native result.
Running this file performs operator checks and writes their small JSON record.
"""
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import sys

sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/core"))
from models import spring as native


@dataclass(frozen=True)
class Closure:
    resting_circumference: float = 2.0
    cells: int = 200

    def __post_init__(self):
        if self.resting_circumference <= 0 or self.cells < 8:
            raise ValueError("Need positive rest circumference and at least eight cells")

    @property
    def material_positions(self):
        return (np.arange(self.cells) + .5) / self.cells

    def initial(self):
        # Fresh uniform source concentrations; no interval checkpoint or cosine.
        return np.r_[np.full(self.cells, .05645),
                     np.full(self.cells, 2.0 - .05645),
                     self.resting_circumference]

    def fields(self, state):
        state = np.asarray(state, dtype=float)
        if state.shape != (2 * self.cells + 1,) or not np.isfinite(state).all():
            raise ValueError("Unexpected or nonfinite state")
        if state[-1] <= 0:
            raise ValueError("Circumference must be positive")
        return state[:self.cells], state[self.cells:-1], float(state[-1])

    def rhs(self, state, cue):
        u, v, circumference = self.fields(state)
        signal = np.broadcast_to(np.asarray(cue, dtype=float), (self.cells,))
        extension = circumference - self.resting_circumference
        # NEW geometry closure, preserving the uniform native length equation.
        speed = 2.0 * (np.mean(native.auth["F"](u)) - .01 * extension)
        reaction = native.auth["f"](u, v, 4.0 + signal, 5, 6, 2,
                                    3.0 + 80.0 * extension)
        lap_u = np.roll(u, 1) - 2.0 * u + np.roll(u, -1)
        lap_v = np.roll(v, 1) - 2.0 * v + np.roll(v, -1)
        inverse_dx2 = (self.cells / circumference) ** 2
        return np.r_[reaction + .01 * lap_u * inverse_dx2 - u * speed / circumference,
                     -reaction + 10.0 * lap_v * inverse_dx2 - v * speed / circumference,
                     speed]

    def amount(self, state):
        u, v, circumference = self.fields(state)
        return float(circumference * np.mean(u + v))


def operator_checks():
    closure = Closure()
    n = closure.cells
    rng = np.random.default_rng(20261017)
    state = np.r_[rng.uniform(.1, 1.7, n), rng.uniform(.3, 1.2, n), 2.045]
    cue = rng.uniform(0, 4, n)
    derivative = closure.rhs(state, cue)
    u, v, circumference = closure.fields(state)
    d_amount = circumference * np.mean(derivative[:n] + derivative[n:-1])
    d_amount += derivative[-1] * np.mean(u + v)
    assert abs(d_amount) < 1e-10
    shift = 37
    rotated = np.r_[np.roll(u, shift), np.roll(v, shift), circumference]
    rotated_derivative = closure.rhs(rotated, np.roll(cue, shift))
    expected_rotated = np.r_[np.roll(derivative[:n], shift),
                             np.roll(derivative[n:-1], shift), derivative[-1]]
    rotation_residual = np.max(np.abs(rotated_derivative - expected_rotated))
    assert rotation_residual < 1e-10
    reflection = np.r_[u[::-1], v[::-1], circumference]
    reflected_derivative = closure.rhs(reflection, cue[::-1])
    expected_reflected = np.r_[derivative[:n][::-1], derivative[n:-1][::-1], derivative[-1]]
    reflection_residual = np.max(np.abs(reflected_derivative - expected_reflected))
    assert reflection_residual < 1e-10
    # Exact agreement of the complete derivative on homogeneous fields. This
    # checks only this special case, not equivalence of nonuniform mechanics.
    uniform = np.r_[np.full(n, 1.15), np.full(n, .72), 2.045]
    periodic_derivative = closure.rhs(uniform, 0.0)
    native_state = np.empty(2 * n + 2)
    native_state[:-2:2], native_state[1:-2:2] = uniform[:n], uniform[n:-1]
    native_state[-2:] = 0., uniform[-1]
    native_derivative = native.rhs(0., native_state, L0=closure.resting_circumference)
    expected_native = np.r_[native_derivative[:-2:2], native_derivative[1:-2:2],
                            native_derivative[-1] - native_derivative[-2]]
    uniform_residual = np.max(np.abs(periodic_derivative - expected_native))
    assert uniform_residual < 1e-12
    # The closure averages the force after the nonlinear Rac switch.
    stepped = np.r_[np.full(n // 4, 1.5), np.full(n - n // 4, .25)]
    mean_force = float(np.mean(native.auth["F"](stepped)))
    force_at_mean = float(native.auth["F"](np.mean(stepped)))
    assert abs(mean_force - force_at_mean) > 1e-4
    initial = closure.initial()
    assert closure.amount(initial) == 2 * closure.resting_circumference
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "Operator checks pass; proposed closure remains subject to review",
        "trajectory_integrations": 0,
        "native_source": "code/core/models/spring.py",
        "native_source_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "proposal_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": {"physical_amount_derivative_residual": float(d_amount),
                   "cyclic_rotation_covariance_residual": float(rotation_residual),
                   "reflection_covariance_residual": float(reflection_residual),
                   "uniform_native_full_derivative_residual": float(uniform_residual),
                   "mean_force_after_switch": mean_force,
                   "force_at_mean_activity_not_used": force_at_mean,
                   "fresh_initial_amount": closure.amount(initial)},
        "not_established": ["existence or stability of a polarized ring state",
                            "equivalence to native spatial mechanics",
                            "any challenger-distance effect"]}


if __name__ == "__main__":
    result = operator_checks()
    output = ROOT / "data/runs/spring_ring_closure_review_20261017/operator_checks.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
