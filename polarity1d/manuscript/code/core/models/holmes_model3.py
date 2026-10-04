from __future__ import annotations
import numpy as np
from dataclasses import replace
from scipy.optimize import root
from polarity1d.contract import Grid1D
from polarity1d.holmes import HolmesModel3
from polarity1d.protocols import GaussianPulse, NoStimulus
from polarity1d.simulation import simulate
EQUATIONS_REVISION="original_PLOS_Eq5_Eq6_factor2_no_PI"
ACTIVATION=3.5
GRADIENT_AMPLITUDE=1.0
TRAINING_AMPLITUDE=0.2
TRAINING_DURATION=40.0
PREPARATION_TIME=2000.0

class HolmesModel3Eq56(HolmesModel3):
    """No-PI Model 3 with the factor 2 printed in original Eqs. (5) and (6).

    Table 1 does not give a unique Model 3 Rac basal activation value. The
    inherited 0.2 value is retained only as the declared low-end search point;
    the revision2 campaign labels every activation choice as tested.
    """

    name = "holmes_model3"

    def metadata(self):
        result = dict(super().metadata())
        source = dict(result["source"])
        source["equations"] = (
            "Original Holmes et al. (2012) Eqs. (1), (5), and (6), no PI feedback; "
            "Eq. (5) Rac activation has denominator 2*(1+(rho/a1)^n), and Eq. (6) "
            "Rho activation has denominator 2*(1+(R/a2)^n)."
        )
        source["parameter_table"] = (
            "Table 1 gives Model 4 Rac activation components and other shared values; "
            "the varied Model 3 basal Rac activation is declared as a search parameter, "
            "not presented as a uniquely published default."
        )
        result["source"] = source
        result["regime"] = {
            "name": "bounded Model 3 no-PI activation screen",
            "rac_activation_status": "candidate value, explicitly tested",
            "equations_revision": EQUATIONS_REVISION,
        }
        return result

    def reaction(self, time, x, state, stimulus):
        del time, x
        rac, rac_i, rho, rho_i = state
        rac_rate = self.rac_activation / (
            2.0 * (1.0 + (rho / self.rho_half_inhibition_of_rac) ** self.hill_n)
        ) + stimulus
        rho_rate = self.rho_activation / (
            2.0 * (1.0 + (rac / self.rac_half_inhibition_of_rho) ** self.hill_n)
        )
        rac_transfer = rac_rate * rac_i / self.rac_total - self.rac_decay * rac
        rho_transfer = rho_rate * rho_i / self.rho_total - self.rho_decay * rho
        return np.vstack((rac_transfer, -rac_transfer, rho_transfer, -rho_transfer))


def homogeneous_fixed_points(rac_activation: float) -> list[dict]:
    """Find distinct homogeneous equilibria and classify their linear stability."""
    model = replace(HolmesModel3Eq56(), rac_activation=float(rac_activation))
    grid = Grid1D(20.0, 5)
    x = grid.x

    def field(active):
        rac, rho = active
        state = np.array([[rac] * 5, [model.rac_total - rac] * 5,
                          [rho] * 5, [model.rho_total - rho] * 5])
        delta = model.reaction(0.0, x, state, np.zeros(5))
        return delta[[0, 2], 0]

    starts = [(r, h) for r in np.linspace(0, model.rac_total, 9)
              for h in np.linspace(0, model.rho_total, 9)]
    found: list[np.ndarray] = []
    for initial in starts:
        sol = root(field, initial)
        if sol.success and np.linalg.norm(field(sol.x), ord=np.inf) < 1e-8:
            if (-1e-8 <= sol.x[0] <= model.rac_total + 1e-8
                    and -1e-8 <= sol.x[1] <= model.rho_total + 1e-8
                    and not any(np.linalg.norm(sol.x - prev) < 1e-5 for prev in found)):
                found.append(sol.x)
    output = []
    for rac, rho in sorted(found, key=lambda z: z[0]):
        eps = 1e-5
        jac = np.column_stack([(field(np.array([rac, rho]) + np.eye(2)[i] * eps)
                                - field(np.array([rac, rho]) - np.eye(2)[i] * eps)) / (2 * eps)
                               for i in range(2)])
        eig = np.linalg.eigvals(jac)
        output.append({"rac_active": float(rac), "rho_active": float(rho),
                       "maximum_nullcline_residual": float(np.linalg.norm(field(np.array([rac, rho])), ord=np.inf)),
                       "eigenvalues": [complex(v).real if abs(complex(v).imag) < 1e-10
                                       else [complex(v).real, complex(v).imag] for v in eig],
                       "stable_homogeneous": bool(np.max(eig.real) < -1e-7)})
    return output


def initial_homogeneous_state(model: HolmesModel3, cells: int) -> np.ndarray:
    """Use the lowest stable no-cue homogeneous root for this activation rate."""
    roots = homogeneous_fixed_points(model.rac_activation)
    stable = [row for row in roots if row["stable_homogeneous"]]
    if not stable:
        raise RuntimeError(f"no stable homogeneous root at Rac activation {model.rac_activation}")
    root_row = min(stable, key=lambda row: row["rac_active"])
    rac, rho = root_row["rac_active"], root_row["rho_active"]
    return np.vstack((np.full(cells, rac), np.full(cells, model.rac_total - rac),
                      np.full(cells, rho), np.full(cells, model.rho_total - rho)))

def _sim(model, state, length, cells, protocol, end, interval=2.5, max_step=10.0):
    state = np.asarray(state, dtype=float)
    if np.min(state) < -1e-8:
        raise ValueError("restart state has a negative concentration")
    return simulate(model, initial_state=np.maximum(state, 0.0), length=float(length),
                    cells=int(cells), protocol=protocol, end_time=float(end),
                    output_interval=float(interval), max_step=float(max_step), backend="bdf")


def _orientation(tr):
    rac = np.maximum(tr.fields["rac_active"], 0.0)
    coord = 2.0 * tr.x / (tr.x[0] + tr.x[-1]) - 1.0
    mass = rac.sum(axis=1)
    return np.divide(rac @ coord, mass, out=np.zeros(len(mass)), where=mass > 1e-12)


def _contrast(tr):
    return np.ptp(np.maximum(tr.fields["rac_active"], 0.0), axis=1)


def _tail(tr, seconds=500.0):
    orientation = _orientation(tr)
    contrast = _contrast(tr)
    keep = tr.time >= tr.time[-1] - seconds
    return {"orientation_mean": float(np.mean(orientation[keep])),
            "orientation_min": float(np.min(orientation[keep])),
            "orientation_max": float(np.max(orientation[keep])),
            "orientation_range": float(np.ptp(orientation[keep])),
            "minimum_active_contrast": float(np.min(contrast[keep]))}


def _low_state(model, length, cells):
    return initial_homogeneous_state(model, cells)


def _trained_left_state(model, length, cells):
    initial = _low_state(model, length, cells)
    pulse = GaussianPulse(amplitude=TRAINING_AMPLITUDE, start=0.0,
                          stop=TRAINING_DURATION, width_fraction=0.1, centers=(0.15,))
    tr = _sim(model, initial, length, cells, pulse, TRAINING_DURATION, interval=2.5)
    state = np.vstack([tr.fields[name][-1] for name in model.state_names])
    free = _sim(model, state, length, cells, NoStimulus(), PREPARATION_TIME, interval=2.5)
    state = np.vstack([free.fields[name][-1] for name in model.state_names])
    return state, _tail(free), "40-s local Rac activation cue (0.2), then 2000-s no-cue settling"

