"""Development ring adapter for three native reaction–diffusion models.

The CLI runs algebraic/operator checks only. It does not integrate trajectories,
prepare a polarized state, submit SLURM jobs, or read endpoint checkpoints.
Future integration must follow review of the exported minimal pilot plan.
"""
from pathlib import Path
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "code/core"), str(ROOT / "code/core/src")]
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
from models.holmes_model3 import HolmesModel3Eq56, initial_homogeneous_state, ACTIVATION
from polarity1d.contract import Grid1D

OUT = ROOT / "figures/length_mechanisms_20261008/dev/ring_distance_prototype_20261017"
KEYS = ("wave_pinning", "goryachev", "holmes_model3")


@dataclass(frozen=True)
class RingGrid:
    circumference: float = 40.0
    cells: int = 400

    def __post_init__(self):
        if not np.isfinite(self.circumference) or self.circumference <= 0:
            raise ValueError("Circumference must be finite and positive")
        if isinstance(self.cells, bool) or not isinstance(self.cells, int) or self.cells < 8:
            raise ValueError("Need an integer number of at least eight ring cells")

    @property
    def dx(self):
        return self.circumference / self.cells

    @property
    def x(self):
        return (np.arange(self.cells) + 0.5) * self.dx

    def arc(self, x, center):
        return (np.asarray(x) - center + self.circumference / 2) % self.circumference - self.circumference / 2

    def distance(self, first, second):
        return np.abs(self.arc(first, second))


def native_model(key):
    if key == "wave_pinning":
        return WavePinning()
    if key == "goryachev":
        return Goryachev()
    if key == "holmes_model3":
        # The active manuscript native.model_for uses this exact regime.
        return replace(HolmesModel3Eq56(), rac_activation=ACTIVATION)
    raise ValueError(f"Not one of the three scoped models: {key}")


def native_rest(key, grid):
    model = native_model(key)
    if key == "holmes_model3":
        return initial_homogeneous_state(model, grid.cells)
    return model.initial_state(Grid1D(grid.circumference, grid.cells))


def periodic_laplacian(state, dx):
    """Cell-average divergence of forward face gradients, including seam face."""
    gradients = (np.roll(state, -1, axis=-1) - state) / dx
    return (gradients - np.roll(gradients, 1, axis=-1)) / dx


def periodic_sparsity(species, cells):
    pattern = sparse.lil_matrix((species * cells, species * cells), dtype=np.int8)
    indices = np.arange(cells)
    for field in range(species):
        rows = field * cells + indices
        for other in range(species):
            pattern[rows, other * cells + indices] = 1
        pattern[rows, field * cells + (indices - 1) % cells] = 1
        pattern[rows, field * cells + (indices + 1) % cells] = 1
    return pattern.tocsr()


def ring_cue(grid, center_cell, total_width, peak):
    """Symmetric raised cosine, native external activation-rate input.

    Discrete rotations preserve the entire cue array exactly, including its
    integrated rate and support. Width is total support, not a half-width.
    Basal external input is zero for all three native models.
    """
    if isinstance(center_cell, bool) or not isinstance(center_cell, (int, np.integer)):
        raise ValueError("Cue centre must be an integer cell offset")
    if not (0 < total_width < grid.circumference) or not np.isfinite(peak) or peak < 0:
        raise ValueError("Need positive sub-circumference width and nonnegative peak")
    distance = grid.distance(grid.x, grid.x[0])
    base = np.where(distance < total_width / 2,
                    peak * 0.5 * (1 + np.cos(2 * np.pi * distance / total_width)), 0.0)
    return np.roll(base, int(center_cell) % grid.cells)


def native_rhs(model, grid, state, stimulus):
    state = np.asarray(state, dtype=float)
    stimulus = np.asarray(stimulus, dtype=float)
    if state.shape != (len(model.state_names), grid.cells) or stimulus.shape != (grid.cells,):
        raise ValueError("State or stimulus does not match ring grid")
    return (np.asarray(model.diffusivities)[:, None] * periodic_laplacian(state, grid.dx)
            + model.reaction(0.0, grid.x, state, stimulus))


def pool_amounts(model, grid, state):
    return np.array([np.sum(np.asarray(state)[list(group)]) * grid.dx for group in model.conserved_groups])


def advance_constant_input(model, grid, initial, stimulus, end_time, sample_times,
                           rtol=1e-9, atol=1e-11, max_step=1.0):
    """Future SLURM phase helper; NOT called by the checks or CLI.

    Each phase has constant native external input. Withdrawal is a separate
    solve starting from the exact last state, avoiding a stepped discontinuity.
    No clipping or mass renormalization is applied. Store the complete fields.
    """
    initial = np.asarray(initial, dtype=float)
    if not np.isfinite(initial).all() or np.min(initial) < 0:
        raise ValueError("Need a finite nonnegative initial state")
    times = np.asarray(sample_times, dtype=float)
    if times[0] != 0 or times[-1] != end_time or np.any(np.diff(times) <= 0):
        raise ValueError("Samples must increase from zero through the exact phase end")
    signal = np.asarray(stimulus, dtype=float).copy()
    species = len(model.state_names)
    solution = solve_ivp(lambda t, y: native_rhs(model, grid, y.reshape(species, grid.cells), signal).ravel(),
                         (0.0, float(end_time)), initial.ravel(), method="BDF", t_eval=times,
                         rtol=rtol, atol=atol, max_step=max_step,
                         jac_sparsity=periodic_sparsity(species, grid.cells))
    if not solution.success:
        raise RuntimeError(solution.message)
    states = solution.y.T.reshape(len(times), species, grid.cells)
    if np.min(states) < -50 * atol:
        raise FloatingPointError("Materially negative periodic solution")
    initial_pools = pool_amounts(model, grid, initial)
    errors = np.array([pool_amounts(model, grid, state) / initial_pools - 1 for state in states])
    if np.max(np.abs(errors)) > 1e-7:
        raise FloatingPointError("Conserved ring pool drift exceeds tolerance")
    return dict(time=times, states=states, pool_relative_error=errors,
                nfev=int(solution.nfev), njev=int(solution.njev), nlu=int(solution.nlu))


def circular_components(mask):
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(len(mask))]
    starts = np.flatnonzero(mask & ~np.roll(mask, 1))
    components = []
    for start in starts:
        component = []
        index = int(start)
        while mask[index]:
            component.append(index)
            index = (index + 1) % len(mask)
        components.append(np.asarray(component))
    return components


def basal_core(primary):
    primary = np.asarray(primary)
    if np.ptp(primary) <= 1e-12:
        raise ValueError("Homogeneous state does not provide a single-site basal core")
    level = float(primary.min() + 0.5 * np.ptp(primary))
    components = circular_components(primary >= level)
    if len(components) != 1:
        raise ValueError("Fresh basal preparation must have exactly one connected active site")
    mask = np.zeros(len(primary), dtype=bool)
    mask[components[0]] = True
    return level, mask


def valid_challenger_locations(grid, primary_basal, offset_cells, width, peak):
    _, core = basal_core(primary_basal)
    incumbent = int(np.argmax(primary_basal))
    records = []
    for offset in offset_cells:
        center = (incumbent + int(offset)) % grid.cells
        cue = ring_cue(grid, center, width, peak)
        overlap = bool(np.any(core & (cue > 0)))
        records.append(dict(offset_cells=int(offset), center_cell=int(center), shortest_arc=float(grid.distance(grid.x[center], grid.x[incumbent])),
                            overlaps_original_core=overlap, eligible=not overlap))
    return records


def translation_compare(current, reference, field_scales):
    """Full-field integer cyclic alignment; diagnostic only, no state remap.

    Compare the unshifted and optimally aligned residuals while recording the
    shift. Low aligned error alone does not qualify a fixed-frame settled
    branch, and this is not a measured translational eigenvalue.
    """
    current, reference = np.asarray(current), np.asarray(reference)
    scales = np.asarray(field_scales)[:, None]
    if current.shape != reference.shape or np.any(scales <= 0) or not np.isfinite(scales).all():
        raise ValueError("Need matching full fields and positive frozen field scales")
    errors = [np.max(np.abs((current - np.roll(reference, shift, axis=-1)) / scales))
              for shift in range(current.shape[-1])]
    shift = int(np.argmin(errors))
    signed = shift if shift <= current.shape[-1] // 2 else shift - current.shape[-1]
    return dict(best_integer_shift_cells=signed, unshifted_scaled_max_residual=float(errors[0]),
                aligned_scaled_max_residual=float(errors[shift]),
                interpretation="Track alignment shift over time; aligned agreement with changing shift indicates motion of an existing pattern, not formation of a second site.")


def spatial_readout(grid, primary, primary_basal):
    """Descriptive fields; never label a moving single site as a new front.

    Frozen basal half-contrast is a visualization/segmentation reference, not
    a universal front-establishment or survival threshold. Keep raw profiles.
    """
    primary = np.asarray(primary)
    level, original = basal_core(primary_basal)
    components = circular_components(primary >= level)
    peaks, properties = find_peaks(np.tile(primary, 3), prominence=0)
    select = (peaks >= grid.cells) & (peaks < 2 * grid.cells)
    peaks = peaks[select] - grid.cells
    prominence = properties["prominences"][select]
    initial_peak = int(np.argmax(primary_basal))
    return dict(frozen_basal_half_contrast=level, connected_regions=len(components),
        component_cells=[component.tolist() for component in components],
        component_overlaps_original_site=[bool(original[component].any()) for component in components],
        peak_cells=peaks.tolist(), peak_values=primary[peaks].tolist(), prominence=prominence.tolist(),
        largest_peak_displacement=float(grid.arc(grid.x[int(np.argmax(primary))], grid.x[initial_peak])),
        interpretation="One moving region permits drift/broadening; distinct simultaneous regions support candidate coexistence. Raw peak/valley histories and later cue-free fate remain required.")


def operator_checks():
    grid = RingGrid()
    rng = np.random.default_rng(20261017)
    results = []
    for key in KEYS:
        model = native_model(key)
        rest = native_rest(key, grid)
        np.testing.assert_array_equal(periodic_laplacian(rest, grid.dx), np.zeros_like(rest))
        basal_signal = np.zeros(grid.cells)
        rest_rate = native_rhs(model, grid, rest, basal_signal)
        if key != "wave_pinning":
            assert np.max(np.abs(rest_rate)) < 1e-8
        state = rest * rng.uniform(0.5, 1.5, size=rest.shape)
        cue = ring_cue(grid, 80, 3.0, {"wave_pinning": 0.1, "goryachev": 6.4, "holmes_model3": 1.0}[key])
        rate = native_rhs(model, grid, state, cue)
        pool_rate = pool_amounts(model, grid, rate)
        scaled_mass_error = float(np.max(np.abs(pool_rate)) / (1 + np.sum(np.abs(rate)) * grid.dx))
        assert scaled_mass_error < 1e-13
        rotation_error = 0.0
        for shift in (1, 17, 200, 399):
            rotated_rate = native_rhs(model, grid, np.roll(state, shift, axis=1), np.roll(cue, shift))
            error = np.max(np.abs(rotated_rate - np.roll(rate, shift, axis=1)))
            rotation_error = max(rotation_error, float(error))
            assert error < 1e-11
        # Confirm that basal external stimulus is zero; the model's own basal
        # activation remains in the native reaction method unchanged.
        local0 = model.reaction(0, grid.x, state, basal_signal)
        locals = model.reaction(0, grid.x, state, cue)
        delta = locals - local0
        expected = cue * state[1] / (model.rac_total if key == "holmes_model3" else 1.0)
        np.testing.assert_allclose(delta[0], expected, rtol=2e-12, atol=2e-14)
        np.testing.assert_allclose(delta[1], -expected, rtol=2e-12, atol=2e-14)
        if key == "holmes_model3":
            np.testing.assert_array_equal(delta[2:], np.zeros_like(delta[2:]))
        # Finite-difference dependency check includes seam couplings and all
        # native local species couplings; no ODE trajectory is integrated.
        small = RingGrid(4, 40)
        y = state[:, :40].copy()
        signal = cue[:40]
        base = native_rhs(model, small, y, signal).ravel()
        pattern = periodic_sparsity(len(model.state_names), small.cells).toarray().astype(bool)
        outside = 0.0
        for column in range(y.size):
            perturbed = y.ravel().copy()
            perturbed[column] += 1e-6
            derivative = (native_rhs(model, small, perturbed.reshape(y.shape), signal).ravel() - base) / 1e-6
            outside = max(outside, float(np.max(np.abs(derivative[~pattern[:, column]]), initial=0)))
        assert outside == 0
        results.append(dict(model=key, species=list(model.state_names), diffusivities=list(model.diffusivities),
            conserved_groups=[list(group) for group in model.conserved_groups],
            native_homogeneous_pool_means=(pool_amounts(model, grid, rest) / grid.circumference).tolist(),
            native_basal_reaction_max_abs=float(np.max(np.abs(rest_rate))),
            scaled_pool_derivative_residual=scaled_mass_error, rotation_covariance_max_abs=rotation_error,
            cue_on_minus_basal_identity=True, finite_difference_dependencies_outside_periodic_sparsity=outside))
    # Independent Fourier eigenmode and seam-face checks of periodic diffusion.
    phase = 2 * np.pi * np.arange(grid.cells) / grid.cells
    wave = np.cos(7 * phase)
    eigenvalue = -4 * np.sin(7 * np.pi / grid.cells) ** 2 / grid.dx ** 2
    spectral_error = float(np.max(np.abs(periodic_laplacian(wave, grid.dx) - eigenvalue * wave)))
    assert spectral_error < 1e-10
    impulse = np.zeros(grid.cells)
    impulse[-1] = 1
    seam = periodic_laplacian(impulse, grid.dx)
    np.testing.assert_allclose(seam[0], 1 / grid.dx ** 2, rtol=0, atol=1e-12)
    np.testing.assert_allclose(seam.sum(), 0, rtol=0, atol=1e-12)
    for center in (0, 80, 120, 160, 200):
        np.testing.assert_array_equal(ring_cue(grid, center, 3, 1), np.roll(ring_cue(grid, 0, 3, 1), center))
    # Synthetic spatial frames test the drift/coexistence diagnostic without
    # integrating any dynamics or claiming biological front establishment.
    basal = 0.1 + np.exp(-(grid.distance(grid.x, grid.x[0]) / 0.7) ** 2)
    drift = np.roll(basal, 80)
    coexist = basal + np.exp(-(grid.distance(grid.x, grid.x[120]) / 0.7) ** 2)
    drift_readout = spatial_readout(grid, drift, basal)
    coexist_readout = spatial_readout(grid, coexist, basal)
    assert drift_readout["connected_regions"] == 1 and len(drift_readout["peak_cells"]) == 1
    assert coexist_readout["connected_regions"] == 2 and len(coexist_readout["peak_cells"]) == 2
    alignment = translation_compare(np.stack([drift, 2 - drift]), np.stack([basal, 2 - basal]), [1, 1])
    assert alignment["best_integer_shift_cells"] == 80
    assert alignment["aligned_scaled_max_residual"] == 0 and alignment["unshifted_scaled_max_residual"] > 0.9
    locations = valid_challenger_locations(grid, basal, [0, 80, 120, 160, 200], 3, 1)
    assert not locations[0]["eligible"] and all(row["eligible"] for row in locations[1:])
    return dict(passed=True, integrated_trajectories=0, slurm_jobs_submitted=0, models=results,
        periodic_fourier_mode_max_abs_error=spectral_error, seam_flux_conservation=True,
        exact_cue_rotation_and_equal_discrete_integral=True,
        synthetic_one_region_drift_vs_two_region_coexistence=True,
        full_field_translation_diagnostic=alignment,
        native_start_note="Wave-pinning's retained u0=0.2683,v0=2 is a rounded homogeneous start, not an exactly solved equilibrium; its nonzero reaction residual is recorded, not removed.",
        actual_periodic_polarized_basal_preparations="Not run; existence, convergence and single-site eligibility unestablished.")


def pilot_plan():
    return dict(status="Proposed for root review; no integration authorized by this prototype step",
        scope=list(KEYS), boundary_condition="Periodic cell-centred ring; unchanged native local kinetics and pool definitions",
        geometry=dict(circumference=40, cells=400, dx=0.1, units="Each model's native length units; no biological cross-model calibration"),
        preparation=dict(initial="Fresh native homogeneous low/rest state, never an endpoint checkpoint",
            cue="Symmetric raised cosine, total support 3, centered on cell 0",
            suggested_peak={"wave_pinning":0.1,"goryachev":0.1,"holmes_model3":1.0},
            tentative_duration_grid=[10,100,1000],
            note="Bounded feasibility choices only; prepare one model at a time, stop at a qualifying state, do not tune by challenger distance.",
            basal="Zero external input, retaining native basal reaction terms",
            settling="Two 1000-unit windows plus an additional 10000-unit cue-free holdout; frozen per-field scales and derivative/range tolerances must be documented before launch.",
            acceptance="Positive finite fields, conserved pools, settled single connected active site, resolved native fields; report every periodic local peak/prominence and valley as well as connected regions. Homogeneous/multisite/drifting nonsettled states remain explicit statuses.",
            translation_diagnostic="Record fixed-frame full-field temporal residual, best integer cyclically aligned residual and its shift over time. Translation is not a polarity reversal; alignment does not by itself qualify settling or establish a neutral eigenvalue.",
            checkpoint="A new periodic-preparation checkpoint with full fields, times, grid, exact hash and source hashes; reused byte-for-byte at every challenger location"),
        challenge=dict(total_support_width=3, candidate_offset_cells=[80,120,160,200],
            candidate_shortest_arc=[8,12,16,20], distance_domain="0 ≤ d ≤ C/2; clockwise/counterclockwise copies are numerical symmetry controls, not independent samples",
            center_rule="Add integer offset to the actual prepared incumbent peak index, modulo N; report actual shortest-arc distance. Every trial uses the same unrotated saved basal array.",
            nonoverlap="Reject cue supports intersecting the measured fixed original incumbent half-contrast core; d=0 is a cue/operator control only.",
            peak={"wave_pinning":0.1,"goryachev":6.4,"holmes_model3":1.0},
            recommendation="Begin with a common maintained-cue window within each model, not a duration-threshold search.",
            decisions_needed=["Maintained-cue diagnostic horizon per model versus finite challenge plus cue-free survival readout.",
                              "Preparation cue duration and settled-field tolerances; the tentative grid is not calibrated on rings.",
                              "Whether C=40 and dx=0.1 resolve each newly prepared site; mesh refinement precedes interpretation."]),
        readout=dict(primary="Full time-resolved primary and inactive field profiles around the ring",
            diagnostics="Track all circular peak positions, prominence/valleys, connected active regions, fixed-frame/aligned full-field changes and incumbent displacement; keep original-site and cue-site regional activity separately. Do not reduce periodic branch structure to its first harmonic.",
            threshold="The frozen basal half-contrast cut is descriptive segmentation only, not a universal front-establishment gate.",
            distinction="One moving site may be incumbent drift; candidate competitor formation requires a distinct simultaneously resolved activity region/peak while the incumbent remains. Overlap/merger and ambiguous broadening remain explicit.",
            outcome="Any later takeover/survival definition requires cue-free holdout and circular alignment so rotation does not masquerade as disappearance or replacement."),
        execution=dict(environment="/home/4473331/.conda/envs/scaling_sandbox/bin/python",
            launch="None. Root reviews this minimal plan before any integration or SLURM submission.",
            proposed_slurm="QoS small; one CPU/model worker, 3 GB, 1 hour initially; no array concurrency throttle; bound wall time and preserve incomplete statuses.",
            output="Fresh data/runs/ring_distance_pilot_<date>/ with frozen sources, separate preparation and challenge archives; no canonical model or retained-run writes."),
        interpretation="A new geometry-control experiment. Periodic reactions/diffusion are feasible algebraically, but polarized-branch existence and distance effects have not been demonstrated.")


def main():
    checks = operator_checks()
    plan = pilot_plan()
    paths = []
    for name, value in (("operator_checks.json", checks), ("minimal_pilot_plan.json", plan)):
        path = dev_output_path(OUT, name)
        path.write_text(json.dumps(value, indent=2) + "\n")
        paths.append(path)
    source_paths = [Path(__file__), Path(__file__).with_name("length_dev_output_layout.py"),
        ROOT / "code/core/models/goryachev.py", ROOT / "code/core/models/wave_pinning.py",
        ROOT / "code/core/models/holmes_model3.py", ROOT / "code/core/src/polarity1d/holmes.py",
        ROOT / "code/core/src/polarity1d/contract.py", ROOT / "code/core/src/polarity1d/solvers.py",
        ROOT / "code/core/native.py", ROOT / "code/core/polarity_establishment.py"]
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer="finished",
        artifact_root="../..", status="Operator-checked development prototype; minimal pilot awaits review",
        integrated_trajectories=0, slurm_jobs_submitted=0,
        source_sha256={str(path.relative_to(ROOT)): digest(path) for path in source_paths},
        sha256={dev_artifact_key(OUT, path): digest(path) for path in paths},
        reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__)))
    dev_output_path(OUT, "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(output=str(OUT), passed=checks["passed"], integrated_trajectories=0,
                          models=checks["models"], periodic_fourier_mode_max_abs_error=checks["periodic_fourier_mode_max_abs_error"]), indent=2))


if __name__ == "__main__":
    main()
