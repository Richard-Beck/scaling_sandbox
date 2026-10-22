"""Immutable matched-initial-physical-arc arm of the stationary-cue experiment.

Freeze and scientific preflight precede a separately authorized submission.
The native worker is exactly the accepted three-size angle-arm worker.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import revision_standard_angles_sizes_20261019 as prior

ROOT = prior.ROOT
DEFAULT = ROOT / "data/runs/revision_initial_arc_sizes_20261019"
OUT = prior.OUT.parent / "standard_arclength_sizes"
ANGLE_RUN = prior.DEFAULT
PYTHON = prior.PYTHON
sha, read, write, now = prior.sha, prior.read, prior.write, prior.now
REFERENCE_FRACTIONS = (1 / 16, 1 / 8, 1 / 4, 1 / 2)


def freeze(run):
    if (run / "manifest.json").exists():
        raise FileExistsError("A completed immutable freeze already exists")
    for folder in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / folder).mkdir(parents=True, exist_ok=True)
    angle_manifest = read(ANGLE_RUN / "manifest.json")
    assert sha(ANGLE_RUN / "manifest.json") == read(ANGLE_RUN / "freeze.json")["manifest_sha256"]
    sources = {}
    for rel, digest in angle_manifest["source_sha256"].items():
        original = ANGLE_RUN / "frozen_sources" / rel
        assert sha(original) == digest
        target = run / "frozen_sources" / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
        sources[rel] = digest
    # The helpers used during freeze must match the immutable implementation.
    for module in (prior, prior.base):
        rel = str(Path(module.__file__).relative_to(ROOT))
        assert sha(module.__file__) == sources[rel]
    rel = Path(__file__).relative_to(ROOT)
    target = run / "frozen_sources" / rel
    shutil.copy2(__file__, target)
    sources[str(rel)] = sha(target)

    summary_path = prior.OUT / "supporting/data/stationary_sizes_summary.csv"
    inventory_path = ANGLE_RUN / "results/condition_inventory.csv"
    summary = pd.read_csv(summary_path)
    inventory = pd.read_csv(inventory_path)
    assert len(summary) == len(inventory) == 240
    assert (inventory.status == "complete").all()
    spring_reference = next(b for b in angle_manifest["basals"] if b["model"] == "spring" and b["size_parameter"] == 1.)
    spring_reference_C = spring_reference["basal_circumference"]
    physical_grids = dict(chemical=[20. * f for f in REFERENCE_FRACTIONS],
                          spring=[spring_reference_C * f for f in REFERENCE_FRACTIONS])
    cases, jobs, contexts, basals, inputs, rejected = [], [], [], [], {}, []
    for basal in angle_manifest["basals"]:
        model, size = basal["model"], basal["size_parameter"]
        spring = model == "spring"
        original = ANGLE_RUN / basal["input"]
        assert sha(original) == basal["input_sha256"]
        target = run / basal["input"]
        shutil.copy2(original, target)
        inputs[basal["input"]] = sha(target)
        basals.append(dict(**basal, angle_arm_basal_path=str(original),
                           angle_arm_basal_sha256=sha(original)))
        with np.load(target) as z:
            state = z["state"].copy()
        Cb, n, q0 = basal["basal_circumference"], basal["cells"], basal["preparation_axis_fraction"]
        q = (np.arange(n) + .5) / n if spring else prior.base.ring.RingGrid(Cb, n).x / Cb
        u, v = (state[:n], state[n:-1]) if spring else (state[0], state[1])
        core = prior.base.core_info(q, u, q0)
        M1 = Cb * np.mean(u * np.exp(2j * np.pi * q))
        offset = basal["initial_M1_offset_from_preparation_deg"]
        template = next(c for c in angle_manifest["cases"] if c["model"] == model and c["size_parameter"] == size)
        width, peak = template["full_physical_support"], template["requested_peak"]
        if spring:
            relevant = ["code/mechanisms/spring_ring_model_20261019.py", "code/mechanisms/spring_ring_closure_20261017.py", "code/core/models/spring.py"]
        else:
            relevant = ["code/mechanisms/ring_distance_prototype_20261017.py", f"code/core/models/{model}.py"]
            relevant += [p for p in sources if p.startswith("code/core/src/polarity1d/")]
        dynamics = {rel: sources[rel] for rel in relevant}
        for index, d0 in enumerate(physical_grids["spring" if spring else "chemical"]):
            # The smallest arm is exactly the accepted angular arm, including floating-point centres.
            smallest = size == (1. if spring else 20.)
            fraction = REFERENCE_FRACTIONS[index] if smallest else d0 / Cb
            assert np.isclose(fraction * Cb, d0, rtol=0, atol=1e-15)
            assert 0 < fraction <= .5
            centre = (q0 + fraction) % 1
            label = f"d_{index + 1}"
            if spring:
                signal, dC = prior.physical_raised_cosine(prior.Closure(size, n), state, centre, peak, width)
                evaluation = "Continuous physical raised cosine at fixed requested material centre; no cell snap"
            else:
                grid = prior.base.ring.RingGrid(Cb, n)
                signal, evaluation = prior.base.chemical_cue(grid, fraction, peak, width)
                dC = np.zeros(n)
            cuefile = run / "inputs" / f"{model}_C{size:g}_{label}_cue.npz"
            np.savez_compressed(cuefile, material_positions=q, initial_cue=signal, initial_cue_dC=dC,
                                centre_fraction=centre, initial_state=state, prescribed_angle_radians=2 * np.pi * fraction,
                                prescribed_initial_arc=d0)
            inputs[str(cuefile.relative_to(run))] = sha(cuefile)
            left, right = core["left"] * Cb, core["right"] * Cb
            overlap = sum(max(0., min(right, d0 + width / 2 + k * Cb) - max(left, d0 - width / 2 + k * Cb)) for k in (-1, 0, 1))
            patch = prior.base.interval_weights(q, centre, width / (2 * Cb))
            total_excess = Cb * np.mean(u - core["background"])
            context = dict(model=model, size_parameter=size, size_label=template["size_label"], circumference=Cb,
                resting_circumference=size if spring else None, cells=n, arm="common_initial_physical_arc",
                placement_label=label, angle_label=label, initial_arc_grid_index=index,
                common_initial_arc=d0, physical_distance_at_onset=d0, fraction=fraction,
                requested_angle_deg=360 * fraction, requested_angle_radians=2 * np.pi * fraction,
                preparation_axis_fraction=q0, preparation_axis_absolute_deg=360 * q0,
                initial_M1_absolute_phase_deg=float(np.rad2deg(np.angle(M1))), initial_M1_offset_from_preparation_deg=offset,
                target_angle_in_P_coordinates_deg=360 * fraction - offset,
                cue_centre_fraction=centre, cue_centre_physical_at_onset=centre * Cb,
                full_physical_support=width, initial_angular_support_deg=360 * width / Cb, requested_peak=peak,
                initial_sampled_cue_peak=float(signal.max()), initial_cue_integral=float(Cb * signal.mean()),
                initial_cue_shape_sha256=hashlib.sha256(np.asarray(signal, dtype="<f8").tobytes()).hexdigest(),
                cue_evaluation=evaluation, original_core_left_physical=left, original_core_right_physical=right,
                original_core_physical_width=right - left, cue_support_core_intersection=overlap,
                overlap_fraction_cue_support=overlap / width, overlap_fraction_original_core=overlap / (right - left),
                initial_overlap=bool(overlap > 1e-12), separate_challenger_interpretation_eligible=bool(overlap <= 1e-12),
                cue_support_contains_preparation_axis=bool(d0 <= width / 2),
                initial_unit_peak_cue_excess_fraction_total=float(Cb * np.mean((signal / peak) * (u - core["background"])) / total_excess),
                initial_geometric_support_excess_fraction_total=float(Cb * ((u - core["background"]) @ patch) / total_excess),
                initial_cue_weighted_active=float(signal @ u / signal.sum()), initial_cue_weighted_inactive=float(signal @ v / signal.sum()),
                input=basal["input"], input_sha256=sha(target), cue_input=str(cuefile.relative_to(run)),
                cue_input_sha256=sha(cuefile), fixed_scales=template["fixed_scales"])
            contexts.append(context)
            for duration in ((.1, 1., 10., 100., 1000.) if spring else (1., 10., 100., 300., 1000.)):
                case = dict(**context, condition_id=len(cases), cue_duration=duration, recovery=13000., native_units_uncalibrated=True)
                candidate = summary[(summary.model == model) & (summary.size_parameter == size) &
                                    np.isclose(summary.fraction, fraction, rtol=0, atol=1e-14) & (summary.cue_duration == duration)]
                proof = None
                if len(candidate):
                    assert len(candidate) == 1
                    selected = candidate.iloc[0].to_dict()
                    try:
                        proof = reuse_proof(case, selected, run, dynamics)
                        origin_manifest = read(Path(proof["source_run"]) / "manifest.json")
                        assert origin_manifest["solver"] == angle_manifest["solver"]
                        proof.update(identical_solver_configuration=True, solver=angle_manifest["solver"],
                                     source_angle_arm_condition_id=int(selected["condition_id"]))
                    except Exception as error:
                        rejected.append(dict(condition_id=case["condition_id"], reason=f"{type(error).__name__}: {error}"))
                if proof is not None:
                    case.update(execution="exact_archive_reuse", reuse=proof,
                                source_angle_arm_condition_id=proof["source_angle_arm_condition_id"])
                else:
                    case.update(execution="new_native_branch", worker_slot=len(jobs))
                    jobs.append(case)
                cases.append(case)
    assert len(cases) == 240 and len(contexts) == 48 and len(basals) == 12
    assert len(jobs) == 115 and not rejected, rejected
    manifest = dict(generated_at=now(), experiment="Matched initial physical arc arm of the stationary-cue experiment",
        stage="frozen_before_independent_preflight", requested_conditions=240, cases=cases, jobs=jobs, basals=basals,
        source_sha256=sources, input_sha256=inputs, rejected_reuse_candidates=rejected,
        reference_angle_arm_manifest=str(ANGLE_RUN / "manifest.json"), reference_angle_arm_manifest_sha256=sha(ANGLE_RUN / "manifest.json"),
        reference_angle_arm_summary=str(summary_path), reference_angle_arm_summary_sha256=sha(summary_path),
        reference_angle_arm_inventory=str(inventory_path), reference_angle_arm_inventory_sha256=sha(inventory_path),
        fixed_chemical_sizes=[20, 30, 40], resting_Spring_sizes=[1, 2, 3], initial_arc_grids=physical_grids,
        spring_reference_basal_circumference=spring_reference_C, reference_fractions=list(REFERENCE_FRACTIONS),
        overlap_cases_retained=True, excluded_for_overlap=0, solver=angle_manifest["solver"],
        cue_full_physical_support=dict(chemical=3., spring=.2), recovery=13000.,
        native_worker_implementation="code/mechanisms/revision_standard_angles_sizes_20261019.py:worker (byte-identical accepted frozen source)",
        angular_reference=angle_manifest["angular_reference"],
        Spring_coordinate_rule="Fixed material centre q_preparation+d0/C_basal. Prescribed arc subsequently equals d0*C(t)/C_basal; no moving cue.",
        interpretation="Arms match initial prescribed cue-centre arc, not instantaneous Spring arc or evolved pole separation. Same native basals, kinetics, channel, physical support, peak, durations and release horizon as the angle arm. All initial overlaps retained.")
    write(run / "manifest.json", manifest)
    pd.DataFrame(contexts).to_csv(run / "results/initial_context.csv", index=False)
    write(run / "results/basal_inventory.json", dict(generated_at=now(), basals=basals))
    counts = dict(generated_at=now(), requested=240, exact_reuse=125, new_native=115,
                  initial_overlap_conditions=sum(c["initial_overlap"] for c in cases), excluded_for_overlap=0,
                  qualified_basals=12, sources=len(sources), inputs=len(inputs),
                  manifest_sha256=sha(run / "manifest.json"), native_jobs_submitted=False)
    write(run / "freeze.json", counts)
    collect(run)
    print(json.dumps(counts, indent=2))


def reuse_proof(case, candidate, run, dynamics):
    """Exact native reuse, including both historical support-metadata schemas."""
    cue, release = Path(candidate["actual_cue_path"]), Path(candidate["actual_recovery_path"])
    assert sha(cue) == candidate["actual_cue_sha256"]
    assert sha(release) == candidate["actual_recovery_sha256"]
    with np.load(run / case["input"]) as z:
        initial = z["state"].copy()
    np.testing.assert_array_equal(prior.base.full_state(cue, 0.), initial)
    np.testing.assert_array_equal(prior.base.full_state(cue, case["cue_duration"]), prior.base.full_state(release, 0.))
    spring = case["model"] == "spring"
    with np.load(cue) as z:
        assert z["time"][-1] == case["cue_duration"]
        if spring:
            key = "cue_centre_fraction" if "cue_centre_fraction" in z else "cue_center_fraction"
            assert float(z[key]) == case["cue_centre_fraction"]
            assert float(z["cue_peak"]) == case["requested_peak"]
            if "full_physical_support" in z:
                assert float(z["full_physical_support"]) == case["full_physical_support"]
            if "resting_circumference" in z:
                assert float(z["resting_circumference"]) == case["resting_circumference"]
        else:
            with np.load(run / case["cue_input"]) as expected:
                np.testing.assert_array_equal(z["stimulus"], expected["initial_cue"])
    with np.load(release) as z:
        assert z["time"][-1] == case["recovery"]
        assert (float(z["cue_peak"]) == 0.) if spring else not np.count_nonzero(z["stimulus"])
    origin = prior.base.source_run(cue)
    for rel, digest in dynamics.items():
        assert sha(origin / "frozen_sources" / rel) == digest, (origin, rel)
    if spring:
        old = read(origin / "manifest.json")
        width = old.get("cue_full_physical_support")
        if isinstance(width, dict):
            width = width["spring"]
        if width is not None:
            assert width == case["full_physical_support"]
        assert any(c.get("model", "spring") == "spring" and
                   c.get("cue_centre_fraction", c.get("cue_center_fraction")) == case["cue_centre_fraction"] and
                   c.get("resting_circumference", c.get("resting_C", c.get("size_parameter"))) == case["resting_circumference"]
                   for c in old.get("cases", old.get("jobs", [])))
    return dict(status="exact_archive_reuse", actual_cue_path=str(cue), actual_cue_sha256=sha(cue),
                actual_recovery_path=str(release), actual_recovery_sha256=sha(release), source_run=str(origin),
                exact_initial_full_state=True, exact_phase_restart=True,
                exact_cue_profile_or_continuous_parameters=True, exact_duration_and_recovery_horizon=True,
                unchanged_native_dynamics_verified=True)


def submit(run):
    assert not (run / "launch.json").exists()
    manifest = read(run / "manifest.json")
    assert sha(run / "manifest.json") == read(run / "freeze.json")["manifest_sha256"]
    for rel, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / rel) == digest
    for rel, digest in manifest["input_sha256"].items():
        assert sha(run / rel) == digest
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = shlex.join([PYTHON, str(runner), "worker", "--run", str(run)]) + ' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=00:30:00",
                 f'--array=0-{len(manifest["jobs"]) - 1}', "--job-name=stationary-initial-arc", "--output=" + str(run / "logs/%A_%a.out"),
                 "--error=" + str(run / "logs/%A_%a.err"), "--wrap", "exec " + command]
    job = subprocess.check_output(arguments, text=True).strip()
    launch = read(run / "freeze.json")
    launch.update(native_jobs_submitted=True, array_job_id=job, submitted_at=now(), sbatch_arguments=arguments)
    write(run / "launch.json", launch)
    print(json.dumps(read(run / "launch.json"), indent=2))


def collect(run):
    # Same inventory schema and status collection as the accepted angle arm.
    prior.collect(run)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["freeze", "submit", "worker", "collect"])
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--slot", type=int)
    arguments = parser.parse_args()
    run = arguments.run.resolve()
    if arguments.action == "freeze":
        freeze(run)
    elif arguments.action == "submit":
        submit(run)
    elif arguments.action == "collect":
        collect(run)
    else:
        prior.worker(run, arguments.slot)
