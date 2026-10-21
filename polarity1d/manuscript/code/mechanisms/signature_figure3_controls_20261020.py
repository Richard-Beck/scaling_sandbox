"""Independent fresh fine-grid basal and mesh/time-step control payloads.

No scheduler submissions occur here. The reviewed main campaign's frozen
native solver and worker are inherited exactly; refined basals are freshly
prepared rather than interpolated from a coarse polarized checkpoint.
"""
from pathlib import Path
import argparse
import copy
import json
import shutil

import numpy as np
import signature_figure3_native_20261020 as native
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_basals_20261020 as basals
import signature_figure3_campaign_20261020 as campaign

ROOT = native.ROOT
DEFAULT = ROOT/"data/runs/signature_figure3_mesh_basals_20261020"
CONTROL_DEFAULT = ROOT/"data/runs/signature_figure3_controls_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def copy_sources(origin, target, hashes):
    for relative, digest in hashes.items():
        source, destination = origin/"frozen_sources"/relative, target/"frozen_sources"/relative
        assert sha(source) == digest
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def freeze_basals(run, primary_run):
    if run.exists():
        raise FileExistsError("Fresh control basal directory required")
    primary = campaign.verify(primary_run)
    original = basals.verify(Path(primary["fresh_basal_manifest"]).parent)
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/name).mkdir(parents=True)
    copy_sources(primary_run, run, primary["source_sha256"])
    jobs, inputs = [], {}
    for key in ("wave_pinning", "goryachev"):
        job = copy.deepcopy(next(j for j in original["jobs"] if j["model"] == key and j["size_parameter"] == 37.8))
        adapter = native.Adapter(key, 37.8, job["cells"]*2)
        state = adapter.initial()
        path = run/"inputs"/f"{key}_C37.8_N{adapter.cells}_fresh_homogeneous.npz"
        np.savez_compressed(path, state=state, material_positions=adapter.q,
                            model=key, size_parameter=37.8, cells=adapter.cells)
        inputs[str(path.relative_to(run))] = sha(path)
        job.update(slot=len(jobs), cells=adapter.cells, reference_density=40, nominal_spacing=.025,
            input=str(path.relative_to(run)), input_sha256=sha(path), fixed_scales=adapter.scales.tolist(),
            native_initial_amounts=adapter.amounts(state).tolist(), role="independent_fresh_fine_mesh_basal")
        jobs.append(job)
    plan = copy.deepcopy(original)
    plan.update(generated_at=now(), experiment="Fresh Figure 3 two-model mesh-control basals",
        jobs=jobs, requested_core_basals=2, source_sha256=primary["source_sha256"], input_sha256=inputs,
        main_campaign=str(primary_run), main_campaign_manifest_sha256=sha(primary_run/"manifest.json"),
        preparation_duration_search=False, coarse_state_interpolation=False)
    write(run/"manifest.json", plan)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    payload = copy.deepcopy(read(Path(primary["fresh_basal_manifest"]).parent/"slurm_payload.json"))
    args = payload["sbatch_arguments"]
    for i, arg in enumerate(args):
        if arg.startswith("--array="):
            args[i] = "--array=0-1"
        elif arg.startswith("--job-name="):
            args[i] = "--job-name=figure3-mesh-basals"
        elif arg.startswith("--output="):
            args[i] = "--output="+str(run/"logs/%A_%a.out")
        elif arg.startswith("--error="):
            args[i] = "--error="+str(run/"logs/%A_%a.err")
    args[-1] = args[-1].replace(str(Path(primary["fresh_basal_manifest"]).parent), str(run))
    payload.update(generated_at=now(), submitted=False, jobs=2,
        manifest_sha256=sha(run/"manifest.json"), sbatch_arguments=args)
    write(run/"slurm_payload.json", payload)
    return dict(generated_at=now(), prepared_fresh_fine_mesh_basals=2, cells_per_basal=1512,
                submitted=False, native_integrations=0, no_cue_horizon=53000.)


def audit_basals(run):
    result = basals.audit(run)
    result.update(requested=2, experiment="Two independently prepared fine-mesh control basals")
    write(run/"results/basal_audit.json", result)
    return result


def compile_controls(run, fine_basal_run, primary_run, worker_run, A1_only=False):
    if run.exists():
        raise FileExistsError("Fresh control campaign directory required")
    primary = campaign.verify(primary_run)
    worker_plan = campaign.verify(worker_run)
    for relative, digest in primary["source_sha256"].items():
        if relative != str(Path(campaign.__file__).relative_to(ROOT)):
            assert worker_plan["source_sha256"][relative] == digest, "Control native dynamics differ from primary producer"
    audit = audit_basals(fine_basal_run)
    if not audit["complete"] or audit["qualified"] != 2:
        raise RuntimeError("Both independent fine-mesh basals must pass qualification and long no-cue control")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results", "plans"):
        (run/name).mkdir(parents=True)
    copy_sources(worker_run, run, worker_plan["source_sha256"])
    sources, inputs, placements, jobs = worker_plan["source_sha256"], {}, [], []
    for variant in ("fine_mesh", "half_max_step_tighter_tolerances"):
        for model in ("wave_pinning", "goryachev"):
            if variant == "fine_mesh":
                record = next(r for r in audit["records"] if r["job"]["model"] == model)
                original = fine_basal_run/record["checkpoint"]
                assert sha(original) == record["checkpoint_sha256"]
                n = record["job"]["cells"]
            else:
                base = next(b for b in primary["basals"] if b["model"] == model and b["size_parameter"] == 37.8)
                original = primary_run/base["input"]
                n = base["cells"]
            target = run/"inputs"/f"{model}_{variant}_exact_basal.npz"
            shutil.copy2(original, target)
            inputs[str(target.relative_to(run))] = sha(target)
            with np.load(target) as z:
                state = z["state"].copy()
            adapter = native.Adapter(model, 37.8, n)
            candidates = geometry.compile_placements(model, adapter.q*37.8, 37.8, adapter.primary(state), adapter.support)
            for gap in (0., 4., 8.):
                placement = copy.deepcopy(next(p for p in candidates if p["placement_status"] == "feasible" and abs(p["edge_gap_10"]-gap) <= 1e-10))
                placement.update(placement_id=len(placements), model=model, size_parameter=37.8, cells=n,
                    cue_peak=native.PEAKS[model], input=str(target.relative_to(run)), input_sha256=sha(target),
                    control_variant=variant, matched_gap=gap)
                placements.append(placement)
                common = dict(placement_id=placement["placement_id"], control=variant,
                    matched_reference_model=model, matched_reference_size=37.8, matched_reference_gap=gap)
                if variant == "half_max_step_tighter_tolerances":
                    common["solver_options"] = dict(max_step=5., rtol=1e-10, atol=1e-12)
                jobs.append(dict(**common, slot=len(jobs), arm="A1"))
                for ratio in (() if A1_only else (.5, 1., 2.)):
                    jobs.append(dict(**common, slot=len(jobs), arm="B", target_ratio=ratio))
    plan = copy.deepcopy(primary)
    plan.update(generated_at=now(), experiment="Figure 3 independent mesh and timestep controls",
        source_sha256=sources, input_sha256=inputs, placements=placements, stages=dict(controls=jobs),
        main_campaign=str(primary_run), main_campaign_manifest_sha256=sha(primary_run/"manifest.json"),
        reviewed_worker_campaign=str(worker_run), reviewed_worker_manifest_sha256=sha(worker_run/"manifest.json"),
        B_included=not A1_only,
        fine_basal_audit=str(fine_basal_run/"results/basal_audit.json"),
        fine_basal_audit_sha256=sha(fine_basal_run/"results/basal_audit.json"),
        control_interpretation="Spatial refinement compares independently qualified native basals at the same exact physical edge gap; timestep refinement restarts the exact coarse basal. Native kinetics and cue support/peak are unchanged.")
    write(run/"manifest.json", plan)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    write(run/"plans/controls.json", dict(generated_at=now(), jobs=jobs))
    write(run/"plans/controls.sha256.json", dict(sha256=sha(run/"plans/controls.json")))
    write(run/"plans/controls_slurm_payload.json", campaign.slurm_payload(run, "controls", jobs))
    return dict(jobs=len(jobs), A1=12, B=0 if A1_only else 36, placements=12, submitted=False,
                independent_fine_basals=2, exact_coarse_basal_timestep_restarts=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze-basals", "audit-basals", "compile"))
    parser.add_argument("--run", type=Path)
    parser.add_argument("--primary-run", type=Path, default=campaign.DEFAULT)
    parser.add_argument("--fine-basal-run", type=Path, default=DEFAULT)
    parser.add_argument("--worker-run", type=Path, default=ROOT/"data/runs/signature_figure3_withdrawals_20261020")
    parser.add_argument("--A1-only", action="store_true")
    args = parser.parse_args()
    if args.action == "freeze-basals":
        result = freeze_basals((args.run or DEFAULT).resolve(), args.primary_run.resolve())
    elif args.action == "audit-basals":
        result = audit_basals((args.run or DEFAULT).resolve())
    else:
        result = compile_controls((args.run or CONTROL_DEFAULT).resolve(), args.fine_basal_run.resolve(), args.primary_run.resolve(), args.worker_run.resolve(), args.A1_only)
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
