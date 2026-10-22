"""Full Figure 3 placement/dependency compiler and native A1/A2/B workers.

Preparation and source freezing are separate from actual scheduler submission.
The compiler requires independently audited fresh basals. A2 is compiled only
after its own A1 handover exists; missing work never becomes right censoring.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import copy
import json
import os
import shlex
import shutil
import time

import numpy as np
import signature_figure3_native_20261020 as native
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_basals_20261020 as basals

ROOT = native.ROOT
DEFAULT = ROOT/"data/runs/signature_figure3_campaign_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def slurm_payload(run, stage, jobs):
    if not jobs:
        return dict(submitted=False, jobs=0, sbatch_arguments=None)
    runner = run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command = "exec "+shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", native.PYTHON, str(runner), "worker",
        "--run", str(run), "--stage", stage])+' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=4G",
        "--time=01:00:00", "--signal=B:USR1@60", f"--array=0-{len(jobs)-1}",
        "--job-name=figure3-"+stage, "--output="+str(run/f"logs/{stage}_%A_%a.out"),
        "--error="+str(run/f"logs/{stage}_%A_%a.err"), "--wrap", command]
    return dict(generated_at=now(), submitted=False, jobs=len(jobs), sbatch_arguments=arguments,
                no_concurrency_throttle=True, resource_reason="One serial native BDF integration per task, complete fields retained")


def freeze_campaign(run, basal_run):
    if run.exists():
        raise FileExistsError("Fresh immutable campaign directory required")
    basal_run = basal_run.resolve()
    basal_plan = basals.verify(basal_run)
    audit = basals.audit(basal_run)
    if not audit["complete"]:
        raise RuntimeError("Fresh basal qualification/control workers are not all terminal")
    if not audit["qualified"]:
        raise RuntimeError("No qualified fresh basal; do not substitute old endpoints")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results", "plans"):
        (run/name).mkdir(parents=True)
    hashes = {}
    # Exact native producer sources are inherited, rather than trusting current
    # native files to still match the basal preparation producer.
    for relative, digest in basal_plan["source_sha256"].items():
        original, target = basal_run/"frozen_sources"/relative, run/"frozen_sources"/relative
        assert sha(original) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
        hashes[relative] = digest
    for source in (Path(__file__), Path(geometry.__file__), Path(native.__file__), Path(basals.__file__)):
        relative = str(source.relative_to(ROOT))
        if relative in hashes and source.name != Path(geometry.__file__).name:
            assert sha(source) == hashes[relative], "Basal producer's native/geometry source has changed: "+relative
        target = run/"frozen_sources"/relative
        shutil.copy2(source, target)
        hashes[relative] = sha(target)
    all_placements, jobs, inputs, accepted, failures = [], [], {}, [], []
    for record in audit["records"]:
        job = record["job"]
        if record["status"] != "qualified":
            failures.append(dict(job=job, status=record["status"],
                source_status=str(basal_run/"workers"/f'{job["slot"]:03d}'/"status.json"),
                largest_size=job["size_parameter"] == native.SIZES[job["model"]][-1],
                retuned=False))
            continue
        original = basal_run/record["checkpoint"]
        assert sha(original) == record["checkpoint_sha256"]
        target = run/"inputs"/original.name.replace("qualified_pre_control_basal", f'{job["model"]}_C{job["size_parameter"]:g}_basal')
        shutil.copy2(original, target)
        inputs[str(target.relative_to(run))] = sha(target)
        with np.load(target) as z:
            state = z["state"].copy()
        adapter = native.Adapter(job["model"], job["size_parameter"], job["cells"])
        circumference = adapter.circumference(state)
        placements = geometry.compile_placements(job["model"], adapter.q*circumference,
            circumference, adapter.primary(state), adapter.support)
        accepted.append(dict(model=job["model"], size_parameter=job["size_parameter"], cells=job["cells"],
            input=str(target.relative_to(run)), input_sha256=sha(target),
            source_checkpoint=str(original), source_checkpoint_sha256=sha(original),
            exact_qualified_pre_control_checkpoint=True, actual_basal_circumference=circumference,
            fixed_scales=adapter.scales.tolist(),
            qualification_audit_sha256=sha(basal_run/"results/basal_audit.json")))
        for placement in placements:
            placement.update(placement_id=len(all_placements), model=job["model"],
                size_parameter=job["size_parameter"], cells=job["cells"],
                input=str(target.relative_to(run)), input_sha256=sha(target),
                cue_peak=native.PEAKS[job["model"]])
            if placement["placement_status"] != "feasible":
                all_placements.append(placement)
                continue
            signal, _ = adapter.cue(state, placement["cue_peak"], placement["cue_center_fraction"])
            placement.update(initial_sampled_cue_peak=float(signal.max()),
                initial_discrete_cue_integral=float(circumference*signal.mean()),
                cue_integral_error_vs_continuous=float(circumference*signal.mean()-placement["cue_peak"]*adapter.support/2))
            all_placements.append(placement)
            jobs.append(dict(slot=len(jobs), arm="A1", placement_id=placement["placement_id"]))
            if placement["classification_10"] == "external":
                for ratio in (.5, 1., 2.):
                    jobs.append(dict(slot=len(jobs), arm="B", placement_id=placement["placement_id"], target_ratio=ratio))
        # A mirrored non-antipodal external challenge tests the actual same
        # checkpoint, with its residual basal asymmetry kept measurable.
        external = [p for p in all_placements if p["model"] == job["model"] and p["size_parameter"] == job["size_parameter"] and p["classification_10"] == "external"]
        non_antipodal = [p for p in external if p["requested_role"] != "antipode"]
        for placement in ([non_antipodal[len(non_antipodal)//2]] if non_antipodal else []) + [external[-1]]:
            jobs.append(dict(slot=len(jobs), arm="A1", placement_id=placement["placement_id"],
                control="mirror_placement", mirrored=True,
                reference_placement_id=placement["placement_id"],
                antipode=placement["requested_role"] == "antipode"))
    B_jobs = [j for j in jobs if j["arm"] == "B"]
    jobs = [j for j in jobs if j["arm"] == "A1"]
    for stage_jobs in (jobs, B_jobs):
        for slot, job in enumerate(stage_jobs):
            job["slot"] = slot
    plan = dict(generated_at=now(), experiment="Figure 3 physical edge-gap competition and resolution",
        source_sha256=hashes, input_sha256=inputs, basals=accepted, failed_basals=failures,
        fresh_basal_manifest=str(basal_run/"manifest.json"), fresh_basal_manifest_sha256=sha(basal_run/"manifest.json"),
        fresh_basal_audit=str(basal_run/"results/basal_audit.json"), fresh_basal_audit_sha256=sha(basal_run/"results/basal_audit.json"),
        placements=all_placements, stages=dict(primary=jobs, B=B_jobs),
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.),
        worker_soft_seconds=3300,
        horizons=dict(chemical_A1_initial=5000., chemical_A1_cap=20000., spring_A1_cap=2000.,
            chemical_A2_tail=13000., spring_A2_tail=1300., chemical_B_tail=20000., spring_B_tail=2000.),
        A2=dict(multipliers=[.5, .75, 1., 1.5, 2.],
            dependency="Each exact non-control A1 event-root handover; unavailable if no observed handover",
            compile_command="compile-a2", no_cap_substitution=True),
        B=dict(target_ratios=[.5, 1., 2.], collapse_fraction=.1,
            withdrawal_event="first of ratio reached or incumbent 90% excess loss; exact BDF dense root",
            floor_and_ROI="Initial basal floor/ROIs for ratio; frozen actual withdrawal floor/regions for instantaneous RHS budgets",
            resolution="One connected region at frozen withdrawal 50%-contrast level, with threshold25/75 sensitivities; cue-associated winner separately measured"),
        controls=dict(mirror=True, antipode=True, accepted_step_mass_positivity=True,
            threshold_sensitivity=[.25, .5, .75],
            mesh_and_timestep=dict(models=["wave_pinning", "goryachev"], size=37.8, gaps=[0., 4., 8.],
                mesh_factor=2, max_step=5., rtol=1e-10, atol=1e-12,
                fresh_fine_mesh_preparations_required=True, compile_command="compile-controls",
                status="pending independent fresh refined-grid basal preparation")),
        main_panels_and_supplements="docs/tasks/briefs/20261020_figure3_requirements.json",
        interpolation="No extrapolation of geometry-only remaps or matched-gap size effects; censored data remain lower bounds",
        global_vector="Native first circular moment with magnitude; heading unavailable at magnitude<=1e-10; never unwrap through cancellation",
        scientific_claims="Prospective; native physical results must settle model-specific distance effects")
    write(run/"manifest.json", plan)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    write(run/"plans/primary.json", dict(generated_at=now(), jobs=jobs))
    write(run/"plans/primary.sha256.json", dict(sha256=sha(run/"plans/primary.json")))
    write(run/"plans/primary_slurm_payload.json", slurm_payload(run, "primary", jobs))
    write(run/"plans/B.json", dict(generated_at=now(), jobs=B_jobs))
    write(run/"plans/B.sha256.json", dict(sha256=sha(run/"plans/B.json")))
    write(run/"plans/B_slurm_payload.json", slurm_payload(run, "B", B_jobs))
    return preflight(run)


def verify(run):
    plan = read(run/"manifest.json")
    assert sha(run/"manifest.json") == read(run/"freeze.json")["manifest_sha256"]
    for relative, digest in plan["source_sha256"].items():
        assert sha(run/"frozen_sources"/relative) == digest
    for relative, digest in plan["input_sha256"].items():
        assert sha(run/relative) == digest
    return plan


def preflight(run):
    plan = verify(run)
    seen, census = set(), {}
    for placement in plan["placements"]:
        if placement["placement_status"] != "feasible":
            assert placement["classification_10"] == "unavailable" and not placement["claim_eligible"]
            continue
        key = placement["model"], placement["size_parameter"]
        adapter = native.Adapter(*key, placement["cells"])
        with np.load(run/placement["input"]) as z:
            state = z["state"].copy()
        C = adapter.circumference(state)
        edge = geometry.contour(adapter.q*C, C, adapter.primary(state), .1)
        measured_gap = placement["center_distance"]-adapter.support/2-edge["right"]
        np.testing.assert_allclose(measured_gap, placement["edge_gap_10"], atol=1e-10)
        assert placement["short_route"] <= C/2+1e-10
        np.testing.assert_allclose(placement["short_route"]+placement["long_route"], C, atol=1e-10)
        if placement["classification_10"] == "external":
            assert placement["edge_gap_10"] >= -1e-10
        identity = (*key, placement["placement_index"])
        assert identity not in seen
        seen.add(identity)
        label = f"{key[0]}_C{key[1]:g}"
        if label not in census:
            census[label] = dict(total=0, external=0, borderline=0, internal=0, antipodes=0)
        census[label]["total"] += 1
        census[label][placement["classification_10"]] += 1
        census[label]["antipodes"] += int(placement["requested_role"] == "antipode")
    assert all(row["antipodes"] == 1 for row in census.values())
    jobs = plan["stages"]["primary"]
    result = dict(generated_at=now(), passed=True, placements=len(plan["placements"]),
        primary_jobs=len(jobs), A1=sum(j["arm"] == "A1" and not j.get("control") for j in jobs),
        B=len(plan["stages"]["B"]), mirror_controls=sum(bool(j.get("control")) for j in jobs),
        deferred_A2_maximum=sum(p["classification_10"] in ("external", "borderline") for p in plan["placements"])*5,
        census=census, failed_basals=plan["failed_basals"],
        infeasible_requested_placements=sum(p["placement_status"] != "feasible" for p in plan["placements"]),
        exact_continuous_geometry_and_antipodes_verified=True, native_integrations=0,
        expected_full_field_tail=True, scientific_conclusions=False)
    write(run/"results/preflight.json", result)
    return result


def first_bracket(t, criterion):
    index = np.flatnonzero(criterion)
    if len(index) == 0:
        return dict(status="right_censored_at_completed_horizon", lower=float(t[-1]), upper=None)
    i = int(index[0])
    return dict(status="already_at_start" if i == 0 else "observed_saved_sample_crossing",
                lower=float(t[max(0, i-1)]), upper=float(t[i]))


def scalar(adapter, state, basal_primary, placement, fraction=.5):
    C = adapter.circumference(state)
    return geometry.region_observables(adapter.q*C, C, adapter.primary(state), basal_primary, placement, fraction)


def current_regions(adapter, state, basal_primary, placement, fraction=.5):
    C = adapter.circumference(state)
    old = placement["initial_10_percent_contour"]["origin"]*C/placement["basal_circumference"]
    cue = placement["cue_center_fraction"]*C
    return geometry.measured_regions(adapter.q*C, C, adapter.primary(state), old, cue, fraction)


def exchange_budget(adapter, state, placement):
    """Exact cue-free RHS in frozen withdrawal ROIs and a frozen floor.

    Chemical weights and background are held fixed. Spring freezes MATERIAL
    ROI fractions and concentration floor; physical weights expand with C, so
    dE/dt contains sum(w*du/dt)+Cdot/C*E. No threshold-motion or min(u)-floor
    derivative enters this distinct instantaneous observable. Symmetric state
    perturbation along the exact RHS verifies that exact convention.
    """
    C = adapter.circumference(state)
    p = adapter.primary(state)
    shape = current_regions(adapter, state, p, placement)
    rhs = adapter.rhs(state)
    if not shape["available"]:
        return dict(available=False, reason=shape["reason"], exact_native_rhs_saved=True)
    background = float(p.min())
    roi_records = []
    for role, index in (("incumbent", shape["incumbent_index"]), ("challenger", shape["challenger_index"])):
        desc = shape["region_descriptors"][index]
        weights = geometry.interval_weights(adapter.q*C, C, desc["left"], 0., desc["width"])
        excess = float((p-background) @ weights)
        field_rate = float(rhs[:adapter.cells] @ weights)
        geometric_rate = float(rhs[-1]/C*excess) if adapter.spring else 0.
        exact = field_rate+geometric_rate
        epsilon = 1e-5/max(1., float(np.max(abs(rhs)/adapter.scales)))
        def value(y):
            scale = adapter.circumference(y)/C if adapter.spring else 1.
            return float((adapter.primary(y)-background) @ (weights*scale))
        difference = (value(state+epsilon*rhs)-value(state-epsilon*rhs))/(2*epsilon)
        error = abs(difference-exact)
        assert error <= 1e-8+1e-5*abs(exact)
        if adapter.spring:
            u, v, _ = adapter.model.fields(state)
            extension = C-adapter.size
            from spring_ring_closure_20261017 import native as spring_native
            reaction = spring_native.auth["f"](u, v, 4., 5, 6, 2, 3.+80.*extension)
            diffusion = .01*(np.roll(u, 1)-2*u+np.roll(u, -1))*(adapter.cells/C)**2
            dilution = -u*rhs[-1]/C
        else:
            fields = state.reshape(-1, adapter.cells)
            reaction = adapter.model.reaction(0., adapter.grid.x, fields, np.zeros(adapter.cells))[0]
            diffusion = adapter.model.diffusivities[0]*native.chemical.periodic_laplacian(fields[0], adapter.grid.dx)
            dilution = np.zeros(adapter.cells)
        chemistry_rate, diffusion_rate, dilution_rate = [float(term @ weights) for term in (reaction, diffusion, dilution)]
        assert abs(chemistry_rate+diffusion_rate+dilution_rate+geometric_rate-exact) <= 1e-9*(1+abs(exact))
        roi_records.append(dict(role=role, frozen_withdrawal_floor=background,
            frozen_material_roi_weights=(weights/C).tolist(), withdrawal_excess=excess,
            exact_rate=exact, chemical_rate=chemistry_rate, diffusion_rate=diffusion_rate,
            concentration_dilution_rate=dilution_rate, physical_weight_dilation_rate=geometric_rate,
            directional_finite_difference=difference, finite_difference_error=error,
            finite_difference_epsilon=epsilon))
    return dict(available=True, convention="Frozen withdrawal floor and material ROIs; physical weights follow Spring circumference; no moving-level derivative",
                exact_cue_free_rhs=True, measured_regions=shape, regions=roi_records)


def winner(adapter, state, placement, level):
    primary = adapter.primary(state)
    groups = geometry.components(primary >= level)
    vector = geometry.activity_vector(primary, adapter.q)
    if len(groups) != 1:
        return dict(winner="unresolved", regions=len(groups), vector=vector)
    cells = groups[0]
    z = np.sum(np.maximum(primary[cells]-level, 0.)*np.exp(2j*np.pi*adapter.q[cells]))
    if abs(z) <= 1e-14:
        return dict(winner="location_unavailable", regions=1, vector=vector)
    center = float((np.angle(z)/(2*np.pi)) % 1)
    old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
    cue = placement["cue_center_fraction"]
    old_distance = abs(float(geometry.arc(center, old, 1.)))
    cue_distance = abs(float(geometry.arc(center, cue, 1.)))
    role = "challenger" if cue_distance < old_distance else "incumbent"
    if abs(cue_distance-old_distance) < 1e-8:
        role = "location_ambiguous"
    return dict(winner=role, regions=1, centroid_fraction=center,
                distance_to_incumbent_fraction=old_distance, distance_to_cue_fraction=cue_distance, vector=vector)


def worker(run, stage, slot):
    native.require_slurm()
    plan = verify(run)
    stage_path = run/f"plans/{stage}.json"
    assert sha(stage_path) == read(run/f"plans/{stage}.sha256.json")["sha256"]
    job = read(stage_path)["jobs"][slot]
    placement = copy.deepcopy(plan["placements"][job["placement_id"]])
    adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
    with np.load(run/placement["input"]) as z:
        initial = z["state"].copy()
    basal_primary = adapter.primary(initial).copy()
    if job.get("mirrored"):
        old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
        placement["cue_center_fraction"] = float((2*old-placement["cue_center_fraction"]) % 1)
    directory = run/"workers"/stage/f"{slot:04d}"
    directory.mkdir(parents=True)
    record = dict(generated_at=now(), producer="running", complete=False, status="starting",
        job=job, placement=placement, slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[],
        phase_continuity="Every withdrawal starts from exact native previous endpoint; no clipping or remapping")
    statuspath = directory/"status.json"
    write(statuspath, record)
    deadline = time.monotonic()+plan["worker_soft_seconds"]
    elapsed = 0.
    state = initial.copy()
    history_t, history_y = [], []

    def phase(kind, duration, peak, event=None):
        nonlocal elapsed, state
        if duration <= 1e-8:
            return None
        record.update(generated_at=now(), status=kind)
        write(statuspath, record)
        result = native.integrate(adapter, state, duration, peak, placement["cue_center_fraction"],
            event=event, deadline=deadline, solver_options=job.get("solver_options"))
        archive = directory/f'{len(record["phases"]):03d}_{kind}.npz'
        metadata = native.save_phase(archive, result, adapter, peak, placement["cue_center_fraction"])
        metadata.update(kind=kind, absolute_start=elapsed, path=str(archive.relative_to(run)))
        record["phases"].append(metadata)
        selection = slice(None) if not history_t else slice(1, None)
        history_t.extend((elapsed+result["time"][selection]).tolist())
        history_y.extend(result["states"][selection])
        elapsed += float(result["time"][-1])
        state = result["states"][-1].copy()
        write(statuspath, record)
        if result["status"] not in ("complete", "event", "event_at_start"):
            raise RuntimeError("Incomplete native phase "+result["status"])
        return result

    try:
        cap = 2000. if adapter.spring else 20000.
        handover_time = None
        if job["arm"] == "A1":
            first_cap = cap if adapter.spring else 5000.
            def handover(y):
                value = scalar(adapter, y, basal_primary, placement)
                return value["challenger_excess"]-value["incumbent_excess"]
            for target in ([first_cap] if adapter.spring else [first_cap, cap]):
                if elapsed >= target-1e-8:
                    continue
                result = phase("continuous_cue", target-elapsed, placement["cue_peak"],
                               event=handover if handover_time is None else None)
                if result and result["event_time"] is not None:
                    handover_time = elapsed
                    phase("continuous_cue_after_handover", target-elapsed, placement["cue_peak"])
                value = scalar(adapter, state, basal_primary, placement)
                frozen_level = float(basal_primary.min()+.5*np.ptp(basal_primary))
                fate = winner(adapter, state, placement, frozen_level)
                if fate["winner"] == "challenger" and value["incumbent_fraction"] <= .1:
                    break
            record.update(handover_time=handover_time,
                handover_status="observed_dense_BDF_event_root" if handover_time is not None else "right_censored_at_completed_horizon",
                handover_lower_bound=elapsed if handover_time is None else None,
                completed_cue_horizon=elapsed,
                final_winner=winner(adapter, state, placement, float(basal_primary.min()+.5*np.ptp(basal_primary))))
        elif job["arm"] == "A2":
            assert sha(Path(job["source_A1_status"])) == job["source_A1_status_sha256"]
            duration = job["cue_duration"]
            phase("pulse", duration, placement["cue_peak"])
            withdrawal = state.copy()
            tail = phase("cue_free_tail", 1300. if adapter.spring else 13000., 0.)
            record.update(cue_duration=duration, withdrawal_time=duration,
                final_winner=winner(adapter, state, placement, float(basal_primary.min()+.5*np.ptp(basal_primary))),
                finite_tail=native.finite_tail(adapter, tail["time"], tail["states"]))
        elif job["arm"] == "B":
            ratio = job["target_ratio"]
            def seed_event(y):
                value = scalar(adapter, y, basal_primary, placement)
                ratio_condition = value["challenger_excess"]-ratio*value["incumbent_excess"]
                collapse_condition = .1*value["basal_incumbent_excess"]-value["incumbent_excess"]
                return max(ratio_condition, collapse_condition)
            result = phase("seed", cap, placement["cue_peak"], event=seed_event)
            value = scalar(adapter, state, basal_primary, placement)
            excess_floor = 1e-10*value["basal_incumbent_excess"]
            signed_excess_gate = value["incumbent_excess"] > excess_floor and value["challenger_excess"] >= -excess_floor
            ratio_reached = signed_excess_gate and value["challenger_excess"] >= ratio*value["incumbent_excess"]-1e-8*value["basal_incumbent_excess"]
            collapsed = value["incumbent_fraction"] <= .1+1e-8
            record["signed_excess_ratio_gate"] = dict(available=bool(signed_excess_gate),
                incumbent_excess=value["incumbent_excess"], challenger_excess=value["challenger_excess"],
                positive_incumbent_excess_floor=excess_floor,
                rule="Incumbent excess >1e-10 basal incumbent excess; challenger excess >=negative roundoff floor. Signed native readouts remain archived.")
            if not result or result["event_time"] is None:
                record.update(seed_status="target_ratio_unreached_at_completed_seed_cap", target_ratio=ratio,
                    seed_lower_bound=elapsed, measured_withdrawal_gap=dict(available=False, reason="cue_not_withdrawn"),
                    B_resolution=dict(status="unavailable_seed_ratio_unreached"))
            else:
                withdrawal = state.copy()
                withdrawal_time = elapsed
                measured_gap = current_regions(adapter, withdrawal, basal_primary, placement)
                eligible = bool(measured_gap["available"] and np.isfinite(measured_gap.get("measured_gap", np.nan))
                                and signed_excess_gate and ratio_reached and not collapsed)
                if collapsed and ratio_reached:
                    seed_status = "ratio_and_incumbent_collapse_coincident"
                elif collapsed:
                    seed_status = "takeover_during_cue"
                elif not signed_excess_gate:
                    seed_status = "signed_excess_ratio_unavailable"
                elif not measured_gap["available"]:
                    seed_status = "not_a_two_region_withdrawal"
                else:
                    seed_status = "ratio_reached_with_two_regions"
                budget = exchange_budget(adapter, withdrawal, placement)
                archive = directory/"withdrawal_exact_rhs.npz"
                np.savez_compressed(archive, state=withdrawal, cue_free_rhs=adapter.rhs(withdrawal),
                    cue_on_rhs=adapter.rhs(withdrawal, placement["cue_peak"], placement["cue_center_fraction"]),
                    material_positions=adapter.q)
                record.update(seed_status=seed_status, B_distance_claim_eligible=eligible,
                    target_ratio=ratio, actual_withdrawal_ratio=value["excess_ratio"],
                    withdrawal_time=withdrawal_time,
                    measured_withdrawal_gap=measured_gap,
                    exact_rhs_archive=str(archive.relative_to(run)), exact_rhs_sha256=sha(archive),
                    instantaneous_cue_free_exchange=budget)
                tail = phase("cue_free_resolution_tail", cap, 0.)
                threshold_results = {}
                for fraction in (.25, .5, .75):
                    level = float(adapter.primary(withdrawal).min()+fraction*np.ptp(adapter.primary(withdrawal)))
                    fates = [winner(adapter, y, placement, level) for y in tail["states"]]
                    resolved = np.asarray([f["regions"] == 1 for f in fates])
                    threshold_results[str(fraction)] = dict(resolution=first_bracket(tail["time"], resolved),
                        final_winner=fates[-1], frozen_withdrawal_level=level,
                        measured_withdrawal_gap=current_regions(adapter, withdrawal, basal_primary, placement, fraction))
                resolution_result = threshold_results["0.5"] if eligible else dict(
                    status="unavailable_for_two_region_distance_claim", reason=seed_status,
                    raw_threshold_50_result=threshold_results["0.5"],
                    note="A single-region withdrawal is not a zero cue-free two-region resolution time")
                record.update(B_resolution=resolution_result, threshold_sensitivity=threshold_results,
                    finite_tail=native.finite_tail(adapter, tail["time"], tail["states"]))
        else:
            raise ValueError("Unknown arm: "+job["arm"])
        t, y = np.asarray(history_t), np.asarray(history_y)
        sensitivity = {}
        arrays = dict(time=t, material_positions=adapter.q)
        for fraction in (.25, .5, .75):
            values = [scalar(adapter, state_y, basal_primary, placement, fraction) for state_y in y]
            old = np.asarray([v["incumbent_fraction"] for v in values])
            challenger = np.asarray([v["challenger_fraction"] for v in values])
            arrays[f"incumbent_fraction_{fraction}"] = old
            arrays[f"challenger_fraction_{fraction}"] = challenger
            sensitivity[str(fraction)] = dict(incumbent_50_percent_excess_loss=first_bracket(t, old <= .5),
                challenger_50_percent_basal_incumbent=first_bracket(t, challenger >= .5),
                handover=first_bracket(t, challenger >= old))
        vectors = [geometry.activity_vector(adapter.primary(state_y), adapter.q) for state_y in y]
        arrays.update(vector_real=np.asarray([v["real"] for v in vectors]),
            vector_imag=np.asarray([v["imag"] for v in vectors]),
            vector_magnitude=np.asarray([v["magnitude"] for v in vectors]),
            vector_heading=np.asarray([v["heading_radians"] if v["heading_radians"] is not None else np.nan for v in vectors]),
            vector_direction_eligible=np.asarray([v["direction_available"] for v in vectors]))
        for floor in (.05, .1, .2):
            arrays[f"vector_magnitude_eligible_{floor}"] = arrays["vector_magnitude"] >= floor
        archive = directory/"scalar_and_vector_readouts.npz"
        np.savez_compressed(archive, **arrays)
        record.update(status="complete", readouts=str(archive.relative_to(run)), readouts_sha256=sha(archive),
            scalar_threshold_sensitivity=sensitivity, completed_native_horizon=elapsed,
            final_native_fields_archive=record["phases"][-1]["path"],
            global_vector_note="No heading unwrap through magnitude<=1e-10; retain cancellation and all raw complex-vector components")
    except Exception as error:
        record.update(status="error", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        record.update(generated_at=now(), producer="finished", complete=True)
        write(statuspath, record)


def compile_a2(run):
    plan = verify(run)
    path = run/"plans/A2.json"
    if path.exists():
        raise FileExistsError("A2 plan is frozen; do not rewrite while producers may run")
    jobs, unavailable, waiting = [], [], []
    for job in plan["stages"]["primary"]:
        if job["arm"] != "A1" or job.get("control"):
            continue
        placement = plan["placements"][job["placement_id"]]
        if placement["classification_10"] not in ("external", "borderline"):
            continue
        source = run/"workers/primary"/f'{job["slot"]:04d}'/"status.json"
        if not source.exists() or not read(source).get("complete"):
            waiting.append(job["placement_id"]); continue
        record = read(source)
        if record["status"] != "complete":
            unavailable.append(dict(placement_id=job["placement_id"], status="A1_native_worker_error", source_status=str(source))); continue
        duration = record.get("handover_time")
        if duration is None or duration <= 0:
            unavailable.append(dict(placement_id=job["placement_id"], status="own_handover_unobserved_or_zero", source_status=str(source))); continue
        for multiplier in (.5, .75, 1., 1.5, 2.):
            jobs.append(dict(slot=len(jobs), arm="A2", placement_id=job["placement_id"],
                multiplier=multiplier, cue_duration=multiplier*duration, own_handover_time=duration,
                source_A1_status=str(source), source_A1_status_sha256=sha(source)))
    if waiting:
        raise RuntimeError(f"A2 dependencies still running/not started at {len(waiting)} placements")
    write(path, dict(generated_at=now(), jobs=jobs, unavailable=unavailable, required_multipliers=[.5, .75, 1., 1.5, 2.]))
    write(run/"plans/A2.sha256.json", dict(sha256=sha(path)))
    write(run/"plans/A2_slurm_payload.json", slurm_payload(run, "A2", jobs))
    return dict(jobs=len(jobs), unavailable_placements=unavailable, submitted=False)


def freeze_B(run, source_run):
    """A separate source version for B while approved A1 producers run."""
    if run.exists():
        raise FileExistsError("Fresh immutable B run required")
    original = verify(source_run)
    for folder in ("inputs", "frozen_sources", "workers", "logs", "results", "plans"):
        (run/folder).mkdir(parents=True)
    plan = copy.deepcopy(original)
    for relative, digest in original["source_sha256"].items():
        old, target = source_run/"frozen_sources"/relative, run/"frozen_sources"/relative
        assert sha(old) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(old, target)
    relative = str(Path(__file__).relative_to(ROOT))
    target = run/"frozen_sources"/relative
    shutil.copy2(__file__, target)
    plan["source_sha256"][relative] = sha(target)
    for relative, digest in original["input_sha256"].items():
        source, target = source_run/relative, run/relative
        assert sha(source) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    jobs = copy.deepcopy(original["stages"]["B"])
    plan.update(generated_at=now(), experiment="Figure 3 qualified two-region cue-free withdrawals",
        source_A1_campaign=str(source_run), source_A1_manifest_sha256=sha(source_run/"manifest.json"),
        stages=dict(B=jobs),
        B_eligibility="Exactly two geometrically associated withdrawal regions with finite measured gap and positive incumbent/signed challenger excess gate. One-region or ambiguous withdrawals are explicit alternate outcomes, not zero-time two-region resolutions.",
        B_coincident_events="Ratio and incumbent collapse coincidence remains explicit; raw signed excesses and exact native fields retained")
    write(run/"manifest.json", plan)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    write(run/"plans/B.json", dict(generated_at=now(), jobs=jobs))
    write(run/"plans/B.sha256.json", dict(sha256=sha(run/"plans/B.json")))
    write(run/"plans/B_slurm_payload.json", slurm_payload(run, "B", jobs))
    verify(run)
    result = dict(generated_at=now(), passed=True, source_A1_untouched=True, jobs=len(jobs),
        placements=sum(p["classification_10"] == "external" for p in plan["placements"]),
        ratios=[.5, 1., 2.], two_region_withdrawal_required=True, submitted=False,
        source_manifest_sha256=sha(source_run/"manifest.json"), manifest_sha256=sha(run/"manifest.json"))
    write(run/"results/preflight.json", result)
    return result


def collect(run):
    plan = verify(run)
    stages, all_records = {}, []
    for name in ("primary", "B", "A2", "controls"):
        path = run/f"plans/{name}.json"
        if not path.exists():
            stages[name] = dict(status="not_compiled", complete=False); continue
        jobs = read(path)["jobs"]
        records = []
        for job in jobs:
            statuspath = run/"workers"/name/f'{job["slot"]:04d}'/"status.json"
            record = read(statuspath) if statuspath.exists() else dict(job=job, complete=False, status="not_started")
            records.append(record)
        stages[name] = dict(jobs=len(jobs), complete=all(r.get("complete") for r in records),
            passed=sum(r["status"] == "complete" for r in records),
            failed=sum(r["status"] == "error" for r in records), records=records)
        all_records.extend(records)
    result = dict(generated_at=now(), stages=stages, complete=all(s["complete"] for s in stages.values()),
        scientific_review="Pending native-field independent audit, controls, complete figures, supplements and manuscript integration")
    write(run/"results/campaign_census.json", result)
    return {k: {kk: vv for kk, vv in v.items() if kk != "records"} for k, v in stages.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "freeze-B", "preflight", "worker", "compile-a2", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--basal-run", type=Path, default=basals.DEFAULT)
    parser.add_argument("--source-run", type=Path, default=DEFAULT)
    parser.add_argument("--stage", default="primary")
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    run = args.run.resolve()
    if args.action == "freeze":
        result = freeze_campaign(run, args.basal_run)
    elif args.action == "freeze-B":
        result = freeze_B(run, args.source_run.resolve())
    elif args.action == "worker":
        result = worker(run, args.stage, args.slot)
    elif args.action == "compile-a2":
        result = compile_a2(run)
    else:
        result = globals()[args.action](run)
    if result is not None:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
