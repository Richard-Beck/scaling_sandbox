"""Complete all97 reviewed active-recovery assays in a fresh frozen campaign.

Reuses exact15 pilot archives; missing82 cases use the identical native assay.
Pilot clocks are fixed per model and are never tuned per size.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import copy
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import all_five_incumbent_recovery_20261017 as assay

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/"data/runs/vision_linear_20261017/recovery"
PILOT = ROOT/"data/runs/all_five_incumbent_recovery_pilot_20261017"
TOPIC = ROOT/"figures/length_mechanisms_20261008/dev/vision_linear_20261017"
METRICS = ROOT/"data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv"
CANONICAL = ROOT/"data/runs/debelly_corrected_three_test_20261008"
sha, read, write, csvwrite = assay.sha, assay.read, assay.write, assay.csvwrite


def semantic_assay_check():
    original = ast.parse((PILOT/"frozen_sources/code/mechanisms/all_five_incumbent_recovery_20261017.py").read_text())
    current = ast.parse(Path(assay.__file__).read_text())
    names = ("initial_geometry", "kick", "integrate", "amounts", "pools", "amount_rates", "evaluate_case", "band_to_dense", "AuditedBanded")
    functions = lambda tree: {node.name: ast.dump(node, include_attributes=False) for node in tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
    a, b = functions(original), functions(current)
    for name in names: assert a[name] == b[name], name
    return list(names)


def cases():
    provenance = read(METRICS.parent/"report_provenance.json")
    assert sha(METRICS) == provenance["summary_csv_sha256"]
    metrics = list(csv.DictReader(METRICS.open())); assert len(metrics) == 97
    pilot = read(PILOT/"manifest.json")
    anchors = {(c["model"], c["size_native"]): c for c in pilot["cases"]}
    entries = []
    for model in assay.MODELS:
        for row in sorted([r for r in metrics if r["model"] == model], key=lambda r: float(r["domain_length"])):
            key = model, float(row["domain_length"])
            source = ROOT/row["checkpoint"]; assert sha(source) == row["checkpoint_sha256"]
            anchor = anchors.get(key)
            entries.append(dict(slot=len(entries), model=model, label=anchor["label"] if anchor else "size_series",
                                size_native=key[1], physical_length=float(row["physical_length"]),
                                core_width_native=float(row["old_pole_width"]), background_native=float(row["background_active"]),
                                source_checkpoint=row["checkpoint"], checkpoint_sha256=row["checkpoint_sha256"],
                                pilot_source_slot=anchor["slot"] if anchor else None))
    return entries


def preflight():
    entries = cases(); checked = []
    for case in entries:
        state, fields, x = assay.load_state(ROOT/case["source_checkpoint"], case["model"])
        geometry = assay.initial_geometry(case, state, fields, x)
        for epsilon in assay.EPSILONS: assay.kick(state, case["model"], geometry["direction"], epsilon)
        checked.append(dict(model=case["model"], size_native=case["size_native"],
                            core_sample_centers=geometry["core_sample_centers"],
                            core_has_fewer_than_8_centers=geometry["core_sample_centers"] < 8,
                            checkpoint_hash_and_kick_checks=True))
    return dict(generated_at=assay.now(), passed=True, cases=len(entries),
                pilot_reuse=sum(c["pilot_source_slot"] is not None for c in entries),
                fresh_cases=sum(c["pilot_source_slot"] is None for c in entries),
                unchanged_reviewed_scientific_functions=semantic_assay_check(),
                dense_linear_solve_fixture=assay.dense_wrapper_fixture(), checks=checked)


def initialize(run):
    if run.exists(): raise ValueError("Use a fresh recovery destination")
    check = preflight(); entries = cases()
    pilot = read(PILOT/"manifest.json"); clocks = read(PILOT/"results/clock_decision.json")
    for folder in ("inputs", "workers", "logs", "results", "frozen_sources"): (run/folder).mkdir(parents=True)
    sources = [Path(__file__), Path(assay.__file__), METRICS, METRICS.parent/"report_provenance.json",
               PILOT/"manifest.json", PILOT/"provenance.json", PILOT/"results/clock_decision.json",
               PILOT/"results/verification.json", CANONICAL/"manifest.json", TOPIC/"DESIGN.md"]
    for folder, pattern in (("code/core", "*.py"), ("code/core", "*.json"),
                            ("code/core/models", "*.py"), ("code/core/src/polarity1d", "*.py")):
        sources.extend(sorted((ROOT/folder).glob(pattern)))
    hashes = {}
    for source in sources:
        relative = source.relative_to(ROOT); target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target); hashes[str(relative)] = sha(target)
    imported = []
    for case in entries:
        source = ROOT/case["source_checkpoint"]; target = run/"inputs"/f'case_{case["slot"]:03d}.npz'
        shutil.copy2(source, target); assert sha(target) == case["checkpoint_sha256"]
        case["input"] = str(target.relative_to(run))
        if case["pilot_source_slot"] is not None:
            old = PILOT/"workers"/f'{case["pilot_source_slot"]:02d}'
            status_path = old/"status.json"
            status = read(status_path)
            if not status.get("complete"):
                status_path = old/"recovered_status.json"; status = read(status_path)
                assert status["original_status_sha256"] == sha(old/"status.json")
            assert status["complete"] and sha(old/"curves.npz") == status["curves_sha256"]
            assert status["case"]["checkpoint_sha256"] == case["checkpoint_sha256"]
            folder = run/"workers"/f'{case["slot"]:02d}'; folder.mkdir()
            shutil.copy2(old/"curves.npz", folder/"curves.npz")
            status = copy.deepcopy(status)
            status.update(case=case, imported_exact_pilot=True, original_pilot_status=str(status_path),
                          original_pilot_status_sha256=sha(status_path), original_pilot_archive=str(old/"curves.npz"),
                          original_pilot_archive_sha256=sha(old/"curves.npz"))
            write(folder/"status.json", status); imported.append(case["slot"])
    manifest = copy.deepcopy(pilot)
    manifest.update(generated_at=assay.now(), experiment="Full97 established-linear common active-deactivation recovery",
                    cases=entries, frozen_source_sha256=hashes, preflight=check,
                    fixed_model_horizons={m["model"]: m["selected_horizon_native"] for m in clocks["models"]},
                    imported_pilot_slots=imported,
                    stopping_rule="Complete and audit all97 at fixed reviewed clocks; report any new clock failure, never retune by size")
    write(run/"manifest.json", manifest); write(run/"provenance.json", dict(generated_at=assay.now(), manifest_sha256=sha(run/"manifest.json")))
    write(run/"results/status.json", dict(generated_at=assay.now(), producer="initialized", complete=False, cases=97, reused=len(imported), fresh=97-len(imported)))
    return manifest


def launch(run):
    manifest = assay.verify(run)
    if (run/"launch.json").exists(): raise ValueError("Already launched")
    missing = [c["slot"] for c in manifest["cases"] if c["slot"] not in manifest["imported_pilot_slots"]]
    frozen = run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command = shlex.join([sys.executable, str(frozen), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    options = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=00:30:00",
               "--job-name=vision-linear-recovery", "--array="+",".join(map(str, missing)), "--chdir="+str(run),
               "--output="+str(run/"logs/recovery_%A_%a.out"), "--error="+str(run/"logs/recovery_%A_%a.err")]
    job = subprocess.check_output(options+["--wrap", command], text=True).strip().split(";")[0]
    write(run/"launch.json", dict(generated_at=assay.now(), job_id=job, tasks=len(missing), reused=15, command=command))
    return job


def collect(run):
    manifest = assay.verify(run); curves, candidates, selected, statuses = [], [], [], []
    replacement_path = run/"results/early_window_replacement.json"
    replacement = read(replacement_path) if replacement_path.exists() else None
    if replacement:
        assert sha(Path(replacement["verification_path"])) == replacement["verification_sha256"]
        assert read(Path(replacement["verification_path"]))["passed"]
        assert sha(Path(replacement["fresh_review_path"])) == replacement["fresh_review_sha256"]
    for case in manifest["cases"]:
        amended = bool(replacement and case["slot"] == replacement["parent_slot"])
        folder = Path(replacement["archive_folder"]) if amended else run/"workers"/f'{case["slot"]:02d}'
        status = read(folder/"status.json")
        assert status["complete"] and status["status"] == "complete"
        assert sha(folder/"curves.npz") == status["curves_sha256"]
        assert status["case"]["checkpoint_sha256"] == case["checkpoint_sha256"]
        evaluation_manifest = read(Path(replacement["amendment_manifest"])) if amended else manifest
        with np.load(folder/"curves.npz") as z:
            a, b = assay.evaluate_case(case, z, evaluation_manifest)
            for row in b: row["measurement_status"] = "measured"
            if amended:
                for unmeasured in evaluation_manifest["unmeasured_later_candidates"]:
                    row = {key: None for key in b[0]}
                    row.update(model=case["model"], size_native=case["size_native"], label=case["label"],
                               horizon_native=unmeasured, measurement_status="unmeasured_after_required_window")
                    b.append(row)
            curves.extend(a); candidates.extend(b)
            horizon = manifest["fixed_model_horizons"][case["model"]]
            gate = next(r for r in b if r["horizon_native"] == horizon)
            detail = next(r for r in a if r["epsilon"] == .005 and abs(r["time_native"]-horizon) < 1e-14)
            initial_amount = float(z["tight_sham_core_amount"][0])
            selected.append(dict(model=case["model"], size_native=case["size_native"], physical_length=case["physical_length"],
                status="complete" if gate["passes"] else "clock_gate_unavailable", horizon_native=horizon,
                restored_fraction=gate["primary_F"], conditional_K_per_native_time=gate["primary_K"],
                empirical_E=gate["primary_E"], coefficient_available=gate["passes"],
                response_resolved_above_3E=gate["primary_change_resolved_above_3E"],
                error_pass=gate["error_pass"], linearity_pass=gate["linearity_pass"], saturation_pass=gate["saturation_pass"],
                prefix_error_max=gate["prefix_error_max"], prefix_saturation_max=gate["prefix_saturation_max"],
                initial_core_damage=detail["initial_core_damage"], initial_damage_fraction_of_core_active=detail["initial_core_damage"]/initial_amount,
                core_deficit=detail["core_deficit"], whole_active_deficit=detail["whole_active_deficit"], outside_core_deficit=detail["outside_core_deficit"],
                sham_core_change=detail["sham_core_change"], probe_core_change=detail["probe_core_change"],
                sham_core_net_rate=detail["sham_core_net_rate"], probe_core_net_rate=detail["probe_core_net_rate"],
                core_sample_centers=status["geometry"]["core_sample_centers"],
                checkpoint=case["source_checkpoint"], checkpoint_sha256=case["checkpoint_sha256"],
                archive=str((folder/"curves.npz").relative_to(run)), archive_sha256=status["curves_sha256"],
                exact_pilot_reuse=case["pilot_source_slot"] is not None,
                required_window_amendment=amended,
                later_candidates_unmeasured=amended))
        statuses.append(status)
    clock_audit = []
    for model in assay.MODELS:
        primary = [r for r in selected if r["model"] == model]
        qualified = all(r["coefficient_available"] for r in primary)
        alternatives = [h for h in manifest["debelly_candidates" if model == "debelly" else "continuous_candidates"]
                        if all(r["passes"] for r in candidates if r["model"] == model and r["horizon_native"] == h)]
        clock_audit.append(dict(model=model, fixed_horizon_native=manifest["fixed_model_horizons"][model],
                                all_size_clock_qualified=qualified, failed_sizes=[r["size_native"] for r in primary if not r["coefficient_available"]],
                                globally_passing_candidates=alternatives,
                                proposed_common_horizon_if_needed=max(alternatives) if not qualified and alternatives else None))
        for row in primary: row["model_clock_qualified"] = qualified
    csvwrite(run/"results/recovery_summary.csv", selected)
    csvwrite(run/"results/recovery_curves.csv", curves); csvwrite(run/"results/recovery_candidate_gates.csv", candidates)
    write(run/"results/clock_audit.json", dict(generated_at=assay.now(), models=clock_audit,
          policy="Fixed reviewed pilot clock, failures visible; any revised global clock requires explicit review"))
    write(run/"results/verification.json", dict(generated_at=assay.now(), passed=True, completed_cases=len(selected),
        reused_pilot_cases=15, newly_integrated_cases=82, collection_source_sha256=sha(Path(__file__)),
        required_window_amendment=replacement,
        tracked_pool_max_relative_error=max(v["max_relative_error"] for s in statuses for group in s["tracked_pool_checks"].values() for v in group),
        zero_kick_all_fields_identity=True, frozen_inputs_and_sources_verified=True,
        all_fixed_clocks_qualified=all(r["all_size_clock_qualified"] for r in clock_audit)))
    write(run/"results/status.json", dict(generated_at=assay.now(), producer="finished", complete=True, cases=97,
        all_fixed_clocks_qualified=all(r["all_size_clock_qualified"] for r in clock_audit)))
    return clock_audit


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("preflight", "initialize", "launch", "worker", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT); parser.add_argument("--slot", type=int)
    args = parser.parse_args(); run = args.run.resolve()
    if args.action == "preflight":
        result = preflight(); write(TOPIC/"supporting/metadata/RECOVERY_PREFLIGHT.json", result); print(json.dumps({k: v for k, v in result.items() if k != "checks"}, indent=2))
    elif args.action == "initialize": print(json.dumps({"cases": len(initialize(run)["cases"]), "run": str(run)}, indent=2))
    elif args.action == "launch": print(launch(run))
    elif args.action == "worker": assay.worker(run, args.slot)
    elif args.action == "collect": print(json.dumps(collect(run), indent=2))


if __name__ == "__main__": main()
