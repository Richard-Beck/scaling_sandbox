"""Scientific arm-selection checks: companion deduplication and no E4 fallback."""
from pathlib import Path
import sys
import gzip
import json

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_regions_stage_20261020 as stage_module
from signature_figure4_regions_stage_20261020 import choose_arm, collect_audit
from signature_figure4_regions_campaign_20261020 import save, sha


def case(identifier, size, f, horizon=32000., kind="noise", sigma=1.):
    return dict(id=identifier, model="wave_pinning", C_label=size, f=f,
        kind=kind, sigma=sigma, seed=123 if kind == "noise" else None,
        native_horizon=horizon, profile_cadence=5., origin="E2",
        scalar_sha256=identifier, profile_sha256=identifier)


def test_common_keeps_dense_transition_inputs_and_new_sizes_without_duplicate_old_seeds():
    shared = case("old", 20., .001)
    dense = case("dense", 20., .02)
    intermediate = case("intermediate", 30., .001)
    rows, _ = choose_arm({"E2": [shared, dense], "E3": [dict(shared), intermediate]}, "common")
    assert {r["id"] for r in rows} == {"old", "dense", "intermediate"}
    with pytest.raises(ValueError, match="exactly"):
        choose_arm({"E3": [intermediate]}, "common")


def test_longest_replaces_selected_conditions_whole_and_keeps_needed_control_horizons():
    common = [case(f"short{i}", 20., .001 + i * .001) for i in range(53)]
    common += [case("zero_short", 20., 0., kind="zero", sigma=0.)]
    long = [case(f"long{i}", 20., .001 + i * .001, 128000.) for i in range(52)]
    long += [case("zero_long", 20., 0., 128000., kind="zero", sigma=0.)]
    rows, _ = choose_arm({"E2": common, "E3": [], "E4": long}, "longest")
    noise = [r for r in rows if r["kind"] == "noise"]
    assert len(noise) == 53 and sum(r["native_horizon"] == 128000. for r in noise) == 52
    assert {r["id"] for r in noise if r["native_horizon"] == 32000.} == {"short52"}
    assert {r["native_horizon"] for r in rows if r["kind"] == "zero"} == {32000., 128000.}
    # Missing a selected arm remains an error, never a short-run substitution.
    with pytest.raises(ValueError, match="52"):
        choose_arm({"E2": common, "E3": [], "E4": long[1:]}, "longest")


def test_sigma_arms_keep_exact_paired_seed_as_distinct_forcing_conditions():
    rows = [case(str(sigma), 20., .03, sigma=sigma) for sigma in (.25, .5, 1., 2.)]
    selected, _ = choose_arm({"E1": rows}, "E1")
    assert len(selected) == 4 and {r["sigma"] for r in selected} == {.25, .5, 1., 2.}
    with pytest.raises(ValueError, match="sigma-one"):
        choose_arm({"E2": rows, "E3": []}, "common")


def test_same_observed_frames_have_identical_archive_bytes_across_analysis_arms(tmp_path):
    value=dict(threshold=.5, result=dict(frames=[dict(time=16000., contrast=1.1234567890123457)]))
    first=tmp_path/"common/observed.json.gz"
    second=tmp_path/"longest/another_filename.json.gz"
    save(first,value);save(second,value)
    assert sha(first)==sha(second)
    with gzip.open(first,"rt") as stream:
        assert json.load(stream)==value
    assert first.read_bytes()[4:8]==b"\0\0\0\0"


def test_independent_master_binds_inventory_all_cases_and_measured_receipts(tmp_path):
    inventory=tmp_path/"inventory.json";save(inventory,dict(cases=[dict(id="native_case")]))
    observed=tmp_path/"observed.json.gz";save(observed,dict(result={}))
    measured=tmp_path/"workers/native_case/result.json"
    save(measured,dict(primary=dict(path=str(observed),sha256=sha(observed),downsample_factor=16)))
    source=Path(stage_module.__file__)
    save(tmp_path/"inputs/worker_freeze.json",dict(inventory_sha256=sha(inventory),source_sha256={
        source.name:sha(source),"signature_figure4_regions_audit_20261020.py":"independent_auditor"}))
    record=dict(case_id="native_case",passed=True,measurement_sha256=sha(measured),
        primary_sha256=sha(observed),views_checked=25,all_view_native_frames=100,frames=10)
    for shard in range(24):
        save(tmp_path/f"results/independent_all_views/shard_{shard:03d}.json",dict(
            passed=True,inventory_sha256=sha(inventory),source_sha256="independent_auditor",
            shard=shard,shards=24,cases=int(shard==0),records=[record] if shard==0 else []))
    assert collect_audit(tmp_path)["cases"]==1
    projection=json.loads((tmp_path/"results/independent_primary_native_frames_audit.json").read_text())
    assert projection["subset_of_all_views"] and projection["inventory_sha256"]==sha(inventory)
    shard_path=tmp_path/"results/independent_all_views/shard_000.json"
    broken=json.loads(shard_path.read_text());broken["inventory_sha256"]="unrelated_previous_run"
    save(shard_path,broken)
    with pytest.raises(ValueError,match="identity or canonical case coverage"):
        collect_audit(tmp_path)
