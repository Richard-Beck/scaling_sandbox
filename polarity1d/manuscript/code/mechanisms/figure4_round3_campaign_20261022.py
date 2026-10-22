"""Compact-only Figure 4 round-3 campaign, with unchanged native equations.

``prepare`` freezes qualified basals, original native sources and storage checks.
``submit`` launches one CPU per independent seed through SLURM. Native fields
exist only in worker memory; the durable output is the declared post-burn-in
float32 descriptor archive and seeds 0/1's bounded float16 illustration.
"""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import time
import traceback

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
sys.dont_write_bytecode = True
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "polarity1d/manuscript").is_dir())
MS = ROOT / "polarity1d/manuscript"
CODE = MS / "code/mechanisms"
DEFAULT = MS / "data/runs/signature_figure4_round3_compact_20261022"
PARENT = MS / "data/runs/revision_noise_strength_replication_scoped_20261019"
INTERMEDIATE = MS / "data/runs/signature_figure4_E3_intermediate_20261020"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS = ("wave_pinning", "holmes_model3", "spring_ring", "goryachev")
FRACTIONS = (.0001, .0003, .001, .003, .005, .01, .02, .03, .05, .1, .3, 1.)
EXTRA = (.004, .007, .015, .025)
REFERENCES = dict(wave_pinning=.1, holmes_model3=1., spring_ring=4., goryachev=6.4)
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def allocated(path):
    """One traversal, hard-linked inode counted once, without following links."""
    seen, total = set(), 0
    for base, _, files in os.walk(path, followlinks=False):
        for name in files:
            p = Path(base) / name
            st = p.lstat()
            key = (st.st_dev, st.st_ino)
            if key not in seen:
                seen.add(key)
                total += st.st_blocks * 512
    return total


def update_storage(out):
    """Maintain one delimited additive record without replacing older audits."""
    out = Path(out).resolve()
    generated = now()
    status = read(out / "status.json")
    prelaunch = read(out / "results/storage_prelaunch.json")
    allocated_bytes = allocated(out)
    retirement_path = ROOT / "reports/figure4_round3_retirement_20261022/execution.json"
    retired = retirement_path.exists() and read(retirement_path).get("complete", False)
    retirement_statement = ("The scoped old Figure4 noise native profiles were deleted after the COMPLETE audit; the separate retirement execution ledger records actual reclaimed storage."
        if retired else "Old Figure4 noise native arrays remain preserved; deletion is gated on the new campaign's audited COMPLETE state.")
    report = dict(generated_utc=generated, campaign=str(out.relative_to(ROOT)),
        allocated_bytes=allocated_bytes, allocated_GiB=allocated_bytes / 2**30,
        status=status.get("status"), complete=status.get("complete", False),
        complete_seeds=status.get("native_seeds_complete", 0),
        initial_seeds=24960, maximum_seeds=99840,
        maximum_estimated_bytes=prelaunch["maximum_estimated_bytes"],
        reserved_bytes=prelaunch["reserved_bytes"],
        shared_filesystem_available_bytes_at_probe=prelaunch["filesystem_available_bytes"],
        personal_quota_available_bytes=None,
        personal_quota_status=prelaunch["personal_quota_status"],
        immutable_manifest_sha256=sha(out / "manifest.json"),
        old_noise_native_deleted=retired,
        measurement="Allocated512byte blocks; hard-linked inode counted once; active-producer traversal, not an atomic snapshot.")
    report_path = ROOT / "reports/figure4_round3_storage_20261022.json"
    write(report_path, report)
    marker = "figure4-round3-compact-storage-20261022"
    begin, end = f"<!-- {marker} begin -->", f"<!-- {marker} end -->"
    record = f"""{begin}

## Figure 4 round-3 compact campaign, 22 October 2026

Measured **{generated}**: **{allocated_bytes:,} allocated bytes
({allocated_bytes / 2**30:.6f} GiB)** under
`{out.relative_to(ROOT)}/`. Status **{status.get('status')}**;
{status.get('native_seeds_complete', 0):,} compact seeds complete, 24,960 initial
and at most 99,840 after precision-driven 96-seed additions. This active-producer
measurement does not replace the older whole-workspace totals above.
[Targeted storage and quota record](reports/figure4_round3_storage_20261022.json).

The scientifically authoritative new records retain post-burn-in float32
activity moments and P = M1/M0, contrast and exact round-1 tracked-region
descriptors, plus bounded float16 illustrations for seeds 0/1 only. Native
profiles exist in memory during extraction and are never written to durable
storage. Initial planning estimated 26.586 decimal GB, maximum 103.823 GB,
with a 207.647 GB reservation (2× margin); the shared filesystem exposed
{prelaunch['filesystem_available_bytes']:,} free bytes at its dated probe.
Personal quota headroom could not be read and is explicitly unknown.
The initial small-QoS arrays were upgraded to xxlarge for independent-seed
throughput; `results/qos_burst_upgrade.json` records the effective scheduler
setting and capacity evidence, while the immutable manifest preserves the
original submission resources. {retirement_statement}

{end}"""
    path = ROOT / "DATA_MANIFEST.md"
    content = path.read_text()
    if begin in content:
        a, b = content.index(begin), content.index(end) + len(end)
        content = content[:a] + record + content[b:]
    else:
        content = content.rstrip() + "\n\n" + record + "\n"
    path.write_text(content)
    retention = f"""{begin}

## Figure 4 round-3 retention exception, authorized 22 October 2026

The user authorized replacement of OLD Figure 4 noise native HPC outputs
**only after the new compact campaign is COMPLETE**. Its current state at
{generated} is **{status.get('status')}**, with
{status.get('native_seeds_complete', 0):,} verified compact seeds.
{retirement_statement} Preserve all preceding native arrays while any
allocated seed or final precision assessment is missing. This exception does
not cover Figures 2/3, exact qualified basals, source/protocol/audit records,
or repository round-2 compact CSVs and figures.

Keep the new `{out.relative_to(ROOT)}/` compact NPZ files, immutable
manifest/frozen sources, exact 20 basal inputs, per-seed integrity and numerical
receipts, CSV index, 2000-resample whole-seed precision reports, allocation
decisions and actual SLURM receipts. Seeds 0/1 retain float16 examples bounded
at 400 spatial bins and 2000 frames. No full native profile history is retained
for new seeds. Unresolvable rates remain explicitly undefined; finite
SE(log k)>0.07 triggers batches of 96 to at most384. Optional largest/weak
extensions require a recorded condition-specific allocation and remain
separate from the original compact evidence.

[Current allocated storage record](reports/figure4_round3_storage_20261022.json).
The immutable initial manifest says small QoS; effective xxlarge throughput
resources are recorded in `results/qos_burst_upgrade.json`. The measured shared
filesystem free space is distinct from unknown personal quota headroom.

{end}"""
    path = ROOT / "RETENTION.md"
    content = path.read_text()
    if begin in content:
        a, b = content.index(begin), content.index(end) + len(end)
        content = content[:a] + retention + content[b:]
    else:
        content = content.rstrip() + "\n\n" + retention + "\n"
    path.write_text(content)
    return report


def storage_check(out, seed_cap=384):
    """The shared-filesystem free value is explicitly distinct from a quota."""
    rows = []
    for command in (["df", "-B1", str(MS)], ["quota", "-sv"], ["quota", "-vs", "-u", os.environ.get("USER", str(os.getuid()))]):
        try:
            p = subprocess.run(command, text=True, capture_output=True, timeout=30)
            rows.append(dict(command=shlex.join(command), returncode=p.returncode,
                             stdout=p.stdout, stderr=p.stderr))
        except (OSError, subprocess.TimeoutExpired) as e:
            rows.append(dict(command=shlex.join(command), error=str(e)))
    fs = os.statvfs(MS)
    available = fs.f_bavail * fs.f_frsize
    # Post-burn-in 3201 frames. Eight tracked components/frame is a deliberately
    # conservative planning scenario, not a scientific cap on observed regions.
    per_seed = 3201 * (18 * 4 + 8 * (6 * 4 + 4 + 1) + 8) + 32 * 1024
    conditions = 4 * 5 * len(FRACTIONS) + 5 * len(EXTRA)
    exemplar = conditions * 2 * (400 * 2000 * 2 + 2000 * 8)
    expected = conditions * seed_cap * per_seed + exemplar
    reserved = int(2 * expected)
    result = dict(generated_utc=now(), probes=rows,
        filesystem_available_bytes=available, personal_quota_available_bytes=None,
        personal_quota_status="Quota lookup did not expose an enforceable personal limit; filesystem availability is shared, not personal quota headroom.",
        planning_regions_per_frame=8, planning_bytes_per_seed=per_seed,
        initial_seeds=96, maximum_seeds=seed_cap, conditions=conditions,
        initial_estimated_bytes=conditions * 96 * per_seed + exemplar,
        maximum_estimated_bytes=expected, reserved_bytes=reserved,
        margin_factor=2., available_over_reservation=available / reserved,
        passed_filesystem_margin=available >= reserved,
        limitation="Planning estimate, not an upper bound; measured allocation and live free space are checked before each submission. No profile histories are written to durable storage.")
    if out:
        write(Path(out) / "results/storage_prelaunch.json", result)
    return result


def prepare(out):
    out = Path(out).resolve()
    if out.exists() or not out.is_relative_to(MS / "data/runs"):
        raise ValueError("Choose a fresh directory under manuscript/data/runs")
    parent, intermediate = read(PARENT / "manifest.json"), read(INTERMEDIATE / "manifest.json")
    basals = {(b["model"], float(b["size"])): (PARENT, b) for b in parent["basals"]}
    basals.update({(b["model"], float(b["size"])): (INTERMEDIATE, b) for b in intermediate["basals"]})
    for directory in ("inputs", "frozen_sources", "workers", "logs", "results", "submissions"):
        (out / directory).mkdir(parents=True)
    sources, inputs = {}, {}
    for relative, digest in parent["source_sha256"].items():
        source = PARENT / relative
        if sha(source) != digest:
            raise ValueError("Original native source changed: " + relative)
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[relative] = digest
    for name in (Path(__file__).name, "signature_figure4_campaign_20261020.py", "signature_figure4_regions_20261020.py"):
        source = CODE / name
        target = out / "frozen_sources" / name
        shutil.copy2(source, target)
        sources[str(target.relative_to(out))] = sha(target)
    brief = ROOT / "docs/tasks/briefs/20261022_figure4_round3.txt"
    shutil.copy2(brief, out / "inputs" / brief.name)
    inputs["inputs/" + brief.name] = sha(brief)
    conditions = []
    basal_records = []
    for model in MODELS:
        sizes = (1., 1.5, 2., 2.5, 3.) if model == "spring_ring" else (20., 30., 40., 50., 60.)
        for size in sizes:
            source_root, original = basals[model, size]
            b = copy.deepcopy(original)
            source = source_root / b["input"]
            if sha(source) != b["input_sha256"]:
                raise ValueError("Qualified basal changed")
            target = out / "inputs" / f"{model}_C{size:g}.npz"
            shutil.copy2(source, target)
            b.update(input=str(target.relative_to(out)), source=str(source), source_manifest=str(source_root / "manifest.json"), source_manifest_sha256=sha(source_root / "manifest.json"))
            basal_records.append(b)
            inputs[b["input"]] = sha(target)
            for fraction in sorted(FRACTIONS + (EXTRA if model == "goryachev" else ())):
                conditions.append(dict(condition_index=len(conditions), condition_id=f"{model}:C{size:g}:f{fraction:g}",
                    model=model, size=size, fraction=fraction, mean=REFERENCES[model] * fraction,
                    sigma=1., kind="noise", platform="spring" if model == "spring_ring" else "chemical",
                    cells=b["cells"], input=b["input"], input_sha256=b["input_sha256"],
                    basal_circumference=b["basal_circumference"]))
    if len(conditions) != 260:
        raise AssertionError("Incorrect full design")
    protocols = copy.deepcopy(parent["native_protocols"])
    manifest = dict(created_utc=now(), campaign="Figure4_round3_compact", schema_version=1,
        conditions=conditions, basals=basal_records, native_protocols=protocols,
        initial_seeds=96, batch_size=96, maximum_seeds=384, precision_se_log_k=.07,
        bootstrap_resamples=2000, burn_in_fraction=.5, saved_cadence=dict(chemical=5., spring=.5),
        optional_extension=dict(maximum_horizon_factor=4, largest_sizes_only=True, weak_fraction_max=.003,
            requires_unresolvable_rate=True, automatically_submitted=False),
        resources=dict(partition="red", qos="small", cpus=1, memory="4G", walltime="02:00:00",
            justification="Independent single-core native integrations; small QoS suits throughput. Scheduler CPU limits control concurrency; no array throttle."),
        source_sha256=sources, input_sha256=inputs,
        compact_schema=dict(
            float32=["time", "circumference", "M0", "M1_real", "M1_imag", "P_real", "P_imag", "P_abs", "R", "contrast", "contrast_ratio", "dominant_heading", "dominant_centroid_physical", "region_heading", "region_centroid_physical", "region_amount", "region_width"],
            integer=["region_count", "dominant_track_id", "region_offset", "region_id"],
            boolean=["region_amount_available", "low_contrast"],
            P="M1/M0; endpoint eligibility uses P_abs>=0.1", R="abs(M1)/abs(M1_initial), legacy diagnostic only",
            region_offset="CSR offsets: region data for frame i use slice offset[i]:offset[i+1]",
            missing="NaN float, track id -1. Undefined region amounts remain NaN; never zero.",
            native_profiles_retained=False, postburn_only=True),
        exemplar_schema=dict(replicates=[0, 1], dtype="float16", maximum_spatial_bins=400,
            maximum_frames=2000, scope="Post-burn-in primary activity; periodic conservative bin average and selected exact native observation times. Illustration only."),
        historical_native_deletion=dict(authorized_after_new_compact_COMPLETE=True, performed=False,
            scope="OLD Figure4 noise native HPC arrays only; preserve round2 compact CSVs/figures and exact basals/source records"))
    write(out / "manifest.json", manifest)
    write(out / "freeze.json", dict(manifest_sha256=sha(out / "manifest.json")))
    write(out / "results/allocation.json", dict(generated_utc=now(), targets={c["condition_id"]: 96 for c in conditions}, precision_evaluated=False, extensions={}))
    check = storage_check(out)
    write(out / "status.json", dict(generated_utc=now(), status="PREPARED", complete=False,
        native_seeds_complete=0, native_seeds_initial=24960, storage_margin_passed=check["passed_filesystem_margin"]))
    worker = out / "worker.sh"
    frozen = out / "frozen_sources" / Path(__file__).name
    worker.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join([PYTHON, "-B", str(frozen), "worker", "--out", str(out)]) + ' --slot "$SLURM_ARRAY_TASK_ID" --batch "$1"\n')
    worker.chmod(0o750)
    return dict(out=str(out), conditions=260, initial_seeds=24960, maximum_seeds=99840,
                manifest_sha256=sha(out / "manifest.json"), storage=check)


def verify(out):
    out = Path(out).resolve()
    manifest = read(out / "manifest.json")
    if sha(out / "manifest.json") != read(out / "freeze.json")["manifest_sha256"]:
        raise ValueError("Immutable manifest changed")
    for mapping in (manifest["source_sha256"], manifest["input_sha256"]):
        for relative, digest in mapping.items():
            if sha(out / relative) != digest:
                raise ValueError("Frozen source/input changed: " + relative)
    return manifest


def native_module(out, platform):
    helper = load_module("f4_r3_unchanged_native_helper", out / "frozen_sources/signature_figure4_campaign_20261020.py")
    return helper.native_module(out, platform), helper


def extract(time_values, profiles, circumference, basal_profile):
    """Use the exact round-1 contour detector and maximum-overlap tracker."""
    import numpy as np
    detector = load_module("f4_r3_exact_round1_regions", Path(__file__).with_name("signature_figure4_regions_20261020.py"))
    t, p = np.asarray(time_values), np.asarray(profiles)
    C = np.broadcast_to(np.asarray(circumference), t.shape)
    basal = np.asarray(basal_profile)
    tracking = detector.track_profiles(t, p, C, basal_contrast=float(np.ptp(basal)), turning_intervals=())
    angle = 2 * np.pi * (np.arange(p.shape[1]) + .5) / p.shape[1]
    M0, M1 = C * p.mean(axis=1), C * (p @ np.exp(1j * angle)) / p.shape[1]
    initial_moment = abs(np.mean(basal * np.exp(1j * angle)))
    # Basal circumference factor cancels only for fixed rings. Caller supplies
    # actual basal circumference when replacing the legacy R array below.
    P = np.divide(M1, M0, out=np.full_like(M1, np.nan + 1j * np.nan), where=M0 > 0)
    fields = dict(time=t, circumference=C, M0=M0, M1_real=M1.real, M1_imag=M1.imag,
        P_real=P.real, P_imag=P.imag, P_abs=abs(P), R=abs(M1) / (C[0] * initial_moment),
        contrast=np.ptp(p, axis=1), contrast_ratio=np.ptp(p, axis=1) / np.ptp(basal))
    frames = tracking["frames"]
    fields.update(region_count=np.array([f["region_count"] for f in frames], dtype="int16"),
        dominant_track_id=np.array([f["dominant_track_id"] if f["dominant_track_id"] is not None else -1 for f in frames], dtype="int32"),
        dominant_heading=np.array([f["dominant_heading"] if f["dominant_heading"] is not None else np.nan for f in frames]),
        low_contrast=np.array([f["low_contrast"] for f in frames], dtype=bool))
    fields["dominant_centroid_physical"] = fields["dominant_heading"] * C / (2 * np.pi)
    rows = [r for f in frames for r in f["regions"]]
    fields["region_offset"] = np.r_[0, np.cumsum(fields["region_count"], dtype="int64")]
    fields["region_id"] = np.array([r["track_id"] for r in rows], dtype="int32")
    fields["region_heading"] = np.array([r["heading"] if r["heading"] is not None else np.nan for r in rows])
    fields["region_centroid_physical"] = np.array([r["heading"] * f["circumference"] / (2 * np.pi) if r["heading"] is not None else np.nan for f in frames for r in f["regions"]])
    fields["region_amount"] = np.array([r["excess"] if r["amount_available"] else np.nan for r in rows])
    fields["region_width"] = np.array([r["width"] for r in rows])
    fields["region_amount_available"] = np.array([r["amount_available"] for r in rows], dtype=bool)
    for key, value in list(fields.items()):
        if np.asarray(value).dtype.kind in "fc":
            fields[key] = np.asarray(value, dtype="float32")
    return fields, tracking["protocol"], tracking["summary"]


def exemplar(fields, profiles):
    import numpy as np
    p = np.asarray(profiles)
    frames = np.unique(np.linspace(0, len(p) - 1, min(2000, len(p)), dtype=int))
    bins = max(i for i in range(1, min(400, p.shape[1]) + 1) if p.shape[1] % i == 0)
    # Equal native-cell groups keep a uniform physical grid without interpolation.
    groups = np.array_split(np.arange(p.shape[1]), bins)
    reduced = np.stack([p[frames][:, group].mean(axis=1) for group in groups], axis=1)
    positions = np.array([(g[0] + g[-1] + 1) / (2 * p.shape[1]) for g in groups], dtype="float32")
    return dict(time=fields["time"][frames], circumference=fields["circumference"][frames],
        primary=reduced.astype("float16"), material_positions=positions,
        native_spatial_cells=np.array(p.shape[1], dtype="int32"),
        native_frame_indices=frames.astype("int32"))


def solve(out, condition, seed, horizon_factor=1, kind="noise", sigma=1.):
    """Native BDF and exact forcing constructors; profile samples stay in RAM."""
    import numpy as np
    from scipy.integrate import BDF
    manifest = read(out / "manifest.json")
    spec = copy.deepcopy(manifest["native_protocols"][condition["platform"]])
    original_horizon = spec["horizon"]
    spec["horizon"] *= horizon_factor
    native, helper = native_module(out, condition["platform"])
    with np.load(out / condition["input"], allow_pickle=False) as z:
        initial = z["state"].copy()
    n = condition["cells"]
    t = np.arange(int(round(spec["horizon"] / spec["scalar_interval"])) + 1) * spec["scalar_interval"]
    deadline = time.monotonic() + (28000 if horizon_factor > 1 else 6900)
    noise = dict(spec["noise"], sigma=float(sigma))
    checks = dict(minimum_concentration=float(initial[:-1].min()) if condition["platform"] == "spring" else float(initial.min()), maximum_pool_relative_error=0.)
    if condition["platform"] == "chemical":
        grid = native.ring.RingGrid(condition["size"], n)
        model = native.ring.native_model(condition["model"])
        forcing = native.NoiseInput(grid, spec["horizon"], seed, condition["mean"], kind=kind, **noise)
        pool = native.ring.pool_amounts(model, grid, initial)
        def inspect(y):
            if not np.isfinite(y).all():
                raise ValueError("Nonfinite native state")
            checks["minimum_concentration"] = min(checks["minimum_concentration"], float(y.min()))
            drift = float(np.max(abs(native.ring.pool_amounts(model, grid, y) / pool - 1)))
            checks["maximum_pool_relative_error"] = max(checks["maximum_pool_relative_error"], drift)
            if checks["minimum_concentration"] < -spec["diagnostics"]["allowed_negative_roundoff"] or drift > spec["diagnostics"]["pool_relative_tolerance"]:
                raise ValueError("Native positivity/pool guard failed")
        rhs = lambda tt, yy: native.ring.native_rhs(model, grid, yy.reshape(initial.shape), forcing(tt)).ravel()
        factory = BDF if horizon_factor == 1 else lambda *a, **kw: helper.OriginalCapBDF(BDF, original_horizon, *a, **kw)
        solver = factory(rhs, 0., initial.ravel(), spec["horizon"],
            jac_sparsity=native.ring.periodic_sparsity(len(model.state_names), n), **spec["solver"])
        profiles = np.empty((len(t), n), dtype="float64")
        profiles[0] = initial[0]
        wanted = 1
        while solver.status == "running":
            if STOP or time.monotonic() > deadline:
                raise TimeoutError("SLURM signal/worker deadline")
            message = solver.step()
            if solver.status == "failed":
                raise RuntimeError(message)
            inspect(solver.y.reshape(initial.shape))
            dense = solver.dense_output()
            while wanted < len(t) and t[wanted] <= solver.t:
                y = dense(t[wanted]).reshape(initial.shape)
                inspect(y)
                profiles[wanted] = y[0]
                wanted += 1
        if wanted != len(t):
            raise RuntimeError("Native solve missed required observation times")
        C = np.full(len(t), condition["size"])
        basal = initial[0]
        checks.update(nfev=int(solver.nfev), njev=int(solver.njev), nlu=int(solver.nlu))
    else:
        closure = native.Closure(condition["size"], n)
        forcing = native.MaterialNoise(n, initial[-1], spec["horizon"], seed, condition["mean"], kind=kind,
            knot=noise["knot"], coupling=noise["coupling"], tau=noise["tau"], ell=noise["ell"], sigma=float(sigma), bridge_base=noise.get("bridge_base"))
        def rhs(tt, y):
            if STOP or time.monotonic() > deadline:
                raise TimeoutError("SLURM signal/worker deadline")
            return closure.rhs(y, forcing(tt))
        solution = helper.spring_native_solve(rhs, lambda tt, y: closure.jacobian(y, forcing(tt)), initial, t, spec["solver"], original_horizon if horizon_factor > 1 else None)
        if not solution.success or len(solution.t) != len(t):
            raise RuntimeError(solution.message)
        y = solution.y.T
        if not np.isfinite(y).all() or y[:, :-1].min() < -5e-10 or y[:, -1].min() <= 0:
            raise ValueError("Native Spring positivity/finite guard failed")
        amount = y[:, -1] * np.mean(y[:, :n] + y[:, n:-1], axis=1)
        drift = float(np.max(abs(amount / amount[0] - 1)))
        if drift > 1e-7:
            raise ValueError("Native Spring physical-pool drift")
        profiles, C, basal = y[:, :n].copy(), y[:, -1].copy(), initial[:n]
        checks.update(minimum_concentration=float(y[:, :-1].min()), maximum_pool_relative_error=drift,
            nfev=int(solution.nfev), njev=int(solution.njev), nlu=int(solution.nlu))
    # Avoid the old audit's large extra midpoint/field arrays: preserve exact
    # latent and kernel hashes plus streaming knot moments.
    knots = forcing.knots
    h = hashlib.sha256(np.asarray(knots, dtype="<f8", order="C").tobytes()).hexdigest() if knots is not None else None
    audit = dict(seed=seed, gaussian_knots_sha256=h, shape=list(knots.shape) if knots is not None else None, sigma=float(sigma),
        mean=condition["mean"], native_reference=REFERENCES[condition["model"]],
        kind="lognormal" if kind == "noise" else kind, formula="mean*exp(sigma*Z-sigma**2/2)" if kind == "noise" else "constant mean", finite_realization_normalization=False,
        gaussian_mean=float(knots.mean()) if knots is not None else None, gaussian_variance=float(knots.var()) if knots is not None else None,
        original_native_protocol=spec, native_kinetics_retuned=False)
    keep = t >= .5 * spec["horizon"]
    return t[keep], profiles[keep], C[keep], basal, checks, audit


def worker(out, slot, batch=0, horizon_factor=1):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Native integration requires a SLURM allocation")
    import numpy as np
    out = Path(out).resolve()
    manifest = verify(out)
    ci, within = divmod(slot, 96)
    replicate = batch * 96 + within
    if not 0 <= ci < len(manifest["conditions"]) or not 0 <= replicate < 384:
        raise ValueError("Seed slot is outside the declared design")
    c = manifest["conditions"][ci]
    if horizon_factor != 1:
        largest = 3. if c["platform"] == "spring" else 60.
        if c["size"] != largest or c["fraction"] > manifest["optional_extension"]["weak_fraction_max"]:
            raise ValueError("Extensions are restricted to largest weak-input conditions")
        allocation = read(out / "results/allocation.json")
        approved = allocation.get("extensions", {}).get(c["condition_id"])
        if approved != horizon_factor:
            raise ValueError("Unresolvable-rate extension requires a recorded allocation decision")
    seed = 202610220000000 + ci * 1000 + replicate
    folder = out / "workers" / f"c{ci:03d}" / f"seed{replicate:03d}"
    if horizon_factor != 1:
        folder = folder / f"extension_x{horizon_factor}"
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "status.json").exists() and read(folder / "status.json").get("complete"):
        return dict(status="already_complete", condition_id=c["condition_id"], replicate=replicate)
    started = time.monotonic()
    status = dict(generated_utc=now(), status="RUNNING", complete=False, condition=c, replicate=replicate,
        seed=seed, horizon_factor=horizon_factor, manifest_sha256=sha(out / "manifest.json"), slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(folder / "status.json", status)
    def stop(*_):
        global STOP
        STOP = True
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGUSR1, stop)
    try:
        t, p, C, basal, checks, audit = solve(out, c, seed, horizon_factor)
        fields, protocol, summary = extract(t, p, C, basal)
        # Legacy R uses the actual qualified basal C, not postburn C[0].
        theta = 2 * np.pi * (np.arange(len(basal)) + .5) / len(basal)
        initial_M1 = c["basal_circumference"] * abs(np.mean(basal * np.exp(1j * theta)))
        fields["R"] = (np.hypot(fields["M1_real"], fields["M1_imag"]) / initial_M1).astype("float32")
        temporary = folder / "compact.tmp.npz"
        np.savez_compressed(temporary, **fields)
        with np.load(temporary, allow_pickle=False) as z:
            if len(z["time"]) != len(t) or z["region_offset"][-1] != len(z["region_id"]):
                raise ValueError("Compact round-trip shape validation failed")
            np.testing.assert_allclose(z["P_abs"], np.hypot(z["P_real"], z["P_imag"]), rtol=3e-6, atol=1e-8)
        temporary.replace(folder / "compact.npz")
        if replicate in (0, 1):
            np.savez_compressed(folder / "exemplar.npz", **exemplar(fields, p))
        write(folder / "noise_audit.json", audit)
        write(folder / "region_audit.json", dict(protocol=protocol, summary=summary,
            centroid="Round-1 excess-weighted contour centroid; physical = heading*circumference/(2*pi)",
            no_full_native_profiles_retained=True))
        status.update(status="COMPLETE", complete=True, validation=checks,
            native_horizon=float(2 * t[0]), analysis_start=float(t[0]), analysis_stop=float(t[-1]),
            compact_sha256=sha(folder / "compact.npz"), exemplar_sha256=sha(folder / "exemplar.npz") if (folder / "exemplar.npz").exists() else None,
            allocated_bytes=allocated(folder))
    except Exception as error:
        status.update(status="FAILED", error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
        raise
    finally:
        status.update(generated_utc=now(), wall_seconds=time.monotonic() - started,
            peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(folder / "status.json", status)
    return dict(condition_id=c["condition_id"], replicate=replicate, status=status["status"], wall_seconds=status["wall_seconds"])


def submit(out, batch=0, precision=None, dry_run=False):
    out = Path(out).resolve()
    m = verify(out)
    state = read(out / "results/allocation.json")
    if not 0 <= batch < 4:
        raise ValueError("Only four 96-seed batches are authorized")
    selected = m["conditions"]
    if batch:
        if precision is None:
            raise ValueError("Adaptive batches require whole-seed precision results")
        report = read(precision)
        rows = {r["condition_id"]: r for r in report["conditions"]}
        selected = []
        for c in m["conditions"]:
            r = rows.get(c["condition_id"])
            if r is None:
                raise ValueError("Precision report omitted a designed condition")
            if r["n_seeds"] < batch * 96:
                continue
            se = r.get("se_log_k")
            uncertain = r.get("reason") == "bootstrap_valid_fraction_below_0.975"
            if (se is not None and math.isfinite(se) and se > .07) or uncertain:
                selected.append(c)
                state["targets"][c["condition_id"]] = (batch + 1) * 96
        state.update(precision_evaluated=True, precision_source=str(Path(precision).resolve()), precision_sha256=sha(precision), generated_utc=now())
    check = storage_check(out)
    if not check["passed_filesystem_margin"]:
        raise RuntimeError("Insufficient filesystem margin for conservative compact campaign reservation")
    records = []
    for model in MODELS:
        chosen = [c for c in selected if c["model"] == model]
        if not chosen:
            continue
        indices = ",".join(f'{c["condition_index"] * 96}-{(c["condition_index"] + 1) * 96 - 1}' for c in chosen)
        receipt_path = out / "submissions" / f"batch{batch}_{model}.json"
        if not dry_run and receipt_path.exists() and read(receipt_path).get("submitted"):
            records.append(read(receipt_path))
            continue
        command = ["sbatch", "--parsable", "--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=1", "--mem=4G", "--time=02:00:00", "--signal=B:USR1@120", "--array=" + indices,
            f"--job-name=F4r3-{model[:3]}-b{batch}", "--output=" + str(out / "logs/%A_%a.out"),
            "--error=" + str(out / "logs/%A_%a.err"), str(out / "worker.sh"), str(batch)]
        receipt = dict(generated_utc=now(), batch=batch, model=model, seeds=len(chosen) * 96, command=shlex.join(command), dry_run=dry_run,
            qos_reason="Burst throughput: independent1CPU4G2h jobs, xxlarge3000CPU/user allowance; red69nodes/4416CPUs and about1TB RAM/node. User permits parallel compute. Initial small-QoS submissions were upgraded separately in qos_burst_upgrade.json.")
        if not dry_run:
            try:
                p = subprocess.run(command, text=True, capture_output=True, timeout=45)
                receipt.update(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr, submitted=p.returncode == 0)
            except (OSError, subprocess.TimeoutExpired) as e:
                receipt.update(returncode=None, error=str(e), submitted=False)
        records.append(receipt)
        write(receipt_path.with_name(receipt_path.stem + "_dryrun.json") if dry_run else receipt_path, receipt)
    if not dry_run:
        write(out / "results/allocation.json", state)
    write(out / "results" / ("submission_dryrun_summary.json" if dry_run else "submission_summary.json"), dict(generated_utc=now(), records=records,
        submitted_seeds=sum(r["seeds"] for r in records if r.get("submitted")),
        richard_command=shlex.join([PYTHON, "-B", str(CODE / Path(__file__).name), "submit", "--out", str(out), "--batch", str(batch)]) + (" --precision " + shlex.quote(str(precision)) if precision else "")))
    return records


def collect(out, precision=None):
    out = Path(out).resolve()
    m = verify(out)
    allocation = read(out / "results/allocation.json")
    if precision is None and allocation.get("precision_source"):
        prior = Path(allocation["precision_source"])
        if prior.exists():
            precision = prior
    rows, failed, missing, running = [], [], 0, 0
    for c in m["conditions"]:
        ci = c["condition_index"]
        for replicate in range(allocation["targets"][c["condition_id"]]):
            folder = out / "workers" / f"c{ci:03d}" / f"seed{replicate:03d}"
            extension = allocation.get("extensions", {}).get(c["condition_id"], 1)
            if extension != 1:
                folder = folder / f"extension_x{extension}"
            path = folder / "status.json"
            if not path.exists():
                missing += 1
                continue
            s = read(path)
            if not s.get("complete"):
                if s.get("status") == "RUNNING":
                    running += 1
                    continue
                failed.append(dict(condition_id=c["condition_id"], replicate=replicate, status=s["status"], error=s.get("error")))
                continue
            compact = folder / "compact.npz"
            if sha(compact) != s["compact_sha256"]:
                raise ValueError("Compact seed archive changed")
            rows.append(dict(condition_id=c["condition_id"], model=c["model"], size=c["size"], fraction=c["fraction"],
                replicate=replicate, seed=s["seed"], compact=str(compact),
                exemplar=str(folder / "exemplar.npz") if (folder / "exemplar.npz").exists() else "",
                status="COMPLETE", horizon=s["native_horizon"], manifest_sha256=sha(out / "manifest.json"), compact_sha256=s["compact_sha256"]))
    path = out / "results/compact_index.csv"
    fieldnames = ["condition_id", "model", "size", "fraction", "replicate", "seed", "compact", "exemplar", "status", "horizon", "manifest_sha256", "compact_sha256"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    converged = False
    unresolved, undefined = [], []
    if precision:
        report = read(precision)
        by_id = {r["condition_id"]: r for r in report["conditions"]}
        for c in m["conditions"]:
            r = by_id.get(c["condition_id"], {})
            target = allocation["targets"][c["condition_id"]]
            if r.get("se_log_k") is None:
                undefined.append(dict(condition_id=c["condition_id"], reason=r.get("reason", "No resolvable rate at a declared lag")))
            needs_precision = (r.get("se_log_k") is not None and r["se_log_k"] > .07) or r.get("reason") == "bootstrap_valid_fraction_below_0.975"
            if r.get("n_seeds", 0) != target or (target < 384 and needs_precision):
                unresolved.append(c["condition_id"])
        converged = len(by_id) == 260 and not unresolved
        allocation.update(precision_evaluated=True, precision_source=str(Path(precision).resolve()), precision_sha256=sha(precision))
        write(out / "results/allocation.json", allocation)
    complete = not missing and not failed and not running and converged
    retirement_execution = ROOT / "reports/figure4_round3_retirement_20261022/execution.json"
    old_deleted = retirement_execution.exists() and read(retirement_execution).get("complete", False)
    status = dict(generated_utc=now(), status="COMPLETE" if complete else "INCOMPLETE", complete=complete,
        native_seeds_complete=len(rows), native_seeds_expected=sum(allocation["targets"].values()),
        missing_seeds=missing, running_seeds=running, failed_seeds=failed, precision_qualified=converged,
        unresolved_precision_conditions=unresolved, allocated_bytes=allocated(out),
        undefined_rate_conditions=undefined,
        native_profiles_retained=False, old_native_deletion_permitted=complete,
        old_native_deletion_performed=old_deleted,
        qualification="COMPLETE requires all allocated seeds plus 2000-resample precision assessment: finite SE(logk)<=.07,384seed cap, or an explicitly undefined rate. Undefined lag/rates are omitted rather than imputed. No old native deletion occurs here.")
    write(out / "status.json", status)
    return status


def extend(out, precision, condition_ids=(), dry_run=False):
    """Optional explicit4x reruns of unresolved largest weak conditions.

    Original compact archives remain beside the extension. The collector uses
    the selected extension for every allocated seed, never a mixed horizon.
    """
    out = Path(out).resolve()
    m = verify(out)
    state = read(out / "results/allocation.json")
    rows = {r["condition_id"]: r for r in read(precision)["conditions"]}
    selected = []
    for c in m["conditions"]:
        r = rows.get(c["condition_id"], {})
        largest = 3. if c["platform"] == "spring" else 60.
        requested = not condition_ids or c["condition_id"] in condition_ids
        if requested and c["size"] == largest and c["fraction"] <= .003 and r.get("se_log_k") is None and r.get("n_seeds", 0) >= 96:
            selected.append(c)
    if condition_ids and {c["condition_id"] for c in selected} != set(condition_ids):
        raise ValueError("Requested extensions include a condition without an unresolved largest/weak rate")
    check = storage_check(out)
    additional = sum(state["targets"][c["condition_id"]] * check["planning_bytes_per_seed"] * 4 for c in selected)
    if check["filesystem_available_bytes"] < 2 * additional + check["reserved_bytes"]:
        raise RuntimeError("Insufficient free-space margin for optional4x compact reruns")
    if not dry_run:
        for c in selected:
            state["extensions"][c["condition_id"]] = 4
        state.update(extension_precision_source=str(Path(precision).resolve()), extension_precision_sha256=sha(precision), generated_utc=now())
        write(out / "results/allocation.json", state)
    records = []
    for c in selected:
        for batch in range(state["targets"][c["condition_id"]] // 96):
            ci = c["condition_index"]
            script = out / "submissions" / f"extension_c{ci:03d}_b{batch}.sh"
            script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join([PYTHON, "-B", str(out / "frozen_sources" / Path(__file__).name), "worker", "--out", str(out), "--batch", str(batch), "--horizon-factor", "4"]) + ' --slot "$SLURM_ARRAY_TASK_ID"\n')
            script.chmod(0o750)
            command = ["sbatch", "--parsable", "--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=1", "--mem=8G", "--time=08:00:00", f"--array={ci*96}-{(ci+1)*96-1}",
                "--job-name=F4r3-extend", "--output=" + str(out / "logs/extension_%A_%a.out"), "--error=" + str(out / "logs/extension_%A_%a.err"), str(script)]
            record = dict(condition_id=c["condition_id"], batch=batch, horizon_factor=4, command=shlex.join(command), dry_run=dry_run)
            if not dry_run:
                p = subprocess.run(command, text=True, capture_output=True, timeout=45)
                record.update(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr, submitted=p.returncode == 0)
            records.append(record)
    write(out / "results" / ("extension_dryrun.json" if dry_run else "extension_submission.json"), dict(generated_utc=now(), records=records, additional_estimated_bytes=additional,
        scope="Optional4x independent full-horizon reruns with exact latent prefix and native solver restart at the original cap; original compact archives preserved."))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("prepare", "storage", "update-storage", "verify", "collect", "submit", "extend", "worker"):
        p = sub.add_parser(action)
        p.add_argument("--out", type=Path, default=DEFAULT)
        if action in ("submit", "collect", "extend"):
            p.add_argument("--precision", type=Path)
        if action in ("worker", "submit"):
            p.add_argument("--batch", type=int, default=0)
        if action in ("submit", "extend"):
            p.add_argument("--dry-run", action="store_true")
        if action == "extend":
            p.add_argument("--condition-id", action="append", default=[])
        if action == "worker":
            p.add_argument("--slot", type=int, required=True)
            p.add_argument("--horizon-factor", type=int, choices=(1, 2, 4), default=1)
    args = parser.parse_args()
    if args.action == "prepare":
        result = prepare(args.out)
    elif args.action == "storage":
        result = storage_check(args.out)
    elif args.action == "update-storage":
        result = update_storage(args.out)
    elif args.action == "verify":
        m = verify(args.out)
        result = dict(passed=True, conditions=len(m["conditions"]), manifest_sha256=sha(args.out / "manifest.json"))
    elif args.action == "worker":
        result = worker(args.out, args.slot, args.batch, args.horizon_factor)
    elif args.action == "submit":
        result = submit(args.out, args.batch, args.precision, args.dry_run)
    elif args.action == "extend":
        if not args.precision:
            raise ValueError("Optional extensions require final precision results")
        result = extend(args.out, args.precision, args.condition_id, args.dry_run)
    else:
        result = collect(args.out, args.precision)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
