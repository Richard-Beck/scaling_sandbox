"""Render the Spring ring development reader from retained native evidence.

This file does not integrate a model or change scientific inputs. The HTML is a
self-contained snapshot: compact raw profiles are embedded, and scheduler state
is queried while rendering rather than inferred from an old status file.
"""
from pathlib import Path
from datetime import datetime, timezone
from html import escape
from io import StringIO
import hashlib
import json
import os
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
PREP = ROOT / "data/runs/spring_ring_preparation_pilot_20261017"
DIST = ROOT / "data/runs/spring_ring_distance_pilot_20261017"
RUN = ROOT / "data/runs/spring_ring_validation_20261019"
OUT = ROOT / "figures/length_mechanisms_20261008/dev/spring_ring_validation_20261019"
EXAMPLE = '''"""Small restart example; run from the repository root inside SLURM.

Output is an API example, not a new frozen scientific campaign.
"""
from pathlib import Path
import argparse
import sys
import numpy as np

root = Path.cwd() / "polarity1d/manuscript"
sys.path.insert(0, str(root / "code/mechanisms"))
from spring_ring_model_20261019 import Closure, CuePhase, run_protocol

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=Path("/tmp/spring_ring_example.npz"))
args = parser.parse_args()
if args.output.exists():
    raise FileExistsError("Choose a fresh --output path")
model = Closure(resting_circumference=2.0, cells=200)
with np.load(root / "data/runs/spring_ring_preparation_pilot_20261017/workers/qualified_basal.npz") as saved:
    state = saved["state"].copy()
incumbent = model.material_positions[np.argmax(state[:model.cells])]
result = run_protocol(model, [
    CuePhase(duration=0.1, peak=4.0, center_fraction=incumbent + 0.4,
             total_support=0.2, output_interval=0.1),
    CuePhase(duration=13000.0, output_interval=5.0),
], initial_state=state)
args.output.parent.mkdir(parents=True, exist_ok=True)
np.savez_compressed(args.output, time=result["time"], states=result["states"],
                    physical_amount=result["physical_amount"],
                    polarity_magnitude=result["polarity_magnitude"])
print("Saved full states:", args.output)
'''


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def link(path, label):
    return '<a href="' + escape(os.path.relpath(path, OUT)) + '">' + escape(label) + '</a>'


def scheduler_snapshot():
    launch = RUN / "launch.json"
    if not launch.exists():
        return {"state": "No launch record available", "queried": False}
    record = read(launch)
    job_id = record.get("array_job_id", record.get("slurm_job_id", record.get("job_id")))
    if not job_id:
        return {"state": "Launch recorded; scheduler identifier unavailable", "queried": False}
    job_id = str(job_id).split(";")[0]
    try:
        result = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%i|%T"],
                                capture_output=True, text=True, timeout=20, check=True)
        rows = [line.strip().split("|") for line in result.stdout.splitlines() if line.strip()]
        if rows:
            return {"state": "; ".join(" ".join(row) for row in rows), "queried": True, "job_id": job_id,
                    "live": True}
        result = subprocess.run(["sacct", "-X", "-n", "-P", "-j", job_id,
                                 "--format=JobIDRaw,State"],
                                capture_output=True, text=True, timeout=20, check=True)
        records = [line.split('|') for line in result.stdout.splitlines() if line.strip()]
        counts = {}
        for row in records:
            state = row[1] if len(row) > 1 else 'unknown'
            counts[state] = counts.get(state, 0) + 1
        description = ', '.join(str(count) + ' ' + state for state, count in sorted(counts.items()))
        return {"state": description or "No scheduler record returned", "queried": True,
                "job_id": job_id, "live": False, "records": records}
    except (OSError, subprocess.SubprocessError) as error:
        return {"state": "Scheduler observation unavailable: " + str(error), "queried": False, "job_id": job_id}


def profile_data():
    manifest = read(DIST / "manifest.json")
    summary = read(DIST / "results/summary.json")
    assert read(DIST / "results/audit.json")["passed"]
    assert read(PREP / "results/audit.json")["passed"]
    with np.load(DIST / "inputs/basal.npz") as archive:
        initial = archive["state"].copy()
    assert sha(DIST / "inputs/basal.npz") == manifest["input_sha256"]
    n = manifest["cells"]
    shift = n // 4 - manifest["incumbent_cell"]
    sources = [DIST / "manifest.json", DIST / "results/summary.json", DIST / "results/audit.json",
               DIST / "inputs/basal.npz", PREP / "results/audit.json", PREP / "workers/status.json"]
    cases = []
    for job in manifest["jobs"]:
        if job["reflection_control"]:
            continue
        record_path = DIST / "workers" / f'{job["slot"]:02d}' / "status.json"
        record = read(record_path)
        assert record["complete"] and record["status"] == "complete"
        sources.append(record_path)
        times, states = [], []
        elapsed = 0.0
        for phase in record["phases"]:
            path = DIST / phase["path"]
            assert sha(path) == phase["sha256"]
            sources.append(path)
            with np.load(path) as archive:
                t, y = archive["time"].copy(), archive["states"].copy()
            np.testing.assert_array_equal(y[0], initial if not states else states[-1])
            take = np.arange(len(t))
            if not times:
                take = take[(t <= 10.0000001) | np.isin(t, [20., 50., 100., 1000.])]
            else:
                take = take[np.isin(t, [.1, 1., 10., 100., 1000., 13000.])]
            times.extend((t[take] + elapsed).tolist())
            states.extend(y[take])
            elapsed += float(t[-1])
        selected = np.asarray(states)
        assert len(times) == len(selected) and np.all(np.diff(times) > 0)
        target = np.asarray(job["target_region_cells"])
        row = next(x for x in summary["jobs"] if x["offset_fraction"] == job["offset_fraction"])
        cases.append({"offset": job["offset_fraction"], "distance": job["shortest_initial_arc"],
                      "target": ((job["center_cell"] + shift) % n + .5) / n,
                      "times": times, "profiles": np.round(np.roll(selected[:, :n], shift, axis=1), 9).tolist(),
                      "circumference": selected[:, -1].tolist(),
                      "target_initial_mean": float(initial[target].mean()),
                      "target_maximum_at_first_sample": float(selected[1, target].max()),
                      "two_components_last": row["sampled_two_half_contrast_component_range"][-1],
                      "two_peaks_last": row["sampled_two_resolved_peak_range"][-1],
                      "final_profile_residual": row["final_profile_max_residual_vs_translated_initial"]})
    return {"cells": n, "initial": np.roll(initial[:n], shift).tolist(),
            "incumbent": (n // 4 + .5) / n, "initial_circumference": float(initial[-1]),
            "threshold": manifest["original_half_contrast_threshold"], "cue_width": .2,
            "cases": cases}, sources


def sci(value):
    if value is None:
        return 'pending'
    if value == 0:
        return '0'
    coefficient, exponent = f'{float(value):.2e}'.split('e')
    return coefficient + ' × 10<sup>' + str(int(exponent)) + '</sup>'


def equation_svg(formula, label, name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['svg.fonttype'] = 'path'
    fig = plt.figure(figsize=(7, .65))
    fig.text(.5, .5, formula, ha='center', va='center', fontsize=19)
    buffer = StringIO()
    fig.savefig(buffer, format='svg', bbox_inches='tight', pad_inches=.08, transparent=True,
                metadata={'Date': None})
    plt.close(fig)
    svg = buffer.getvalue()
    path = OUT / 'supporting/vectors' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg)
    svg = svg[svg.index('<svg'):]
    return svg.replace('<svg ', '<svg role="img" aria-label="' + escape(label) + '" ', 1)


def validation_content(summary, scheduler, audit=None):
    if summary is None:
        return ('<p>New validation results have not yet been collected. The retained preparation and '
                'three-position result below are available; their successful audits do not establish '
                'spatial convergence.</p><p class="muted">Scheduler at reader generation: '
                + escape(scheduler["state"]) + '</p>')
    status = summary.get("status", summary.get("producer", "Summary available"))
    complete = summary.get("complete", False)
    items = []
    for key in ("verdict", "interpretation", "conclusion", "scope"):
        if isinstance(summary.get(key), str):
            items.append('<p>' + escape(summary[key]) + '</p>')
    items.insert(0, '<p><strong>' + escape(str(status)) + '</strong> · '
                 + ('collected output complete' if complete else 'partial output') + '</p>')
    jobs = summary.get('jobs', [])
    finished = [job for job in jobs if job.get('complete')]
    items.append('<p>' + str(len(finished)) + ' of ' + str(len(jobs))
                 + ' workers have finished. '
                 + ('Independent saved-state audit passes.' if audit and audit.get('passed') else
                    'Final independent saved-state audit is pending.') + '</p>')
    fresh = next((job for job in jobs if job.get('job', {}).get('name') == 'fresh_200'), {})
    basal = next((phase for phase in fresh.get('phases', []) if phase['kind'] == 'basal'), None)
    if basal:
        items.append('<p>A fresh preparation reaches a single basal peak at C = '
                     + f'{basal["final"]["circumference"]:.6f}'
                     + '. Its maximum scaled full-field/circumference range in the last 1,000 units is '
                     + sci(basal['final_1000_scaled_range']) + '.</p>')
    comparisons = summary.get('comparisons', {})
    chosen = [('grid_200_vs_400', '200 versus 400 cells'),
              ('fresh_tolerance', '100× tighter tolerances: fresh preparation'),
              ('retained_radau_vs_retained', 'Radau versus BDF: retained checkpoint'),
              ('rotation_73_vs_retained', 'Shift incumbent and cue by 73 cells')]
    table = []
    for key, label in chosen:
        phases = comparisons.get(key, {})
        largest = max((v['max_chemical_field_difference'] for v in phases.values()), default=None)
        final = max((v['final_chemical_field_difference'] for k, v in phases.items()
                     if k in ('basal', 'recovery')), default=None)
        table.append('<tr><td>' + escape(label) + '</td><td>'
                     + sci(largest) + '</td><td>' + sci(final) + '</td></tr>')
    items.append('<div class="tablewrap"><table class="choices"><thead><tr><th>Bounded check</th>'
                 '<th>Largest sampled field difference</th><th>Final cue-free field difference</th></tr></thead>'
                 '<tbody>' + ''.join(table) + '</tbody></table></div>')
    items.append('<p class="caption">Differences use complete active and inactive fields in raw concentration units. '
                 'The first numerical column takes the largest difference across all compared phases; '
                 'the final column takes the larger difference at the available basal/recovery endpoints. '
                 'Fresh preparations are independent at each mesh, followed by the same offset-0.4 challenge. '
                 'The full record retains every phase and circumference difference. '
                 'Two meshes test sensitivity at this point; they do not prove spatial convergence. '
                 'The fresh tolerance row shows bounded transient sensitivity despite matching final fields.</p>')
    perturb = next((job for job in jobs if job.get('job', {}).get('name') == 'conserved_perturbation'), {})
    if perturb.get('complete'):
        items.append('<p>One small amount-conserving chemical perturbation returns to a single stationary '
                     'basal cap. Its final activity direction shifts by '
                     + f'{perturb["final"]["vector_direction_fraction"]:.3g}'
                     + ' of a turn. This one perturbation tests recovery of the prepared pattern; it '
                     'does not establish stability against every disturbance.</p>')
    pulses = [job for name in ('finite_0p1', 'finite_10') for job in jobs
              if job.get('job', {}).get('name') == name and job.get('complete')]
    if pulses:
        items.append('<h3>Withdrawing a short cue changes the retained site</h3>'
                     '<p>The same incumbent is challenged at offset 0.4, with the same peak 4 and '
                     'physical support 0.2. Only pulse duration changes; the final state is measured '
                     '13,000 native units after withdrawal.</p>')
        rows = []
        for job in pulses:
            final = job['final']
            rows.append('<tr>' + ''.join('<td>' + escape(value) + '</td>' for value in (
                f'{job["job"]["cue_duration"]:g}', str(final['resolved_peaks']),
                f'{final["vector_direction_fraction"]:.5f}', f'{final["vector_magnitude"]:.5f}',
                'passes' if job.get('tail_stationary') else 'fails')) + '</tr>')
        items.append('<div class="tablewrap"><table class="choices"><thead><tr><th>Cue duration</th>'
                     '<th>Final maxima</th><th>Direction: turns from original</th>'
                     '<th>Moment magnitude</th><th>Tail stationarity</th></tr></thead><tbody>'
                     + ''.join(rows) + '</tbody></table></div>')
        items.append('<p class="caption">The original incumbent is at 0 turns and the cue at 0.4. '
                     'The 0.1-unit pulse leaves a cap near the original site; the 10-unit pulse leaves '
                     'a cap near the challenger. Small phase shifts remain, so neither is an exact '
                     'copy of a cap pinned to the specified centre. These two doses demonstrate a '
                     'withdrawal-dependent contrast; they do not locate a minimum successful duration. '
                     + ('The final independent audit covers these trajectories.' if audit and audit.get('passed') else
                        'Values are preliminary worker output until the independent audit is available.') + '</p>')
    items.append('<p class="muted">Summary generated: ' + escape(str(summary.get("generated_at", "unrecorded")))
                 + '. Scheduler at reader generation: ' + escape(scheduler["state"]) + '</p>')
    items.append('<details><summary>Full validation record</summary><pre>'
                 + escape(json.dumps(summary, indent=2)) + '</pre></details>')
    return ''.join(items)


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Spring on a ring — development evidence</title><style>
:root{color-scheme:light;--ink:#19354b;--blue:#137c9b;--orange:#c46220;--green:#246d55;--muted:#526774;--paper:#fff;--line:#d9e3e8}
*{box-sizing:border-box}body{margin:0;background:#eef3f5;color:var(--ink);font:16px/1.55 system-ui,sans-serif}main{max-width:1080px;margin:auto;background:var(--paper);padding:36px 44px 55px}h1{font-size:2rem;line-height:1.2;margin:8px 0 14px}h2{font-size:1.32rem;line-height:1.3;margin:34px 0 12px}h3{font-size:1rem;margin:0 0 5px}p{margin:9px 0 15px}a{color:#086c91;text-decoration-thickness:1px;text-underline-offset:3px}small,.muted{color:var(--muted);font-size:.88rem}.eyebrow{font-size:.78rem;letter-spacing:.09em;text-transform:uppercase;color:var(--blue)}.lead{font-size:1.12rem;max-width:900px}.badge{display:inline-block;padding:4px 10px;border-radius:14px;background:#e5f2ee;color:var(--green);font-size:.8rem}.notice{padding:14px 18px;border-left:4px solid var(--orange);background:#fff6ec;margin:20px 0}.geometry{background:#f5f9fb;border:1px solid var(--line);border-radius:8px;padding:16px;display:grid;grid-template-columns:1fr 1.3fr;gap:28px}.geometry svg{width:100%;height:155px}.geometry p{font-size:.91rem;margin-bottom:0}.equation{display:grid;grid-template-columns:1fr auto;align-items:center;margin:24px 0;font-size:1.15rem;text-align:center}.equation math{margin:auto;font-size:1.15rem}.equation span{font-size:.9rem}.choices{width:100%;border-collapse:collapse;margin:15px 0;font-size:.92rem}th{text-align:left;font-weight:600;background:#f1f6f8}td,th{padding:10px 12px;border-bottom:1px solid var(--line);vertical-align:top}.tablewrap{overflow-x:auto}.viewer{border:1px solid var(--line);border-radius:9px;padding:20px;margin-top:18px}.controls{display:flex;align-items:center;gap:13px;flex-wrap:wrap}.controls label{font-size:.92rem}select,button{font:inherit;border:1px solid #bacbd4;background:white;color:var(--ink);border-radius:5px;padding:5px 9px}button{cursor:pointer}button:hover{background:#e8f3f6}input[type=range]{accent-color:var(--blue);width:min(570px,100%)}.slider{display:flex;align-items:center;gap:14px;margin:15px 0 6px}.slider output{min-width:165px;font-size:.9rem}.plots{display:grid;grid-template-columns:minmax(0,1fr) 235px;gap:8px;align-items:center}#profile{width:100%;height:auto;overflow:visible}#ring{width:220px;height:220px}.readouts{display:flex;gap:22px;flex-wrap:wrap;background:#f4f8fa;border-radius:5px;padding:10px 12px;font-size:.89rem}.legend{display:flex;gap:18px;flex-wrap:wrap;font-size:.83rem;color:var(--muted);margin:7px 0}.key{display:inline-block;width:22px;height:3px;vertical-align:middle;margin-right:5px}.caption{font-size:.88rem;color:var(--muted);margin-top:12px}.tick{fill:#526774;font-size:12px}.axis{stroke:#9eb2bd;stroke-width:1}.grid{stroke:#e4ecef;stroke-width:1}details{margin:15px 0}summary{cursor:pointer;color:var(--blue)}pre{white-space:pre-wrap;word-break:break-word;background:#f5f8fa;padding:14px;font-size:.8rem;max-height:420px;overflow:auto}.sources{font-size:.88rem}.sources a{display:inline-block;margin:3px 14px 3px 0}footer{border-top:1px solid var(--line);margin-top:32px;padding-top:15px}.statusbox{border:1px solid var(--line);border-radius:7px;padding:14px 18px}.statusbox p:last-child{margin-bottom:0}@media(max-width:720px){main{padding:24px 18px}.geometry{grid-template-columns:1fr;gap:22px}.plots{grid-template-columns:1fr}#ring{margin:auto;display:block}.viewer{padding:12px}h1{font-size:1.7rem}.choices{font-size:.82rem}td,th{padding:8px}}
</style><style>.equation{grid-template-columns:minmax(0,1fr) auto}.equation svg{width:100%;max-width:600px;height:auto;display:block;margin:auto}.slider input{min-width:0;flex:1;width:0}.slider output{min-width:145px}@media(max-width:720px){.equation{font-size:.9rem}.slider output{min-width:135px;font-size:.8rem}}</style></head><body><main>
<div class="eyebrow">Spring · periodic domain · development reader</div>
<h1>A polarized Spring analogue works on a ring</h1>
<p class="lead">Periodic chemical fields coupled to a moving circumference retain one active Rac region after cue withdrawal. A new challenger appears promptly at each tested position; the near case bridges to the incumbent, while farther sites coexist as the old peak decays.</p>
<span class="badge">Retained preparation and distance audits pass</span>
<p class="muted">Reader generated @@GENERATED@@. Retained trajectories generated 17 October 2026. This page is a generated snapshot; reload after regeneration to view new results. No automatic browser refresh.</p>
<div class="notice"><strong>The modeling choice is explicit.</strong> Removing the original endpoints requires a new mechanical rule. This reader evaluates the <em>mean-force global-extension analogue</em>. On 19 October the user explicitly agreed: “Yes—carry over the existing stretching, dilution and diffusion terms.” The source comparison below explains that carryover. Force averaging remains the declared circular adaptation.</div>

<h2>What changes when the interval closes?</h2>
<div class="geometry">
<div><h3>Original interval</h3><svg viewBox="0 0 370 155" role="img" aria-label="An interval with forces acting at its two endpoints">
<defs><marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#c46220"/></marker></defs>
<line x1="65" y1="77" x2="305" y2="77" stroke="#137c9b" stroke-width="9"/><circle cx="65" cy="77" r="9" fill="#19354b"/><circle cx="305" cy="77" r="9" fill="#19354b"/>
<line x1="60" y1="77" x2="17" y2="77" stroke="#c46220" stroke-width="3" marker-end="url(#arrow)"/><line x1="310" y1="77" x2="353" y2="77" stroke="#c46220" stroke-width="3" marker-end="url(#arrow)"/>
<text x="64" y="42" text-anchor="middle" fill="#19354b" font-size="17">F(u₋)</text><text x="305" y="42" text-anchor="middle" fill="#19354b" font-size="17">F(u₊)</text><text x="185" y="116" text-anchor="middle" fill="#526774" font-size="15">two moving ends · length L(t)</text></svg>
<p>Endpoint Rac sets two outward forces. Their sum drives extension against the spring. Chemical flux vanishes at each end.</p></div>
<div><h3>Declared circular analogue</h3><svg viewBox="0 0 440 155" role="img" aria-label="A periodic chemical ring whose distributed mean force controls one global circumference">
<circle cx="107" cy="78" r="51" fill="none" stroke="#137c9b" stroke-width="9"/>
<g stroke="#c46220" stroke-width="2.5" marker-end="url(#arrow)"><path d="M107 22v-17"/><path d="M107 134v17"/><path d="M51 78H31"/><path d="M163 78h20"/><path d="M67 38L53 24"/><path d="M147 38l14-14"/><path d="M67 118l-14 14"/><path d="M147 118l14 14"/></g>
<text x="107" y="76" text-anchor="middle" font-size="16" fill="#19354b">C(t)</text><text x="107" y="97" text-anchor="middle" font-size="12" fill="#526774">dynamic</text>
<text x="214" y="53" font-size="15" fill="#19354b">Periodic active + inactive Rac</text><text x="214" y="80" font-size="15" fill="#19354b">Mean of local force switches</text><text x="214" y="107" font-size="15" fill="#19354b">One uniform extension</text></svg>
<p>Every position has the same mechanical status. Spatially averaged forcing drives a single circumference; the diagram represents this scalar closure, with no local deformation or force-propagation field.</p></div>
</div>
<div class="equation">@@CLOSURE_EQUATION@@<span>(1)</span></div>
<p>The average is taken <em>after</em> the nonlinear Rac force switch. Local reaction, diffusion, force-switch and extension-feedback coefficients retain their native values. Concentrations dilute as the ring expands. Periodic transport conserves the physical amount of active plus inactive Rac.</p>
<div class="tablewrap"><table class="choices"><thead><tr><th>Choice</th><th>Meaning for this result</th></tr></thead><tbody>
<tr><td>One dynamic circumference</td><td>All trials share C = 2.04418 initially; subsequent circumference may differ. Extension still changes global Rac deactivation and dilution.</td></tr>
<tr><td>Resting circumference C₀ = 2</td><td>A native-unit feasibility point. It supplies no micron conversion or unique interval-length-to-perimeter mapping.</td></tr>
<tr><td>Periodic cue</td><td>Peak 4; full physical support 0.2; centre fixed in material coordinates. It is a full circular patch, with a different total dose from an endpoint half-patch.</td></tr>
<tr><td>Global mechanical feedback</td><td>Rearranging an unchanged Rac distribution preserves its mean force. Separation can alter chemical interactions and thereby the global extension; the closure contains no local mechanical transmission time.</td></tr>
</tbody></table></div>

<details open><summary>Why dilution is retained: it is already the interval convention</summary>
<p>The native interval adapter computes the fractional expansion rate from the difference of its endpoint velocities, divided by the current length. It applies that same rate to <em>every</em> active and inactive concentration through the terms −u L̇/L and −v L̇/L, and scales both material-coordinate diffusion terms by 1/L². The equations therefore use uniform affine deformation even when only one physical endpoint moves.</p>
<p>The ring carries over those existing terms with circumference C in place of length L: −u Ċ/C, −v Ċ/C and diffusion coefficients D/C². Its scalar expansion law changes from endpoint sampling to the mean local force. Adding periodic transport without those dilution terms would change the native geometric convention and would fail to conserve physical Rac amount on the expanding domain.</p>
<p>In this preparation C rises from 2 to 2.04418, an increase of 2.21%. Expansion alone would reduce a concentration with fixed species amount by 2.16%; the actual active/inactive split also changes through reactions. Material-coordinate diffusion coefficients decrease by 4.28%, while the native extension-dependent Rac deactivation coefficient increases from 3 to 6.534. These are distinct geometric and feedback effects; the physical diffusivities remain unchanged.</p>
<p>During the fresh long-cue challenge at offset 0.4, C instead falls from 2.04418 to 2.04157. Its net geometric contribution is a 0.128% concentration increase, followed by the inverse expansion during recovery. “Dilution” therefore includes signed concentration changes: contraction concentrates material while expansion dilutes it. @@GEOMETRY_LINKS@@</p>
<p class="caption">Source: @@DILUTION_SOURCE@@. The preserved dilution is not an added ring-only mechanism. A separate diagnostic can measure its contribution along the coupled trajectory; these percentages are endpoint geometry comparisons, not a causal partition of polarization.</p>
</details>

<h2>See the two-site competition in the raw fields</h2>
<p>Each trial restarts exactly the same qualified incumbent. Its original half-contrast core occupies 95 of 200 cells; the nearer cue starts outside that core but samples a higher Rac flank. Changing location changes both arc lengths and the initial local chemistry.</p>
<div class="viewer">
<div class="controls"><label for="case">Challenger offset</label><select id="case"><option value="0">0.3 C — near</option><option value="1">0.4 C — intermediate</option><option value="2">0.5 C — opposite</option></select><button id="previous" type="button" aria-label="Previous sampled time">←</button><button id="next" type="button" aria-label="Next sampled time">→</button><button id="at02" type="button">t = 0.2</button><button id="at5" type="button">t = 5</button><button id="final" type="button">Final basal</button></div>
<div class="slider"><input id="time" type="range" min="0" max="100" value="2" step="1" aria-label="Saved time index"><output id="clock" for="time"></output></div>
<div class="plots"><svg id="profile" viewBox="0 0 700 330" role="img" aria-label="Raw active Rac concentration across the connected periodic domain"></svg><svg id="ring" viewBox="0 0 235 235" role="img" aria-label="The same raw activity shown around a connected ring"></svg></div>
<div class="legend"><span><i class="key" style="background:#137c9b"></i>active Rac u(q,t)</span><span><i class="key" style="background:#aab8c0"></i>same initial incumbent</span><span><i class="key" style="background:#c46220"></i>cue location/support</span><span>Dashed: frozen initial half-contrast level</span></div>
<div class="readouts" aria-live="polite"><span id="size"></span><span id="phase"></span><span id="moment"></span></div>
<p class="caption">Raw concentration uses one common scale. Both plot edges connect. Display coordinates shift the same initial incumbent to q ≈ 0.25; simulation states remain unchanged. The small ring shows identical values on a common color scale. The circular activity moment measures activity organization, with no cell-velocity rule. The time slider traverses saved samples, densely every 0.1 unit through t = 10 and sparsely thereafter.</p>
<details><summary>Activity moment definition</summary><p class="caption">The first circular moment uses raw active Rac, including its basal background. Activity direction is arg(P), with arbitrary material-angle origin; magnitude is |P|. Direction is unavailable at vanishing magnitude. Multiple separated sites can reduce magnitude by cancellation, so the full profile remains necessary. The drawn arrow has fixed length and shows direction only; the readout reports magnitude.</p><div class="equation">@@MOMENT_EQUATION@@<span>(2)</span></div></details>
</div>
<div class="tablewrap"><table class="choices"><thead><tr><th>Offset d/C</th><th>Initial shortest arc</th><th>Initial target mean</th><th>Last sample: two threshold regions</th><th>Last sample: two resolved maxima</th></tr></thead><tbody>@@ROWS@@</tbody></table></div>
<p>The near activity regions connect by t = 0.2, while two maxima remain resolved through t = 0.8. Farther cases retain two separated regions until t = 6.6 or 6.9; the old maximum decays close to its starting position. A resolved maximum requires circular prominence ≥ 0.05. These are sampled classifications, with 0.1-unit sampling in this window.</p>
<p>All target-region maxima exceed activity one at the first saved time, t = 0.1. This is prompt activation at the available resolution. It is neither a takeover criterion nor evidence for identical activation times. All three long cues ultimately leave the same basal profile translated to the target: full chemical-field residuals are below 9.5 × 10⁻¹¹ after withdrawal.</p>

<h2>What has been verified?</h2>
<div class="statusbox"><h3>Retained 17 October evidence</h3><p>The 10-unit preparation cue is followed by 13,000 cue-free units. Both chemical fields and circumference pass successive settling windows and a 10,000-unit holdout. The retained three-position assay uses a 1,000-unit cue plus 13,000 cue-free units; phase continuity, complete fields, amount conservation and final stationarity pass. The reflected control agrees in the full fields within 1.9 × 10⁻¹⁰.</p><p class="muted">The original preparation recorded a SciPy internal-difference warning. Its saved fields and final derivative are finite and its saved-state audit passes. Fresh replay, resolution and solver checks below address separate numerical questions.</p></div>
<h3 style="margin-top:22px">19 October validation snapshot</h3><div class="statusbox" id="validation">@@VALIDATION@@</div>
<details><summary>Run the reusable ring adapter</summary><p>The development adapter exposes the unchanged ring equations, their complete analytic Jacobian, periodic physical-width cues, exact phase restarts and full activity-moment histories. @@API_LINKS@@</p><p class="caption">The separate API smoke verifies the latest helper, including a zero-active pool and exact staged restart with warnings promoted to errors. It is distinct from the nine-worker scientific validation and does not replace its frozen source record.</p><p>From the repository root, submit the @@EXAMPLE_LINK@@ as a one-CPU SLURM job:</p><pre>@@COMMAND@@</pre><p class="caption">This example restarts the retained basal state, applies the 0.1-unit challenger and recovers without cue. It writes full native states to /tmp/spring_ring_example.npz and refuses to overwrite an existing output. A figure campaign should additionally freeze its sources and record inputs.</p><details><summary>Example Python source</summary><pre>@@EXAMPLE_SOURCE@@</pre></details></details>

<h2>Scientific scope</h2>
<p>This establishes a workable development ring under the stated force-averaging choice, and a location-dependent <em>route through competition</em> at one resting size. The 1,000-unit maintained cue does not measure a minimum pulse duration. The original pilot contains one mesh; new resolution checks must be read at their reported scope. A qualified stationary checkpoint does not prove mathematical stability against every perturbation.</p>
<p>WGD-associated enlargement and increased migration persistence motivate the study. These ring readouts concern local Rac organization and site competition. Observed migration summaries, conditional OU velocity-memory estimates and polarity-site lifetimes remain different measurements. A size-dependent or migration-persistence conclusion requires additional evidence.</p>
<details><summary>Equations, assumptions and exact source evidence</summary><div class="sources">@@SOURCES@@</div><p class="caption">Compact viewer profiles are rounded to 9 decimal places for display; all source checks use retained native arrays. Provenance records bind this reader to the exact input archives and generation snapshot.</p></details>
<footer><p class="sources">@@FOOTER@@</p><p class="muted">Development evidence for review. This reader does not promote Spring into the accepted three-model main ring comparison.</p></footer>
</main><script id="spring-data" type="application/json">@@DATA@@</script><script>
'use strict';
const D=JSON.parse(document.getElementById('spring-data').textContent), E=id=>document.getElementById(id);
const ns='http://www.w3.org/2000/svg';
function tag(name,attrs={},text=''){const e=document.createElementNS(ns,name);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text)e.textContent=text;return e;}
function color(u){const f=Math.max(0,Math.min(1,(u-.65)/1.2));return `rgb(${Math.round(228-215*f)},${Math.round(240-116*f)},${Math.round(244-89*f)})`;}
function draw(){
 const c=D.cases[+E('case').value];E('time').max=c.times.length-1;let index=Math.min(+E('time').value,c.times.length-1);E('time').value=index;
 const t=c.times[index],u=c.profiles[index],C=c.circumference[index],W=700,H=330,L=57,R=18,T=19,B=52,ymin=.65,ymax=1.85;
 const x=q=>L+q*(W-L-R),y=v=>T+(ymax-v)/(ymax-ymin)*(H-T-B),plot=E('profile');plot.replaceChildren();
 for(let v=.7;v<=1.801;v+=.2){plot.append(tag('line',{x1:L,x2:W-R,y1:y(v),y2:y(v),class:'grid'}));plot.append(tag('text',{x:L-9,y:y(v)+4,'text-anchor':'end',class:'tick'},v.toFixed(1)));}
 for(let q=0;q<=1.001;q+=.2){plot.append(tag('line',{x1:x(q),x2:x(q),y1:H-B,y2:H-B+5,class:'axis'}));plot.append(tag('text',{x:x(q),y:H-B+22,'text-anchor':'middle',class:'tick'},q.toFixed(1)));}
 const width=D.cue_width/C,shade=tag('rect',{x:x(c.target-width/2),y:T,width:x(width)-L,height:H-B-T,fill:'#fff0df'});plot.append(shade);
 plot.append(tag('line',{x1:L,x2:W-R,y1:y(D.threshold),y2:y(D.threshold),stroke:'#a16b47','stroke-dasharray':'5 4'}));
 plot.append(tag('line',{x1:x(D.incumbent),x2:x(D.incumbent),y1:T,y2:H-B,stroke:'#aab8c0','stroke-dasharray':'2 4'}));
 function path(values){const edge=(values[0]+values[values.length-1])/2;return `M${x(0)},${y(edge)} `+values.map((v,i)=>`L${x((i+.5)/D.cells)},${y(v)}`).join(' ')+` L${x(1)},${y(edge)}`;}
 plot.append(tag('path',{d:path(D.initial),fill:'none',stroke:'#aab8c0','stroke-width':2}));plot.append(tag('path',{d:path(u),fill:'none',stroke:'#137c9b','stroke-width':2.7}));
 plot.append(tag('line',{x1:L,x2:L,y1:T,y2:H-B,class:'axis'}));plot.append(tag('line',{x1:L,x2:W-R,y1:H-B,y2:H-B,class:'axis'}));
 plot.append(tag('text',{x:(W+L-R)/2,y:H-8,'text-anchor':'middle',class:'tick'},'Material position q — the edges connect'));
 plot.append(tag('text',{x:16,y:(H-B+T)/2,'text-anchor':'middle',transform:`rotate(-90 16 ${(H-B+T)/2})`,class:'tick'},'Raw active Rac'));
 plot.append(tag('text',{x:x(c.target),y:14,'text-anchor':'middle',fill:'#a34f15','font-size':12},'target'));
 const ring=E('ring');ring.replaceChildren();const cx=117.5,cy=111,r=77;
 for(let i=0;i<D.cells;i++){const a=2*Math.PI*i/D.cells-Math.PI/2,b=2*Math.PI*(i+1)/D.cells-Math.PI/2;ring.append(tag('path',{d:`M${cx+r*Math.cos(a)},${cy+r*Math.sin(a)} A${r} ${r} 0 0 1 ${cx+r*Math.cos(b)},${cy+r*Math.sin(b)}`,fill:'none',stroke:color(u[i]),'stroke-width':17}));}
 for(const [q,col,label] of [[D.incumbent,'#738b98','initial'],[c.target,'#c46220','target']]){let a=2*Math.PI*q-Math.PI/2;ring.append(tag('circle',{cx:cx+96*Math.cos(a),cy:cy+96*Math.sin(a),r:4,fill:col}));}
 const re=u.reduce((s,v,i)=>s+v*Math.cos(2*Math.PI*(i+.5)/D.cells),0),im=u.reduce((s,v,i)=>s+v*Math.sin(2*Math.PI*(i+.5)/D.cells),0),total=u.reduce((a,b)=>a+b,0),m=Math.hypot(re,im)/total,angle=Math.atan2(im,re),dx=52*Math.cos(angle-Math.PI/2),dy=52*Math.sin(angle-Math.PI/2);
 ring.append(tag('line',{x1:cx,y1:cy,x2:cx+dx,y2:cy+dy,stroke:'#19354b','stroke-width':2}));ring.append(tag('circle',{cx:cx+dx,cy:cy+dy,r:3,fill:'#19354b'}));ring.append(tag('text',{x:cx,y:219,'text-anchor':'middle','font-size':12,fill:'#526774'},'Arrow: activity direction'));
 E('clock').textContent=`t = ${Number(t.toFixed(2))} native units`;E('size').textContent=`C = ${C.toFixed(6)}`;E('phase').textContent=t===0?'Exact shared initial state':t<1000?'Cue on':t===1000?'Cue withdrawal':'Cue-free continuation';E('moment').textContent=`Activity moment magnitude = ${m.toFixed(4)}`;
}
E('case').addEventListener('change',draw);E('time').addEventListener('input',draw);
E('previous').onclick=()=>{E('time').value=Math.max(0,+E('time').value-1);draw();};E('next').onclick=()=>{E('time').value=Math.min(+E('time').max,+E('time').value+1);draw();};
function at(t){const times=D.cases[+E('case').value].times;E('time').value=times.findIndex(v=>Math.abs(v-t)<1e-8);draw();}
E('at02').onclick=()=>at(.2);E('at5').onclick=()=>at(5);E('final').onclick=()=>at(14000);draw();
</script></body></html>'''


def main():
    generated = datetime.now(timezone.utc).isoformat()
    data, sources = profile_data()
    validation_path = RUN / "results/summary.json"
    validation = read(validation_path) if validation_path.exists() else None
    audit_path = RUN / "results/audit.json"
    audit = read(audit_path) if audit_path.exists() else None
    if validation is not None:
        sources.append(validation_path)
    if audit is not None:
        sources.append(audit_path)
    scheduler = scheduler_snapshot()
    for path in (RUN / "manifest.json", RUN / "launch.json", ROOT / "code/core/models/spring.py",
                 ROOT / "code/mechanisms/spring_ring_model_20261019.py"):
        if path.exists():
            sources.append(path)
    rows = ''.join('<tr>' + ''.join('<td>' + escape(value) + '</td>' for value in (
        f'{c["offset"]:.1f}', f'{c["distance"]:.4f}', f'{c["target_initial_mean"]:.4f}',
        f't = {c["two_components_last"]:.1f}', f't = {c["two_peaks_last"]:.1f}')) + '</tr>'
        for c in data["cases"])
    source_links = [link(ROOT / "SPRING_RING_CLOSURE_20261017.md", "Declared closure and equations"),
                    link(ROOT / "SPRING_RING_DISTANCE_PILOT_20261017.md", "Retained distance result"),
                    link(REPO / "current manuscript figs/reviews/spring_distance_r1.md", "Independent scientific review"),
                    link(PREP / "results/audit.json", "Preparation audit"),
                    link(DIST / "results/audit.json", "Distance audit"),
                    link(DIST / "manifest.json", "Frozen stimulus and solver definition"),
                    link(DIST / "workers", "Retained full native fields"),
                    link(OUT / "supporting/data/spring_ring_reader_data.json", "Compact viewer data"),
                    link(OUT / "supporting/metadata/spring_ring_reader.manifest.json", "Reader provenance")]
    if validation_path.exists():
        source_links.append(link(validation_path, "New validation summary"))
    if audit_path.exists():
        source_links.append(link(audit_path, "New validation audit"))
    for relative, label in [('results/independent_reviewer_audit.json', 'Independent reviewer audit'),
                            ('results/geometry_dilution_20261019/summary.json', 'Signed geometry/dilution measurements')]:
        path = RUN / relative
        if path.exists():
            source_links.append(link(path, label))
            sources.append(path)
    if RUN.exists():
        source_links.append(link(RUN, "New validation source/run record"))
    geometry_dir = RUN / 'results/geometry_dilution_20261019'
    geometry_links = ' · '.join(link(geometry_dir / name, label) for name, label in (
        ('summary.json', 'Signed geometry measurements'), ('timecourses.csv', 'Native signed-rate time courses'))
        if (geometry_dir / name).exists())
    api_smoke = ROOT / 'data/runs/spring_ring_api_smoke_20261019/results/api_smoke.json'
    api_links = link(ROOT / 'code/mechanisms/spring_ring_model_20261019.py', 'Reusable model source')
    if api_smoke.exists():
        api_links += ' · ' + link(api_smoke, 'Latest API smoke record')
        sources.append(api_smoke)
    for path in (geometry_dir / 'timecourses.csv',):
        if path.exists():
            sources.append(path)
    example_path = OUT / 'supporting/data/spring_ring_example.py'
    example_path.parent.mkdir(parents=True, exist_ok=True)
    example_path.write_text(EXAMPLE)
    command = ("sbatch --qos=small --ntasks=1 --cpus-per-task=1 --mem=3G --time=00:10:00 \\\n"
               "  --wrap='env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "
               "/home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(example_path.relative_to(REPO)) + "'")
    footer = link(ROOT / "SPRING_RING_STATUS_20261019.md", "Concise status note") + link(
        ROOT / "figures/length_mechanisms_20261008/dev/README.md", "All development evidence")
    encoded = json.dumps(data, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
    closure_svg = equation_svg(r'$\frac{dC}{dt}=2\left[\left\langle F(u)\right\rangle_q-k(C-C_0)\right],\quad k=0.01$',
                               'dC/dt equals twice the mean force minus spring resistance; k equals 0.01',
                               'spring_ring_closure_equation.svg')
    moment_svg = equation_svg(r'$P=\frac{\sum_j u_j\exp(2\pi i q_j)}{\sum_j u_j}$',
                              'P equals the first circular activity moment normalized by total active Rac',
                              'spring_ring_activity_moment_equation.svg')
    html = HTML
    for name, value in {"GENERATED": escape(generated), "ROWS": rows,
                        "VALIDATION": validation_content(validation, scheduler, audit),
                        "SOURCES": ''.join(source_links), "FOOTER": footer,
                        "DILUTION_SOURCE": link(ROOT / 'code/core/models/spring.py', 'Native interval rhs and moving-domain adapter'),
                        "GEOMETRY_LINKS": geometry_links,
                        "API_LINKS": api_links,
                        "EXAMPLE_LINK": link(example_path, 'downloadable restart example'),
                        "EXAMPLE_SOURCE": escape(EXAMPLE), "COMMAND": escape(command),
                        "DATA": encoded, "CLOSURE_EQUATION": closure_svg,
                        "MOMENT_EQUATION": moment_svg}.items():
        html = html.replace('@@' + name + '@@', value)
    assert '@@' not in html
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'index.html').write_text(html)
    write(OUT / 'supporting/data/spring_ring_reader_data.json', data)
    (OUT / 'RESULTS.md').write_text(
        '# Spring ring development reader\n\n'
        '[Open the scientific reader](index.html) for the declared closure, the raw interactive '
        'three-position competition profiles and the latest collected validation snapshot.\n\n'
        f'Generated {generated}. This page is a snapshot; it has no timed refresh. '
        'The retained trajectories are unchanged. On 19 October the user explicitly agreed to '
        'carry over the existing stretching, dilution and diffusion terms. Historical agent '
        'review of the separately declared mean-force adaptation is recorded separately.\n')
    source_hashes = {str(path.relative_to(REPO)): sha(path) for path in sources}
    source_hashes[str(Path(__file__).relative_to(REPO))] = sha(__file__)
    artifacts = ['index.html', 'supporting/data/spring_ring_reader_data.json', 'RESULTS.md',
                 'supporting/vectors/spring_ring_closure_equation.svg',
                 'supporting/vectors/spring_ring_activity_moment_equation.svg',
                 'supporting/data/spring_ring_example.py']
    write(OUT / 'supporting/metadata/spring_ring_reader.manifest.json', {
        'generated_at': generated, 'artifact_root': '../..',
        'producer': 'finished', 'status': 'Reader generated from saved native fields',
        'trajectory_integrations': 0, 'scheduler_snapshot': scheduler,
        'validation_summary_available': validation is not None,
        'sha256': {p: sha(OUT / p) for p in artifacts}, 'source_sha256': source_hashes,
        'display_precision': 'Nine decimal places; scientific checks use full saved arrays',
        'user_conversion_direction': '19 October explicit user direction: Yes—carry over the existing stretching, dilution and diffusion terms. This approves the native geometric carryover; mean-force closure remains the declared circular adaptation.'})
    print(json.dumps({'reader': str(OUT / 'index.html'), 'generated_at': generated,
                      'source_files': len(source_hashes), 'viewer_samples_per_case': len(data['cases'][0]['times']),
                      'scheduler': scheduler}, indent=2))


if __name__ == '__main__':
    main()
