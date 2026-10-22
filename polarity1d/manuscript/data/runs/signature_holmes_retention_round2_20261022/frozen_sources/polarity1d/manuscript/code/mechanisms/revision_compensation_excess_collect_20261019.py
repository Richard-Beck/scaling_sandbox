"""Compact collection of background-aware derived responses; no integration."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import os
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'data/runs/revision_compensation_excess_baseline_20261019'
OUTPUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/five_model_compensation_20261019/excess'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def rows(p):return list(csv.DictReader(p.open()))
def write(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def csvwrite(p,a):
    keys=list(dict.fromkeys(k for r in a for k in r))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(a)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,default=DEFAULT);parser.add_argument('--out',type=Path,default=OUTPUT)
    args=parser.parse_args();run=args.run.resolve();out=args.out.resolve();m=read(run/'manifest.json')
    initials=[];endpoints=[];events=[];summaries=[];sources={str(run/'manifest.json'):sha(run/'manifest.json'),str(Path(__file__)):sha(Path(__file__))}
    for model in m['models']:
        folder=run/'workers'/model;p=folder/'result.json'
        if not p.exists():
            summaries.append(dict(model=model,state=read(folder/'status.json')['state'] if (folder/'status.json').exists() else 'not_started'));continue
        r=read(p);assert r['manifest_sha256']==sha(run/'manifest.json')
        sources[str(p)]=sha(p)
        for name,target in [('initial_core_definitions.csv',initials),('pulse_endpoints.csv',endpoints),('sampled_loss_events.csv',events)]:
            path=folder/name;assert sha(path)==r['derived_artifact_hashes'][str(path.relative_to(run))]
            sources[str(path)]=sha(path);target.extend(rows(path))
        summaries.append(dict(model=model,state='complete',completed_case_groups=r['completed_case_groups'],
            maximum_relative_organization_replay_error=r['maximum_relative_organization_replay_error'],endpoint_rows=r['endpoint_rows'],
            event_rows=r['event_rows'],wall_seconds=r['wall_seconds'],peak_RSS_KiB=r['peak_RSS_KiB']))
    # Floor/ROI sensitivity is descriptive, without new outcome gates.
    sensitivity=[]
    for r in endpoints:
        base={k:r[k] for k in ['case_slot','model','length','width','native_peak_multiplier','phase','pulse_duration']}
        a=float(r['organization_core_current_minimum_excess_fraction_last'])
        b=float(r['organization_core_current_sample_quantile_0p05_excess_fraction_last'])
        c=float(r['organization_core_frozen_initial_minimum_excess_fraction_last'])
        d=float(r['worker_core_current_minimum_excess_fraction_last'])
        sensitivity.append(dict(**base,organization_raw_fraction=float(r['organization_core_raw_fraction_last']),
            organization_current_minimum_excess_fraction=a,organization_current_quantile_excess_fraction=b,
            organization_frozen_minimum_excess_fraction=c,worker_current_minimum_excess_fraction=d,
            current_minimum_minus_quantile=a-b,current_minimum_minus_frozen=a-c,fractional_minus_whole_cell_ROI=a-d,
            scalar_full_polarity_success=r['scalar_full_polarity_success'],native_field_tail_pass=r['native_field_tail_pass']))
    now=datetime.now(timezone.utc).isoformat();stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    generation=out/'supporting/generations'/stamp;(generation/'data').mkdir(parents=True);(generation/'metadata').mkdir()
    tables={'initial_core_definitions.csv':initials,'pulse_endpoints.csv':endpoints,'sampled_loss_events.csv':events,'endpoint_floor_ROI_sensitivity.csv':sensitivity,'worker_summary.csv':summaries}
    for name,rr in tables.items():csvwrite(generation/'data'/name,rr)
    done=sum(r['state']=='complete' for r in summaries)
    status=dict(generated_at=now,producer='finished' if done==5 else 'partial_derived_analysis',complete=done==5,
        completed_model_workers=done,total_model_workers=5,completed_initial_core_rows=len(initials),endpoint_rows=len(endpoints),event_rows=len(events),
        no_native_integrations=True,independent_derived_review='pending')
    write(generation/'metadata/latest_status.json',status)
    write(generation/'metadata/validation.json',dict(generated_at=now,source_hashes=sources,
        artifact_hashes={str(p.relative_to(generation)):sha(p) for p in (generation/'data').glob('*.csv')},
        scientific_scope='Named localized-signal measures and floor/ROI sensitivity; not molecular depletion or a replacement gate.'))
    text='''# Fixed-initial-core organization and challenge response

This analysis keeps native raw regional activity and separately measures localized signal above explicitly named profile floors. It does not replace the scalar or full-field outcome flags.

The retained organization figure defines its basal background as the initial profile minimum. Its initial half-height edge is linearly interpolated; amount uses a fractional finite-volume edge cell, or piecewise-linear trapezoid integration on Spring's physical domain. The experiment worker's frozen whole-cell set is also retained and separately named. deBelly's activity integrals remain nominal-coordinate quantities.

For each fixed initial core, positive excess is integrated above (i) the current full-profile minimum, (ii) the current native-sample fifth percentile, and (iii) the frozen basal minimum. Every rule uses its own initial denominator. The current floor describes a concentration reference; it is not a measured molecular pool. The fifth percentile is an unweighted native-grid concentration quantile, not a physical-volume percentile.

Spring's current-floor excess amount with evolving physical weights is invariant under pure uniform dilution; its frozen basal concentration subtraction is not. Thus organization loss, regional redistribution, global concentration shifts and retained replacement must remain distinct.

Entry and retained-entry times are saved-sample brackets in native phase-local clocks. Cue-on windows stop at the actual pulse withdrawal. A finite failure to cross a threshold is not an infinite latency or an optimized minimum-duration bound.
'''
    text+=f'\nGenerated {now}; **{status["producer"]}**, {done}/5 model workers.\n\n'
    text+='[Native-amplitude endpoints](data/pulse_endpoints.csv) · [Floor/ROI sensitivity](data/endpoint_floor_ROI_sensitivity.csv) · [Sampled loss brackets](data/sampled_loss_events.csv) · [Initial definitions](data/initial_core_definitions.csv) · [Source proof](metadata/validation.json).\n'
    (generation/'README.md').write_text(text)
    temp=out/'current.tmp'
    if temp.is_symlink():temp.unlink()
    temp.symlink_to(os.path.relpath(generation,out));temp.replace(out/'current')
    (out/'README.md').write_text('# Background-aware challenge readouts\n\n[Latest evidence, definitions and sensitivity tables](current/README.md).\n\n'+f'Generated {now}; {done}/5 derived model workers complete. Native workers and gates are unchanged.\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
