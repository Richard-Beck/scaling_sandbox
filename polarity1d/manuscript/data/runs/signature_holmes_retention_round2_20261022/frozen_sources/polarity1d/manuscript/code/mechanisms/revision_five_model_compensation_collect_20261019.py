"""Atomic pilot collection and measured cost forecast; no native integration."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/revision_five_model_compensation_pilot_v2_20261019'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/five_model_compensation_20261019/pilot'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):return json.loads(p.read_text())


def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp.json')
    tmp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');tmp.replace(p)


def csvwrite(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True);keys=list(dict.fromkeys(k for r in rows for k in r));tmp=p.with_suffix('.tmp.csv')
    with tmp.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    tmp.replace(p)


def allocated(path):return sum(p.stat().st_blocks*512 for p in path.rglob('*') if p.is_file())


def main():
    for name in ['supporting/data','supporting/metadata']:(OUT/name).mkdir(parents=True,exist_ok=True)
    manifest=read(RUN/'manifest.json');provenance=read(RUN/'provenance.json');digest=sha(RUN/'manifest.json')
    assert digest==provenance['manifest_sha256']
    for rel,h in provenance['frozen_sources'].items():assert sha(RUN/'frozen_sources'/rel)==h,rel
    for rel,h in provenance['frozen_inputs'].items():assert sha(RUN/rel)==h,rel
    launch=read(RUN/'launch.json');job=str(launch.get('array_job_id',launch.get('job_id','21202256')))
    command=['sacct','-j',job,'--noheader','--parsable2','--format=JobID%30,State,ExitCode,ElapsedRaw']
    raw=subprocess.check_output(command,text=True)
    scheduler={r.split('|')[0].strip():r.split('|')[1:4] for r in raw.strip().splitlines() if r}
    tasks={k:v for k,v in scheduler.items() if re.fullmatch(re.escape(job)+r'_\d+',k)}
    terminal=bool(len(tasks)==10 and all(v[0] not in ['PENDING','RUNNING','CONFIGURING','COMPLETING'] for v in tasks.values()))
    summaries=[];branches=[];prefixes=[];sources={str((RUN/'manifest.json').relative_to(ROOT)):digest};phase_costs=[]
    for c in manifest['cases']:
        folder=RUN/'workers'/f'{c["slot"]:03d}';result=folder/'result.json'
        row=dict(slot=c['slot'],model=c['model'],length=c['length'],state='not_started',
                 checkpoint_sha256=c['checkpoint_sha256'])
        if not result.exists():
            if (folder/'error.json').exists():row['state']='error';row['error_path']=str((folder/'error.json').relative_to(ROOT))
            elif (folder/'status.json').exists():row.update(state=read(folder/'status.json')['state'])
            summaries.append(row);continue
        r=read(result);assert r['case']==c and r['manifest_sha256']==digest
        assert sha(folder/'cue.npz')==r['cue_archive_sha256']
        sources[str(result.relative_to(ROOT))]=sha(result)
        sources[str((folder/'cue.npz').relative_to(ROOT))]=r['cue_archive_sha256']
        row.update(state='complete',wall_seconds=r['wall_seconds'],peak_resident_memory_KiB=r['peak_resident_memory_KiB'],
            allocated_worker_bytes=allocated(folder),native_on_wall_seconds=r['native_on_wall_seconds'],
            on_raw_field_bytes=r['on_raw_field_bytes'],on_allocated_archive_bytes=(folder/'cue.npz').stat().st_blocks*512,
            maximum_ON_prefix_scaled_difference=max(p['maximum_scaled_state_difference'] for p in r['prefixes']),
            zero_control_valid=r['withdrawals'][0]['no_cue_preparation_valid'])
        for p in r['prefixes']:prefixes.append(dict(slot=c['slot'],model=c['model'],length=c['length'],**p))
        refs=[]
        for w in r['withdrawals']:
            for path,hkey in [('archive','archive_sha256'),('restart','restart_sha256')]:
                assert sha(RUN/w[path])==w[hkey];sources[str((RUN/w[path]).relative_to(ROOT))]=w[hkey]
            ref=w['direct_native_OFF_reference']
            if ref:
                assert sha(RUN/ref['final_state'])==ref['final_state_sha256']
                sources[str((RUN/ref['final_state']).relative_to(ROOT))]=ref['final_state_sha256'];refs.append(ref)
            b=dict(slot=c['slot'],model=c['model'],length=c['length'],pulse_duration=w['duration'],
                orientation_success=w['orientation_success'],full_polarity_success=w['full_polarity_success'],field_tail_pass=w['field_tail_pass'],
                final_orientation=w['final_orientation'],field_tail_max_scaled=w['field_tail_max_scaled'],
                saved_initial_max_absolute_roundoff=w['saved_initial_max_absolute_roundoff'],
                saved_initial_roundoff_bound=w['saved_initial_roundoff_bound'],
                native_off_wall_seconds=w['native_off_wall_seconds'],off_raw_field_bytes=w['off_raw_field_bytes'],
                allocated_off_archive_bytes=(RUN/w['archive']).stat().st_blocks*512,
                direct_OFF_reference_present=bool(ref),direct_OFF_max_final_scaled_difference=ref['maximum_final_scaled_difference'] if ref else '',
                direct_OFF_all_flags_agree=ref['all_scalar_and_field_tail_flags_agree'] if ref else '',
                actual_archive=str((RUN/w['archive']).relative_to(ROOT)),actual_archive_sha256=w['archive_sha256'])
            branches.append(b)
        row.update(direct_OFF_references=len(refs),maximum_direct_OFF_final_scaled_difference=max(q['maximum_final_scaled_difference'] for q in refs),
                   all_direct_OFF_flags_agree=all(q['all_scalar_and_field_tail_flags_agree'] for q in refs))
        summaries.append(row)
    forecasts=[]
    for model in manifest['models']:
        ss=[r for r in summaries if r['model']==model and r['state']=='complete'];bb=[r for r in branches if r['model']==model]
        if not ss:continue
        lengths=sum(c['model']==model for c in manifest['length_series_inventory'])
        # 3 baseline-width ON histories + four extra middle-width strengths.
        on_groups=7*lengths;off_branches=35*lengths;zero_controls=lengths
        normal=[b for b in bb if b['pulse_duration']>0];zero=[b for b in bb if b['pulse_duration']==0]
        mean=lambda v:sum(v)/len(v)
        cpu=on_groups*mean([s['native_on_wall_seconds'] for s in ss])+off_branches*mean([b['native_off_wall_seconds'] for b in normal])+zero_controls*mean([b['native_off_wall_seconds'] for b in zero])
        compressed=on_groups*mean([s['on_allocated_archive_bytes'] for s in ss])+off_branches*mean([b['allocated_off_archive_bytes'] for b in normal])+zero_controls*mean([b['allocated_off_archive_bytes'] for b in zero])
        rawbytes=on_groups*mean([s['on_raw_field_bytes'] for s in ss])+(off_branches+zero_controls)*mean([b['off_raw_field_bytes'] for b in bb])
        forecasts.append(dict(model=model,qualified_lengths=lengths,pilot_complete_workers=len(ss),
            full_constant_ON_groups=on_groups,full_positive_withdrawals=off_branches,full_zero_controls=zero_controls,
            baseline_pilot_mean_native_CPU_hours=cpu/3600,baseline_pilot_mean_allocated_archive_bytes=compressed,
            baseline_pilot_single_CPU_elapsed_hours_including_output=compressed/mean([s['allocated_worker_bytes'] for s in ss])*mean([s['wall_seconds'] for s in ss])/3600,
            raw_full_field_payload_bytes=rawbytes,max_pilot_peak_RSS_KiB=max(s['peak_resident_memory_KiB'] for s in ss),
            forecast_limit='Baseline short/long pilot extrapolation; widths/strengths may alter stiffness/compression. Excludes refinements and extra high-dose numerical companions.'))
    now=datetime.now(timezone.utc).isoformat();done=sum(r['state']=='complete' for r in summaries);errors=sum(r['state']=='error' for r in summaries)
    producer='finished' if done==10 else 'terminal_array_incomplete' if terminal else 'execution_errors' if errors else 'running' if any(r['state']!='not_started' for r in summaries) else 'queued'
    status=dict(generated_at=now,producer=producer,complete=done==10,completed_workers=done,total_workers=10,errors=errors,
        terminal_scheduler_observed=terminal,array_job_id=job,pilot_protocol_rows=60,completed_branch_rows=len(branches),
        actual_run_allocated_bytes=allocated(RUN),full_strength_levels=[.03,.1,.3,1.,3.],full_planned_rows=3492,
        numerical_fidelity_independent_terminal_review='pending')
    for name,table in [('worker_summary.csv',summaries),('withdrawal_outcomes.csv',branches),('prefix_fidelity.csv',prefixes),('measured_forecast.csv',forecasts)]:csvwrite(OUT/'supporting/data'/name,table)
    write(RUN/'results/latest_status.json',status);write(OUT/'supporting/metadata/latest_status.json',status)
    write(OUT/'supporting/metadata/scheduler.json',dict(generated_at=now,command=command,raw=raw))
    sources[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    artifacts={str(p.relative_to(OUT)):sha(p) for p in (OUT/'supporting/data').glob('*.csv')}
    write(OUT/'supporting/metadata/validation.json',dict(generated_at=now,source_hashes=sources,artifact_hashes=artifacts))
    (OUT/'README.md').write_text('# Five-model native platform pilot\n\n'+f'Generated {now}; producer **{producer}**. Completed workers {done}/10 and withdrawal rows {len(branches)}/60. No native rerun or integration is performed by this collector.\n\n'+
        '[Worker fidelity/resources](supporting/data/worker_summary.csv) · [Actual withdrawal outcomes](supporting/data/withdrawal_outcomes.csv) · [Native ON-prefix comparisons](supporting/data/prefix_fidelity.csv) · [Measured full-stage forecast](supporting/data/measured_forecast.csv) · [Status](supporting/metadata/latest_status.json)\n\n'+
        'The planned full series contains all 97 accepted qualified basals, three baseline widths and five middle-width strength levels, with five coarse pulse withdrawals and one no-cue control per basal. Forecasts extrapolate baseline short/long measurements; high-strength numerical checks and refined transitions remain additional work. No tested success among sampled pulses does not exclude an unobserved success island. Independent terminal fidelity review is pending.\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
