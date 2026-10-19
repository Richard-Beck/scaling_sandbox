"""Atomic actual-result report for the expanded linear cue experiments.

Read-only with respect to native campaigns. Completed groups are checked once
and cached by immutable result/input/report identities. No integration occurs.
"""
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime, timezone
import argparse
import collections
import csv
import hashlib
import html
import importlib
import json
import os
import re
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'data/runs/revision_five_model_compensation_strength_v2_20261019'
OUTPUT = ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/five_model_compensation_20261019/strength'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()


def read(path): return json.loads(path.read_text())


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def csvwrite(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def scheduler(run, workers):
    launch = read(run/'launch.json') if (run/'launch.json').exists() else {}
    job = str(launch.get('array_job_id', launch.get('job_id', '')))
    if not job: return dict(job_id=None, terminal=False, tasks={}, raw='', query_error=None)
    command = ['sacct', '-j', job, '--noheader', '--parsable2', '--format=JobID%32,State,ExitCode,ElapsedRaw']
    try:
        raw = subprocess.check_output(command, text=True, stderr=subprocess.STDOUT)
        tasks = {}
        for line in raw.splitlines():
            fields = line.split('|'); match = re.fullmatch(re.escape(job)+r'_(\d+)', fields[0].strip())
            if match: tasks[int(match[1])] = dict(state=fields[1], exit_code=fields[2], elapsed_seconds=fields[3])
        live = {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED', 'RESIZING'}
        terminal = len(tasks) == workers and all(t['state'].split()[0] not in live for t in tasks.values())
        return dict(job_id=job, terminal=terminal, tasks=tasks, raw=raw, query_error=None, command=command)
    except subprocess.CalledProcessError as exc:
        return dict(job_id=job, terminal=False, tasks={}, raw=exc.output, query_error=str(exc), command=command)


def arrays(path):
    with np.load(path) as z: return {k: z[k].copy() for k in z.files}


def derive(archive, case, z, engine):
    a = arrays(archive); fields = {k[6:]: v for k, v in a.items() if k.startswith('field_')}
    active = fields[engine.p.active_key(case['model'])]
    tr = SimpleNamespace(time=a['time'], fields=fields, active=active,
                         orientation=a['orientation'], contrast=a['contrast'])
    return engine.regional(tr, case, z, 'derived')[0]


def checked(path, expected, sources):
    actual = sha(path); assert actual == expected, str(path)
    sources[str(path.relative_to(ROOT))] = actual


def compact_events(a, case, z, engine, phase, duration):
    out = []
    time = a['time']; fraction = a['original_core_active_fraction']
    for loss in [.1, .25, .5]:
        out.append(dict(phase=phase, pulse_duration=duration, observable='fixed_original_core_active_fraction',
                        threshold=1-loss, **engine.event(time, fraction <= 1-loss)))
    unavailable = bool(np.any(z['initial_core_mask'] & z['initial_cue_footprint']))
    desc = dict(phase=phase, pulse_duration=duration,
                observable='separated_challenger_peak_above_initial_halfcontrast',
                threshold=float(z['initial_halfcontrast_threshold']))
    if unavailable: out.append(dict(**desc, status='unavailable_initial_cue_core_overlap', lower=None, upper=None))
    else:
        good = (a['actual_cue_region_active_peak'] >= z['initial_halfcontrast_threshold']) & (a['actual_cue_original_core_overlap_fraction'] == 0)
        out.append(dict(**desc, **engine.event(time, good)))
    for row in out:
        row.update(observation_start=float(time[0]), observation_end=float(time[-1]),
                   actual_min_sample_interval=float(np.min(np.diff(time))), actual_max_sample_interval=float(np.max(np.diff(time))),
                   timing='phase-local native time; sampled-entry bracket, not confidence interval')
    return out


def process(case, run, out, engine, report_sha, digest):
    reuse = case.get('exact_reuse'); folder = run/'workers'/f'{case["slot"]:04d}'
    result_path = ROOT/reuse['source_result'] if reuse else folder/'result.json'
    if not result_path.exists(): return None
    result_sha = sha(result_path)
    if reuse:
        assert result_sha == reuse['source_result_sha256']
        assert sha(ROOT/reuse['numerical_screen']) == reuse['numerical_screen_sha256']
    key = hashlib.sha256((result_sha+case['regional_input_sha256']+report_sha+digest).encode()).hexdigest()
    cache = out/'supporting/data/case_cache'/f'{case["slot"]:04d}_{key}.json'
    if cache.exists(): return read(cache)
    result = read(result_path)
    if not reuse: assert result['case'] == case and result['manifest_sha256'] == digest
    z = arrays(run/case['regional_input']); sources = {str(result_path.relative_to(ROOT)): result_sha}
    checked(run/case['regional_input'], case['regional_input_sha256'], sources)
    if reuse:
        sources[reuse['numerical_screen']] = reuse['numerical_screen_sha256']
        sources[reuse['source_manifest']] = reuse['source_manifest_sha256']
        assert sha(ROOT/reuse['source_manifest']) == reuse['source_manifest_sha256']
        withdrawals = reuse['withdrawals']
        on_spec = None if case['case_kind']=='zero_control' else dict(archive=reuse['cue_archive'], archive_sha256=reuse['cue_archive_sha256'])
        base = ROOT
    else:
        withdrawals = result['withdrawals']; on_spec = result['on']; base = run
    on = None; on_path = None
    if on_spec:
        on_path = base/on_spec['archive']; checked(on_path, on_spec['archive_sha256'], sources)
        if 'regional' in on_spec:
            regional_path = base/on_spec['regional']; checked(regional_path, on_spec['regional_sha256'], sources)
            on = arrays(regional_path)
        else:
            on = derive(on_path, case, z, engine)
            regional_path = out/'supporting/data/case_cache'/f'{case["slot"]:04d}_reused_on.npz'
            regional_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(regional_path, **on)
    rows = []; events = []; traces = []
    for w in withdrawals:
        duration = float(w['duration']); archive = base/w['archive']
        checked(archive, w['archive_sha256'], sources); checked(base/w['restart'], w['restart_sha256'], sources)
        if 'regional' in w:
            regional = base/w['regional']; checked(regional, w['regional_sha256'], sources); off = arrays(regional)
        else:
            off = derive(archive, case, z, engine)
            regional = out/'supporting/data/case_cache'/f'{case["slot"]:04d}_reused_off_{duration:g}.npz'
            regional.parent.mkdir(parents=True, exist_ok=True); np.savez_compressed(regional, **off)
        row = dict(pulse_duration=duration, source_kind='exact_validated_high_strength_native_reuse' if reuse else 'new_native_full_fields',
                   actual_off_archive=str(archive.relative_to(ROOT)), actual_off_archive_sha256=w['archive_sha256'],
                   actual_restart=str((base/w['restart']).relative_to(ROOT)), actual_restart_sha256=w['restart_sha256'],
                   off_regional_readouts=str(regional.relative_to(ROOT)),
                   on_archive=str(on_path.relative_to(ROOT)) if on_path else '',
                   off_observation_end=float(off['time'][-1]),
                   original_scalar_orientation_success=w['orientation_success'],
                   original_scalar_full_polarity_success=w['full_polarity_success'],
                   native_field_tail_pass=w['field_tail_pass'], native_field_tail_max_scaled=w['field_tail_max_scaled'],
                   final_orientation=w['final_orientation'],
                   maximum_relative_polarity_error_last_500=w['maximum_relative_polarity_error_last_500'],
                   no_cue_preparation_valid=w['no_cue_preparation_valid'],
                   saved_initial_max_absolute_roundoff=w['saved_initial_max_absolute_roundoff'])
        for name, values in off.items():
            if name != 'time': row[name+'_at_withdrawal'] = float(values[0]); row[name+'_final'] = float(values[-1])
        if on is not None:
            mask = on['time'] <= duration
            assert on['time'][mask][-1] == duration
            window = {k: a[mask] for k, a in on.items()}
            events.extend(compact_events(window, case, z, engine, 'cue_on_truncated_at_actual_withdrawal', duration))
            row['minimum_original_core_fraction_during_cue'] = float(window['original_core_active_fraction'].min())
            row['maximum_actual_cue_region_peak_during_cue'] = float(window['actual_cue_region_active_peak'].max())
        events.extend(compact_events(off, case, z, engine, 'after_withdrawal', duration))
        rows.append(row)
        traces.append(dict(pulse_duration=duration, phase='after_withdrawal', path=str(regional.relative_to(ROOT)),
                           actual_native_archive=str(archive.relative_to(ROOT)), sample_count=len(off['time'])))
    if on is not None:
        traces.append(dict(pulse_duration='shared_constant_ON_prefix', phase='cue_on', path=str(regional_path.relative_to(ROOT)),
                           actual_native_archive=str(on_path.relative_to(ROOT)), sample_count=len(on['time'])))
    data = dict(rows=rows, events=events, traces=traces, source_hashes=sources,
                result_path=str(result_path.relative_to(ROOT)), result_sha256=result_sha,
                native_worker_wall_seconds=result.get('wall_seconds'), peak_RSS_KiB=result.get('peak_resident_memory_KiB'))
    write(cache, data); return data


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run', type=Path, default=DEFAULT); parser.add_argument('--out', type=Path, default=OUTPUT)
    args = parser.parse_args(); run=args.run.resolve(); out=args.out.resolve()
    manifest=read(run/'manifest.json'); provenance=read(run/'provenance.json'); digest=sha(run/'manifest.json')
    assert digest == provenance['manifest_sha256']
    sources = {str((run/'manifest.json').relative_to(ROOT)):digest, str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}
    for rel, expected in provenance['frozen_sources'].items(): checked(run/'frozen_sources'/rel, expected, sources)
    for rel, expected in provenance['frozen_inputs'].items(): checked(run/rel, expected, sources)
    sys.path.insert(0, str(run/'frozen_sources/code/mechanisms'))
    engine=importlib.import_module('revision_five_model_compensation_campaign_20261019')
    schedule=scheduler(run, len(manifest['jobs'])); task_by_slot={j['case_slot']:j['task'] for j in manifest['jobs']}
    groups=[]; outcomes=[]; events=[]; traces=[]; errors=[]
    for c in manifest['cases']:
        common=dict(c['initial_geometry'])
        common.update(case_slot=c['slot'], model=c['model'], length=c['length'], width=c['width'],
                    native_peak_multiplier=c['native_peak_multiplier'], requested_native_peak=c['response_amplitude'],
                    case_kind=c['case_kind'], checkpoint_sha256=c['checkpoint_sha256'],
                    assay_coordinate_system=c['assay_coordinate_system'])
        # Geometry carries model/length/checkpoint too; all values must agree.
        folder=run/'workers'/f'{c["slot"]:04d}'; task=task_by_slot.get(c['slot'])
        try: data=process(c,run,out,engine,sources[str(Path(__file__).relative_to(ROOT))],digest)
        except Exception as exc:
            data=None; errors.append(dict(case_slot=c['slot'], error=repr(exc))); state='collector_validation_error'
        else:
            state='complete' if data else 'error' if (folder/'error.json').exists() else read(folder/'status.json')['state'] if (folder/'status.json').exists() else 'not_started'
        group=dict(**common, state=state, native_task=task,
                   scheduler_state=schedule['tasks'].get(task,{}).get('state'),
                   exact_native_reuse=bool(c.get('exact_reuse')))
        if data:
            sources.update(data['source_hashes']); group.update(result_path=data['result_path'], result_sha256=data['result_sha256'],
                native_worker_wall_seconds=data['native_worker_wall_seconds'], peak_RSS_KiB=data['peak_RSS_KiB'])
            for row in data['rows']: outcomes.append(dict(**common, state=state, **row))
            events.extend(dict(case_slot=c['slot'],model=c['model'],length=c['length'],width=c['width'],native_peak_multiplier=c['native_peak_multiplier'],**row) for row in data['events'])
            traces.extend(dict(case_slot=c['slot'],model=c['model'],length=c['length'],width=c['width'],**row) for row in data['traces'])
        else:
            for duration in c['pulse_durations'] if c['case_kind']=='constant_on' else [0.]: outcomes.append(dict(**common,state=state,pulse_duration=duration))
        groups.append(group)
    assert len(outcomes)==manifest['logical_protocol_rows']
    done=sum(g['state']=='complete' for g in groups); completed=sum(r['state']=='complete' for r in outcomes)
    failed=sum(g['state'] in ['error','collector_validation_error'] for g in groups)
    producer='finished' if done==len(groups) else 'terminal_array_incomplete' if schedule['terminal'] else 'execution_errors' if failed else 'running' if schedule['job_id'] else 'frozen_not_submitted'
    now=datetime.now(timezone.utc).isoformat(); stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    generation=out/'supporting/generations'/stamp
    (generation/'data').mkdir(parents=True); (generation/'metadata').mkdir()
    for name, rows in [('case_groups.csv',groups),('pulse_outcomes.csv',outcomes),('sampled_events.csv',events),('regional_trace_sources.csv',traces)]: csvwrite(generation/'data'/name,rows)
    status=dict(generated_at=now,producer=producer,complete=done==len(groups),completed_case_groups=done,total_case_groups=len(groups),
                completed_protocol_rows=completed,total_protocol_rows=len(outcomes),errors=failed,
                array_job_id=schedule['job_id'],terminal_scheduler_observed=schedule['terminal'],
                scheduler_counts=dict(collections.Counter(t['state'] for t in schedule['tasks'].values())),
                exact_high_strength_reuse_groups=manifest['exact_reuse_case_groups'],manifest_sha256=digest,
                numerical_review='Native pilot validated; new full-stage independent terminal audit pending.')
    write(generation/'metadata/latest_status.json',status); write(generation/'metadata/scheduler.json',schedule)
    write(generation/'metadata/collection_errors.json',errors)
    validation=dict(generated_at=now,source_hashes=sources,
                    artifact_hashes={str(p.relative_to(generation)):sha(p) for p in (generation/'data').glob('*.csv')},
                    caching='Actual result and archive hashes checked on first collection; immutable cached records keyed by result/input/report/manifest SHA.',
                    no_native_integration=True)
    write(generation/'metadata/validation.json',validation)
    body=[]
    for model in manifest['models']:
        rr=[r for r in outcomes if r['model']==model]; ready=[r for r in rr if r['state']=='complete']
        positive=[r for r in ready if r['pulse_duration']>0]
        body.append(f'<tr><td>{html.escape(model)}</td><td>{len(ready)}/{len(rr)}</td><td>{sum(r["original_scalar_full_polarity_success"] for r in positive)}/{len(positive)}</td><td>{sum(r["native_field_tail_pass"] for r in ready)}/{len(ready)}</td></tr>')
    display=[]
    for r in outcomes:
        if r['state']!='complete' or r['pulse_duration']==0:continue
        f=lambda k:format(r[k],'.4g')
        display.append('<tr data-model="'+r['model']+'"><td>'+html.escape(r['model'])+'</td><td>'+f('length')+'</td><td>'+f('width')+'</td><td>'+f('pulse_duration')+'</td><td>'+f('original_core_active_fraction_at_withdrawal')+'</td><td>'+f('actual_cue_region_active_peak_at_withdrawal')+'</td><td>'+str(r['original_scalar_full_polarity_success'])+'</td><td>'+str(r['native_field_tail_pass'])+'</td><td>'+('overlap' if not r['separated_challenger_readout_available'] else 'separate')+'</td></tr>')
    page='''<!doctype html><meta charset="utf-8"><title>Five-model cue compensation: actual results</title><style>body{font:16px system-ui;max-width:1250px;margin:28px auto;padding:0 18px;line-height:1.45}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #ddd;padding:6px;text-align:left}thead{position:sticky;top:0;background:#fff}p{max-width:1100px}.state{background:#eef4fa;padding:12px}a{color:#145c91}</style>'''
    page+=f'<h1>Five-model cue compensation: actual results</h1><p class="state">Generated {now}. Producer: <b>{producer}</b>. Completed {completed}/{len(outcomes)} requested protocols ({done}/{len(groups)} case groups). Array {schedule["job_id"]}. No timed refresh; reopen or manually refresh after a new collection.</p>'
    page+='<p>All 97 exact qualified basals are tested at the middle native width with four additional peak factors (.03, .1, .3 and 3) and five pulse withdrawals. The separately retained baseline adds 1× and two other widths plus no-cue controls. The tables retain pending/error rows, actual source identities and initial overlap. Requested input is model-specific; equal peak factors are not matched optical doses. Broader cues also change integrated input.</p>'
    page+='<p><a href="data/pulse_outcomes.csv">All requested pulse outcomes</a> · <a href="data/sampled_events.csv">Sampled regional entry/retained-entry brackets</a> · <a href="data/case_groups.csv">Case groups and geometry</a> · <a href="data/regional_trace_sources.csv">Exact regional/native trace sources</a> · <a href="metadata/latest_status.json">Status</a> · <a href="metadata/validation.json">Source hashes</a></p>'
    page+='<table><thead><tr><th>Model</th><th>Completed protocols</th><th>Positive pulses: scalar full-strength pass</th><th>Field-tail pass (all completed)</th></tr></thead><tbody>'+''.join(body)+'</tbody></table>'
    page+='<p>Original-core activity is the amount in the frozen initial halfcontrast component; it is regional depletion, not identified molecular damage. Cue-associated peak is a native concentration and may include incumbent activity when the initial cue footprint overlaps. Spring regional amounts use evolving physical weights; deBelly Rac integrals use nominal coordinates, with physical length separately recorded.</p>'
    page+='<p>Cue-on event searches are truncated at each actual withdrawal. After-withdrawal times start at zero. Brackets reflect saved native sampling, not confidence intervals. A failed coarse duration scan means <b>no tested success among sampled durations up to the maximum</b>; it neither proves monotonicity nor excludes success islands between two failed tests. Original scalar gates and native field-tail stability remain separate from regional measures.</p>'
    page+='<h2>Completed positive pulse outcomes</h2><label>Model <select id="model"><option value="all">All</option>'+''.join(f'<option>{html.escape(m)}</option>' for m in manifest['models'])+'</select></label><table><thead><tr><th>Model</th><th>Native length</th><th>Width</th><th>Pulse</th><th>Old-core fraction at withdrawal</th><th>Cue-ROI native peak at withdrawal</th><th>Scalar full strength</th><th>Field tail</th><th>Initial cue/core</th></tr></thead><tbody id="rows">'+''.join(display)+'</tbody></table><script>document.getElementById("model").addEventListener("change",e=>{for(const r of document.querySelectorAll("#rows tr"))r.hidden=e.target.value!=="all"&&r.dataset.model!==e.target.value;});</script>'
    (generation/'index.html').write_text(page)
    target=os.path.relpath(generation,out);temp=out/'current.tmp'
    if temp.is_symlink():temp.unlink()
    temp.symlink_to(target);temp.replace(out/'current')
    (out/'README.md').write_text('# Expanded five-model cue experiment: live actual results\n\n[Open the latest actual-results report](current/index.html). The report shows its generation time and producer state; it does not refresh automatically.\n\n'+f'Latest snapshot: {completed}/{len(outcomes)} protocols complete; **{producer}**, generated {now}.\n\n'+
        'This is a coarse native length/width/pulse scan, not a final figure or a monotone minimum-duration estimate. [All requested outcomes](current/data/pulse_outcomes.csv), [sampled regional event brackets](current/data/sampled_events.csv), [exact trace sources](current/data/regional_trace_sources.csv) and [source hashes](current/metadata/validation.json) remain available. No new native integration is performed by this report.\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
