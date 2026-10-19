"""Two targeted OFF continuations from exact native DB T500 saved states."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import copy
import json
import resource
import shutil
import time
import traceback
import numpy as np
import revision_five_model_compensation_20261019 as p
import revision_five_model_compensation_campaign_20261019 as campaign

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'
DEFAULT=ROOT/'data/runs/revision_debelly_paired_withdrawal_20261019'


def initialize(out):
    assert not out.exists();m=p.read(BASE/'manifest.json');base_prov=p.read(BASE/'provenance.json')
    for rel,h in p.n.sources().items():assert base_prov['frozen_sources'][rel]==h
    for name in ['inputs','frozen_sources','workers','logs']:(out/name).mkdir(parents=True,exist_ok=True)
    inputs={};cases=[]
    for old in m['cases']:
        if old['model']!='debelly' or old['length']!=39. or old['width'] not in [1.5,5.] or old['case_kind']!='constant_on':continue
        result_path=BASE/'workers'/f'{old["slot"]:04d}'/'result.json';r=p.read(result_path);assert r['case']==old
        archive=BASE/r['on']['archive'];assert p.n.sha(archive)==r['on']['archive_sha256']
        with np.load(archive) as a:
            ii=np.flatnonzero(a['time']==500.);assert len(ii)==1
            restart={name:a['field_'+name][int(ii[0])].copy() for name in p.n.DB_NAMES}
        c=copy.deepcopy(old);slot=len(cases);c.update(slot=slot,source_baseline_case_slot=old['slot'],pulse_duration=500.,off_duration=10000.)
        for field,name in [('frozen_checkpoint','basal'),('regional_input','region')]:
            src=BASE/old[field];dest=out/'inputs'/f'{slot:02d}_{name}.npz';shutil.copy2(src,dest)
            inputs[str(dest.relative_to(out))]=p.n.sha(dest);c[field]=str(dest.relative_to(out))
        assert p.n.sha(out/c['regional_input'])==c['regional_input_sha256']
        rp=out/'inputs'/f'{slot:02d}_restart_T500.npz';p.save_state(rp,'debelly',restart);inputs[str(rp.relative_to(out))]=p.n.sha(rp)
        c.update(frozen_restart=str(rp.relative_to(out)),frozen_restart_sha256=p.n.sha(rp),
            source_ON_archive=str(archive.relative_to(ROOT)),source_ON_archive_sha256=r['on']['archive_sha256'],
            source_result=str(result_path.relative_to(ROOT)),source_result_sha256=p.n.sha(result_path),source_saved_time=500.)
        cases.append(c)
    assert len(cases)==2
    hashes=p.n.sources()
    for path in [Path(__file__),ROOT/'code/mechanisms/revision_five_model_compensation_20261019.py',ROOT/'code/mechanisms/revision_five_model_compensation_campaign_20261019.py']:
        hashes[str(path.relative_to(ROOT))]=p.n.sha(path)
    for rel in hashes:
        dest=out/'frozen_sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
    control=next(c for c in m['cases'] if c['model']=='debelly' and c['length']==39 and c['case_kind']=='zero_control')
    controlpath=ROOT/control['exact_reuse']['source_result'] if control['exact_reuse'] else BASE/'workers'/f'{control["slot"]:04d}'/'result.json'
    cr=p.read(controlpath);assert next(w for w in cr['withdrawals'] if w['duration']==0.)['no_cue_preparation_valid']
    manifest=dict(created_at=p.n.now(),cases=cases,models=['debelly'],native_ON_integrations=0,new_native_OFF_branches=2,
        cue_widths_nominal=[1.5,5.],requested_peak=3.6,withdrawal_time=500.,off_duration=10000.,
        source_baseline_manifest_sha256=p.n.sha(BASE/'manifest.json'),source_no_cue_control=str(controlpath.relative_to(ROOT)),source_no_cue_control_sha256=p.n.sha(controlpath),
        purpose='Test whether one common pulse links observed different weakening latencies to different retained outcomes. No optimized minimum or monotonic duration inference.',
        scientific_native_engine='Accepted frozen platform sampled_run/gates, independent endpoints, unchanged half-unit stagger, all11states restored.',
        resource_forecast='Two one-CPU workers,2GiB/10min each; approximately150MB allocated output,nonbinding.')
    p.write(out/'manifest.json',manifest);p.write(out/'provenance.json',dict(manifest_sha256=p.n.sha(out/'manifest.json'),frozen_sources=hashes,frozen_inputs=inputs))
    shell=out/'worker.sh';shell.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(out/'frozen_sources/code/mechanisms'/Path(__file__).name)+' worker --out '+str(out)+' --task "$SLURM_ARRAY_TASK_ID"\n');shell.chmod(0o750)
    cmd='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:10:00 --array=0-1 --job-name=F2-DB-T500-pair --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(shell)
    (out/'submission_command.txt').write_text(cmd+'\n');print(p.n.sha(out/'manifest.json'));print(cmd)


def worker(out,task):
    start=time.monotonic();m=p.verify(out);c=m['cases'][task];folder=out/'workers'/f'{task:02d}';folder.mkdir(exist_ok=True);assert not (folder/'result.json').exists()
    _,prep=p.n.read_checkpoint(out/c['frozen_checkpoint'],'debelly')
    with np.load(out/c['frozen_restart']) as a:restart={name:a['state_'+name].copy() for name in p.n.DB_NAMES}
    with np.load(out/c['regional_input']) as a:z={k:a[k].copy() for k in a.files}
    p.write(folder/'status.json',dict(state='off_continuation',generated_at=p.n.now()))
    off=p.sampled_run(c,restart,p.e.withdrawal(c),p.sample_times('debelly',10000.,off=True))
    first=p.state_at(off,'debelly',0);roundoff={name:float(np.max(abs(first[name]-restart[name]))) for name in p.n.DB_NAMES}
    assert max(roundoff.values())==0
    path=folder/'off.npz';p.save_phase(path,off);a,events=campaign.regional(off,c,z,'after_withdrawal');reg=folder/'regional.npz';np.savez_compressed(reg,**a)
    p.write(folder/'result.json',dict(case=c,manifest_sha256=p.n.sha(out/'manifest.json'),generated_at=p.n.now(),
        off_archive=str(path.relative_to(out)),off_archive_sha256=p.n.sha(path),regional_archive=str(reg.relative_to(out)),regional_archive_sha256=p.n.sha(reg),
        initial_restart_max_absolute_errors=roundoff,events=events,diagnostics=p.gates(off,c,prep),
        regional_endpoint={name:float(v[-1]) for name,v in a.items() if name!='time'},wall_seconds=time.monotonic()-start,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    p.write(folder/'status.json',dict(state='complete',generated_at=p.n.now()))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','worker']);parser.add_argument('--out',type=Path,default=DEFAULT);parser.add_argument('--task',type=int)
    a=parser.parse_args()
    if a.action=='init':initialize(a.out.resolve())
    else:
        try:worker(a.out.resolve(),a.task)
        except Exception:
            folder=a.out/'workers'/f'{a.task:02d}';folder.mkdir(exist_ok=True);p.write(folder/'error.json',dict(traceback=traceback.format_exc()));raise
