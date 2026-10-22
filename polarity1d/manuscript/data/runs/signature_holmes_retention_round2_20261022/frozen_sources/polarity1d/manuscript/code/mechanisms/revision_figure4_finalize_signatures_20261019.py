"""Guarded read-only F4 final analysis after terminal replication collection.

Never submits native work, alters source inputs or relabels numerical evidence.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
RUN=MS/'data/runs/revision_noise_strength_replication_scoped_20261019'
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/signature_assembly_20261019'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=OUT);a=parser.parse_args();out=a.out
    assert os.environ.get('SLURM_JOB_ID'),'Full read-only assembly runs on SLURM'
    metadata=out/'supporting/metadata';frozen=metadata/'final_frozen_sources';closure=read(metadata/'FINAL_ANALYSIS_SOURCE_CLOSURE.json')
    for name,digest in closure['source_sha256'].items():assert sha(frozen/name)==digest
    assert sha(RUN/'manifest.json')==closure['native_manifest_sha256']
    assert sha(RUN/'results/STABLE_ESTIMAND_ACCEPTANCE.json')==closure['scope_acceptance_sha256']
    status=read(RUN/'results/partial/status.json')
    assert status['complete'] and status['producer']=='terminal' and status['failures']==0
    assert status['scientific_cases_available']==936 and status['new_completed_outputs']==360
    assert not (metadata/'FINAL_POSTPROCESS_RECEIPT.json').exists(),'Do not automatically regenerate a delivered package'
    python=sys.executable
    stages=[
        ('replicated_bootstrap','assembly.py',['--run',str(RUN),'--out',str(out)]),
        ('native_spatial_windows','spatial.py',['--run',str(RUN),'--expected-cases','936','--out',str(out)]),
        ('native_examples','local.py',['--out',str(out)]),
        ('all_temporal_windows','temporal.py',['--out',str(out),'--source',str(out/'supporting/data/memory_all_conditions_windows_lags.csv')]),
        ('all_organization_windows','organization.py',['--out',str(out),'--conditions',str(out/'supporting/data/condition_window_summaries.csv')]),
        ('late_quarter_sensitivity','quarters.py',['--out',str(out)]),
        ('report_legend_package','report.py',['--out',str(out)])]
    progress=dict(started_at=now(),native_integrations=0,slurm_job_id=os.environ['SLURM_JOB_ID'],stages=[],status='running',native_generation=status['generated_at'])
    write(metadata/'FINAL_POSTPROCESS_STATUS.json',progress)
    try:
        for label,source,args in stages:
            started=now();print(json.dumps(dict(stage=label,started_at=started)),flush=True)
            log=metadata/('final_'+label+'.log')
            with log.open('w') as handle:
                subprocess.run([python,str(frozen/source),*args],cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT,check=True)
            progress['stages'].append(dict(stage=label,started_at=started,finished_at=now(),source_sha256=sha(frozen/source),log=str(log.relative_to(out))))
            write(metadata/'FINAL_POSTPROCESS_STATUS.json',progress)
        progress.update(status='complete',finished_at=now(),independent_scientific_review='pending')
        write(metadata/'FINAL_POSTPROCESS_STATUS.json',progress)
        files={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.suffix not in ('.out','.err','.log') and p.name!='FINAL_POSTPROCESS_RECEIPT.json'}
        write(metadata/'FINAL_POSTPROCESS_RECEIPT.json',dict(**progress,outputs=files,source_closure_sha256=sha(metadata/'FINAL_ANALYSIS_SOURCE_CLOSURE.json')))
        print(json.dumps(dict(status='complete',out=str(out),finished_at=now(),native_integrations=0)),flush=True)
    except BaseException as exc:
        progress.update(status='failed',stopped_at=now(),error=f'{type(exc).__name__}: {exc}')
        write(metadata/'FINAL_POSTPROCESS_STATUS.json',progress);raise


if __name__=='__main__':main()
