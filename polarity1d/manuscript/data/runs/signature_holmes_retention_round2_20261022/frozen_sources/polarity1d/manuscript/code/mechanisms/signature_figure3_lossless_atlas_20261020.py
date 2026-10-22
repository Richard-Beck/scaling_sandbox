"""Lossless, pixel-audited portable images for the complete native atlas."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil
from datetime import datetime,timezone
import numpy as np
from PIL import Image

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/"polarity1d/manuscript"
OUT=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance"
DEFAULT=ROOT/"data/runs/signature_figure3_lossless_atlas_20261020"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def freeze(run,out):
    if run.exists():raise FileExistsError("A fresh immutable conversion run is required")
    for directory in ['workers','results','logs']:(run/directory).mkdir(parents=True)
    shutil.copy2(__file__,run/"frozen_converter.py")
    index=out/"supporting/data/native_history_index.json"
    cases=read(index)["cases"]
    jobs=[dict(case_id=c["case_id"],png=c["png"],png_sha256=c["png_sha256"],
        webp=str(out/"supporting/viewers/native_histories_webp"/(c["case_id"]+".webp"))) for c in cases]
    chunks=[list(range(i,min(i+20,len(jobs)))) for i in range(0,len(jobs),20)]
    write(run/"manifest.json",dict(jobs=jobs,chunks=chunks,source_index_sha256=sha(index),
        source_sha256=sha(run/"frozen_converter.py"),format="WebP lossless, method6, exact RGBA pixel equality against immutable PNG"))
    command="exec "+shlex.join(['env','PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',
        '/home/4473331/.conda/envs/scaling_sandbox/bin/python',str(run/"frozen_converter.py"),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    payload=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G','--time=00:30:00',f'--array=0-{len(chunks)-1}',
        '--job-name=figure3-lossless-atlas','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),'--wrap',command]
    write(run/'slurm_payload.json',dict(jobs=len(chunks),images=len(jobs),submitted=False,sbatch_arguments=payload))
    return dict(images=len(jobs),jobs=len(chunks))


def worker(run,slot):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Complete image conversion uses SLURM')
    plan=read(run/'manifest.json');assert sha(Path(__file__))==plan['source_sha256']
    rows=[]
    for index in plan['chunks'][slot]:
        job=plan['jobs'][index];source,target=Path(job['png']),Path(job['webp'])
        assert sha(source)==job['png_sha256'];target.parent.mkdir(exist_ok=True)
        with Image.open(source) as image:
            reference=np.asarray(image.convert('RGBA'))
            image.save(target,format='WEBP',lossless=True,method=6,exact=True)
        with Image.open(target) as image:converted=np.asarray(image.convert('RGBA'))
        assert np.array_equal(reference,converted),'Lossless pixel audit failed: '+job['case_id']
        rows.append(dict(**job,webp_sha256=sha(target),pixel_equal=True,png_bytes=source.stat().st_size,webp_bytes=target.stat().st_size))
    write(run/'workers'/f'{slot:04d}.json',dict(complete=True,records=rows))


def collect(run):
    plan=read(run/'manifest.json');rows=[]
    for slot,indices in enumerate(plan['chunks']):
        record=read(run/'workers'/f'{slot:04d}.json');assert record['complete']
        assert len(indices)==len(record['records'])
        for index,row in zip(indices,record['records']):
            assert plan['jobs'][index]['case_id']==row['case_id'] and row['pixel_equal']
            assert sha(Path(row['webp']))==row['webp_sha256']
            rows.append(row)
    assert len(rows)==len(plan['jobs'])==len({r['case_id'] for r in rows})
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),complete=True,pixel_equal_images=len(rows),
        png_bytes=sum(r['png_bytes'] for r in rows),webp_bytes=sum(r['webp_bytes'] for r in rows),
        source_manifest_sha256=sha(run/'manifest.json'),converter_sha256=plan['source_sha256'],records=rows)
    write(run/'results/lossless_pixel_audit.json',report)
    write(OUT/'supporting/metadata/lossless_pixel_audit.json',report)
    return {key:value for key,value in report.items() if key!='records'}


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','worker','collect']);p.add_argument('--run',type=Path,default=DEFAULT)
    p.add_argument('--output',type=Path,default=OUT);p.add_argument('--slot',type=int);args=p.parse_args()
    result=freeze(args.run.resolve(),args.output.resolve()) if args.action=='freeze' else worker(args.run.resolve(),args.slot) if args.action=='worker' else collect(args.run.resolve())
    if result is not None:print(json.dumps(result))


if __name__=='__main__':main()
