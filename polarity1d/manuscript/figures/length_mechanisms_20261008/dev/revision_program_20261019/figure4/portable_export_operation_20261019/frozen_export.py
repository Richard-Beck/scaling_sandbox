"""Exact compact F4 redraw exports; native archives remain replay lineage.

Field CSV.gz chunks round-trip every plotted float64 value. Only declared
condition medoids and the four main mean controls are exported as fields.
"""
from pathlib import Path
from datetime import datetime,timezone
import gzip
import hashlib
import importlib.util
import io
import json
import os
import shutil
import time
import numpy as np
import pandas as pd

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript';FIG=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
SOURCE=FIG/'signature_assembly_20261019';OUT=FIG/'portable_reconstruction_20261019'
RUN=MS/'data/runs/revision_noise_strength_replication_scoped_20261019'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def main():
    assert os.environ.get('SLURM_JOB_ID'),'Whole atlas float export runs on SLURM'
    start=time.monotonic();core=json.loads((SOURCE/'supporting/metadata/FINAL_POSTPROCESS_RECEIPT.json').read_text());assert core['status']=='complete'
    for name in ('fields','tables','metadata','helpers','selection_cache/results/partial'):(OUT/name).mkdir(parents=True,exist_ok=True)
    assert not (OUT/'metadata/EXPORT_RECEIPT.json').exists()
    for p in (SOURCE/'supporting/data').glob('*.csv'):shutil.copy2(p,OUT/'tables'/p.name)
    for name in ('native_generation_snapshot.json','BOOTSTRAP_PROTOCOL.json','FINAL_POSTPROCESS_RECEIPT.json','FINAL_REPORT_RECEIPT.json'):
        shutil.copy2(SOURCE/'supporting/metadata'/name,OUT/'metadata'/name)
    cachepath=RUN/'results/partial/case_cache.json';cache=json.loads(cachepath.read_text());assert len(cache)==936
    selected=pd.read_csv(SOURCE/'supporting/data/all_condition_native_atlas_selection.csv').case_id.tolist()
    context=pd.read_csv(SOURCE/'supporting/data/local_example_context.csv')
    for _,row in context.iterrows():
        match=[k for k,v in cache.items() if v['profile_sha256']==row.mean_profile_sha256];assert len(match)==1;selected+=match
    selected=list(dict.fromkeys(selected));assert len(selected)==112
    compact={}
    needed=('window','model','platform','size','fraction','kind','seed','case_id','mean_R','mean_A','mean_Q','mean_contrast_ratio','qualified_multiregion_fraction','weak_R_fraction','low_contrast_fraction')
    for case_id,case in cache.items():
        d=next(d for d in case['descriptors'] if d['window']=='last_half')
        compact[case_id]={k:case[k] for k in ('profile','trace','profile_sha256','trace_sha256')}
        compact[case_id].update(descriptors=[{k:d[k] for k in needed}],portable_case_id=case_id)
    write(OUT/'selection_cache/results/partial/case_cache.json',compact)
    shutil.copy2(SOURCE/'supporting/metadata/native_generation_snapshot.json',OUT/'selection_cache/results/partial/status.json')
    exported={};values=0
    for ordinal,case_id in enumerate(selected):
        case=cache[case_id]
        for key in ('profile','trace'):assert sha(case[key])==case[key+'_sha256']
        with np.load(case['profile']) as z:
            t=z['time'];keep=np.r_[0,np.flatnonzero(t>=t[-1]/2)]
            times=t[keep];u=z['primary'][keep];C=np.broadcast_to(z['circumference'],len(t))[keep]
            coordinate=z['material_positions'] if 'material_positions' in z else z['x']
        with np.load(case['trace']) as z:
            ix=np.searchsorted(z['time'],times);np.testing.assert_array_equal(z['time'][ix],times)
            R=np.hypot(z['P_real'][ix],z['P_imag'][ix]);A=z['M0'][ix]/z['M0'][0];Q=z['q1'][ix]/z['q1'][0]
            regions=z['current_half_contrast_regions'][ix];contrast=z['contrast'][ix]/z['contrast'][0]
        matrix=np.column_stack((times,C,R,A,Q,regions,contrast,u));name=f'field_{ordinal:03d}.csv.gz';target=OUT/'fields'/name
        columns=['time','C','R','A','Q','regions','contrast_ratio']+[f'u_{j:04d}' for j in range(u.shape[1])]
        with target.open('wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as zipped:
                with io.TextIOWrapper(zipped,encoding='utf-8',newline='\n') as text:
                    np.savetxt(text,matrix,fmt='%.17g',delimiter=',',header=','.join(columns),comments='')
        with gzip.open(target,'rt') as handle:roundtrip=np.loadtxt(handle,delimiter=',',skiprows=1)
        np.testing.assert_array_equal(roundtrip,matrix);values+=matrix.size
        exported[case_id]=dict(path='fields/'+name,sha256=sha(target),shape=list(matrix.shape),coordinate=coordinate.tolist(),
            profile_cadence=40 if case['descriptors'][0]['platform']=='chemical' else 4,
            selection='Exact basal row plus every saved second-half profile; original native amplitudes; no decimation',
            lineage={k:case[k] for k in ('profile','trace','profile_sha256','trace_sha256')})
        if ordinal%12==0:print(json.dumps(dict(fields=ordinal+1,elapsed_seconds=time.monotonic()-start)),flush=True)
    write(OUT/'metadata/FIELD_INDEX.json',exported)
    frozen=SOURCE/'supporting/metadata/final_frozen_sources'
    for name in ('local.py','temporal.py','organization.py','quarters.py','assembly.py'):
        shutil.copy2(frozen/name,OUT/'helpers'/name)
    for source,name in [('revision_figure4_compact_main_20261019.py','compact_main.py'),('revision_figure4_portable_redraw_20261019.py','redraw.py')]:
        shutil.copy2(MS/'code/mechanisms'/source,OUT/'helpers'/name)
    expected={p.name:sha(p) for p in SOURCE.glob('*.png')};compactmain=FIG/'main_plate_20261019/figure4_main.png'
    if compactmain.exists():expected['figure4_main.png']=sha(compactmain)
    write(OUT/'metadata/EXPECTED_PNG_HASHES.json',expected)
    receipt=dict(generated_at=datetime.now(timezone.utc).isoformat(),native_integrations=0,source_sha256=sha(__file__),
        source_core_receipt_sha256=sha(SOURCE/'supporting/metadata/FINAL_POSTPROCESS_RECEIPT.json'),
        exact_fields=112,condition_medoid_fields=108,main_mean_control_fields=4,float64_values_roundtrip_verified=values,
        cache_cases=936,roundtrip_maximum_error=0.,
        redraw_dependencies='Only fields/*.csv.gz, tables/*.csv, metadata/*.json, helpers/*.py and selection_cache JSON. Absolute native paths in lineage are not opened by redraw.',
        scientific_replay='Recomputing all native moments/region geometry requires preserved local native archives; redraw is a distinct verified closure.',
        allocated_bytes=sum(p.stat().st_blocks*512 for p in OUT.rglob('*') if p.is_file()),
        outputs={str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file()})
    write(OUT/'metadata/EXPORT_RECEIPT.json',receipt);print(json.dumps({k:v for k,v in receipt.items() if k!='outputs'},indent=2))


if __name__=='__main__':main()
