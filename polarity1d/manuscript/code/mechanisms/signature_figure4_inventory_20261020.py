"""Canonical read-only E0 inputs; no simulation or hypothesis analysis."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
PACKAGE=ROOT/'polarity1d/manuscript'
SCOPED=PACKAGE/'data/runs/revision_noise_strength_replication_scoped_20261019'
REFINEMENT=PACKAGE/'data/runs/revision_gory_strength_transition_20261019'


def read(path):return json.loads(Path(path).read_text())
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def basal_profile(path,model,C_label):
    with np.load(path,allow_pickle=False) as archive:
        state=archive['state'].copy()
        if model=='spring_ring':
            q=archive['material_positions'].copy();n=len(q)
            if state.shape!=(2*n+1,):raise ValueError('Spring ring requires blocked active[N], inactive[N], circumference state')
            active=state[:n];C=float(state[-1])
        else:
            active=state[0];C=float(C_label);q=archive['x'].copy()/C
    if active.ndim!=1 or len(active)!=len(q) or C<=0 or not np.isfinite(active).all():
        raise ValueError('Invalid exact basal profile')
    return active,q,C


def prepare(out):
    out=Path(out).resolve()
    if out.exists() or not out.is_relative_to(PACKAGE/'data/runs'):
        raise ValueError('Choose a fresh E0 directory beneath manuscript/data/runs')
    (out/'inputs').mkdir(parents=True);(out/'results').mkdir()
    scoped=read(SCOPED/'manifest.json');cache_path=SCOPED/'results/partial/case_cache.json'
    cache=read(cache_path);refined=read(REFINEMENT/'manifest.json');basals={}
    for row in scoped['basals']:
        model=row['model'];label=float(row['size']);checkpoint=SCOPED/row['input']
        if sha(checkpoint)!=row['input_sha256']:raise ValueError('Exact basal checkpoint changed')
        active,q,C=basal_profile(checkpoint,model,label)
        destination=out/'inputs/basals'/f'{model}_C{label:g}.npz';destination.parent.mkdir(exist_ok=True)
        np.savez_compressed(destination,primary=active,material_positions=q,circumference=C,
            source_checkpoint_sha256=np.array(row['input_sha256']))
        basals[(model,label)]=dict(basal_checkpoint_path=str(checkpoint),basal_checkpoint_sha256=sha(checkpoint),
            basal_profile_path=str(destination),basal_profile_sha256=sha(destination),
            basal_contrast=float(np.ptp(active)),actualC0=C,native_cells=len(active))
    candidates=[]
    for identifier,item in cache.items():
        meta=item['descriptors'][0]
        candidates.append(dict(id=identifier,model=meta['model'],C_label=float(meta['size']),
            f=float(meta['fraction']),kind='zero' if meta['kind']=='basal' else meta['kind'],
            seed=meta['seed'],scalar_path=item['trace'],scalar_sha256=item['trace_sha256'],
            profile_path=item['profile'],profile_sha256=item['profile_sha256'],origin='core'))
    for row in refined['jobs']:
        worker=REFINEMENT/'workers'/f"{row['slot']:04d}";status=read(worker/'status.json')
        if not status['complete']:raise ValueError('Refinement native worker is incomplete')
        candidates.append(dict(id=f"refinement_{row['slot']:04d}",model=row['model'],C_label=float(row['size']),
            f=float(row['fraction']),kind=row['kind'],seed=row['seed'],
            scalar_path=str(worker/'trace.npz'),scalar_sha256=status['output_sha256']['trace.npz'],
            profile_path=str(worker/'primary.npz'),profile_sha256=status['output_sha256']['primary.npz'],origin='refinement'))
    records=[];identities={};paths={};duplicates=[]
    for row in candidates:
        scalar=Path(row['scalar_path']).resolve();profile=Path(row['profile_path']).resolve()
        key=(row['model'],row['C_label'],row['f'],row['kind'],row['seed'] if row['kind']=='noise' else None)
        if key in identities or scalar in paths:
            previous=identities.get(key,paths.get(scalar))
            if previous['scalar_sha256']!=row['scalar_sha256']:
                raise ValueError('A claimed independent realization has conflicting native outputs')
            duplicates.append(dict(duplicate_id=row['id'],retained_id=previous['id']));continue
        if not scalar.is_file() or not profile.is_file():raise FileNotFoundError('Canonical native scalar/profile is missing')
        status_path=scalar.parent/'status.json';status=read(status_path)
        if status.get('complete') is not True:raise ValueError('Core native worker is incomplete')
        row.update(scalar_path=str(scalar),profile_path=str(profile),
            sigma=1. if row['kind']=='noise' else 0.,
            source_status_path=str(status_path),source_status_sha256=sha(status_path),
            source_manifest_sha256=status.get('manifest_sha256'),
            scalar_cadence=.5 if row['model']=='spring_ring' else 5.,
            profile_cadence=4. if row['model']=='spring_ring' else 40.,
            native_horizon=3200. if row['model']=='spring_ring' else 32000.,
            **basals[(row['model'],row['C_label'])])
        records.append(row);identities[key]=row;paths[scalar]=row
    records.sort(key=lambda r:(r['model'],r['C_label'],r['kind'],r['f'],-1 if r['seed'] is None else r['seed']))
    basal_initial_checks=[]
    for row in records:
        if row['kind']!='zero':continue
        with np.load(row['basal_profile_path'],allow_pickle=False) as z:active=z['primary'].copy()
        with np.load(row['profile_path'],allow_pickle=False) as z:native=z['primary'][0].copy()
        if not np.array_equal(active,native):raise ValueError('Compact basal must exactly match the native first primary profile')
        basal_initial_checks.append(dict(model=row['model'],C_label=row['C_label'],exact_native_initial_primary=True))
    counts=dict(Counter(row['kind'] for row in records))
    if counts!={'noise':840,'mean':110,'zero':12}:
        raise ValueError('E0 independent native case census differs from the declared 936+26 authority')
    sources={str(p):sha(p) for p in (SCOPED/'manifest.json',cache_path,REFINEMENT/'manifest.json',
        SCOPED/scoped['reuse']['index'])}
    for name,path in [('scoped_manifest.json',SCOPED/'manifest.json'),('source_case_cache.json',cache_path),
                      ('refinement_manifest.json',REFINEMENT/'manifest.json')]:shutil.copy2(path,out/'inputs'/name)
    write(out/'inventory.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=sources,inventory_source_sha256=sha(Path(__file__)),
        cases=records,counts=counts,total_independent_cases=len(records),duplicate_reuses=duplicates,
        basal_native_initial_checks=basal_initial_checks,
        raw_native_file_rehash_complete=False,
        scientific_status='Canonical expected hashes from native output receipts; full rehash is a separate gate.',
        region_qualification='Current contrast >=0.1 times this exact basal contrast; raw half-contrast components.',
        R_definition='abs(P), P=M1/M1(initial), as retained by the existing native pipeline; not q1 coherence.',
        no_native_integrations=True))
    with (out/'inventory.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    write(out/'freeze.json',dict(inventory_sha256=sha(out/'inventory.json'),csv_sha256=sha(out/'inventory.csv')))
    return dict(output=str(out),counts=counts,total_independent_cases=len(records),duplicate_reuses=len(duplicates),
        inventory_sha256=sha(out/'inventory.json'),raw_native_file_rehash_complete=False)


def verify(out):
    out=Path(out);inventory=read(out/'inventory.json')
    if sha(out/'inventory.json')!=read(out/'freeze.json')['inventory_sha256']:
        raise ValueError('Canonical E0 case inventory changed')
    checked={}
    for row in inventory['cases']:
        for role in ('scalar','profile','basal_checkpoint','basal_profile','source_status'):
            path=Path(row[role+'_path']);expected=row[role+'_sha256']
            if path not in checked:checked[path]=sha(path)
            if checked[path]!=expected:raise ValueError('Native E0 file changed: '+str(path))
    result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,
        inventory_sha256=sha(out/'inventory.json'),unique_files_rehashed=len(checked),
        cases=len(inventory['cases']),counts=inventory['counts'],no_native_integrations=True)
    write(out/'results/native_identity_verification.json',result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for command in ('prepare','verify'):sub.add_parser(command).add_argument('--out',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(prepare(args.out) if args.command=='prepare' else verify(args.out),indent=2))
