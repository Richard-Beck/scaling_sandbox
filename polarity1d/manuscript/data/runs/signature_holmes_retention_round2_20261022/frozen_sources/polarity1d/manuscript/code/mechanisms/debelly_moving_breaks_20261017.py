#!/usr/bin/env python3
"""One fixed-grid moving-connectivity-defect prototype; preparation only."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode=True
import numpy as np

PROJECT=next(p for p in Path(__file__).resolve().parents
             if (p/'polarity1d/manuscript/AGENTS.md').is_file())
FIXED_RUN=PROJECT/'polarity1d/manuscript/data/runs/debelly_rac_continuity_20261017'
RUN=PROJECT/'polarity1d/manuscript/data/runs/debelly_moving_breaks_complete_ledgers_20261017'
DEV=PROJECT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_moving_breaks_20261017'
DEP=Path(__file__).resolve().parent
if not (DEP/'ring_core.py').exists():DEP=FIXED_RUN/'frozen_sources'
sys.path.insert(0,str(DEP))
from ring_core import FIELDS, RingCore, face_trace
from gap_events import geometry, inventory


def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference():
    path=DEP/('reference_pilot.py' if (DEP/'reference_pilot.py').exists() else 'pilot.py')
    spec=importlib.util.spec_from_file_location('moving_break_reference',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    return module


REF=reference()


def circular_distance(a,b,C):return abs((a-b+.5*C)%C-.5*C)


def nearest_face(position,edges,previous=None):
    C=float(edges[-1]);x=float(position%C)
    distances=abs(edges[:-1]-x);distances=np.minimum(distances,C-distances)
    tolerance=64*np.finfo(float).eps*max(1.,C)
    candidates=np.flatnonzero(distances<=distances.min()+tolerance).tolist()
    if previous is not None and previous in candidates:return int(previous)
    # A strict cell-center maximum lies exactly between two native faces.
    # Remaining ties choose the clockwise face in this declared orientation.
    return min(candidates,key=lambda k:float((edges[k]-x)%C))


def moving_transition(widths,state,time,threshold,next_id):
    widths=np.asarray(widths,dtype=float)
    if not np.allclose(widths,widths[0],rtol=0.,atol=1e-12):
        raise RuntimeError('Moving prototype requires the unchanged uniform mesh')
    edges,_=geometry(widths);C=float(edges[-1]);old_fields={q:state[q] for q in FIELDS}
    mass0=float(np.dot(widths,state['b']+state['c']))
    state={**state,'cuts':deepcopy(state.get('cuts',{}))}
    records=[]
    for face,gap in list(state['cuts'].items()):
        sensor=face_trace(state['R'],widths,face)
        if sensor<threshold:
            records.append(dict(event='close',time=float(time),id=int(gap['id']),
                x=float(gap['x']),Rac_trace=sensor,controller=deepcopy(gap)))
            del state['cuts'][face]
    info=inventory(widths,state,threshold)
    by_id={int(gap['id']):(face,gap) for face,gap in state['cuts'].items()}
    proposals=[];suppressed=[]
    for component in info['components']:
        candidates=sorted(component['peaks'],key=lambda p:(-p['maximum'],p['position']))
        identifiers=component['existing_gap_ids']
        if len(identifiers)>1:
            raise RuntimeError('Moving-break merger has multiple surviving breaks; no deletion/identity law invented')
        if not candidates:continue
        old_face,old_gap=by_id[identifiers[0]] if identifiers else (None,None)
        tolerance=64*np.finfo(float).eps*max(1.,abs(candidates[0]['maximum']))
        tied=[p for p in candidates if candidates[0]['maximum']-p['maximum']<=tolerance]
        if len(tied)>1:
            if old_gap is None:
                raise RuntimeError('Equal distinct highest maxima in gap-free region; no birth identity chosen')
            distances=[circular_distance(p['position'],old_gap['x'],C) for p in tied]
            winners=[p for p,d in zip(tied,distances) if d<=min(distances)+tolerance]
            if len(winners)>1:
                raise RuntimeError('Equal highest maxima have ambiguous moving-gap association')
            peak=winners[0]
        else:peak=candidates[0]
        face=nearest_face(peak['position'],edges,old_face)
        proposals.append((face,old_face,old_gap,peak))
        suppressed.extend(dict(position=p['position'],value=p['maximum'],
            reason='one break per occupied above-threshold region') for p in candidates if p is not peak)
    destinations=[face for face,_,_,_ in proposals]
    if len(destinations)!=len(set(destinations)):
        raise RuntimeError('Two moving regions nominate the same face; no collision rule invented')
    original_cuts=deepcopy(state['cuts']);newcuts=deepcopy(state['cuts'])
    move_sham_cuts=deepcopy(state['cuts'])
    for face,old_face,gap,peak in proposals:
        x=float(edges[face]);quantization=circular_distance(x,peak['position'],C)
        if gap is None:
            if face in newcuts:raise RuntimeError('Birth destination is occupied')
            value=face_trace(state['u'],widths,face)
            newgap=dict(id=int(next_id),x=x,birth_x=x,birth_time=float(time),
                Uleft=value,Uright=value,raw_left=value,raw_right=value)
            newcuts[face]=newgap;move_sham_cuts[face]=deepcopy(newgap)
            records.append(dict(event='open',time=float(time),id=int(next_id),x=x,
                detected_peak_x=float(peak['position']),detected_peak_Rac=float(peak['maximum']),
                inserted_face_Rac=face_trace(state['R'],widths,face),
                quantization_error=quantization,initialized_U=value,removed_material=0.,cells_after=len(widths)))
            next_id+=1
        elif face!=old_face:
            if face in original_cuts and original_cuts[face]['id']!=gap['id']:
                raise RuntimeError('Relocation target is occupied; no break merger rule invented')
            carried=deepcopy(gap);carried['x']=x;carried.setdefault('birth_x',gap['x'])
            del newcuts[old_face];newcuts[face]=carried
            background=face_trace(state['u'],widths,face)
            records.append(dict(event='move',time=float(time),id=int(gap['id']),
                old_face=int(old_face),new_face=int(face),old_x=float(gap['x']),x=x,
                detected_peak_x=float(peak['position']),detected_peak_Rac=float(peak['maximum']),
                old_gap_Rac=face_trace(state['R'],widths,old_face),
                new_gap_Rac=face_trace(state['R'],widths,face),quantization_error=quantization,
                signed_reference_shift=float((x-gap['x']+.5*C)%C-.5*C),
                old_background_u=face_trace(state['u'],widths,old_face),new_background_u=background,
                carried_controller=deepcopy(gap),
                carried_center_minus_new_background=.5*(gap['Uleft']+gap['Uright'])-background,
                carried_flank_jump=gap['Uleft']-gap['Uright'],removed_material=0.))
    state['cuts']=newcuts
    mass1=float(np.dot(widths,state['b']+state['c']))
    if any(not np.array_equal(state[q],old_fields[q]) for q in FIELDS):
        raise RuntimeError('Topology event projected a stored field')
    return widths,state,next_id,dict(events=records,distinct_maximum_suppressions=suppressed,
        material_impulse=mass1-mass0,move_sham_cuts=move_sham_cuts,prebirth_inventory=info)


class MovingEngine(REF.Engine):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.motion_ledgers=[];self.first_move_fields=None;self.worst_move_fields=None
        self.largest_move_difference=-1.

    def step(self,location,cue_duration,dt):
        before=self.state
        self.widths,self.state,self.next_id,event=moving_transition(
            self.widths,self.state,self.time,self.cfg['activity_threshold'],self.next_id)
        self.events.extend(event['events'])
        for row in event['distinct_maximum_suppressions']:
            self.suppressions.append(dict(time=self.time,**row))
        self.max_errors['event_material_impulse']=max(self.max_errors['event_material_impulse'],abs(event['material_impulse']))
        external=self.cue(location,cue_duration)
        new,ledger=self.core.step(self.state,dt,external,target_mass=self.cfg['total_MCA'])
        moved=[row for row in event['events'] if row['event']=='move']
        if moved:
            sham,sham_ledger=self.core.step(self.state,dt,external,target_mass=self.cfg['total_MCA'],
                                          cuts=event['move_sham_cuts'])
            difference={q:float(np.max(abs(new[q]-sham[q]))) for q in FIELDS}
            score=max(difference[q]/self.cfg['scales'][q] for q in FIELDS)
            motion=dict(time=self.time,step=dt,events=deepcopy(moved),conditional_one_step_field_difference=difference,
                maximum_scaled_one_step_difference=score,
                actual_budget={k:float(v) for k,v in ledger.items() if isinstance(v,(int,float,np.number))},
                frozen_position_budget={k:float(v) for k,v in sham_ledger.items() if isinstance(v,(int,float,np.number))},
                actual_cut_ledgers=deepcopy(ledger['cuts']),frozen_cut_ledgers=deepcopy(sham_ledger['cuts']))
            self.motion_ledgers.append(motion)
            snapshot=dict(time=np.asarray(self.time),step=np.asarray(dt),after_time=np.asarray(self.time+dt),widths=self.widths.copy(),
                **{f'before_{q}':before[q].copy() for q in FIELDS},
                **{f'actual_{q}':new[q].copy() for q in FIELDS},
                **{f'frozen_{q}':sham[q].copy() for q in FIELDS})
            if self.first_move_fields is None:self.first_move_fields=snapshot
            if score>self.largest_move_difference:self.largest_move_difference=score;self.worst_move_fields=snapshot
        self.state=new
        for name in self.max_errors:
            if name in ledger:self.max_errors[name]=max(self.max_errors[name],abs(float(ledger[name])))
        for name in self.integrated:self.integrated[name]+=float(ledger[name])
        for row in event['events']:
            if row['event']=='open':
                face=next(k for k,v in self.state['cuts'].items() if v['id']==row['id'])
                self.first_step_by_event[str(row['id'])]=deepcopy(ledger['cuts'][face])
        self.time+=dt;self.last_ledger=ledger

    def dump(self,folder,status):
        result=super().dump(folder,status)
        write(folder/'motion_ledgers.json',self.motion_ledgers)
        for name,snapshot in [('first_move_fields',self.first_move_fields),('worst_move_fields',self.worst_move_fields)]:
            if snapshot is not None:np.savez_compressed(folder/(name+'.npz'),**snapshot)
        result.update(prototype='moving native-face connectivity defect',
            event_counts={name:sum(row['event']==name for row in self.events) for name in ['open','close','move']},
            largest_scaled_one_step_relocation_difference=self.largest_move_difference)
        write(folder/'summary.json',result)
        return result


def initialize(run):
    if run.exists():raise RuntimeError('Refusing to overwrite an existing moving-break run')
    source_manifest=json.loads((FIXED_RUN/'manifest.json').read_text())
    for path,digest in source_manifest['sha256'].items():
        if sha(FIXED_RUN/path)!=digest:raise RuntimeError('Fixed comparator input hash changed: '+path)
    run.mkdir(parents=True);source=run/'frozen_sources';source.mkdir();(run/'logs').mkdir()
    cfg=json.loads((FIXED_RUN/'inputs.json').read_text())
    cfg.update(gap_rule='one moving native-face break per occupied Rac>0.8 region; follow highest maximum; close local trace below0.8; carry controller unchanged',
        motion_rule='nearest existing face; face-distance ties prefer previous location then clockwise; no remeshing',
        peak_tie_rule='occupied region: uniquely nearest equal-high maximum to old break; unoccupied or unresolved association tie stops',
        controller_rule='carry Uleft,Uright,raw_left,raw_right unchanged on relocation; no stored field projection',
        reference_fixed_run=str(FIXED_RUN),tasks=[],prospective_holdout_extension_allowed=False)
    write(run/'inputs.json',cfg)
    for name in ['ring_core.py','peak_detector.py','gap_events.py','debelly.py','polarity1d/contract.py','polarity1d/protocols.py']:
        dest=source/name;dest.parent.mkdir(exist_ok=True,parents=True)
        shutil.copy2(FIXED_RUN/'frozen_sources'/name,dest)
    shutil.copy2(FIXED_RUN/'frozen_sources/pilot.py',source/'reference_pilot.py')
    shutil.copy2(Path(__file__).resolve(),source/'moving_pilot.py')
    paths=['inputs.json',*[str(p.relative_to(run)) for p in source.rglob('*.py')]]
    write(run/'manifest.json',dict(sha256={p:sha(run/p) for p in paths},
        fixed_run_manifest_sha256=sha(FIXED_RUN/'manifest.json'),fixed_control_run=str(FIXED_RUN)))
    print(run)


def submit(run):
    import shlex
    command=['/home/4473331/.conda/envs/scaling_sandbox/bin/python',str(run/'frozen_sources/moving_pilot.py'),
             'prepare','--run',str(run)]
    wrapped='env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '+shlex.join(command)
    args=['sbatch','--parsable','--qos=small','--cpus-per-task=1','--mem=3G','--time=00:30:00',
        '--job-name=deB_moving_break','--output='+str(run/'logs/prepare_%j.out'),'--wrap',wrapped]
    result=subprocess.run(args,text=True,capture_output=True,check=True)
    job=result.stdout.strip().split(';')[0];write(run/'prepare_job.json',dict(job_id=job,command=args));print(job)


def prepare(run):
    REF.Engine=MovingEngine
    REF.prepare(run)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['initialize','submit','prepare'])
    parser.add_argument('--run',type=Path,default=RUN)
    args=parser.parse_args();run=args.run.resolve()
    if args.mode=='initialize':initialize(run)
    elif args.mode=='submit':submit(run)
    else:prepare(run)
