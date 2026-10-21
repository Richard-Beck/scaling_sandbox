"""Additive literal Figure 4 A3 coordinates; accepted analyses stay immutable.

The primary coordinate is unwrapped physical heading times actual C/(2*pi),
anchored at the native material origin. The initial branch in every qualified
section is [0,2*pi). Circumference changes therefore contribute to this
coordinate. Existing dilation-removed and track-reanchored estimates remain
separate sensitivities. This module does not estimate persistence or H1--H3.
Primary tracks split at unresolved antipodal edges without removing frames;
the two tied route choices are separate deterministic sensitivities.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import signature_figure4_analysis_20261020 as A
import signature_figure4_e0_20261020 as E

METHODS=('heading_literal_actualC','profile_peak_literal_actualC')
POLICIES=('split','plus_pi','minus_pi')
ANTIPODAL_ATOL=1e-12


def sections_with_cuts(valid,cuts):
    for start,stop in A.sections(valid):
        edges=[start]+[int(i) for i in np.flatnonzero(cuts) if start<i<stop]+[stop]
        yield from zip(edges[:-1],edges[1:])


def route_data(signal,C,valid,policy='split'):
    signal=np.asarray(signal,complex);C=np.asarray(C,float);valid=np.asarray(valid,bool)
    if signal.shape!=C.shape or signal.shape!=valid.shape or np.any(C<=0) or not np.isfinite(C).all() or not np.isfinite(signal[valid]).all() or np.any(abs(signal[valid])==0):raise ValueError('Matching finite signals and positive actual circumference required')
    if policy not in POLICIES:raise ValueError('An explicit antipodal route policy is required')
    phase=np.angle(signal);steps=(np.diff(phase)+np.pi)%(2*np.pi)-np.pi
    ties=np.r_[False,valid[1:]&valid[:-1]&(abs(abs(steps)-np.pi)<=ANTIPODAL_ATOL)]
    cuts=ties if policy=='split' else np.zeros(len(valid),bool)
    if policy!='split':steps=np.where(ties[1:],np.pi if policy=='plus_pi' else -np.pi,steps)
    values=np.full(len(C),np.nan)
    for start,stop in sections_with_cuts(valid,cuts):
        theta=np.r_[phase[start]%(2*np.pi),phase[start]%(2*np.pi)+np.cumsum(steps[start:stop-1])]
        values[start:stop]=theta*C[start:stop]/(2*np.pi)
    return values,cuts,ties


def literal_coordinate(signal,C,valid,policy='split'):
    return route_data(signal,C,valid,policy)[0]


def split_msd(values,valid,cuts,max_lag):
    nums=np.zeros(max_lag);counts=np.zeros(max_lag,dtype=np.int64)
    for start,stop in sections_with_cuts(valid,cuts):
        n=min(max_lag,stop-start-1)
        if n<1:continue
        # Center the completed x coordinate inside A.section_msd, preserving
        # the literal actual-C dilation term and both frames next to a cut.
        sn,sc,_=A.section_msd(values[start:stop],np.ones(stop-start,bool),n)
        nums[:n]+=sn;counts[:n]+=sc
    return nums,counts,np.maximum(len(values)-np.arange(1,max_lag+1),0)


def pair_counts(valid,cuts,max_lag):
    sizes=[stop-start for start,stop in sections_with_cuts(valid,cuts)]
    hist=np.bincount(sizes,minlength=len(valid)+1);indices=np.arange(len(hist))
    tail_n=np.cumsum(hist[::-1])[::-1];tail_size=np.cumsum((hist*indices)[::-1])[::-1];lag=np.arange(1,max_lag+1)
    return tail_size[lag+1]-lag*tail_n[lag+1]


def native_tracks(row):
    for key in ('scalar','profile','basal_profile','source_status'):E.require_hash(row[key+'_path'],row[key+'_sha256'])
    with np.load(row['scalar_path'],allow_pickle=False) as z:
        t=z['time'].astype(float);P=z['P_real']+1j*z['P_imag'];M=z['M1_real']+1j*z['M1_imag']
        C=np.broadcast_to(z['circumference'],t.shape).copy();contrast=z['contrast'].copy();regions=z['current_half_contrast_regions'].copy()
    if not np.allclose(P/P[0],M/M[0],rtol=1e-10,atol=1e-10):raise ValueError('Native physical moment and normalized heading disagree')
    if not np.allclose(np.diff(t),row['scalar_cadence'],rtol=1e-9,atol=1e-8) or not np.isclose(t[-1],row['native_horizon']):raise ValueError('Incomplete scalar observation schedule')
    with np.load(row['profile_path'],allow_pickle=False) as z:
        pt=z['time'].astype(float);pc=np.broadcast_to(z['circumference'],pt.shape).copy();primary=z['primary']
        q=z['material_positions'] if 'material_positions' in z else z['x']/pc[0]
        peak=np.mod(q[np.argmax(primary,axis=1)],1.);pcontrast=np.ptp(primary,axis=1)
    pi=np.rint((pt-t[0])/row['scalar_cadence']).astype(int)
    if not np.allclose(t[pi],pt,rtol=0,atol=1e-7) or not np.allclose(np.diff(pt),row['profile_cadence'],rtol=1e-9,atol=1e-8):raise ValueError('Peak profiles must align exactly with native scalar observations')
    valid=(abs(P)>=.1)&(regions==1)&(contrast>=.1*row['basal_contrast'])
    pv=(abs(P[pi])>=.1)&(regions[pi]==1)&(pcontrast>=.1*row['basal_contrast'])
    return {'heading_literal_actualC':(t,M,C,valid),'profile_peak_literal_actualC':(pt,np.exp(2j*np.pi*peak),pc,pv)}


def checked(root):
    m=E.read_json(root/'manifest.json')
    for p,digest in m['source_sha256'].items():E.require_hash(p,digest)
    for p,digest in m['input_receipts'].items():E.require_hash(p,digest)
    if str(Path(__file__).resolve()) not in m['source_sha256']:raise ValueError('Execute the actual frozen literal-site source')
    return m


def accepted_children(root):
    status_path=root/'results/analysis/status.json';status=E.read_json(status_path)
    if not status.get('passed') or not status.get('complete'):raise ValueError('Only a completed accepted stage may supply literal-site inputs')
    receipts={str(status_path):E.digest(status_path)}
    for name,digest in status['outputs'].items():E.require_hash(status_path.parent/name,digest)
    if (root/'results/arms.json').exists():
        arms_path=root/'results/arms.json';arms=E.read_json(arms_path)
        if not arms.get('verified'):raise ValueError('Stage native identity verification is required')
        receipts[str(arms_path)]=E.digest(arms_path)
        children=[(Path(r['path']),r['sigma'],r['analysis_manifest_sha256']) for r in arms['arms']]
        if len(children)!=len(status['arms']):raise ValueError('Every declared sigma arm is required')
        for (child,sigma,_),row in zip(children,status['arms']):
            if sigma!=row['sigma']:raise ValueError('Stage arm identity changed')
            E.require_hash(child/'results/analysis/status.json',row['status_sha256'])
            E.require_hash(child/'results/analysis/independent_estimand_audit.json',row['audit_sha256'])
    else:children=[(root,1.,E.digest(root/'analysis_manifest.json'))]
    groups=[]
    for child,sigma,manifest_sha in children:
        p=child/'analysis_manifest.json';E.require_hash(p,manifest_sha);m=E.read_json(p)
        audit=child/'results/analysis/independent_estimand_audit.json';a=E.read_json(audit)
        if not a.get('passed') or not a.get('complete') or a['analysis_manifest_sha256']!=manifest_sha or a['condition_groups']!=len(m['groups']):raise ValueError('The complete independent parent analysis audit is required')
        if a['inventory_sha256']!=m['inventory_sha256']:raise ValueError('The independent parent audit belongs to a different inventory')
        E.require_hash(child/'inventory.json',m['inventory_sha256'])
        receipts.update({str(p):manifest_sha,str(audit):E.digest(audit)})
        for source,digest in m['source_sha256'].items():E.require_hash(source,digest)
        for g in m['groups']:
            identity=hashlib.sha256(json.dumps(sorted((r['scalar_sha256'],r['profile_sha256'],r['basal_profile_sha256']) for r in g['noise'])).encode()).hexdigest()
            groups.append(dict(g,source_group_task=g['task'],source_root=str(root),source_stage=status.get('stage','E0'),sigma=sigma,native_group_sha256=identity))
    return groups,receipts


def prepare(root,parent,rule):
    if root.exists():raise FileExistsError('Use a fresh additive literal-site run root')
    definition=E.read_json(rule)
    if definition.get('root_approved') is not True or definition.get('primary',{}).get('method')!=METHODS[0] or definition.get('antipodal',{}).get('primary_policy')!='split':raise ValueError('The reviewed literal A3 antipodal-splitting definition is required')
    groups,receipts=accepted_children(parent);receipts[str(rule)]=E.digest(rule);root.mkdir(parents=True)
    for name in ('frozen_sources','workers','audits','logs','results'): (root/name).mkdir()
    sources=Path(__file__).parent
    for path in (Path(__file__),Path(A.__file__),Path(E.__file__),sources/'signature_figure4_inventory_20261020.py',sources/'signature_figure4_site_audit_20261020.py'):shutil.copy2(path,root/'frozen_sources'/path.name)
    for i,g in enumerate(groups):g['task']=i
    E.write_json(root/'manifest.json',dict(created_utc=E.stamp(),groups=groups,parent=str(parent),input_receipts=receipts,source_sha256={str(p):E.digest(p) for p in (root/'frozen_sources').glob('*.py')},
        bootstrap_resamples=2000,views=E.VIEWS,methods=METHODS,antipodal_policies=POLICIES,antipodal_atol_radians=ANTIPODAL_ATOL,source_version=2,scope='Additive literal A3 only; no native integration or H1-H3 change',coordinate='arg(M1) at material q=0; start each qualified section in[0,2pi), unwrap inside section, x=theta*C(t)/(2pi); primary splits labels at unresolved antipodal edges, preserving every native frame',
        qualification='R=abs(P)>=.1, exactly one native raw half-contrast component, contrast>=.1 exact native basal contrast; profile peak uses its own profile contrast'))
    return dict(prepared=True,groups=len(groups),manifest_sha256=E.digest(root/'manifest.json'))


def worker(root,task):
    m=checked(root);g=m['groups'][task];folder=root/f'workers/{task:04d}'
    if folder.exists():raise FileExistsError('Completed literal-site workers are immutable')
    folder.mkdir();per={v:{method:{policy:[] for policy in POLICIES} for method in METHODS} for v in E.VIEWS};cases=[];ambiguities=[]
    for row in g['noise']:
        tracks=native_tracks(row);cases.append(dict(id=row['id'],scalar_sha256=row['scalar_sha256'],profile_sha256=row['profile_sha256']))
        for view in E.VIEWS:
            for method,(t,signal,C,valid) in tracks.items():
                ix=E.select(t,view,row['native_horizon']);mask=valid[ix];lag=(len(ix)-1)//2
                _,cuts,ties=route_data(signal[ix],C[ix],mask,'split');base=pair_counts(mask,np.zeros(len(ix),bool),lag);unambiguous=pair_counts(mask,cuts,lag)
                ambiguities.append(dict(id=row['id'],view=view,method=method,cadence=float(t[ix[1]]-t[ix[0]]),qualified_frames=int(mask.sum()),total_frames=len(ix),qualified_adjacent_steps=int(np.sum(mask[1:]&mask[:-1])),antipodal_steps=int(ties.sum()),tracks_before_split=len(A.sections(mask)),tracks_after_split=sum(1 for _ in sections_with_cuts(mask,cuts))))
                for policy in POLICIES:
                    coord,route_cuts,_=route_data(signal[ix],C[ix],mask,policy);nums,counts,possible=split_msd(coord,mask,route_cuts,lag)
                    per[view][method][policy].append((nums,counts,possible,float(t[ix[1]]-t[ix[0]]),int(mask.sum()),len(ix),int(ties.sum()),base,unambiguous))
    key=g.get('bootstrap_pairing_key',f"{g['model']}:{g['C_label']}:{g['f']}")
    seed=2026102004+int(hashlib.sha256(key.encode()).hexdigest()[:8],16);w=A.bootstrap_weights(len(cases),2000,seed)
    arrays={};qualifications=[]
    for view in E.VIEWS:
        for method in METHODS:
            for policy in POLICIES:
                rows=per[view][method][policy];nums=np.array([r[0] for r in rows]);counts=np.array([r[1] for r in rows]);possible=np.sum([r[2] for r in rows],axis=0)
                base=np.sum([r[7] for r in rows],axis=0);unambiguous=np.sum([r[8] for r in rows],axis=0)
                msd=A.ratio(nums.sum(axis=0),counts.sum(axis=0));boot=A.ratio(w@nums,w@counts);ci=A.finite_quantiles(boot);prefix=view+'__'+method+'__'+policy+'__'
                for name,value in dict(lag=np.arange(1,len(msd)+1)*rows[0][3],MSD=msd,MSD_ci_lower=ci[0],MSD_ci_upper=ci[1],coverage=A.ratio(counts.sum(axis=0),possible),eligible_pairs=counts.sum(axis=0),possible_pairs=possible,
                        eligible_pairs_without_ambiguity_split=base,coverage_without_ambiguity_split=A.ratio(base,possible),unambiguous_pairs=unambiguous,pairs_lost_to_ambiguity_split=base-unambiguous,coverage_lost_to_ambiguity_split=A.ratio(base-unambiguous,possible),route_assumed_pairs=counts.sum(axis=0)-unambiguous).items():arrays[prefix+name]=value
                qualifications.append(dict(view=view,method=method,antipodal_policy=policy,cadence=rows[0][3],qualified_frames=sum(r[4] for r in rows),total_frames=sum(r[5] for r in rows),antipodal_steps=sum(r[6] for r in rows),seeds_with_antipodal_steps=sum(r[6]>0 for r in rows)))
    np.savez_compressed(folder/'arrays.npz',**arrays)
    E.write_json(folder/'result.json',dict(generated_utc=E.stamp(),passed=True,task=task,manifest_sha256=E.digest(root/'manifest.json'),arrays_sha256=E.digest(folder/'arrays.npz'),whole_seeds=len(cases),bootstrap_seed=seed,bootstrap_resamples=2000,cases=cases,qualifications=qualifications,case_ambiguities=ambiguities,all_qualified_native_frames_preserved=True))
    return dict(passed=True,task=task)


def collect(root):
    m=checked(root);out=root/'results/site_curves.csv.gz'
    if out.exists():raise FileExistsError('Do not overwrite a completed literal-site export')
    fields=['source_stage','source_root','model','C_label','f','sigma','native_horizon','native_group_sha256','view','method','antipodal_policy','cadence','n_seeds','bootstrap_resamples','antipodal_steps','seeds_with_antipodal_steps','qualified_frames','lag','MSD','MSD_ci_lower','MSD_ci_upper','coverage','eligible_pairs','possible_pairs','eligible_pairs_without_ambiguity_split','coverage_without_ambiguity_split','unambiguous_pairs','pairs_lost_to_ambiguity_split','coverage_lost_to_ambiguity_split','route_assumed_pairs','primary_A3','coordinate','profile_cadence_limited']
    receipts=[];count=0;case_rows=[]
    with out.open('wb') as binary,gzip.GzipFile(filename='',mode='wb',fileobj=binary,mtime=0) as compressed:
        import io
        with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as text:
            writer=csv.DictWriter(text,fieldnames=fields);writer.writeheader()
            for g in m['groups']:
                folder=root/f"workers/{g['task']:04d}";result=E.read_json(folder/'result.json');audit=root/f"audits/{g['task']:04d}.json";a=E.read_json(audit)
                if not result['passed'] or result['manifest_sha256']!=E.digest(root/'manifest.json') or not a.get('passed') or a['result_sha256']!=E.digest(folder/'result.json'):raise ValueError('All literal-site worker and independent audits must pass')
                E.require_hash(folder/'arrays.npz',result['arrays_sha256'])
                receipts.append(dict(task=g['task'],result_sha256=E.digest(folder/'result.json'),audit_sha256=E.digest(audit)))
                qualifiers={(r['view'],r['method'],r['antipodal_policy']):r for r in result['qualifications']}
                case_rows.extend(dict(r,source_stage=g['source_stage'],source_root=g['source_root'],model=g['model'],C_label=g['C_label'],f=g['f'],sigma=g['sigma'],native_horizon=g['noise'][0]['native_horizon'],native_group_sha256=g['native_group_sha256']) for r in result['case_ambiguities'])
                with np.load(folder/'arrays.npz') as z:
                    for view in E.VIEWS:
                        for method in METHODS:
                            for policy in POLICIES:
                                prefix=view+'__'+method+'__'+policy+'__';names=('lag','MSD','MSD_ci_lower','MSD_ci_upper','coverage','eligible_pairs','possible_pairs','eligible_pairs_without_ambiguity_split','coverage_without_ambiguity_split','unambiguous_pairs','pairs_lost_to_ambiguity_split','coverage_lost_to_ambiguity_split','route_assumed_pairs');values={k:z[prefix+k] for k in names}
                                qualification=qualifiers[(view,method,policy)]
                                context={k:g[k] for k in ('source_stage','source_root','model','C_label','f','sigma','native_group_sha256')};context.update(native_horizon=g['noise'][0]['native_horizon'],view=view,method=method,antipodal_policy=policy,cadence=float(values['lag'][0]),n_seeds=len(g['noise']),bootstrap_resamples=2000,
                                    antipodal_steps=qualification['antipodal_steps'],seeds_with_antipodal_steps=qualification['seeds_with_antipodal_steps'],qualified_frames=qualification['qualified_frames'],
                                    primary_A3=method==METHODS[0] and policy=='split',coordinate='fixed native material origin; section initial angle in[0,2pi); literal theta*C(t)/(2pi)',profile_cadence_limited=method==METHODS[1])
                                for j in range(len(values['lag'])):
                                    writer.writerow(dict(context,**{k:format(float(values[k][j]),'.17g') for k in names}));count+=1
    ambiguity=root/'results/case_ambiguities.csv.gz'
    with ambiguity.open('wb') as binary,gzip.GzipFile(filename='',mode='wb',fileobj=binary,mtime=0) as compressed:
        import io
        with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as text:
            writer=csv.DictWriter(text,fieldnames=sorted({k for r in case_rows for k in r}));writer.writeheader();writer.writerows(case_rows)
    E.write_json(root/'results/status.json',dict(generated_utc=E.stamp(),passed=True,complete=True,source_version=2,manifest_sha256=E.digest(root/'manifest.json'),condition_groups=len(m['groups']),rows=count,methods=METHODS,antipodal_policies=POLICIES,antipodal_atol_radians=ANTIPODAL_ATOL,
        outputs={p.name:E.digest(p) for p in (out,ambiguity)},receipts=receipts,all_qualified_native_frames_preserved=True,route_choices_are_not_independent_seeds=True,no_native_integrations=True,no_H1_H2_H3_estimand_change=True,Dx_interpretation='Existing Dtheta*(mean actual C/(2pi))^2 is a separate conversion, not a slope fit to literal coordinate MSD'))
    return dict(passed=True,complete=True,groups=len(m['groups']),rows=count)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=('prepare','worker','collect'));p.add_argument('--out',type=Path,required=True);p.add_argument('--parent',type=Path);p.add_argument('--rule',type=Path);p.add_argument('--task',type=int);args=p.parse_args();root=args.out.resolve()
    result=prepare(root,args.parent.resolve(),args.rule.resolve()) if args.command=='prepare' else worker(root,args.task) if args.command=='worker' else collect(root)
    print(json.dumps(result,indent=2))
