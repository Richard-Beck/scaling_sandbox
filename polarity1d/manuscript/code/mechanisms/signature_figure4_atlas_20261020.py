"""Frozen native-field S11 atlas with declared medoid/random selection.

No native integration, spatial resampling or main-example selection. Spring
mesh-face interpolation changes boundary geometry only, never field values.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import argparse
import hashlib
import html
import json
import os
import resource
import shlex
import shutil
import time
from urllib.parse import quote
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/'polarity1d/manuscript'
E0=ROOT/'data/runs/signature_figure4_e0_corrected_basals_20261020'
DEFAULT=ROOT/'data/runs/signature_figure4_native_atlas_E0_20261020'
RULES=REPO/'docs/tasks/briefs/20261020_figure4_s11_decisions.json'
FEATURES=['mean_relative_contrast','mean_basal_normalized_P','P_coverage_01','single_fraction','multiple_fraction','low_contrast_fraction','mean_qualified_count']
NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring='Spring',spring_ring='Spring')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def physical_faces(time,circumference,cells):
    t,C=np.asarray(time,float),np.asarray(circumference,float)
    if t.ndim!=1 or len(t)<2 or C.shape!=t.shape or np.any(np.diff(t)<=0) or np.any(C<=0):
        raise ValueError('Increasing native times and aligned positive actual circumferences required')
    # Include native times as faces too, so a circumference extremum at a
    # saved frame is not erased by midpoint averaging of the domain edge.
    edges=np.sort(np.r_[t,.5*(t[:-1]+t[1:])])
    rows=np.repeat(np.arange(len(t)),2)[1:-1]
    # Geometry only: repeated strips retain the exact nearest saved field.
    Cedges=np.interp(edges,t,C)
    x=Cedges[:,None]*np.linspace(0.,1.,cells+1)[None,:]
    return x,np.broadcast_to(edges[:,None],x.shape),rows


def selection_features(time,contrast,P,regions,basal_contrast,horizon):
    t,c,p,n=[np.asarray(v) for v in [time,contrast,P,regions]]
    if not (t.shape==c.shape==p.shape==n.shape) or np.any(np.diff(t)<=0) or abs(t[-1]-horizon)>1e-8:
        raise ValueError('Complete aligned native scalar histories required')
    if not np.isfinite(c).all() or not np.isfinite(p).all() or basal_contrast<=0:
        raise ValueError('Finite scalar fields and exact basal contrast required')
    starts=np.flatnonzero(abs(t-.5*horizon)<1e-8)
    if len(starts)!=1:raise ValueError('Half-horizon must be an exact retained scalar sample')
    k=int(starts[0]);t,c,p,n=t[k:],c[k:],abs(p[k:]),n[k:]
    dt=np.diff(t);duration=t[-1]-t[0]
    weighted=lambda value:float(value[:-1]@dt/duration)
    mean=lambda value:float((.5*(value[:-1]+value[1:]))@dt/duration)
    qualified=c>=.1*basal_contrast
    result=dict(mean_relative_contrast=mean(c/basal_contrast),mean_basal_normalized_P=mean(p),
        P_coverage_01=weighted(p>=.1),single_fraction=weighted(qualified&(n==1)),
        multiple_fraction=weighted(qualified&(n>=2)),low_contrast_fraction=weighted(~qualified),
        mean_qualified_count=weighted(np.where(qualified,n,0)),
        seed_forced_fraction=weighted((p<.1)|~qualified))
    return result


def choose(features,case_ids,condition,master_seed):
    matrix=np.asarray(features,float)
    if matrix.shape!=(len(case_ids),len(FEATURES)) or not np.isfinite(matrix).all() or len(set(case_ids))!=len(case_ids):
        raise ValueError('Unique complete independent cases and finite declared features required')
    iqr=np.quantile(matrix,.75,axis=0)-np.quantile(matrix,.25,axis=0)
    extent=np.ptp(matrix,axis=0);scale=np.where(iqr>0,iqr,extent)
    # Subtract physical features first: narrow IQRs must not amplify a
    # subtraction of large, separately normalized means.
    delta=abs(matrix[:,None,:]-matrix[None,:,:])
    distance=np.divide(delta,scale,out=np.zeros_like(delta),where=scale>0).sum(axis=2)
    scores=distance.sum(axis=1)
    best=float(scores.min());tied=np.flatnonzero(abs(scores-best)<=1e-12*max(1.,float(abs(scores).max())))
    medoid=min(tied,key=lambda i:case_ids[i]).item()
    random_keys=[hashlib.sha256((str(master_seed)+'|'+condition+'|'+identifier).encode()).hexdigest() for identifier in case_ids]
    remaining=[i for i in range(len(case_ids)) if i!=medoid]
    random=min(remaining,key=lambda i:random_keys[i]) if remaining else None
    return dict(medoid=medoid,random=random,scales=scale.tolist(),distance=distance.tolist(),scores=scores.tolist(),random_keys=random_keys)


def freeze(run,inventory,verification,condition_metrics,rules,stage,memory_gib=4,qos='small'):
    if run.exists():raise FileExistsError('Atlas uses a fresh immutable runroot')
    source,proof=read(inventory),read(verification)
    assert proof['passed'] and proof['inventory_sha256']==sha(inventory)
    assert proof['cases']==len(source['cases'])
    for name in ['inputs','workers','results','logs']:(run/name).mkdir(parents=True)
    shutil.copy2(__file__,run/'inputs/frozen_atlas.py');shutil.copy2(rules,run/'inputs/selection_rules.json')
    shutil.copy2(inventory,run/'inputs/canonical_inventory.json');shutil.copy2(verification,run/'inputs/native_identity_verification.json')
    conditions=pd.read_csv(condition_metrics)
    conditions=conditions[conditions.view=='main']
    entries=[]
    for case in source['cases']:
        case=dict(case);case['experiment']=case.get('experiment',case.get('origin',stage))
        case['atlas_context']=stage
        if stage!='E0' and 'independent_realization_id' not in case:
            raise ValueError('Expanded canonical inventories must explicitly declare independent-realization/reuse IDs')
        case['independent_realization_id']=case.get('independent_realization_id',case['id'])
        case['sigma']=float(case.get('sigma',1. if case['kind']=='noise' else 0.))
        if case['kind']=='noise':
            matches=conditions[(conditions.model==case['model'])&np.isclose(conditions.C_label,case['C_label'],rtol=1e-12,atol=1e-12)&np.isclose(conditions.f,case['f'],rtol=1e-10,atol=0.)]
            if 'sigma' in matches and matches.sigma.nunique()>1:
                matches=matches[np.isclose(matches.sigma,case['sigma'])]
            if len(matches)!=1:raise ValueError('Each noisy native group requires one canonical forced-state record')
            value=matches.iloc[0].get('forced_combined',matches.iloc[0].get('forced'))
            if pd.isna(value):case['condition_forced']=None
            elif isinstance(value,(bool,np.bool_)):case['condition_forced']=bool(value)
            else:raise ValueError('Forcing classification must be explicit True/False or unavailable')
        else:case['condition_forced']=None
        entries.append(case)
    # Densest verified observation of the same declared realization wins
    # within an arm. Other arms remain separate provenance, not extra seeds.
    preferred={};duplicates=[]
    for case in entries:
        identity=(case['model'],case['C_label'],case['f'],case['kind'],case['sigma'],stage,case['native_horizon'],case['independent_realization_id'])
        if identity not in preferred:preferred[identity]=case;continue
        old=preferred[identity]
        winner=min([old,case],key=lambda c:(c['profile_cadence'],not c.get('canonical_preferred',False),c['id']))
        preferred[identity]=winner;duplicates.append(dict(retained=winner['id'],observation_companions=[old['id'],case['id']],independent_realization_id=case['independent_realization_id']))
    groups=defaultdict(list)
    fields=['model','C_label','f','kind','sigma','atlas_context','native_horizon','profile_cadence']
    for case in preferred.values():groups[tuple(case[k] for k in fields)].append(case)
    jobs=[dict(slot=slot,condition=dict(zip(fields,key),experiment=stage),cases=sorted(cases,key=lambda c:c['id'])) for slot,(key,cases) in enumerate(sorted(groups.items()))]
    maximum_quads=max((2*int(round(j['condition']['native_horizon']/j['condition']['profile_cadence'])))*max(c['native_cells'] for c in j['cases'])*min(2,len(j['cases'])) for j in jobs)
    estimated_memory_gib=1.+maximum_quads*200./1024**3
    if memory_gib<estimated_memory_gib:
        raise ValueError(f'Dense native atlas needs a resource pilot and at least {estimated_memory_gib:.1f} GiB estimated plot memory; declared {memory_gib} GiB')
    plan=dict(created_utc=datetime.now(timezone.utc).isoformat(),jobs=jobs,stage=stage,
        code_sha256=sha(run/'inputs/frozen_atlas.py'),rules_sha256=sha(run/'inputs/selection_rules.json'),
        inventory_sha256=sha(inventory),verification_sha256=sha(verification),condition_metrics_sha256=sha(condition_metrics),
        independent_native_cases=len(preferred),source_case_records=len(entries),deduplicated_observation_companions=duplicates,
        stochastic_realizations=sum(c['kind']=='noise' for c in preferred.values()),
        deterministic_controls=sum(c['kind']!='noise' for c in preferred.values()),
        candidate_count_scope='Physical runs; deterministic controls contribute zero independent stochastic seeds.',
        main_example_selection=False,native_integrations=0,
        resources=dict(memory_gib=memory_gib,qos=qos,largest_native_quadrilaterals=maximum_quads,estimated_plot_memory_gib=estimated_memory_gib,estimate_bytes_per_quadrilateral=200,requires_representative_resource_pilot=True))
    write(run/'manifest.json',plan)
    command='exec '+shlex.join(['env','PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
        '/home/4473331/.conda/envs/scaling_sandbox/bin/python',str(run/'inputs/frozen_atlas.py'),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments=['sbatch','--parsable',f'--qos={qos}','--ntasks=1','--cpus-per-task=1',f'--mem={memory_gib}G','--time=01:00:00',f'--array=0-{len(jobs)-1}',
        '--job-name=figure4-native-atlas','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),'--wrap',command]
    write(run/'slurm_payload.json',dict(sbatch_arguments=arguments,jobs=len(jobs),submitted=False))
    return dict(groups=len(jobs),independent_native_cases=len(preferred),noise_groups=sum(j['condition']['kind']=='noise' for j in jobs))


def load_native(case):
    for key in ['profile','scalar','basal_profile','source_status']:
        assert sha(Path(case[key+'_path']))==case[key+'_sha256']
    with np.load(case['profile_path']) as z:
        t=z['time'].copy();p=z['primary'].copy();C=np.broadcast_to(z['circumference'],t.shape).copy()
        q=z['material_positions'].copy() if 'material_positions' in z else z['x'].copy()/float(C[0])
    assert p.shape==(len(t),case['native_cells'])
    np.testing.assert_allclose(q,(np.arange(len(q))+.5)/len(q),rtol=0,atol=1e-12)
    assert abs(t[-1]-case['native_horizon'])<1e-8 and np.isfinite(p).all() and np.all(C>0)
    np.testing.assert_allclose(np.diff(t),case['profile_cadence'],rtol=0,atol=1e-8)
    with np.load(case['scalar_path']) as z:
        st=z['time'].copy();P=z['P_real']+1j*z['P_imag'];input_mean=z['local_input_mean'].copy()
    return dict(time=t,primary=p,C=C,q=q,scalar_time=st,P=P,input_mean=input_mean)


def plot(job,selected,run):
    cases=[job['cases'][selected['medoid']]]
    roles=['Summary-feature medoid']
    if selected['random'] is not None:cases.append(job['cases'][selected['random']]);roles.append('Deterministic distinct random')
    movies=[load_native(c) for c in cases]
    cols=len(cases);fig,axes=plt.subplots(3,cols,figsize=(14 if cols==2 else 8,12),squeeze=False,gridspec_kw=dict(height_ratios=[1.8,1.05,.65]))
    high=max(float(m['primary'].max()) for m in movies);snapshot_records=[]
    for column,(case,movie,role) in enumerate(zip(cases,movies,roles)):
        t,p,C,q=[movie[k] for k in ['time','primary','C','q']];ax=axes[0,column]
        xfaces,tfaces,rows=physical_faces(t,C,len(q))
        image=ax.pcolormesh(xfaces,tfaces,p[rows],shading='flat',vmin=0.,vmax=high,cmap='viridis',rasterized=True)
        ax.axhspan(0.,.5*t[-1],color='white',alpha=.20);ax.axhline(.5*t[-1],color='white',ls='--',lw=1.)
        ax.invert_yaxis();ax.set_xlabel('Physical arc from fixed material origin');ax.set_ylabel('Native time')
        forced=case['condition_forced'];tag=('forced condition' if forced else 'nonforced condition' if forced is False else 'forcing status unavailable') if case['kind']=='noise' else 'deterministic control'
        ax.set_title(f"{role} · seed {case['seed']}\n{tag} · profile cadence {case['profile_cadence']:g}",fontsize=11)
        fig.colorbar(image,ax=ax,label='Native active Rac' if case['model']=='holmes_model3' else 'Native active concentration')
        snapshots=[];ax=axes[1,column]
        for fraction,color in zip([0.,.5,.75,1.],['#aaaaaa','#3974a0','#c58428','#76518c']):
            index=int(np.argmin(abs(t-fraction*t[-1])))
            ax.plot(q*C[index],p[index],color=color,lw=1.,label=f't={t[index]:g}, C={C[index]:.5g}')
            snapshots.append(dict(requested_horizon_fraction=fraction,native_saved_time=float(t[index]),native_frame=index,actual_circumference=float(C[index])))
        ax.set_xlabel('Actual physical arc at each saved snapshot');ax.set_ylabel('Native active concentration');ax.legend(fontsize=8,ncol=2)
        ax=axes[2,column];st,P=movie['scalar_time'],movie['P'];ax.plot(st,abs(P),color='#267f75',lw=.7)
        ax.axhline(.1,color='#aaaaaa',ls=':',lw=.7);ax.axvspan(0.,.5*st[-1],color='#dddddd',alpha=.5)
        ax.set_xlabel('Native time');ax.set_ylabel('|P| (basal M₁ normalized)')
        noise=ax.twinx();noise.plot(st,movie['input_mean'],color='#956c43',lw=.4,alpha=.6);noise.set_ylabel('Native spatial mean input',fontsize=8)
        snapshot_records.append(dict(case_id=case['id'],role=role,snapshots=snapshots,native_frames=len(t),native_cells=len(q),
            minimum_actual_circumference=float(C.min()),maximum_actual_circumference=float(C.max())))
    condition=job['condition'];model=condition['model'];C=condition['C_label']
    fig.suptitle(f"S11 · {NAMES.get(model,model)} · size label {C:g} · f={condition['f']:g} · {condition['kind']}\n"
        f"context {condition['experiment']} · σ={condition['sigma']:g} · horizon={condition['native_horizon']:g} · "+
        (f"n={len(job['cases'])} independent noisy realizations" if condition['kind']=='noise' else f"{len(job['cases'])} deterministic "+('control' if len(job['cases'])==1 else 'controls')),fontsize=13)
    caption=('Native concentrations and every retained field frame are plotted without spatial resampling.\n'
        'First half shaded: excluded from summary selection.\n'
        'Spring mesh-face interpolation changes only boundary geometry; snapshots use actual saved C(t).\n'
        'P uses basal M₁ normalization, distinct from instantaneous coherence.\n'
        'Supplemental selection and controls; no main example is chosen.\n'
        'Native NPZ/status links remain local/HPC data.') if cols==1 else (
        'Native concentrations and every retained field frame are plotted without spatial resampling. First half shaded: excluded from summary selection.\n'
        'Spring mesh-face interpolation applies only to boundary geometry; snapshots use their actual saved C(t). P uses basal M₁ normalization, not instantaneous coherence.\n'
        'Selection and controls are supplemental; no main example is chosen. Native NPZ/status links remain local/HPC data.')
    fig.text(.5,.015,caption,ha='center',fontsize=7 if cols==1 else 8)
    fig.tight_layout(rect=(0,.075,1,.94))
    target=run/'results/export/png'/f"S11_{condition['experiment']}_{job['slot']:04d}.png";target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target,dpi=170);plt.close(fig)
    webp=run/'results/export/webp'/(target.stem+'.webp');webp.parent.mkdir(exist_ok=True)
    with Image.open(target) as im:
        reference=np.asarray(im.convert('RGBA'));im.save(webp,format='WEBP',lossless=True,method=6,exact=True)
    with Image.open(webp) as im:assert np.array_equal(reference,np.asarray(im.convert('RGBA')))
    return dict(png=str(target),png_sha256=sha(target),webp=str(webp),webp_sha256=sha(webp),pixel_equal=True,
        displayed=snapshot_records,native_fields_interpolated=False,physical_mesh_geometry_interpolation_only=True)


def worker(run,slot):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native atlas production uses SLURM')
    started=time.monotonic()
    plan=read(run/'manifest.json');assert sha(Path(__file__))==plan['code_sha256']
    rules=read(run/'inputs/selection_rules.json');assert sha(run/'inputs/selection_rules.json')==plan['rules_sha256']
    job=plan['jobs'][slot];features=[]
    for case in job['cases']:
        assert sha(Path(case['scalar_path']))==case['scalar_sha256']
        with np.load(case['scalar_path']) as z:
            row=selection_features(z['time'],z['contrast'],z['P_real']+1j*z['P_imag'],z['current_half_contrast_regions'],case['basal_contrast'],case['native_horizon'])
        features.append(row)
    selected=choose([[row[k] for k in FEATURES] for row in features],[case['id'] for case in job['cases']],json.dumps(job['condition'],sort_keys=True),rules['master_seed'])
    artifact=plot(job,selected,run)
    write(run/'workers'/f'{slot:04d}.json',dict(complete=True,condition=job['condition'],cases=job['cases'],features=features,
        feature_names=FEATURES,selection=selected,artifact=artifact,code_sha256=plan['code_sha256'],rules_sha256=plan['rules_sha256'],
        resource=dict(wall_seconds=time.monotonic()-started,peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            SLURM_JOB_ID=os.environ['SLURM_JOB_ID'],SLURM_ARRAY_TASK_ID=os.environ.get('SLURM_ARRAY_TASK_ID'))))


def collect(run):
    plan=read(run/'manifest.json');records=[];feature_rows=[]
    assert sha(Path(__file__))==plan['code_sha256'] and sha(run/'inputs/selection_rules.json')==plan['rules_sha256']
    for job in plan['jobs']:
        record=read(run/'workers'/f"{job['slot']:04d}.json");assert record['complete'] and record['condition']==job['condition']
        a=record['artifact'];assert a['pixel_equal'] and sha(Path(a['webp']))==a['webp_sha256'] and sha(Path(a['png']))==a['png_sha256']
        with Image.open(a['png']) as original,Image.open(a['webp']) as portable:
            assert np.array_equal(np.asarray(original.convert('RGBA')),np.asarray(portable.convert('RGBA')))
        selected=record['selection'];matrix=[[row[k] for k in FEATURES] for row in record['features']]
        recomputed=choose(matrix,[case['id'] for case in job['cases']],json.dumps(job['condition'],sort_keys=True),read(run/'inputs/selection_rules.json')['master_seed'])
        assert recomputed==selected
        for i,(case,values) in enumerate(zip(job['cases'],record['features'])):
            role='medoid' if i==selected['medoid'] else 'random' if i==selected['random'] else 'not_selected'
            feature_rows.append(dict(**job['condition'],case_id=case['id'],independent_realization_id=case['independent_realization_id'],seed=case['seed'],role=role,
                source_observation_id=case.get('source_observation_id',case['id']),origin=case.get('origin',case['experiment']),
                physical_run_id=case.get('physical_run_id',case['independent_realization_id']),latent_Gaussian_cohort_id=case.get('latent_Gaussian_cohort_id'),
                medoid_distance_sum=selected['scores'][i],random_sha256=selected['random_keys'][i],**values))
        records.append(record)
    out=run/'results/export';pd.DataFrame(feature_rows).to_csv(out/'candidate_features_and_selection.csv',index=False)
    shutil.copy2(run/'inputs/selection_rules.json',out/'selection_rules.json')
    shutil.copy2(run/'inputs/frozen_atlas.py',out/'frozen_atlas.py')
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),complete=True,passed=True,groups=len(records),
        independent_native_cases=plan['independent_native_cases'],selected_displays=sum(len(r['artifact']['displayed']) for r in records),
        lossless_pixel_equal_groups=len(records),manifest_sha256=sha(run/'manifest.json'),selection_rules_sha256=plan['rules_sha256'],
        native_integrations=0,main_example_selection=False,stage=plan['stage'],sourcecode_sha256=plan['code_sha256'])
    summary.update(stochastic_realizations=plan.get('stochastic_realizations'),deterministic_controls=plan.get('deterministic_controls'),
        candidate_count_scope='Physical runs, with deterministic controls separate from independent stochastic realizations.')
    write(out/'native_atlas_index.json',dict(summary=summary,groups=records))
    css='body{font:17px/1.5 system-ui;margin:30px auto;max-width:1300px;padding:0 20px;color:#14232b}img{width:100%;border:1px solid #ccc}select,input{font-size:16px;padding:8px;margin:8px}code{font-size:12px;overflow-wrap:anywhere}'
    entries=[]
    rel=lambda path:quote(os.path.relpath(Path(path),out),safe='/._-')
    for record in records:
        c=record['condition'];a=record['artifact'];title=f"{c['experiment']} · {NAMES.get(c['model'],c['model'])} · size {c['C_label']:g} · f {c['f']:g} · {c['kind']} · σ {c['sigma']:g} · horizon {c['native_horizon']:g} · cadence {c['profile_cadence']:g}"
        files=[]
        for display in a['displayed']:
            case=next(case for case in record['cases'] if case['id']==display['case_id'])
            for kind in ['profile','scalar','basal_profile','source_status']:
                files.append(dict(label=display['role']+' · '+case['id']+' · '+kind+' (local/HPC data)',url=rel(case[kind+'_path']),sha256=case[kind+'_sha256']))
        entries.append(dict(title=title,image=rel(a['webp']),model=c['model'],kind=c['kind'],files=files))
    document='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure4 S11 native atlas</title><style>'+css+'</style><h1>S11 · Native activity atlas</h1><p>Generated '+html.escape(summary['generated_utc'])+'. All producers stopped; no timed refresh. '+str(summary['groups'])+' complete groups. Every displayed WebP is lossless and pixel-identical to its retained PNG.</p><p>All native active-field frames and concentrations are shown. Physical arcs use actual circumference, including Spring. The first half is shaded and excluded from summary-feature selection. Four snapshots use native saved samples nearest 0, 0.5, 0.75 and 1 horizon; actual times and circumferences are printed. P uses basal-first-moment normalization. Native NPZ/status links are retained local/HPC archives, not implied to be included in a compact checkout.</p><p>Each group has the declared summary-feature medoid and a deterministic distinct random case; singleton controls have one display. No main example is selected. Every σ, horizon, stage and cadence remains separate; declared observation companions do not add independent cases. <a href="selection_rules.json">Predeclared selection rules</a> · <a href="frozen_atlas.py">Frozen renderer and selection source</a> · <a href="candidate_features_and_selection.csv">Every candidate/selection score</a> · <a href="native_atlas_index.json">Complete source and geometry record</a></p><input id="search" placeholder="Model, size, f, kind or stage"><select id="group"></select><h2 id="title"></h2><img id="image" alt="Native physical-arc kymographs and saved profiles"><div id="sources"></div><script>const groups='+json.dumps(entries,separators=(',',':')).replace('</',r'<\/')+';const byId=id=>document.getElementById(id);let selected=[];function filter(){const q=byId("search").value.toLowerCase();selected=groups.filter(g=>g.title.toLowerCase().includes(q));byId("group").replaceChildren(...selected.map((g,i)=>{const o=document.createElement("option");o.value=i;o.textContent=g.title;return o}));show()}function show(){const g=selected[byId("group").value];if(!g)return;byId("title").textContent=g.title;byId("image").src=g.image;const list=document.createElement("ul");for(const f of g.files){const li=document.createElement("li"),a=document.createElement("a"),code=document.createElement("code");a.href=f.url;a.textContent=f.label;code.textContent=" · SHA256 "+f.sha256;li.append(a,code);list.append(li)}byId("sources").replaceChildren(list)}byId("search").addEventListener("input",filter);byId("group").addEventListener("change",show);filter();</script></html>'
    (out/'index.html').write_text(document)
    write(out/'complete.json',dict(summary,outputs={name:sha(out/name) for name in ['index.html','native_atlas_index.json','candidate_features_and_selection.csv','selection_rules.json','frozen_atlas.py']}))
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','worker','collect']);p.add_argument('--run',type=Path,default=DEFAULT)
    p.add_argument('--inventory',type=Path,default=E0/'inventory.json');p.add_argument('--verification',type=Path,default=E0/'results/native_identity_verification.json')
    p.add_argument('--condition-metrics',type=Path,default=E0/'results/analysis/condition_metrics.csv');p.add_argument('--rules',type=Path,default=RULES)
    p.add_argument('--stage',default='E0');p.add_argument('--slot',type=int);p.add_argument('--memory-gib',type=int,default=4);p.add_argument('--qos',default='small');args=p.parse_args();run=args.run.resolve()
    result=freeze(run,args.inventory.resolve(),args.verification.resolve(),args.condition_metrics.resolve(),args.rules.resolve(),args.stage,args.memory_gib,args.qos) if args.action=='freeze' else worker(run,args.slot) if args.action=='worker' else collect(run)
    if result is not None:print(json.dumps(result))


if __name__=='__main__':main()
