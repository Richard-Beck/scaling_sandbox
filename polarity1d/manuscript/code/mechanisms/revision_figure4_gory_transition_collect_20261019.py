"""Read-only terminal collection of the separate26-case Gory refinement."""
from pathlib import Path
from datetime import datetime,timezone
import gzip,hashlib,importlib.util,json,os,sys
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript';RUN=MS/'data/runs/revision_gory_strength_transition_20261019'
REPL=MS/'data/runs/revision_noise_strength_replication_scoped_20261019'
FIG=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
OUT=FIG/'gory_transition_refinement_20261019'
COLLECTOR=MS/'data/runs/revision_noise_strength_broad_20261019/results/monitor/collector_v2.py'
COLLECTOR_SHA='79747ebe6d831571e302a1fec528b7c6f58b223d7c9d7a8ea827b22a134fdd83'
ASSEMBLY=FIG/'signature_assembly_20261019/supporting/metadata/final_frozen_sources/assembly.py'
ASSEMBLY_SHA='7c1ef9929f8fbbda798e02c97ee7db29d152c229082a5d7242767f72fd2283fe'
FEATURES=('mean_R','mean_A','mean_Q','mean_contrast_ratio','qualified_multiregion_fraction','weak_R_fraction','low_contrast_fraction')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    assert os.environ.get('SLURM_JOB_ID')
    for name in ('data','metadata','vectors'):(OUT/'supporting'/name).mkdir(parents=True,exist_ok=True)
    data=OUT/'supporting/data';meta=OUT/'supporting/metadata';m=json.loads((RUN/'manifest.json').read_text())
    assert sha(RUN/'manifest.json')=='739cd95964c9dcb8818676e82f844fff995fcaff3aec08059c1d742ec99281d0'
    assert sha(COLLECTOR)==COLLECTOR_SHA and sha(ASSEMBLY)==ASSEMBLY_SHA
    reader=load('accepted_gory_readouts',COLLECTOR);analysis=load('accepted_gory_joint_bootstrap',ASSEMBLY)
    old=json.loads((REPL/'results/partial/case_cache.json').read_text());cache={};inventory=[]
    for record in m['retained_context']:
        c=old[record['case_id']]
        for key in ('trace','profile'):assert sha(c[key])==record[key+'_sha256']
        cache[record['case_id']]=c
    for job in m['jobs']:
        directory=RUN/'workers'/f'{job["slot"]:04d}';p=directory/'status.json'
        status=json.loads(p.read_text()) if p.exists() else {}
        inventory.append(dict(slot=job['slot'],kind=job['kind'],fraction=job['fraction'],seed=job['seed'],complete=status.get('complete',False),status=status.get('status','missing'),elapsed_seconds=status.get('elapsed_seconds'),error=status.get('error')))
        if not status.get('complete'):continue
        assert status['job']==job and status['producer']=='terminal' and status['manifest_sha256']==sha(RUN/'manifest.json')
        for name,digest in status['output_sha256'].items():assert sha(directory/name)==digest
        base=dict(case_id=f'transition_{job["slot"]:04d}',model='goryachev',platform='chemical',size=60.,fraction=job['fraction'],kind=job['kind'],seed=job['seed'],legacy=False)
        c=reader.summarize_case(base,directory/'trace.npz',directory/'primary.npz')
        c.update(trace=str(directory/'trace.npz'),profile=str(directory/'primary.npz'),trace_sha256=sha(directory/'trace.npz'),profile_sha256=sha(directory/'primary.npz'))
        cache[base['case_id']]=reader.clean(c)
    pd.DataFrame(inventory).to_csv(data/'worker_inventory.csv',index=False)
    complete=all(r['complete'] for r in inventory)
    status=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='terminal' if complete else 'terminal_with_incomplete_outputs',new_complete=sum(r['complete'] for r in inventory),new_target=26,available_cases=len(cache),target_cases=53,native_integrations=0)
    write(meta/'COLLECTION_STATUS.json',status)
    if not complete:
        (OUT/'README.md').write_text('# Gory transition refinement\n\nTerminal collection found incomplete native outputs. See supporting/metadata/COLLECTION_STATUS.json and supporting/data/worker_inventory.csv. No resubmission or complete-result inference was made.\n');return
    assert len(cache)==53
    descriptors=pd.DataFrame([d for c in cache.values() for d in c['descriptors']]);memory=pd.DataFrame([d for c in cache.values() for d in c['memory']])
    descriptors.to_csv(data/'case_window_descriptors.csv',index=False);memory.to_csv(data/'case_window_memory.csv',index=False)
    keys=['model','size','fraction','kind','window'];condition=descriptors.groupby(keys,sort=True)[list(analysis.METRICS)].mean().reset_index();condition.to_csv(data/'condition_window_summaries.csv',index=False)
    groups=[]
    for key,g in memory.groupby(keys+['lag'],sort=True):
        n=g.pairs_eligible.sum();total=g.pairs_available.sum();support=n>=100 and n/total>=.01
        groups.append(dict(zip(keys+['lag'],key),direction_cosine_reported=g.direction_cosine_sum.sum()/n if support else np.nan,direction_coverage=n/total,raw_vector_correlation=g.raw_product_sum.sum()/np.sqrt(g.energy_a.sum()*g.energy_b.sum())))
    pd.DataFrame(groups).to_csv(data/'memory_all_conditions_windows_lags.csv',index=False)
    ci,cm,error=analysis.analyze(OUT);assert error<1e-12
    native=[];selection=[]
    for f in (.01,.015,.02,.03):
        candidates=[(key,c,next(d for d in c['descriptors'] if d['window']=='last_half')) for key,c in cache.items() if c['descriptors'][0]['kind']=='noise' and c['descriptors'][0]['fraction']==f]
        candidates.sort(key=lambda v:v[2]['seed']);assert len(candidates)==12
        X=np.array([[d[k] for k in FEATURES] for _,_,d in candidates]);span=np.ptp(X,axis=0);span[span<1e-12]=1
        scores=np.sum(np.sum(abs(X[:,None]-X[None,:])/span,axis=2),axis=1);chosen=np.flatnonzero(np.isclose(scores,scores.min(),rtol=0,atol=1e-12))[0]
        case_id,c,d=candidates[int(chosen)];mean,=[v for v in cache.values() if v['descriptors'][0]['kind']=='mean' and v['descriptors'][0]['fraction']==f]
        with np.load(c['profile']) as z:t=z['time'];u=z['primary'];x=z['x']
        with np.load(mean['profile']) as z:mu=z['primary']
        snapshot=int(np.flatnonzero(t==18000)[0]);take=(t>=16000)&(t<=20000)
        native.append(dict(f=f,t=t[take],u=u[take],x=x,basal=u[0],snapshot=u[snapshot],mean=mu[snapshot]))
        selection.append(dict(fraction=f,case_id=case_id,seed=d['seed'],snapshot_time=18000.,profile=c['profile'],profile_sha256=c['profile_sha256'],mean_profile=mean['profile'],mean_profile_sha256=mean['profile_sha256']))
        with gzip.open(data/f'native_medoid_f{f:g}.csv.gz','wt') as handle:
            np.savetxt(handle,np.column_stack((np.r_[0,t[take]],np.vstack((u[0],u[take])))),fmt='%.17g',delimiter=',',header='time,'+','.join(f'u_{j:04d}' for j in range(u.shape[1])),comments='')
        pd.DataFrame(dict(x=x,basal=u[0],noisy=u[snapshot],matched_mean=mu[snapshot])).to_csv(data/f'profiles_f{f:g}.csv',index=False)
    pd.DataFrame(selection).to_csv(data/'native_medoid_selection.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
    with PdfPages(OUT/'gory_transition_review.pdf') as pdf:
        fig,axes=plt.subplots(2,2,figsize=(10.5,8.4));fig.subplots_adjust(left=.1,right=.9,bottom=.15,top=.82,hspace=.4,wspace=.4)
        fig.suptitle('Goryachev C60: targeted refinement of the sampled transition',x=.08,ha='left',fontsize=16)
        fig.text(.08,.9,'f=.015/.02 are new; .01/.03 are exact retained anchors. Twelve noise seeds and a matched mean control per strength.',fontsize=9.5)
        d=ci[(ci.kind=='noise')&(ci.window=='last_half')].sort_values('fraction');means=ci[(ci.kind=='mean')&(ci.window=='last_half')].sort_values('fraction');m1=cm[(cm.kind=='noise')&(cm.window=='last_half')&(cm.lag==8000)].sort_values('fraction')
        for ax,metric,label,color,style in [(axes[0,0],'qualified_multiregion_fraction','Qualified multiple regions','#a13958','-'),(axes[0,0],'low_contrast_fraction','Low contrast','#a13958','--'),(axes[0,1],'mean_contrast_ratio','Noisy contrast/basal','#4055a8','-'),(axes[1,0],'mean_R','R','#4055a8','-'),(axes[1,0],'mean_A','Active amount A','#168b84','-'),(axes[1,0],'mean_Q','Coherence Q','#b87c27','-')]:
            ax.plot(d.fraction,d[metric],color=color,ls=style,marker='o',label=label);ax.vlines(d.fraction,d[metric+'_ci_low'],d[metric+'_ci_high'],color=color,lw=.8)
        axes[0,1].plot(means.fraction,means.mean_contrast_ratio,'s--',color='#168b84',label='Matched mean control')
        axes[1,1].plot(m1.fraction,m1.direction_cosine_reported,'o-',color='#4055a8',label='Conditional heading cosine');axes[1,1].vlines(m1.fraction,m1.direction_cosine_reported_ci_low,m1.direction_cosine_reported_ci_high,color='#4055a8',lw=.8)
        axes[1,1].plot(m1.fraction,m1.direction_coverage,'s--',color='#168b84',label='Usable-pair fraction')
        for ax in axes.flat:ax.set_xticks([.01,.015,.02,.03]);ax.set_xlabel('Input factor f');ax.spines[['top','right']].set_visible(False);ax.legend(fontsize=8,frameon=False);ax.grid(axis='y',alpha=.2)
        axes[0,0].set_ylim(-.03,1.03);axes[1,1].set_ylim(-1.03,1.03)
        axes[0,0].set_ylabel('Fraction of all scalar times');axes[0,1].set_ylabel('Contrast relative to own basal');axes[1,0].set_ylabel('Paired amount/coherence readouts');axes[1,1].set_ylabel('Cosine / fraction')
        fig.text(.08,.07,'Second-half window16000–32000; 2000 whole-seed pointwise intervals. Direction lag8000 requires endpoint R≥0.1.',fontsize=9)
        fig.text(.08,.045,'No phase boundary, competition, continuous lifetime or migration law inferred. Positive mean and fluctuations scale together.',fontsize=9)
        fig.savefig(OUT/'gory_transition_readouts.png',dpi=170);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(12.4,7.5),gridspec_kw={'height_ratios':[1.2,1.]});fig.subplots_adjust(left=.07,right=.93,bottom=.15,top=.8,wspace=.35,hspace=.5)
        fig.suptitle('Goryachev transition: native medoids and fixed-time profiles',x=.07,ha='left',fontsize=16)
        fig.text(.07,.89,'Same zero-based native activity scale; fixed window16000–20000 and snapshot18000. Medoids use the unchanged seven-descriptor rule.',fontsize=9)
        maximum=max(float(v['basal'].max()) for v in native)
        for col,v in enumerate(native):
            ax=axes[0,col];im=ax.imshow(v['u'].T,origin='lower',aspect='auto',extent=[16000,20000,0,60],vmin=0,vmax=maximum,cmap='viridis');ax.axvline(18000,color='white',ls='--',lw=.8);ax.set_title(f'f={v["f"]:g}',loc='left');ax.set_xlabel('Native time');ax.set_xticks([16000,20000])
            pr=axes[1,col];pr.plot(v['x'],v['basal'],color='#9aa2a6',label='Exact basal');pr.plot(v['x'],v['mean'],'--',color='#168b84',label='Matched mean');pr.plot(v['x'],v['snapshot'],color='#a13958',label='Noisy snapshot');pr.set_ylim(0,maximum*1.04);pr.set_xlabel('Physical arc');pr.set_xlim(0,60)
            zoom=pr.inset_axes([.36,.38,.59,.5]);zoom.plot(v['x'],v['mean'],'--',color='#168b84',lw=.8);zoom.plot(v['x'],v['snapshot'],color='#a13958',lw=.8);zoom.set_title('Native-scale zoom',fontsize=6.5);zoom.tick_params(labelsize=6)
            for item in (ax,pr):item.tick_params(labelsize=8);item.spines[['top','right']].set_visible(False)
        axes[0,0].set_ylabel('Physical arc');axes[1,0].set_ylabel('Native activity')
        fig.legend(*axes[1,0].get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.065,.86),ncol=3,fontsize=8,frameon=False)
        cbar=fig.colorbar(im,ax=list(axes[0]),fraction=.014,pad=.02);cbar.set_label('Native activity',fontsize=8)
        fig.text(.07,.07,'White lines mark the exact native profiles below. Insets enlarge native noisy/mean amplitudes, without normalization.',fontsize=9)
        fig.text(.07,.045,'The incumbent crosses the 0/60 periodic seam: its two plotted edge portions are one circular region.',fontsize=9)
        fig.savefig(OUT/'gory_transition_native.png',dpi=170);pdf.savefig(fig);plt.close(fig)
    receipt=dict(**status,source_sha256=sha(__file__),accepted_readout_source_sha256=COLLECTOR_SHA,accepted_bootstrap_source_sha256=ASSEMBLY_SHA,
        source_manifest_sha256=sha(RUN/'manifest.json'),unchanged936_main=True,maximum_point_reconciliation=error,independent_review='pending',
        outputs={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.suffix in ('.csv','.png','.pdf','.gz')})
    write(meta/'COLLECTION_RECEIPT.json',receipt)
    (OUT/'README.md').write_text(f'''# Goryachev C60 transition refinement\n\n[Two-page review](gory_transition_review.pdf) · [Readouts](gory_transition_readouts.png) · [Native histories/profiles](gory_transition_native.png) · [Source/collection receipt](supporting/metadata/COLLECTION_RECEIPT.json)\n\nData generated {status['generated_at']}; native producer terminal. All26 new cases completed:24 noisy paths and2 matched-mean controls. With exact retained .01/.03 and no-input contexts,53 cases are represented. The completed936-case main/supplement closure remains unchanged.\n\nNew f=.015/.02 refine the unsampled(.01,.03) interval; they do not exclude still narrower intervals. All native protocols and the qualified basal are identical. Summaries use explicit full/half windows and independent-seed uncertainty. Read occupancy together with contrast, R/A/Q, matched means and conditional heading coverage. Multiple regions alone do not establish competition or a phase. Native profiles are saved every40units; sampled lifetimes/establishment are not inferred. The local medoid/snapshot selection is fixed by the original descriptor/time rules, not a search for extra peaks.\n\nAll case and condition/CI tables, exact selection/source hashes, compact kymographs and actual profiles are under supporting/data. Native fields remain local for replay. No timed browser refresh is used.\n''')
    print(json.dumps({k:v for k,v in receipt.items() if k!='outputs'},indent=2))
if __name__=='__main__':main()
