"""Select by analytical score ONLY, freeze parameters, then check tracks."""
import html,json,base64
from datetime import datetime,timezone
import numpy as np
from fit import *

def table(rows,keys):
    def cell(v):return f'{v:.5g}' if isinstance(v,float) else html.escape(str(v))
    return '<table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in keys)+'</tr>'+''.join('<tr>'+''.join('<td>'+cell(r.get(k,''))+'</td>' for k in keys)+'</tr>' for r in rows)+'</table>'

def publish():
    from quality_figures import generate
    from parameter_figures import generate as generate_parameters
    generate()
    generate_parameters()
    fits=[json.loads(p.read_text()) for p in sorted(OUT.glob('fit_*.json'))]
    sections=[]
    for scenario in SCENARIOS:
        candidates=[f for f in fits if f['scenario']==scenario]
        if not candidates:continue
        best=min(candidates,key=lambda f:f['total_score']);vpath=OUT/f'validation_{scenario}.json'
        sections.append('<details><summary style="font-size:20px;margin-top:24px;cursor:pointer">'+SCENARIO_LABELS[scenario]+' — detailed tables</summary>')
        sections.append(f'<h2>{scenario}: {len(candidates)}/8 starts finished</h2><p>Best analytical objective {best["total_score"]:.5g}; data {best["data_score"]:.5g}; pooling {best["pooling_penalty"]:.5g}; persistence prior {best["persistence_prior_penalty"]:.5g}.</p>')
        sections.append(table([dict(parameter=k,**v) for k,v in best['hyperparameters'].items()],['parameter','population_baseline','population_wgd_fold','fixed_log_baseline_sd','fixed_log_effect_sd']))
        pars,drift,_=decode(best['z'])
        sections.append(table([dict(study=s,tau_minus=pars[i,0,0],tau_plus=pars[i,1,0],sigma_minus=pars[i,0,1],sigma_plus=pars[i,1,1],tau_fold=pars[i,1,0]/pars[i,0,0],sigma_fold=pars[i,1,1]/pars[i,0,1]) for i,s in enumerate(STUDIES)],['study','tau_minus','tau_plus','tau_fold','sigma_minus','sigma_plus','sigma_fold']))
        sections.append(f'<p>Hosny signed drift: reference {drift[0]:.5g}, WGD {drift[1]:.5g} µm/h.</p>')
        sections.append(table(best['rows'],['study','group','metric','quantile','obs','predicted','residual']))
        if vpath.exists():
            val=json.loads(vpath.read_text())
            sections.append('<h3>Frozen-parameter track checks</h3><p>16,384 tracks per condition in each of two independent scrambles. No optimization uses these results. Quartile index: 0=25th, 1=median, 2=75th.</p>')
            sections.append(table(val['comparison'],['study','group','metric','quantile','obs','analytic','tracks_seed1','tracks_seed2','role']))
        sections.append(table(candidates,['start','total_score','data_score','pooling_penalty','persistence_prior_penalty','success','nfev','optimality','seconds']))
        sections.append('</details>')
    now=datetime.now(timezone.utc).isoformat()
    doc='''<!doctype html><meta charset="utf-8"><title>Gaussian hierarchical WGD fit</title><style>body{font:16px system-ui;max-width:1400px;margin:30px auto;padding:20px;color:#203040}table{border-collapse:collapse;font-size:13px;margin:20px 0}td,th{padding:7px;border:1px solid #ccd;text-align:right}th{background:#edf2f8}h2{margin-top:40px}</style><h1>Gaussian hierarchical WGD fit</h1>'''
    doc+=f'<p>Generated {now}. Completed optimization starts: {len(fits)}/{8*len(SCENARIOS)}. Refresh manually to load updated files.</p>'
    doc+='''<p>Five studies; Hosny regular hypoxic N/P baseline only. All optimization predictions are analytical approximations. Each study has reference persistence and noise strength and its own WGD fold changes. Log baselines and log fold changes are shrunk toward separately fitted common means. Primary fixed log standard deviations are log(10) and log(2); strong/weak halve/double both. These are regularization choices, not inferred population variances or confidence intervals.</p><p>Persistence is in hours, bounded between 10⁻⁶ and 12 for both groups, with effects in either direction. Noise strength is in µm/h^(3/2). Hosny alone has a signed drift baseline and WGD change, unpooled because only one directed assay is available.</p><p>Mallin mean directness is approximated by mean endpoint/mean path. Hosny fits only median path, directness and signed directional migration: median path ≈ mean Rice step length × number of steps; median directness ≈ Rice endpoint median/mean path; median signed migration ≈ mean parallel displacement/mean path. Path variability and numerator–denominator dependence are omitted in optimization. Hosny quartiles are reserved for checks. Other studies retain their analytical speed/MSD/path predictions.</p><p>Score = squared data residuals + squared pooling deviations + squared standardized log-persistence deviations from the prior center. Scraped ranges set log tolerances; Hosny residuals use observed IQRs. No likelihood or posterior precision is claimed. Eight independent starts per pooling choice use bounded trust-region least squares. Selection is based exclusively on the analytical objective; subsequent track checks never trigger refitting.</p>'''
    doc+='<p><b>Persistence prior:</b> the same soft log-scale penalty applies to all ten study/condition persistence times. Primary median 10 minutes, log SD log(10)/1.96 (untruncated central 95%: 1–100 minutes). Sensitivity medians 5 and 20 minutes retain this width, using primary pooling; strong/weak pooling retain the 10-minute prior. Existing numerical bounds remain. This is a Gaussian penalty in log-persistence coordinates, not a density-mode calculation in raw hours. No separate noise prior was added.</p>'
    completed=len(list(OUT.glob('validation_*.json')))
    doc+=f'<p>Track validation: {completed}/{len(SCENARIOS)} settings complete. '+('Campaign complete; producer stopped.' if completed==len(SCENARIOS) else 'Validation pending or in progress.')+'</p>'
    doc+='<h2>Parameters: group means and individual experiments</h2><p>Black diamonds show geometric group means, consistent with the fitted log-scale hierarchy. Colored points show each experiment; connecting lines pair WGD− and WGD+. Persistence is displayed in minutes. The WGD-effect plots show the corresponding fold changes. Hosny drift is unpooled and has no across-study mean. These are fitted point estimates, not uncertainty intervals.</p>'
    for scenario in SCENARIOS:
        if scenario!='primary':doc+=f'<details><summary style="font-size:20px;cursor:pointer">{SCENARIO_LABELS[scenario]} — parameter plots</summary>'
        for prefix in ['parameters','parameter_effects']:
            path=OUT/f'{prefix}_{scenario}.png'
            if path.exists():doc+=f'<img src="data:image/png;base64,{base64.b64encode(path.read_bytes()).decode()}" style="width:100%" alt="{prefix}: {SCENARIO_LABELS[scenario]}">'
        if scenario!='primary':doc+='</details>'
    doc+='<h2>Observed measurements versus predictions</h2><p>Start with the primary fit below. Dark circles show observations, orange diamonds analytical predictions, and blue squares the frozen-parameter track checks. Observed intervals are scraped reading ranges except in Hosny, where bars show the middle 50% of cells. They are not confidence intervals. Panels identify fitted measurements and held-out checks.</p>'
    for scenario in SCENARIOS:
        if scenario!='primary':doc+=f'<details><summary style="font-size:20px;cursor:pointer">{SCENARIO_LABELS[scenario]} — measurement plots</summary>'
        for study in STUDIES:
            name=f'quality_{scenario}_{study}'
            if (OUT/(name+'.png')).exists():
                doc+=f'<figure style="margin:24px 0"><img src="data:image/png;base64,{base64.b64encode((OUT/(name+".png")).read_bytes()).decode()}" style="width:100%;max-width:1400px" alt="{study}: observations and model predictions, {SCENARIO_LABELS[scenario]}"></figure>'
        if scenario!='primary':doc+='</details>'
    if (OUT/'tracks_primary.png').exists():doc+='<h2>Example tracks: primary fit</h2><img src="data:image/png;base64,'+base64.b64encode((OUT/"tracks_primary.png").read_bytes()).decode()+'" style="width:100%" alt="Example tracks for five studies and both WGD conditions">'
    doc+=''.join(sections)
    tmp=OUT/'report.tmp';tmp.write_text(doc);tmp.replace(OUT/'report.html')

def sampled_predictions(z,n,seed):
    out=predict(z,n,seed)
    pars,_,_=decode(z)
    for i,study in enumerate(STUDIES):
        if study=='Hosny_baseline':continue
        cad,T=SCHEDULE[study]
        for g in (0,1):
            tau,sigma=pars[i,g]
            inc=unit_increments(float(tau),cad,T,n,seed+100*i)*sigma
            path=np.linalg.norm(inc,axis=2).sum(axis=1)
            end=np.linalg.norm(inc.sum(axis=1),axis=1)
            p=out[study][g]
            p.update(speed=float(path.mean()/T),path_length=float(path.mean()),endpoint=float(end.mean()),directness=float((end/path).mean()),early_msd=float((inc[:,0,:]**2).sum(axis=1).mean()))
            # Exact marginal displacement sample at the specified MSD lag.
            v=float(displacement_variance(11.1111111111,np.array([1/tau]),np.array([sigma*sigma*tau/2])))
            normals=normal_bank(1,n,seed+100*i+51).reshape(n,2)
            p['msd']=float(v*(normals**2).sum(axis=1).mean())
    return out

def validate():
    fits=[json.loads(p.read_text()) for p in sorted(OUT.glob('fit_*.json'))]
    assert len(fits)==8*len(SCENARIOS),'Incomplete campaign'
    for scenario in SCENARIOS:
        best=min((f for f in fits if f['scenario']==scenario),key=lambda f:f['total_score'])
        z=best['z'];a=analytic(z)
        sampled=[sampled_predictions(z,16384,seed) for seed in (74101,92317)]
        comparisons=[]
        for checks in (False,True):
            _,ra=data_residuals(a,checks)
            rr=[data_residuals(p,checks)[1] for p in sampled]
            for j,row in enumerate(ra):
                held=checks or (row['study']=='Hosny_baseline' and row.get('quantile')!=1)
                # Constant-path approximation provides no path quartile prediction.
                ap=None if row['study']=='Hosny_baseline' and row['metric']=='path_length' and row.get('quantile')!=1 else row['predicted']
                comparisons.append(dict(study=row['study'],group=row['group'],metric=row['metric'],quantile=row.get('quantile',''),obs=row['obs'],analytic=ap,tracks_seed1=rr[0][j]['predicted'],tracks_seed2=rr[1][j]['predicted'],role='held out' if held else 'fit'))
        atomic(OUT/f'validation_{scenario}.json',dict(selected_start=best['start'],z=z,comparison=comparisons,n=16384,seeds=[74101,92317],generated=datetime.now(timezone.utc).isoformat()))
        # Save example tracks with unaltered fitted parameters, including Hosny drift.
        pars,drifts,_=decode(z);tracks={}
        for i,s in enumerate(STUDIES):
            cad,T=SCHEDULE[s]
            for g in (0,1):
                inc=unit_increments(float(pars[i,g,0]),cad,T,32,17001+i)*pars[i,g,1]
                if s=='Hosny_baseline':inc[:,:,1]+=drifts[g]*cad
                tracks[f'{s}_{g}']=np.concatenate([np.zeros((32,1,2)),np.cumsum(inc,axis=1)],axis=1)
        np.savez_compressed(OUT/f'tracks_{scenario}.npz',**tracks)
        publish()
    publish()
if __name__=='__main__':
    import sys
    if '--publish-only' in sys.argv:publish()
    else:validate()
