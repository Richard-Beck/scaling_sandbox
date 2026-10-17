"""Independent corrected preparations, reused identities and full native trial audit."""
from pathlib import Path
import copy,gzip,json
import numpy as np
import audit_reversal as a

def preparations(run,package):
    p=run/'preparation';m=a.read(p/'manifest.json');v=a.read(p/'provenance.json')
    a.require(a.sha(p/'manifest.json')==v['manifest_sha256'],'Preparation manifest changed')
    a.require(m['center_boundary_motion'] is False and m['boundary_method']=='independent_endpoints','Wrong preparation mechanics')
    records={}
    for c in m['cases']:
        result=p/'work/cases'/f'case_{c["case_index"]:04d}/result.json';data=a.read(result);ts=data['trials'];s=data['summary']
        valid=[]
        for t in ts:
            z=np.load(p/t['trace']);o=z['release_orientation'];contrast=z['release_contrast'];time=z['release_time']
            a.close(o[-1],t['final_orientation'],'Preparation final moment');a.close(contrast[-1],t['final_contrast'],'Preparation final contrast')
            for k in z.files:a.require(np.all(np.isfinite(z[k])),'Nonfinite preparation '+k)
            a.require(a.sha(p/t['checkpoint'])==t['checkpoint_sha256'],'Preparation checkpoint changed')
            if t['verified']:
                checks=t['convergence']['checks'];a.require(checks[-1]['kind']=='holdout' and checks[-1]['verified'],'Missing holdout')
                a.require(checks[-1]['maximum_scaled_excursion']<=m['steady_state']['holdout_tolerance'],'Holdout excursion')
                a.require(checks[-1]['steady'],'Holdout final window')
                if contrast[-1]>=.05:valid.append(t)
        best=max((t['cue_aligned_strength'] for t in valid),default=max((t['cue_aligned_strength'] for t in ts if t['verified']),default=0))
        passed=[t for t in valid if best>=.05 and t['cue_aligned_strength']>=.99*best and t['polarity_arrival'] is not None]
        selected=min(passed,key=lambda t:t['duration']) if passed else None
        a.close(s['best_observed_strength'],best,'Strongest preparation')
        a.require((s['status']=='polarized')==(selected is not None),'Preparation eligibility selection')
        a.require(not s['duration_boundary_unresolved'] and not s['upper_optimum_unresolved'],'Unresolved preparation')
        if selected:
            a.close(s['cue_duration'],selected['duration'],'Shortest preparation');a.require(s['selected_checkpoint']==selected['checkpoint'],'Selected preparation identity')
        records[c['length']]=(s,p)
    return records

def inputs(run,m,package):
    v=a.read(run/'provenance.json');a.require(a.sha(run/'manifest.json')==v['manifest_sha256'],'Manifest changed')
    for name,digest in v['sources'].items():a.require(a.sha(package/name)==digest and a.sha(run/'frozen_sources'/name)==digest,'Frozen source '+name)
    a.require(m['experiment']=='three_geometry_corrected_debelly' and m['center_boundary_motion'] is False,'Wrong corrected campaign')
    a.require(set(m['models'])==a.MODELS and set(m['excluded_models'])=={'otsuji','legi'},'Model scope')
    a.require(m['recovery_horizon']==10000 and m['sample_interval']==2.5 and m['gate_window']==500,'Clocks')
    original=a.read(package/m['previous_campaign']/'manifest.json');old={c['case_index']:c for c in original['cases']}
    prep=preparations(run,package);initial={};cases=m['cases']
    a.require(len(cases)==330 and [c['case_index'] for c in cases]==list(range(330)),'Coverage')
    for c in cases:
        o=old[c['case_index']]
        for k in ['model','length','condition','durations','amplitude','response_amplitude','response_factor','width','width_mode','cells','fixed_scales']:
            a.require(c.get(k)==o.get(k),'Unchanged input '+k)
        if c['condition']=='fixed_slope':a.require(c['gradient_slope']==o['gradient_slope'],'Slope')
        if c['model']!='debelly':a.require(c==o,'Untouched model case changed')
        else:
            s,p=prep[c['length']];a.require(c['center_boundary_motion'] is False,'Wrong boundary rule')
            a.require(c['eligible']==(s['status']=='polarized'),'Eligibility differs')
            if c['eligible']:
                a.require(package/c['checkpoint']==p/s['selected_checkpoint'],'Old checkpoint reused')
                a.require(c['checkpoint_sha256']==s['selected_checkpoint_sha256'],'Selected checkpoint differs')
        if c['eligible']:
            cp=package/c['checkpoint'];a.require(a.sha(cp)==c['checkpoint_sha256'],'Checkpoint digest')
            with np.load(cp) as z:
                a.close(z['orientation'][-1],c['initial_orientation'],'Initial moment');a.close(z['contrast'][-1],c['initial_contrast'],'Initial contrast')
                a.require(c['initial_orientation']<-.05 and c['initial_contrast']>=.05,'Eligibility')
                initial[c['case_index']]=a.checkpoint_fields(z,c)
    for path in (run/'work/chunks').glob('r0_*.json'):
        data=a.read(path);reuse=data.get('reuse_source')
        if not reuse:continue
        c=cases[data['case_index']];a.require(c['model']!='debelly','Deprecated deBelly simulation imported')
        origin=package/reuse['campaign']/reuse['trace'];a.require(a.sha(origin)==reuse['archive_sha256']==data['trace_sha256'],'Imported archive changed')
        original_chunk=a.read(origin.with_suffix('.json'));a.require(original_chunk['checkpoint_sha256']==c['checkpoint_sha256'],'Imported checkpoint differs')
        oldrows={t['duration']:t for t in original_chunk['trials']}
        for row in data['trials']:
            lhs=copy.deepcopy(row);rhs=copy.deepcopy(oldrows[row['duration']]);lhs.pop('trace');rhs.pop('trace');a.require(lhs==rhs,'Imported metadata differs')
    return cases,initial

a.audit_inputs=inputs
if __name__=='__main__':a.main()
