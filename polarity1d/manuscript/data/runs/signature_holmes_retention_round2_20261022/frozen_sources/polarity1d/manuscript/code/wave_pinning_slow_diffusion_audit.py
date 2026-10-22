"""Independent variant inputs, full native trial gates and endpoint audit."""
import copy
import numpy as np
import audit_reversal as a


def inputs(run,manifest,package):
    provenance=a.read(run/'provenance.json')
    a.require(a.sha(run/'manifest.json')==provenance['manifest_sha256'],'Manifest changed')
    for name,digest in provenance['sources'].items():
        a.require(a.sha(package/name)==digest==a.sha(run/'frozen_sources'/name),'Frozen source changed: '+name)
    a.require(manifest['models']==['wave_pinning'],'Variant model scope')
    a.require(manifest['parameters']==dict(du=.01,dv=1.),'Diffusion parameters')
    a.require(manifest['experiment']=='wave_pinning_tenfold_reduced_diffusion','Experiment identity')
    a.require(manifest['recovery_horizon']==10000. and manifest['sample_interval']==2.5 and manifest['gate_window']==500.,'Observation clocks')
    plan=a.read(run/'preparation_plan.json')
    a.require(a.sha(run/'preparation_plan.json')==manifest['preparation_plan_sha256'],'Preparation plan changed')
    source=a.read(package/'data/runs/three_test_reversal_20261006/manifest.json')['cases']
    original={c['length']:c for c in source if c['model']=='wave_pinning' and c['condition']=='fixed_relative'}
    slopes={c['length']:c['gradient_slope'] for c in source if c['model']=='wave_pinning' and c['condition']=='fixed_slope'}
    cases=manifest['cases'];initial={};shared={}
    a.require(len(cases)==3*len(original),'Variant length coverage')
    a.require({(c['length'],c['condition']) for c in cases}=={(L,g) for L in original for g in ['fixed_physical','fixed_relative','fixed_slope']},'Geometry coverage')
    a.require([c['case_index'] for c in cases]==list(range(len(cases))),'Index coverage')
    for c in cases:
        old=original[c['length']]
        for key in ['model','cells','amplitude','amplitude_factor','durations','fixed_scales']:
            a.require(c[key]==old[key],'Changed retained input: '+key)
        a.require(c['variant_parameters']==manifest['parameters'],'Case diffusion differs')
        a.require(c['original_checkpoint']==old['checkpoint'] and c['original_checkpoint_sha256']==old['checkpoint_sha256'],'Wrong original branch')
        a.close(c['response_amplitude'],.1,'Patch amplitude')
        if c['condition']=='fixed_slope':a.close(c['gradient_slope'],slopes[c['length']],'Slope changed')
        else:
            a.require(c['width_mode']==('physical' if c['condition']=='fixed_physical' else 'relative'),'Width mode')
            a.close(c['width'],3. if c['condition']=='fixed_physical' else .2*c['length'],'Width changed')
        prep=a.read(package/c['variant_preparation'])
        a.require(c['eligible']==prep['eligible'],'Preparation eligibility differs')
        if not old['eligible']:
            a.require(not c['eligible'],'Unqualified original branch substituted')
            continue
        a.require(a.sha(package/old['checkpoint'])==old['checkpoint_sha256'],'Original checkpoint changed')
        a.require(prep['verified'],'Unverified variant preparation')
        checkpoint=package/c['checkpoint']
        a.require(a.sha(checkpoint)==c['checkpoint_sha256']==prep['checkpoint_sha256'],'Variant checkpoint changed')
        if c['length'] in shared:a.require(shared[c['length']]==c['checkpoint_sha256'],'Geometries do not share exact start')
        shared[c['length']]=c['checkpoint_sha256']
        with np.load(checkpoint) as z:
            a.close(float(z['orientation'][-1]),c['initial_orientation'],'Initial scalar')
            a.close(float(z['contrast'][-1]),c['initial_contrast'],'Initial contrast')
            a.require(c['eligible']==bool(float(z['orientation'][-1])<-.05 and float(z['contrast'][-1])>=.05),'Eligibility')
            fields=a.checkpoint_fields(z,c)
            if manifest['preparation_mode']=='settle':
                t=z['time'];rule=plan['steady_state']
                a.require(t[-1]>=rule['holdout']+2*rule['chunk'],'Insufficient preparation holdout')
                hold=t>=t[-1]-rule['holdout'];last=t>=t[-1]-rule['chunk']
                for name,scale in c['fixed_scales'].items():
                    y=z['field_'+name]
                    # Compare holdout to its initial recorded state; independent of runner checks.
                    excursion=np.max(abs(y[hold]-y[hold][0]))/scale
                    a.require(excursion<=rule['holdout_tolerance']+1e-10,'Preparation holdout excursion')
                    change=np.max(np.ptp(y[last],axis=0))/scale
                    drift=np.max(abs(y[-1]-y[-2]))/(t[-1]-t[-2])*rule['chunk']/scale
                    a.require(max(change,drift)<=rule['tolerance']+1e-10,'Preparation final window not steady')
                with np.load(package/old['checkpoint']) as previous:
                    # Native finite volumes have equal cell weights.
                    mass0=previous['state'].sum()*c['length']/c['cells']
                mass=z['state'].sum()*c['length']/c['cells']
                a.require(np.isclose(mass,mass0,rtol=1e-5),'Variant preparation mass drift')
            if c['eligible']:initial[c['case_index']]=fields
    return cases,initial


a.audit_inputs=inputs
if __name__=='__main__':a.main()
