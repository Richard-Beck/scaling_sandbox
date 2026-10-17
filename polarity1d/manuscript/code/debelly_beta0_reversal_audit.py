"""Independent beta=0 preparation, gate, coverage and shortest-pulse audit."""
from pathlib import Path
import numpy as np
import audit_reversal as a
read=a.read;sha=a.sha;require=a.require;close=a.close

def inputs(run,m,package):
    p=read(run/'provenance.json')
    require(sha(run/'manifest.json')==p['manifest_sha256'],'Manifest changed')
    for name,digest in p['sources'].items():require(sha(package/name)==digest==sha(run/'frozen_sources'/name),'Frozen source changed: '+name)
    require(m['experiment']=='debelly_beta0_shared_steady_reversal' and m['models']==['debelly'] and m['center_boundary_motion'] is False,'Variant identity')
    require(m['recovery_horizon']==10000 and m['sample_interval']==2.5 and m['gate_window']==500,'Observation clocks')
    plan=read(run/'preparation_plan.json');require(sha(run/'preparation_plan.json')==m['preparation_plan_sha256'],'Plan changed')
    original={c['length']:c for c in read(package/'data/runs/debelly_corrected_three_test_20261008/manifest.json')['cases'] if c['model']=='debelly' and c['condition']=='fixed_relative'}
    require(len(m['cases'])==3*len(original),'Length coverage')
    require({(c['length'],c['condition']) for c in m['cases']}=={(L,g) for L in original for g in ['fixed_physical','fixed_relative','fixed_slope']},'Cue coverage')
    initial={};shared={}
    for c in m['cases']:
        old=original[c['length']];prep=read(package/c['variant_preparation'])
        for key in ['model','cells','amplitude','durations','fixed_scales']:require(c[key]==old[key],'Changed retained input: '+key)
        require(c['original_checkpoint']==old['checkpoint'] and c['original_checkpoint_sha256']==old['checkpoint_sha256'],'Wrong native branch')
        require(c['eligible']==prep['eligible'],'Eligibility differs')
        require(c['orientation_cutoff']==(.05 if m['eligibility_policy']=='original' else .005),'Orientation cutoff differs')
        close(c['response_amplitude'],3.6,'Patch peak')
        if c['condition']=='fixed_slope':close(c['gradient_slope'],.24,'Gradient slope')
        else:
            require(c['width_mode']==('physical' if c['condition']=='fixed_physical' else 'relative'),'Width mode')
            close(c['width'],3. if c['condition']=='fixed_physical' else .2*c['length'],'Patch width')
        if not old['eligible']:
            require(not c['eligible'],'Unqualified native branch substituted');continue
        require(prep['verified'],'Unsettled beta0 state')
        cp=package/c['checkpoint'];require(sha(cp)==prep['checkpoint_sha256']==c['checkpoint_sha256'],'Checkpoint changed')
        if c['length'] in shared:require(shared[c['length']]==c['checkpoint_sha256'],'Three cues use different starts')
        shared[c['length']]=c['checkpoint_sha256']
        with np.load(cp) as z,np.load(package/old['checkpoint']) as oldz:
            t=z['time'];rule=plan['steady_state'];require(t[-1]>=rule['holdout']+2*rule['chunk'],'Insufficient preparation')
            hold=t>=t[-1]-rule['holdout'];last=t>=t[-1]-rule['chunk']
            require(np.array_equal(z['field_delta_beta_tension'],np.zeros_like(z['field_tension'])),'Beta0 violated')
            require(np.array_equal(z['field_rho_activation'],np.ones_like(z['field_rho_activation'])),'Rho basal activation changed')
            import sys
            sys.path.insert(0,str(package/'code/core'))
            import native as n
            from dataclasses import replace
            parameters=replace(n.db._parameters(),center_boundary_motion=False)
            R,_=n.db._homogeneous_rac_rho(parameters,.5);close(1.3*R,.899061648345644,'Protrusion threshold changed')
            native_alpha,_=n.db._feedback_terms(parameters,z['field_mca_bound'],z['field_tension'])
            require(np.array_equal(z['field_delta_alpha_mca'],native_alpha),'MCA feedback changed')
            for name in n.DB_NAMES:require(np.array_equal(z['field_'+name][0],oldz['state_'+name]),'Preparation restart changed '+name)
            fields=a.checkpoint_fields(z,c)
            for name,scale in c['fixed_scales'].items():
                if name=='physical_length':y=c['length']+z['field_boundary_separation']
                else:y=z['field_'+name]
                if name=='displacement':y=y-y.mean(axis=1)[:,None]
                excursion=float(abs(y[hold]-y[hold][0]).max()/scale)
                variation=float(np.ptp(y[last],axis=0).max()/scale)
                drift=float(abs(y[-1]-y[-2]).max()/(t[-1]-t[-2])*rule['chunk']/scale)
                require(excursion<=rule['holdout_tolerance']+1e-10 and max(variation,drift)<=rule['tolerance']+1e-10,'Preparation stability failed: '+name)
            ori=float(z['orientation'][-1]);contrast=float(z['contrast'][-1])
            close(ori,c['initial_orientation'],'Initial polarity');close(contrast,c['initial_contrast'],'Initial contrast')
            require(c['eligible']==bool(ori < -c['orientation_cutoff'] and contrast>=.05),'Variant eligibility')
            require(c['original_cutoff_eligible']==bool(ori<-.05 and contrast>=.05),'Original cutoff flag')
            mass=(z['field_mca_bound']+z['field_mca_unbound']).sum(axis=1)
            require(abs(mass/mass[0]-1).max()<1e-10,'MCA conservation')
            if c['eligible']:initial[c['case_index']]=fields
    return m['cases'],initial

def audit_trial(row, arrays, case, initial):
    duration = row['duration']
    require(np.isfinite(duration) and duration >= 0., 'Invalid pulse duration')
    t = arrays['recovery_time']
    o = arrays['recovery_orientation']
    contrast = arrays['recovery_contrast']
    require(len(t) == len(o) == len(contrast) == 4001, 'Recovery sample coverage incomplete')
    require(np.array_equal(t, np.arange(4001) * 2.5), 'Recovery sampling times changed')
    require(np.isfinite(o).all() and np.isfinite(contrast).all(), 'Nonfinite scalar histories')
    tail = t >= 9500.
    require(int(tail.sum()) == 201, 'Incomplete final-500 sample set')
    y, q = o[tail], contrast[tail]
    orientation = bool(np.min(y) > case['orientation_cutoff'] and np.min(q) >= .05 and np.ptp(y) < .005)
    error = float(np.max(np.abs(y / abs(case['initial_orientation']) - 1.)))
    full = bool(orientation and error <= .01)
    require(row['orientation_success'] == orientation, 'Orientation gate differs from arrays')
    require(row['full_polarity_success'] == full, 'Scalar-strength gate differs from arrays')
    close(row['maximum_relative_polarity_error_last_500'], error, 'Strength error differs')
    close(row['minimum_contrast_last_500'], float(np.min(q)), 'Tail contrast differs')
    close(row['orientation_range_last_500'], float(np.ptp(y)), 'Orientation range differs')
    close(row['final_orientation'], float(y[-1]), 'Final orientation differs')
    pulse_t = arrays['pulse_time']
    require(len(pulse_t) >= 1 and pulse_t[0] == 0. and np.isclose(pulse_t[-1], duration),
            'Pulse timing differs')
    close(float(arrays['pulse_orientation'][0]), case['initial_orientation'],
          'Pulse does not restart original orientation')
    close(float(arrays['pulse_contrast'][0]), case['initial_contrast'],
          'Pulse does not restart original contrast')
    maximum, excursion = 0., 0.
    fields = {name[len('fieldtail_'):-len('_min')] for name in arrays
              if name.startswith('fieldtail_') and name.endswith('_min')}
    require(fields == set(initial), 'Monitored native-field coverage changed')
    for name in fields:
        scale = case['fixed_scales'][name]
        lower, upper, previous, last = [arrays[f'fieldtail_{name}_{key}']
                                      for key in ('min','max','penultimate','last')]
        require(lower.shape == upper.shape == previous.shape == last.shape == initial[name].shape,
                'Field-tail shape changed')
        require(all(np.isfinite(v).all() for v in (lower,upper,previous,last)), 'Nonfinite field statistic')
        require(np.all(lower <= upper) and np.all(previous >= lower) and np.all(previous <= upper)
                and np.all(last >= lower) and np.all(last <= upper), 'Inconsistent tail extrema')
        variation = float(np.max(upper - lower)) / scale
        drift = float(np.max(np.abs(last - previous))) / 2.5 * 500. / scale
        maximum = max(maximum, variation, drift)
        close(row['field_tail_metrics'][name]['range_scaled'], variation, 'Field range differs')
        close(row['field_tail_metrics'][name]['projected_drift_500'], drift, 'Field drift differs')
        if duration == 0.:
            a, b, start = [arrays[f'no_cue_{name}_{key}'] for key in ('min','max','initial')]
            require(np.array_equal(start, initial[name]), 'No-cue control initial field differs from checkpoint')
            require(a.shape == b.shape == start.shape and np.isfinite(a).all() and np.isfinite(b).all()
                    and np.all(a <= b), 'Invalid no-cue whole-history extrema')
            excursion = max(excursion, float(max(np.max(np.abs(a-start)), np.max(np.abs(b-start)))) / scale)
    close(row['field_tail_max_scaled'], maximum, 'Field flag maximum differs')
    require(row['field_tail_pass'] == bool(maximum <= .001), 'Field-tail flag differs')
    final = [v for name,v in arrays.items() if name == 'finalstate' or name.startswith('finalstate_')]
    require(bool(final) and all(np.isfinite(v).all() for v in final), 'Missing/nonfinite final native state')
    if duration == 0.:
        close(row['no_cue_maximum_scaled_excursion'], excursion, 'No-cue full-history excursion differs')
        left = bool(np.max(y) < -case['orientation_cutoff'] and np.min(q) >= .05 and np.ptp(y) < .005)
        valid = bool(excursion <= .01 and not orientation and left)
        require(row['no_cue_left_tail_valid'] == left, 'No-cue left-state gate differs')
        require(row['no_cue_preparation_valid'] == valid, 'No-cue preparation validity differs')


_original_trial=a.audit_trial

def checked_trial(row,arrays,case,initial):
    audit_trial(row,arrays,case,initial)
    tail=arrays['recovery_time']>=9500
    o=arrays['recovery_orientation'][tail];q=arrays['recovery_contrast'][tail]
    require(row['original_orientation_success']==bool(o.min()>.05 and q.min()>=.05 and np.ptp(o)<.005),'Original outcome flag')
    checks=row['feedback_audits'];require(len(checks)==(1 if row['duration']==0 else 2),'Missing full-history intervention checks')
    for check in checks:
        require(check['initial_full_state_exact'] and check['delta_beta_zero_exact'] and check['native_delta_alpha_exact'] and check['mca_mass_error']<1e-10,'Intervention or state identity failed')
a.audit_inputs=inputs
a.audit_trial=checked_trial
if __name__=='__main__':a.main()
