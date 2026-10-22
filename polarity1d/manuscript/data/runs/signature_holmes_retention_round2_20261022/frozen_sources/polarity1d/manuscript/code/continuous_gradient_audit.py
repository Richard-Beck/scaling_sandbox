"""Independently recompute continuous-input readouts from native fields."""
from pathlib import Path
import argparse,json,sys
import numpy as np
CODE=Path(__file__).resolve().parent
sys.path.insert(0,str(CODE))
import continuous_gradient as c
n=c.n


def audit(out):
    manifest=c.verify(out);checks=[]
    for case in manifest['cases']:
        if not case['eligible']:continue
        slot=case['case_index'];key=case['model'];dest=out/'work/cases'/f'case_{slot:03d}'
        row=c.read(dest/'result.json');archive=dest/'trajectory.npz'
        assert n.sha(archive)==row['archive_sha256'],slot
        checkpoint=n.ROOT/case['checkpoint']
        assert n.sha(checkpoint)==case['checkpoint_sha256'],slot
        _,prep=n.read_checkpoint(checkpoint,key)
        with np.load(archive) as z:
            time=z['time'];active=z['active'];recorded=z['orientation'];x=z['x']
            assert time[0]==0 and np.all(np.diff(time)>0),slot
            assert active.shape==(len(time),201),slot
            assert np.isfinite(active).all(),slot
            coord=np.linspace(-1,1,201) if key=='spring' else (2*(np.arange(201)+.5)/201-1)
            if key=='spring':
                orientation=np.trapezoid(active*coord,dx=1/200,axis=1)/np.trapezoid(active,dx=1/200,axis=1)
            else:orientation=np.sum(active*coord,axis=1)/np.sum(active,axis=1)
            assert np.max(abs(orientation-recorded))<1e-12,(slot,'orientation')
            assert np.max(abs(np.ptp(active,axis=1)-z['contrast']))<1e-12,slot
            assert np.max(abs(active[0]-prep.active[-1]))<1e-12,(slot,'exact initial active')
            target=-float(prep.orientation[-1]);assert target>0.05
            assert abs(target-case['target_orientation'])<1e-12,slot
            hit=np.where(orientation>=target)[0]
            if len(hit):
                i=int(hit[0]);assert i>0
                assert row['status']=='success' and row['crossing_index']==i,slot
                assert row['duration']==float(time[i]) and row['lower_time']==float(time[i-1]),slot
                assert np.all(orientation[:i]<target),slot
                estimate=time[i-1]+(target-orientation[i-1])/(orientation[i]-orientation[i-1])*(time[i]-time[i-1])
                assert abs(estimate-row['interpolated_time'])<1e-8*max(1,estimate),slot
            else:
                assert row['status']=='not_reached_within_horizon' and row['duration'] is None,slot
                assert time[-1]==manifest['maximum_application_duration'],slot
            fields={name[6:]:z[name] for name in z.files if name.startswith('field_')}
            for name,values in fields.items():assert np.isfinite(values).all(),(slot,name)
            for name in prep.fields:
                if name in fields and (key!='debelly' or name in n.DB_NAMES):
                    assert np.max(abs(fields[name][0]-prep.fields[name][-1]))<1e-10,(slot,name,'initial profile')
            if key=='spring':masses=[np.trapezoid(fields['active']+fields['inactive'],dx=1/200,axis=1)*fields['length'].reshape(len(time),-1)[:,0]]
            elif key=='debelly':masses=[np.sum(fields['mca_bound']+fields['mca_unbound'],axis=1)]
            else:
                model=n.model_for(key)
                masses=[sum(np.sum(fields[model.state_names[j]],axis=1) for j in group) for group in model.conserved_groups]
            error=max(float(np.max(abs(m/m[0]-1))) for m in masses)
            assert error<1e-5,(slot,'mass closure',error)
            assert case['gradient_slope']==manifest['slope_by_model'][key],slot
            cue=c.g.response_cue(case);assert cue.stop is None
            for at in [0.,1.,manifest['maximum_application_duration']]:
                assert np.allclose(cue(at,x,case['length']),case.get('background_stimulus',0.)+case['gradient_slope']*x,rtol=0,atol=0),slot
            checks.append(dict(case_index=slot,model=key,length=case['length'],status=row['status'],
                               maximum_orientation_reconstruction_error=float(np.max(abs(orientation-recorded))),
                               conserved_material_relative_error=error,checkpoint_verified=True,
                               first_crossing_verified=True,constant_gradient_verified=True))
    assert len(checks)==97
    result=dict(generated_at=n.now(),passed=True,complete=True,checked_cases=len(checks),
                preparation_ineligible=13,checks=checks,
                method='Independent normalized first moment from all retained active/Rac fields, exact original profile/hash checks, first crossing and censor horizon recomputation, native conserved material and constant cue audit.')
    n.write(out/'results/independent_audit.json',result)
    print(json.dumps(dict(audit_passed=True,checked_cases=len(checks))))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=c.DEFAULT);audit(p.parse_args().out.resolve())
