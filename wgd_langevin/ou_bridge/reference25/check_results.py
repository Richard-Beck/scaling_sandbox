"""Basic mathematical and numerical correctness checks for the review."""
import json
import numpy as np
import campaign as c


def main():
    assert len(c.SPATIAL_JOBS)==90
    x=(np.arange(10000)+.5)/10000*25
    for model,case in c.SPATIAL_JOBS:
        for seg in case['segments']:
            t=(seg[0]+seg[1])/2
            y=c.Protocol(case)(t,x,25)
            actual=y[5000:].mean()-y[:5000].mean()
            predicted=c.forcing(t/case['unit'],c.normalized_segments(case),case['width'])
            assert abs(actual-predicted)<1e-8,(model,case['name'])
    t=np.linspace(0,10,101)
    s=np.array([[0.,10.,0.,0.,1.,0.]])
    theta=np.array([2.,1.,0.,.3,.3,.4,0.,.1])
    p=c.integrate(theta,t,s,.1,0)[:,0]
    assert max(abs(p-.08*(1-np.exp(-t/2))))<1e-8
    x=np.linspace(0,1,126)
    p=2*c.halfmass(np.array([np.ones(126),1+x,2-x]))
    assert abs(p[0])<1e-12 and abs(p[1]+p[2])<1e-12
    rows=[]
    for model,family in c.FIT_JOBS:
        path=c.OUT/'fits'/f'{model}_{family}.json'
        record=json.loads(path.read_text())
        errors=[]
        with np.load(path.with_suffix('.npz')) as stored:
            for case in c.load_cases(model):
                refined=c.integrate(np.array(record['theta']),case['t'],case['spec'],
                    case['width'],c.FAMILIES.index(family),2)[:,0]
                errors.append(float(np.max(abs(refined-stored[case['name']]))))
        error=max(errors)
        rows.append(dict(model=model,family=family,max_step_halving_difference=error,
                         passed=error<.002))
        print(model,family,error,flush=True)
    c.design.save_json(c.OUT/'numerical_check.json',dict(generated=c.now(),
        projection_check=True,analytic_ou_check=True,reflection_check=True,
        tolerance=.002,fits=rows))
    assert all(r['passed'] for r in rows),'Reduced integration check failed'


if __name__=='__main__': main()
