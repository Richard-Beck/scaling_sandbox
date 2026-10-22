"""Resolve the wave-pinning peak-rate dip using positive-sink accuracy replays."""
import json
import numpy as np
import inactive_sink_loss as s


def main():
    out=s.OUT;manifest=s.verify(out);records=[]
    for job in manifest['jobs']:
        if job['case']['model']!='wave_pinning' or job['case']['length']!=39.:continue
        L=39.;directory=out/'workers'/f'{job["slot"]:03d}'
        with np.load(s.n.ROOT/job['case']['checkpoint']) as z:original=z['state'].copy()
        for i in [2,3,4]:
            k=manifest['strengths'][i]
            model=s.SinkModel('wave_pinning',k,job['sink_left'],job['sink_right'],L,original.shape[1])
            tr=s.simulate_reaction_diffusion(model,initial_state=np.vstack((original,np.zeros_like(original[0]))),
                protocol=s.n.ConstantInput(0.),length=L,cells=original.shape[1],
                times=np.r_[0.,np.geomspace(1e-6,s.HORIZON,4801)],backend='bdf',rtol=2e-10,atol=1e-12,max_step=10.)
            mask=s.overlap_weights(0,job['pole_edge'],L,original.shape[1])>0
            peak=np.max(tr.fields['active'][:,mask],axis=1)
            frac=peak/peak[0]
            t90,t70=[s.crossing(tr.time,frac,q) for q in [.9,.7]]
            rate=.2/(t70-t90)
            with np.load(directory/f'k_{i:02d}.npz') as z:
                peak0=z['peak_height'];frac0=peak0/peak0[0]
                a,b=[s.crossing(z['time'],frac0,q) for q in [.9,.7]]
                nominal=.2/(b-a)
            refinement=s.r.read(out/'peak_refinement'/f'slot_{job["slot"]:03d}.json')['records']
            record=next((r for r in refinement if r['strength']==k),None)
            if record is not None:nominal=record['peak_fractional_loss_rate_10_to_30']
            error=abs(rate/nominal-1.)
            records.append(dict(geometry=job['geometry'],model='wave_pinning',length=L,strength=k,
                rate=rate,reference_rate=nominal,relative_error=error,
                t10_peak_loss=t90,t30_peak_loss=t70,passed=True))
            print(json.dumps(records[-1]),flush=True)
            assert error<.01
    s.r.write(out/'results/wave_peak_dip_accuracy.json',dict(generated_at=s.n.now(),passed=True,
        definition='Six positive-sink accuracy replays at L=39; strengths 0.01, 0.0316228, 0.1 and both geometries. 100x tighter relative tolerance, 10x smaller maximum step, 4x denser sampling. Original spatial checkpoint mesh.',cases=records))


if __name__=='__main__':main()
