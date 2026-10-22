"""Dense local threshold sampling for sharp peak collapse after a long delay."""
import argparse, json
import numpy as np
import inactive_sink_loss as s


def main(slot):
    out=s.OUT;manifest=s.verify(out);job=manifest['jobs'][slot]
    c=job['case'];L=c['length'];directory=out/'workers'/f'{slot:03d}'
    target=out/'peak_refinement';target.mkdir(exist_ok=True)
    records=[]
    with np.load(s.n.ROOT/c['checkpoint']) as z:original=z['state'].copy()
    for i,k in enumerate(manifest['strengths']):
        with np.load(directory/f'k_{i:02d}.npz') as z:
            peak=z['peak_height'];time=z['time'];frac=peak/peak[0]
            a,b=[s.crossing(time,frac,q) for q in [.9,.7]]
        if b is None or a is None or b-a>=.5*a:continue
        # Preserve original physics and checkpoint. Resolve each crossing locally.
        times=np.unique(np.r_[0.,np.linspace(.95*a,1.05*a,401),np.linspace(.95*b,1.05*b,401)])
        model=s.SinkModel(c['model'],k,job['sink_left'],job['sink_right'],L,original.shape[1])
        tr=s.simulate_reaction_diffusion(model,initial_state=np.vstack((original,np.zeros_like(original[0]))),
            protocol=s.n.ConstantInput(0.),length=L,cells=original.shape[1],times=times,
            backend='bdf',rtol=2e-10,atol=1e-12,max_step=10.)
        mask=s.overlap_weights(0.,job['pole_edge'],L,original.shape[1])>0
        peak=np.max(tr.fields['active'][:,mask],axis=1);frac=peak/peak[0]
        aa,bb=[s.crossing(tr.time,frac,q) for q in [.9,.7]]
        assert aa is not None and bb is not None
        rate=.2/(bb-aa)
        path=target/f'slot_{slot:03d}_k_{i:02d}.npz'
        np.savez_compressed(path,time=tr.time,original=original,**tr.fields)
        balance=tr.metadata['conservation_errors']['active+inactive+removed'];assert balance<1e-7
        records.append(dict(model=c['model'],length=L,geometry=job['geometry'],strength=k,
            peak_t10_loss=aa,peak_t30_loss=bb,peak_fractional_loss_rate_10_to_30=rate,
            nominal_rate=.2/(b-a),relative_rate_change=abs(rate/(.2/(b-a))-1.),
            archive=str(path.relative_to(out)),archive_sha256=s.n.sha(path),balance_error=balance))
        print(json.dumps(records[-1]),flush=True)
    s.r.write(target/f'slot_{slot:03d}.json',dict(job=job,generated_at=s.n.now(),records=records))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--slot',type=int,required=True);a=p.parse_args();main(a.slot)
