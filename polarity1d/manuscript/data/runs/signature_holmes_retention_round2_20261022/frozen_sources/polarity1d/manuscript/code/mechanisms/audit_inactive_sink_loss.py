"""Independent equation/budget audit and selected positive-sink accuracy replays."""
from pathlib import Path
import json, sys
import numpy as np
import inactive_sink_loss as s


def main():
    out=s.OUT; manifest=s.verify(out); records=[]
    for job in manifest['jobs']:
        if job['case']['length'] not in [9.,39.]:continue
        c=job['case'];L=c['length'];directory=out/'workers'/f'{job["slot"]:03d}'
        with np.load(directory/'k_06.npz') as z:
            original=z['original'].copy();u=z['active'];v=z['inactive'];time=z['time'];dx=L/u.shape[1]
            edges=np.arange(u.shape[1])*dx
            weights=np.clip(job['pole_edge']-edges,0.,dx)
            sink=np.maximum(0.,np.minimum(edges+dx,job['sink_right'])-np.maximum(edges,job['sink_left']))
            # Equations written independently of the sink wrapper and metric helper.
            if c['model']=='goryachev':reaction=(.33/101*u*u+.67/101*u)*v-.01733*u
            else:reaction=v*(.067+u*u/(1.+u*u))-u
            lap=np.empty_like(u);lap[:,0]=(u[:,1]-u[:,0])/dx**2;lap[:,-1]=(u[:,-2]-u[:,-1])/dx**2
            lap[:,1:-1]=(u[:,:-2]-2*u[:,1:-1]+u[:,2:])/dx**2
            diffusion=.0025 if c['model']=='goryachev' else .1
            np.testing.assert_allclose(z['core_rate'],(reaction+diffusion*lap)@weights,rtol=1e-11,atol=1e-10)
            np.testing.assert_allclose(z['removal_flux'],v@sink,rtol=1e-12,atol=1e-12)
            available=(u+v).sum(axis=1)*dx
            np.testing.assert_allclose(available+z['removed'].sum(axis=1)*dx,original.sum()*dx,rtol=1e-7,atol=1e-8)
            # Dense logarithmic quadrature is an independent approximate flux check.
            integral=np.sum(.5*(z['removal_flux'][1:]+z['removal_flux'][:-1])*np.diff(time))
            flux_error=abs(integral/z['removed_total'][-1]-1.)
            assert flux_error<.002
        model=s.SinkModel(c['model'],1.,job['sink_left'],job['sink_right'],L,original.shape[1])
        tr=s.simulate_reaction_diffusion(model,initial_state=np.vstack((original,np.zeros_like(original[0]))),
            protocol=s.n.ConstantInput(0.),length=L,cells=original.shape[1],
            times=np.r_[0.,np.geomspace(1e-6,s.HORIZON,2401)],backend='bdf',rtol=2e-10,atol=1e-12,max_step=25.)
        row,_=s.measurements(dict(time=tr.time,**tr.fields),job,1.)
        nominal=json.loads((directory/'k_06.json').read_text())['metrics']
        errors={key:abs(row[key]/nominal[key]-1.) for key in ['t10_loss','t50_loss','fractional_loss_rate_10_to_30']}
        assert max(errors.values())<.01
        record=dict(model=c['model'],length=L,geometry=job['geometry'],strength=1.,
                    relative_errors=errors,integrated_flux_error=float(flux_error),passed=True)
        records.append(record);print(json.dumps(record),flush=True)
    s.r.write(out/'results/numerical_accuracy.json',dict(generated_at=s.n.now(),passed=True,
        definition='Eight positive-sink replays: 100x tighter relative tolerance, 4x smaller maximum step, 2x denser time samples. Independent native reaction, active-diffusion, sink-flux and material-budget reconstruction.',
        cases=records,spatial_mesh='Unchanged exact checkpoint mesh; no mesh-convergence claim.'))


if __name__=='__main__':main()
