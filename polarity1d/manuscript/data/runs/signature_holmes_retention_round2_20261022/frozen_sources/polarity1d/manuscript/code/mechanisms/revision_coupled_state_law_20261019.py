"""Evaluate the frozen native Goryachev state law; no trajectory integration."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/revision_width_fullfield_gory_retry_20261019'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/coupled_state_law_20261019'
sys.path[:0]=[str(RUN/'frozen_sources/code/core'),str(RUN/'frozen_sources/code/core/src')]
from models.goryachev import Goryachev
from polarity1d.solvers import laplacian_no_flux


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    for name in ('supporting/data','supporting/metadata'):(OUT/name).mkdir(parents=True,exist_ok=True)
    manifest=json.loads((RUN/'manifest.json').read_text())
    provenance=json.loads((RUN/'provenance.json').read_text())
    sources={str((RUN/'manifest.json').relative_to(ROOT)):sha(RUN/'manifest.json')}
    for path in ['frozen_sources/code/core/models/goryachev.py','frozen_sources/code/core/src/polarity1d/solvers.py']:
        assert sha(RUN/path)==provenance['frozen_sources'][path.removeprefix('frozen_sources/')]
        sources[str((RUN/path).relative_to(ROOT))]=sha(RUN/path)
    model=Goryachev();rows=[]
    for case in manifest['cases']:
        folder=RUN/'workers'/f'{case["slot"]:02d}'
        result=json.loads((folder/'result.json').read_text())
        archive=folder/'recovery.npz'
        assert sha(archive)==result['recovery_sha256']
        sources[str(archive.relative_to(ROOT))]=sha(archive)
        maskpath=RUN/case['regional_input'];assert sha(maskpath)==case['regional_input_sha256']
        sources[str(maskpath.relative_to(ROOT))]=sha(maskpath)
        with np.load(maskpath) as zz:core=zz['incumbent_core_mask'].copy()
        with np.load(archive) as z:
            t,x=z['time'],z['x'];u,v=z['field_active'],z['field_inactive']
            dx=case['length']/case['cells']
            for elapsed in [0.,100.]:
                ids=np.flatnonzero(t==elapsed);assert len(ids)==1
                i=int(ids[0]);state=np.vstack((u[i],v[i]))
                reaction=model.reaction(elapsed,x,state,np.zeros_like(x))
                diffusion=np.asarray(model.diffusivities)[:,None]*laplacian_no_flux(state,dx)
                rhs=reaction+diffusion
                quadratic=model.a*model.ec*u[i]**2*v[i]
                linear=model.b*model.ec*u[i]*v[i]
                deactivation=model.c*u[i]
                np.testing.assert_allclose(reaction[0],quadratic+linear-deactivation,atol=2e-16,rtol=2e-14)
                for name,mask in [('distal_centre_selected_ROI',x>=case['length']/2),('frozen_initial_core',core)]:
                    indices=np.flatnonzero(mask);assert len(indices)>0 and np.all(np.diff(indices)==1)
                    left,right=int(indices[0]),int(indices[-1]);weights=mask.astype(float)*dx
                    native_diffusion=float(np.dot(diffusion[0],weights))
                    left_flux=0. if left==0 else model.du*(u[i,left-1]-u[i,left])/dx
                    right_flux=0. if right==len(x)-1 else model.du*(u[i,right+1]-u[i,right])/dx
                    face_flux=left_flux+right_flux
                    recruitment=float(np.dot(quadratic+linear,weights))
                    loss=float(np.dot(deactivation,weights))
                    direct=float(np.dot(rhs[0],weights))
                    residual=direct-(recruitment-loss+face_flux)
                    assert abs(native_diffusion-face_flux)<1e-13
                    assert abs(residual)<1e-13
                    lo=max(0,i-1);hi=min(len(t)-1,i+1)
                    secant=float(np.dot(u[hi]-u[lo],weights)/(t[hi]-t[lo]))
                    rows.append(dict(length=case['length'],physical_width=case['width'],pulse_duration=case['cue_duration'],
                        time_after_withdrawal=elapsed,region=name,region_first_cell_index=left,region_last_cell_index=right,
                        region_cell_count=len(indices),dx=dx,region_left_face=left*dx,region_right_face=(right+1)*dx,
                        active_amount=float(np.dot(u[i],weights)),active_peak=float(u[i,mask].max()),
                        inactive_mean=float(v[i,mask].mean()),endogenous_recruitment=recruitment,
                        quadratic_recruitment=float(np.dot(quadratic,weights)),linear_recruitment=float(np.dot(linear,weights)),
                        deactivation=loss,net_reaction=recruitment-loss,integrated_native_diffusion=native_diffusion,
                        left_face_diffusive_inflow=left_flux,right_face_diffusive_inflow=right_flux,
                        integrated_native_active_RHS=direct,regional_RHS_balance_residual=residual,
                        saved_secant_left_time=float(t[lo]),saved_secant_right_time=float(t[hi]),
                        saved_secant_dt=float(t[hi]-t[lo]),saved_active_amount_secant=secant,
                        source_archive=str(archive.relative_to(ROOT)),source_archive_sha256=sha(archive)))
    csvpath=OUT/'supporting/data/native_regional_state_law.csv'
    with csvpath.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    bykey={(r['length'],r['physical_width'],r['time_after_withdrawal'],r['region']):r for r in rows}
    narrative=['# Native maintenance at the observed states','',
      'This is a direct evaluation of the frozen Goryachev equation on the already saved full fields, with the cue set to zero after withdrawal. It adds no integration or new perturbation. All rates are native activity amount per native time; activity concentration and amount remain distinct.','',
      'The native transfer is', '',
      r'$$J=(a e_c u^2+b e_c u+S)v-cu.\tag{1}$$','',
      'Recruitment below means endogenous recruitment only (the first two terms); stimulation is absent. Diffusion uses the exact frozen conservative cell-centred no-flux Laplacian. Regional amount is the sum of active concentration times the fixed native cell width. The distal ROI contains cells whose centres satisfy $x_i \\ge L/2$, including the entire centre cell for these odd $N=201$ meshes; it therefore begins at the left face of that centre cell, slightly before the ideal geometric half. The CSV records every endpoint index, cell count, $\\Delta x$ and face coordinate. The frozen original-core mask is evaluated separately.','',
      '| L | Cue width | Off time | Distal active amount | Endogenous recruitment | Deactivation | Net reaction | Native diffusion |',
      '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['region']=='distal_centre_selected_ROI':
            narrative.append(f"| {r['length']:g} | {r['physical_width']:g} | {r['time_after_withdrawal']:g} | {r['active_amount']:.5g} | {r['endogenous_recruitment']:.5g} | {r['deactivation']:.5g} | {r['net_reaction']:.5g} | {r['integrated_native_diffusion']:.3g} |")
    narrative+=['',
      'At L9 withdrawal, the broader response has more distal active amount but less endogenous recruitment and more deactivation. Both responses initially lose distal activity. At 100 cue-free units the narrow response has positive net endogenous reaction, while the broad response remains negative. This makes the concentration-dependent native maintenance account concrete at those states. The small regional active-diffusion terms are instantaneous contributions across this particular ROI boundary; they do not establish that transport or inactive redistribution is unimportant for the preceding or subsequent competition.','',
      'The regional balance independently equals endogenous recruitment minus deactivation plus the exact sum of incoming face fluxes. Maximum absolute residual is '+f"{max(abs(r['regional_RHS_balance_residual']) for r in rows):.3g}"+' native amount/time. The CSV also records actual neighbouring-frame secants and their intervals. In particular the forward 0.25-unit secant at withdrawal is a finite-time observation, not the instantaneous derivative; no exact agreement with it is assumed.','',
      'The L25 pair is contextual, with its distinct pulse 815 and width contrast 3 / 5; it is not a controlled change in length alone. Both L25 distal regions still have slightly negative net reaction at 100 cue-free units, although the broader challenge ultimately wins. These two snapshots therefore do not explain that width-fate difference or define a universal maintenance threshold. The table identifies native production/loss at each saved state without allocating the causal explanation of eventual replacement. The inherited $N=201$ spatial resolution has not undergone an independent convergence study; the native input and clocks remain uncalibrated to an optical actuator or reporter.','',
      '[Exact regional values and provenance](supporting/data/native_regional_state_law.csv) · [Source/artifact checks](supporting/metadata/validation.json)','']
    (OUT/'README.md').write_text('\n'.join(narrative))
    sources[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    artifacts={str(p.relative_to(OUT)):sha(p) for p in [csvpath,OUT/'README.md']}
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(source_hashes=sources,artifact_hashes=artifacts,
        rows=len(rows),integration_performed=False,max_regional_RHS_balance_residual=max(abs(r['regional_RHS_balance_residual']) for r in rows)),indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),rows=len(rows)),indent=2))


if __name__=='__main__':main()
