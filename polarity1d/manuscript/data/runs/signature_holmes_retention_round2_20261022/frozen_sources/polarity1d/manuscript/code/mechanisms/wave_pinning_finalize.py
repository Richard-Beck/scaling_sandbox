"""Close the Wave-pinning mechanism artifact inventory and replay audit."""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import wave_pinning_mechanism as w
n=w.n
FIG=w.MANUSCRIPT/'figures/fixed_slope_mechanisms/wave_pinning'

PAGES=[
 ('wave_pinning_spatial','Three spatial routes to persistent reversal',
  'Increasing length strengthens the ramp and starts right activation faster, but leaves a more complicated profile to resolve after withdrawal.',
  'All columns use the exact relative_20pct preparation, unchanged native physics and g=0.00666667. Native lengths5,13,39 sample the decreasing arm, minimum and increasing arm. The upper row shows the steady starting active profile and the physical ramp; the middle row separates on-cue snapshots from the lower basal-recovery snapshots. Right-edge active concentration first exceeds the descriptive level0.6 at30.3,5.6,1.3 respectively; this is an activation marker, not the response clock. AtL13, two boundary caps compete after removal. AtL39, broad activation collapses into an interior competitor and a right cap; the selected orientation pulse39 passes without replacing the original scalar strength.'),
 ('wave_pinning_decisions','Why just-too-short pulses fail',
  'Failure means loss of polarity at short length, restoration of the old winner at the minimum, and insufficient right bias at the fixed horizon in the long cell.',
  'Each passing pulse is paired with its immediately lower tested failing duration:30.75/30.5 atL5,13.875/13.75 atL13 and39/38.75 atL39. Heatmaps show native active protein after withdrawal, using logarithmic time spacing; the curves display the orientation readout separately. The official gate uses all201 samples in recovery9,500–10,000. The short failure becomes homogeneous, and the middle failure returns to the left cap despite weak positive orientation at cutoff. The long failure is still slightly below0.05 at10,000, but crosses that level during unchanged native continuation:0.04950 at10,000 becomes0.05823 at50,000. The extended trace documents later evolution; it is not an added threshold search or a newly audited201-sample endpoint.'),
 ('wave_pinning_mechanism','Weak short drive versus long pool redistribution',
  'A short-cell cue-strength intervention and an identical-cutoff inactive-mobility intervention support different bottlenecks on the two arms.',
  'PanelA retains both scalar duration curves and native-field flags (open markers fail the independent0.1% field-tail diagnostic). B shows the strongly unequal inactive pool in the long cell at cutoff. C holdsL5 and pulse13.875 fixed: raising slope2.6-fold to match theL13 peak changes failure to full-strength reversal under otherwise native physics. D shows long-cell full-strength success depends on later competitor clearance. E changes only inactive diffusion after the exactL39,pulse39 native cutoff;10→100 gives a single right cap by about100 recovery units. No-cue controls preserve the polarized branch, with only0.668% scalar-strength shift. F shows the native competitor, its clearance under faster inactive redistribution, and qualitative agreement with a bounded401-cell native-parameter check. Mobility controls alter diagnostic recovery physics; native continuation establishes the horizon effect without altering parameters.'),
 ('wave_pinning_full_strength','Long-cell physical cleanup and the full-strength gate',
  'A quarter-unit pulse difference shifts interior-domain collapse across the official final-500 window.',
  'BothL39 pulses49.5 and49.25 first recover into an interior active region plus a right boundary cap. The passing49.5 pulse clears the competitor before the official9,500–10,000 recovery window. The49.25 pulse clears it during that window: its final orientation is already0.43238, but the tail range0.05120 fails the persistence gate and the maximum initial-strength error0.1155 fails the1% criterion. It ultimately has the native single right cap by50,000. Thus the reported threshold selects cleanup within a fixed observation window, not a measured asymptotic basin boundary. A401-cell check preserves the physical sequence and the failed native-window outcome, while exact near-threshold timing changes.'),
]

def finalize():
    audit=[]
    for k in range(8):
        p=w.OUT/f'trajectory_{k:02d}.npz';record=json.loads(p.with_suffix('.json').read_text());job=record['job'];c=job['case'];D=job['duration']
        if n.sha(p)!=record['archive_sha256']:raise RuntimeError('archive hash changed')
        with np.load(p) as z:
            t=z['time'];tail=(t>=D+9500-1e-8)&(t<=D+10000+1e-8)
            if np.count_nonzero(tail)!=201:raise RuntimeError('tail coverage mismatch')
            u=z['active'][tail];v=z['inactive'][tail];x=z['x'];center=2*x/(x[0]+x[-1])-1
            orientation=(np.maximum(u,0)@center)/np.maximum(u,0).sum(axis=1)
            contrast=np.ptp(u,axis=1)
            ordinary=bool(orientation.min()>.05 and contrast.min()>=.05 and np.ptp(orientation)<.005)
            error=float(np.max(abs(orientation/abs(c['initial_orientation'])-1)))
            full=bool(ordinary and error<=.01)
            maximum=0.
            for name,a in [('active',u),('inactive',v)]:
                maximum=max(maximum,float(np.max(np.ptp(a,axis=0)))/c['fixed_scales'][name],float(np.max(abs(a[-1]-a[-2])))/2.5*500/c['fixed_scales'][name])
            field=bool(maximum<=.001)
            source=job['original_row']
            if (ordinary,full,field)!=(source['orientation_success'],source['full_polarity_success'],source['field_tail_pass']):raise RuntimeError('independent gate mismatch')
            if abs(float(orientation[-1])-source['final_orientation'])>1e-9:raise RuntimeError('orientation mismatch')
            mass=(z['active']+z['inactive']).sum(axis=1)*c['length']/c['cells'];masserror=float(np.max(abs(mass/mass[0]-1)))
            if masserror>2e-13:raise RuntimeError('mass conservation')
            audit.append(dict(slot=k,length=c['length'],duration=D,tail_samples=201,orientation_success=ordinary,full_polarity_success=full,field_tail_pass=field,mass_relative_error=masserror,checkpoint_sha256=c['checkpoint_sha256']))
    result=dict(generated_at=n.now(),complete=True,passed=True,method='Separate raw native-field reconstruction of every formal tail gate; exact201 samples; archive/checkpoint hashes and conserved totalmass.',trajectories=audit)
    n.write(w.OUT/'reconstruction_audit.json',result)
    pages=[]
    for name,title,claim,caption in PAGES:
        pages.append(dict(title=title,main_claim=claim,caption=caption,pdf=str(FIG/f'{name}.pdf'),png=str(FIG/f'{name}.png')))
    assets=dict(model='wave_pinning',generated_at=n.now(),complete=True,pages=pages,notes=str(w.OUT/'mechanism_notes.txt'),summary=str(w.OUT/'mechanism_summary.json'),audit=str(w.OUT/'reconstruction_audit.json'))
    n.write(FIG/'assets.json',assets)
    sourcepaths=list(Path(__file__).parent.glob('wave_pinning*.py'))
    sourcepaths+=[w.MANUSCRIPT/'code/core/native.py',w.MANUSCRIPT/'code/core/models/wave_pinning.py',w.MANUSCRIPT/'code/core/src/polarity1d/solvers.py',w.MANUSCRIPT/'code/gradient_reversal.py']
    outputpaths=list(FIG.glob('wave_pinning*.*'))+[FIG/'assets.json',w.OUT/'mechanism_notes.txt',w.OUT/'mechanism_summary.json',w.OUT/'reconstruction_audit.json']
    nativeinputs=[w.CAMPAIGN/'manifest.json',w.CAMPAIGN/'results/endpoints.csv',w.CAMPAIGN/'results/records.json.gz']
    nativeinputs+=[w.MANUSCRIPT/c['job']['case']['checkpoint'] for c in [json.loads((w.OUT/f'trajectory_{k:02d}.json').read_text()) for k in [0,2,4]]]
    provenance=dict(generated_at=n.now(),slurm_jobs=dict(native_reconstruction='20724893',mobility_controls='20724968',cue_controls='20724990',mesh_controls='20724996'),
                    sources={str(p.relative_to(w.MANUSCRIPT)):n.sha(p) for p in sourcepaths},inputs={str(p.relative_to(w.MANUSCRIPT)):n.sha(p) for p in nativeinputs},
                    outputs={str(p.relative_to(w.MANUSCRIPT)):n.sha(p) for p in outputpaths},archive_records=[str(p.relative_to(w.MANUSCRIPT)) for p in sorted(w.OUT.glob('*control_*.json'))]+[str((w.OUT/f'trajectory_{k:02d}.json').relative_to(w.MANUSCRIPT)) for k in range(8)])
    n.write(w.OUT/'final_provenance.json',provenance)
    print(json.dumps(dict(assets=str(FIG/'assets.json'),pages=len(pages),audit_passed=result['passed']),indent=2))

if __name__=='__main__':finalize()
