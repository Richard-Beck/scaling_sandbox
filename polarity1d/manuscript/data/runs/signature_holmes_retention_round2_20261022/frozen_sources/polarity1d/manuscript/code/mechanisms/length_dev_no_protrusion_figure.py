"""Independent native-array audit and development figure for threshold=Inf."""
from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace, asdict
import argparse, csv, json, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code/core'),str(Path(__file__).parent)]
import native as n
import settling
import length_dev_no_protrusion as campaign
import length_dev_output_layout as layout
OUT=ROOT/'figures/length_mechanisms_20261008/dev/no_protrusion'
def read(p):return json.loads(Path(p).read_text())
def output(name):return layout.dev_output_path(OUT,name)
def key(p):return str(p.relative_to(OUT))

def audit(run):
    manifest=read(run/'manifest.json');prov=read(run/'provenance.json')
    assert n.sha(run/'manifest.json')==prov['manifest_sha256']
    for name,digest in prov['sources'].items():assert n.sha(run/'frozen_sources'/name)==digest
    native=(run/'frozen_sources/code/core/models/debelly.py').read_text()
    baseline_provenance=read(campaign.BASE/'provenance.json')
    for name,digest in prov['sources'].items():
        if name.startswith('code/core/'):
            assert baseline_provenance['sources'][name]==digest,name
    variant=run/'frozen_sources/code/mechanisms/debelly_no_protrusion_native.py'
    assert variant.read_text()==native.replace(prov['source_change']['before'],prov['source_change']['after'])
    assert n.sha(variant)==prov['variant_sha256']
    missing=[j['slot'] for j in manifest['jobs'] if not (run/'workers'/f"{j['slot']:04d}.json").exists()]
    if missing:raise RuntimeError(f'Missing {len(missing)} workers: {missing}')
    records=[];checks=[];max_mass_error=0.;max_recurrence_error=0.
    p=replace(n.db._parameters(),center_boundary_motion=False)
    for job in manifest['jobs']:
        prefix=run/'workers'/f"{job['slot']:04d}";r=read(prefix.with_suffix('.json'))
        assert r['job']==job and r['parameters']==asdict(p)
        assert n.sha(prefix.with_suffix('.npz'))==r['trace_sha256']
        state,tr=n.read_checkpoint(prefix.with_suffix('.npz'),'debelly')
        assert all(np.isfinite(v).all() for v in tr.fields.values())
        rac=tr.fields['rac'];frac=(np.arange(rac.shape[1])+.5)/rac.shape[1]
        orientation=np.sum(rac*(2*frac-1),axis=1)/np.sum(rac,axis=1)
        contrast=np.ptp(rac,axis=1)
        np.testing.assert_allclose(orientation,tr.orientation,rtol=0,atol=2e-15)
        np.testing.assert_allclose(contrast,tr.contrast,rtol=0,atol=2e-15)
        assert abs(r['final_orientation']-orientation[-1])<2e-15
        assert abs(r['final_contrast']-contrast[-1])<2e-15
        for k in n.DB_NAMES:np.testing.assert_array_equal(state[k],tr.fields[k][-1])
        mass=np.sum(tr.fields['mca_bound']+tr.fields['mca_unbound'],axis=1)
        masserror=float(np.max(abs(mass/mass[0]-1)));assert masserror<1e-10
        max_mass_error=max(max_mass_error,masserror)
        # Native early-release sampling takes .25-s steps; subsequent steps are .5 s.
        exponent=np.minimum(tr.time,10)/.25+np.maximum(tr.time-10,0)/.5
        for side in ['left','right']:
            name='raw_'+side+'_displacement';u=tr.fields[name]
            expected=u[0]*p.boundary_relaxation**exponent[:,None]
            residual=float(np.max(abs(u-expected)));assert residual<1e-11
            max_recurrence_error=max(max_recurrence_error,residual)
            np.testing.assert_array_equal(u,tr.fields[side+'_boundary_displacement'])
        if job['kind']=='native_branch':
            cp=run/'checkpoints'/f"L{job['length']:g}.npz"
            assert n.sha(cp)==job['checkpoint_sha256']
            original,_=n.read_checkpoint(cp,'debelly')
            for k in n.DB_NAMES:np.testing.assert_array_equal(original[k],tr.fields[k][0])
        if job['kind']=='from_rest' and job['duration']>0:
            path=prefix.with_suffix('.pulse.npz');assert n.sha(path)==r['pulse_sha256']
            pulse_state,pulse=n.read_checkpoint(path,'debelly')
            assert abs(pulse.time[-1]-job['duration'])<1e-9
            for k in n.DB_NAMES:np.testing.assert_array_equal(pulse_state[k],tr.fields[k][0])
            for side in ['left','right']:assert np.max(abs(pulse.fields['raw_'+side+'_displacement']))==0
        # Recompute every full-field convergence window and each holdout excursion.
        rule=manifest['steady_state'];scales=job['case']['fixed_scales'];previous=10.;reverified=False
        for c in r['convergence']['checks']:
            end=c['time'];start=end-rule['chunk'] if c['kind']=='convergence' else previous
            mask=(tr.time>=start-1e-8)&(tr.time<=end+1e-8)
            segment=SimpleNamespace(time=tr.time[mask],fields={k:v[mask] for k,v in tr.fields.items()})
            metric=settling.diagnostics(segment,'debelly',job['length'],scales,rule)
            assert metric['steady']==c['steady']
            for name,m in metric['metrics'].items():
                for label,value in m.items():assert abs(value-c['metrics'][name][label])<1e-10
            if c['kind']=='holdout':
                fields=settling.monitored(segment,'debelly',job['length'])
                excursion=float(settling.deviations(fields,{k:v[0] for k,v in fields.items()},scales).max())
                assert abs(excursion-c['maximum_scaled_excursion'])<1e-10
                reverified=bool(excursion<=rule['holdout_tolerance'] and metric['steady'])
                assert reverified==c['verified']
            previous=end
        assert reverified==r['convergence']['verified']
        records.append(r);checks.append(dict(slot=job['slot'],passed=True,verified=reverified,mca_mass_relative_error=masserror))
        if (job['slot']+1)%25==0:print('Audited',job['slot']+1,'/',len(manifest['jobs']),flush=True)
        del tr,segment,rac
    # Independent implementation of the intervention: zero speed gives the same dynamics.
    cn,_,_=campaign.configure(run);ds=cn.ds
    cp=run/'checkpoints/L39.npz';state,_=cn.read_checkpoint(cp,'debelly')
    ablated=cn.run('debelly',39.,state,cn.ConstantInput(0),10.,interval=.5)
    cn.ds=cn.db;cn.model_for=lambda k:cn.db.DeBelly(replace(cn.db._parameters(),center_boundary_motion=False,protrusion_speed=0.))
    zero_speed=cn.run('debelly',39.,state,cn.ConstantInput(0),10.,interval=.5)
    for name in ablated.fields:np.testing.assert_array_equal(ablated.fields[name],zero_speed.fields[name])
    rows=[];selected={}
    for b in manifest['native_baseline']:
        assert n.sha(ROOT/b['result'])==b['result_sha256']
        L=b['length'];rs=[r for r in records if r['job']['length']==L]
        rest=[r for r in rs if r['job']['kind']=='from_rest' and r['convergence']['verified']]
        assert len([r for r in rs if r['job']['kind']=='from_rest'])==18
        # Highest observed left moment, independently of cue-duration optimization.
        best=max(rest,key=lambda r:-r['final_orientation']) if rest else None
        branch=next((r for r in rs if r['job']['kind']=='native_branch'),None)
        seed=next(r for r in rs if r['job']['kind']=='small_seed')
        row=dict(length_um=L,native_eligible=b['summary']['status']=='polarized',
                 native_best_strength=b['summary']['best_observed_strength'],
                 native_selected_strength=b['summary']['selected_final_polarity'],
                 native_selected_contrast=b['summary']['selected_final_contrast'],
                 off_rest_strength=max(0.,-best['final_orientation']) if best else None,
                 off_rest_contrast=best['final_contrast'] if best else None,
                 off_rest_eligible=bool(best and -best['final_orientation']>=.05 and best['final_contrast']>=.05),
                 off_rest_slot=best['job']['slot'] if best else None,
                 off_rest_witness_duration_s=best['job']['duration'] if best else None,
                 off_rest_verified_trials=len(rest),off_rest_unsettled_trials=18-len(rest),
                 off_branch_strength=-branch['final_orientation'] if branch else None,
                 off_branch_contrast=branch['final_contrast'] if branch else None,
                 off_branch_verified=branch['convergence']['verified'] if branch else None,
                 off_branch_eligible=bool(branch and branch['convergence']['verified'] and -branch['final_orientation']>=.05 and branch['final_contrast']>=.05),
                 off_branch_slot=branch['job']['slot'] if branch else None,
                 small_seed_abs_strength=abs(seed['final_orientation']),small_seed_contrast=seed['final_contrast'],
                 small_seed_finite_window_verified=seed['convergence']['verified'])
        rows.append(row)
    result=dict(generated_at=n.now(),passed=True,complete=True,trials=len(records),checks=checks,
                only_one_threshold_line_changed=True,all_other_parameters_identical=True,
                zero_speed_independent_control_bit_identical=True,all_checkpoint_and_phase_identities_passed=True,
                all_full_field_convergence_checks_recomputed=True,maximum_mca_mass_relative_error=max_mass_error,
                maximum_endpoint_relaxation_residual_um=max_recurrence_error,
                verified_trials=sum(r['convergence']['verified'] for r in records),
                unsettled_slots=[r['job']['slot'] for r in records if not r['convergence']['verified']],
                native_polarizable_lengths=[r['length_um'] for r in rows if r['native_eligible']],
                no_protrusion_polarizable_lengths=[r['length_um'] for r in rows if r['off_rest_eligible']],
                retained_native_branch_lengths=[r['length_um'] for r in rows if r['off_branch_eligible']],
                input_provenance_sha256=n.sha(run/'provenance.json'))
    n.write(run/'results/independent_audit.json',result)
    n.csvwrite(run/'results/length_summary.csv',rows)
    n.write(run/'results/length_summary.json',rows)
    n.write(run/'results/status.json',dict(generated_at=n.now(),producer='finished',complete=True,passed=True,total=len(records),verified=result['verified_trials'],unsettled_slots=result['unsettled_slots']))
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))

def render(run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    audit_result=read(run/'results/independent_audit.json');assert audit_result['passed'] and audit_result['complete']
    manifest=read(run/'manifest.json');rows=read(run/'results/length_summary.json')
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    blue='#0072B2';orange='#D55E00';green='#009E73';gray='#aaaaaa'
    fig=plt.figure(figsize=(13.8,8.2),layout='constrained')
    gs=fig.add_gridspec(3,3,height_ratios=[.85,1.25,1.3])
    a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,:2]);c=fig.add_subplot(gs[1,2])
    profiles=[fig.add_subplot(gs[2,i]) for i in range(3)]
    lengths=np.array([r['length_um'] for r in rows])
    eligibility=np.array([[r['native_eligible'] for r in rows],[r['off_rest_eligible'] for r in rows],[r['off_branch_eligible'] if r['off_branch_slot'] is not None else -1 for r in rows]],float)
    # Uniformly spaced columns retain the full irregular original length grid.
    a.imshow(eligibility,aspect='auto',cmap=ListedColormap(['#eeeeee','#bdbdbd',green]),vmin=-1,vmax=1,interpolation='nearest')
    a.set_xticks(np.arange(len(rows)),[f'{L:g}' for L in lengths]);a.set_yticks([0,1,2],['Native: from rest','Protrusion off: from rest','Protrusion off: native start'])
    a.set_xlabel('Nominal cell length (µm; original sampled grid)')
    a.set_title('A   Polarizable lengths are unchanged on the sampled grid',loc='left',fontweight='bold')
    for i,row in enumerate(eligibility):
        for j,value in enumerate(row):a.text(j,i,'✓' if value==1 else ('–' if value==-1 else '×'),ha='center',va='center',color='white' if value==1 else '#444444',fontsize=11)
    a.text(1.0,1.05,'✓ eligible   × below gate   – no native polarized start',ha='right',transform=a.transAxes,fontsize=8)
    native=[r['native_selected_strength'] for r in rows]
    off=[r['off_rest_strength'] for r in rows]
    b.plot(lengths,native,'o-',color=blue,ms=4,label='Native selected basal state')
    b.plot(lengths,off,'s-',color=orange,ms=4,label='Protrusion off: best verified from rest')
    br=[r for r in rows if r['off_branch_slot'] is not None]
    b.plot([r['length_um'] for r in br],[r['off_branch_strength'] for r in br],'o',color=green,mfc='none',ms=8,lw=1.5,label='Protrusion off: continued native state')
    b.axhline(.05,color='#666666',ls=':',label='Standard polarity gate = 0.05')
    b.set(xlabel='Nominal cell length (µm)',ylabel='Left Rac polarity (−normalized first moment)',ylim=(-.004,.115),xlim=(0,40))
    b.set_title('B   Polarity persists; its length dependence changes',loc='left',fontweight='bold')
    b.legend(frameon=False,fontsize=8,loc='upper left');b.grid(alpha=.12)
    c.plot(lengths,[r['native_selected_contrast'] or 0. for r in rows],'o-',color=blue,ms=4)
    c.plot(lengths,[r['off_rest_contrast'] for r in rows],'s-',color=orange,ms=4)
    c.axhline(.05,color='#666666',ls=':')
    c.set(xlabel='Nominal cell length (µm)',ylabel='Rac contrast (max − min)',xlim=(0,40),ylim=(0,None))
    c.set_title('C   Clear spatial Rac asymmetry',loc='left',fontweight='bold');c.grid(alpha=.12)
    data={};sources={}
    for ax,L,letter in zip(profiles,[18.,23.,39.],'DEF'):
        row=next(r for r in rows if r['length_um']==L)
        base=next(r for r in manifest['native_baseline'] if r['length']==L)
        _,native=n.read_checkpoint(ROOT/base['checkpoint'],'debelly')
        offpath=run/'workers'/f"{row['off_rest_slot']:04d}.npz"
        branchpath=run/'workers'/f"{row['off_branch_slot']:04d}.npz"
        _,offtr=n.read_checkpoint(offpath,'debelly');_,branchtr=n.read_checkpoint(branchpath,'debelly')
        for tr,color,label,style in [(native,blue,'Native','-'),(offtr,orange,'Off: from rest','-'),(branchtr,green,'Off: native start','--')]:
            ax.plot(tr.x,tr.fields['rac'][-1],color=color,ls=style,label=label,lw=2)
        ax.set(xlabel='Nominal position (µm)',ylabel='Basal Rac concentration',xlim=(0,L),ylim=(.54,.94))
        ax.set_title(f'{letter}   L = {L:g} µm',loc='left',fontweight='bold');ax.grid(alpha=.12)
        if L==18:ax.legend(frameon=False,fontsize=8)
        for name,tr in [('native',native),('off_rest',offtr),('off_branch',branchtr)]:
            data[f'L{L:g}_{name}_x']=tr.x;data[f'L{L:g}_{name}_rac']=tr.fields['rac'][-1]
        for path in [ROOT/base['checkpoint'],offpath,branchpath]:sources[str(path.relative_to(ROOT))]=n.sha(path)
    fig.suptitle('deBelly polarization with Rac-driven protrusion disabled',fontsize=16,fontweight='bold')
    fig.get_layout_engine().set(rect=(0,.065,1,.935))
    fig.text(.015,.018,'Only the endpoint Rac threshold changes: 1.3 × basal Rac → ∞. All remaining mechanochemical feedbacks are retained.\nBasal endpoints follow cue withdrawal and full-field convergence/holdout checks. Eligibility: left polarity ≥ 0.05 and Rac contrast ≥ 0.05.',fontsize=9)
    exports={}
    for ext in ['png','pdf','svg']:
        path=output('no_protrusion_polarizability.'+ext);fig.savefig(path,dpi=220,facecolor='white');exports[key(path)]=n.sha(path)
    plt.close(fig)
    for source,name in [(run/'results/length_summary.csv','no_protrusion_length_summary.csv'),(run/'results/independent_audit.json','no_protrusion_independent_audit.json')]:
        import shutil
        path=output(name);shutil.copy2(source,path);exports[key(path)]=n.sha(path)
    path=output('no_protrusion_profiles.npz');np.savez_compressed(path,**data);exports[key(path)]=n.sha(path)
    caption=('Protrusion-disabled polarization in the corrected independent-endpoint deBelly model. The frozen native solver differs by exactly one line: the shared endpoint Rac threshold is infinity instead of 1.3 times homogeneous basal Rac, making both protrusion commands exactly zero. All other parameters, endogenous Rac/Rho reactions, MCA feedback, Rho active stress, tension feedback, mechanics, transport and endpoint relaxation remain native. A: all 22 original lengths (1,3,5,7,9,11,12,13,15,17,18,19,21,23,25,27,29,31,33,35,37,39 µm). The 12 eligible lengths are the same with and without protrusion: 18,19,21,23,25,27,29,31,33,35,37,39 µm. Every eligible native branch also remains eligible when protrusion is disabled. B: native selected basal state, strongest observed verified left-polarized endpoint from the ablated preparation grid, and continuation of the exact native polarized checkpoint. This compares endpoint polarity, not minimum preparation times or a reversal test. C: max-minus-min Rac contrast; native ineligible values without a selected endpoint are shown at zero. D–F: full Rac profiles at L=18,23,39 in nominal coordinates; retained branch and de novo results overlap closely. Each preparation independently starts from homogeneous native rest, receives a nominal peak-3.6 left raised-cosine patch of support 20% of length, then returns to cue-free basal input. The 18 original durations are zero plus 17 logarithmic samples from 0.5 to 2000 s; no new duration threshold or strongest-endpoint refinement is claimed. A further small conserved MCA seed at every length is an auxiliary diagnostic only; finite-window settling of a tiny perturbation is not proof of homogeneous linear stability. There are 430 integrations: 396 from-rest trials, 12 native-branch continuations, and 22 seed controls. Unsettled trials stay explicit in the data and do not invalidate separate verified positive witnesses at the same length. No eligible endpoint was found at sampled lengths up to 17 µm; this bounded scan does not locate an exact continuous-length bifurcation or exclude every possible initial state. All displayed ablated endpoints meet the unchanged full-field convergence protocol: two 1000-s windows at 1e-4 frozen scale tolerance followed by a 10000-s holdout at 1e-3, with a 60000-s cap. Independent audits recompute all moments, contrasts, full-field convergence windows and holdout excursions, MCA conservation, phase/checkpoint identities and zero-command endpoint relaxation. A ten-second independent zero-protrusion-speed implementation is bit-identical to the infinity-threshold variant. '+f"Completed {audit_result['generated_at']}; {audit_result['verified_trials']}/430 integrations verified, unsettled slots {audit_result['unsettled_slots']}.")
    path=output('no_protrusion_polarizability_caption.txt');path.write_text(caption+'\n');exports[key(path)]=n.sha(path)
    for path in [run/'manifest.json',run/'provenance.json',run/'results/independent_audit.json',run/'results/length_summary.json',Path(__file__)]:sources[str(path.relative_to(ROOT))]=n.sha(path)
    metadata=dict(artifact_root='../..',generated_at=n.now(),producer='finished',legend=caption,sha256=exports,source_sha256=sources,
                  native_baseline='debelly_corrected_three_test_20261008',campaign=str(run.relative_to(ROOT)),independent_audit_passed=True,
                  reproduction=[f'{campaign.PYTHON} {Path(__file__)} audit --run {run}',f'{campaign.PYTHON} {Path(__file__)} render --run {run}'])
    n.write(output('no_protrusion_polarizability.manifest.json'),metadata)
    (run/'README.md').write_text('# Protrusion-disabled polarization\n\n'+caption+'\n\nFrozen scientific sources and the one-line threshold variant are under `frozen_sources/`. All full native release histories, pulse histories and worker records remain in `workers/`. `results/` contains audited length summaries and completion status.\n')
    print(output('no_protrusion_polarizability.png'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['audit','render']);ap.add_argument('--run',type=Path,default=campaign.DEFAULT)
    a=ap.parse_args();audit(a.run.resolve()) if a.action=='audit' else render(a.run.resolve())
