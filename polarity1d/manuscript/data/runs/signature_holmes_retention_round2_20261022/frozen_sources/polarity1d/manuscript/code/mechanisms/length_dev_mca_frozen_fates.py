"""Frozen-snapshot MCA first-exit diagnostics on the Figure X backgrounds.

Each row holds that snapshot's transport coefficients fixed indefinitely.
Tagged molecules still diffuse, bind, detach, and leave by outward advection.
The earlier future-protocol companion is preserved in the same topic folder.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import json
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import solve_banded

import length_dev_mca_molecular_fates as base

ROOT, DEV = base.ROOT, base.DEV
RUN = ROOT / 'data/runs/debelly_mca_frozen_fates_20261010_v4'
OUT = DEV / 'mca_molecular_fates'
STEM = 'MCA_molecular_fates_frozen'
CONFIGS = [('main',1601),('coarse',801),('refined',3201)]
LONG_RESIDENCE = 1e8


def output(name):
    return base.output(name)


def frozen_rates(bg, index, cells, zero_membrane=False):
    """Freeze saved cortex velocity and trailing-interval membrane velocity."""
    faces = np.linspace(0,bg.L,cells+1)
    edge = bg.edges[index]
    bface = np.interp(faces,np.r_[0,bg.x,bg.L],np.r_[edge[0],bg.f['velocity'][index],edge[1]])
    if index == 0 or zero_membrane:
        uface = np.zeros_like(bface)
    else:
        uface = np.interp(faces,np.r_[0,bg.x,bg.L],np.r_[bg.me[index-1,0],bg.mv[index-1],bg.me[index-1,1]])
    return base.generator(bface,uface,bg.L/cells,bg.p)


def pinned_solve(q,left,right):
    """Stable first-exit solve even when escape is below matrix precision.

    Write the backward solution as z+c, with z pinned to zero at the state
    having the greatest absorption rate. Solve for z and the scaled constant
    c*max(absorption) separately, retaining tiny absorption explicitly.
    Unlike a direct inverse of -Q, this does not lose 1e-18 escape rates when
    they are added to O(1..10000) transport/exchange diagonal entries.
    """
    killing = left+right
    scale = float(killing.max())
    if scale == 0:
        return np.column_stack([np.full(len(left),np.inf),np.full((len(left),2),np.nan)]),dict(closed=True)
    anchor = int(killing.argmax())
    a = -q
    pinned = a.copy()
    rhs = np.column_stack([np.ones(len(left)),left,right])
    augmented = np.column_stack([rhs,killing/scale])
    augmented[anchor] = 0
    for offset in range(-2,3):
        col = anchor+offset
        if 0 <= col < len(left):
            pinned[2-offset,col] = 0
    pinned[2,anchor] = 1
    solved = solve_banded((2,2),pinned,augmented,check_finite=False)
    # Refine the pin solve using extended-precision residuals. This is cheap
    # for the banded operator and resolves cancellation in steep flow fields.
    wide_pinned = pinned.astype(np.longdouble)
    for _ in range(2):
        correction_rhs = augmented.astype(np.longdouble)-base.apply_banded(wide_pinned,solved.astype(np.longdouble))
        solved += solve_banded((2,2),pinned,np.asarray(correction_rhs,dtype=float),check_finite=False)
    solved = solved.astype(np.longdouble)
    row = sum(np.longdouble(a[2-offset,anchor+offset])*solved[anchor+offset]
              for offset in range(-2,3) if 0 <= anchor+offset < len(left))
    constant_scaled = (rhs[anchor]-row[:3])/(killing[anchor]/scale-row[3])
    varying = solved[:,:3]-solved[:,3,None]*constant_scaled
    result = np.asarray(varying+constant_scaled/scale,dtype=float)
    # Evaluate the equation before recombining the huge constant, avoiding
    # catastrophic cancellation in Q*T when residence is billions of years.
    residual = base.apply_banded(a.astype(np.longdouble),varying)+killing[:,None].astype(np.longdouble)*(constant_scaled/scale)-rhs
    residual_scale = base.apply_banded(abs(a).astype(np.longdouble),abs(varying))+killing[:,None].astype(np.longdouble)*abs(constant_scaled/scale)+abs(rhs)
    # Unit floor handles identically-zero probability equations and harmless
    # cancellation of a constant probability. The mean-time RHS is exactly 1.
    residual_scale = np.maximum(residual_scale,1.)
    errors = dict(closed=False,stable_equation_max_absolute_residual=float(abs(residual).max()),
                  stable_equation_max_componentwise_backward_error=float((abs(residual)/residual_scale).max()),
                  probability_sum_max_error=float(abs(result[:,1:].sum(axis=1)-1).max()),
                  minimum_probability=float(result[:,1:].min()),maximum_probability=float(result[:,1:].max()),
                  minimum_mean_residence_s=float(result[:,0].min()),maximum_mean_residence_s=float(result[:,0].max()))
    assert errors['stable_equation_max_componentwise_backward_error'] < 1e-11
    assert errors['probability_sum_max_error'] < 1e-7
    assert errors['minimum_probability'] > -1e-7 and errors['maximum_probability'] < 1+1e-7
    assert errors['minimum_mean_residence_s'] > 0
    return result,errors


def selfcheck():
    checks = base.selfcheck()
    p = dict(mca_on=5.,mca_off=1.4,mca_bound_diffusion=.3,mca_unbound_diffusion=.3)
    n,dx = 11,.2
    maximum_relative_error = 0.
    for speed in [.1,1e-8,1e-18]:
        faces = np.zeros(n+1); faces[0],faces[-1] = -speed,speed
        q,l,r = base.generator(faces,faces,dx,p)
        result,audit = pinned_solve(q,l,r)
        i = np.arange(n)
        expected = n*dx/(2*speed)+i*(n-1-i)*dx**2/(2*.3)
        error = float((abs(result[:,0]-np.repeat(expected,2))/np.repeat(expected,2)).max())
        maximum_relative_error = max(maximum_relative_error,error)
        assert error < 1e-10
        exact_probability = (i*dx+.3/speed)/((n-1)*dx+2*.3/speed)
        np.testing.assert_allclose(result[:,2],np.repeat(exact_probability,2),atol=1e-11)
        if speed == .1:
            direct = base.stationary(q,l,r)[0]
            np.testing.assert_allclose(result,direct,rtol=1e-10,atol=1e-10)
    q,l,r = base.generator(np.zeros(n+1),np.zeros(n+1),dx,p)
    closed,audit = pinned_solve(q,l,r)
    assert audit['closed'] and np.isinf(closed[:,0]).all() and np.isnan(closed[:,1:]).all()
    checks.update(pinned_solve_analytic_tiny_escape=True,analytic_tiny_escape_max_relative_error=maximum_relative_error,
                  exact_zero_escape_returns_infinite_time_and_undefined_exit_odds=True)
    return checks


def launch():
    assert not RUN.exists(),'Preserve producer evidence; use a fresh campaign'
    checks = selfcheck()
    jobs,parameters,sources = base.inputs()
    for directory in ['workers','logs','frozen_sources']:
        (RUN/directory).mkdir(parents=True,exist_ok=True)
    for p in [Path(__file__),Path(base.__file__)]:
        frozen = RUN/'frozen_sources'/p.name
        shutil.copy2(p,frozen)
        sources[str(frozen.relative_to(ROOT))] = base.sha(frozen)
    base.write(RUN/'manifest.json',dict(generated_at=base.now(),jobs=jobs,parameters=parameters,source_sha256=sources,
                                      configurations=CONFIGS,selfcheck=checks,
                                      conditioning='Hold snapshot transport, binding and boundary coefficients fixed indefinitely; no future protocol.'))
    base.write(RUN/'status.json',dict(generated_at=base.now(),producer='running',complete=False))
    command = f'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; {base.PYTHON} {RUN / "frozen_sources" / Path(__file__).name} worker --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G',
                                   '--time=00:10:00','--array=0-2','--job-name=mca-frozen',
                                   '--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err'),
                                   '--wrap',command],text=True).strip()
    base.write(RUN/'launch.json',dict(generated_at=base.now(),array=job))
    print(job)


def worker(slot):
    start = time.monotonic()
    m = base.read(RUN/'manifest.json')
    for name,digest in m['source_sha256'].items():
        assert base.sha(ROOT/name) == digest,name
    bg = base.Background(m['jobs'][slot],m['parameters'])
    ids = np.searchsorted(bg.t,bg.display)
    weights = bg.weights(bg.display,len(bg.x))
    arrays = dict(time=bg.display,x=bg.x,duration=np.array(bg.duration),bound_fraction=weights,
                  cortical_edge_velocity=bg.edges[ids],
                  membrane_edge_velocity=np.vstack([np.zeros(2) if i == 0 else bg.me[i-1] for i in ids]))
    audit_rows = []
    for name,cells in m['configurations']:
        fields = []
        for index in ids:
            result,audit = pinned_solve(*frozen_rates(bg,index,cells))
            fields.append(base.sample_states(bg,result,cells))
            audit_rows.append(dict(configuration=name,time_after_onset_s=float(bg.t[index]),**audit))
        arrays[name+'_conditioned'] = np.array(fields)
        arrays[name+'_mixed'] = base.mix(np.array(fields),weights)
        print(bg.L,name,'finished',time.monotonic()-start,flush=True)
    # Freeze displacement literally (zero membrane velocity) as a sensitivity
    # control; primary freezes the actual instantaneous transport coefficients.
    fields = [base.sample_states(bg,pinned_solve(*frozen_rates(bg,i,1601,zero_membrane=True))[0],1601) for i in ids]
    arrays['zero_membrane_mixed'] = base.mix(np.array(fields),weights)
    path = RUN/'workers'/f'{slot:03d}.npz'
    np.savez_compressed(path,**arrays)
    base.write(path.with_suffix('.json'),dict(generated_at=base.now(),archive_sha256=base.sha(path),length_um=bg.L,
                                            wall_seconds=time.monotonic()-start,audits=audit_rows,
                                            native_field_tail_pass=m['jobs'][slot]['field_tail_pass']))
    print(bg.L,'complete',time.monotonic()-start,flush=True)


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm,Normalize
    m = base.read(RUN/'manifest.json')
    for name,digest in m['source_sha256'].items():
        assert base.sha(ROOT/name) == digest,name
    checks = selfcheck()
    arrays,traces,audits,rows = {},{},[],[]
    for slot,job in enumerate(m['jobs']):
        path = RUN/'workers'/f'{slot:03d}.npz'; meta = base.read(path.with_suffix('.json'))
        assert base.sha(path) == meta['archive_sha256']
        with np.load(path) as z:
            a = {k:z[k].copy() for k in z.files}
        L,t,x = job['length'],a['time'],a['x']; value = a['main_mixed']
        refined = base.comparison(value,a['refined_mixed']);coarse = base.comparison(value,a['coarse_mixed'])
        assert refined['mean_time_max_relative_difference'] < .03
        assert refined['right_probability_max_absolute_difference'] < .005
        pre = np.arange(-60,0,.5)
        display = np.concatenate([np.broadcast_to(value[0],(len(pre),len(x),3)),value])
        traces[L] = dict(time=np.r_[pre,t],x=x,value=display,duration=float(a['duration']))
        arrays.update({f'L{L}_{key}':v for key,v in a.items()})
        arrays[f'L{L}_display_time'] = np.r_[pre,t]
        arrays[f'L{L}_display_mean_residence_s'] = display[...,0]
        arrays[f'L{L}_display_p_left'] = display[...,1]
        arrays[f'L{L}_display_p_right'] = display[...,2]
        positive = [r for r in meta['audits'] if not r['closed']]
        audit = dict(length_um=L,tracer_1601_vs_3201=refined,tracer_801_vs_1601=coarse,
                     maximum_stable_equation_residual=max(r['stable_equation_max_absolute_residual'] for r in positive),
                     maximum_componentwise_backward_error=max(r['stable_equation_max_componentwise_backward_error'] for r in positive),
                     maximum_probability_sum_error=max(r['probability_sum_max_error'] for r in positive),
                     closed_snapshot_count=sum(r['closed'] for r in meta['audits'] if r['configuration']=='main'),
                     mean_residence_range_s=[float(value[...,0].min()),float(value[...,0].max())],
                     almost_closed_pixel_fraction=float((value[...,0]>=LONG_RESIDENCE).mean()),
                     zero_membrane_control=base.comparison(value,a['zero_membrane_mixed']),
                     native_field_tail_pass=job['field_tail_pass'])
        audits.append(audit)
        for when in [0,60,float(a['duration']),600,1200]:
            index = int(np.argmin(abs(t-when)))
            for fraction in [.05,.25,.5,.75,.95]:
                rows.append(dict(length_um=L,time_after_onset_s=float(t[index]),position_fraction=fraction,
                                 mean_remaining_residence_s=float(np.interp(L*fraction,x,value[index,:,0])),
                                 p_left=float(np.interp(L*fraction,x,value[index,:,1])),
                                 p_right=float(np.interp(L*fraction,x,value[index,:,2]))))
    np.savez_compressed(output(STEM+'_fields.npz'),**arrays)
    with output(STEM+'_measurements.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none',
                         'hatch.color':'#777777','hatch.linewidth':.45})
    fig = plt.figure(figsize=(12,8.1))
    gs = fig.add_gridspec(2,4,width_ratios=[23,31,39,2.2],left=.105,right=.91,bottom=.19,top=.865,wspace=.18,hspace=.30)
    for row in range(2):
        norm = LogNorm(1,LONG_RESIDENCE) if row == 0 else Normalize(0,1)
        for col,L in enumerate(base.LENGTHS):
            tr = traces[L];ax = fig.add_subplot(gs[row,col])
            values = tr['value'][...,0 if row == 0 else 2]
            mesh = ax.pcolormesh(tr['x'],tr['time'],values,shading='nearest',norm=norm,
                                 cmap='viridis' if row == 0 else 'RdBu_r',rasterized=True)
            long = tr['value'][...,0]>=LONG_RESIDENCE
            if long.any():
                ax.contourf(tr['x'],tr['time'],long.astype(float),levels=[.5,1.5],colors='none',hatches=['///'],alpha=.25)
            ax.set_yscale('symlog',linthresh=120,linscale=.65)
            ax.set(xlim=(0,L),ylim=(1200,-60),xticks=[0,L/2,L],yticks=[-60,0,60,120,300,600,1200])
            ax.set_yticklabels(['−60','0','60','120','300','600','1200']);ax.tick_params(labelsize=9,length=3)
            if col:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(('Frozen mean residence' if row == 0 else 'Frozen first-exit side')+'\nSnapshot time from cue onset (s)',labelpad=8)
            ax.axhline(0,color='white',lw=2);ax.axhline(0,color='#222222',lw=.65)
            ax.axhline(60,color='white',ls=':',lw=1.2)
            ax.axhline(tr['duration'],color='white',ls='--',lw=2);ax.axhline(tr['duration'],color='#222222',ls='--',lw=.7)
            ax.text(.025,.975,chr(65+3*row+col),transform=ax.transAxes,va='top',weight='bold',fontsize=11,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.6))
            if row == 0:ax.set_title(f'L = {L} µm',fontsize=12,pad=12)
            else:ax.set_xlabel('Reference position x (µm)')
        cb = fig.colorbar(mesh,cax=fig.add_subplot(gs[row,3]),extend='max' if row == 0 else 'neither')
        cb.ax.tick_params(labelsize=9)
        if row == 0:
            cb.set_ticks([1,1e2,1e4,1e6,1e8]);cb.set_ticklabels(['1','100','10⁴','10⁶','10⁸'])
            cb.set_label('Mean time to first exit (s; log scale)',labelpad=8)
        else:
            cb.set_ticks([0,.25,.5,.75,1]);cb.set_ticklabels(['0 · LHS','.25','.50','.75','1 · RHS'])
            cb.set_label('Probability of RHS first exit',labelpad=8)
    fig.suptitle('MCA residence and exit fate in frozen snapshots',fontsize=15,y=.975)
    fig.text(.105,.93,'Each row holds that snapshot’s transport fields fixed; molecules continue to diffuse, bind and detach.',fontsize=10,color='#444444')
    fig.text(.105,.108,'Solid: cue onset     Dotted: 60 s     Dashed: cue withdrawal     −60–0 s: repeated stationary basal snapshot',fontsize=9)
    fig.text(.105,.075,'Hatching: frozen mean ≥10⁸ s (almost closed). Times above 10⁸ s share the top color; uncapped values are exported.',fontsize=9)
    fig.text(.105,.042,'Zero diffusive exit at both ends. Each prediction uses only its own snapshot, with no later cue or recovery.',fontsize=9)
    for ext in ['png','pdf','svg']:
        fig.savefig(output(STEM+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    caption = (
        'Frozen-snapshot MCA residence and first-exit fate, companion to Figure X. Columns show L=23,31,39 µm and the exact same '
        'corrected relative-20% trials with pulses 232.244898,285.280612,809.081633 s. Each position/time asks what would happen to '
        'a molecule if all tracer transport, binding and endpoint coefficients at that snapshot were held fixed indefinitely. '
        'It does not use subsequent cue or recovery fields and needs no future-tail extrapolation. A–C show mean remaining time '
        'until first domain exit; the common color scale is logarithmic from 1 to 10^8 s with an upper extension. Values exceeding '
        '10^8 s share the top color and are hatched in both rows; actual uncapped values are retained in NPZ/CSV. D–F show RHS '
        'first-exit probability; LHS probability is the complement. Starting bound/unbound state is weighted by the local native '
        'MCA pool fractions; results conditioned on either state are exported. Snapshot fields and output positions match Figure X. '
        'The −60–0-s strip repeats the initial stationary snapshot; unlike the earlier future-protocol companion, its conditioning '
        'is identical to the calculation at t=0. Time increases downwards, linear through 120 s then logarithmic through 1200 s. '
        'Cortical velocities are the saved native field, with endpoint commands reconstructed from the same Rac/tension snapshot. '
        'Membrane velocities are trailing native output-interval displacement derivatives, including endpoint velocities, and are '
        'held at their instantaneous values. Initial basal membrane velocity is zero. Freezing these prescribed velocities is a '
        'transport-coefficient diagnostic, not a self-consistent moving-boundary continuation. An additional zero-membrane-velocity '
        'control is exported for a literal static-displacement interpretation. Diffusivities are 0.001/0.3 µm²/s and off/on rates '
        '1.4/5 s^-1 for bound/unbound MCA. Diffusion reflects at both endpoints; outward advection alone removes tags. Inward '
        'advective input and uniform replenishment contain no tags; native split-exchange numerical residuals do not create tags. '
        'No identity after first exit is assumed. Almost closed frozen snapshots can yield enormous hypothetical lifetimes: tiny '
        'smooth-switch endpoint commands and numerical membrane velocities then strongly affect escape odds. Hatching identifies '
        'this regime rather than implying such lifetimes occur during the actual evolving reversal. A row-pinned banded Poisson '
        'solve keeps tiny absorption rates explicit, avoiding loss of those rates against transport/exchange matrix diagonals. '
        'Analytic checks include finite and exactly zero escape, symmetric diffusion with 10^-18 boundary speeds, and agreement '
        'with direct solutions when well-conditioned. Primary tracer mesh is 1601 cells with exponentially fitted conservative '
        'interior fluxes; 801/3201-cell comparisons are evaluated at the original 201 native cell centers. Native background '
        'dynamics and scalar-success/field-tail flags are unchanged. The future-protocol companion is retained separately.')
    output(STEM+'_caption.txt').write_text(caption+'\n')
    audit = dict(generated_at=base.now(),passed=True,selfcheck=checks,checks=audits,
                 background_identity_with_Figure_X=True,future_fields_used=False,new_native_integrations=0,
                 lifetime_display_upper_s=LONG_RESIDENCE,raw_lifetimes_capped=False,
                 membrane_velocity='Frozen trailing-interval instantaneous estimate; zero-membrane control also exported')
    base.write(output(STEM+'_audit.json'),audit)
    html = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>Frozen MCA fates</title>'
            '<style>body{max-width:1250px;margin:32px auto;padding:0 24px;font:16px/1.55 system-ui}img{width:100%;height:auto}</style>'
            '<h1>MCA residence and exit fate with fields frozen at each snapshot</h1><p>Generated '+base.now()+' · Producer finished</p>'
            '<a href="../vectors/'+STEM+'.pdf"><img src="../../'+STEM+'.png" alt="Frozen MCA residence and exit probability"></a>'
            '<p><a href="../vectors/'+STEM+'.pdf">PDF</a> · <a href="../vectors/'+STEM+'.svg">SVG</a> · '
            '<a href="../data/'+STEM+'_fields.npz">Fields and uncapped times</a> · <a href="../data/'+STEM+'_measurements.csv">Measurements</a> · '
            '<a href="../metadata/'+STEM+'_audit.json">Numerical audit</a> · <a href="MCA_molecular_fates.html">Earlier future-protocol calculation</a></p>'
            '<p>'+caption+'</p></html>')
    output(STEM+'.html').write_text(html)
    sources = dict(m['source_sha256']);sources[str(Path(__file__).relative_to(ROOT))] = base.sha(Path(__file__))
    sources[str(Path(base.__file__).relative_to(ROOT))] = base.sha(Path(base.__file__))
    for slot in range(3):
        for ext in ['.npz','.json']:
            p = RUN/'workers'/f'{slot:03d}{ext}';sources[str(p.relative_to(ROOT))] = base.sha(p)
    paths = [output(STEM+ext) for ext in ['.png','.pdf','.svg','.html','_fields.npz','_measurements.csv','_caption.txt','_audit.json']]
    base.write(output(STEM+'.manifest.json'),dict(artifact_root='../..',generated_at=base.now(),producer='finished',passed=True,
               legend=caption,description=caption,sha256={str(p.relative_to(OUT)):base.sha(p) for p in paths},source_sha256=sources,
               audit=audit,configurations=CONFIGS,future_fields_used=False,
               array_definitions=dict(channels=['mean_remaining_residence_s','p_left_first_exit','p_right_first_exit'],
                                      conditioned='time × interleaved [bound_i,unbound_i] × channels',
                                      mixed='time × original 201 native cell centers × channels; initial local pool mixture',
                                      zero_membrane_mixed='Freeze snapshot cortex with membrane velocity forced to zero',
                                      infinite_residence='Exact zero escape; probabilities undefined (NaN)'),
               reproduction=[base.PYTHON+' '+str(Path(__file__))+' render']))
    base.write(RUN/'status.json',dict(generated_at=base.now(),producer='finished',complete=True,audit_passed=True))
    print(json.dumps(dict(output=str(output(STEM+'.png')),checks=audits,passed=True),indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command',choices=['check','launch','worker','render'])
    parser.add_argument('--slot',type=int)
    args = parser.parse_args()
    if args.command == 'check':print(json.dumps(selfcheck(),indent=2))
    elif args.command == 'launch':launch()
    elif args.command == 'worker':worker(args.slot)
    else:render()
