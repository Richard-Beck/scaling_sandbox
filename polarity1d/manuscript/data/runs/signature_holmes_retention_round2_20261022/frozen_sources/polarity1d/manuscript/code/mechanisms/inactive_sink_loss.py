"""Localized inactive-pool removal from exact basal checkpoints; no cues/shams."""
from pathlib import Path
import argparse, json, os, shlex, shutil, subprocess, sys, traceback
import numpy as np

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))
import three_test_reversal as t
n, r = t.n, t.r
from polarity1d.solvers import simulate_reaction_diffusion, laplacian_no_flux

OUT = n.DATA / 'runs/inactive_sink_loss_20261008'
FIG = n.ROOT / 'figures/inactive_sink_loss'
STRENGTHS = np.logspace(-3, 3, 13).tolist()
HORIZON = 20000.


def overlap_weights(left, right, length, cells):
    dx = length / cells
    edges = np.arange(cells) * dx
    return np.maximum(0., np.minimum(edges + dx, right) - np.maximum(edges, left))


def pole_edge(u, length):
    x = (np.arange(len(u)) + .5) * length / len(u)
    threshold = u.min() + .5 * np.ptp(u)
    peak = int(np.argmax(u))
    below = np.flatnonzero(u[peak:] < threshold)
    if not len(below): return float(length)
    j = peak + int(below[0])
    return float(np.interp(threshold, [u[j], u[j-1]], [x[j], x[j-1]]))


class SinkModel:
    """Third field records removed material without returning it to chemistry."""
    state_names = ('active', 'inactive', 'removed')
    observable_name = 'active'
    nonnegative = True
    conserved_groups = ((0, 1, 2),)

    def __init__(self, key, strength, left, right, length, cells):
        self.base = n.make_model(key)
        self.name = self.base.name + '_inactive_sink'
        self.preferred_step = self.base.preferred_step
        self.diffusivities = (*self.base.diffusivities, 0.)
        self.strength = strength
        self.mask = overlap_weights(left, right, length, cells) / (length / cells)

    def reaction(self, time, x, state, stimulus):
        native = self.base.reaction(time, x, state[:2], stimulus)
        removal = self.strength * self.mask * state[1]
        return np.vstack((native[0], native[1] - removal, removal))

    def metadata(self):
        return dict(self.base.metadata(), sink_strength=self.strength,
                    sink_cell_fractions=self.mask.tolist(), removed_field='inert bookkeeping')


def submit(out, action, array=None, dependency=None):
    cmd = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ' + shlex.join(
        [sys.executable, str(Path(__file__)), action, '--out', str(out)])
    if array: cmd += ' --slot "$SLURM_ARRAY_TASK_ID"'
    opts = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
            '--mem=3G', '--time=02:00:00', '--job-name=inactive-sink-' + action,
            '--chdir=' + str(n.ROOT), '--output=' + str(out/'logs'/f'{action}_%A_%a.out'),
            '--error=' + str(out/'logs'/f'{action}_%A_%a.err')]
    if array: opts += ['--array=' + array]
    if dependency: opts += ['--dependency=afterany:' + dependency]
    return subprocess.check_output(opts + ['--wrap', cmd], text=True).strip().split(';')[0]


def initialize(out):
    if out.exists(): raise ValueError('Choose a fresh output directory')
    assert r.read(t.DEFAULT/'results/status.json')['complete']
    cases = sorted([c for c in r.read(t.DEFAULT/'manifest.json')['cases']
                    if c['model'] in ['wave_pinning', 'goryachev'] and c['condition']=='fixed_relative'],
                   key=lambda c: (c['model'], c['length']))
    jobs, exclusions = [], []
    for c in cases:
        if not c['eligible']:
            exclusions.append(dict(model=c['model'], length=c['length'], reason='Preparation ineligible'))
            continue
        with np.load(n.ROOT/c['checkpoint']) as z: edge = pole_edge(z['state'][0], c['length'])
        for geometry in ['opposite', 'fixed_gap']:
            left = c['length'] - 1. if geometry=='opposite' else edge + 3.
            right = left + 1.
            if right > c['length'] + 1e-10 or left < edge:
                exclusions.append(dict(model=c['model'], length=c['length'], geometry=geometry,
                                       reason='Full width-1 sink beyond pole does not fit'))
                continue
            jobs.append(dict(slot=len(jobs), case=c, geometry=geometry, pole_edge=edge,
                             sink_left=left, sink_right=right))
    for folder in ['workers', 'logs', 'results', 'frozen_sources']:
        (out/folder).mkdir(parents=True, exist_ok=True)
    manifest = dict(created_at=n.now(), jobs=jobs, exclusions=exclusions,
        strengths=STRENGTHS, horizon=HORIZON, sink_width=1., fixed_gap=3.,
        definition='Persistent sink: inactive removal k*v within fixed width 1. Opposite sink [L-1,L]. Fixed-gap sink starts 3 native length units beyond initial half-height incumbent edge. Both positions remain fixed. Native reactions/diffusion; zero added stimulus. Removed material is inert bookkeeping with no feedback.',
        measurements='Active amount and peak height in fixed original half-height core, relative to their own initial values. Loss rate = 0.2/(t70-t90), the mean initial-fraction loss rate between 10% and 30% active-amount loss. t90 is loss delay; t50 is half-loss time. Crossings interpolate only adjacent recorded samples; right censor at 20000. Width is the contiguous left domain above the fixed initial half-height threshold.',
        solver=dict(backend='bdf', rtol=2e-8, atol=1e-11, max_step=100.),
        sampling='0 plus 1201 logarithmic samples from 1e-6 to 20000 native time units.',
        no_sink_controls=False, added_cues=False)
    r.write(out/'manifest.json', manifest)
    sources = n.sources()
    for p in [Path(__file__), CODE/'three_test_reversal.py', t.DEFAULT/'manifest.json']:
        sources[str(p.relative_to(n.ROOT))] = n.sha(p)
    for name in sources:
        dest = out/'frozen_sources'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(n.ROOT/name, dest)
    r.write(out/'provenance.json', dict(sources=sources, manifest_sha256=n.sha(out/'manifest.json')))
    r.write(out/'results/status.json', dict(generated_at=n.now(), complete=False, producer='running',
                                           jobs=len(jobs), trials=len(jobs)*len(STRENGTHS)))
    (out/'README.md').write_text('# Inactive-pool sink experiment\n\n'+manifest['definition']+'\n\n'+manifest['measurements']+'\n\nNo no-sink simulations and no added cues. All positive strengths use exact retained basal checkpoints. Original evidence is immutable.\n')
    array = submit(out, 'worker', f'0-{len(jobs)-1}')
    final = submit(out, 'finish', dependency=array)
    r.write(out/'launch.json', dict(array=array, finalizer=final, generated_at=n.now()))
    print(json.dumps(dict(jobs=len(jobs), trials=len(jobs)*len(STRENGTHS), array=array, finalizer=final)))


def verify(out):
    p = r.read(out/'provenance.json')
    assert n.sha(out/'manifest.json') == p['manifest_sha256']
    for name, digest in p['sources'].items():
        assert n.sha(n.ROOT/name)==digest, 'Changed source: '+name
        assert n.sha(out/'frozen_sources'/name)==digest
    return r.read(out/'manifest.json')


def crossing(time, value, level):
    hits = np.flatnonzero(value <= level)
    if not len(hits): return None
    j = int(hits[0])
    if j==0: return float(time[0])
    return float(time[j-1] + (level-value[j-1]) / (value[j]-value[j-1]) * (time[j]-time[j-1]))


def measurements(z, job, strength):
    u, v, time = z['active'], z['inactive'], z['time']
    L = job['case']['length']; dx = L/u.shape[1]
    w = overlap_weights(0., job['pole_edge'], L, u.shape[1])
    mask = w > 0
    amount = u @ w
    peak = np.max(u[:, mask], axis=1)
    norm = amount/amount[0]
    t90, t70, t50 = [crossing(time, norm, q) for q in [.9, .7, .5]]
    height50 = crossing(time, peak/peak[0], .5)
    threshold = u[0].min() + .5*np.ptp(u[0])
    widths = []
    x = (np.arange(u.shape[1])+.5)*dx
    for field in u:
        hits = np.flatnonzero(field < threshold)
        j = int(hits[0]) if len(hits) else len(field)
        widths.append(0. if j==0 else L if j==len(field) else float(np.interp(threshold,[field[j],field[j-1]],[x[j],x[j-1]])))
    widths = np.asarray(widths)
    base = n.make_model(job['case']['model'])
    reaction = base.reaction(0., np.tile(x, len(time)),
        np.vstack((u.ravel(), v.ravel())), np.zeros(u.size))[0].reshape(u.shape)
    rate = (reaction + base.diffusivities[0]*laplacian_no_flux(u, dx)) @ w
    sw = overlap_weights(job['sink_left'], job['sink_right'], L, len(x))
    flux = strength*(v@sw)
    removed = z['removed'].sum(axis=1)*dx
    row = dict(model=job['case']['model'], length=L, geometry=job['geometry'], strength=strength,
        pole_edge=job['pole_edge'], sink_left=job['sink_left'], sink_right=job['sink_right'],
        distance_from_pole_edge=job['sink_left']-job['pole_edge'], initial_active_amount=float(amount[0]),
        initial_peak_height=float(peak[0]), initial_inactive_concentration=float(v[0].mean()),
        initial_removal_flux=float(flux[0]), final_removed_amount=float(removed[-1]),
        t10_loss=t90, t30_loss=t70, t50_loss=t50, peak_t50_loss=height50,
        fractional_loss_rate_10_to_30=(.2/(t70-t90) if t70 is not None and t90 is not None and t70>t90 else None),
        final_active_fraction=float(norm[-1]),
        peak_fraction_at_amount_half=(float(np.interp(t50,time,peak/peak[0])) if t50 is not None else None),
        width_fraction_at_amount_half=(float(np.interp(t50,time,widths/widths[0])) if t50 is not None else None),
        final_peak_fraction=float(peak[-1]/peak[0]),
        censor_horizon=HORIZON)
    arrays = dict(core_amount=amount, peak_height=peak, width=widths, core_rate=rate,
                  removal_flux=flux, removed_total=removed)
    return row, arrays


def worker(out, slot):
    m = verify(out); job = m['jobs'][slot]; c = job['case']; L=c['length']
    assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
    with np.load(n.ROOT/c['checkpoint']) as z: original=z['state'].copy()
    directory = out/'workers'/f'{slot:03d}'; directory.mkdir(exist_ok=True)
    times = np.r_[0., np.geomspace(1e-6, HORIZON, 1201)]
    for i, k in enumerate(m['strengths']):
        model = SinkModel(c['model'], k, job['sink_left'], job['sink_right'], L, original.shape[1])
        state = np.vstack((original, np.zeros_like(original[0])))
        tr = simulate_reaction_diffusion(model, initial_state=state, protocol=n.ConstantInput(0.),
            length=L, cells=original.shape[1], times=times, **m['solver'])
        arrays = dict(time=tr.time, x=tr.x, original=original, **tr.fields)
        row, measured = measurements(arrays, job, k); arrays.update(measured)
        total = (arrays['active']+arrays['inactive']+arrays['removed']).sum(axis=1)
        error = float(np.max(abs(total/total[0]-1.)))
        assert error<1e-7 and min(arrays[a].min() for a in model.state_names)>-1e-8
        path = directory/f'k_{i:02d}.npz'; np.savez_compressed(path, **arrays)
        r.write(path.with_suffix('.json'), dict(generated_at=n.now(), job=job, strength=k,
            metrics=row, solver=tr.metadata, mass_balance_error=error, archive_sha256=n.sha(path)))
        print(json.dumps(dict(slot=slot, k=k, t50=row['t50_loss'], seconds=tr.metadata['wall_seconds'])), flush=True)
    r.write(directory/'complete.json', dict(generated_at=n.now(), trials=len(m['strengths'])))


def finish(out):
    m=verify(out); rows=[]; checks=[]
    for job in m['jobs']:
        directory=out/'workers'/f'{job["slot"]:03d}'
        assert r.read(directory/'complete.json')['trials']==len(m['strengths'])
        c=job['case']; assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
        with np.load(n.ROOT/c['checkpoint']) as cp: original=cp['state'].copy()
        for i,k in enumerate(m['strengths']):
            path=directory/f'k_{i:02d}.npz'; meta=r.read(path.with_suffix('.json'))
            assert meta['job']==job and meta['strength']==k and n.sha(path)==meta['archive_sha256']
            with np.load(path) as z:
                np.testing.assert_array_equal(z['original'],original)
                np.testing.assert_array_equal(z['active'][0],original[0]); np.testing.assert_array_equal(z['inactive'][0],original[1])
                assert np.isfinite(z['active']).all() and np.isfinite(z['inactive']).all()
                total=(z['active']+z['inactive']+z['removed']).sum(axis=1)
                balance=float(np.max(abs(total/total[0]-1.)))
                assert balance<1e-7
                assert min(z[a].min() for a in ['active','inactive','removed'])>-1e-8
                assert np.min(np.diff(z['removed_total']))>-1e-7
                row, arrays=measurements(z,job,k)
                for name,val in arrays.items(): np.testing.assert_allclose(z[name],val,rtol=1e-12,atol=1e-12)
                assert row==meta['metrics']
                # Independent spatial integral and global loss identity.
                dx=c['length']/len(z['x'])
                np.testing.assert_allclose((z['active']+z['inactive']).sum(axis=1)*dx+z['removed_total'],
                    original.sum()*dx, rtol=1e-7,atol=1e-8)
                row['archive']=str(path.relative_to(out)); rows.append(row)
                checks.append(dict(archive=row['archive'],passed=True,mass_error=balance))
    n.csvwrite(out/'results/measurements.csv',rows)
    r.write(out/'results/audit.json',dict(generated_at=n.now(),passed=True,trials=len(rows),checks=checks,
        checks_description='Exact checkpoint identity/hashes; frozen source hashes; complete positive-strength coverage; finite nonnegative fields; native-plus-removed material balance; monotone removed amount; recomputed spatial observables and first-crossing rates.'))
    render(out,m,rows)
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=True,producer='finished',
        jobs=len(m['jobs']),trials=len(rows),audit_passed=True,
        censored_half_loss=sum(row['t50_loss'] is None for row in rows)))
    print(json.dumps(dict(complete=True,trials=len(rows))))


def render(out,m,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.backends.backend_pdf import PdfPages
    figdir = FIG if out==OUT.resolve() else out/'figures'
    figdir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    cmap=plt.get_cmap('viridis'); norm=Normalize(3,39)
    columns=[(model,geometry) for model in ['wave_pinning','goryachev'] for geometry in ['opposite','fixed_gap']]
    titles={'wave_pinning':'Wave-pinning','goryachev':'Goryachev'}
    with PdfPages(figdir/'inactive_sink_loss.pdf') as pdf:
        fig,axes=plt.subplots(3,4,figsize=(17,11),layout='constrained')
        specs=[('fractional_loss_rate_10_to_30','Active-amount loss rate\n(initial fraction / native time)'),
               ('t10_loss','Time to 10% active-amount loss'),('t50_loss','Time to 50% active-amount loss')]
        for col,(model,geo) in enumerate(columns):
            rs=[v for v in rows if v['model']==model and v['geometry']==geo]
            for i,(key,label) in enumerate(specs):
                ax=axes[i,col]
                for L in sorted(set(v['length'] for v in rs)):
                    series=sorted([v for v in rs if v['length']==L],key=lambda v:v['strength'])
                    x=[v['strength'] for v in series]; y=[v[key] if v[key] is not None else np.nan for v in series]
                    ax.plot(x,y,'o-',ms=2.6,lw=1,color=cmap(norm(L)))
                    if i>0:
                        missing=[v['strength'] for v in series if v[key] is None]
                        ax.scatter(missing,[HORIZON]*len(missing),marker='^',s=15,facecolors='none',edgecolors=[cmap(norm(L))])
                ax.set(xscale='log',yscale='log',xlabel='Sink strength k (1 / native time)')
                ax.grid(alpha=.15)
                if i>0: ax.axhline(HORIZON,color='.6',ls=':',lw=.8)
                if col==0: ax.set_ylabel(label)
            axes[0,col].set_title(titles[model]+'\n'+('Opposite end' if geo=='opposite' else 'Fixed gap from initial pole'))
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,label='Cell length L',shrink=.8)
        fig.suptitle('Inactive-pool sinks: incumbent loss rate and delay',fontsize=18)
        fig.supxlabel('Width 1; fixed gap 3 from initial half-height edge. Continuous removal; no added cues.\nRate measured over 10–30% loss. Open triangles: threshold not reached by 20,000; missing rates remain gaps.',fontsize=10)
        fig.savefig(figdir/'inactive_sink_loss.png',dpi=180);fig.savefig(figdir/'inactive_sink_loss.svg');pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(17,7),layout='constrained')
        for col,(model,geo) in enumerate(columns):
            for i,(key,label) in enumerate([('peak_fraction_at_amount_half','Peak height / initial height\nat 50% active-amount loss'),('width_fraction_at_amount_half','Domain width / initial width\nat 50% active-amount loss')]):
                ax=axes[i,col]
                for L in sorted(set(v['length'] for v in rows if v['model']==model and v['geometry']==geo)):
                    series=sorted([v for v in rows if v['model']==model and v['geometry']==geo and v['length']==L],key=lambda v:v['strength'])
                    ax.plot([v['strength'] for v in series],[v[key] if v[key] is not None else np.nan for v in series],'o-',ms=3,lw=1,color=cmap(norm(L)))
                ax.set(xscale='log',xlabel='Sink strength k',ylim=(0,1.1));ax.grid(alpha=.15)
                if col==0:ax.set_ylabel(label)
            axes[0,col].set_title(titles[model]+'\n'+('Opposite end' if geo=='opposite' else 'Fixed gap'))
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,label='Cell length L',shrink=.8)
        fig.suptitle('How the incumbent erodes: height and occupied width',fontsize=17)
        fig.supxlabel('Active amount and peak height use the fixed initial half-height core. Width uses the fixed INITIAL half-height threshold.\nGaps indicate 50% amount loss was not reached; width can reach zero before active amount vanishes.',fontsize=10)
        fig.savefig(figdir/'incumbent_morphology.png',dpi=180);pdf.savefig(fig);plt.close(fig)
        # Representative dynamics: k=1, lengths 9/21/39, same strength in each panel.
        fig,axes=plt.subplots(3,4,figsize=(17,10),layout='constrained')
        for col,(model,geo) in enumerate(columns):
            for L in [9.,21.,39.]:
                row=next(v for v in rows if v['model']==model and v['geometry']==geo and v['length']==L and v['strength']==1.)
                with np.load(out/row['archive']) as z:
                    for i,(val,label) in enumerate([(z['core_amount']/z['core_amount'][0],'Active amount / initial amount'),(z['peak_height']/z['peak_height'][0],'Peak height / initial height'),(z['width']/z['width'][0],'Domain width / initial width')]):
                        axes[i,col].plot(z['time'],val,label=f'L={L:g}',color=cmap(norm(L)),lw=1.7)
                        axes[i,col].set(xscale='log',xlim=(.001,HORIZON),xlabel='Time since sink onset',ylim=(0,1.05));axes[i,col].grid(alpha=.15)
                        if col==0:axes[i,col].set_ylabel(label)
            axes[0,col].set_title(titles[model]+'\n'+('Opposite end' if geo=='opposite' else 'Fixed gap'));axes[0,col].legend(frameon=False)
        fig.suptitle('Incumbent trajectories at sink strength k = 1',fontsize=17)
        fig.savefig(figdir/'incumbent_trajectories.png',dpi=180);pdf.savefig(fig);plt.close(fig)
    stamp=n.now(); csvlink=os.path.relpath(out/'results/measurements.csv',figdir)
    (figdir/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Inactive sink and incumbent loss</title></head><body style="font:17px system-ui;max-width:1600px;margin:30px auto;padding:0 24px"><h1>Inactive-pool sinks and incumbent loss</h1><p>Generated '+stamp+'; producer finished. '+str(len(rows))+' sink-only trials; checkpoint, source and material-budget audits passed.</p><p>'+m['definition']+'</p><p>'+m['measurements']+'</p><p><a href="inactive_sink_loss.pdf">Three-page PDF</a> · <a href="'+csvlink+'">Measurements CSV</a></p>'+''.join('<img style="width:100%" src="'+name+'">' for name in ['inactive_sink_loss.png','incumbent_morphology.png','incumbent_trajectories.png'])+'<p>Positive strengths only. No no-sink simulations or reversal cues. Sink removal reduces available active+inactive material. The inert removed field checks the budget and never feeds back. Fixed gap is measured from the initial pole edge, not its peak or moving edge. Oversized or overlapping sink geometries are excluded explicitly in the manifest. This diagnostic measures incumbent erosion, not reversal.</p></body></html>')


def self_check():
    # Algebra/geometry checks only; no no-sink simulations.
    for key in ['wave_pinning','goryachev']:
        L=39.;cells=201;x=(np.arange(cells)+.5)*L/cells
        model=SinkModel(key,2.,7.123,8.123,L,cells)
        state=np.vstack((np.full(cells,1.),np.full(cells,.2),np.zeros(cells)))
        f=model.reaction(0.,x,state,np.zeros(cells));native=model.base.reaction(0.,x,state[:2],np.zeros(cells))
        np.testing.assert_array_equal(f[0],native[0])
        np.testing.assert_allclose(f.sum(axis=0),0.,atol=1e-15)
        np.testing.assert_allclose(model.mask.sum()*L/cells,1.,atol=1e-14)
        np.testing.assert_allclose(f[2].sum()*L/cells,.4,atol=1e-14)
    assert crossing(np.array([0.,1.,2.]),np.array([1.,.8,.6]),.9)==.5
    assert crossing(np.array([0.,1.]),np.array([1.,.8]),.5) is None
    print('Sink-only algebra, unit physical width, removal budget and threshold interpolation checks passed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker','finish','self-check'])
    p.add_argument('--out',type=Path,default=OUT);p.add_argument('--slot',type=int)
    a=p.parse_args();out=a.out.resolve()
    try:
        if a.action=='self-check':self_check()
        elif a.action=='init':initialize(out)
        elif a.action=='worker':worker(out,a.slot)
        else:finish(out)
    except Exception:
        if (out/'results').exists():r.write(out/'results'/f'{a.action}_{a.slot}.error.json',dict(generated_at=n.now(),error=traceback.format_exc()))
        raise
