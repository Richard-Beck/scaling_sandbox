"""Figure 3 round 3: retained half-heading readouts and four frozen relaxations.

freeze/worker are the only integration workflow. collect and render reuse fields;
they never launch a simulation or replace a censoring bound by an observation.
"""
from pathlib import Path
import argparse
import copy
import json
import os
import shlex
import shutil
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True
ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MAN = ROOT/'polarity1d/manuscript'
PRIMARY = MAN/'data/runs/signature_figure3_campaign_20261020'
# Identical executable native model, solver and parameters as the retained A1.
FROZEN = PRIMARY/'frozen_sources/code'
sys.path.insert(0, str(FROZEN/'mechanisms'))
sys.path.insert(0, str(FROZEN))
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import signature_figure3_native_20261020 as native
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_regions_20261020 as regions

DEV = MAN/'figures/length_mechanisms_20261008/dev'
SOURCE = DEV/'revision_program_20261020/figure3_distance/supporting/data'
OUT = DEV/'revision_round3_20261022/figure3'
RUN = MAN/'data/runs/figure3_round3_relaxation_v2_20261022'
NAMES = dict(wave_pinning='Wave-pinning', goryachev='Goryachev', holmes_model3='Holmes', spring='Spring')
COLORS = dict(wave_pinning='#2878b5', goryachev='#c83b3b', holmes_model3='#36944c', spring='#b78b16')
STYLES = ('-', '--', ':')
EXAMPLES = [('i', 30., -1.), ('ii', 30., 4.), ('iii', 47.6, -1.), ('iv', 47.6, 4.)]
MAIN_EXAMPLES = EXAMPLES
SINGLE_C_LAYOUT = False
EXTRA_SUPPLEMENT_SOURCES = []
DATA = OUT/'supporting/data'
META = OUT/'supporting/metadata'
write, read, sha = native.write, native.read, native.sha


def source_tables():
    h = pd.read_csv(SOURCE/'A1_heading_target_time_brackets.csv')
    a = pd.read_csv(SOURCE/'A1_complete_region_evidence.csv')
    return h, a


def freeze():
    if (RUN/'manifest.json').exists():
        return read(RUN/'manifest.json')
    h, a = source_tables()
    for d in ('inputs', 'workers', 'logs', 'frozen_sources'):
        (RUN/d).mkdir(parents=True, exist_ok=True)
    jobs = []
    for label, size, gap in EXAMPLES:
        row = a[(a.model == 'wave_pinning') & (a.size_parameter == size) &
                np.isclose(a.gap, gap) & a.control.isna()].iloc[0]
        status = read(Path(row.native_status))
        p = status['placement']
        original = PRIMARY/p['input']
        assert sha(original) == p['input_sha256']
        dest = RUN/'inputs'/original.name
        shutil.copy2(original, dest)
        jobs.append(dict(label=label, slot=int(row.slot), placement=p,
            source_status=str(row.native_status), source_status_sha256=sha(Path(row.native_status)),
            input=str(dest.relative_to(RUN)), input_sha256=sha(dest)))
    script = RUN/'frozen_sources'/Path(__file__).name
    shutil.copy2(__file__, script)
    # These readout modules are absent from the original producer freeze.
    # Preserve their current detector code without modifying that producer.
    for module in (regions,):
        shutil.copy2(module.__file__, RUN/'frozen_sources'/Path(module.__file__).name)
    original = read(PRIMARY/'manifest.json')
    plan = dict(generated_at=native.now(), jobs=jobs, native_protocols=4, retuned=False,
        native_source_sha256=original['source_sha256'], primary=str(PRIMARY),
        script_sha256=sha(script), cue_on_horizon=20000., cue_off_fraction=.95,
        relaxation_duration_factor=4., solver=original['solver'],
        stationarity=dict(active_scale=2.2683, window=1000., scaled_tolerance=1e-4,
            active_derivative_tolerance=2.2683e-7,
            source='signature_figure3_basals_20261020.py:213-215; same scaled RHS gate, restricted to native active species as requested'),
        readout_source=str(SOURCE/'A1_heading_target_time_brackets.csv'),
        readout_sha256=sha(SOURCE/'A1_heading_target_time_brackets.csv'),
        reason='Exact 95% event states are not archived. Replay four identical cue-on phases from exact basal states; stop at event or 20,000, then relax.',
        no_new_readout_campaign=True)
    write(RUN/'manifest.json', plan)
    command = shlex.join([native.PYTHON, '-B', str(script), 'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
        '--mem=3G', '--time=01:00:00', '--signal=B:USR1@60', '--array=0-3',
        '--job-name=fig3-r3-relax', '--chdir='+str(ROOT),
        '--output='+str(RUN/'logs/%A_%a.out'), '--error='+str(RUN/'logs/%A_%a.err'), '--wrap', command]
    write(RUN/'submission.json', dict(sbatch_arguments=args, no_concurrency_throttle=True,
        reason='Four independent serial BDF protocols; small QoS; original native solver settings.'))
    return dict(manifest=str(RUN/'manifest.json'), submission=args)


def worker(slot):
    native.require_slurm()
    plan = read(RUN/'manifest.json')
    assert sha(Path(__file__)) == plan['script_sha256']
    for rel, digest in plan['native_source_sha256'].items():
        assert sha(PRIMARY/'frozen_sources'/rel) == digest
    job = plan['jobs'][slot]
    assert sha(RUN/job['input']) == job['input_sha256']
    p = job['placement']
    ad = native.Adapter('wave_pinning', p['size_parameter'], p['cells'])
    with np.load(RUN/job['input']) as z:
        initial = z['state'].copy()
    phi0 = np.angle(np.sum(ad.primary(initial)*np.exp(2j*np.pi*ad.q)))
    target = float(geometry.arc(2*np.pi*p['cue_center_fraction'], phi0, 2*np.pi))
    def off_event(y):
        z = np.sum(ad.primary(y)*np.exp(2j*np.pi*ad.q))
        return float(np.angle(z*np.exp(-1j*(phi0+.95*target))))
    directory = RUN/'workers'/job['label']
    directory.mkdir(parents=True, exist_ok=True)
    status = dict(generated_at=native.now(), producer='running', job=job, complete=False,
                  slurm_job_id=os.environ['SLURM_JOB_ID'])
    write(directory/'status.json', status)
    cue = native.integrate(ad, initial, 20000., p['cue_peak'], p['cue_center_fraction'],
        event=off_event, solver_options=plan['solver'], deadline=time.monotonic()+3200)
    on = native.save_phase(directory/'cue_on.npz', cue, ad, p['cue_peak'], p['cue_center_fraction'])
    assert cue['status'] in ('complete', 'event'), cue['status']
    duration = float(cue['time'][-1])
    tolerance = plan['stationarity']['active_derivative_tolerance']
    def stationary(y):
        return tolerance-float(np.max(abs(ad.rhs(y)[:ad.cells])))
    off = native.integrate(ad, cue['states'][-1], 4*duration,
        event=stationary, solver_options=plan['solver'], deadline=time.monotonic()+3200)
    tail = native.save_phase(directory/'cue_off.npz', off, ad)
    assert off['status'] in ('complete', 'event', 'event_at_start'), off['status']
    final_rate = float(np.max(abs(ad.rhs(off['states'][-1])[:ad.cells])))
    status.update(generated_at=native.now(), producer='finished', complete=True,
        cue_on=on, cue_off=tail, cue_on_duration=duration,
        withdrawal_condition='95_percent_heading' if cue['event_time'] is not None else '20000_horizon',
        relaxation_duration=float(off['time'][-1]),
        stop_condition='stationarity' if off['event_time'] is not None else '4x_cue_on_duration',
        stationarity_tolerance=tolerance, final_max_active_derivative=final_rate,
        target_angle=float(target), final_cue_on_heading_fraction=float(np.angle(
            np.sum(ad.primary(cue['states'][-1])*np.exp(2j*np.pi*ad.q))*np.exp(-1j*phi0))/target),
        exact_withdrawal_restart=bool(np.array_equal(cue['states'][-1], off['states'][0])))
    write(directory/'status.json', status)
    print(json.dumps({k:status[k] for k in ('complete','cue_on_duration','withdrawal_condition',
        'relaxation_duration','stop_condition','final_max_active_derivative')}))


def full_history(status):
    times, fields, paths = [], [], []
    for phase in status['phases']:
        path = PRIMARY/phase['path']
        assert sha(path) == phase['sha256']
        with np.load(path) as z:
            t, y = z['time'].copy()+phase['absolute_start'], z['states'].copy()
        use = slice(1, None) if times else slice(None)
        times.extend(t[use]); fields.extend(y[use]); paths.append(str(path))
    return np.asarray(times), np.asarray(fields), paths


def crossing_state(t, y, cells, target, phi0):
    """Linear saved-field interpolant, only within an existing valid bracket."""
    phase = np.exp(2j*np.pi*(np.arange(cells)+.5)/cells)
    v = y[:, :cells]@phase
    angle = np.angle(v*np.exp(-1j*phi0))
    idx = np.flatnonzero((angle[:-1] < target) & (angle[1:] >= target))[0]
    f = brentq(lambda x: np.angle(((1-x)*v[idx]+x*v[idx+1])*np.exp(-1j*phi0))-target, 0, 1)
    return float((1-f)*t[idx]+f*t[idx+1]), (1-f)*y[idx]+f*y[idx+1], int(idx), float(f)


def collect():
    DATA.mkdir(parents=True, exist_ok=True); META.mkdir(parents=True, exist_ok=True)
    h, a = source_tables()
    selected = h[(h.target_fraction == .5) & (h.coherence_floor == .05) & h.control.isna()]
    rows, provenance, profile_rows, profile_meta = [], [], [], []
    for row in selected.itertuples():
        record = dict(slot=int(row.slot), model=row.model, size_parameter=row.size_parameter,
            gap=0. if abs(row.gap) < 1e-9 else row.gap, status=row.status,
            target_angle=row.target_angle*2, half_angle=row.target_angle,
            lower=None if pd.isna(row.lower) else row.lower,
            upper=None if pd.isna(row.upper) else row.upper,
            t50=None, speed=None, speed_lower=None, speed_upper=None,
            heading_undefined='gap' in row.status, no_turn=row.target_angle <= 1e-8,
            topology='crossing_unavailable', saddle=None)
        source = a[a.slot == row.slot].iloc[0]
        assert sha(Path(source.native_status)) == source.native_status_sha256
        status = read(Path(source.native_status)); p = status['placement']
        if row.status == 'observed_saved_valid_segment_crossing' and not row.prior_heading_gap:
            t, y, paths = full_history(status)
            ad = native.Adapter(row.model, row.size_parameter, p['cells'])
            phi0 = np.angle(np.sum(ad.primary(y[0])*np.exp(2j*np.pi*ad.q)))
            found = (t >= row.lower-1e-9) & (t <= row.upper+1e-9)
            tb, yb = t[found], y[found]
            assert len(tb) >= 2
            tc, yc, index, f = crossing_state(tb, yb, ad.cells, row.target_angle, phi0)
            old = p['initial_10_percent_contour']['origin']/p['basal_circumference']
            topo = regions.readout(ad.primary(yc), ad.circumference(yc), old, p['cue_center_fraction'])
            end_topo = [regions.readout(ad.primary(z), ad.circumference(z), old,
                                       p['cue_center_fraction'], detailed=False)['two_region_available'] for z in yb]
            record.update(t50=tc, speed=row.target_angle/tc,
                speed_lower=row.target_angle/row.upper, speed_upper=row.target_angle/row.lower,
                topology=topo['state_class'], saddle=topo['two_region_available'],
                topology_bracket_consistent=len(set(end_topo+[topo['two_region_available']])) == 1,
                interpolation_fraction=f, time_error_width=row.upper-row.lower)
            provenance.append(dict(slot=int(row.slot), source_status=str(source.native_status),
                source_status_sha256=sha(Path(source.native_status)), native_archives=paths,
                interpolated_topology=topo, endpoint_two_region=end_topo))
        elif row.status == 'right_censored_saved_heading_at_completed_horizon':
            record.update(speed_upper=row.target_angle/row.lower)
        rows.append(record)
        print('collected', row.model, row.size_parameter, row.slot, flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA/'half_heading_readouts.csv', index=False)
    write(META/'topology_provenance.json', provenance)
    regimes, bands = [], []
    for (model, size), g in frame.groupby(['model', 'size_parameter']):
        g = g[~g.no_turn].sort_values('gap')
        saddles = g[g.saddle == True]
        regimes.append(dict(model=model, size_parameter=size,
            first_saddle_gap=None if saddles.empty else float(saddles.iloc[0].gap),
            first_saddle_slot=None if saddles.empty else int(saddles.iloc[0].slot),
            crossings_scored=int(g.t50.notna().sum()),
            censored_or_undefined_crossings=int(g.t50.isna().sum())))
        x = g.gap.to_numpy()
        for j, r in enumerate(g.itertuples()):
            if r.heading_undefined:
                left = float((x[j-1]+x[j])/2) if j else float(x[j])
                right = float((x[j]+x[j+1])/2) if j+1 < len(x) else float(x[j])
                bands.append(dict(model=model, size_parameter=size, slot=r.slot,
                    left=left, right=right, gap=r.gap, status=r.status))
    pd.DataFrame(regimes).to_csv(DATA/'regime_lines.csv', index=False)
    pd.DataFrame(bands).to_csv(DATA/'undefined_heading_bands.csv', index=False)
    # A uses the retained cue-on profiles, with only the final row supplied by
    # the four new protocols. All snapshot reconstruction is explicit.
    for label, size, gap in EXAMPLES:
        r = a[(a.model == 'wave_pinning') & (a.size_parameter == size) & np.isclose(a.gap, gap) & a.control.isna()].iloc[0]
        status = read(Path(r.native_status)); p = status['placement']
        t, y, _ = full_history(status)
        ad = native.Adapter('wave_pinning', size, p['cells'])
        phi0 = np.angle(np.sum(ad.primary(y[0])*np.exp(2j*np.pi*ad.q)))
        target = float(geometry.arc(2*np.pi*p['cue_center_fraction'], phi0, 2*np.pi))
        snapshots = [(0, 0., y[0], 'saved_initial')]
        for fraction in (.5, .9):
            tc, yc, _, _ = crossing_state(t, y, ad.cells, fraction*target, phi0)
            snapshots.append((fraction, tc, yc, 'linear_saved_field_interpolation'))
        relaxed = read(RUN/'workers'/label/'status.json')
        assert relaxed['complete'] and relaxed['exact_withdrawal_restart']
        with np.load(RUN/'workers'/label/'cue_off.npz') as z:
            final = z['states'][-1].copy()
        snapshots.append((1., relaxed['cue_on_duration']+relaxed['relaxation_duration'], final, 'new_native_relaxation'))
        old = p['initial_10_percent_contour']['origin']/p['basal_circumference']
        for stage, (fraction, tc, yc, origin) in enumerate(snapshots):
            v = geometry.activity_vector(ad.primary(yc), ad.q)
            topo = regions.readout(ad.primary(yc), size, old, p['cue_center_fraction'])
            # Show the short-route (inter-region) saddle, rather than marking
            # the opposite resting-floor basin boundary as an active saddle.
            valleys = topo.get('valley_descriptors', []) if topo['two_region_available'] else []
            saddle_positions = [float(geometry.arc((v['cell_coordinate']+.5)/ad.cells,
                                                   old, 1.)*size) for v in valleys]
            active_saddles = [s for s in saddle_positions if 0 < s < p['center_distance']]
            meta = dict(label=label, stage=stage, fraction=fraction, elapsed=tc, origin=origin,
                size=size, gap=gap, target_angle=target, heading=float(geometry.arc(v['heading_radians']/(2*np.pi), old, 1.)*size),
                heading_fraction=float(geometry.arc(v['heading_radians'], phi0, 2*np.pi)/target),
                cue_left=p['cue_near_edge_relative'], cue_right=p['cue_far_edge_relative'],
                saddles=active_saddles, topology=topo['state_class'])
            profile_meta.append(meta)
            for x, u in zip(geometry.arc(ad.q, old, 1.)*size, ad.primary(yc)):
                profile_rows.append(dict(label=label, stage=stage, x=float(x), active=float(u)))
        print('profiles', label, flush=True)
    pd.DataFrame(profile_rows).to_csv(DATA/'activity_profiles.csv.gz', index=False, compression='gzip')
    write(DATA/'activity_profile_metadata.json', profile_meta)
    write(META/'source_sha256.json', {str((SOURCE/name).relative_to(ROOT)):sha(SOURCE/name)
        for name in ('A1_heading_target_time_brackets.csv', 'A1_complete_region_evidence.csv')})
    write(META/'collection.json', dict(generated_at=native.now(), placements=len(frame),
        observed_crossings=int(frame.t50.notna().sum()), heading_undefined=int(frame.heading_undefined.sum()),
        zero_target=int(frame.no_turn.sum()), bracket_topology_disagreements=int((frame.topology_bracket_consistent == False).sum()),
        native_integrations=0, interpolation='Linear native saved-field interpolant at heading root within the source bracket; no interpolation across direction gaps',
        topology='Existing two-region periodic watershed detector; native prominence .05; no placement-class proxy',
        coherence_floor=.05))


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import to_rgb
    from pypdf import PdfReader, PdfWriter
    plt.rcParams.update({'font.size':8.5, 'axes.labelsize':8.5, 'axes.titlesize':9,
        'xtick.labelsize':7.5, 'ytick.labelsize':7.5, 'font.family':'DejaVu Sans',
        'pdf.fonttype':42, 'svg.fonttype':'none', 'axes.spines.top':False, 'axes.spines.right':False})
    def color(model, i):
        c = np.array(to_rgb(COLORS[model])); return tuple(c*(.62+.19*i)+(1-(.62+.19*i)))
    def clean(ax):
        ax.tick_params(length=3, pad=2); ax.grid(False)
    profiles = pd.read_csv(DATA/'activity_profiles.csv.gz')
    metadata = read(DATA/'activity_profile_metadata.json')
    frame = pd.read_csv(DATA/'half_heading_readouts.csv')
    regimes = pd.read_csv(DATA/'regime_lines.csv'); bands = pd.read_csv(DATA/'undefined_heading_bands.csv')
    if SINGLE_C_LAYOUT:
        fig = plt.figure(figsize=(16.2, 8.6))
        agrid = fig.add_gridspec(4, 2, left=.052, right=.348, bottom=.16,
            top=.86, wspace=.28, hspace=.20)
        bgrid = fig.add_gridspec(1, 4, left=.414, right=.972, bottom=.59,
            top=.86, wspace=.42)
        cgrid = fig.add_gridspec(1, 4, left=.414, right=.972, bottom=.16,
            top=.43, wspace=.42)
    else:
        fig = plt.figure(figsize=(11.7, 13.6))
        outer = fig.add_gridspec(3, 1, height_ratios=[6.2, 2.35, 2.35],
            hspace=.34, top=.918, bottom=.112, left=.082, right=.965)
        agrid = outer[0].subgridspec(4, 4, wspace=.27, hspace=.18)
        bgrid=outer[1].subgridspec(1,4,wspace=.38)
        cgrid=outer[2].subgridspec(1,4,wspace=.38)
    a_axes=[]
    for j, (label, size, gap) in enumerate(MAIN_EXAMPLES):
        column=[]
        ymax=profiles[profiles.label==label].active.max()*1.1
        for stage in range(4):
            ax=fig.add_subplot(agrid[stage,j]); clean(ax); column.append(ax)
            m=next(x for x in metadata if x['label']==label and x['stage']==stage)
            p=profiles[(profiles.label==label)&(profiles.stage==stage)].sort_values('x')
            ax.axvspan(m['cue_left'],m['cue_right'],color='.87',zorder=0)
            ax.plot(p.x,p.active,color=color('wave_pinning',0 if size==30 else 2),lw=1.4)
            ax.set_ylim(0,ymax); ax.set_xlim(-size/2,size/2)
            for saddle in m['saddles']:
                u=float(np.interp(saddle,p.x,p.active))
                ax.plot(saddle,u,'D',mfc='white',mec='black',ms=4,zorder=5)
            ax.plot(m['heading'],0,'|',color='black',ms=9,mew=1.5,clip_on=False,zorder=6)
            if stage==0:
                placement_title=(f'({label})  gap {gap:g}' if SINGLE_C_LAYOUT else f'({label})  C = {size:g}, gap {gap:g}')
                ax.set_title(placement_title+'\n'+('straddling' if gap<0 else 'separated'),pad=7)
            if stage<3: ax.tick_params(labelbottom=False)
            else: ax.set_xlabel('Arc position (native length)')
            ax.text(.97,.89,('cue off' if stage==3 else f'{int(m["fraction"]*100)}%'),
                transform=ax.transAxes,ha='right',va='top',fontsize=8,color='black')
            if j==0: ax.set_ylabel('Active species, u')
        a_axes.append(column)
    apos=a_axes[0][0].get_position()
    if SINGLE_C_LAYOUT:
        fig.text(.015,.955,'A',fontsize=13,fontweight='bold',color='black')
        fig.text(.052,.955,'Activity profiles · Wave-pinning, C = 30',fontsize=10,color='black')
    else:
        fig.text(.025,apos.y1+.070,'A',fontsize=13,fontweight='bold',color='black')
        fig.text(.082,apos.y1+.070,'Activity profiles · Wave-pinning',fontsize=10,color='black')
    # Symbols beside their own panel, never in a detached side box.
    handles=[Patch(facecolor='.87',edgecolor='none',label='cue support (fixed; also shown after cue-off)'),
        Line2D([],[],marker='D',mfc='white',mec='black',ls='none',ms=4,label='inter-region saddle'),
        Line2D([],[],marker='|',color='black',ls='none',ms=9,label='heading on arc axis')]
    if SINGLE_C_LAYOUT:
        fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.048,.018),
            ncol=1,frameon=False,fontsize=7.6,handlelength=1.2,borderaxespad=0)
    else:
        fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.077,apos.y1+.050),
            ncol=3,frameon=False,fontsize=7.6,handlelength=1.2,columnspacing=1.8,borderaxespad=0)
    baxes=[]; caxes=[]
    render_records=[]
    for j, model in enumerate(NAMES):
        b=fig.add_subplot(bgrid[0,j]); c=fig.add_subplot(cgrid[0,j],sharex=b)
        baxes.append(b); caxes.append(c)
        sizes=sorted(frame[frame.model==model].size_parameter.unique())
        b.set_title(NAMES[model],color='black',pad=9)
        for ax in (b,c):
            clean(ax); ax.set_yscale('log')
            for r in bands[bands.model==model].itertuples():
                ax.axvspan(r.left,r.right,color='.92',zorder=-3,lw=0)
        for i,size in enumerate(sizes):
            g=frame[(frame.model==model)&(frame.size_parameter==size)&(~frame.no_turn)].sort_values('gap')
            col=color(model,i)
            rg=regimes[(regimes.model==model)&(regimes.size_parameter==size)].iloc[0]
            for ax in (b,c):
                if pd.notna(rg.first_saddle_gap):ax.axvline(rg.first_saddle_gap,color=col,ls=STYLES[i],lw=1.8-.4*i,zorder=-1)
            for ax,key,shape in ((b,'t50','^'),(c,'speed','v')):
                values=g[key].to_numpy()
                ax.plot(g.gap,values,color=col,ls=STYLES[i],lw=1.45,marker='o',ms=2.7,zorder=3)
                censored=g[g.status=='right_censored_saved_heading_at_completed_horizon']
                for r in censored.itertuples():
                    bound=r.lower if key=='t50' else r.speed_upper
                    ax.plot(r.gap,bound,marker=shape,mfc='white',mec=col,ls='none',ms=6,zorder=5)
                valid=g[g[key].notna()]
                if len(valid):
                    r=valid.iloc[-1]
                    # Direct labels use point offsets; manual small separation
                    # retains each label's exact terminal datum as its anchor.
                    dy={0:-9,1:1,2:10}[i]
                    ax.annotate(f'C={size:g}',(r.gap,r[key]),xytext=(5,dy),textcoords='offset points',
                        color=col,fontsize=7.2,va='center',ha='left',clip_on=False)
            if model=='wave_pinning':
                for label, sz, gap in MAIN_EXAMPLES:
                    if sz!=size:continue
                    r=g[np.isclose(g.gap,gap)].iloc[0]
                    for ax,key in ((b,'t50'),(c,'speed')):
                        ax.plot(r.gap,r[key],marker='o',ms=6.5,mfc=col,mec='black',mew=.6,zorder=7)
                        offsets={'i':(-1,-11),'ii':(0,-13),'iii':(-1,9),'iv':(0,10)}
                        ax.annotate(label,(r.gap,r[key]),xytext=offsets[label],textcoords='offset points',
                            fontsize=8,color='black',fontweight='bold',ha='center',zorder=8)
            render_records.append(dict(model=model,size_parameter=float(size),
                observed_points=int(g.t50.notna().sum()), censored_symbols=int(len(censored)),
                gaps_with_no_connected_curve=g[g.t50.isna()].gap.tolist(),
                regime_line=None if pd.isna(rg.first_saddle_gap) else float(rg.first_saddle_gap)))
        b.tick_params(labelbottom=False)
        c.set_xlabel('Cue edge gap'+('\n(Spring native length)' if model=='spring' else '\n(native length)'))
        g=frame[(frame.model==model)&(~frame.no_turn)]
        span=g.gap.max()-g.gap.min()
        b.set_xlim(g.gap.min()-.025*span,g.gap.max()+.13*span)
        # Headroom keeps direct size labels clear of the plot boundaries.
        for ax in (b,c):
            lo,hi=ax.get_ylim(); ax.set_ylim(lo/1.2,hi*1.45)
        if j==0:
            b.set_ylabel(r'Time to 50%, $t_{50}$'+'\n(native time)')
            c.set_ylabel('Mean angular speed\n(rad / native time)')
    for label,axes,title in [('B',baxes,'Time to 50% of cue-direction heading'),
                             ('C',caxes,'Mean angular speed to 50%')]:
        pos=axes[0].get_position()
        label_x, title_x = (.377,.414) if SINGLE_C_LAYOUT else (.025,.082)
        label_y = (.955 if label=='B' else .505) if SINGLE_C_LAYOUT else pos.y1+.032
        fig.text(label_x,label_y,label,fontsize=13,fontweight='bold',color='black')
        fig.text(title_x,label_y,title,fontsize=10,color='black')
    placement_labels='i–ii' if SINGLE_C_LAYOUT else 'i–iv'
    symbols=[Patch(facecolor='.92',edgecolor='none',label='grey band: heading undefined at ≥1 size'),
        Line2D([],[],color='.25',ls=':',label='vertical line: first saddle at 50% (same C style)'),
        Line2D([],[],marker='^',mfc='white',mec='black',ls='none',label='time: lower bound'),
        Line2D([],[],marker='v',mfc='white',mec='black',ls='none',label='speed: upper bound'),
        Line2D([],[],marker='o',mfc=COLORS['wave_pinning'],mec='black',ls='none',label=placement_labels+': profile placements in A')]
    legend_x = .414 if SINGLE_C_LAYOUT else .08
    fig.legend(handles=[symbols[k] for k in (0,1,2,4)],loc='upper left',
        bbox_to_anchor=(legend_x,baxes[0].get_position().y0-(.023 if SINGLE_C_LAYOUT else .008)),ncol=2,
        frameon=False,fontsize=7.2,columnspacing=2,handlelength=1.8,borderaxespad=0)
    c_symbols=[symbols[k] for k in (0,1,3,4)] if SINGLE_C_LAYOUT else symbols
    fig.legend(handles=c_symbols,loc='lower left',bbox_to_anchor=(legend_x,.018),ncol=2 if SINGLE_C_LAYOUT else 3,
        frameon=False,fontsize=7.2 if SINGLE_C_LAYOUT else 7.6,columnspacing=1.5,handlelength=1.8)
    for ext in ('pdf','png','svg'):
        fig.savefig(OUT/f'figure3_round3.{ext}',dpi=220)
    plt.close(fig)
    write(META/'render_coordinates.json',render_records)
    if SINGLE_C_LAYOUT:
        boxes=lambda axes:[list(ax.get_position().bounds) for ax in axes]
        write(META/'layout.json',dict(generated_at=native.now(),figure_size_inches=[16.2,8.6],
            profile_sizes=[30.],profile_labels=[x[0] for x in MAIN_EXAMPLES],
            profile_axes=boxes([ax for column in a_axes for ax in column]),
            time_axes=boxes(baxes),speed_axes=boxes(caxes),native_integrations=0))
    # Add the requested 90/99 timing pages, then retain every preceding page.
    added=OUT/'heading_time_supplement.pdf'
    h,_=source_tables()
    with PdfPages(added) as pdf:
        for fraction in (.9,.99):
            f,axes=plt.subplots(2,2,figsize=(10,7.5))
            f.suptitle(f'Supplement: time to {fraction:.0%} of target heading; continuous cue',fontsize=12)
            for ax,(model,name) in zip(axes.flat,NAMES.items()):
                g=h[(h.model==model)&h.control.isna()&(h.target_fraction==fraction)&(h.coherence_floor==.05)&(h.target_angle>1e-8)]
                for i,(size,sg) in enumerate(g.groupby('size_parameter')):
                    sg=sg.sort_values('gap'); y=(sg.lower+sg.upper)/2
                    y[sg.prior_heading_gap]=np.nan
                    ax.plot(sg.gap,y,color=color(model,i),ls=STYLES[i],marker='o',ms=3,label=f'C={size:g}')
                    for r in sg[sg.status=='right_censored_saved_heading_at_completed_horizon'].itertuples():
                        ax.plot(r.gap,r.lower,marker='^',mfc='white',mec=color(model,i),ls='none')
                    for r in sg[sg.status.str.contains('gap')].itertuples():
                        ax.axvspan(r.gap-.04*(.1 if model=='spring' else 1),r.gap+.04*(.1 if model=='spring' else 1),color='.9')
                ax.set_title(name,color='black');ax.set_yscale('log');ax.set_xlabel('Cue edge gap (native length)');ax.set_ylabel('Native time');ax.legend(frameon=False,fontsize=8)
            f.text(.08,.025,'Open up-triangles: lower bounds; grey: undefined heading. Midpoints of saved brackets; quality gate R ≥ 0.05. No cap is a measured time.',fontsize=8)
            f.tight_layout(rect=[0,.055,1,.95]);pdf.savefig(f);plt.close(f)
    writer=PdfWriter();preserved=[]
    prior=DEV/'revision_round2_20261022/figure3'
    for path in (added,*EXTRA_SUPPLEMENT_SOURCES,prior/'figure3_main.pdf',prior/'figure3_supplement.pdf'):
        writer.append(path)
        preserved.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),pages=len(PdfReader(path).pages)))
    writer.write(OUT/'figure3_round3_supplement.pdf')
    write(META/'supplement_preservation.json',dict(generated_at=native.now(),sources=preserved,
        total_pages=len(PdfReader(OUT/'figure3_round3_supplement.pdf').pages),
        complete_round2_main_and_supplement_retained=True))
    documents()


def documents():
    statuses=[read(RUN/'workers'/label/'status.json') for label,_,_ in MAIN_EXAMPLES]
    cue_rule=[]
    for (label,_,_),s in zip(MAIN_EXAMPLES,statuses):
        withdrawal='95% heading' if s['withdrawal_condition']=='95_percent_heading' else 'the 20,000 horizon (95% not reached)'
        stop='stationarity' if s['stop_condition']=='stationarity' else 'the four-times-duration limit'
        cue_rule.append(f'{label}: cue off at {withdrawal}, t = {s["cue_on_duration"]:.1f}; relaxed for {s["relaxation_duration"]:.1f}, stopped at {stop}')
    title="Cue placement relative to the incumbent’s edge switches reorientation between sliding the incumbent and takeover by a new region; enlargement slows both."
    legend=[
        'A. Native Wave-pinning active species u versus arc position centred on the incumbent’s initial centre, for (i) C = 30, gap −1; (ii) C = 30, gap 4; (iii) C = 47.6, gap −1; and (iv) C = 47.6, gap 4. Rows show 0%, 50% and 90% of the target heading under continuous cue, then the final cue-free profile. Grey shading marks the fixed cue support in every row; open diamonds mark an inter-region saddle where resolved, and black ticks on the arc axis mark heading. Each column shares its y-scale.',
        'B. Time to 50% of the angle from initial activity-vector heading to the cue direction versus initial cue edge gap (near edge of cue support minus the incumbent’s 10%-contrast boundary). C. Mean angular speed over that turn, equal to half the target angle divided by t₅₀. Time and length use each model’s native units; speed is radians per native time. Spring retains its native mechanical units. Lighter solid, intermediate dashed and darker dotted curves denote increasing C, labelled at each line end. Enlarged filled i–iv markers identify A in both Wave-pinning plots. Both readouts are evaluated with the cue continuously on; points interpolate the saved 50% crossing brackets.',
        'Vertical lines mark the first tested gap with two resolved regions and a saddle at the 50% crossing for each C, using that size’s line style; the existing periodic watershed detector is applied to the native saved-field interpolant. Light grey bands mark sampled near-antipode placements with undefined heading at one or more sizes (activity-vector coherence R < 0.05 or a cancellation interval); that size’s curve breaks there. Open up-triangles mark censored time lower bounds, and matching open down-triangles mark speed upper bounds; horizon caps are not measured values.',
        'For the relaxation row only, cue-off occurs at 95% heading or at 20,000 if 95% is not reached, restarting from the full native event state after an identical cue-on replay. No-cue relaxation stops when max |∂u/∂t| reaches the existing tolerance 2.2683 × 10⁻⁷ (native activity per native time; the original 10⁻⁴ scaled-RHS tolerance over 1,000 time units), or after four times the cue-on duration. '+'; '.join(cue_rule)+'.',
        'deBelly is excluded from rings.',
        'Cellular measurement motivated by this assay: photoactivation at varying offsets from the leading-edge boundary in WGD+ versus WGD− cells, with time to half-turn as the readout.'
    ]
    if SINGLE_C_LAYOUT:
        legend[0]=legend[0].replace('for (i) C = 30, gap −1; (ii) C = 30, gap 4; (iii) C = 47.6, gap −1; and (iv) C = 47.6, gap 4.',
            'at C = 30 for (i) gap −1 (straddling) and (ii) gap 4 (separated).')
        legend[1]=legend[1].replace('filled i–iv markers','filled i–ii markers')
    (OUT/'legend.md').write_text('**Figure 3. '+title+'**\n\n'+'\n\n'.join(legend)+'\n')
    rel=lambda p:str(p.relative_to(ROOT))
    section=dict(id='figure3',title=title,status='round3_complete',legend=legend,
        role='Continuous-cue ring reorientation across cue placements and native sizes.',
        narrative=['Straddling Wave-pinning examples slide the incumbent, whereas separated examples develop a second activity region and saddle. Both examples take longer to reach half of the cue-direction target in the enlarged ring.'],
        measurements=['Photoactivate at varying offsets from the leading-edge boundary in WGD+ versus WGD− cells; measure time to half-turn.'],
        limitations=['Native model clocks are not calibrated across models. Topology at the 50% crossing uses saved-field interpolation; the C = 47.6, gap 0.5 Wave-pinning bracket changes topology between its endpoints. Undefined heading cannot supply a timing bound. Four fixed-parameter relaxation protocols are the only new simulations.'],
        assets=[dict(path=rel(OUT/'figure3_round3.png'),native_pdf_page=rel(OUT/'figure3_round3.pdf'),
            page_size=[16.2,8.6] if SINGLE_C_LAYOUT else [11.7,13.6],minimum_printed_width_inches=16.2 if SINGLE_C_LAYOUT else 11.7,page_title='')],
        main_pdf=rel(OUT/'figure3_round3.pdf'),supplement_pdf=rel(OUT/'figure3_round3_supplement.pdf'),
        reader=rel(OUT/'index.html'),sources=[dict(label='50% readouts and topology',path=rel(DATA/'half_heading_readouts.csv')),
            dict(label='Four native relaxation protocols',path=rel(RUN/'manifest.json'))])
    if SINGLE_C_LAYOUT:
        section['narrative']=['At C = 30, the straddling Wave-pinning example slides the incumbent, whereas the separated example develops a second activity region and saddle. B and C retain all three sizes and show slower half-heading reorientation in enlarged rings.']
    write(OUT/'legend.json',section)
    generated=native.now()
    collection=read(META/'collection.json')
    body='<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 3 round 3</title><style>body{max-width:1150px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}img{width:100%}a{color:#166aa2}</style>'
    body+=f'<h1>Figure 3 — round 3</h1><p>Data generated {generated}; all four native relaxation producers finished.</p><p><a href="figure3_round3.pdf">Main PDF</a> · <a href="figure3_round3_supplement.pdf">Complete supplement, including every round-2 page</a> · <a href="legend.md">Text legend</a> · <a href="supporting/data/half_heading_readouts.csv">50% source coordinates and topology</a></p><img src="figure3_round3.png" alt="Activity profiles and half-heading time and speed versus cue gap">'
    import html
    body+='<h2>'+html.escape(title)+'</h2>'+''.join('<p>'+html.escape(p)+'</p>' for p in legend)
    body+=f'<p><a href="README.md">Methods, source provenance and reproducibility</a> · <a href="supporting/metadata/collection.json">Collection audit</a></p></html>'
    (OUT/'index.html').write_text(body)
    (OUT/'README.md').write_text(f'''# Figure 3 round 3

[Main PDF](figure3_round3.pdf) · [legend](legend.md) · [complete supplement](figure3_round3_supplement.pdf) · [reader](index.html).

Generated {generated}; four relaxation protocols complete. Original evidence is immutable. No kinetic retuning or new response campaign. Readout source: `{rel(SOURCE/'A1_heading_target_time_brackets.csv')}`, target_fraction = 0.5, coherence_floor = 0.05, non-control placements. Coincident cue/heading axes have no nonzero turn and are not log-axis observations. All other placements remain in the coordinate CSV, including undefined headings.

Time coordinates are roots of the linear native-field interpolant within the original saved brackets, rather than new integrations. Complete bracket endpoints and reciprocal speed intervals are exported. Topology is classified by the existing two-region watershed/saddle detector at that same interpolated state; endpoint agreement is checked. Undefined-heading cells extend halfway to adjacent tested gaps, separately for each size. Grey bands show their union across sizes; only the affected size's curve is interrupted. Regime lines report the first tested saddle at 50%, not an inferred continuous bifurcation. No-saddle sizes have no line. The 50% source contains no positive-target censored observations at the selected gate; triangle meanings remain defined for reproducibility. There are {collection['bracket_topology_disagreements']} crossings whose bracket endpoints disagree on two-region topology.

A's cue-on profiles reuse the original full fields, reconstructed at 50/90% inside saved brackets. Only its relaxation row is new. Four cue-on replays start from the exact immutable basal input, using the original frozen model and BDF settings (rtol 10⁻⁹, atol 10⁻¹¹, max step 10). They stop at the dense 95% heading event or 20,000. Cue-off starts from exactly that full state. Existing active-field scaled-RHS tolerance is 10⁻⁴ over 1,000 time units and active scale 2.2683, hence max |∂u/∂t| ≤ 2.2683 × 10⁻⁷. The original all-species basal qualification gate is distinct from this requested active-species relaxation stop. Full cue-on/cue-off states, accepted-step mass/positivity checks and stop conditions are in `{rel(RUN)}`.

Every round-2 main and supplement page is appended to the new supplement, after 90/99%-heading pages. Original files are untouched. Native clocks have no cross-model calibration; observed migration and conditional OU persistence remain separate evidence from these simulated polarity readouts.

Rebuild using the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py freeze
# Submit the exact sbatch_arguments in the generated submission.json once.
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py collect
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py render
```
''')
    if SINGLE_C_LAYOUT:
        note='\nThis layout revision places the eight C = 30 profiles on the left and B above C on the right. Only i–ii are highlighted in the main. The previous four-placement main and its legend are appended to the supplement, retaining the C = 47.6 profiles and their horizon qualification. All original source CSVs and four completed relaxation protocols are reused; no simulations were run for this layout revision.\n'
        with (OUT/'README.md').open('a') as stream:stream.write(note)
    write(META/'build_receipt.json',dict(generated_at=generated,completed=True,
        main_sha256=sha(OUT/'figure3_round3.pdf'),new_protocols=4,retuned=False,
        all_round2_preserved=True,collection=collection,
        layout_revision_native_integrations=0,main_profile_sizes=[30.] if SINGLE_C_LAYOUT else [30.,47.6]))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['freeze','worker','collect','render'])
    parser.add_argument('--slot',type=int)
    args=parser.parse_args()
    if args.command=='freeze':print(json.dumps(freeze(),indent=2))
    elif args.command=='worker':worker(args.slot)
    elif args.command=='collect':collect()
    else:render()
