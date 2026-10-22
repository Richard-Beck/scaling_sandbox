"""Revised deBelly Figures 8, 10 and 11 from retained native evidence.

No integrations or scientific input changes. New measurements identify spatial
flow, the Rho maximum, and MCA progress through withdrawal and recovery. The
short-length stability panel replaces the redundant old Figure 11.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D

import length_dev_figure_x_overview as overview
import length_dev_figure_z_mca_distance_speed as distance

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'revision_20261013/debelly_sequence'
OLD = REPORT / 'revision_20261010/debelly'
SOURCES = {}
BLUE, ORANGE, GREEN, PURPLE, INK, GRAY = '#0072B2', '#D55E00', '#009E73', '#762A83', '#252525', '#777777'
COLORS = {18: '#777777', 23: BLUE, 31: PURPLE, 39: ORANGE}
GEOS = distance.GEOS
GEO_COLORS = dict(zip(GEOS, [GRAY, BLUE, ORANGE]))
GEO_NAMES = dict(zip(GEOS, ['Fixed 3 µm', 'Relative 20%', 'Slope 0.24']))
CHAPTERS = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def rows(path):
    with source(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def csvout(name, records):
    path = OUT / name
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    outside = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > bounds.width + 1 or box.y1 > bounds.height + 1:
                outside.append(artist.get_text())
    assert not outside, outside
    paths = []
    for ext in ['png', 'pdf', 'svg']:
        path = OUT / (stem + '.' + ext)
        fig.savefig(path, dpi=200, facecolor='white')
        paths.append(path)
    plt.close(fig)
    return paths


def finish(stem, number, title, paragraphs, caption, methods, audit, artifacts, extra_links=()):
    cp = OUT / (stem + '_caption.txt'); cp.write_text(caption + '\n')
    mp = OUT / (stem + '_methods.txt'); mp.write_text(methods + '\n')
    ap = OUT / (stem + '_audit.json')
    audit.update(passed=True, new_integrations=0, plot_text_within_page=True)
    ap.write_text(json.dumps(audit, indent=2) + '\n')
    artifacts = list(artifacts) + [cp, mp, ap]
    meta = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer='finished', passed=True,
                renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(__file__),
                artifacts={p.name: sha(p) for p in artifacts}, sources=dict(SOURCES),
                sha256={p.name: sha(p) for p in artifacts}, source_sha256=dict(SOURCES),
                caption=caption, audit=audit, new_integrations=0,
                reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_revision_20261013_sequence.py')
    manifest = OUT / (stem + '.manifest.json')
    manifest.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + '\n')
    rel = lambda p: str(p.relative_to(REPORT))
    CHAPTERS.append(dict(id=stem.lower(), figure_label=str(number), title=title,
                         image=rel(OUT / (stem + '.png')), paragraphs=paragraphs, caption=caption,
                         manifest=rel(manifest), sources=[[rel(mp), 'Full methods and measurement definitions'],
                         [rel(ap), 'Numerical checks and source identities']] + list(extra_links)))


def style(ax, letter=None):
    ax.grid(alpha=.14, linewidth=.6)
    if letter:
        ax.text(-.13, 1.07, letter, transform=ax.transAxes, weight='bold', fontsize=13)


def interp_field(tr, key, time):
    a = tr['fields'][key]
    return np.array([np.interp(time, tr['time'], a[:, i]) for i in range(a.shape[1])])


def figure8():
    traces, original_measurements, checks, arrays = overview.load()
    SOURCES.update(overview.SOURCES)
    measurements = []
    for length, tr in traces.items():
        for t in [0., 5., 30., 60., tr['duration'], tr['duration'] + 30., tr['duration'] + 60.]:
            f = {key: interp_field(tr, key, t) for key in overview.FIELDS}
            measurements.append(dict(length_um=length, time_after_onset_s=t,
                cue_duration_s=tr['duration'], time_from_withdrawal_s=t-tr['duration'],
                rac_old_left=float(f['rac'][0]), rac_challenged_right=float(f['rac'][-1]),
                rho_old_left=float(f['rho'][0]), rho_challenged_right=float(f['rho'][-1]),
                rho_maximum_x_over_L=float(tr['x'][np.argmax(f['rho'])] / length),
                rightward_cortical_flow_fraction=float(np.mean(f['velocity'] > 0)),
                mean_cortical_flow_um_s=float(f['velocity'].mean()),
                bound_mca_centroid_um=float(f['mca_bound'] @ tr['x'] / f['mca_bound'].sum())))
    fig = plt.figure(figsize=(12.6, 13.4))
    gs = fig.add_gridspec(5, 4, width_ratios=[23, 31, 39, 2.2], left=.11,
                         right=.93, bottom=.065, top=.962, hspace=.23, wspace=.18)
    for row, (key, label, unit) in enumerate(zip(overview.FIELDS, overview.LABELS, overview.UNITS)):
        visible = [tr['fields'][key][tr['time'] <= overview.LIMIT] for tr in traces.values()]
        lo, hi = min(float(a.min()) for a in visible), max(float(a.max()) for a in visible)
        if key == 'rac':
            norm, cmap = TwoSlopeNorm(vmin=lo, vcenter=1., vmax=hi), 'viridis'
        elif key in ['velocity', 'tension']:
            norm, cmap = Normalize(-max(abs(lo), abs(hi)), max(abs(lo), abs(hi))), 'RdBu_r'
        else:
            norm, cmap = Normalize(lo, hi), 'viridis'
        for col, length in enumerate(overview.LENGTHS):
            ax = fig.add_subplot(gs[row, col]); tr = traces[length]
            mesh = ax.pcolormesh(tr['x'], tr['time'], tr['fields'][key], shading='nearest',
                                norm=norm, cmap=cmap, rasterized=True)
            ax.set_yscale('symlog', linthresh=120, linscale=.65)
            ax.set(xlim=(0, length), ylim=(1200, -60), xticks=[0, length/2, length],
                   yticks=[-60, 0, 60, 120, 300, 600, 1200])
            ax.set_yticklabels(['−60', '0', '60', '120', '300', '600', '1200'])
            ax.tick_params(labelsize=10, length=3)
            ax.axhline(0, color='white', linewidth=1.6)
            ax.axhline(60, color='white', linestyle=':', linewidth=1.8)
            ax.axhline(tr['duration'], color='white', linestyle='--', linewidth=2)
            ax.axhline(tr['duration'], color='#222222', linestyle='--', linewidth=.7)
            if col:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(label + '\nTime from cue onset (s)', fontsize=11, labelpad=8)
                ax.text(-.28, 1.03, 'ABCDE'[row], transform=ax.transAxes, fontsize=13, weight='bold')
            if row == 0: ax.set_title(f'L = {length} µm', fontsize=12, pad=12)
            if row == 4: ax.set_xlabel('Reference position x (µm)', fontsize=11)
        cb = fig.colorbar(mesh, cax=fig.add_subplot(gs[row, 3]))
        if key == 'rac': cb.set_ticks([.6, .8, 1., 2., 4.])
        cb.ax.tick_params(labelsize=10, length=2)
        cb.set_label(unit, fontsize=10, labelpad=6)
        cb.ax.yaxis.get_offset_text().set_fontsize(10)
    stem = 'Figure_8_reversal_overview'
    artifacts = save(fig, stem)
    artifacts += [csvout(stem + '_spatial_response.csv', measurements),
                  csvout(stem + '_diagnostics.csv', original_measurements)]
    p = OUT / (stem + '_fields.npz'); np.savez_compressed(p, **arrays); artifacts.append(p)
    caption = ('Figure 8. Rac activation and cortical-flow reversal precede most MCA redistribution, with a distinct response at 39 µm. '
        'Rows A–E show Rac, Rho, bound MCA, cortical velocity and membrane tension at 23, 31 and 39 µm. '
        'Solid lines mark cue onset, dotted lines 60 s, and dashed lines withdrawal of the shortest tested 1%-strength-passing '
        'relative-20% pulses (232.244898, 285.280612 and 809.081633 s; peak 3.6). '
        'At 60 s, Rho maxima lie inside the cell at x/L = 0.619, 0.634 and 0.749. '
        'Rightward flow occupies 100%, 99.5% and 6.0% of the 201 native cells at that time. '
        'All columns share one absolute color scale for each field. The Rac scale places half its palette below 1.0 and half above; '
        '1.0 is a display breakpoint. Time increases downward, linearly through 120 s and logarithmically afterward. '
        'The pre-cue strip repeats the exact basal checkpoint for visibility. Positions are nominal reference coordinates and column widths scale with length. '
        'No fields are reset at withdrawal; cue removal changes their subsequent dynamics. Scalar Rac strength is evaluated after 10000 s of basal recovery, '
        'beyond the displayed interval; separate native-field flags remain in the exported diagnostics.')
    paragraphs = [
        'A finite Rac cue over the rightmost 20% of the cell challenges the left-polarized basal deBelly state. The raised-cosine patch has the same peak, 3.6, at every length. In this model, Rac and Rho locally inhibit one another. Lower bound membrane–cortex attachment (MCA) increases Rac activation, Rac drives endpoint protrusion, tension raises Rho activation, and Rho gradients drive cortical motion. The cue therefore changes chemical activity and the motion that redistributes MCA.',
        'By 5 s, Rac at the challenged right end has risen from 0.651 to 1.071, 0.623 to 1.027 and 0.562 to 0.928 at 23, 31 and 39 µm (Fig. 8A). Rho falls there while its maximum moves into the interior (Fig. 8B): at 60 s the maximum is at 62%, 63% and 75% of cell length. Bound MCA remains biased toward the original right rear at that time (Fig. 8C); its centroid has moved only 0.237, 0.199 and 0.216 µm toward the left.',
        'Cortical flow responds very differently at the longest length (Fig. 8D). Five seconds after onset, rightward flow occupies 87.1%, 68.2% and 8.0% of the 23-, 31- and 39-µm cells; at 60 s those fractions are 100%, 99.5% and 6.0%. Most of the 39-µm cortex continues moving toward the old left pole even though Rac has risen at the right end. Its spatial mean flow crosses zero near 485 s, and its signed tension moment near 530 s (Fig. 8E); these are distinct measures of the delayed mechanical response.',
        'Withdrawal substantially reduces cue-supported Rac. In the next 30 s, right-end Rac falls from 4.29 to 1.95, 4.47 to 2.07 and 4.64 to 2.19 across the three lengths. MCA continues redistributing during that relaxation: its centroid moves another 0.099, 0.114 and 0.029 µm leftward over the same interval. The cue-duration readout therefore precedes completion of the basal MCA rearrangement. Figure 9 resolves the local processes that produce it.'
    ]
    methods = (caption + '\n\nFields are read without resampling for the kymographs. Descriptive values at exact stated times use linear '
        'interpolation of saved native fields. Positive-flow fraction counts cells with v > 0, each with equal nominal finite-volume width. '
        'Rho maximum is the nominal center of its largest-density cell. Endpoint values are first/last native cells. '
        'Bound MCA centroid is sum(x*B)/sum(B); its shift is positive toward the left. Flow sign crossing refers to spatial mean velocity. '
        'Tension sign crossing refers to sum((2x/L-1)*T)/sum(abs(T)), with a reversed sign for convenience in the original export. '
        'Neither identifies complete spatial reversal. All source archives, starts, withdrawal continuity, cue selection and scalar gates '
        'are verified by the original loader; field flags are not used to exclude curves. No new native step is taken.')
    finish(stem, 8, 'At 60 s, flow is rightward across almost all of the 23–31-µm cells and mostly leftward at 39 µm',
           paragraphs, caption, methods, dict(original_checks=checks, measured_times_s=[0,5,30,60,'D','D+30','D+60'],
           local_spatial_response_measured=True, rho_maximum_not_equated_to_first_moment=True), artifacts)


def load_distance():
    data, checks = distance.measure(); SOURCES.update(distance.SOURCES)
    manifest = read(distance.RUN / 'manifest.json')
    traces = {}
    for row in data:
        job = next(j for j in manifest['jobs'] if j['case']['length'] == row['length_um'] and j['case']['condition'] == row['condition'])
        p = ROOT / job['source'] if job['source'] else distance.RUN / f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        with np.load(source(p)) as z:
            c = distance.centroid(z['recovery_field_mca_bound'], z['x'])
            full = row['initial_distance_to_reflected_state_um']
            cp = distance.centroid(z['pulse_field_mca_bound'], z['x'])
            cue_progress = (cp[0]-cp)/full
            maximum = int(np.argmax(cue_progress))
            row['maximum_cue_fraction_of_reflected_centroid_distance'] = float(cue_progress[maximum])
            row['time_of_maximum_cue_progress_s'] = float(z['pulse_time'][maximum])
            row['cue_retreat_after_maximum_fraction'] = float(cue_progress[maximum]-cue_progress[-1])
            row['final_recovery_centroid_um'] = float(c[-1])
            row['final_recovery_fraction_of_reflected_distance'] = (row['initial_centroid_um'] - float(c[-1])) / full
            row['withdrawal_to_final_recovery_shift_um'] = row['withdrawal_centroid_um'] - float(c[-1])
            row['fraction_of_actual_onset_to_final_shift_during_cue'] = row['centroid_shift_during_cue_um'] / (row['initial_centroid_um'] - float(c[-1]))
            row['final_centroid_error_from_reflected_um'] = float(c[-1]) - row['reflected_target_centroid_um']
            row['final500_centroid_range_um'] = float(np.ptp(c[z['recovery_time'] >= 9500]))
            if row['condition'] == 'fixed_relative' and row['length_um'] in [18,23,31,39]:
                traces[int(row['length_um'])] = dict(time=np.r_[z['pulse_time'], row['cue_duration_s'] + z['recovery_time'][1:]],
                    progress=np.r_[(cp[0]-cp)/full, (cp[0]-c[1:])/full], duration=row['cue_duration_s'])
            np.testing.assert_allclose(row['cue_average_speed_um_s'] * row['cue_duration_s'], row['centroid_shift_during_cue_um'], atol=1e-14)
    return data, checks, traces


def figure10():
    data, checks, traces = load_distance()
    linear = read(ROOT / 'data/runs/debelly_corrected_zero_crossing_20261008/linear.json')
    assert read(ROOT / 'data/runs/debelly_corrected_zero_crossing_20261008/independent_audit.json')['passed']
    linear_provenance = read(ROOT / 'data/runs/debelly_corrected_zero_crossing_20261008/provenance.json')
    assert sha(source(ROOT / 'code/core/models/debelly.py')) == linear_provenance['model_sha256']
    near_root = ROOT / 'data/runs/debelly_mca_threshold_profiles_20261010'
    near = rows(near_root / 'results/length_summary.csv')
    near_audit = read(near_root / 'results/audit.json'); assert near_audit['passed']
    near_manifest = read(near_root / 'manifest.json')
    for name,digest in near_manifest['sources'].items():
        assert sha(source(near_root / 'frozen_sources' / name)) == digest
    for record in near_audit['inputs']:
        assert record['convergence']['verified']
        assert sha(source(near_root / record['file'])) == record['sha256']
    flat = [float(r['length_um']) for r in near if float(r['bound_mca_range']) < 1e-4]
    polarized = [float(r['length_um']) for r in near if float(r['bound_mca_range']) >= 1e-4]
    loss = [max(flat), min(polarized)]; assert loss == [17.9,17.925]
    lx = np.array([r['length'] for r in linear]); ly = np.array([r['rate'] for r in linear])
    i = np.flatnonzero(ly[:-1]*ly[1:] < 0); assert len(i) == 1
    onset = float(np.interp(0, ly[i[0]:i[0]+2], lx[i[0]:i[0]+2]))
    fig, axes = plt.subplots(2,3,figsize=(13.2,9.5))
    fig.subplots_adjust(left=.08,right=.98,bottom=.09,top=.95,hspace=.43,wspace=.36)
    a,b,c,d,e,f = axes.flat
    for ax, letter in zip(axes.flat,'ABCDEF'): style(ax,letter)
    for ax in [a,b,c,d,e]:
        ax.set(xlim=(17,40),xticks=[18,23,31,39],xlabel='Cell length (µm)')
        ax.axvspan(*loss,color=GRAY,alpha=.45,zorder=0)
        ax.axvspan(loss[1],onset,color='#E69F00',alpha=.13,zorder=0)
        ax.axvline(onset,color='#B27600',linewidth=1,linestyle=':')
    for geo in GEOS:
        rr=sorted([r for r in data if r['condition']==geo],key=lambda r:r['length_um'])
        L=np.array([r['length_um'] for r in rr])
        a.plot(L,[r['cue_duration_s'] for r in rr],color=GEO_COLORS[geo],marker='o',ms=4,lw=1.7,label=GEO_NAMES[geo])
        b.plot(L,[100*r['cue_fraction_of_reflected_centroid_distance'] for r in rr],color=GEO_COLORS[geo],marker='o',ms=4,lw=1.7)
        # Open versus filled markers preserve the separate native-field diagnostic.
        for ax,key,mult in [(a,'cue_duration_s',1),(b,'cue_fraction_of_reflected_centroid_distance',100)]:
            for r in rr:
                ax.plot(r['length_um'],mult*r[key],marker='o',ms=4,mec=GEO_COLORS[geo],
                    mfc=GEO_COLORS[geo] if r['native_field_tail_pass'] else 'white',linestyle='none',zorder=4)
    a.set(ylabel='Shortest passing cue (s)',ylim=(0,870));a.legend(frameon=False,fontsize=10,loc='upper center')
    b.axhline(50,color=PURPLE,ls=':',lw=1.3);b.axhline(100,color=INK,ls='--',lw=1.2)
    b.set(ylabel='Centroid shift during cue\n(% of reflected-state distance)',ylim=(0,108))
    c.plot(lx,1000*ly,color=INK,lw=2);c.axhline(0,color=GRAY,lw=.8)
    c.set(ylabel='Homogeneous growth rate\n(10⁻³ s⁻¹)',ylim=(-.40,.80))
    c.text(.98,.98,f'Zero near {onset:.2f} µm',transform=c.transAxes,ha='right',va='top',fontsize=10)
    rr=sorted([r for r in data if r['condition']=='fixed_relative'],key=lambda r:r['length_um'])
    L=np.array([r['length_um'] for r in rr])
    for key,color,ls,label in [
        ('initial_distance_to_reflected_state_um',INK,'--','Initial → reflected'),
        ('initial_distance_to_midpoint_um',PURPLE,':','Initial → midpoint'),
        ('centroid_shift_during_cue_um',BLUE,'-','During cue'),
        ('withdrawal_to_final_recovery_shift_um',GREEN,'-.','After withdrawal')]:
        d.plot(L,[r[key] for r in rr],color=color,ls=ls,marker='o',ms=3.5,lw=1.7,label=label)
    d.set(ylabel='Relative-20% centroid shift (µm)',ylim=(0,5.3));d.legend(frameon=False,fontsize=10,loc='upper left')
    for key,color,ls,label in [('early_mean_speed_um_s',GREEN,'-','First 60 s'),('late_mean_speed_um_s',ORANGE,'-','Last 60 s'),
                            ('cue_average_speed_um_s',BLUE,'-','Whole cue')]:
        e.plot(L,[1000*r[key] for r in rr],color=color,ls=ls,marker='o',ms=3.5,lw=1.7,label=label)
    e.axhline(0,color=GRAY,lw=.8);e.set(ylabel='Relative-20% centroid speed\n(10⁻³ µm/s)');e.legend(frameon=False,fontsize=10,loc='upper left')
    for length,tr in traces.items():
        f.plot(tr['time'],100*tr['progress'],color=COLORS[length],lw=1.6,label=f'{length} µm')
        f.plot(tr['duration'],100*np.interp(tr['duration'],tr['time'],tr['progress']),marker='o',mfc='white',mec=COLORS[length],ms=6)
    f.axhline(50,color=PURPLE,ls=':',lw=1.2);f.axhline(100,color=INK,ls='--',lw=1.1)
    f.set_xscale('symlog',linthresh=120,linscale=.7)
    f.set(xlim=(0,11000),ylim=(-5,115),xticks=[0,120,1000,10000],
          xlabel='Time from cue onset (s)',ylabel='Relative-20% centroid progress\n(% of reflected-state distance)')
    f.set_xticklabels(['0','120','1000','10000']);f.legend(frameon=False,fontsize=10,loc='lower right',ncol=2)
    stem='Figure_10_MCA_progress_and_regimes';artifacts=save(fig,stem)
    table=csvout(stem+'.csv',data);artifacts.append(table)
    artifacts.append(csvout(stem+'_linear_growth.csv',linear))
    long=[r for r in rr if r['length_um']>=21]
    plateau=[r['cue_fraction_of_reflected_centroid_distance'] for r in long]
    actual=[r['fraction_of_actual_onset_to_final_shift_during_cue'] for r in long]
    # The requested alternative is retained as a companion because the two
    # axes share the exact shift/D identity and do not independently explain D.
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.5),sharey=True)
    fig.subplots_adjust(left=.075,right=.89,bottom=.20,top=.88,wspace=.24)
    norm=Normalize(18,39)
    for ax,geo in zip(axes,GEOS):
        rrgeo=sorted([r for r in data if r['condition']==geo],key=lambda r:r['length_um'])
        ax.plot([r['centroid_shift_during_cue_um'] for r in rrgeo],[1000*r['cue_average_speed_um_s'] for r in rrgeo],color='#bbbbbb',lw=.8,zorder=1)
        scatter=ax.scatter([r['centroid_shift_during_cue_um'] for r in rrgeo],[1000*r['cue_average_speed_um_s'] for r in rrgeo],
            c=[r['length_um'] for r in rrgeo],cmap='viridis',norm=norm,s=40,zorder=2)
        for r in rrgeo:
            if r['length_um'] in [18,21,31,39]:
                ax.annotate(str(r['length_um']),(r['centroid_shift_during_cue_um'],1000*r['cue_average_speed_um_s']),
                            xytext=(4,5),textcoords='offset points',fontsize=10)
        ax.set_title(GEO_NAMES[geo],fontsize=12);ax.set_xlabel('Net cue shift (µm)');style(ax)
        ax.margins(x=.17,y=.18)
    axes[0].set_ylabel('Whole-cue centroid speed\n(10⁻³ µm/s)')
    cb=fig.colorbar(scatter,cax=fig.add_axes([.915,.22,.017,.61]));cb.set_label('Cell length (µm)');cb.ax.tick_params(labelsize=10)
    optional=save(fig,'Figure_10_optional_shift_speed');artifacts.extend(optional)
    caption=('Figure 10. Short-cell reversal, MCA progress and later length effects occupy different regimes. '
        'A) Shortest tested 1%-strength-passing cue duration under all three geometries. B) Net bound-MCA centroid shift at withdrawal '
        'as a percentage of its initial separation from the reflected basal centroid: 50% reaches the cell midpoint and 100% the reflected centroid. '
        'C) Homogeneous-state asymmetric growth rate; its sign changes near 19.53 µm. Orange shading in A–E spans the interval between '
        'sampled polarized-branch loss at 17.9–17.925 µm (gray) and homogeneous instability. '
        'D) Relative-20% distances before, during and after the cue. E) Relative-20% signed centroid speeds over the first 60 s, final 60 s and whole cue. '
        'F) Relative-20% progress through pulse and 10000-s basal recovery at 18, 23, 31 and 39 µm; open circles mark withdrawal. '
        'Its time axis is linear through 120 s and logarithmic afterward. '
        'At 18 µm the cue supplies 91.5% of the reflected-centroid shift; at 19 µm 53.5%; at 21–39 µm 32.7–40.5%. '
        'The full relative-20% recovery ends at 99.55–100.01% across the twelve lengths. '
        'Filled/open markers in A–B retain passing/failing native-field tail flags separately from scalar success. '
        'All geometries share the exact same corrected basal checkpoint at each length. Whole-cue speed equals net cue shift divided by the selected duration; '
        'it is a descriptive summary. Progress is measured at the selected Rac-gated duration, rather than at an independently tested MCA threshold. '
        'The optional companion shows the shift–speed view colored by length.')
    paragraphs=[
        'The duration curves have a pronounced short-length decline and further nonmonotonic structure (Fig. 10A). Relative-20% duration falls from 472 s at 18 µm to 226 s at 21 µm, rises to 285 s at 31 µm, dips to 275 s at 33 µm, then reaches 809 s at 39 µm. Fixed-width duration rises to 351 s at 27 µm, falls to 308 s at 31 µm and remains 310 s at 33 µm. Constant-slope duration decreases to 150 s at 23 µm before rising to 293 s at 39 µm.',
        'The short-cell decline coincides with a distinct polarity regime (Fig. 10C). The retained polarized branch disappears between 17.9 and 17.925 µm, whereas small asymmetric perturbations of the homogeneous state decay until about 19.53 µm. Thus 18–19-µm cells can retain a finite cue-established pole while homogeneity remains linearly stable; above 19.53 µm, small asymmetry grows. The later dips at 27–33 µm occur entirely above this transition.',
        'MCA is much closer to its reversed basal centroid at withdrawal in the shortest cell (Fig. 10B,D). A relative-20% pulse supplies 91.5% of the initial-to-reflected shift at 18 µm and 53.5% at 19 µm. At 21–39 µm it supplies 32.7–40.5%, leaving the centroid on the original side of the midpoint. This restricted range is approximately a fixed fraction of the initial centroid separation. Using the actual final recovery centroid instead changes the range only to 32.8–40.5%.',
        'The remaining redistribution occurs after cue withdrawal (Fig. 10F). At 23 and 39 µm, the relative-20% cue shifts the centroid by 1.040 and 1.982 µm; basal recovery adds 1.678 and 2.996 µm. The 39-µm centroid advances farther earlier in the cue, then retreats before withdrawal. Whole-cue speed rises from 4.48 × 10^-3 µm/s at 23 µm to 5.51 × 10^-3 at 33 µm, then falls to 2.45 × 10^-3 at 39 µm (Fig. 10E). These averages describe achieved progress over unequal selected durations. Figure 11 explains the starting mechanical state, and Figure 12 examines the long-cell exception during reversal.'
    ]
    methods=(caption+'\n\nAll 36 pulses are read by the existing native-history selection and gate audit. '
        'The mass-normalized centroid is c=sum(x*B)/sum(B), using the 201 equal nominal finite volumes. '
        'The reflected reference is L-c(0); reflected distance is 2c(0)-L. Cue progress is [c(0)-c(D)]/[2c(0)-L]. '
        'After-withdrawal shift is c(D)-c(D+10000). Final recovery values are recorded endpoints, not declared full-field steady states. '
        'Final-500 centroid range and reflected-target error are exported. Early/late speeds use signed differences over [0,60] and [D-60,D], '
        'with linear interpolation at exact boundaries. Every selected D is at least 120 s. Whole-cue speed is [c(0)-c(D)]/D, verified algebraically. '
        'The growth-rate curve is the independently audited infinitesimal tangent calculation, separate from finite branch continuation. '
        'Its zero is linearly interpolated between the neighboring sampled signs. The branch bracket uses the audited near-threshold '
        'continuation with convergence and holdout checks and a bound-MCA range cutoff of 1e-4. '
        'Below the bracket, this continuation supplies no retained polarized branch; it does not prove absence of every possible state. '
        'The optional scatter shares algebraically related variables and is not an independent explanation or fit of duration.')
    finish(stem,10,'MCA progress at withdrawal changes sharply across the short-cell regime',paragraphs,caption,methods,
           dict(native_history_checks=checks, homogeneous_instability_um=onset, sampled_branch_loss_um=loss,
                relative20_progress_range_L21_39=plateau and [min(plateau),max(plateau)],
                relative20_actual_final_progress_range_L21_39=[min(actual),max(actual)],
                cue_speed_identity_verified=True, final_recovery_not_assumed_full_field_steady=True),artifacts,
           [[str((OUT/'Figure_10_optional_shift_speed.pdf').relative_to(REPORT)),'Optional shift–speed view colored by length'],
            [str(table.relative_to(REPORT)),'All 36 centroid, duration and recovery measurements']])


def audited_old(stem):
    meta=read(OLD/(stem+'.manifest.json'))
    for p,digest in meta['artifacts'].items():
        assert sha(source(OLD/p))==digest,p
    for p,digest in meta['sources'].items():
        assert sha(source(ROOT/p))==digest,p
    assert meta['audit']['passed']
    return meta


def figure11():
    old='Figure_12_steady_state_mechanism';meta=audited_old(old)
    control_rows=rows(ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/tension_command_measurements.csv')
    control_summary=[]
    for length,expected_tension in [(23,2.762),(39,1.214)]:
        r=next(r for r in control_rows if float(r['length_um'])==length)
        np.testing.assert_allclose(float(r['equal_command_um_per_s']),-.04,atol=1e-15)
        np.testing.assert_allclose(float(r['equal_endpoint_displacement_um']),-1.,atol=1e-14)
        np.testing.assert_allclose(float(r['equal_front_tension_pN']),expected_tension,atol=5e-4)
        control_summary.append(dict(length_um=length,front_tension_pN=float(r['equal_front_tension_pN']),
            equal_command_um_s=float(r['equal_command_um_per_s']),
            equal_displacement_um=float(r['equal_endpoint_displacement_um']),
            frozen_profiles='Each length retains its own distinct native Rho and MCA profiles'))
    names=['endpoint_metrics','spatial_profiles','spatial_flux_components','spatial_flow_components','transport_flow']
    raw={name:rows(OLD/(old+'_'+name+'.csv')) for name in names}
    endpoint=[{k:float(v) for k,v in r.items() if k not in ['checkpoint','checkpoint_sha256']} for r in raw['endpoint_metrics']]
    threshold=endpoint[0]['rac_switch_threshold']
    flow_error=max(abs(sum(float(r[k+'_cortical_velocity_um_per_s']) for k in ['endpoint','rho','membrane'])-
                       float(r['total_cortical_velocity_um_per_s'])) for r in raw['spatial_flow_components'])
    flux_error=max(abs(sum(float(r[k]) for k in ['bound_advection','bound_diffusion','unbound_advection','unbound_diffusion'])-
                      float(r['total_flux'])) for r in raw['spatial_flux_components'])
    assert flow_error<1e-11 and flux_error<1e-12
    gradient_rows, gradient_summary = [], []
    for length in [19,21,25,31]:
        ff = [r for r in raw['spatial_flow_components'] if float(r['length_um'])==length]
        pp = [r for r in raw['spatial_profiles'] if float(r['length_um'])==length]
        assert len(ff)==len(pp)==201
        velocity = np.array([float(r['total_cortical_velocity_um_per_s']) for r in ff])
        bound = np.array([float(r['mca_bound']) for r in pp])
        vg, bg = np.gradient(velocity,length/201), np.gradient(bound,length/201)
        if length in [25,31]:
            assert np.all(vg>0) and np.all(bg>0) and np.all(velocity<0)
        gradient_summary.append(dict(length_um=length, minimum_velocity_gradient_s_inverse=float(vg.min()),
            maximum_velocity_gradient_s_inverse=float(vg.max()), positive_velocity_gradient_cell_fraction=float(np.mean(vg>0)),
            rightward_total_flow_cell_fraction=float(np.mean(velocity>0))))
        gradient_rows.extend(dict(length_um=length,cell_index=j,x_nominal_um=float(ff[j]['x_nominal_um']),
            displayed_velocity_gradient_s_inverse=float(vg[j]), basal_bound_mca_gradient_conc_per_um=float(bg[j]))
            for j in range(201))
    fig=plt.figure(figsize=(13.2,12.5))
    a=fig.add_axes([.075,.71,.245,.22]);d=fig.add_axes([.075,.385,.245,.22])
    bg=fig.add_gridspec(1,4,left=.425,right=.98,bottom=.71,top=.93,wspace=.29)
    cg=fig.add_gridspec(1,4,left=.425,right=.98,bottom=.385,top=.605,wspace=.29)
    bottom=fig.add_gridspec(1,4,left=.075,right=.98,bottom=.075,top=.265,wspace=.48)
    e,f,g,h=[fig.add_subplot(bottom[0,i]) for i in range(4)]
    for ax,letter in [(a,'A'),(d,'D'),(e,'E'),(f,'F'),(g,'G'),(h,'H')]:style(ax,letter)
    spatial=raw['spatial_profiles']
    for length,color,ls in [(17,'#aaaaaa','--'),(18,INK,':'),(19,GREEN,'-.'),(23,BLUE,'-'),(31,PURPLE,'-'),(39,ORANGE,'-')]:
        rr=[r for r in spatial if float(r['length_um'])==length];assert len(rr)==201
        a.plot([float(r['x_over_L']) for r in rr],[float(r['mca_bound']) for r in rr],color=color,ls=ls,lw=1.8,label=f'{length} µm')
    a.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Bound MCA (native density)',ylim=(.043,.115))
    a.legend(ncol=2,frameon=False,fontsize=10,loc='upper left',columnspacing=.7,handlelength=1.5)
    flux_specs=[('bound_advection','Bound advection',ORANGE,'-'),('bound_diffusion','Bound diffusion',GREEN,'-'),
        ('unbound_advection','Unbound advection','#CC79A7',':'),('unbound_diffusion','Unbound diffusion',BLUE,'-'),('total_flux','Total flux',INK,'--')]
    flow_specs=[('endpoint','Endpoint',ORANGE,'-'),('rho','Rho stress',BLUE,'-'),('membrane','Membrane',GRAY,'-'),('total','Total',INK,'--')]
    for index,length in enumerate([19,21,25,31]):
        sub=bg[0,index].subgridspec(2,1,height_ratios=[2.1,1],hspace=.13)
        main=fig.add_subplot(sub[0]);detail=fig.add_subplot(sub[1],sharex=main)
        rr=[r for r in raw['spatial_flux_components'] if float(r['length_um'])==length];assert len(rr)==202
        x=np.array([float(r['x_nominal_um']) for r in rr])
        for key,_,color,ls in flux_specs:main.plot(x,[1000*float(r[key]) for r in rr],color=color,ls=ls,lw=1.5)
        detail.plot(x,[1e6*float(r['bound_diffusion']) for r in rr],color=GREEN,lw=1.5)
        main.set(xlim=(0,length),ylim=(-2.15,.38));main.tick_params(labelbottom=False)
        detail.set(xlim=(0,length),ylim=(-3.8,.25),xlabel='x (µm)',xticks=[0,length])
        main.set_title(f'{length} µm',fontsize=11,pad=8)
        for ax in [main,detail]:style(ax);ax.axhline(0,color='#aaaaaa',lw=.65)
        if index==0:
            main.text(-.45,1.07,'B',transform=main.transAxes,weight='bold',fontsize=13)
            main.set_ylabel('MCA flux\n(10⁻³ conc. µm/s)');detail.set_ylabel('Bound diff.\n(10⁻⁶)')
        else:main.tick_params(labelleft=False);detail.tick_params(labelleft=False)
        cortex=fig.add_subplot(cg[0,index])
        rr=[r for r in raw['spatial_flow_components'] if float(r['length_um'])==length];assert len(rr)==201
        x=np.array([float(r['x_nominal_um']) for r in rr])
        for key,_,color,ls in flow_specs:cortex.plot(x,[float(r[key+'_cortical_velocity_um_per_s']) for r in rr],color=color,ls=ls,lw=1.7)
        cortex.set(xlim=(0,length),ylim=(-.045,.016),xlabel='x (µm)',xticks=[0,length]);cortex.set_title(f'{length} µm',fontsize=11,pad=8)
        style(cortex);cortex.axhline(0,color='#aaaaaa',lw=.65)
        if index==0:
            cortex.text(-.45,1.07,'C',transform=cortex.transAxes,weight='bold',fontsize=13)
            cortex.set_ylabel('Cortical velocity (µm/s)')
        else:cortex.tick_params(labelleft=False)
    handles=lambda specs:[Line2D([],[],color=color,ls=ls,lw=1.8,label=label) for _,label,color,ls in specs]
    fig.legend(handles=handles(flux_specs),loc='upper center',bbox_to_anchor=(.70,.985),ncol=3,frameon=False,fontsize=10,columnspacing=1)
    fig.legend(handles=handles(flow_specs),loc='upper center',bbox_to_anchor=(.70,.66),ncol=4,frameon=False,fontsize=10,columnspacing=.7)
    width_summary=[]
    for length,color in [(23,BLUE),(31,PURPLE),(39,ORANGE)]:
        rr=[r for r in spatial if float(r['length_um'])==length]
        q=np.array([float(r['x_over_L']) for r in rr]);rac=np.array([float(r['rac']) for r in rr]);rho=np.array([float(r['rho']) for r in rr])
        d.plot(q,rac,color=color,lw=1.9);d.plot(q,rho,color=color,lw=1.9,ls='--')
        h.plot(q,[float(r['tension']) for r in rr],color=color,lw=1.9,label=f'{length} µm')
        width_summary.append(dict(length_um=length,rac_halfcontrast_fraction=float(np.mean(rac>=.5*(rac.min()+rac.max()))),
                                  rho_halfcontrast_fraction=float(np.mean(rho>=.5*(rho.min()+rho.max())))))
    d.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Rac / Rho (native conc.)',ylim=(.54,.94))
    d.legend(handles=[Line2D([],[],color=INK,lw=1.8,label='Rac'),Line2D([],[],color=INK,ls='--',lw=1.8,label='Rho')],frameon=False,ncol=2,loc='upper right',fontsize=10)
    h.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Membrane tension (pN)',ylim=(-.45,10));h.axhline(0,color=GRAY,lw=.65);h.legend(frameon=False,fontsize=10,loc='upper left')
    lengths=[r['length_um'] for r in endpoint]
    def line(ax,key,color,ls='-',label=None):ax.plot(lengths,[r[key] for r in endpoint],color=color,ls=ls,lw=1.8,marker='o',ms=3.5,label=label)
    for ax in [e,f,g]:ax.set(xlim=(17,40),xticks=[18,31,39],xlabel='Length (µm)')
    line(e,'left_rac',BLUE);e.axhline(threshold,color='#9c2525',ls=':',lw=1.2,label='Switch midpoint')
    e.set(ylabel='Front Rac (native conc.)',ylim=(.844,.908));e.legend(frameon=False,loc='lower right',fontsize=10)
    line(f,'derived_left_command_outward_um_per_s',BLUE,label='Front');line(f,'derived_right_command_outward_um_per_s',GRAY,'--','Rear')
    f.set(ylabel='Outward command (µm/s)',ylim=(-.004,.105));f.legend(frameon=False,fontsize=10,loc='upper left')
    line(g,'left_tension_pN',BLUE,label='Front');line(g,'right_tension_pN',GRAY,'--','Rear')
    g.set(ylabel='Endpoint tension (pN)',ylim=(-.45,10));g.axhline(0,color=GRAY,lw=.65);g.legend(frameon=False,fontsize=10,loc='upper left')
    stem='Figure_11_steady_state_mechanism';artifacts=save(fig,stem)
    for name,rr in raw.items():artifacts.append(csvout(stem+'_'+name+'.csv',rr))
    artifacts.append(csvout(stem+'_normalized_halfcontrast_width.csv',width_summary))
    artifacts.append(csvout(stem+'_spatial_gradients.csv',gradient_rows))
    byL={r['length_um']:r for r in endpoint};r23,r39=byL[23],byL[39]
    command_ratio=r39['derived_left_command_outward_um_per_s']/r23['derived_left_command_outward_um_per_s']
    caption=('Figure 11. The basal MCA distribution is maintained by endpoint-driven flow, opposing Rho stress and distributed MCA turnover. '
        'A) Bound-MCA profiles at 17–39 µm. B) Full spatial bound/unbound advection and diffusion fluxes at 19, 21, 25 and 31 µm; '
        'lower axes display bound diffusion separately in 10^-6 rather than 10^-3 flux units. C) Full cortical-flow source components at the same lengths. '
        'D) Rac (solid) and Rho (dashed) profiles at 23, 31 and 39 µm. E) Rac in the old left-front cell and the smooth protrusion-switch midpoint, 0.899061648346. '
        'F) Native outward front and rear endpoint commands. G) Front and rear membrane tension. H) Full tension profiles. '
        'Front Rac rises from 0.88893 to 0.90220 between 23 and 39 µm while the native front command rises 6.20-fold. '
        'Tension and the Rac switch both affect that command. Endpoint forcing is leftward and Rho-driven flow rightward at every displayed length. '
        'Rho forcing produces rightward total flow over 87.1% of the cells at 19 µm and 58.7% at 21 µm; total flow is leftward throughout at 25 and 31 µm. '
        'Both diffusion fluxes point leftward. At 19 µm, interior advection opposes diffusion; at 25 and 31 µm, bound advection and diffusion point the same way. '
        'Local binding and uniform replacement complete the concentration balance shown in Figure 9. '
        'A and D–H use exact corrected basal checkpoints; B–C reuse audited exports after one additional native 0.5-s basal step to resolve membrane velocity. '
        'Spatial scales are native nominal coordinates. Positive flux and flow are rightward. Flow components are linear source superpositions at fixed state. '
        'No new integrations or parameter changes are made.')
    paragraphs=[
        'The pre-cue MCA profiles coexist with continuous spatial transport and turnover (Fig. 11A,B). At 19 µm, rightward interior bound-MCA advection opposes leftward diffusion. At 25 and 31 µm, bound advection points leftward throughout, in the same direction as both diffusion fluxes. The separate lower axes resolve bound diffusion using 10^-6 rather than 10^-3 flux units. Left-end outflow removes material; local exchange and uniform replacement replenish bound MCA. A nearly stationary density profile can therefore carry substantial MCA flux.',
        'Endpoint forcing supplies a leftward cortical contribution and Rho stress a rightward contribution at all four displayed lengths (Fig. 11C). At 19 µm, Rho overcomes endpoint forcing over 87.1% of the cells; at 21 µm, over 58.7%. At 25 and 31 µm, endpoint forcing wins throughout. Flow becomes less negative toward the right, so its positive gradient dilutes bound MCA throughout these two cells. Transport of the existing density profile counters that dilution: leftward flow brings in denser MCA from the right. Local exchange and uniform replacement complete the balance.',
        'The starting Rac and Rho regions broaden in absolute units as the cell lengthens (Fig. 11D). Front Rac changes only from 0.88893 at 23 µm to 0.90220 at 39 µm, but this crosses the steep protrusion-switch region centered at 0.89906 (Fig. 11E). The native front command rises 6.20-fold over that interval (Fig. 11F). The command also depends on local tension, which directly attenuates protrusion.',
        'Front tension turns downward at longer lengths while rear tension rises (Fig. 11G,H). Tension also promotes Rho activation, linking the mechanical load back to the Rac–Rho profiles. The stronger native endpoint command and changed tension are coupled outputs of this feedback. An equal-command control retains a length-dependent mechanical response: imposing the same left command of −0.040 µm/s and displacement of −1.00 µm gives front tensions of 2.762 and 1.214 pN at 23 and 39 µm. Each length’s native Rho and MCA profiles are held fixed, so geometry and inherited spatial coefficients remain combined. This control does not independently evolve the feedback.'
    ]
    methods=(caption+'\n\nThe original spatial checkpoint, face-flux and instantaneous mechanical-source exports are hash-verified against '
        'their prior manifest and native sources, then read without numerical integration. Every panel retains all 201 cells or 202 faces. '
        'Cortical components are endpoint forcing, Rho active-stress gradient and membrane drag with the same bound-MCA friction. '
        'Their pointwise sum equals total flow; four signed face fluxes sum to total MCA flux. Bound advection uses cortical face velocity '
        'and native upwind density, unbound advection uses membrane velocity; the two diffusion coefficients are native. '
        'The native front command is 0.15*[1+tanh((Rac-front-0.899061648346)/0.01)]/2/[1+T-front^2/10], in µm/s, '
        'with the reflected sign at the right endpoint. The checkpoint fields determine every command. '
        'Equal-command context reuses the separately audited control linked below; it changes commands while holding each length’s own '
        'native Rho and MCA profiles fixed. Those profiles differ across lengths, so this comparison combines geometry and inherited '
        'spatial coefficients rather than isolating geometry alone; no feedback pathway is independently evolved. '
        'The normalized halfcontrast-width CSV reports the fraction of native cells above (minimum+maximum)/2 for each field. '
        'These are descriptive widths, distinct from growth or branch thresholds in Figure 10. '
        'The spatial-gradient CSV uses centered differences between native cell centers, with one-sided differences at the ends, '
        'to describe the displayed flow and basal MCA profiles. At 25 µm the velocity gradient ranges from 2.458e-4 to 1.971e-3 s^-1; '
        'at 31 µm from 5.089e-4 to 2.768e-3 s^-1. It is positive at every cell at both lengths, while velocity is negative and '
        'the bound-MCA gradient positive. In the continuous decomposition, -B*dv/dx therefore dilutes B and -v*dB/dx opposes it. '
        'These gradient diagnostics describe the profiles; the exact finite-volume bound-MCA step budget is shown separately in Figure 9.')
    finish(stem,11,'A steep Rac protrusion switch strengthens the basal cortical drive as length increases',paragraphs,caption,methods,
        dict(original_audit=meta['audit'],full_spatial_profiles_retained=True,flow_superposition_max_error_um_s=flow_error,
             total_flux_sum_max_error=flux_error,front_command_ratio_39_over_23=command_ratio,normalized_width_summary=width_summary,
             displayed_gradient_summary=gradient_summary,equal_command_control=control_summary),artifacts,
        [['narrative_20261009/mechanics/debelly_resting_equal_commands.pdf','Equal-command mechanical control'],
         ['../../data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/tension_command_measurements.csv','Equal-command tension measurements']])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source(Path(__file__));source(Path(overview.__file__));source(Path(distance.__file__))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'xtick.labelsize':10,
        'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    for fn in [figure8,figure10,figure11]:fn();print(fn.__name__,flush=True)
    for chapter in CHAPTERS:
        meta=read(REPORT/chapter['manifest'])
        for name,digest in meta['artifacts'].items():assert sha(OUT/name)==digest
        for name,digest in meta['sources'].items():assert sha(ROOT/name)==digest
    (OUT/'chapters.json').write_text(json.dumps(CHAPTERS,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(passed=True,chapters=str(OUT/'chapters.json')),indent=2))


if __name__=='__main__':main()
