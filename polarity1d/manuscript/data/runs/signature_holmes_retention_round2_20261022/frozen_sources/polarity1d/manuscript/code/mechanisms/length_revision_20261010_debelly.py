"""Revised Figures 8--12: preserve the approved deBelly visual narrative.

Read-only rendering from audited corrected-model checkpoints and compact exports.
No native integration, original development output, or report writes.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import importlib.util
import inspect
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
CODE = Path(__file__).parent
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'revision_20261010/debelly'
SOURCES = {}
STEMS = ['Figure_8_reversal_overview', 'Figure_9_cortical_motion_mca',
         'Figure_10_MCA_distance_speed', 'Figure_11_MCA_starting_states',
         'Figure_12_steady_state_mechanism']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(register(path).read_text())


def rows(path):
    with register(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def csvout(path, data):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def module(name):
    path = register(CODE / (name + '.py'))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    outside = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > bounds.width+1 or box.y1 > bounds.height+1:
                outside.append(artist.get_text())
    assert not outside, outside
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / (stem + '.' + ext), dpi=220, facecolor='white')
    plt.close(fig)


def finish(stem, caption, audit, extra=()):
    methods = OUT / (stem + '_methods.txt')
    methods.write_text(caption + '\n')
    caption = CAPTIONS[stem]
    cp = OUT / (stem + '_caption.txt')
    cp.write_text(caption + '\n')
    ap = OUT / (stem + '_audit.json')
    audit.update(passed=True, new_integrations=0, plot_text_within_page=True)
    ap.write_text(json.dumps(audit, indent=2) + '\n')
    files = [OUT / (stem + '.' + ext) for ext in ['png', 'pdf', 'svg']] + [cp, methods, ap] + list(extra)
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
        passed=True, renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(__file__),
        artifacts={p.name: sha(p) for p in files}, sources=dict(SOURCES),
        sha256={p.name: sha(p) for p in files}, source_sha256=dict(SOURCES),
        caption=caption, legend=caption, audit=audit, new_integrations=0,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_revision_20261010_debelly.py')
    (OUT / (stem + '.manifest.json')).write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')


def figure_x():
    m = module('length_dev_figure_x_overview')
    traces, measurements, checks, arrays = m.load()
    SOURCES.update(m.SOURCES)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
                         'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig = plt.figure(figsize=(11.8, 12.8))
    gs = fig.add_gridspec(5, 4, width_ratios=[23,31,39,1.9], left=.10,
                         right=.93, bottom=.065, top=.965, hspace=.24, wspace=.18)
    axes = {}
    for row, (key, label, unit) in enumerate(zip(m.FIELDS, m.LABELS, m.UNITS)):
        visible = [tr['fields'][key][tr['time'] <= m.LIMIT] for tr in traces.values()]
        lo, hi = min(float(a.min()) for a in visible), max(float(a.max()) for a in visible)
        if key == 'rac':
            norm, cmap = TwoSlopeNorm(vmin=lo, vcenter=1., vmax=hi), 'viridis'
        elif key in ['velocity', 'tension']:
            norm, cmap = Normalize(-max(abs(lo),abs(hi)), max(abs(lo),abs(hi))), 'RdBu_r'
        else:
            norm, cmap = Normalize(lo, hi), 'viridis'
        for col, length in enumerate(m.LENGTHS):
            ax = fig.add_subplot(gs[row,col]); tr = traces[length]
            mesh = ax.pcolormesh(tr['x'], tr['time'], tr['fields'][key],
                                shading='nearest', cmap=cmap, norm=norm, rasterized=True)
            ax.set_yscale('symlog', linthresh=120, linscale=.65)
            ax.set(xlim=(0,length), ylim=(m.LIMIT,-m.PRE_CUE),
                   xticks=[0,length/2,length], yticks=[-60,0,60,120,300,600,1200])
            ax.set_yticklabels(['−60','0','60','120','300','600','1200'])
            ax.tick_params(labelsize=9, direction='out', length=3)
            if col: ax.tick_params(labelleft=False)
            else: ax.set_ylabel(label+'\nTime from cue onset (s)', fontsize=10, labelpad=8)
            ax.axhline(0,color='white',lw=1.5);ax.axhline(0,color='#444444',lw=.6)
            ax.axhline(60,color='white',ls=':',lw=1.5,alpha=.9)
            ax.axhline(tr['duration'],color='white',ls='--',lw=2)
            ax.axhline(tr['duration'],color='#222222',ls='--',lw=.7)
            if row==0: ax.set_title(f'L = {length} µm',fontsize=12,pad=12)
            if row==4: ax.set_xlabel('Reference position x (µm)',fontsize=10)
            axes[row,length] = ax
        cb = fig.colorbar(mesh,cax=fig.add_subplot(gs[row,3]))
        if key=='rac': cb.set_ticks([.6,.8,1.,2.,4.])
        cb.ax.tick_params(labelsize=8,length=2); cb.set_label(unit,fontsize=9,labelpad=6)
        cb.ax.yaxis.get_offset_text().set_fontsize(8)
    for key, letter, xy, offset in [((0,31),'a',(29,30),(-35,-10)),
        ((1,31),'b',(26,60),(-35,-10)),((2,31),'c',(24,60),(-25,-18)),
        ((3,39),'d',(29,485),(-25,28)),((4,39),'e',(10,530),(25,28))]:
        m.marker(axes[key],letter,xy,offset)
    stem=STEMS[0];save(fig,stem)
    table=OUT/(stem+'_measurements.csv');csvout(table,measurements)
    fields=OUT/(stem+'_fields.npz');np.savez_compressed(fields,**arrays)
    caption=('Figure 8. MCA retains the old spatial bias after the fast chemical response to reversal. '
        'a) Rac gains rapidly at the challenged right end. b) Rho reorients toward the old left front. '
        'c) At the common 60-s dotted line, bound MCA remains biased toward the original right rear at all three lengths. '
        'd) The 39-µm cell delays whole-cortex flow reversal until about 485 s. e) Its tension bias reverses near 530 s. '
        'Columns show 23, 31 and 39 µm; rows show Rac, Rho, bound MCA, cortical flow and membrane tension. '
        'The solid horizontal line marks cue onset and dashed lines mark withdrawal of the shortest tested 1%-strength-passing '
        'relative-20% pulses, lasting 232.244898, 285.280612 and 809.081633 s. Patch peak is 3.6. '
        'Rac first moments reverse in approximately 7–12 s and Rho in 26–50 s. MCA redistribution takes hundreds of seconds. '
        'The 39-µm mechanical delay is distinct from its early local chemical response. '
        'Each trajectory starts from its exact qualifying corrected independent-endpoint basal checkpoint; the −60-to-0-s strip repeats '
        'that checkpoint for visibility. Fields are continuous at withdrawal. Time increases downward, linearly through 120 s and logarithmically '
        'thereafter. Column widths scale with nominal length, and positions are undeformed reference coordinates. '
        'Each field shares one absolute color scale across lengths. The Rac color scale allocates half its palette below 1.0 and half above '
        '1.0 to resolve basal polarity while retaining the complete cue-induced range; the breakpoint is a display choice. Velocity and tension '
        'are signed with zero-centered scales. MCA can cross and recross the midpoint; a first crossing is not complete redistribution. '
        'All three pulses pass the final-500-s scalar Rac strength gate after 10000 s of basal recovery, beyond the displayed interval; '
        'the separate native-field tail flags are retained in the measurements.')
    finish(stem,caption,dict(checks=checks),[table,fields])


def audited_exports(folder, stem):
    meta=read(folder/'supporting/metadata'/(stem+'.manifest.json'))
    for path,digest in meta['sha256'].items():
        assert sha(folder/path)==digest,(folder/path)
    for path,digest in meta.get('source_sha256',{}).items():
        assert sha(ROOT/path)==digest, path
        SOURCES[path]=digest
    return meta


def figure_y():
    m=module('length_dev_figure_y_mass_balance')
    old='figure_Y_cortical_motion_mca_mass_balance'
    folder=REPORT/'dev/resting_state'
    meta=audited_exports(folder,old); assert meta['audit']['passed']
    datasets={key:rows(folder/'supporting/data'/(old+'_'+key+'.csv'))
              for key in ['profiles','fluxes','drawing_nodes','mass_balance','mass_balance_rates']}
    snapshots=[]
    for i,check in enumerate(meta['audit']['snapshots']):
        label=check['label']; data={key:[r for r in rr if r['snapshot']==label] for key,rr in datasets.items()}
        array=lambda key,field:np.array([float(r[field]) for r in data[key]])
        f={key:array('profiles',key) for key in ['rac','rho','mca_bound','mca_unbound','tension','displacement','velocity']}
        s=dict(label=label,time=float(data['profiles'][0]['time_after_onset_s']),f=f,audit=check,
            face_x=array('fluxes','reference_face_position_um'),
            advective=array('fluxes','bound_advective_flux'),diffusive=array('fluxes','bound_diffusive_flux'),
            unbound_advective=array('fluxes','unbound_advective_flux'),unbound_diffusive=array('fluxes','unbound_diffusive_flux'),
            total_flux=array('fluxes','total_mca_spatial_flux'),membrane_x=array('drawing_nodes','plotted_membrane_position_um'),
            cortex_x=array('drawing_nodes','plotted_cortex_marker_position_um'),
            components={k:array('profiles',k+'_velocity_um_per_s') for k in ['endpoint','rho','membrane']},
            budget={k:array('mass_balance',k) for k in ['transport','replenishment','exchange','net']})
        s['boundary']=s['advective'][[0,-1]]/f['mca_bound'][[0,-1]]
        assert np.max(abs(sum(s['components'].values())-f['velocity']))<1e-5
        assert np.max(abs(sum(s['budget'][k] for k in ['transport','replenishment','exchange'])-s['budget']['net']))<1e-12
        snapshots.append(s)
    x=array('profiles','reference_position_um'); nodes=array('drawing_nodes','reference_position_um')
    # Reuse the approved drawing and layout, reconstructing every plotted input
    # from its hash-verified compact export rather than rerunning the native model.
    source=inspect.getsource(m.main)
    draw=source[source.index("    plt.rcParams.update"):source.index("    caption=(")]
    draw='\n'.join(line[4:] if line.startswith('    ') else line for line in draw.splitlines())
    draw=draw.replace("top=.850", "top=.900")
    draw=draw.replace("fig.suptitle('Figure Y | Mechanics and MCA transport through cue withdrawal',x=.08,ha='left',y=.984,fontsize=21)", "")
    draw=draw.replace("fig.text(.08,.960,f'L = {LENGTH} µm · '+clock_description,fontsize=11,color='#555555')", "")
    draw=draw.replace("fig.text(.08,.941,'Left: deformed mechanics and native profiles. Middle: cortical flow. Right: bound + unbound MCA fluxes and local mass balance.',fontsize=11,color='#555555')", "")
    draw=draw.replace("bbox_to_anchor=(.075,.927)","bbox_to_anchor=(.075,.989)")
    draw=draw.replace(".896,.22,.007", ".944,.22,.007")
    draw=draw.replace("fig.text(.08,.018,'Generated '+generated+' · source simulations complete · current report unchanged',fontsize=8,color='#777777')", "")
    draw=draw.replace("ax.text(.04, .96, 'Positive = rightward', transform=ax.transAxes, fontsize=9, va='top', color='#555555')", "")
    draw=draw.replace("ax.text(.04,.96,'Positive = rightward',transform=ax.transAxes,fontsize=9,va='top',color='#555555')", "")
    draw=draw.replace("        rates = s['audit']['mass_budget']\n", "")
    start=draw.index("    ax.text(.5,-.13,'Next 0.5 s:")
    stop=draw.index("    budget_ax =",start)
    draw=draw[:start]+draw[stop:]
    draw=draw.replace("ax.text(.015,.02,'Membrane displacement and 10-s relative slip ×2.5',transform=ax.transAxes,fontsize=8,color='#555555')", "")
    env=dict(m.__dict__);env.update(snapshots=snapshots,x=x,nodes=nodes,p=meta['parameters'])
    exec(compile(draw,str(CODE/'length_dev_figure_y_mass_balance.py')+':retained_visual_layout','exec'),env)
    fig=env['fig'];stem=STEMS[1];save(fig,stem)
    caption=('Figure 9. Cortical motion sweeps MCA out of the domain, while uniform replacement maintains the total pool. '
        'A) Basal mechanics and Rac, Rho, bound-MCA and tension profiles in the 31-µm cell. '
        'B) Endpoint forcing, Rho stress gradients and membrane motion sum to cortical flow. '
        'C) Bound and unbound MCA fluxes and the exact local total-MCA balance. D–F) The same views 30 s into the constant-slope cue. '
        'G–I) The inherited state at withdrawal, 218.354592 s after onset. J–L) The same views 115.316713 s after withdrawal, '
        '60 s after old-front Rac crosses the smooth protrusion-switch midpoint. '
        'The shortest tested 1%-strength-passing cue has constant physical slope 0.24. Bound MCA moves with the cortex; unbound MCA moves with '
        'the membrane. Both pools diffuse. Boundary outflow removes material; the implemented spatially uniform mass correction redistributes '
        'replacement across the cell. Outflow therefore does not imply a pile-up at the downstream boundary. '
        'Springs carry membrane tension, purple axial dashpots represent cortical viscosity, and green MCA linkers couple the layers through drag. '
        'Geometry and 10-s local slip are magnified 2.5-fold; profiles retain native positions and values. Cortex markers show local 10-s motion '
        'relative to the membrane, rather than permanent particles. Flow arrows use 4 µm per 0.1 µm/s. '
        'The cortical-flow decomposition is instantaneous source superposition, not separately evolved pathway removal. Positive flow and flux are rightward. '
        'Upper flux plots show displayed-state face fluxes; lower budgets show the next native 0.5-s implicit step. Orange transport, green uniform '
        'replacement and gray sequential-exchange residual sum to the measured dashed-black total change at every cell. The gray term is a discretization '
        'residual, rather than a biological source. Pool-wide integrated sources balance boundary loss. Budget vertical limits vary by row. '
        'All plotted inputs reproduce the original audited exports without new integration. Fields are continuous at withdrawal. '
        'This pulse passes the scalar strength gate; its separate native-field tail flag is false.')
    extras=[]
    for key,rr in datasets.items():
        path=OUT/(stem+'_'+key+'.csv');csvout(path,rr);extras.append(path)
    finish(stem,caption,dict(original_audit=meta['audit'],compact_exports_hash_verified=True,
        cortical_superposition_and_local_budget_rechecked=True),extras)


def figure_z():
    m=module('length_dev_figure_z_mca_distance_speed');data,checks=m.measure();SOURCES.update(m.SOURCES)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,3,figsize=(13.2,7.3),sharex=True,sharey='row')
    fig.subplots_adjust(left=.078,right=.985,top=.85,bottom=.16,hspace=.37,wspace=.14)
    distance=[('initial_distance_to_reflected_state_um','#555555','--','s','Initial → opposite polarized centroid'),
        ('initial_distance_to_midpoint_um','#762A83',':','^','Initial → midpoint'),
        ('centroid_shift_during_cue_um','#0072B2','-','o','Net shift during cue')]
    speed=[('early_mean_speed_um_s','#009E73','-','^','Early: first 60 s'),
        ('late_mean_speed_um_s','#D55E00','-','s','Late: final 60 s'),
        ('cue_average_speed_um_s','#0072B2','-','o','Whole cue')]
    for col,geo in enumerate(m.GEOS):
        rr=sorted([r for r in data if r['condition']==geo],key=lambda r:r['length_um']);lengths=np.array([r['length_um'] for r in rr])
        ref=next(r for r in rr if r['length_um']==23)
        for row,specs in enumerate([distance,speed]):
            ax=axes[row,col];ax.axhline(0,color='#aaaaaa',lw=.7);ax.grid(axis='y',alpha=.16)
            if col==0:ax.text(-.17,1.11,'AB'[row],transform=ax.transAxes,fontsize=13,weight='bold')
            for key,color,ls,marker,label in specs:
                ax.plot(lengths,[r[key]*(1000 if row else 1) for r in rr],color=color,ls=ls,marker=marker,ms=4.3,lw=1.7)
            ax.set(xlim=(17,40),xticks=[18,23,27,31,35,39])
            if row==0:ax.set_title(m.TITLES[col],fontsize=12,pad=15);ax.set_ylim(bottom=0)
            else:
                ax.plot(lengths,ref['cue_average_speed_um_s']*lengths/23*1000,color='#777777',ls=(0,(5,3)),lw=1.2,zorder=1)
                ax.set_xlabel('Cell length (µm)')
    axes[0,0].set_ylabel('Centroid distance / net shift (µm)')
    axes[1,0].set_ylabel('Speed toward new rear\n(10⁻³ µm/s)')
    handles=lambda specs:[Line2D([],[],color=c,ls=ls,marker=mark,ms=4,label=label) for _,c,ls,mark,label in specs]
    fig.legend(handles=handles(distance),loc='upper center',bbox_to_anchor=(.53,.982),ncol=3,frameon=False,fontsize=9.3)
    fig.legend(handles=handles(speed)+[Line2D([],[],color='#777777',ls=(0,(5,3)),label='Proportional to length; anchored at 23 µm')],
               loc='lower center',bbox_to_anchor=(.53,.02),ncol=2,frameon=False,fontsize=9.3)
    stem=STEMS[2];save(fig,stem);table=OUT/(stem+'.csv');csvout(table,data)
    summaries={}
    for geo in m.GEOS:
        short,long=[next(r for r in data if r['condition']==geo and r['length_um']==length) for length in [23,39]]
        summaries[geo]={key+'_ratio_39_over_23':long[key]/short[key] for key in ['initial_distance_to_reflected_state_um','early_mean_speed_um_s','late_mean_speed_um_s','cue_average_speed_um_s']}
    caption=('Figure 10. MCA must redistribute farther in larger cells, and its cue-average centroid speed does not keep pace. '
        'A) Across fixed 3-µm, relative-20% and constant-slope cues, compare net bound-MCA centroid shift during the cue (blue), '
        'initial distance to the midpoint (purple), and initial distance to the exactly reflected polarized centroid (gray). '
        'B) Compare signed speeds over the first 60 s (green), last 60 s before withdrawal (orange), and the entire cue (blue). '
        'The dashed gray benchmark scales speed in proportion to length, anchored to each geometry’s whole-cue value at 23 µm. '
        'From 23 to 39 µm, the initial reflected-centroid distance grows 1.83-fold; whole-cue speed changes by 0.90-, 0.55- and '
        '1.09-fold for fixed, relative and gradient cues. Intermediate lengths are nonmonotonic, and gradient-driven early speed '
        'rises 1.66-fold, close to the 1.70-fold increase in length. The whole-cue result supports the larger redistribution burden '
        'without claiming that every instantaneous speed decreases with length. '
        'All twelve eligible lengths use their exact shared corrected basal checkpoints and shortest tested 1%-strength-passing pulses. '
        'The late and whole-cue windows therefore have different durations across lengths. Centroids are mass-normalized bound-MCA '
        'distributions in nominal reference coordinates; exchange and changing bound amount contribute to their movement. These are '
        'distribution speeds, rather than tracked molecular speeds or cortical velocities. Initial benchmark distances are not remaining '
        'distances at withdrawal. Reflected MCA is a geometric reference; the Rac scalar gate does not require MCA to match it. '
        'Field-tail flags and complete measurements remain in the CSV; no length is excluded because of its field-tail flag.')
    finish(stem,caption,dict(checks=checks,summary=summaries),[table])


def figure_alpha():
    m=module('length_dev_figure_alpha')
    m.OUT=OUT;m.STEM=STEMS[4];m.output=lambda name:OUT/name
    body=inspect.getsource(m.main).replace("    fig.text(.065,.975,'Figure α',fontsize=14,weight='bold',va='top')", "")
    exec(compile(body,str(CODE/'length_dev_figure_alpha.py')+':retained_visual_layout','exec'),m.__dict__)
    m.main();SOURCES.update(m.SOURCES)
    stem=STEMS[4];old=json.loads((OUT/(stem+'.manifest.json')).read_text()); audit=old['audit']
    caption=('Figure 12. Length amplifies the native protrusive drive that organizes steady MCA and membrane tension. '
        'A) Bound-MCA profiles from 17 to 39 µm. B) Full spatial MCA-flux components at 19, 21, 25 and 31 µm; the lower axes '
        'magnify bound diffusion by 1000-fold. C) Full spatial endpoint, Rho-gradient and membrane-motion contributions to cortical flow '
        'at the same lengths. D) Rac (solid) and Rho (dashed) occupy approximately fixed fractions of the domain at 23, 31 and 39 µm. '
        'E) Front Rac approaches the steep protrusion switch and first exceeds its 0.899061648346 midpoint at the sampled 37-µm length. '
        'F) The native outward front command rises strongly while the rear command stays small. G) Front tension turns down at long lengths '
        'as rear tension rises. H) Full spatial tension profiles expose the growing rear load. '
        'Small changes in endpoint Rac produce large changes in the smooth command law. Increasing endpoint drive strengthens MCA advection '
        'against opposing diffusion and Rho-driven flow. Greater length distributes endpoint displacement over more membrane, allowing '
        'stronger native cortical drive; the resulting tension is concentrated toward the original unpolarized rear. '
        'The command is 6.20-fold greater at 39 than at 23 µm; the reflected MCA-centroid distance is 1.83-fold greater. '
        'A and D–H use exact qualifying corrected basal checkpoints. B–C reuse audited spatial exports after one additional native 0.5-s '
        'no-cue step, resolving membrane velocity; this revision adds no integration. Flux and velocity are positive rightward. B and C '
        'share vertical scales across lengths. Bound MCA uses cortical velocity and unbound MCA uses membrane velocity. Local exchange '
        'and uniform replacement lie outside the spatial-flux sum. Flow components are instantaneous source superpositions, not separately '
        'evolved pathway ablations. The positions are nominal coordinates; concentrations retain native units. Figure 11 separately '
        'locates the sampled polarized-branch loss and homogeneous instability, neither of which is defined by equal flux magnitudes.')
    extras=[p for p in OUT.glob(stem+'*.csv')]
    finish(stem,caption,audit,extras)


def figure_bridge():
    folder=REPORT/'dev/mca_threshold_profiles'
    audited_exports(folder,'mca_threshold_profiles')
    profile=rows(folder/'supporting/data/mca_threshold_spatial_profiles.csv')
    length_summary=rows(folder/'supporting/data/mca_threshold_length_summary.csv')
    growth=read(ROOT/'data/runs/debelly_corrected_zero_crossing_20261008/linear.json')
    assert read(ROOT/'data/runs/debelly_corrected_zero_crossing_20261008/independent_audit.json')['passed']
    transport=rows(OUT/(STEMS[4]+'_transport_flow.csv'))
    linear_x=np.array([r['length'] for r in growth]);linear_y=np.array([r['rate'] for r in growth])
    i=np.flatnonzero(linear_y[:-1]*linear_y[1:]<0);assert len(i)==1
    onset=float(np.interp(0,linear_y[i[0]:i[0]+2],linear_x[i[0]:i[0]+2]))
    flat=[float(r['length_um']) for r in length_summary if float(r['bound_mca_range'])<1e-4]
    asymmetric=[float(r['length_um']) for r in length_summary if float(r['bound_mca_range'])>=1e-4]
    loss=[max(flat),min(asymmetric)];assert loss==[17.9,17.925]
    fig=plt.figure(figsize=(13.2,7.4))
    a=fig.add_axes([.07,.18,.43,.70]);b=fig.add_axes([.61,.59,.35,.29]);c=fig.add_axes([.61,.15,.35,.29])
    for ax,letter in [(a,'A'),(b,'B'),(c,'C')]:
        ax.text(-.12,1.06,letter,transform=ax.transAxes,fontsize=13,weight='bold');ax.grid(alpha=.13)
    colors={17:'#aaaaaa',18:'#252525',19:'#009E73',21:'#56B4E9',23:'#0072B2',31:'#762A83',39:'#D55E00'}
    for length,color in colors.items():
        rr=[r for r in profile if float(r['length_um'])==length];assert len(rr)==201
        a.plot([float(r['x_fraction']) for r in rr],[float(r['bound_mca']) for r in rr],color=color,
               ls='--' if length==17 else '-',lw=2,label=f'{length} µm')
    a.axhline(.068,color='#9D245D',ls=':',lw=1.2,label='Rac-feedback midpoint')
    a.set(xlim=(0,1),xticks=[0,.5,1],xlabel='Reference position x/L',ylabel='Bound MCA (native density)')
    a.legend(frameon=False,ncol=3,loc='upper left',fontsize=9)
    # Regimes come from independently audited branch continuation and linear
    # stability, not from flux equality. Bands intentionally retain that distinction.
    for ax in [b,c]:
        ax.axvspan(17,loss[0],color='#dddddd',alpha=.5)
        ax.axvspan(loss[1],onset,color='#E69F00',alpha=.14)
        ax.axvspan(onset,40,color='#009E73',alpha=.09)
        ax.axvspan(*loss,color='#777777',alpha=.5)
        ax.axvline(onset,color='#009E73',lw=1.1,ls=':')
        ax.set(xlim=(17,40),xticks=[18,23,31,39],xlabel='Length (µm)')
    b.plot(linear_x,1000*linear_y,color='#252525',lw=2)
    b.axhline(0,color='#777777',lw=.8)
    b.set(ylabel='Homogeneous growth rate\n(10⁻³ s⁻¹)',ylim=(-.40,.80))
    rr=sorted(transport,key=lambda r:float(r['length_um']))
    for key,label,color in [('rms_bound_advective_flux','Bound advection','#D55E00'),
                           ('rms_combined_diffusive_flux','Bound + unbound diffusion','#0072B2')]:
        c.plot([float(r['length_um']) for r in rr],[1000*float(r[key]) for r in rr],color=color,lw=2,label=label)
    c.set(ylabel='Interior RMS MCA flux\n(10⁻³ conc. µm/s)')
    c.legend(frameon=False,loc='upper left',fontsize=9)
    # The narrow bistable band is labeled by letters so prose can explain it
    # clearly without stacking unreadable wording over two microns of x range.
    for xx,label,color in [(17.45,'i','#888888'),(18.7,'ii','#B27600'),(30,'iii','#007A59')]:
        b.text(xx,1.04,label,transform=b.get_xaxis_transform(),ha='center',fontsize=10,weight='bold',color=color)
    stem=STEMS[3];save(fig,stem)
    caption=(f'Figure 11. MCA polarity persists below the length at which homogeneous cells become spontaneously polarizable. '
        'A) Verified bound-MCA profiles span the nearly homogeneous 17-µm state and polarized 18–39-µm states; the dotted magenta '
        'line marks the MCA midpoint in the Rac-feedback law. B) The independently computed growth rate of an infinitesimal asymmetric '
        f'perturbation changes sign near {onset:.2f} µm. Region i has no retained polarized branch in this continuation; region ii '
        'supports a finite polarized state while the homogeneous state is linearly stable; region iii has a linearly unstable homogeneous '
        'state. The gray interval marks the sampled polarized-branch loss between 17.9 and 17.925 µm. The orange band therefore describes '
        'finite perturbation requirements rather than an absence of polarized states. C) Interior RMS bound-MCA advection and the combined '
        'bound-plus-unbound diffusive flux at the verified steady states. Their opposing contributions are comparable near the short-cell '
        'regime; advection becomes much larger as length and protrusive drive increase. Flux magnitudes measure the transport balance; '
        'their equality does not define either boundary in B. '
        'Original qualifying checkpoints supply the retained 18–39-µm branch; fine lengths are continued from the exact corrected 18-µm '
        'checkpoint. Convergence and holdout audits pass. C uses the native steady velocity, upwind bound-MCA face donors and the two '
        'native diffusion constants, excluding boundary faces when calculating RMS. Bound and unbound diffusion are added as signed '
        'face fluxes before calculating their RMS. The homogeneous growth calculation and polarized-branch continuation are separate '
        'measurements. Uniform replacement and local exchange complete the MCA balance as shown in Figure 9. The new figure reuses '
        'audited profiles and calculations without new integration.')
    lp=OUT/(stem+'_linear_growth.csv');csvout(lp,growth)
    pp=OUT/(stem+'_profiles.csv');csvout(pp,profile)
    tp=OUT/(stem+'_flux_summary.csv');csvout(tp,transport)
    finish(stem,caption,dict(sampled_branch_loss_um=loss,homogeneous_instability_um=onset,
        regimes_from_independent_growth_and_continuation=True,flux_equality_not_used_to_locate_bifurcations=True),[lp,pp,tp])


def main():
    OUT.mkdir(parents=True,exist_ok=True);register(Path(__file__))
    for fn in [figure_x,figure_y,figure_z,figure_alpha,figure_bridge]:
        fn();print(fn.__name__,flush=True)
    manifests=[read(OUT/(stem+'.manifest.json')) for stem in STEMS]
    for stem,meta in zip(STEMS,manifests):
        for name,digest in meta['artifacts'].items():assert sha(OUT/name)==digest
        for name,digest in meta['sources'].items():assert sha(ROOT/name)==digest
    (OUT/'NARRATIVE.md').write_text(NARRATIVE)
    print(json.dumps(dict(passed=True,figures=STEMS,generated_at=datetime.now(timezone.utc).isoformat()),indent=2))


CAPTIONS = {
STEMS[0]: ('Figure 8. MCA retains the old spatial bias after the fast chemical response to reversal. '
    'a) Rapid Rac gain at the challenged right end. b) Rho reorientation toward the old left front. '
    'c) At 60 s, bound MCA remains biased toward the original right rear at every length. '
    'd) Whole-cortex flow reversal is delayed to about 485 s at 39 µm. e) Its tension bias reverses near 530 s. '
    'Columns show 23, 31 and 39 µm; rows show Rac, Rho, bound MCA, cortical flow and tension. '
    'Solid lines mark cue onset, dotted lines the common 60-s comparison, and dashed lines withdrawal of the shortest tested '
    '1%-strength-passing relative-20% pulses (232.244898, 285.280612 and 809.081633 s; peak 3.6). '
    'Rac first moments reverse in approximately 7–12 s and Rho in 26–50 s; MCA redistribution takes hundreds of seconds. '
    'The 39-µm mechanical delay follows an early local chemical response. The pre-cue strip repeats the exact corrected basal checkpoint '
    'for visibility. Time increases downward, linearly through 120 s and logarithmically thereafter. Column widths scale with nominal '
    'length; each field shares one absolute color scale. Rac uses a piecewise-linear mapping with half the palette below 1.0 and half '
    'above. Fields remain continuous at withdrawal. All three pulses pass the Rac scalar gate after 10000 s of recovery; field-tail flags '
    'are retained separately. MCA can cross and recross the midpoint.'),
STEMS[1]: ('Figure 9. Cortical motion sweeps MCA out of the domain, while uniform replacement maintains the total pool. '
    'A) Basal mechanics and Rac, Rho, bound-MCA and tension profiles in a 31-µm cell. B) Endpoint, Rho-gradient and membrane-motion '
    'contributions sum to cortical flow. C) MCA fluxes and the local bound-plus-unbound balance. D–F) The same views 30 s into the cue. '
    'G–I) The inherited state at withdrawal, 218.354592 s after onset. J–L) The same views 115.316713 s after withdrawal, 60 s after '
    'old-front Rac crosses the smooth protrusion-switch midpoint. This is the shortest tested 1%-strength-passing constant-slope pulse '
    '(slope 0.24). Bound MCA moves with the cortex and unbound MCA with the membrane; both diffuse. Endpoint outflow removes MCA, '
    'and replacement is uniform across the cell rather than accumulated at the downstream boundary. Springs carry tension; purple '
    'dashpots represent cortical viscosity; green MCA linkers carry drag. Geometry and 10-s local slip are magnified 2.5-fold; profiles '
    'retain native values. Flow arrows use 4 µm per 0.1 µm/s. Upper flux plots show instantaneous diagnostics. Lower budgets show the '
    'next native 0.5-s step: transport, uniform replacement and the numerical sequential-exchange residual sum to the measured total '
    'change. Positive flow and flux are rightward. Flow decomposition is instantaneous source superposition. The pulse passes the scalar '
    'strength gate; its separate field-tail flag is false.'),
STEMS[2]: ('Figure 10. MCA must redistribute farther in larger cells, and its cue-average centroid speed does not keep pace. '
    'A) Across fixed 3-µm, relative-20% and constant-slope cues, compare net bound-MCA centroid shift during the cue (blue), initial '
    'distance to the midpoint (purple), and initial distance to the reflected polarized centroid (gray). B) Signed speeds over the '
    'first 60 s (green), last 60 s before withdrawal (orange), and whole cue (blue). Dashed gray benchmarks scale speed in proportion '
    'to length, anchored to each geometry’s whole-cue value at 23 µm. From 23 to 39 µm, the reflected-centroid distance grows 1.83-fold; '
    'whole-cue speed changes 0.90-, 0.55- and 1.09-fold for fixed, relative and gradient cues. Intermediate lengths are nonmonotonic; '
    'gradient early speed rises 1.66-fold, close to the 1.70-fold length increase. All twelve eligible lengths use exact shared corrected '
    'basal checkpoints and shortest tested 1%-strength-passing pulses. Late and whole-cue windows therefore differ across lengths. '
    'Centroids describe mass-normalized bound-MCA distributions in nominal coordinates; exchange and changing bound amount contribute '
    'to their movement. Initial benchmark distances are not remaining distances at withdrawal, and the reflected MCA state is a '
    'geometric reference rather than a Rac-gate requirement. Field-tail flags remain separate.'),
STEMS[3]: ('Figure 11. MCA polarity persists below the length at which homogeneous cells become spontaneously polarizable. '
    'A) Verified bound-MCA profiles from the nearly homogeneous 17-µm state through the polarized 18–39-µm states. Magenta marks the '
    'MCA midpoint in the Rac-feedback law. B) The homogeneous state’s asymmetric growth rate changes sign near 19.53 µm. Region i '
    'has no retained polarized branch in the continuation; region ii supports a polarized state while homogeneous cells resist '
    'infinitesimal perturbations; region iii has a linearly unstable homogeneous state. Gray marks sampled polarized-branch loss '
    'between 17.9 and 17.925 µm. C) Interior RMS bound-MCA advection and combined bound-plus-unbound diffusion at verified steady '
    'states. Their opposing contributions are comparable near the short-cell regime; advection grows strongly at larger lengths. '
    'Flux magnitudes illustrate transport balance; their equality does not define either boundary in B. Fine lengths continue from '
    'the exact corrected 18-µm checkpoint, and convergence and holdout audits pass. RMS excludes boundary faces and adds the two signed '
    'diffusive fluxes before calculating their magnitude. The growth-rate and branch-continuation calculations locate distinct transitions. '
    'Uniform replacement and local exchange complete the MCA balance as shown in Figure 9.'),
STEMS[4]: ('Figure 12. Length amplifies the native protrusive drive that organizes steady MCA and membrane tension. '
    'A) Bound-MCA profiles from 17 to 39 µm. B) Spatial MCA-flux components at 19, 21, 25 and 31 µm; lower axes magnify bound diffusion '
    '1000-fold. C) Spatial endpoint, Rho-gradient and membrane-motion contributions to cortical flow. D) Rac (solid) and Rho (dashed) '
    'profiles occupy approximately fixed fractions of the 23-, 31- and 39-µm domains. E) Front Rac approaches the steep protrusion switch, '
    'first exceeding its 0.899061648346 midpoint at the sampled 37-µm length. F) Native outward front command rises while rear command '
    'stays small. G) Front tension turns down at long lengths as rear tension rises. H) Full tension profiles show the growing rear load. '
    'The command grows 6.20-fold from 23 to 39 µm, against a 1.83-fold rise in reflected MCA-centroid distance. More membrane distributes '
    'endpoint displacement, allowing stronger native cortical drive; the resulting tension grows toward the original unpolarized rear. '
    'A and D–H use exact corrected basal checkpoints. B–C reuse audited exports after one native 0.5-s no-cue step and share vertical '
    'scales across lengths. Positive flux and flow are rightward. Local exchange and uniform replacement lie outside the flux sum. '
    'Flow components are instantaneous source superpositions. Figure 11 separately locates the polarized-branch loss and homogeneous '
    'instability; equal flux magnitudes define neither transition.')
}


NARRATIVE='''# Revised deBelly sequence: report prose and evidence

Figure 8 opens with the visible reversal. The challenged right end gains Rac within seconds (a), and Rho turns toward the left end (b). At the shared 60-s line, bound MCA still carries the original right-rear bias in every cell (c). MCA is the slow spatial memory. The local chemical response begins quickly at 39 µm too, but its whole-cortex flow and tension redistribution are delayed to approximately 485 and 530 s (d–e). This long-cell exception motivates examining the mechanics rather than assuming a uniform speed.

Figure 9 explains the mechanism that moves MCA. Follow A–C from the resting cell to D–F during the cue, G–I at withdrawal, and J–L after old-front Rac falls through the protrusion-switch midpoint. MCA linkers couple membrane springs to the viscous cortex. Endpoint commands dominate the displayed flow; Rho gradients and membrane motion alter its spatial pattern. Bound MCA is swept along the cortex and out through an endpoint. The solver then replaces lost MCA uniformly across the domain. The uniform source and exchange balance can oppose local clearance even when advective flux is large; this is redistribution rather than downstream boundary accumulation.

Figure 10 compares the growing spatial burden with the achieved distribution speed. From 23 to 39 µm the MCA distance to the reflected polarized centroid increases 1.83-fold (A). Whole-cue centroid speed changes 0.90-fold for a fixed patch, 0.55-fold for a relative patch, and 1.09-fold for a gradient (B). MCA must redistribute farther without a matching increase in average speed. The gradient early speed is a useful exception, rising 1.66-fold; later and whole-cue behavior therefore matter. These signed centroid speeds describe changing distributions, including exchange, rather than tracked molecular velocities.

Figure 11 explains the inherited MCA profiles. The polarized branch already exists at 18 µm, although a small perturbation of the homogeneous state decays until approximately 19.53 µm (A–B). Between the branch-loss bracket of 17.9–17.925 µm and the homogeneous-instability point, polarization requires a finite perturbation. The short polarized profiles arise while advective and diffusive MCA fluxes are comparable (C); longer cells strengthen advection. These flux magnitudes illustrate transport organization, while the separately measured branch and growth-rate calculations define polarizability.

Figure 12 follows that increased advection back to endpoint mechanics. MCA profiles and full spatial fluxes show the growing cortical transport (A–C). Rac and Rho occupy approximately fixed fractions of the cell (D), but the small rise in front Rac reaches the steep endpoint switch (E). The outward command rises 6.20-fold from 23 to 39 µm (F). More membrane spreads the endpoint displacement, allowing the stronger native command while front tension turns down (G). The resulting spatial tension grows toward the old unpolarized rear (H). The starting MCA distribution is therefore an output of coupled chemical and mechanical feedback, rather than an arbitrary initial profile.
'''


if __name__=='__main__':
    main()
