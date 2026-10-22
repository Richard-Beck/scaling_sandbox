"""Draft-aligned Figures 9 and 10 from retained corrected native trajectories.

No integration is executed. The audited spatial curves and spring/dashpot views
are retained, with mechanics preceding bound-MCA budgets. Common-time outbound
flux and endpoint changes are separated from length-dependent cue durations.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import length_revision_20261014_figure9 as spatial

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
REV = REPORT / 'revision_20261014_draft_alignment'
OUT = REV / 'figures'
DRAFT = REPORT / 'dev/mca_endpoint_change_drafts_20261014'
BUDGET = ROOT / 'data/runs/debelly_mca_speed_decomposition_20261010'
PRIOR_SEQUENCE = REPORT / 'revision_20261014/visual_palette/debelly'
TIMES = [30., 60., 120., 180.]
TIME_COLORS = dict(zip(TIMES, plt.get_cmap('viridis')(np.linspace(.12, .85, 4))))
LENGTH_COLORS = {23:'#332288', 31:'#882255', 35:'#44AA99', 37:'#88CCEE', 39:'#999933'}
GEO_COLORS = {'fixed_physical':'#0072B2', 'fixed_relative':'#D55E00', 'fixed_slope':'#8B3F82'}
GEO_NAMES = {'fixed_physical':'Fixed 3 µm', 'fixed_relative':'Relative 20%', 'fixed_slope':'Slope 0.24'}
SOURCES = {}


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


def csvwrite(name, records):
    path = OUT/name
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    return path


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def text_bounds(fig):
    fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    outside=[]
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_visible() and obj.get_text():
            bb=obj.get_window_extent(renderer)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.width+1 or bb.y1 > fig.bbox.height+1:
                outside.append(obj.get_text())
    assert not outside, outside


def save(fig, stem):
    text_bounds(fig)
    result=[]
    for ext in ['png','pdf','svg']:
        path=OUT/(stem+'.'+ext); fig.savefig(path, dpi=220, facecolor='white'); result.append(path)
    plt.close(fig)
    return result


def panel(ax, letter):
    ax.text(-.025,1.045,letter,transform=ax.transAxes,weight='bold',fontsize=13,va='bottom')
    ax.grid(axis='y',alpha=.15,lw=.65); ax.tick_params(labelsize=10)


def load_endpoint():
    manifest=read(DRAFT/'supporting/metadata/endpoint_bound_MCA.manifest.json')
    assert manifest['passed'] and manifest['producer']=='finished'
    for name,digest in manifest['sha256'].items():
        assert sha(source(DRAFT/name))==digest, name
    for name,digest in manifest['source_sha256'].items():
        assert sha(source(ROOT/name))==digest, name
    changes=rows(DRAFT/'supporting/data/endpoint_bound_MCA_change_summary.csv')
    outflux=rows(DRAFT/'supporting/data/outward_advected_bound_MCA_summary.csv')
    time_rows=rows(DRAFT/'supporting/data/endpoint_bound_MCA_native_time_series.csv')
    manifest=read(BUDGET/'manifest.json')
    assert manifest['geometry']=='fixed_relative'
    assert read(BUDGET/'status.json')['audit_passed']
    traces={}; basal=[]; audits=[]; endpoints=[]
    for slot,L in enumerate(manifest['lengths_um']):
        path=BUDGET/'workers'/f'{slot:03d}.npz'; meta=read(path.with_suffix('.json'))
        assert meta['passed'] and meta['parameters']['center_boundary_motion'] is False
        assert sha(source(path))==meta['npz_sha256']
        cp=source(ROOT/meta['checks'][0]['checkpoint'])
        with np.load(path) as z:
            t=z[f'L{L}_time'].copy(); b=z[f'L{L}_bound'][:,[0,-1]].copy()
            rac=z[f'L{L}_rac'][:,[0,-1]].copy(); x=z[f'L{L}_x'].copy()
        with np.load(cp) as z:
            np.testing.assert_array_equal(b[0],z['state_mca_bound'][[0,-1]])
            np.testing.assert_array_equal(rac[0],z['state_rac'][[0,-1]])
        exported=[r for r in time_rows if int(r['length_um'])==L]
        np.testing.assert_array_equal(t,[float(r['elapsed_time_s']) for r in exported])
        np.testing.assert_array_equal(b[:,0],[float(r['left_bound_MCA']) for r in exported])
        np.testing.assert_array_equal(b[:,1],[float(r['right_bound_MCA']) for r in exported])
        D=float(meta['checks'][0]['cue_duration_s'])
        assert t[0]==0 and t[-1]==D and D>180
        traces[L]=dict(time=t,bound=b,rac=rac,duration=D,x=x,field_tail=meta['checks'][0]['native_field_tail_pass'])
        basal.append(dict(length_um=L,left_bound_MCA=float(b[0,0]),right_bound_MCA=float(b[0,1]),
            left_Rac=float(rac[0,0]),right_Rac=float(rac[0,1]),left_sample_x_um=float(x[0]),right_sample_x_um=float(x[-1]),
            cue_duration_s=D,native_field_tail_pass=traces[L]['field_tail']))
        for i in range(len(t)):
            endpoints.append(dict(length_um=L,time_after_onset_s=float(t[i]),cue_duration_s=D,
                old_left_bound_MCA=float(b[i,0]),new_right_bound_MCA=float(b[i,1]),
                old_left_Rac=float(rac[i,0]),new_right_Rac=float(rac[i,1])))
        audits.append(dict(length_um=L,archive_hash_verified=True,exact_checkpoint_verified=True,
                           endpoint_draft_CSV_matches_native_arrays_exactly=True,common_times_before_withdrawal=True))
    duration_manifest=read(PRIOR_SEQUENCE/'Figure_10_MCA_progress_and_regimes.manifest.json')
    duration_path=source(PRIOR_SEQUENCE/'Figure_10_MCA_progress_and_regimes.csv')
    assert duration_manifest['passed']
    assert sha(duration_path)==duration_manifest['artifacts'][duration_path.name]
    duration_rows=rows(duration_path)
    return changes,outflux,traces,basal,endpoints,duration_rows,audits


def figure9(x,p,snapshots,trace,outflux):
    """Reposition retained artists; replace only the old global-pool panel.

    This preserves every chemical, flow-source and local-budget x/y array and
    every data-derived mechanics drawing. The original renderer saves into the
    isolated output directory through a single controlled save hook.
    """
    original_save=matplotlib.figure.Figure.savefig
    modified=[]; artist_audit=[]; axes_layout=[]
    def hook(fig,*args,**kwargs):
        if not modified:
            modified.append(True)
            h=float(fig.get_size_inches()[1]); assert h==16.7 and len(fig.axes)==15
            saved_arrays=[[(np.asarray(z.get_xdata()).copy(),np.asarray(z.get_ydata()).copy()) for z in ax.lines]
                          for ax in fig.axes[:12]]
            spatial_axes={}
            for ax in fig.axes[:12]:
                texts=[z for z in ax.texts if z.get_text() in list('ABCDEFGHIJKL')]
                assert len(texts)==1
                spatial_axes[texts[0].get_text()]=(ax,texts[0])
            for old,new in zip('GHI','JKL'):
                ax,label=spatial_axes[old]; bb=ax.get_position(); ax.set_position([bb.x0,3.60/h,bb.width,bb.height]); label.set_text(new)
            for old,new in zip('JKL','GHI'):
                ax,label=spatial_axes[old]; bb=ax.get_position(); ax.set_position([bb.x0,8.35/h,bb.width,bb.height]); label.set_text(new)
            # Keys follow their own row, with physical plotting sizes retained.
            fig.legends[1].set_bbox_to_anchor((.535,7.07/h))
            fig.legends[2].set_bbox_to_anchor((.535,8.29/h))
            for ax in fig.axes[9:12]:ax.set_title('')
            for ax in fig.axes[12:14]:
                bb=ax.get_position(); ax.set_position([bb.x0,7.68/h,bb.width,bb.height])
            assert len(fig.texts)==1
            fig.texts[0].set_position((.535,7.08/h))
            ax=fig.axes[-1]; ax.clear(); bb=ax.get_position()
            ax.set_position([bb.x0,.65/h,bb.width,bb.height]); spatial.panel(ax,'M')
            for t in TIMES:
                selected=[r for r in outflux if r['sample']==f'{int(t)}_s']
                ax.plot([float(r['length_um']) for r in selected],[float(r['both_ends_bound_loss_amount']) for r in selected],
                        color=TIME_COLORS[t],marker='o',ms=4,lw=1.9,label=f'{int(t)} s')
            ax.set(xlim=(17,40),xticks=[18,23,27,31,35,39],ylim=(0,None),xlabel='Nominal cell length (µm)',
                   ylabel='Cumulative outward bound MCA\n(native concentration × µm)')
            ax.legend(title='Common time from cue onset',frameon=False,ncol=4,loc='lower center',
                      bbox_to_anchor=(.5,1.04),fontsize=10,title_fontsize=10)
            for old,ax in zip(saved_arrays,fig.axes[:12]):
                assert len(old)==len(ax.lines)
                for (xx,yy),line in zip(old,ax.lines):
                    np.testing.assert_array_equal(xx,line.get_xdata()); np.testing.assert_array_equal(yy,line.get_ydata())
            for ax in fig.axes:
                for line in ax.lines:
                    artist_audit.append(dict(label=line.get_label(),color=matplotlib.colors.to_hex(line.get_color()),
                        x_sha256=hashlib.sha256(np.asarray(line.get_xdata(),dtype=float).tobytes()).hexdigest(),
                        y_sha256=hashlib.sha256(np.asarray(line.get_ydata(),dtype=float).tobytes()).hexdigest()))
                labels=[z.get_text() for z in ax.texts if z.get_text() in list('ABCDEFGHIJKLM')]
                if labels: axes_layout.append(spatial.layout_record(ax,labels[0],fig))
            text_bounds(fig)
        return original_save(fig,*args,**kwargs)
    spatial.OUT=OUT
    matplotlib.figure.Figure.savefig=hook
    try:
        spatial.render(x,p,snapshots,trace)
    finally:
        matplotlib.figure.Figure.savefig=original_save
    return dict(passed=True,original_spatial_and_mechanical_arrays_unchanged=True,
        figure_width_in=13.1,figure_height_in=16.7,axes=axes_layout,curves=artist_audit,
        common_times_s=TIMES,plot_text_within_page=True,new_integrations=0)


def figure10(changes,traces,basal,durations):
    fig,axes=plt.subplots(3,2,figsize=(13.1,13.0))
    fig.subplots_adjust(left=.105,right=.985,bottom=.062,top=.943,hspace=.67,wspace=.28)
    for ax,label in zip(axes.flat,'ABCDEF'): panel(ax,label)
    ax=axes[0,0]
    for geo in GEO_NAMES:
        selected=[r for r in durations if r['condition']==geo]
        for style,flag in [('-',True),('none',False)]:
            rr=[r for r in selected if (r['native_field_tail_pass']=='True')==flag]
            if style=='-':
                ax.plot([float(r['length_um']) for r in selected],[float(r['cue_duration_s']) for r in selected],
                        color=GEO_COLORS[geo],lw=1.8,label=GEO_NAMES[geo])
            ax.plot([float(r['length_um']) for r in rr],[float(r['cue_duration_s']) for r in rr],ls='',
                    marker='o',ms=5,mfc=GEO_COLORS[geo] if flag else 'white',mec=GEO_COLORS[geo],color=GEO_COLORS[geo])
    ax.set(xlabel='Nominal cell length (µm)',ylabel='Shortest tested strength-passing cue (s)',xlim=(17,40),ylim=(0,None),
           xticks=[18,23,27,31,35,39])
    ax.legend(frameon=False,fontsize=10,loc='lower left',bbox_to_anchor=(0,1.04),ncol=3)
    ax=axes[0,1]
    for key,color,label in [('left_bound_MCA','#0072B2','Old left pole'),('right_bound_MCA','#D55E00','Right cue end')]:
        ax.plot([r['length_um'] for r in basal],[r[key] for r in basal],color=color,marker='o',ms=4,lw=1.8,label=label)
    ax.set(xlabel='Nominal cell length (µm)',ylabel='Pre-cue bound MCA (native concentration)',xlim=(17,40),xticks=[18,23,27,31,35,39])
    ax.legend(frameon=False,fontsize=10,loc='lower left',bbox_to_anchor=(0,1.04),ncol=2)
    endpoint_limits=[]
    for ax,side,title in [(axes[1,0],'left','Old left pole'),(axes[1,1],'right','Right cue end')]:
        for t in TIMES:
            rr=[r for r in changes if r['sample']==f'{int(t)}_s']
            values=[float(r[f'{side}_net_change_conc']) for r in rr];endpoint_limits.extend(values)
            ax.plot([float(r['length_um']) for r in rr],values,color=TIME_COLORS[t],marker='o',ms=4,lw=1.8,label=f'{int(t)} s')
        ax.axhline(0,color='#999999',lw=.8)
        ax.set(xlim=(17,40),xticks=[18,23,27,31,35,39],xlabel='Nominal cell length (µm)',
               ylabel='Bound-MCA change from pre-cue\n(native concentration)',title=title)
        ax.legend(title='Common time from cue onset',frameon=False,ncol=4,loc='lower left',
                  bbox_to_anchor=(0,1.16),fontsize=10,title_fontsize=10,columnspacing=1.)
    pad=.08*(max(endpoint_limits)-min(endpoint_limits))
    for ax in axes[1]:ax.set_ylim(min(endpoint_limits)-pad,max(endpoint_limits)+pad)
    threshold=.899061648345644
    for ax,key,label in [(axes[2,0],'bound','Old-end bound MCA (native concentration)'),
                         (axes[2,1],'rac','Old-end Rac (native concentration)')]:
        for L,color in LENGTH_COLORS.items():
            tr=traces[L]
            ax.plot(tr['time'],tr[key][:,0],color=color,lw=1.9,label=f'{L} µm')
            ax.plot(tr['duration'],tr[key][-1,0],marker='s',ms=5,color=color,ls='')
        ax.set(xlim=(0,850),xticks=[0,180,400,600,800],xlabel='Time after cue onset (s)',ylabel=label)
        ax.legend(title='Cell length',frameon=False,ncol=3,loc='lower left',bbox_to_anchor=(0,1.035),
                  fontsize=10,title_fontsize=10,columnspacing=1.2)
    axes[2,1].axhline(threshold,color='#777777',ls=':',lw=1.1)
    fig.text(.545,.018,'E–F: right relative-20% cue; squares mark withdrawal. Dotted Rac line: protrusion-switch midpoint.',
             ha='center',fontsize=10)
    audit=dict(passed=True,new_integrations=0,common_times_s=TIMES,endpoint_change_panels_share_concentration_limits=True,
        all_twelve_eligible_lengths_displayed_in_A_to_D=True,trajectory_lengths_um=list(LENGTH_COLORS),
        trajectory_curves_end_at_own_withdrawal=True,protrusion_midpoint_not_reversal_gate=True,
        separate_native_field_flags_in_A=True,figure_width_in=13.1,figure_height_in=13.,plot_text_within_page=True)
    return save(fig,'Figure_10_endpoint_MCA_and_incumbent'),audit


def finish(stem,number,caption,methods,audit,files):
    cp=OUT/(stem+'_caption.txt');cp.write_text(caption+'\n')
    mp=OUT/(stem+'_methods.txt');mp.write_text(methods+'\n')
    ap=OUT/(stem+'_audit.json');write(ap,audit)
    files=list(files)+[cp,mp,ap]
    manifest=OUT/(stem+'.manifest.json')
    stamp=datetime.now(timezone.utc).isoformat()
    record=dict(generated_at=stamp,producer='finished',passed=True,new_integrations=0,artifact_root='.',
        renderer=str(Path(__file__).relative_to(ROOT)),renderer_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in files},sha256={p.name:sha(p) for p in files},
        sources=dict(SOURCES),source_sha256=dict(SOURCES),caption=caption,audit=audit,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__)))
    write(manifest,record)
    rel=lambda p:str(p.relative_to(REPORT))
    return dict(figure_label=str(number),image=rel(OUT/(stem+'.png')),manifest=rel(manifest),caption=caption,
        sources=[[rel(mp),'Native-data methods'],[rel(ap),'Source identity and layout audit'],
                 [rel(OUT/(stem+'.pdf')),'Vector figure']])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    protected={name:sha(REPORT/name) for name in ['length_mechanisms.pdf','index.html']}
    source(Path(__file__));source(Path(spatial.__file__))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                        'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    changes,outflux,traces,basal,endpoints,durations,endpoint_audit=load_endpoint()
    source(REPORT/'revision_20261014/visual_palette/semantic_palette.json')
    crossing_path=source(REV/'endpoint_descriptive_crossings.csv')
    endpoint_files=[csvwrite('endpoint_basal_state.csv',basal),csvwrite('endpoint_native_cue_trajectories.csv',endpoints),
                    csvwrite('endpoint_common_time_changes.csv',changes),csvwrite('outward_bound_MCA_common_time_amount.csv',outflux),
                    csvwrite('all_geometry_strength_duration.csv',durations)]
    x,p,snapshots,trace,native_audit,chapter,drawrows,drawsummary=spatial.load()
    SOURCES.update(spatial.SOURCES)
    stem=spatial.STEM
    copied=[]
    for suffix in ['_profiles_and_rates.csv','_face_fluxes.csv','_regional_rates.csv','_global_pool_rates.csv','_snapshot_diagnostics.csv','_parameters.json']:
        path=OUT/(stem+suffix);shutil.copy2(spatial.PRIOR/(stem+suffix),path);copied.append(path)
    copied += [csvwrite(stem+'_drawing_nodes.csv',drawrows),csvwrite(stem+'_drawing_snapshots.csv',drawsummary)]
    audit9=figure9(x,p,snapshots,trace,outflux);audit9['retained_native_audit']=native_audit
    audit9['endpoint_data_audit']=endpoint_audit
    caption9=('Figure 9. MCA redistribution in a 31-µm cell under the right relative-20% cue, '
        'using the shortest tested canonical 1%-scalar-strength-passing pulse (285.280612 s). '
        'A–C: bound MCA before the cue, at 30 s and at 180 s. D–F: cortical velocity and its endpoint, Rho-gradient and membrane-motion contributions. '
        'Positive velocity is rightward; endpoint forcing acts at both Rac-active ends. G–I: membrane springs, viscous cortical dashpots and MCA couplers at the same times. '
        'Gray nodes mark the reference configuration. Membrane displacement and local 10-s cortical slip relative to membrane are magnified ×2.5; these are fixed-reference illustrations, not tracked material particles. '
        'Spring color shows signed tension; coupler width/color depicts bound-MCA drag. Black arrows show recorded flow, purple arrows endpoint commands, and pink signs uniform correction to both MCA pools. '
        'J–L: exact next-step contributions to local bound MCA. Profile advection, flow-gradient dilution/concentration, bound diffusion, net binding and uniform replacement sum to net change. '
        'Rate axes are symmetric logarithmic outside ±1 × 10⁻⁴ native concentration/s, exposing the endpoint cells. '
        'M: cumulative bound MCA advected out through both endpoint faces at common 30-, 60-, 120- and 180-s times, across all twelve eligible lengths. '
        'The amount is native concentration × µm; it is outward throughput and can be replenished, rather than the concentration loss at either endpoint. All common times precede withdrawal. '
        'The spatial profiles and process arrays are unchanged from the audited native L31 trajectory.')
    methods9=('The corrected independent-endpoint model, exact basal checkpoints, parameters, native finite-volume face conventions, 0.5-s staggering and 201-cell grid are retained. '
        'No integration is performed. Original spatial plots and the data-derived mechanics artists are repositioned with every plotted x/y array checked for exact identity. '
        'The mechanics row precedes the bound-MCA process row. Profile advection and flow-gradient terms retain the native donor convention; they are diagnostic source allocations at the coupled state. '
        'Endpoint outflux uses implicit pre-uniform-correction bound donors and commands used in that step; exact duration-weighted native fluxes are summed. '
        'The displayed cortical flow is from the completed preceding mechanical solve, whose source allocation uses the preceding Rho, completed-step bound drag, membrane velocity and endpoint commands. '
        'The local bound-MCA budget is the exact next 0.5-s step: interior transport uses that stored cortical flow, while endpoint commands use the current displayed Rac/tension. '
        'The ordinary mechanical-flow/bound-budget offset is therefore up to one native 0.5-s step and is intentionally retained. '
        'No stored interior velocity or trapezoidal donor approximation substitutes for endpoint commands. The previous total-pool/global split-residual panel remains available in the linked prior audit and companion; it is not used to explain local bound-MCA concentrations. '
        'M uses equal elapsed times, and deliberately excludes own-withdrawal averages. Color meanings and spring/drag scales retain the current semantic palette. '
        'Source-array identity, hash-verified endpoint exports, exact basal states and text bounds are audited. Producer finished; native data are unchanged.')
    files9=[OUT/(stem+'.'+ext) for ext in ['png','pdf','svg']]+copied+[endpoint_files[3]]
    override9=finish(stem,9,caption9,methods9,audit9,files9)
    files10,audit10=figure10(changes,traces,basal,durations);audit10['endpoint_data_audit']=endpoint_audit
    caption10=('Figure 10. Basal endpoint MCA and incumbent-pole dynamics across length. '
        'A: shortest tested cue durations passing the canonical 1% scalar-strength gate at all twelve eligible lengths, for fixed 3-µm, relative-20% and constant-slope 0.24 cues. '
        'All scalar successes are retained; filled/open markers show passing/failing separate native-field tail checks. '
        'B: bound MCA at the old left pole and challenged right end before the cue. '
        'C–D: observed changes in the first/last native-cell bound MCA at common 30-, 60-, 120- and 180-s times; both panels share concentration limits. '
        'These local changes include all bound-MCA processes and are distinct from endpoint outward flux. '
        'E–F: absolute old-end bound MCA and Rac during the relative-20% cue at 23, 31, 35, 37 and 39 µm. Curves end at each selected cue withdrawal, marked by a square. '
        'The dotted Rac line is the smooth protrusion-switch midpoint, 0.899062 native concentration, rather than a polarity-reversal gate or a universal MCA threshold. '
        'All endpoint observations use the native cells nearest the boundaries, x=L/402 and x=L−L/402. '
        'Short-length and longest-length exceptions are shown without fitting a common length law.')
    methods10=('All endpoint fields are read from retained native relative-20% trajectories in debelly_mca_speed_decomposition_20261010. '
        'Native arrays agree exactly with the audited endpoint-draft CSV and the qualifying basal checkpoint. '
        'Common elapsed times 30, 60, 120 and 180 s precede every selected withdrawal (226–809 s); no interpolation is required for these samples. '
        'Concentration change is B(t)−B(0), including every native term. E–F show native cue samples only; no recovery or extrapolation is displayed. '
        'Panel A reuses the current hash-registered three-geometry cue-duration table; marker fill preserves field-tail flags independently of scalar success. '
        'No horizontal MCA concentration threshold is inferred. The Rac reference is the native algebraic protrusion-switch midpoint and does not establish complete incumbent-pole disruption. '
        'No parameter, selection gate, checkpoint, integration or figure 1–7 asset changes are made.')
    override10=finish('Figure_10_endpoint_MCA_and_incumbent',10,caption10,methods10,audit10,files10+endpoint_files)
    override10['sources'].append([str(crossing_path.relative_to(REPORT)),
        'Descriptive incumbent-Rac decline clocks; separate from reversal gates'])
    write(OUT/'figure_overrides.json',{'9':override9,'10':override10})
    write(OUT/'generation_audit.json',dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        protected_report_unchanged={name:sha(REPORT/name)==digest for name,digest in protected.items()},
        new_native_integrations=0,figures=['9','10'],figure9_panel_map={
            'A–C':'bound-MCA profiles','D–F':'cortical velocity sources','G–I':'membrane/cortex mechanics',
            'J–L':'bound-MCA local rates','M':'common-time cumulative outward bound-MCA amount'},
        figure10_panel_map={'A':'three-geometry cue durations','B':'pre-cue endpoint MCA',
            'C':'old-end net bound-MCA change','D':'new-end net bound-MCA change','E':'old-end bound MCA during cue','F':'old-end Rac during cue'}))
    assert all(sha(REPORT/name)==digest for name,digest in protected.items())
    print(json.dumps(dict(passed=True,output=str(OUT),overrides=str(OUT/'figure_overrides.json')),indent=2))


if __name__=='__main__':main()
