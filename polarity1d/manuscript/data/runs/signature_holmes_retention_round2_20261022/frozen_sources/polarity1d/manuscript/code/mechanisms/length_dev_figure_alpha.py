"""Figure alpha: native MCA transport, protrusion and resting tension.

Read-only analysis of corrected, audited checkpoints and spatial exports. No integrations. Run with the existing scaling_sandbox environment.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from length_dev_basal_front_tension_decomposition import solve
from length_dev_output_layout import dev_output_path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/deBelly_redraft'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
MCA = ROOT / 'data/runs/debelly_mca_threshold_profiles_20261010'
ZERO = ROOT / 'data/runs/debelly_corrected_zero_crossing_20261008'
STEM = 'Figure_alpha_steady_state_mechanism'
COLORS = {23.: '#0072B2', 31.: '#762a83', 39.: '#D55E00'}
ORANGE, BLUE, GREEN, GRAY, INK = '#D55E00', '#0072B2', '#009E73', '#888888', '#252525'
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
    with source(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def output(name):
    return dev_output_path(OUT, name)


def csvout(name, records):
    path = output(name)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def state(path, digest):
    assert sha(source(path)) == digest, path
    with np.load(path) as archive:
        fields = {k[6:]: archive[k].copy() for k in archive.files if k.startswith('state_')}
    assert fields and all(np.isfinite(v).all() for v in fields.values())
    return fields


def main():
    source(Path(__file__))
    source(Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py'))
    source(Path(__file__).with_name('length_dev_output_layout.py'))
    original_report = {name: sha(REPORT / name) for name in ['index.html', 'length_mechanisms.pdf']}
    campaign = read(BASE / 'manifest.json')
    assert not campaign['center_boundary_motion']
    assert read(BASE / 'results/independent_audit.json')['passed']
    meta = read(REPORT / 'dev/resting_state/supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    p = meta['parameters']
    assert not p['center_boundary_motion']
    endpoint_path = REPORT / 'dev/resting_state/supporting/data/basal_endpoint_metrics_vs_length.csv'
    assert sha(endpoint_path) == meta['sha256']['supporting/data/basal_endpoint_metrics_vs_length.csv']
    raw_endpoint = rows(endpoint_path)
    endpoint = [{k: float(v) for k, v in row.items() if k not in ['checkpoint', 'checkpoint_sha256']}
                for row in raw_endpoint]
    native = {}
    for raw in raw_endpoint:
        length = float(raw['length_um'])
        siblings = [c for c in campaign['cases'] if c['model'] == 'debelly' and c['length'] == length]
        assert len(siblings) == 3 and all(c['eligible'] and c['checkpoint_sha256'] == raw['checkpoint_sha256'] for c in siblings)
        native[length] = state(ROOT / raw['checkpoint'], raw['checkpoint_sha256'])
    threshold = endpoint[0]['rac_switch_threshold']

    mca_manifest = read(MCA / 'manifest.json')
    mca_audit = read(MCA / 'results/audit.json')
    assert mca_audit['passed'] and mca_audit['all_profiles_verified']
    for name, digest in mca_manifest['sources'].items():
        assert sha(source(MCA / 'frozen_sources' / name)) == digest
    family = dict(native)
    for inp in mca_audit['inputs']:
        assert inp['convergence']['verified']
        f = state(MCA / inp['file'], inp['sha256'])
        length = float(inp['length'])
        if length in family:
            for name in f:
                np.testing.assert_array_equal(f[name], family[length][name])
        family[length] = f
    near_summary = rows(MCA / 'results/length_summary.csv')
    flat = [float(r['length_um']) for r in near_summary if float(r['bound_mca_range']) < 1e-4]
    polarized = [float(r['length_um']) for r in near_summary if float(r['bound_mca_range']) >= 1e-4]
    loss_bracket = [max(flat), min(polarized)]
    assert loss_bracket == [17.9, 17.925]

    assert read(ZERO / 'independent_audit.json')['passed']
    zero_provenance = read(ZERO / 'provenance.json')
    assert not zero_provenance['center_boundary_motion']
    assert sha(source(ROOT / 'code/core/models/debelly.py')) == zero_provenance['model_sha256']
    linear = read(ZERO / 'linear.json')
    lx = np.array([r['length'] for r in linear])
    lr = np.array([r['rate'] for r in linear])
    crossing = np.flatnonzero(lr[:-1] * lr[1:] < 0)
    assert len(crossing) == 1
    i = crossing[0]
    onset = float(np.interp(0, lr[i:i+2], lx[i:i+2]))
    seeds = rows(ZERO / 'spontaneous_polarizability.csv')
    seed_rate_errors = []
    for r in seeds:
        length, eps = float(r['length']), float(r['seed_amplitude'])
        assert sha(source(ZERO / f'seed_L{length:g}_e{eps:g}.npz')) == r['file_sha256']
        np.testing.assert_allclose(float(r['tangent_rate']), np.interp(length, lx, lr), atol=1e-15)
        seed_rate_errors.append(abs(float(r['tangent_rate']) - float(r['measured_growth_rate'])))
    assert max(seed_rate_errors) < 1e-7

    transport, spatial, audits = [], [], []
    for length, f in sorted(family.items()):
        size = len(f['rac'])
        dx = length / size
        x = (np.arange(size) + .5) * dx
        friction = p['linker_friction'] * np.maximum(f['mca_bound'], 1e-10)
        switch = .5 * (1 + np.tanh((f['rac'][[0, -1]] - threshold) / p['protrusion_switch_scale']))
        commands = p['protrusion_speed'] * switch / (1 + f['tension'][[0, -1]]**2 / p['protrusion_tension_scale'])
        commands[1] *= -1
        zero = np.zeros(size)
        flow_end = solve(friction, p['cortical_viscosity'], zero, dx, *commands)
        flow_rho = solve(friction, p['cortical_viscosity'], p['active_stress'] * p['rho_to_rac_coupling'] * np.gradient(f['rho'], dx), dx)
        velocity_error = float(np.max(abs(flow_end + flow_rho - f['velocity'])))
        # At verified native rest the membrane is stationary. Error includes
        # finite convergence tolerance and the native staggered mechanical solve.
        assert velocity_error < 1e-5, (length, velocity_error)
        face_v = .5 * (f['velocity'][1:] + f['velocity'][:-1])
        adv = face_v * np.where(face_v >= 0, f['mca_bound'][:-1], f['mca_bound'][1:])
        bdiff = -p['mca_bound_diffusion'] * np.diff(f['mca_bound']) / dx
        udiff = -p['mca_unbound_diffusion'] * np.diff(f['mca_unbound']) / dx
        rms = lambda y: float(np.sqrt(np.mean(y**2)))
        stretch = p['membrane_stiffness'] * (f['right_boundary_displacement'][0] - f['left_boundary_displacement'][0]) / length
        transport.append(dict(length_um=length, rms_bound_advective_flux=rms(adv),
            rms_bound_diffusive_flux=rms(bdiff), rms_unbound_diffusive_flux=rms(udiff),
            rms_combined_diffusive_flux=rms(bdiff + udiff), advection_to_combined_diffusion_rms=rms(adv)/rms(bdiff+udiff),
            mean_endpoint_flow_um_per_s=float(flow_end.mean()), mean_rho_flow_um_per_s=float(flow_rho.mean()),
            mean_native_cortical_flow_um_per_s=float(f['velocity'].mean()),
            stationary_flow_reconstruction_max_error_um_per_s=velocity_error,
            affine_endpoint_stretch_tension_pN=float(stretch),
            bound_mca_centroid_um=float(f['mca_bound'] @ x / f['mca_bound'].sum()),
            reflected_bound_mca_centroid_separation_um=float(abs(length - 2 * f['mca_bound'] @ x / f['mca_bound'].sum()))))
        audits.append(dict(length_um=length, checkpoint_finite=True, stationary_flow_max_error_um_per_s=velocity_error))
        spatial.extend(dict(length_um=length, cell_index=j, x_over_L=float(xx/length),
            x_nominal_um=float(xx), endpoint_flow_um_per_s=float(flow_end[j]), rho_flow_um_per_s=float(flow_rho[j]),
            **{k: float(f[k][j]) for k in ['rac', 'rho', 'mca_bound', 'mca_unbound', 'velocity', 'tension']})
            for j, xx in enumerate(x))

    for r in endpoint:
        f = native[r['length_um']]
        np.testing.assert_array_equal([f['rac'][0], f['tension'][0], f['tension'][-1]],
                                     [r['left_rac'], r['left_tension_pN'], r['right_tension_pN']])
        speed = .15 * .5 * (1 + np.tanh((f['rac'][0] - threshold)/.01)) / (1+f['tension'][0]**2/10)
        np.testing.assert_allclose(speed, r['derived_left_command_outward_um_per_s'], atol=1e-15)

    # Cross-check the new summary against independently audited Figure Y-style
    # face-flux exports, which use the end of one additional 0.5-s native step.
    flux_dir = REPORT / 'dev/mca_threshold_profiles/supporting'
    flux_meta = read(flux_dir / 'metadata/mca_threshold_flux_components.manifest.json')
    assert flux_meta['audit_passed']
    flux_path = flux_dir / 'data/mca_threshold_flux_components.csv'
    flux_key = str(flux_path.relative_to(REPORT / 'dev/mca_threshold_profiles'))
    assert sha(flux_path) == flux_meta['sha256'][flux_key]
    prior_flux = rows(flux_path)
    flux_crosscheck_errors = []
    for length in [19.,21.,25.,31.]:
        interior = [r for r in prior_flux if float(r['length_um']) == length and 0 < int(r['face_index']) < 201]
        assert len(interior) == 200
        row = next(r for r in transport if r['length_um'] == length)
        for old,new in [('bound_advection','rms_bound_advective_flux'),
                        ('bound_diffusion','rms_bound_diffusive_flux'),
                        ('unbound_diffusion','rms_unbound_diffusive_flux')]:
            rms = float(np.sqrt(np.mean([float(r[old])**2 for r in interior])))
            flux_crosscheck_errors.append(abs(rms-row[new]))
    assert max(flux_crosscheck_errors) < 1e-10

    flow_path = flux_dir / 'data/mca_threshold_flux_components_cortical_flow.csv'
    flow_key = str(flow_path.relative_to(REPORT / 'dev/mca_threshold_profiles'))
    assert sha(flow_path) == flux_meta['sha256'][flow_key]
    prior_flow = rows(flow_path)
    profile_lengths = [19.,21.,25.,31.]
    flow_superposition_error = max(abs(sum(float(r[k+'_cortical_velocity_um_per_s'])
        for k in ['endpoint','rho','membrane'])-float(r['total_cortical_velocity_um_per_s'])) for r in prior_flow)
    assert flow_superposition_error < 1e-11

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
        'axes.spines.top':False,'axes.spines.right':False,
        'legend.fontsize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(16,11.8))
    a=fig.add_axes([.065,.685,.255,.225])
    d=fig.add_axes([.065,.365,.255,.225])
    bg=fig.add_gridspec(1,4,left=.405,right=.983,bottom=.685,top=.91,wspace=.28)
    cg=fig.add_gridspec(1,4,left=.405,right=.983,bottom=.365,top=.59,wspace=.28)
    bottom=fig.add_gridspec(1,4,left=.065,right=.983,bottom=.075,top=.265,wspace=.35)
    e,f,g,h=[fig.add_subplot(bottom[0,i]) for i in range(4)]
    def panel(ax,letter):
        ax.text(0,1.07,letter,transform=ax.transAxes,fontsize=13,weight='bold')
    def style(ax):
        ax.grid(alpha=.12,lw=.6)
    for ax,letter in [(a,'A'),(d,'D'),(e,'E'),(f,'F'),(g,'G'),(h,'H')]:
        panel(ax,letter);style(ax)
    for length,color,ls in [(17.,'#aaaaaa','--'),(18.,INK,':'),(19.,GREEN,'-.'),
                            *[(L,c,'-') for L,c in COLORS.items()]]:
        field=family[length];q=(np.arange(len(field['rac']))+.5)/len(field['rac'])
        a.plot(q,field['mca_bound'],color=color,ls=ls,lw=1.9,label=f'{length:g} µm')
    a.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Bound MCA (native density)',ylim=(.043,.115))
    a.legend(ncol=3,frameon=False,loc='upper left',columnspacing=1.1,handlelength=1.8)
    flux_specs=[('bound_advection','Bound advection',ORANGE,'-'),
        ('bound_diffusion','Bound diffusion',GREEN,'-'),
        ('unbound_advection','Unbound advection','#CC79A7',':'),
        ('unbound_diffusion','Unbound diffusion',BLUE,'-'),
        ('total_flux','Total flux',INK,'--')]
    flow_specs=[('endpoint','Endpoint forcing',ORANGE,'-'),
        ('rho','Rho stress gradient',BLUE,'-'),
        ('membrane','Membrane motion',GRAY,'-'),
        ('total','Total flow',INK,'--')]
    for index,length in enumerate(profile_lengths):
        sub=bg[0,index].subgridspec(2,1,height_ratios=[2.2,1],hspace=.12)
        main=fig.add_subplot(sub[0]);detail=fig.add_subplot(sub[1],sharex=main)
        data=[r for r in prior_flux if float(r['length_um'])==length]
        assert len(data)==202
        x=np.array([float(r['x_nominal_um']) for r in data])
        for key,label,color,ls in flux_specs:
            main.plot(x,[1000*float(r[key]) for r in data],color=color,ls=ls,lw=1.6)
        detail.plot(x,[1e6*float(r['bound_diffusion']) for r in data],color=GREEN,lw=1.5)
        main.set(xlim=(0,length),ylim=(-2.15,.38))
        detail.set(xlim=(0,length),ylim=(-3.8,.25),xlabel='x (µm)',xticks=[0,round(length/2),length])
        main.tick_params(labelbottom=False)
        main.set_title(f'{length:g} µm',fontsize=10,pad=8)
        for ax in [main,detail]:
            ax.axhline(0,color='#aaaaaa',lw=.65,zorder=0);style(ax)
        if index==0:
            panel(main,'B');main.set_ylabel('MCA flux\n(10⁻³ conc. µm/s)')
            detail.set_ylabel('Bound diff.\n(10⁻⁶)')
        else:
            main.tick_params(labelleft=False);detail.tick_params(labelleft=False)
        cortex=fig.add_subplot(cg[0,index])
        data=[r for r in prior_flow if float(r['length_um'])==length]
        assert len(data)==201
        x=np.array([float(r['x_nominal_um']) for r in data])
        for key,label,color,ls in flow_specs:
            cortex.plot(x,[float(r[key+'_cortical_velocity_um_per_s']) for r in data],color=color,ls=ls,lw=1.7)
        cortex.set(xlim=(0,length),ylim=(-.045,.016),xlabel='x (µm)',xticks=[0,round(length/2),length])
        cortex.set_title(f'{length:g} µm',fontsize=10,pad=8)
        cortex.axhline(0,color='#aaaaaa',lw=.65,zorder=0);style(cortex)
        if index==0:
            panel(cortex,'C');cortex.set_ylabel('Cortical velocity (µm/s)')
        else:cortex.tick_params(labelleft=False)
    def handles(specs):
        return [Line2D([],[],color=color,ls=ls,lw=1.8,label=label) for _,label,color,ls in specs]
    fig.legend(handles=handles(flux_specs),loc='upper center',bbox_to_anchor=(.694,.978),
        ncol=3,frameon=False,columnspacing=1.7,fontsize=9)
    fig.legend(handles=handles(flow_specs),loc='upper center',bbox_to_anchor=(.694,.656),
        ncol=4,frameon=False,columnspacing=1.4,fontsize=9)
    for length,color in COLORS.items():
        field=native[length];q=(np.arange(len(field['rac']))+.5)/len(field['rac'])
        d.plot(q,field['rac'],color=color,lw=1.9)
        d.plot(q,field['rho'],color=color,lw=1.9,ls='--')
        h.plot(q,field['tension'],color=color,lw=1.9,label=f'{length:g} µm')
    d.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Rac / Rho (native conc.)',ylim=(.54,.94))
    d.legend(handles=[Line2D([],[],color=INK,lw=1.9,label='Rac'),Line2D([],[],color=INK,ls='--',lw=1.9,label='Rho')],
        frameon=False,loc='upper right',ncol=2,fontsize=9)
    h.set(xlim=(0,1),xticks=[0,.5,1],xlabel='x/L',ylabel='Membrane tension (pN)',ylim=(-.45,10))
    h.axhline(0,color=GRAY,lw=.65);h.legend(frameon=False,loc='upper left')
    lengths=np.array([r['length_um'] for r in endpoint])
    def line(ax,key,color,ls='-',label=None):
        ax.plot(lengths,[r[key] for r in endpoint],color=color,ls=ls,lw=1.8,marker='o',ms=3.5,label=label)
    for ax in [e,f,g]:
        ax.set(xlim=(17,40),xticks=[18,23,31,39],xlabel='Length (µm)')
    line(e,'left_rac',BLUE)
    e.axhline(threshold,color='#9c2525',ls=':',lw=1.2,label='Switch midpoint')
    e.set(ylabel='Front Rac (native conc.)',ylim=(.844,.908))
    e.legend(frameon=False,loc='lower right',fontsize=8.5)
    line(f,'derived_left_command_outward_um_per_s',BLUE,label='Front')
    line(f,'derived_right_command_outward_um_per_s',GRAY,'--','Rear')
    f.set(ylabel='Outward command (µm/s)',ylim=(-.004,.105))
    f.legend(frameon=False,loc='upper left')
    line(g,'left_tension_pN',BLUE,label='Front')
    line(g,'right_tension_pN',GRAY,'--','Rear')
    g.set(ylabel='Endpoint tension (pN)',ylim=(-.45,10))
    g.axhline(0,color=GRAY,lw=.65);g.legend(frameon=False,loc='upper left')
    fig.text(.065,.975,'Figure α',fontsize=14,weight='bold',va='top')
    by_length={r['length_um']:r for r in transport}
    distance_ratio=by_length[39.]['reflected_bound_mca_centroid_separation_um']/by_length[23.]['reflected_bound_mca_centroid_separation_um']
    command_ratio=endpoint[-1]['derived_left_command_outward_um_per_s']/next(r for r in endpoint if r['length_um']==23)['derived_left_command_outward_um_per_s']
    caption=f'''Figure alpha. Native basal MCA transport, protrusion and tension.

A. Bound-MCA profiles at 17, 18, 19, 23, 31 and 39 µm. Only position is divided by nominal length; densities retain native units. The 17-µm profile is a nearly homogeneous continuation of the polarized 18-µm state. Blue, purple and orange denote 23, 31 and 39 µm throughout A, D and H.
B. Full spatial MCA face-flux components at 19, 21, 25 and 31 µm. Bound advection uses cortical velocity; unbound advection uses membrane velocity. Bound and unbound diffusion use their respective native diffusivities. The lower axes magnify bound diffusion by 1000 relative to the upper axes. Vertical scales are shared across lengths for each quantity. Positive flux is rightward. These are instantaneous spatial fluxes; local exchange and uniform mass restoration are outside their sum. Figure Y explains the boundary sweep-out and replenishment.
C. Full spatial cortical-velocity components at the same four lengths: endpoint forcing, Rho stress gradient, membrane-motion drag and native total. Their pointwise sum is independently verified. Vertical scales are shared. Positive velocity is rightward. These are linear source superpositions at the recorded native state, not independently evolved pathway ablations.
D. Rac (solid) and Rho (dashed) profiles at 23, 31 and 39 µm on a common native concentration axis.
E. Rac at the old left front across all twelve qualified original lengths. The dotted line marks the protrusion-switch midpoint, {threshold:.10f}. The first sampled length above it is 37 µm; the switch remains smooth below and above the midpoint.
F. Actual native outward endpoint commands, evaluated at each checkpoint. Front command is blue; rear command is gray dashed. The native endpoint law and model parameters are unchanged.
G. Native membrane tension in the first and last cells, representing the front and rear.
H. Full membrane-tension profiles at 23, 31 and 39 µm.

A and D–H use exact corrected basal checkpoints. B–C reuse the independently audited spatial exports sampled after one additional 0.5-s native no-cue step, which resolves membrane motion; no new integration was performed for this revision. The same shared corrected preparations underlie all three reversal assays. Original convergence and holdout checks pass. Positions and lengths follow the native nominal-coordinate convention. The sampled continuation loss at 17.9–17.925 µm and homogeneous small-seed onset near {onset:.2f} µm remain supporting context; neither is inferred from a spatial-flux magnitude crossing.
'''
    files = []
    for ext in ['png','pdf','svg']:
        path = output(STEM+'.'+ext); fig.savefig(path,dpi=220,facecolor='white');files.append(path)
    # Verify all plot annotations are inside the saved page, including inset text.
    fig.canvas.draw(); renderer=fig.canvas.get_renderer(); bounds=fig.bbox
    outside=[]
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box=artist.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > bounds.width+1 or box.y1 > bounds.height+1:
                outside.append(artist.get_text())
    assert not outside, outside
    plt.close(fig)
    files.append(csvout(STEM+'_transport_flow.csv',transport))
    files.append(csvout(STEM+'_spatial_profiles.csv',spatial))
    files.append(csvout(STEM+'_endpoint_metrics.csv',raw_endpoint))
    files.append(csvout(STEM+'_small_seed_linear_rates.csv',linear))
    files.append(csvout(STEM+'_spatial_flux_components.csv',prior_flux))
    files.append(csvout(STEM+'_spatial_flow_components.csv',prior_flow))
    path=output(STEM+'_caption.txt');path.write_text(caption);files.append(path)
    audit = dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),
        original_checkpoints_hash_verified=True,shared_three_geometry_starts_verified=True,
        all_profiles_verified=True,spatial_profile_lengths_um=profile_lengths,
        independent_spatial_flow_superposition_max_error_um_per_s=flow_superposition_error,
        spatial_flux_and_flow_exports_hash_verified=True,
        plot_text_within_page=True,new_integrations=0,small_seed_growth_max_abs_rate_error_per_s=max(seed_rate_errors),
        small_seed_onset_um=onset,sampled_continuation_loss_bracket_um=loss_bracket,
        maximum_stationary_flow_reconstruction_error_um_per_s=max(a['stationary_flow_max_error_um_per_s'] for a in audits),
        independent_prior_face_flux_rms_max_error=max(flux_crosscheck_errors),
        reflected_centroid_distance_ratio_39_over_23=distance_ratio,front_command_ratio_39_over_23=command_ratio,
        original_report_unchanged=all(sha(REPORT/name)==digest for name,digest in original_report.items()),checks=audits)
    assert audit['original_report_unchanged']
    path=output(STEM+'_audit.json');path.write_text(json.dumps(audit,indent=2)+'\n');files.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',artifact_root='../..',
        caption=caption,source_sha256=SOURCES,sha256={str(f.relative_to(OUT)):sha(f) for f in files},
        renderer_sha256=sha(__file__),revision=2,new_integrations=0,audit=audit,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_dev_figure_alpha.py')
    output(STEM+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(figure=str(output(STEM+'.png')),audit=audit),indent=2))


if __name__ == '__main__':
    main()
