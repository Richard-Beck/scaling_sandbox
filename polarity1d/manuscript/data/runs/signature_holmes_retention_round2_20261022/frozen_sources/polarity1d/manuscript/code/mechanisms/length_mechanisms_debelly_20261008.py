"""Focused report views from independently audited corrected deBelly arrays.

No integration or model changes. Outputs preserve the original mechanism figures.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / 'data/runs'
ZERO = RUNS / 'debelly_corrected_zero_crossing_20261008'
FLOW = RUNS / 'debelly_corrected_advection_20261008'
OUT = RUNS / 'length_mechanisms_20261008/debelly'
FIG = ROOT / 'figures/length_mechanisms_20261008'
LENGTHS = np.arange(23., 40., 2.)
COLORS = ['#0072B2', '#E69F00', '#009E73']
USED = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    USED[str(path.relative_to(ROOT))] = sha(path)
    return json.loads(path.read_text())


def archive(path):
    USED[str(path.relative_to(ROOT))] = sha(path)
    return np.load(path)


def polarity(field, x, length):
    return np.sum(field * (2 * x / length - 1), axis=-1) / np.sum(field, axis=-1)


def save(fig, name, caption):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(FIG / f'{name}.{ext}', dpi=200)
    (FIG / f'{name}_caption.txt').write_text(caption + '\n')
    plt.close(fig)


def style(ax):
    ax.grid(alpha=.15)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    for run in [ZERO, FLOW]:
        assert read(run / 'independent_audit.json')['passed']
    source_budget = read(RUNS / 'debelly_corrected_advection_revision_20261008/figure_audit.json')
    assert source_budget['max_polarity_budget_residual'] < 1e-11
    for run in [ZERO, FLOW]:
        provenance = read(run / 'provenance.json')
        if 'sources' in provenance:
            for path, digest in provenance['sources'].items():
                assert sha(path) == digest
        else:
            assert provenance['center_boundary_motion'] is False
            assert sha(ROOT / 'code/core/models/debelly.py') == provenance['model_sha256']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    summaries = read(ZERO / 'withdrawal_summary.json')
    withdrawal_rows = []
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 8.7))
    fig.subplots_adjust(left=.075, right=.975, top=.88, bottom=.14,
                        wspace=.33, hspace=.52)
    fig.suptitle('deBelly: MCA attachment bias provides a useful withdrawal landmark',
                 fontsize=17, weight='bold', y=.98)
    fig.text(.075, .925, 'After the cue ends, residual Rac/Rho activity keeps shifting MCA; strong polarity then rebuilds in the amplifying regime.', fontsize=10.5, color='#43525B')
    for col, length in enumerate([18., 23., 39.]):
        ax = axes[0, col]
        for offset, color in zip([-100., -50., 0.], COLORS):
            path = ZERO / f'withdrawal_L{length:g}_offset{offset:g}.npz'
            with archive(path) as z:
                assert 'preparation_ineligible' not in z.files
                t = z['recovery_time']
                pm = polarity(z['recovery_field_mca_bound'], z['x'], length)
                pr = polarity(z['recovery_field_rac'], z['x'], length)
                label = f'{abs(offset):g} s before zero' if offset else 'At MCA zero'
                ax.plot(t, pm, color=color, lw=1.9, label=label)
                if offset == 0:
                    axes[1, 0].plot(t, pr, color=['#8D9499', '#0072B2', '#009E73'][col],
                                    lw=2, label=f'L = {length:g} µm')
                row = next(a for a in summaries if a['length'] == length and a['offset'] == offset)
                assert sha(path) == row['archive_sha256']
                assert abs(row['final_rac'] - pr[-1]) < 1e-12
                assert abs(row['final_mca'] - pm[-1]) < 1e-12
                withdrawal_rows.append(dict(length=length, offset=offset,
                                            final_mca=float(pm[-1]), final_rac=float(pr[-1])))
        ax.axhline(0, color='#666', lw=.8)
        ax.set_xscale('symlog', linthresh=100)
        ax.set_xlim(0, 10000)
        ax.set_ylim(-.15, .055)
        ax.set_xticks([0, 100, 1000, 10000], ['0', '100', '1,000', '10,000'])
        ax.set_xlabel('Time after cue withdrawal (s)')
        ax.set_title(f'{chr(65+col)}   L = {length:g} µm', loc='left', weight='bold')
        outcome = 'Weak polarity at 10,000 s' if length == 18 else 'Strong reversed polarity retained'
        ax.text(.04, .06, outcome, transform=ax.transAxes, fontsize=9,
                color='#43525B')
        style(ax)
    axes[0, 0].set_ylabel('MCA attachment bias\npositive = old pole; negative = new pole')
    axes[0, 1].legend(frameon=False, fontsize=8, loc='upper right')
    ax = axes[1, 0]
    ax.set_xscale('symlog', linthresh=100)
    ax.set_xlim(0, 10000)
    ax.set_xticks([0, 100, 1000, 10000], ['0', '100', '1,000', '10,000'])
    ax.axhline(.05, color='#999', ls=':', lw=1)
    ax.set_ylim(0, .55)
    ax.set_ylabel('Rightward Rac polarity')
    ax.set_xlabel('Time after withdrawal at MCA zero (s)')
    ax.set_title('D   Rac settles after cue removal', loc='left', weight='bold')
    ax.legend(frameon=False, fontsize=8)
    linear = read(ZERO / 'linear.json')
    ll = np.array([a['length'] for a in linear])
    rates = np.array([a['rate'] for a in linear])
    i = np.flatnonzero(rates[:-1]*rates[1:] < 0)[0]
    critical = float(np.interp(0, rates[i:i+2], ll[i:i+2]))
    ax = axes[1, 1]
    ax.plot(ll, rates*1000, color='#65418A', lw=2)
    ax.axhline(0, color='#777', lw=.8)
    ax.axvline(critical, color='#777', ls=':', lw=1)
    ax.set_title('E   Short cells amplify polarity less readily', loc='left', weight='bold', fontsize=10)
    ax.set_xlabel('Cell length (µm)')
    ax.set_ylabel('Growth of small perturbations (10⁻³ s⁻¹)')
    ax.text(.98, .05, f'Spontaneous onset ≈ {critical:.2f} µm', transform=ax.transAxes, ha='right', fontsize=9)
    seedrows = []
    with (ZERO / 'spontaneous_polarizability.csv').open() as f:
        seeds = list(csv.DictReader(f))
    USED[str((ZERO/'spontaneous_polarizability.csv').relative_to(ROOT))] = sha(ZERO/'spontaneous_polarizability.csv')
    ax = axes[1, 2]
    for eps, ls, marker in [(1e-5, '-', 'o'), (1e-6, '--', 's')]:
        rows = [a for a in seeds if float(a['seed_amplitude']) == eps]
        ax.plot([float(a['length']) for a in rows],
                [float(a['final_abs_rac_polarity']) for a in rows],
                ls=ls, marker=marker, ms=4, color='#0072B2',
                label=f'MCA seed {eps:g}')
        seedrows += rows
        for row in rows:
            path = ZERO / f"seed_L{float(row['length']):g}_e{eps:g}.npz"
            assert sha(path) == row['file_sha256']
            USED[str(path.relative_to(ROOT))] = row['file_sha256']
    ax.axvspan(critical, 23, color='#E69F00', alpha=.13)
    ax.axhline(.05, color='#999', ls=':', lw=1)
    ax.set_title('F   Amplification is slow near onset', loc='left', weight='bold')
    ax.set_xlabel('Cell length (µm)')
    ax.set_ylabel('|Rac polarity| after 40,000 s')
    ax.legend(frameon=False, fontsize=8)
    for ax in axes[1]:
        style(ax)
    caption = ('Independent local endpoint commands; exact corrected basal starts. Continuous gradient slope 0.24. '
               'A–D: all three lengths have qualifying polarized initial states. MCA zero works well as a conservative '
               'withdrawal landmark at L=23 and 39; at L=18 the new bias remains too weak at 10,000 s. '
               'E–F: no cue, homogeneous basal rest with small zero-mean MCA perturbations. Short-cell branch '
               'establishment and spontaneous amplification are different questions. Time axes compress after 100 s.')
    fig.text(.075, .06, 'All starts corrected and polarized. L=18 has weaker post-cue amplification.', fontsize=9, color='#43525B')
    save(fig, 'debelly_memory_landmark', caption)

    native = {}
    measurement = []
    stage_rows = []
    cases = read(FLOW/'cases.json')
    for length in LENGTHS:
        summary = read(FLOW/f'L{length:g}_summary.json')
        file = FLOW/f'L{length:g}_native.npz'
        assert sha(file) == summary['native_sha256']
        case = next(a for a in cases if a['length'] == length)
        checkpoint = ROOT / case['checkpoint']
        assert sha(checkpoint) == case['checkpoint_sha256'] == summary['checkpoint_sha256']
        with archive(file) as z:
            pm = polarity(z['field_mca_bound'], z['x'], length)
            time = z['time'].copy()
            j = np.flatnonzero(pm <= 0)[0]
            assert time[j] == summary['crossing_time']
            native[length] = (time, pm)
            shift = float(length / 2 * (pm[0] - pm[j]))
            measurement.append(dict(length=float(length), time_to_zero=float(time[j]),
                                    initial_mca_polarity=float(pm[0]), nominal_centroid_shift=shift,
                                    mean_net_centroid_speed=shift/float(time[j]),
                                    checkpoint_sha256=case['checkpoint_sha256']))
        stage_rows += summary['rows']
    lookup = {(a['length'], a['stage']): a for a in stage_rows}
    times = np.array([a['time_to_zero'] for a in measurement])
    coef = np.polyfit(LENGTHS, times, 1)
    # Independent closed-form OLS agrees with NumPy fit.
    slope = np.sum((LENGTHS-LENGTHS.mean())*(times-times.mean()))/np.sum((LENGTHS-LENGTHS.mean())**2)
    intercept = times.mean()-slope*LENGTHS.mean()
    np.testing.assert_allclose(coef, [slope, intercept], atol=1e-10)
    r2 = float(1-np.sum((times-np.polyval(coef,LENGTHS))**2)/np.sum((times-times.mean())**2))
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 8.7))
    fig.subplots_adjust(left=.075, right=.975, top=.88, bottom=.16,
                        wspace=.34, hspace=.53)
    fig.suptitle('deBelly: a longer MCA shift and slower late progress delay memory reversal', fontsize=16, weight='bold', y=.98)
    fig.text(.075, .925, 'Length changes cortical flow; the new bias must travel farther, while opposing endpoint flow slows its late progress.', fontsize=10.5, color='#43525B')
    ax = axes[0, 0]
    ax.plot(LENGTHS, np.polyval(coef, LENGTHS), color='#AAA', lw=1.5, ls='--', label='Descriptive linear fit')
    ax.plot(LENGTHS, times, color='#65418A', marker='o', lw=1.7, label='Audited MCA-zero times')
    ax.set_ylabel('Time to MCA zero (s)')
    ax.set_title('A   Approximately linear timing', loc='left', weight='bold')
    ax.text(.04, .92, f'R² = {r2:.3f} · slope {slope:.1f} s/µm', transform=ax.transAxes, va='top', fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax = axes[0, 1]
    ax.plot(LENGTHS, [a['nominal_centroid_shift'] for a in measurement], color='#0072B2', marker='o', lw=1.8)
    ax.set_ylabel('Nominal MCA centroid shift to zero (µm)')
    ax.set_title('B   The required shift grows', loc='left', weight='bold')
    ax = axes[0, 2]
    ax.plot(LENGTHS, [a['mean_net_centroid_speed']*1000 for a in measurement], color='#009E73', marker='o', lw=1.8)
    ax.set_ylabel('Mean net centroid speed (10⁻³ µm/s)')
    ax.set_ylim(0, .006*1000)
    ax.set_title('C   Mean speed does not keep pace', loc='left', weight='bold')
    for ax in axes[0]:
        ax.set_xlabel('Cell length (µm)')
        style(ax)
    # Static mechanical control: retain local Rho forcing and extend its low tail.
    sys.path.insert(0, str(ROOT/'code/core'))
    import native as n
    p = n.model_for('debelly').parameters
    b0 = .1*5/6.4
    friction = p.linker_friction*b0
    profilefile = FLOW/'L23_profiles.npz'
    assert sha(profilefile) == read(FLOW/'L23_summary.json')['profiles_sha256']
    with archive(profilefile) as z:
        refx = z['zero_x_fraction']*23
        refrho = z['zero_rho_concentration']
    def controlled(length, cells=2001):
        dx = length/cells
        x = (np.arange(cells)+.5)*dx
        rho = np.interp(x,refx,refrho,left=refrho[0],right=refrho[-1])
        force = p.active_stress*n.ds._gradient(rho,dx)
        v = n.ds._dirichlet_elliptic(np.full(cells,friction),p.cortical_viscosity,force,dx,0.,0.)
        return x, rho, -v
    ax = axes[1,0]
    checks = []
    for length, color in zip([23.,39.,60.,100.], ['#65418A','#009E73','#E69F00','#B55D60']):
        x,rho,v = controlled(length)
        ax.plot(np.r_[0,x,length],np.r_[0,v,0],color=color,lw=1.7,label=f'L={length:g} µm')
        ax.scatter(length,0,color=color,s=20,zorder=3)
        xf,rf,vf = controlled(length,4001)
        error = float(abs(np.interp(5,x,v)/np.interp(5,xf,vf)-1))
        assert error < .001
        checks.append(dict(length=length,relative_mesh_change=error))
    ax.axvspan(0,23,color='#B4BEC5',alpha=.12)
    ax.set_xlim(0,100)
    ax.set_ylim(0,.018)
    ax.set_xlabel('Nominal distance from left end (µm)')
    ax.set_ylabel('Rho-only speed toward left (µm/s)')
    ax.set_title('D   Same Rho profile; longer return path',loc='left',weight='bold',fontsize=10)
    ax.legend(frameon=False,fontsize=8,ncol=2,loc='upper left')
    ax.text(.03,.05,'Rho profile held fixed\nVelocity = 0 at both ends',transform=ax.transAxes,fontsize=8)
    ax.axvline(5,color='#777',ls=':',lw=.8)
    style(ax)
    ax = axes[1,1]
    control = RUNS/'debelly_corrected_advection_revision_20261008/fixed_reference_rho_extension.csv'
    USED[str(control.relative_to(ROOT))] = sha(control)
    with control.open() as f:
        controls = list(csv.DictReader(f))
    for keepfriction, color, label in [('True','#0072B2','With MCA friction'),('False','#999','Without friction')]:
        rows = [a for a in controls if a['friction']==keepfriction]
        ax.plot([float(a['length']) for a in rows],[float(a['speed_at_x5']) for a in rows],
                color=color,lw=1.8,ls='-' if keepfriction=='True' else '--',label=label)
    ax.set_xlim(23,150)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Cell length (µm)')
    ax.set_ylabel('Fixed-profile speed at x=5 µm (µm/s)')
    ax.set_title('E   Geometry strengthens local Rho flow',loc='left',weight='bold',fontsize=10)
    ax.legend(frameon=False,fontsize=8,loc='lower right')
    ax.text(.04,.94,f'MCA screening length ≈ {np.sqrt(p.cortical_viscosity/friction):.1f} µm',transform=ax.transAxes,va='top',fontsize=8)
    style(ax)
    ax = axes[1,2]
    for key, label, color in [('boundary','Endpoint advection','#D55E00'), ('rho','Rho advection','#0072B2')]:
        ax.plot(LENGTHS,[lookup[length,'zero']['rate_'+key]*1e4 for length in LENGTHS], color=color, lw=1.8, label=label)
    other = [sum(lookup[length,'zero']['rate_'+key] for key in ['membrane','diffusion','exchange','replenishment']) for length in LENGTHS]
    ax.plot(LENGTHS,np.array(other)*1e4,color='#92999E',lw=1.8,label='Other terms')
    ax.plot(LENGTHS,[lookup[length,'zero']['reorientation_rate']*1e4 for length in LENGTHS], color='#202020', ls='--', lw=2,label='Net rate')
    assert max(abs(a['budget_residual']) for a in stage_rows) < 1e-11
    ax.axhline(0,color='#777',lw=.8)
    ax.set_title('F   Endpoint flow opposes late progress',loc='left',weight='bold',fontsize=10)
    ax.set_ylabel('MCA reorientation rate (10⁻⁴ s⁻¹)')
    ax.set_xlabel('Cell length (µm)')
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    style(ax)
    caption = ('Corrected independent-endpoint trajectories; continuous gradient slope 0.24, L=23–39 µm. '
               'A: first recorded nonpositive bound-MCA polarity; step 0.5 s. B–C: centroid shift and time-averaged '
               'net centroid speed computed from the same trajectories. The centroid is a nominal-coordinate MCA '
               'distribution measure, not a tracked material parcel. D–E: a static mechanical control holds the L=23 '
               'local Rho profile fixed, extending only its terminal low-Rho region. Uniform basal MCA friction and '
               'zero boundary/membrane forcing isolate domain geometry. A farther endpoint allows a more gradual '
               'velocity return; friction limits this gain. F: exact instantaneous MCA-rate budget at zero in native '
               'trajectories. Rho helps reorientation; endpoint advection opposes it. Exchange and other '
               'terms are retained in the net rate. These zero times are separate from minimum reversal-pulse durations.')
    fig.text(.075,.07,'MCA-zero time = centroid shift / mean net centroid speed. Nominal coordinates; distributions also change by exchange and diffusion.',fontsize=9,color='#43525B')
    save(fig,'debelly_flow_length',caption)
    with (OUT/'length_measurements.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=measurement[0].keys());writer.writeheader();writer.writerows(measurement)
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',passed=True,
               corrected_independent_endpoint_only=True,no_new_simulations=True,
               source_sha256=USED,renderer_sha256=sha(Path(__file__)),
               descriptive_fit=dict(slope_s_per_um=float(slope),intercept_s=float(intercept),r_squared=r2,
                                    closed_form_verified=True,length_min=23,length_max=39),
               centroid_measurement_definition='L/2*(P_initial-P_first_nonpositive); mean net speed=shift/time',
               native_budget_max_residual=max(abs(a['budget_residual']) for a in stage_rows),
               controlled_profile_mesh_checks=checks,
               withdrawal_final_moments=withdrawal_rows,
               spontaneous_linear_onset_um=critical,
               figure_sha256={p.name:sha(p) for name in ['debelly_memory_landmark','debelly_flow_length']
                              for p in FIG.glob(name+'.*')})
    (OUT/'figure_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','withdrawal_final_moments','figure_sha256']},indent=2))


if __name__ == '__main__':
    main()
