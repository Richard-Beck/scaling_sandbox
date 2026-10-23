"""Re-letter Figures 1–2 from frozen outline_v11 data, without new analysis."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / 'polarity1d/manuscript/code'
sys.path[:0] = [str(CODE), str(CODE / 'mechanisms')]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from pypdf import PdfReader, PdfWriter
import figure_builder as foundation
import three_test_figures as three
import revision_round3_figure2_20261022 as previous

DEV = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
OUT = DEV / 'revision_round6_20261023'
BASELINE = json.loads((ROOT / 'manuscript_outline/outline_v11.json').read_text())
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def csvout(folder, name, rows):
    rows = list(rows)
    with (folder / 'supporting/data' / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def readcsv(path):
    SOURCES[relative(path)] = sha(path)
    with path.open() as stream:
        return list(csv.DictReader(stream))


def setup(number):
    folder = OUT / f'figure{number}'
    (folder / 'supporting/data').mkdir(parents=True, exist_ok=True)
    for category in ('metadata', 'vectors', 'viewers'):
        (folder / 'supporting' / category).mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.labelsize': 8, 'axes.titlesize': 8.5, 'xtick.labelsize': 7,
        'ytick.labelsize': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.spines.top': False, 'axes.spines.right': False,
        'text.color': 'black', 'axes.titlecolor': 'black', 'axes.labelcolor': 'black'})
    return folder


def panel(fig, letter, left, top):
    fig.text(left, top, letter, fontsize=12, weight='bold', color='black', va='bottom')


def export(fig, folder, number, section, legend, remap):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            bounds = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            if bounds.x0 < -.003 or bounds.x1 > 1.003 or bounds.y0 < -.003 or bounds.y1 > 1.003:
                clipped.append({'text': artist.get_text(), 'bounds': list(bounds.bounds)})
    write(folder / 'supporting/metadata/display_audit.json', {'one_quantity_per_letter': remap,
        'clipped_text': clipped, 'page_inches': list(fig.get_size_inches())})
    assert not clipped, clipped
    for suffix in ('pdf', 'png', 'svg'):
        fig.savefig(folder / f'figure{number}_round6.{suffix}', dpi=220, facecolor='white')
    page_size = list(fig.get_size_inches())
    plt.close(fig)
    stem = folder / f'figure{number}_round6'
    section.update(status='round6_layout_only', legend=legend,
        assets=[{'path': relative(stem.with_suffix('.png')),
                 'native_pdf_page': relative(stem.with_suffix('.pdf')),
                 'page_size': page_size, 'minimum_printed_width_inches': page_size[0],
                 'page_title': ''}],
        main_pdf=relative(stem.with_suffix('.pdf')),
        supplement_pdf=relative(folder / f'figure{number}_round6_supplement.pdf'),
        reader=relative(folder / 'scientific_readme.md'))
    section['main_source_data'] = [{'path': relative(p), 'sha256': sha(p)}
                                   for p in sorted((folder / 'supporting/data').glob('*.csv'))]
    write(folder / 'legend.json', section)
    (folder / 'legend.md').write_text(section['title'] + '\n\n' + '\n\n'.join(legend) + '\n')
    receipt = {'generated_utc': datetime.now(timezone.utc).isoformat(),
        'baseline_outline_pdf': 'manuscript_outline/outline_v11.pdf',
        'baseline_outline_pdf_sha256': sha(ROOT / 'manuscript_outline/outline_v11.pdf'),
        'new_analysis': False, 'new_simulation': False, 'data_unchanged': True,
        'input_sha256': SOURCES, 'panel_remap': remap,
        'output_sha256': {p.name: sha(p) for p in sorted(folder.glob('*.pdf'))}}
    write(folder / 'supporting/metadata/provenance.json', receipt)
    (folder / 'scientific_readme.md').write_text('''# Unexpected

This layout revision produces no new scientific result. Numerical observations,
eligibility rules and qualifications are those of the starting outline_v11 PDF.

Every distinct plotted quantity has its own letter. Facets repeat the same
quantity. Previous main artwork and complete supplements remain in the supplement.
The plotted coordinates are in `supporting/data/`; source hashes and letter mapping are in
`supporting/metadata/provenance.json`. The source data and model kinetics are unchanged.

Rebuild from the repository root with the existing scaling_sandbox environment:

```sh
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B ''' + relative(folder / 'rebuild.py') + '''
```
''')


def figure1():
    folder = setup(1)
    prep, duration, arrival, prep_path = foundation.preparation_example(ROOT / 'polarity1d/manuscript/data')
    SOURCES[relative(prep_path)] = sha(prep_path)
    correlation_path = ROOT / 'polarity1d/manuscript/data/derived/length_correlations.csv'
    correlations = readcsv(correlation_path)
    shutil.copy2(correlation_path, folder / 'supporting/data/E_length_correlations.csv')
    fig = plt.figure(figsize=(11.8, 8.7))
    grid = fig.add_gridspec(1, 3, left=.075, right=.97, top=.92, bottom=.73, wspace=.4)
    foundation.measured_row(fig, grid, 0, 'A', '', prep, duration, arrival)
    # The historical helper adds one shared A. Replace it with independent A–C.
    for artist in list(fig.texts):
        artist.remove()
    for ax, letter in zip(fig.axes[:3], 'ABC'):
        panel(fig, letter, ax.get_position().x0-.045, .945)
    fig.axes[1].set_xlabel('Position / length')
    fig.axes[1].get_legend().set_bbox_to_anchor((.5, -.25), transform=fig.axes[1].transAxes)
    x, times, active = prep['x'], prep['time'], prep['active']
    q = x / 15
    stop = float(prep['orientation_time'][-1])
    cue_times = np.unique(np.r_[np.linspace(0, duration, 21), duration,
        np.geomspace(max(duration, .001), stop, 201)])
    csvout(folder, 'A_preparation_cue.csv',
        ({'time': float(t), 'position_over_length': float(v),
          'cue_over_peak': float(c)}
         for t in cue_times for v, c in zip(q, foundation.cosine(q, .2, 'left') if t < duration else np.zeros_like(q))))
    indices = sorted(set(int(np.argmin(np.abs(times - t))) for t in [0, duration, duration+30, duration+250, stop]))
    minimum, span = active[indices].min(), np.ptp(active[indices])
    csvout(folder, 'B_active_profiles.csv',
        ({'time': float(times[k]), 'position_over_length': float(v),
          'native_active': float(raw), 'scaled_active': float((raw-minimum)/span)}
         for k in indices for v, raw in zip(q, active[k])))
    csvout(folder, 'C_orientation.csv',
        ({'time': float(t), 'orientation': float(v)} for t, v in zip(prep['orientation_time'], prep['orientation'])))
    panel(fig, 'D', .03, .625)
    curve_rows = []
    position = np.linspace(0, 1, 401)
    for j, title in enumerate(three.LABELS):
        ax = fig.add_axes([.075+.31*j, .42, .26, .175])
        for length, color, style in [(15, '#333333', '-'), (30, '#73A6B8', '--')]:
            signal = foundation.cosine(position, 3/length if j == 0 else .2, 'right') if j < 2 else length/15*position
            ax.plot(position, signal, color=color, ls=style, lw=1.4, label=f'L = {length}')
            curve_rows.extend({'geometry': three.TESTS[j], 'length': length,
                'position_over_length': float(v), 'cue_over_reference_peak': float(c)} for v, c in zip(position, signal))
        ax.set(xlim=(0, 1), ylim=(0, 2.1 if j == 2 else 1.1), xticks=[0, .5, 1], xlabel='Position / length')
        ax.set_title(title, color='black')
        if j == 0:
            ax.set_ylabel('Cue / reference peak')
            ax.legend(frameon=False, fontsize=7)
    csvout(folder, 'D_challenge_geometries.csv', curve_rows)
    panel(fig, 'E', .03, .325)
    lookup = {(r['model'], r['condition'], r['definition']): r for r in correlations}
    matrix = np.array([[float(lookup[m, g, gate]['spearman_rho']) for g in three.TESTS for gate, *_ in foundation.GATES] for m in foundation.MODELS])
    ax = fig.add_axes([.20, .085, .77, .235])
    image = ax.imshow(matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i, model in enumerate(foundation.MODELS):
        for j, condition in enumerate(three.TESTS):
            for k, (gate, *_) in enumerate(foundation.GATES):
                row = lookup[model, condition, gate]
                value = float(row['spearman_rho'])
                ax.text(2*j+k, i, f'{value:+.2f}\nn={row["n_correlation"]}',
                        ha='center', va='center', fontsize=7, color='white' if abs(value) > .65 else 'black')
    ax.set_yticks(range(5), [foundation.LABELS[m] for m in foundation.MODELS])
    ax.set_xticks(range(6), ['Orientation', 'Strength']*3)
    ax.tick_params(length=0)
    for j, title in enumerate(three.LABELS):
        ax.text((2*j+1)/6, 1.08, title, transform=ax.transAxes, ha='center', fontsize=8)
    for border in [1.5, 3.5]:
        ax.axvline(border, color='white', lw=2)
    cax = fig.add_axes([.20, .028, .25, .013])
    bar = fig.colorbar(image, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])
    bar.ax.tick_params(labelsize=7)
    cax.text(1.15, .5, 'Spearman correlation', transform=cax.transAxes, va='center', fontsize=8)
    section = next(s.copy() for s in BASELINE['sections'] if s['id'] == 'figure1')
    old = section['legend']
    legend = [
        'A, Recorded Wave-pinning preparation cue at length 15: a finite left patch, withdrawn at 0.773 native time units. White to black is zero to peak input; coloured vertical lines mark the profile samples in B.',
        'B, Native active profiles from the same preparation, scaled using the common minimum and span of the displayed snapshots. Colours identify recorded times.',
        'C, Orientation from the same preparation. Grey marks basal verification after establishment arrival at about 116; the cue has already ended. The matched time scale in A and C separates input duration from autonomous development.',
        old[1].replace('B,', 'D,', 1), old[2].replace('C,', 'E,', 1), *old[3:]]
    legend[-1] = ('deBelly uses corrected endpoint commands. Spring retains native length and time units. '
                  'Cross-model amplitudes and clocks are not calibrated biochemical doses or cellular times.')
    supplement = PdfWriter()
    for p in (ROOT / 'polarity1d/manuscript/figures/main_polarity.pdf', ROOT / 'polarity1d/manuscript/figures/supplement_width_response.pdf'):
        supplement.append(p)
        SOURCES[relative(p)] = sha(p)
    supplement.write(folder / 'figure1_round6_supplement.pdf')
    export(fig, folder, 1, section, legend,
        {'old A': ['A cue kymograph', 'B active profiles', 'C orientation'], 'old B': 'D challenge shapes', 'old C': 'E correlation heatmap'})


def figure2():
    folder = setup(2)
    oldfolder = DEV / 'revision_round3_20261022/figure2'
    names = ('A_incumbent_organization.csv', 'B_initiation.csv', 'C_replacement_times_all450.csv',
             'D_gap1_vs_pulse2_all180.csv', 'E_L39_challenger_share_traces.csv')
    datasets = []
    for name in names:
        source = oldfolder / 'supporting/data' / name
        rows = readcsv(source)
        shutil.copy2(source, folder / 'supporting/data' / name)
        datasets.append(rows)
    org, initiation, widths, endpoints, traces = datasets
    # Parsing copies the existing coordinates; no summaries are recomputed.
    for rows in datasets:
        for row in rows:
            for key, value in list(row.items()):
                if value not in ('', 'True', 'False'):
                    try:
                        row[key] = float(value)
                    except ValueError:
                        pass
    for row in endpoints:
        row['below_activity_gate'] = row['below_activity_gate'] == 'True'
    for row in widths:
        for key in ('lower_bound', 'time_estimate', 'time_upper'):
            if row[key] == '': row[key] = None
    fig = plt.figure(figsize=(11.8, 16.8))
    def facets(y, height, count=5, left=.075, right=.975, gap=.032):
        width = (right-left-gap*(count-1))/count
        return [fig.add_axes([left+j*(width+gap), y, width, height]) for j in range(count)]
    def models(axes):
        for ax, model in zip(axes, previous.ORDER):
            ax.set_title(previous.NAMES[model], color='black', pad=5)
    fig.legend([Line2D([], [], color=previous.COLORS[m], lw=1.5) for m in previous.ORDER],
        [previous.NAMES[m] for m in previous.ORDER], loc='upper center', bbox_to_anchor=(.53, .998),
        ncol=5, frameon=False, fontsize=8)
    axes = facets(.897, .073, count=3, gap=.065)
    for ax, metric, label, letter in zip(axes, ('excess_amount', 'width', 'peak'), ('Amount', 'Extent', 'Peak'), 'ABC'):
        panel(fig, letter, ax.get_position().x0-.045, .977)
        for model in previous.ORDER:
            rows = sorted([r for r in org if r['model'] == model], key=lambda r:r['size_native'])
            ax.plot([r['size_over_reference'] for r in rows], [r[metric+'_over_reference'] for r in rows], color=previous.COLORS[model], lw=1.2)
        ax.set(xlabel='Length / reference', ylabel=label+' / reference')
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=3, prune='both'))
        previous.clean(ax)
    panel(fig, 'D', .03, .863)
    axes = facets(.785, .067)
    models(axes)
    for ax, model in zip(axes, previous.ORDER):
        rows = [r for r in initiation if r['model'] == model]
        if model == 'debelly':
            ax.plot(previous.LENGTHS, [.5]*6, 'v', ms=4, mfc='white', mec=previous.COLORS[model])
            for length in previous.LENGTHS:
                ax.annotate('', xy=(length,.05), xytext=(length,.46), arrowprops=dict(arrowstyle='->',lw=.7,color=previous.COLORS[model]))
            ax.set(ylim=(0,.6), yticks=[0,.5], ylabel='Native time')
        else:
            values = [r['ratio_estimate'] for r in rows]
            ax.errorbar([r['length'] for r in rows], values,
                yerr=[[v-r['ratio_lower'] for r,v in zip(rows,values)], [r['ratio_upper']-v for r,v in zip(rows,values)]],
                color=previous.COLORS[model], lw=1.1, capsize=2)
            ax.set(yscale='log', ylim=(.75,4.1), yticks=[1,2,4], yticklabels=['1','2','4'])
            ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.axhline(1,color='.7',ls=':',lw=.6)
        ax.set(xlim=(17,41),xticks=[19,29,39],xlabel='Length L')
        previous.clean(ax)
    axes[0].set_ylabel('Initiation / L19')
    panel(fig, 'E', .03, .752)
    fig.legend([Line2D([],[],color='.15',lw=1,ls=s) for s in previous.STYLES],
        [f'w = {w:g}' for w in previous.WIDTHS],loc='upper center',bbox_to_anchor=(.55,.738),ncol=5,frameon=False,fontsize=8)
    for i, area in enumerate(previous.AREAS):
        axes = facets(.654-i*.075, .052)
        if i == 0: models(axes)
        for ax, model in zip(axes, previous.ORDER):
            previous.width_axes(ax,widths,model,area,i)
        axes[0].text(.03,.94,f'Area {area:g}×',transform=axes[0].transAxes,va='top',fontsize=7)
    fig.text(.015,.59,'Replacement time (native units)',rotation=90,va='center',fontsize=8)
    panel(fig, 'F', .03, .467)
    axes=facets(.367,.087)
    models(axes)
    for j,(ax,model) in enumerate(zip(axes,previous.ORDER)):
        others=[r for r in endpoints if r['model']!=model]
        ax.scatter([r['activity_kept'] for r in others],[r['headroom_gain'] for r in others],s=10,c='.84',linewidths=0)
        for row in [r for r in endpoints if r['model']==model]:
            light=.58*(39-row['length'])/20
            color=tuple((1-light)*c+light for c in to_rgb(previous.COLORS[model]))
            ax.scatter(row['activity_kept'],row['headroom_gain'],s=18,
                facecolors='white' if row['below_activity_gate'] else [color],edgecolors=[color],linewidths=.8 if row['below_activity_gate'] else .2)
        ax.set(xlim=(-.055,1.06),ylim=(-.05,1.02),xticks=[0,.5,1],yticks=[0,.5,1],xlabel='Gap-end / pulse 1')
        ax.axhline(0,lw=.5,color='.65');ax.axvline(0,lw=.5,color='.65')
        if j:ax.tick_params(labelleft=False)
        previous.clean(ax)
    axes[0].set_ylabel('Pulse-2 gain / headroom')
    fig.text(.075,.337,'Open circles: below activity gate at x = 0. Size increases from light (L19) to dark (L39).',fontsize=7)
    panel(fig, 'G', .03, .323)
    axes=facets(.239,.071)
    classes = dict(goryachev='lost', wave_pinning='kept', holmes_model3='kept',
                   debelly='lost then commits', spring='lost then commits')
    for ax, model in zip(axes, previous.ORDER):
        ax.set_title(previous.NAMES[model] + '\n' + classes[model], color='black', fontsize=7.5, pad=3)
    for j,(ax,model) in enumerate(zip(axes,previous.ORDER)):
        rows=[r for r in traces if r['model']==model]
        ax.plot([r['pulse_time'] for r in rows],[r['challenger_share'] for r in rows],color=previous.COLORS[model],lw=1.1)
        for pulse in range(4):ax.axvspan(pulse*5,pulse*5+1,color=previous.COLORS[model],alpha=.12,lw=0)
        ax.set(ylim=(-.02,1.04),yticks=[0,.5,1],xlabel='Time / pulse',xlim=(0,max(r['pulse_time'] for r in rows)))
        ax.set_xticks([0,16,32] if model=='debelly' else [0,8,16])
        if j:ax.tick_params(labelleft=False)
        ax.text(.98,.04,'L39',ha='right',transform=ax.transAxes,fontsize=7)
        previous.clean(ax)
    axes[0].set_ylabel('Challenger share')
    panel(fig,'H',.03,.203)
    panel(fig,'I',.03,.098)
    upper=facets(.133,.042,count=3,gap=.065)
    lower=facets(.043,.042,count=3,gap=.065)
    schematic_rows=[]
    for j,(ax,bx,name) in enumerate(zip(upper,lower,('Initiation-limited','Saturating activation','Slow non-activity state'))):
        previous.schematic(ax,j)
        previous.schematic(bx,j,width=True)
        counterpart = ('Goryachev', 'Wave-pinning / Holmes', 'deBelly / Spring')[j]
        ax.set_title(name + '\n' + counterpart,color='black',fontsize=7.5)
        ax.set_xlabel('Time (schematic units)',fontsize=7)
        ax.set_xticks([0,4,8])
        for panel_name,axis in [('H',ax),('I',bx)]:
            for line_index,line in enumerate(axis.lines):
                schematic_rows.extend({'panel':panel_name,'response_class':name,'group':'WGD−' if line_index==0 else 'WGD+',
                    'x':float(x),'y':float(y)} for x,y in zip(line.get_xdata(),line.get_ydata()))
    upper[0].set_ylabel('Biosensor\n(schematic)')
    lower[0].set_ylabel('Replacement time\n(schematic)')
    fig.legend([Line2D([],[],color='.42',lw=1),Line2D([],[],color='.08',lw=1,ls=(0,(3,1.6)))],['WGD−','WGD+'],loc='upper right',bbox_to_anchor=(.98,.216),ncol=2,frameon=False,fontsize=7)
    fig.legend([Line2D([],[],color='.42',lw=1),Line2D([],[],color='.08',lw=1,ls=(0,(3,1.6)))],['WGD−','WGD+'],loc='upper right',bbox_to_anchor=(.98,.112),ncol=2,frameon=False,fontsize=7)
    csvout(folder,'H_I_prospective_schematics.csv',schematic_rows)
    section=next(s.copy() for s in BASELINE['sections'] if s['id']=='figure2')
    section['narrative'] = [
        'At reference cue area 1, replacement is approximately 5–12-fold slower from L19 to L39 in every model wherever both lengths resolve. Enlarged cells resist replacement through distinct initiation, cue-width and pulse-memory signatures.',
        *section['narrative'][1:]]
    old=section['legend']
    legend=[old[0],
        'A, Incumbent excess amount normalized to the model reference, across the same 97 retained states.',
        'B, Incumbent extent normalized to the model reference, for those same states.',
        'C, Incumbent peak normalized to the model reference, for those same states.',
        old[2].replace('B,','D,',1),old[3].replace('C,','E,',1),
        old[4].replace('D,','F,',1),old[5].replace('E,','G,',1),
        'H, Prospective biosensor histories for the same pulse train in WGD− and enlarged WGD+ cells. Initiation-limited resetting (Goryachev), saturating activation with gap retention (Wave-pinning/Holmes), and progress retained after activity loss (deBelly/Spring) give distinct dimensionless schematic histories.',
        'I, Prospective equal-dose width responses for those classes: weak width rescue, strong width rescue, and weaker width rescue. Solid grey: WGD−; dashed black: WGD+. These are qualitative assay predictions.']
    # Repair remaining panel references in retained section prose.
    replacements={' in D':' in F','D codes':'F codes'}
    for key in ('limitations','equations'):
        def remap(value):
            if isinstance(value,str):
                for a,b in replacements.items():value=value.replace(a,b)
            elif isinstance(value,list):value=[remap(v) for v in value]
            elif isinstance(value,dict):value={k:remap(v) for k,v in value.items()}
            return value
        if key in section:section[key]=remap(section[key])
    supplement=PdfWriter()
    for p in (oldfolder/'figure2_round3.pdf',oldfolder/'figure2_round3_supplement.pdf'):
        supplement.append(p);SOURCES[relative(p)]=sha(p)
    supplement.write(folder/'figure2_round6_supplement.pdf')
    export(fig,folder,2,section,legend,{'old A':['A amount','B extent','C peak'],
        'old B':'D initiation','old C':'E replacement','old D':'F retention/gain',
        'old E':'G histories','old F':['H biosensor schematics','I width schematics']})


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--figure',type=int,choices=[1,2],required=True)
    args=parser.parse_args()
    (figure1 if args.figure==1 else figure2)()


if __name__=='__main__':
    main()
