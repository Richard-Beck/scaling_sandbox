"""A compact, hash-checked view of the audited sink plateau; no integration."""
from pathlib import Path
import csv
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'data/runs/inactive_sink_loss_20261008'
OUT = ROOT / 'figures/length_mechanisms_20261008'
DATA = ROOT / 'data/runs/length_mechanisms_20261008/sink'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source = RUN / 'results/peak_and_amount_measurements.csv'
    provenance = json.loads((RUN / 'results/peak_render_provenance.json').read_text())
    assert sha(source) == provenance['measurement_sha256'] and provenance['complete']
    for filename in ['audit.json', 'numerical_accuracy.json', 'wave_peak_dip_accuracy.json']:
        assert json.loads((RUN / 'results' / filename).read_text())['passed']
    all_rows = list(csv.DictReader(source.open()))
    rows = [r for r in all_rows if float(r['strength']) == 1000 and float(r['length']) >= 9]
    assert len(rows) == 72
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    with (DATA / 'plateau_measurements.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(10.6, 7.5))
    grid = fig.add_gridspec(3, 2, height_ratios=[.52, 1, 1], left=.12,
                           right=.96, top=.94, bottom=.16, hspace=.52, wspace=.33)
    geometry = fig.add_subplot(grid[0, :])
    geometry.set(xlim=(0, 40), ylim=(-.5, 2.3)); geometry.axis('off')
    for y, sink_x, text in [(1.35, 32, 'Opposite end: pole–sink distance grows with length'),
                             (.1, 12, 'Fixed gap: 3 beyond the initial pole edge')]:
        geometry.plot([0, 33], [y, y], color='#bac4ca', lw=8, solid_capstyle='round')
        geometry.add_patch(Rectangle((0, y-.18), 8, .36, color='#32729b', zorder=3))
        geometry.add_patch(Rectangle((sink_x, y-.25), 1, .50, color='#d46c31', zorder=4))
        geometry.annotate('', xy=(sink_x, y+.43), xytext=(8, y+.43),
                          arrowprops={'arrowstyle':'<->', 'color':'#6f7b84', 'lw':1})
        geometry.text(35, y, 'inactive\nsink', color='#b65423', va='center', fontsize=9)
        geometry.text(0, y+.64, text, fontsize=10)
    geometry.text(3, -.48, 'incumbent', ha='center', fontsize=9, color='#32729b')
    colors = {'opposite':'#286b96', 'fixed_gap':'#d46c31'}
    folds = {}
    for col, model in enumerate(['wave_pinning', 'goryachev']):
        for row, key in enumerate(['peak_fractional_loss_rate_10_to_30', 'fractional_loss_rate_10_to_30']):
            ax = fig.add_subplot(grid[row+1, col])
            for geo, ls in [('opposite', '-'), ('fixed_gap', '--')]:
                selected = sorted([r for r in rows if r['model']==model and r['geometry']==geo], key=lambda r:float(r['length']))
                x = np.array([float(r['length']) for r in selected])
                y = np.array([float(r[key]) for r in selected])
                assert len(x)==18 and np.isfinite(y).all() and (y>0).all()
                label = 'Opposite end' if geo=='opposite' else 'Fixed gap 3'
                ax.plot(x, y, ls, marker='o', ms=3.2, lw=1.7, color=colors[geo], label=label)
                folds[f'{model}/{geo}/{key}'] = float(y[0]/y[-1])
            ax.set(yscale='log', xlim=(8, 40), xticks=[9, 15, 21, 27, 33, 39])
            ax.grid(alpha=.18)
            ax.text(-.14, 1.02, 'ABCD'[row*2+col], transform=ax.transAxes, weight='bold', fontsize=12)
            if row==0:
                ax.set_title('Wave-pinning' if col==0 else 'Goryachev', fontsize=12, weight='bold')
                ax.legend(frameon=False, fontsize=9)
            else:
                ax.set_xlabel('Cell length L (native units)')
                fold=folds[f'{model}/fixed_gap/{key}']
                ax.text(.05, .08, f'Fixed-gap slowdown: {fold:.1f}×', transform=ax.transAxes, fontsize=9,
                        color=colors['fixed_gap'])
            if col==0:
                ax.set_ylabel(('Peak-height' if row==0 else 'Active-amount')+' loss rate\n(initial fraction / native time)')
    fig.text(.12, .038, 'Strong sink (k = 1000); width 1. Rates average 10–30% loss in the fixed original pole core.\n'
             'Native kinetics; no reversal cue. Diagram illustrates geometry and is not to scale.', fontsize=9, color='#52616c')
    outputs=[]
    for ext in ['png', 'pdf', 'svg']:
        p=OUT/f'sink_distance_and_resistance.{ext}'
        fig.savefig(p, dpi=220, facecolor='white'); outputs.append(p)
    plt.close(fig)
    (DATA/'audit.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),
        producer='finished', passed=True, source=str(source.relative_to(ROOT)), source_sha256=sha(source),
        source_generated_at=provenance['generated_at'], source_trials=len(all_rows), selected_rows=len(rows),
        fixed_gap=3, sink_strength=1000, folds_L9_to_L39=folds,
        renderer_sha256=sha(Path(__file__)), outputs={str(p.relative_to(ROOT)):sha(p) for p in outputs}), indent=2)+'\n')
    print(json.dumps({'figure':str(outputs[0]), 'selected_rows':len(rows), 'folds':folds}, indent=2))


if __name__=='__main__':
    main()
