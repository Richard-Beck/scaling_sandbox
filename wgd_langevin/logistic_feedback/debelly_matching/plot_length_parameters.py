"""Plot independent length-specific fits at fixed De Belly feedback (1x)."""
from pathlib import Path
import argparse
import hashlib
import json
from datetime import datetime, timezone

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE/'results_self_activation'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['self_activation', 'ungated'], default='self_activation')
    args = parser.parse_args()
    out = HERE/('results_'+args.model)
    source = out/'parameter_mapping.csv'
    frame = pd.read_csv(source)
    frame = frame[frame.variant == 'beta_1x'].sort_values('length_um')
    assert frame.length_um.tolist() == [17., 21., 25., 29., 33., 39.]
    fields = [('tau_p', 'Rac relaxation', 's'), ('tau_v', 'MCA relaxation', 's'),
              ('alpha', 'Rac → MCA coupling', '1/s'), ('beta', 'MCA → Rac coupling', '1/s'),
              ('ps', 'Rac saturation scale', 'centered polarity'),
              ('vs', 'MCA saturation scale', 'centered polarity'),
              ('a', 'Self-activation amplitude', 'polarity/s'),
              ('K', 'Self-activation half-saturation', 'centered polarity'),
              ('cue_gain', 'Projected cue gain', 'polarity / (input · s)')]
    fig, axes = plt.subplots(3, 3, figsize=(13, 10), layout='constrained')
    for ax, (key, title, unit) in zip(axes.flat, fields):
        if key not in frame:
            ax.axis('off')
            ax.text(.5, .5, f'{key}: absent in this model', ha='center', va='center',
                    transform=ax.transAxes, color='#777777', fontsize=12)
            continue
        polarized = frame[frame.length_um >= 21]
        unpolarized = frame[frame.length_um == 17]
        ax.plot(polarized.length_um, polarized[key], '-o', color='#0072b2', ms=5,
                label='Polarized reference cases')
        ax.scatter(unpolarized.length_um, unpolarized[key], marker='o', s=45,
                   facecolors='none', edgecolors='#777777', label='17 µm: unpolarized reference')
        ref = frame[frame.length_um == 29]
        ax.scatter(ref.length_um, ref[key], marker='*', s=150, color='#d55e00',
                   edgecolors='white', linewidth=.5, zorder=4, label='29 µm: multistart reference')
        ax.set(title=f'{key}: {title}', xlabel='Length (µm)', ylabel=unit)
        ax.set_xticks(frame.length_um)
        ax.grid(alpha=.18)
        if key in ['alpha', 'beta', 'a', 'cue_gain']:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=False)
        ax.set_ylim(0, 1.12*frame[key].max())
        if args.model == 'ungated' and key in ['ps', 'vs']:
            ax.axhline(2., color='#777777', ls=':', lw=1)
            ax.text(.97, .93, 'Fit upper bound = 2', ha='right', va='top',
                    transform=ax.transAxes, fontsize=8, color='#777777')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=3, fontsize=9)
    model_label = 'self-activation model, n = 2' if args.model == 'self_activation' else 'ungated model, no polarity self-activation'
    count = 'nine' if args.model == 'self_activation' else 'seven'
    fig.suptitle(f'Length-specific parameters: {model_label}, De Belly feedback = 1x\n'
                 f'All {count} parameters refitted jointly; other lengths use one start. Lines guide the eye, not a scaling law.', fontsize=13)
    fig.savefig(out/'parameters_vs_length.png', dpi=200)
    fig.savefig(out/'parameters_vs_length.pdf')
    plt.close(fig)
    frame.to_csv(out/'parameters_vs_length.csv', index=False)
    rows = []
    for key, title, unit in fields:
        if key not in frame:
            continue
        low = float(frame.loc[frame.length_um == 21, key].iloc[0])
        high = float(frame.loc[frame.length_um == 39, key].iloc[0])
        rows.append(dict(parameter=key, value_21um=low, value_39um=high,
                         percent_change=100*(high/low-1)))
    pd.DataFrame(rows).to_csv(out/'parameter_endpoint_changes.csv', index=False)
    (out/'parameters_vs_length_provenance.json').write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(), variant='beta_1x', model=args.model,
        data_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note='Independent descriptive estimates; no confidence intervals. 17 um is unpolarized.'), indent=2)+'\n')
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
