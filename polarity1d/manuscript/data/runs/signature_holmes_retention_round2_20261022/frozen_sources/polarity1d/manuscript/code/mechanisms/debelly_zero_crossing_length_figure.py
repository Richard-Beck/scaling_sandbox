"""Render corrected MCA crossing times from independently audited native arrays."""
from pathlib import Path
import csv
from datetime import datetime, timezone
import hashlib
import json
import textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'data/runs/debelly_corrected_advection_20261008'
FIG = ROOT / 'figures/deBelly_mechanism'
STEM = 'debelly_zero_crossing_vs_length'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    audit_path = SOURCE / 'independent_audit.json'
    audit = json.loads(audit_path.read_text())
    assert audit['passed'] and audit['scientific_sources_unchanged']
    cases = json.loads((SOURCE / 'cases.json').read_text())
    archives = {a['length']:a for a in audit['archives']}
    rows = []
    for case in sorted(cases, key=lambda c:c['length']):
        L = case['length']
        assert case['eligible']
        path = SOURCE / f'L{L:g}_native.npz'
        summary = json.loads((SOURCE / f'L{L:g}_summary.json').read_text())
        sha = digest(path)
        assert sha == summary['native_sha256'] == archives[L]['sha256']
        assert case['checkpoint_sha256'] == summary['checkpoint_sha256'] == archives[L]['checkpoint_sha256']
        with np.load(path) as z:
            t,x,b = z['time'],z['x'],z['field_mca_bound']
            assert np.all(np.isfinite(b)) and np.all(b >= 0)
            p = np.sum(b*(2*x/L-1),axis=1)/np.sum(b,axis=1)
            assert p[0] > 0 and np.all(np.diff(t) > 0)
            i = int(np.flatnonzero(p <= 0)[0])
            assert i > 0 and p[i-1] > 0
            bracket = [float(t[i-1]),float(t[i])]
            assert bracket == summary['crossing_bracket']
            assert bracket[1] == summary['crossing_time'] and bracket[1]-bracket[0] == .5
            interpolated = float(np.interp(0,p[i-1:i+1][::-1],t[i-1:i+1][::-1]))
        rows.append(dict(length_um=L,first_nonpositive_time_s=bracket[1],bracket_start_s=bracket[0],
                         bracket_end_s=bracket[1],interpolated_zero_time_s=interpolated,
                         native_sha256=sha,checkpoint_sha256=case['checkpoint_sha256']))
    assert len(rows) == audit['completed_lengths'] == 9
    FIG.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax = plt.subplots(figsize=(7.6,5.8))
    lengths = np.array([a['length_um'] for a in rows])
    times = np.array([a['first_nonpositive_time_s'] for a in rows])
    lower = np.array([a['bracket_start_s'] for a in rows])
    ax.errorbar(lengths,times,yerr=np.vstack([times-lower,np.zeros(len(times))]),
                color='#0072B2',marker='o',ms=6,lw=2,capsize=3)
    ax.set(xlabel='Cell length (µm)',ylabel='Time to bound-MCA zero crossing (s)',xticks=lengths,xlim=(22,40))
    ax.grid(alpha=.18)
    ax.set_title('deBelly: time to MCA zero crossing',loc='left',weight='bold',pad=32)
    ax.text(0,1.025,'Continuous rightward gradient · fixed slope 0.24',transform=ax.transAxes,
            fontsize=11,color='#43525B')
    ax.annotate(f'{times[0]:g} s',(lengths[0],times[0]),xytext=(10,10),textcoords='offset points',fontsize=10)
    ax.annotate(f'{times[-1]:g} s',(lengths[-1],times[-1]),xytext=(-10,-19),textcoords='offset points',ha='right',fontsize=10)
    caption = ('First recorded nonpositive normalized bound-MCA first moment, measured from cue onset. '
               'All nine lengths start from exact corrected relative-20% patch/basal-polarized checkpoints. '
               'The gradient remains on; its peak increases with length at fixed slope 0.24. '
               'Bars show the preceding 0.5 s crossing bracket and are smaller than the markers; connecting lines guide the eye. '
               'This memory landmark does not establish retained reversal after withdrawal. '
               'Independent local endpoint commands; 201 cells, 0.5 s native step. '
               'Only the audited 23–39 µm mechanism trajectories are included.')
    fig.subplots_adjust(left=.14,right=.96,top=.82,bottom=.31)
    fig.text(.14,.205,'\n'.join(textwrap.wrap(caption,100)),fontsize=8.5,va='top',linespacing=1.4)
    for ext in ['png','pdf','svg']:
        fig.savefig(FIG / f'{STEM}.{ext}',dpi=220)
    plt.close(fig)
    with (FIG / f'{STEM}.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    stamp = datetime.now(timezone.utc).isoformat()
    (FIG / f'{STEM}_caption.txt').write_text(caption+f'\nGenerated {stamp}; simulation producer finished.\n')
    provenance = dict(generated_at=stamp,producer='finished',passed=True,
                      definition='First recorded bound-MCA polarity <= 0 under continuous cue',
                      native_crossing_brackets_independently_verified=True,
                      source=str(SOURCE.relative_to(ROOT)),source_audit_sha256=digest(audit_path),
                      renderer_sha256=digest(Path(__file__)),measurements=rows,
                      artifact_sha256={ext:digest(FIG / f'{STEM}.{ext}') for ext in ['png','pdf','svg','csv']})
    (FIG / f'{STEM}_audit.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps([{'length_um':a['length_um'],'time_s':a['first_nonpositive_time_s']} for a in rows]))

if __name__ == '__main__':
    main()
