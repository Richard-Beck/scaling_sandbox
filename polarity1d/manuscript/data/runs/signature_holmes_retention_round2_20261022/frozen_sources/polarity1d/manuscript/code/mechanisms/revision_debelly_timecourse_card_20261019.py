"""Native DB width-pair card; compact CSV-only redraw, no integrations."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'
EXCESS=ROOT/'data/runs/revision_compensation_excess_baseline_20261019'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/debelly_width_timecourse_card_20261019'
SOURCES={}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def rows(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return list(csv.DictReader(p.open()))
def csvwrite(p,rr):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)


def extract():
    m=read(BASE/'manifest.json');em=read(EXCESS/'manifest.json');assert em['native_manifest_sha256']==sha(BASE/'manifest.json')
    data=[];protocols=[]
    for c in m['cases']:
        if c['model']!='debelly' or c['length'] not in [18.,39.] or c['width'] not in [1.5,5.] or c['case_kind']!='constant_on':continue
        result=read(BASE/'workers'/f'{c["slot"]:04d}'/'result.json');assert result['case']==c
        derived=read(EXCESS/'workers/debelly/result.json');off=next(w for w in result['withdrawals'] if w['duration']==2000.)
        protocols.append(dict(length=c['length'],width=c['width'],requested_native_peak=c['response_amplitude'],
            native_pulse_duration=2000.,native_off_duration=10000.,scalar_orientation_success=off['orientation_success'],
            scalar_full_strength_success=off['full_polarity_success'],native_field_tail_pass=off['field_tail_pass']))
        # ON plus the exact 2000-pulse OFF branch; derived phase indices are explicit.
        j=next(i for i,w in enumerate(result['withdrawals']) if w['duration']==2000.)+1
        for phase,item,index in [('cue_on',result['on'],0),('after_withdrawal',off,j)]:
            path=BASE/item['archive'];assert sha(path)==item['archive_sha256'];SOURCES[str(path.relative_to(ROOT))]=item['archive_sha256']
            dp=EXCESS/'workers/debelly'/f'{c["slot"]:04d}_{index:02d}_{phase}.npz'
            assert sha(dp)==derived['derived_artifact_hashes'][str(dp.relative_to(EXCESS))];SOURCES[str(dp.relative_to(ROOT))]=sha(dp)
            with np.load(path) as z:
                t=z['time'].copy();rac=z['field_rac'][:,-1].copy();tension=z['field_tension'][:,-1].copy();disp=z['field_right_boundary_displacement'][:,0].copy()
            with np.load(dp) as z:
                np.testing.assert_array_equal(t,z['time'])
                values={k:z[k].copy() for k in ['organization_core_current_minimum_excess_fraction','organization_core_current_sample_quantile_0p05_excess_fraction','organization_core_frozen_initial_minimum_excess_fraction','organization_core_raw_fraction','organization_core_current_minimum_excess_amount']}
            for k,time in enumerate(t):
                data.append(dict(length=c['length'],width=c['width'],phase=phase,native_time=float(time),challenged_edge_rac=float(rac[k]),
                    challenged_edge_tension=float(tension[k]),right_endpoint_displacement=float(disp[k]),
                    conditional_tension_factor=float(1/(1+tension[k]**2/10.)),
                    **{name:float(v[k]) for name,v in values.items()}))
    assert len(protocols)==4
    return data,protocols


def card(data,length):
    plt.rcParams.update({'font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,aa=plt.subplots(1,3,figsize=(7.2,3.45));fig.subplots_adjust(left=.085,right=.96,bottom=.25,top=.67,wspace=.53)
    specs=[('challenged_edge_rac','Challenged-edge Rac','Native concentration'),('right_endpoint_displacement','Right endpoint response','Displacement (µm)'),('organization_core_current_minimum_excess_fraction','Localized incumbent','Excess signal / initial')]
    for ax,(key,title,ylabel) in zip(aa,specs):
        for width,color in [(1.5,'#444444'),(5.,'#9961a9')]:
            rr=[r for r in data if r['length']==length and r['width']==width and (r['phase']=='cue_on' or r['native_time']>0)]
            times=[r['native_time']+(2000. if r['phase']=='after_withdrawal' else 0.) for r in rr]
            ax.plot(times,[r[key] for r in rr],color=color,lw=1.7,label=f'Width {width:g}')
        ax.set(title=title,xlabel='Native elapsed time (s)',ylabel=ylabel,xlim=(0,12000),xscale='symlog');ax.set_xscale('symlog',linthresh=500.)
        ax.set_xticks([0,500,2000,12000]);ax.set_xticklabels(['0','500','2000','12000']);ax.axvline(2000,color='#888',ls=':',lw=.8)
        ax.grid(alpha=.15);ax.tick_params(labelsize=8)
    aa[-1].axhline(.5,color='#888',lw=.8,ls='--');aa[-1].set_ylim(-.03,1.05)
    fig.suptitle('More local activation need not mean earlier incumbent loss',fontsize=11,y=.98)
    fig.text(.5,.865,f'Corrected linear deBelly L{length:g}; cue off at 2000 s, followed by 10000 s off',ha='center',fontsize=8.5)
    fig.legend(*aa[0].get_legend_handles_labels(),loc='upper center',ncol=2,bbox_to_anchor=(.50,.835),frameon=False)
    fig.text(.5,.045,'Time axis is linear through 500 s, logarithmic thereafter; dotted line marks withdrawal.',ha='center',fontsize=8)
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--extract',action='store_true');args=parser.parse_args()
    for n in ['supporting/data','supporting/vectors','supporting/metadata']:(OUT/n).mkdir(parents=True,exist_ok=True)
    table=OUT/'supporting/data/native_timecourses.csv';protocol=OUT/'supporting/data/native_protocols.csv'
    if args.extract:
        data,pp=extract();csvwrite(table,data);csvwrite(protocol,pp)
        (OUT/'supporting/metadata/native_evidence_sources.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    else:
        data=rows(table);pp=rows(protocol)
        for r in data:
            for k in r:
                if k!='phase':r[k]=float(r[k])
    for length,name in [(39.,'debelly_L39_width_card'),(18.,'debelly_L18_width_support')]:
        fig=card(data,length);fig.savefig(OUT/(name+'.png'),dpi=220)
        for ext in ['pdf','svg']:fig.savefig(OUT/'supporting/vectors'/(name+'.'+ext))
        plt.close(fig)
    caption='''# Native corrected-deBelly width pair

Exact native fields at L39 compare nominal right-hand cue widths 1.5 and 5, both requested peak3.6, with identical initial corrected independent-endpoint state. The card shows the complete2000-second pulse and10000 cue-free seconds. Its time axis is linear through500seconds and logarithmic thereafter; the dotted line marks withdrawal. The compact CSV retains the complete ON/OFF native sampling and native amplitudes. L18 is exported as a supporting contrast.

The wider input produces more challenged-edge Rac but delays the early endpoint response and the first50% loss of localized initial-core excess. The wider case later catches up and exceeds the narrow endpoint displacement; it is not uniformly less mobile. At L39,t100, edgeRac is4.0137 versus2.7684, endpoint displacement0.5087 versus1.1319µm, and edge tension7.2824 versus3.0185. In the native command law, the tension denominator is1+tension²/10. Both Rac switches are saturated at this saved state, giving a conditional broad/narrow command ratio≈0.303. That command is distinct from the actual coupled membrane velocity and from accumulated endpoint displacement, which includes the accepted per-update boundary relaxation.

Localized incumbent signal integrates positive Rac above the current full-profile minimum in the fixed, fractional initial organization core; deBelly uses nominal chemical-solver coordinates. Both ROIs and background sensitivity remain in the comprehensive source tables. Primary50% loss entries are(370,372.5] versus(447.5,450]seconds. These are sampled observation brackets, not confidence intervals. Both2000-second width cases pass the original retained scalar-strength and native field-tail criteria; those criteria are separate from regional organization loss and do not by themselves prove whole-field replacement.

The actuator remains fixed in nominal coordinates; a wider input also changes integrated input. Physical cue-sample extrema and boundary displacement are exported separately in the native snapshot source. The native Rac production depends on MCA and Rho, while mechanical response depends on Rac and tension. The observed coupled sequence and saved-state conditional law explain why more Rac need not imply faster regional reorganization. No evolved mediation experiment, unique molecular explanation or quantitative causal share of latency is claimed. N=201 and the original0.5-second staggered update are retained without a new convergence claim.
'''
    (OUT/'CAPTION.md').write_text(caption)
    (OUT/'README.md').write_text('# Native explanatory pair: corrected linear deBelly\n\n![L39 native timecourse](debelly_L39_width_card.png)\n\n[PDF](supporting/vectors/debelly_L39_width_card.pdf) · [Caption](CAPTION.md) · [Complete compact native timecourses](supporting/data/native_timecourses.csv) · [L18 contrast](debelly_L18_width_support.png).\n\nThe display is a candidate explanatory card, not an isolated mechanical ablation or a final figure composition.\n')
    closure=dict(command='python polarity1d/manuscript/code/mechanisms/revision_debelly_timecourse_card_20261019.py',working_directory='repository root',native_archives_required_for_redraw=False,
        source_hashes={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},compact_inputs={str(p.relative_to(ROOT)):sha(p) for p in [table,protocol]})
    (OUT/'supporting/metadata/redraw_closure.json').write_text(json.dumps(closure,indent=2)+'\n')
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(source_hashes=SOURCES,redraw_source_sha256=sha(Path(__file__)),timecourse_rows=len(data),
        artifact_hashes={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='validation.json'}),indent=2)+'\n')
    print(str(OUT/'debelly_L39_width_card.png'))


if __name__=='__main__':main()
