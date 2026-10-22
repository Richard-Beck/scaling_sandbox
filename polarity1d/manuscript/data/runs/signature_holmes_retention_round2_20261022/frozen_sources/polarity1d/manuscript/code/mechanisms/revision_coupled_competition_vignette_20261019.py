"""Result-led, native-field vignette of coupled activation and retention."""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import revision_width_fullfield_review_20261019 as native

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/coupled_competition_vignette_20261019'
PRIOR=OUT.parent/'finite_response_candidate_20261019/figure2_organization.png'
COLORS={'narrow':'#2171b5','broad':'#983ca0'}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def foundation():
    fig=plt.figure(figsize=(15,9.4))
    fig.text(.06,.945,'Larger incumbents: more material, organized differently',fontsize=19,weight='bold')
    fig.text(.06,.90,'Prior all-five evidence retained unchanged: excess amount, physical extent and native peak concentration.',fontsize=11)
    ax=fig.add_axes([.025,.34,.95,.51]);ax.imshow(plt.imread(PRIOR));ax.axis('off')
    fig.text(.065,.265,'Replacing an incumbent requires more than activating a new region',fontsize=17,weight='bold')
    fig.text(.065,.207,'A wider cue can recruit more distal activity and weaken the old region more, yet leave a challenger that fails after withdrawal.',fontsize=11.5)
    fig.text(.065,.153,'The next sheet shows this directly in Goryachev L9: same basal, same requested peak, same pulse; only physical width changes.',fontsize=11)
    fig.text(.065,.104,'This illustrates coupled competition. Larger-cell and WP examples provide context; it is not a universal width rule or a unique molecular signature.',fontsize=10)
    fig.text(.065,.055,'The organization comparison is observational. The challenge changes both spatial support and integrated input; it does not isolate incumbent height or width causally.',fontsize=9)
    return fig


def vignette(pair):
    fig=plt.figure(figsize=(15,9.4))
    fig.text(.06,.945,'Greater incumbent weakening can accompany challenger failure',fontsize=18,weight='bold')
    fig.text(.06,.902,'Goryachev L9: exact shared basal, requested native peak 6.4, pulse 74, widths 1.8 / 3; then 10000 cue-free units.',fontsize=11)
    heat_axes=[fig.add_axes([.065,.43,.235,.36]),fig.add_axes([.365,.43,.235,.36])]
    image=None
    for ax,d in zip(heat_axes,pair):
        c=d['case'];on,off=d['phases']['cue_on'],d['phases']['off']
        ids=np.unique(np.r_[np.arange(0,len(on['time']),5),len(on['time'])-1])
        offids=np.flatnonzero(off['time']<=1000.)[::10]
        offids=offids[off['time'][offids]>0]
        time=np.r_[on['time'][ids],74.+off['time'][offids]]
        active=np.vstack([on['field_active'][ids],off['field_active'][offids]])
        image=ax.pcolormesh(on['x'],time,active,shading='nearest',cmap='magma',vmin=0,vmax=28.,rasterized=True)
        ax.set_yscale('symlog',linthresh=74.)
        ax.set_ylim(1074,0);ax.set_xlim(0,9)
        ax.set_yticks([0,37,74,174,1074],['0','37','74: off','174','1074'])
        ax.axhline(74,color='white',ls='--',lw=1)
        ax.plot([9-c['width'],9-c['width']],[0,74],color='cyan',ls=':',lw=1.4)
        ax.set_xlabel('Native interval position',fontsize=9)
        ax.set_ylabel('Time from cue onset (native)',fontsize=9)
        ax.set_title(f'{"A" if c["width_label"]=="narrow" else "B"}  w = {c["width"]:g}: {"new region retained" if c["width_label"]=="narrow" else "incumbent returns"}',loc='left',fontsize=11,weight='bold')
        ax.tick_params(labelsize=8)
    cb=fig.colorbar(image,cax=fig.add_axes([.615,.43,.012,.36]));cb.set_label('Native active concentration',fontsize=9);cb.ax.tick_params(labelsize=8)
    fig.text(.065,.371,'Actual saved fields; time is linear through 74 and logarithmic later. Cyan marks the physical cue-support edge during stimulation.',fontsize=8.5)

    ax=fig.add_axes([.735,.43,.215,.36])
    for d in pair:
        c=d['case'];core=d['masks']['incumbent_core_mask'];roi=d['masks']['cue_footprint_mask']
        for phase,offset in (('cue_on',0.),('off',74.)):
            p=d['phases'][phase];keep=p['time']<= (74 if phase=='cue_on' else 1000)
            tt=p['time'][keep]+offset;u=p['field_active'][keep]
            ax.plot(tt,u[:,core].max(axis=1),color=COLORS[c['width_label']],ls='--',lw=1.4)
            ax.plot(tt,u[:,roi].max(axis=1),color=COLORS[c['width_label']],lw=1.7)
    ax.set_xscale('symlog',linthresh=74);ax.set_xlim(0,1074);ax.set_ylim(0,29)
    ax.set_xticks([0,74,174,1074],['0','74','174','1074'])
    ax.axvline(74,color='.5',lw=1)
    ax.set_title('C  Both challengers activate',loc='left',fontsize=11,weight='bold')
    ax.set_xlabel('Time from cue onset (native)',fontsize=9);ax.set_ylabel('Native active peak',fontsize=9)
    ax.legend(handles=[Line2D([],[],color=COLORS['narrow'],label='w1.8'),Line2D([],[],color=COLORS['broad'],label='w3'),
        Line2D([],[],color='black',label='cue ROI'),Line2D([],[],color='black',ls='--',label='original core')],fontsize=7.5,frameon=False)
    ax.tick_params(labelsize=8)

    ax=fig.add_axes([.075,.145,.21,.18])
    fractions=[d['summary']['incumbent_fraction_at_withdrawal'] for d in pair]
    ax.bar([0,1],fractions,color=[COLORS['narrow'],COLORS['broad']],width=.6)
    ax.set_ylim(0,1);ax.set_xticks([0,1],['w1.8','w3']);ax.set_ylabel('Original-core active amount / initial',fontsize=8.5)
    ax.set_title('D  More original-region loss',loc='left',fontsize=10.5,weight='bold')
    for i,value in enumerate(fractions):ax.text(i,value+.04,f'{value:.3f}',ha='center',fontsize=9)
    ax.tick_params(labelsize=8)
    ax=fig.add_axes([.38,.145,.21,.18])
    amounts=[d['summary']['fixed_right_half_active_amount_at_withdrawal'] for d in pair]
    ax.bar([0,1],amounts,color=[COLORS['narrow'],COLORS['broad']],width=.6)
    ax.set_ylim(0,16);ax.set_xticks([0,1],['w1.8','w3']);ax.set_ylabel('Native active amount in fixed right half',fontsize=8.5)
    ax.set_title('E  More distal activity; lower peak',loc='left',fontsize=10.5,weight='bold')
    for i,value in enumerate(amounts):ax.text(i,value+.5,f'{value:.3f}',ha='center',fontsize=9)
    ax.tick_params(labelsize=8)
    ax=fig.add_axes([.70,.135,.26,.20]);ax.axis('off')
    ax.text(0,1,'Native mechanism and hypothesis',weight='bold',fontsize=10.5,va='top')
    ax.text(0,.73,r'$J=(a e_c u^2+b e_c u+S)v-cu$',fontsize=14,va='top')
    ax.text(0,.46,'Stimulation recruits from inactive v.\nWithdrawal removes S; nonlinear local\nmaintenance and transport continue.\nDistributed activation may fail to sustain\na concentrated challenger in this case.',fontsize=9,va='top',linespacing=1.35)
    fig.text(.065,.055,'Supported: the broader response forms with more distal active amount, then disappears. The shared-pool/local-maintenance account is a hypothesis, not a unique causal attribution.',fontsize=8.5)
    fig.text(.065,.028,'Larger Goryachev and WP width contrasts differ; this is not a universal peak or depletion threshold. Native clocks/input are uncalibrated; N = 201; no spatial-convergence claim.',fontsize=8.5)
    return fig


def main():
    for folder in ('supporting/data','supporting/metadata','supporting/vectors','supporting/viewers'):(OUT/folder).mkdir(parents=True,exist_ok=True)
    cases,sources=native.read_cases();rows=[]
    for d in cases:
        r=native.metrics(d);c=d['case'];u=d['phases']['cue_on']['field_active'][-1];x=d['masks']['x'];core=d['masks']['incumbent_core_mask']
        r.update(fixed_right_half_active_amount_at_withdrawal=float(u[x>=c['length']/2].sum()*c['length']/c['cells']),
                 original_core_peak_at_withdrawal=float(u[core].max()),
                 actual_cue_footprint_active_amount_at_withdrawal=float(u[d['masks']['cue_footprint_mask']].sum()*c['length']/c['cells']))
        d['summary']=r;rows.append(r)
    with (OUT/'supporting/data/coupled_width_readouts.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    pair=sorted([d for d in cases if d['case']['model']=='goryachev' and d['case']['length']==9],key=lambda d:d['case']['width'])
    assert len(pair)==2
    files=[]
    with PdfPages(OUT/'supporting/vectors/figure2_coupled_competition_candidate.pdf') as book:
        for stem,fig in [('figure2_incumbent_foundation',foundation()),('figure2_coupled_competition',vignette(pair))]:
            png=OUT/(stem+'.png');fig.savefig(png,dpi=160);files.append(png)
            svg=OUT/'supporting/vectors'/(stem+'.svg');fig.savefig(svg);files.append(svg)
            book.savefig(fig);plt.close(fig)
    (OUT/'CLAIM.md').write_text('''# One claim, one observable contrast

**Claim supported by the fields:** weakening the original region and generating substantial distal activity do not guarantee that a challenger persists after withdrawal. At Goryachev L9, a broader equal-peak cue produces greater original-core activity loss and more active amount in the same fixed distal half, yet a lower challenger peak and eventual incumbent return. The narrower cue retains the new pole. Both challengers genuinely form; the failed case is not an absence of stimulation response.

**Characterized challenge and observable:** exact shared corrected basal, native requested peak 6.4, physical widths 1.8 / 3, pulse 74, then 10000 cue-free units. Image simultaneous old-region activity, distal activity distribution and their subsequent retention. Native active concentration and physical cell-centre integrals are used; no reporter kinetics or optical-actuator equivalence is assumed. The actual cue footprint changes with width; the main amount comparison uses the identical fixed right-half ROI for both cases.

**Exact supporting observation:** at withdrawal, original-core active fraction is 0.449 / 0.412 for narrow / broad cues; distal-half active amount is 12.428 / 13.466; challenger peak is 12.025 / 7.433. After release, both challenger peaks initially fall; the narrow challenger subsequently grows into the retained pole, while the broad challenger disappears and the original side recovers. All first-stage finite native field-tail tests pass. Actual full-field histories supply this evidence independently from the original scalar-strength gates.

**Mechanistic hypothesis:** the native stimulus × inactive channel redistributes a shared conserved pool. A broader response can contain more active material while being less concentrated; after stimulation ceases, the native concentration-dependent feedback, deactivation and transport determine whether it is maintained. The Goryachev activation law contains quadratic/linear active-dependent recruitment from inactive substrate. These fields are consistent with that coupled account but do not uniquely identify its causal allocation. Larger Goryachev L25 behaves differently: broadening 3 → 5 at pulse 815 lowers the challenger peak but increases distal amount and now favors retention. WP broadening also favors retention in the selected cases. No universal peak or original-region-depletion threshold is asserted.

**How this adds to the size argument:** the unchanged all-five foundation shows that larger cells carry more incumbent material, with WP mainly widening and Goryachev mainly increasing height. The width challenge exposes an additional difficulty beyond lowering the incumbent's activity: a newly activated distribution must persist in competition after cue withdrawal. The size cohorts use different pulses and width intervals, so the comparison is not an isolated causal intervention on size or incumbent height. The existing corrected pulse brackets show the width-effect direction across neighboring lengths, while this sheet explains the selected finite histories.

**Operating-point limit:** prior weakening-time dose sensitivity differs mainly in degree and may reflect tested dose range/substrate depletion. A separate bounded amplitude extension probes below that range. Fractional peak matching is not matched effective actuation. No new all-five factorial or unique mechanism signature is implied.
''')
    (OUT/'README.md').write_text('# Coupled competition: compact Figure 2 review candidate\n\n[Two-sheet PDF](supporting/vectors/figure2_coupled_competition_candidate.pdf) · [One claim / observable / evidence / hypothesis](CLAIM.md) · [Exact eight-case context](supporting/data/coupled_width_readouts.csv)\n\n![Prior all-five foundation](figure2_incumbent_foundation.png)\n\n![Actual coupled width contrast](figure2_coupled_competition.png)\n\nThis is a separate result-led candidate for author steering, not an accepted final layout. Prior artifacts remain unchanged. Heatmaps show selected actual saved frames through 1000 units after withdrawal; retained outcome and source fields extend through 10000 cue-free units. Kymographs are linear through withdrawal at 74 and logarithmic later. Source fields were saved every 0.1 unit during stimulation and 0.25 unit after withdrawal; displayed heatmaps use 0.5 / 2.5-unit saved-frame intervals.\n')
    sources[str(PRIOR.relative_to(ROOT))]=sha(PRIOR);sources[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    artifacts={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='validation.json'}
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(source_hashes=sources,artifact_hashes=artifacts,
        native_context_cases=8,main_cases=2,initial_cue_duration=74.,displayed_off_horizon=1000.,retained_off_horizon=10000.),indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),pages=2),indent=2))


if __name__=='__main__':main()
