"""Full F4 organization/amount/coherence/width observed-window supplement.

Every model, size and strength is retained; native numerical evidence unchanged.
Instantaneous widths remain conditional descriptive measurements, not lifetimes.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
FIG = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
OUT = FIG / 'postburn_signature_candidate_20261019'
CONDITIONS = FIG / 'broad_strength_campaign_20261019/review_complete_scientific/supporting/data/condition_window_snapshot.csv'
MODELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev', holmes_model3='Holmes 3', spring_ring='Spring ring')
GRID = np.array([.0001, .0003, .001, .003, .01, .03, .1, .3, 1.])
COLORS = ('#2676a6', '#28866b', '#c77722')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    parser.add_argument('--conditions', type=Path, default=CONDITIONS)
    parser.add_argument('--spatial', type=Path)
    args = parser.parse_args(); out = args.out
    spatial_path = args.spatial or out / 'supporting/data/postburn_spatial_condition_windows.csv'
    scalar = pd.read_csv(args.conditions); spatial = pd.read_csv(spatial_path)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'svg.fonttype':'none', 'pdf.fonttype':42})
    artists = []; pages = []
    with PdfPages(out / 'supporting/vectors/all_models_sizes_strengths_organization_windows.pdf') as pdf:
        for window, window_label in (('full','Whole observed run'), ('first_half','First half'), ('last_half','Second half')):
            for sheet, title in (
                ('organization','Organization, contrast and global magnitude'),
                ('amount','Active amount, coherence and spatial harmonics'),
                ('width','Instantaneous region widths with qualification coverage')):
                fig, axes = plt.subplots(4,3,figsize=(13.3,14.1))
                fig.subplots_adjust(left=.078,right=.98,bottom=.14,top=.84,hspace=.65,wspace=.29)
                fig.suptitle(title + ': all sizes and strengths',x=.078,y=.975,ha='left',fontsize=17)
                observed_limits={'full':'chemicals0–32000; Spring0–3200',
                    'first_half':'chemicals0–16000; Spring0–1600',
                    'last_half':'chemicals16000–32000; Spring1600–3200'}[window]
                fig.text(.078,.941,window_label + ' · ' + observed_limits + '. Half-window midpoint shared.',fontsize=10.3)
                fig.text(.078,.914,'Size colors: blue smallest, teal middle, orange largest. Squares n=12; circles n=2. Whiskers: observed seed range, not confidence intervals.',fontsize=9.7)
                key = ('solid: raw multi / noisy low contrast / noisy R     dashed: contrast-qualified multi / matched mean-control low contrast / mean-control R' if sheet=='organization'
                    else 'solid: active amount A / raw q₁ / noisy contrast     dashed: coherence Q / q₂ / matched mean-control contrast' if sheet=='amount'
                    else 'solid: qualified multi occupancy / largest width     dashed: qualified-frame coverage / second-largest width')
                fig.text(.078,.881,key,fontsize=9.1)
                for row,(model,label) in enumerate(MODELS.items()):
                    for size,color in zip((1.,2.,3.) if model=='spring_ring' else (20.,40.,60.),COLORS):
                        n = scalar[(scalar.model==model)&(scalar['size']==size)&(scalar.kind=='noise')&(scalar.window==window)].set_index('fraction').reindex(GRID)
                        mean = scalar[(scalar.model==model)&(scalar['size']==size)&(scalar.kind=='mean')&(scalar.window==window)].set_index('fraction').reindex(GRID)
                        s = spatial[(spatial.model==model)&(spatial['size']==size)&(spatial.kind=='noise')&(spatial.window==window)].set_index('fraction').reindex(GRID)
                        if sheet=='organization':
                            specifications=(((n,'multiregion_fraction','-'),(n,'qualified_multiregion_fraction','--')),
                                            ((n,'low_contrast_fraction','-'),(mean,'low_contrast_fraction','--')),
                                            ((n,'mean_R','-'),(mean,'mean_R','--')))
                        elif sheet=='amount':
                            specifications=(((n,'mean_A','-'),(n,'mean_Q','--')),
                                            ((n,'mean_raw_q1','-'),(n,'mean_q2','--')),
                                            ((n,'mean_contrast_ratio','-'),(mean,'mean_contrast_ratio','--')))
                        else:
                            specifications=(((s,'qualified_multiregion_fraction_profile_cadence','-'),(s,'qualified_profile_fraction','--')),
                                            ((s,'mean_largest_physical','-'),(s,'mean_second_physical','--')),
                                            ((s,'mean_largest_fraction','-'),(s,'mean_second_fraction','--')))
                        for column,specs in enumerate(specifications):
                            ax=axes[row,column]
                            for table,metric,style in specs:
                                y=table[metric].to_numpy(float);ax.plot(GRID,y,color=color,ls=style,lw=1.45)
                                noisy=table is not mean
                                if noisy:
                                    ax.vlines(GRID,table[metric+'_seed_min'],table[metric+'_seed_max'],color=color,lw=.7,alpha=.4)
                                    counts=table['seeds'] if sheet=='width' else table['available_cases']
                                    contributors=table[metric+'_contributing_seeds'] if metric+'_contributing_seeds' in table else counts
                                    for f,value,count,contributing in zip(GRID,y,counts,contributors):
                                        marker='^' if contributing<count else 's' if count==12 else 'o'
                                        if np.isfinite(value):ax.plot(f,value,color=color,marker=marker,mfc=color if count==12 else 'white',ms=3.5,ls='none')
                                artists.append(dict(sheet=sheet,window=window,model=model,size=size,metric=metric,style=style,
                                    values=[float(v) if np.isfinite(v) else None for v in y]))
                    for column,ax in enumerate(axes[row]):
                        ax.set_xscale('log');ax.set_xlim(.00007,1.4);ax.set_xticks([.0001,.001,.01,.1,1.],['10⁻⁴','10⁻³','10⁻²','10⁻¹','1'])
                        ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8.5);ax.grid(axis='y',color='#e6eaec',lw=.7)
                        ax.set_xlabel('Whole-signal strength factor f',fontsize=9.3)
                        if sheet=='organization':
                            titlecol=('Raw / qualified multiple regions','Low contrast / matched mean control','First-moment magnitude R')[column]
                            if column<2:ax.set_ylim(-.03,1.03)
                            else:ax.set_ylim(bottom=0)
                            ylabel='Fraction of scalar times' if column<2 else 'Mean R'
                        elif sheet=='amount':
                            titlecol=('Amount A / coherence Q','Raw first / second harmonic coherence','Contrast / matched mean control')[column]
                            if column==1:ax.set_ylim(-.03,1.03)
                            elif column==2:ax.set_yscale('symlog',linthresh=.1);ax.set_ylim(bottom=0)
                            else:ax.set_ylim(bottom=0)
                            ylabel=('Normalized quantity','Raw q₁ / q₂','Contrast / basal (symlog)')[column]
                        else:
                            titlecol=('Qualified multi / qualified-frame coverage','Conditional physical region widths','Conditional width / circumference')[column]
                            if column!=1:ax.set_ylim(-.03,1.03)
                            else:ax.set_ylim(bottom=0)
                            ylabel=('Fraction of profile times','Native physical arc width','Fraction of circumference')[column]
                        ax.set_ylabel(ylabel,fontsize=9)
                        ax.set_title(f'{label}: '+titlecol,loc='left',fontsize=9.5)
                    handles=[Line2D([],[],color=c,lw=1.5,label=f'{"C₀" if model=="spring_ring" else "C"}={v:g}') for c,v in zip(COLORS,(1,2,3) if model=='spring_ring' else (20,40,60))]
                    axes[row,0].legend(handles=handles,frameon=False,fontsize=7.5,ncol=3,loc='upper left',bbox_to_anchor=(-.015,1.26),handlelength=1.4,columnspacing=1.1)
                if sheet=='width':
                    fig.text(.078,.09,'Widths join periodic seam-crossing regions: cell count/N × actual C(t). Threshold: current min+half contrast; qualification: contrast≥10% basal. Profile cadence40/4.',fontsize=9.1)
                    fig.text(.078,.062,'Largest: one per qualified frame; second: one per qualified multiregion frame. Equal available-seed means; △ fewer contributing seeds than condition n (exact counts exported).',fontsize=9.1)
                    fig.text(.078,.034,'Missing conditional widths stay unavailable. Widths are descriptive geometry, outside the ten-metric numerical gate; no region births, establishment, lifetimes or travel inferred.',fontsize=9.1)
                else:
                    fig.text(.078,.09,'All scalar times remain in occupancy denominators; cadence5 chemicals/0.5 Spring. Threshold regions may coexist with weak or strong resultants; these are not exclusive phases.',fontsize=9.1)
                    fig.text(.078,.062,'R=AQ exactly per frame; mean R is mean of the paired product, not product of means. q₁=|M₁|/M₀ and q₂=|M₂|/M₀ retain basal coherence differences.',fontsize=9.1)
                    fig.text(.078,.034,'Mean controls match declared ensemble input means, without per-path normalization. Native clocks/laws/units differ; full/half differences do not certify stationarity.',fontsize=9.1)
                stem=f'full_organization_{window}_{sheet}'
                fig.savefig(out/'supporting/viewers'/(stem+'.png'),dpi=140);pdf.savefig(fig);plt.close(fig);pages.append(stem)
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),source=str(Path(__file__).relative_to(ROOT)),source_sha256=sha(__file__),
        scalar_source=str(args.conditions.relative_to(ROOT)),scalar_source_sha256=sha(args.conditions),
        spatial_source=str(spatial_path.relative_to(ROOT)),spatial_source_sha256=sha(spatial_path),
        pdf_pages=9,pages=pages,native_integrations=0,artist_count=len(artists),artists=artists,
        width_weighting='Conditional frame mean then equal available seeds; actual profile frames and contributors exported',
        noise_and_mean_controls=True,all_models_sizes_nine_strengths_three_windows=True)
    (out/'supporting/metadata/FULL_ORGANIZATION_RENDER_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artists'},indent=2))


if __name__=='__main__':main()
