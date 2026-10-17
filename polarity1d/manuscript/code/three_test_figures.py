"""Canonical five-model x three-geometry figures, after complete trial audit."""
from pathlib import Path
import argparse, csv, gzip, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
import figure_builder as f

ROOT=Path(__file__).resolve().parents[1]
TESTS=['fixed_physical','fixed_relative','fixed_slope']
LABELS=['Fixed width 3','Relative width 20%','Constant physical slope']
f.CONDITIONS=TESTS
f.CONDITION_LABELS=LABELS
CORRELATION_RULE = ('All audited scalar successes, irrespective of the separate native-field tail flag; '
                    'at least four distinct lengths and nonconstant durations. '
                    'Censored and preparation-ineligible readouts are excluded, never imputed. '
                    'Ties receive average ranks; no significance or effect-size cutoff.')


def correlation_summary(rows):
    result = []
    for model in f.MODELS:
        for condition in TESTS:
            for gate, _, _, _ in f.GATES:
                chosen = [r for r in rows if (r['model'], r['condition'], r['definition']) ==
                          (model, condition, gate)]
                assert chosen, (model, condition, gate)
                passed = [r for r in chosen if r['status'] == 'success']
                x = np.array([r['length'] for r in passed])
                y = np.array([r['endpoint'] for r in passed])
                assert np.isfinite(x).all() and np.isfinite(y).all()
                assert len(np.unique(x)) == len(x), 'Duplicate successful lengths'
                row = dict(model=model, condition=condition, definition=gate,
                           n_total=len(chosen), n_success=len(passed), n_correlation=len(passed),
                           n_field_tail_passed=sum(r['selected_field_tail_pass'] is True for r in passed),
                           n_field_tail_failed=sum(r['selected_field_tail_pass'] is False for r in passed),
                           n_field_tail_unknown=sum(r['selected_field_tail_pass'] is None for r in passed),
                           n_censored=sum(r['status'] == 'censored' for r in chosen),
                           n_ineligible=sum(r['status'] == 'ineligible' for r in chosen),
                           n_preparation_invalid=sum(r['status'] == 'preparation_invalid' for r in chosen),
                           length_min=float(x.min()) if len(x) else None,
                           length_max=float(x.max()) if len(x) else None,
                           correlation_status='insufficient_successes', spearman_rho=None)
                if len(x) >= 4:
                    row['correlation_status'] = 'constant_duration' if np.ptp(y) == 0 else 'available'
                    if np.ptp(y) > 0:
                        row['spearman_rho'] = float(spearmanr(x, y).statistic)
                        assert np.isfinite(row['spearman_rho']) and abs(row['spearman_rho']) <= 1
                result.append(row)
    return result


def correlation_check(rows, correlations):
    """Independent pairwise midranks and centered dot product, including ties."""
    count = 0
    for summary in correlations:
        if summary['correlation_status'] != 'available':
            continue
        chosen = [r for r in rows if all(r[k] == summary[k] for k in
                  ['model', 'condition', 'definition']) and r['status'] == 'success']
        ranks = []
        for key in ['length', 'endpoint']:
            values = np.array([r[key] for r in chosen])
            rank = np.array([1 + np.sum(values < v) + .5 * (np.sum(values == v) - 1)
                             for v in values])
            ranks.append(rank - rank.mean())
        a, b = ranks
        expected = np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b))
        np.testing.assert_allclose(summary['spearman_rho'], expected, rtol=1e-12, atol=1e-12)
        assert summary['n_correlation'] == len(chosen)
        count += 1
    return dict(passed=True, independently_verified_correlations=count,
                method='Pairwise average ranks and centered dot product, independent of scipy.stats.spearmanr.',
                inclusion_rule=CORRELATION_RULE)


def regression_check(rows,fits):
    count=0
    for fit in fits:
        if fit['fit_status']!='fit':continue
        chosen=[r for r in rows if all(r[k]==fit[k] for k in ['model','condition','definition']) and r['status']=='success' and r['selected_field_tail_pass']]
        x=np.array([r['length'] for r in chosen]);y=np.array([r['endpoint'] for r in chosen])
        slope=np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
        intercept=y.mean()-slope*x.mean();mid=(x.min()+x.max())/2
        np.testing.assert_allclose([slope,intercept,slope*mid/(slope*mid+intercept)],
                                   [fit['slope'],fit['intercept'],fit['normalized_slope']],rtol=1e-8,atol=1e-11)
        if np.ptp(y)>0:
            rsq=1-np.sum((y-slope*x-intercept)**2)/np.sum((y-y.mean())**2)
            np.testing.assert_allclose(rsq,fit['r_squared'],rtol=1e-8,atol=1e-11)
        count+=1
    return dict(passed=True,independently_verified_fits=count,method='Closed-form centered ordinary least squares, independent of plotted lstsq coefficients.')


def main_figure(rows,correlations,manifest,output):
    prep,duration,arrival,prep_path=f.preparation_example(ROOT/'data')
    # This historical example starts from a different preparation; only panel A
    # is reused. Panel B explicitly shows exact shared-checkpoint response cues.
    fig=plt.figure(figsize=(11,7.6))
    grid=fig.add_gridspec(1,3,left=.075,right=.97,top=.94,bottom=.68,wspace=.36)
    for x,label in zip([.075,.385,.695],['Preparation cue','Native active profiles','Measured polarity']):fig.text(x,.967,label,weight='bold',fontsize=9)
    f.measured_row(fig,grid,0,'A','',prep,duration,arrival)
    fig.axes[1].get_legend().set_bbox_to_anchor((.5,-.15),transform=fig.axes[1].transAxes)
    fig.text(.035,.565,'B',weight='bold',fontsize=12)
    q=np.linspace(0,1,401)
    for j,title in enumerate(LABELS):
        ax=fig.add_axes([.075+.31*j,.395,.27,.16])
        for L,color,ls in [(15,'#333333','-'),(30,'#73A6B8','--')]:
            signal=f.cosine(q,3/L if j==0 else .2,end='right') if j<2 else (L/15)*q
            ax.plot(q,signal,color=color,ls=ls,lw=1.4,label=f'L={L}')
        ax.set(xlim=(0,1),ylim=(0,2.1 if j==2 else 1.1),xticks=[0,.5,1],xlabel='Position x / L')
        ax.set_title(title,fontsize=10,loc='left');ax.tick_params(labelsize=7)
        if j==0:ax.set_ylabel('Cue / reference peak',fontsize=8);ax.legend(frameon=False,fontsize=7)
    fig.text(.035,.288,'C',weight='bold',fontsize=12)
    lookup={(r['model'],r['condition'],r['definition']):r for r in correlations}
    matrix=np.array([[lookup[m,g,gate]['spearman_rho'] if lookup[m,g,gate]['spearman_rho'] is not None else np.nan for g in TESTS for gate,_,_,_ in f.GATES] for m in f.MODELS])
    ax=fig.add_axes([.22,.095,.75,.18]);cmap=plt.colormaps['RdBu_r'].copy();cmap.set_bad('#eeeeee')
    im=ax.imshow(np.ma.masked_invalid(matrix),cmap=cmap,vmin=-1,vmax=1,aspect='auto')
    for i,m in enumerate(f.MODELS):
        for j,g in enumerate(TESTS):
            for k,(gate,_,_,_) in enumerate(f.GATES):
                summary=lookup[m,g,gate];v=summary['spearman_rho']
                label=(f'{v:+.2f}' if v is not None else 'Undefined' if summary['correlation_status']=='constant_duration' else 'Unavailable')+f'\nn={summary["n_correlation"]}'
                ax.text(2*j+k,i,label,ha='center',va='center',fontsize=7,color='white' if v is not None and abs(v)>.65 else '#222222')
    ax.set_yticks(range(5),[f.LABELS[m] for m in f.MODELS],fontsize=8)
    ax.set_xticks(range(6),['Orientation','Strength']*3,fontsize=7);ax.tick_params(length=0)
    for j,g in enumerate(LABELS):ax.text((2*j+1)/6,1.07,g,transform=ax.transAxes,ha='center',fontsize=8)
    for border in [1.5,3.5]:ax.axvline(border,color='white',lw=2)
    cax=fig.add_axes([.22,.043,.25,.012])
    colorbar=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,-.5,0,.5,1])
    colorbar.ax.tick_params(labelsize=6)
    cax.text(1.1,.5,'Spearman correlation',transform=cax.transAxes,va='center',fontsize=8)
    f.save(fig,output,'main_polarity')
    return prep_path


def supplement(rows,output,stamp,manifest):
    fig,axes=plt.subplots(5,3,figsize=(14.5,15.5))
    for i,m in enumerate(f.MODELS):
        for j,g in enumerate(TESTS):
            ax=axes[i,j];chosen=sorted([r for r in rows if r['model']==m and r['condition']==g],key=lambda r:r['length'])
            valid=[r['endpoint'] for r in chosen if r['status']=='success'];caps=[r['maximum_tested_duration'] for r in chosen if r['maximum_tested_duration'] is not None];maximum=max(caps)
            for gate,color,marker,_ in f.GATES:
                rs=[r for r in chosen if r['definition']==gate]
                ax.plot([r['length'] for r in rs],[r['endpoint'] if r['status']=='success' else np.nan for r in rs],color=color,ls='-' if gate=='orientation' else '--',lw=1.3)
                for r in rs:
                    if r['status']=='success':ax.plot(r['length'],r['endpoint'],marker,ms=4,color=color,mfc=color if r['selected_field_tail_pass'] else 'white')
                    elif r['status']=='censored':ax.plot(r['length'],r['maximum_tested_duration'],'^',ms=4.5,mfc='white',mec=color)
                    elif gate=='orientation':ax.plot(r['length'],-.13,'x' if r['status']=='ineligible' else '+',transform=ax.get_xaxis_transform(),clip_on=False,color='#999999',ms=5)
            x=[r['length'] for r in chosen];delta=max(x)-min(x);ax.set_xlim(min(x)-.03*delta,max(x)+.03*delta);ax.set_yscale('log')
            cap=maximum if any(r['status']=='censored' for r in chosen) else 0
            ax.set_ylim(min(valid)/1.5 if valid else maximum/100,max(max(valid) if valid else maximum,cap)*1.6)
            if i==0:ax.set_title(LABELS[j]+ ('\nSpring: native width 0.2' if j==0 else ''),fontsize=10)
            if j==0:ax.set_ylabel(f.LABELS[m]+'\nShortest passing tested cue duration',fontsize=9)
            ax.set_xlabel('Native resting length L0' if m=='spring' else 'Nominal cell length (native units)',fontsize=8,labelpad=20);ax.tick_params(labelsize=7);ax.grid(axis='y',color='#eeeeee',lw=.6)
            ax.text(.98,.035,f'{sum(r["status"]=="ineligible" for r in chosen if r["definition"]=="orientation")} ineligible',transform=ax.transAxes,ha='right',fontsize=6.5,color='#666666')
    fig.suptitle('Supplement 1 · Length versus finite cue duration',fontsize=15,y=.98)
    fig.text(.5,.955,'Five models × three tests · exact shared relative-20% basal checkpoint · unchanged scalar gates and recovery',ha='center',fontsize=9)
    if manifest.get('debelly_boundary_method')=='independent_endpoints':
        fig.text(.5,.942,'deBelly uses fresh corrected preparations and independent endpoint commands',ha='center',fontsize=8,color='#43525B')
    handles=[Line2D([0],[0],color=c,marker=m,ls='-' if g=='orientation' else '--',label=t) for g,c,m,t in f.GATES]
    handles += [Line2D([0],[0],color='#777777',marker='^',mfc='white',ls='none',label='No passing pulse within tested range'),Line2D([0],[0],color='#999999',marker='x',ls='none',label='Ineligible preparation (below axis)'),Line2D([0],[0],color='#777777',marker='o',mfc='white',ls='none',label='Scalar success; native-field flag fails')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.935),ncol=3,fontsize=8,frameon=False)
    factor=manifest['goryachev_selection']['selected_factor']
    fig.text(.5,.035,f'Both patches use a constant peak across all lengths; Goryachev K={factor:g}, peak={.1*factor:g}. Gradients retain one constant physical slope per model.\n'
             'Each positive duration restarts the original checkpoint. Full strength adds the 1% reversed initial-polarity band to the orientation gate.\n'
             'Open markers retain failing field-tail flags. Censoring marks explicit tested limits; missing endpoints break lines. Native times and scales differ.\n'
             f'Data generated {stamp}; producer finished. Every imported and fresh trial passed independent native-array audit.',ha='center',fontsize=7.5,linespacing=1.6)
    fig.subplots_adjust(left=.09,right=.985,top=.865,bottom=.135,hspace=.5,wspace=.22)
    f.save(fig,output,'supplement_width_response',dpi=180)


def render(campaign):
    results=campaign/'results';status=json.loads((results/'status.json').read_text());v=json.loads((results/'validation.json').read_text())
    assert status['complete'] and v['passed'] and v['independent_audit_passed']
    assert f.digest(results/'independent_audit.json')==v['independent_audit_sha256']
    rows=f.read_endpoint_rows(results/'endpoints.csv');assert len(rows)==660
    manifest=json.loads((campaign/'manifest.json').read_text());fits=f.fit_summary(rows)
    check=regression_check(rows,fits);(results/'fit_validation.json').write_text(json.dumps(check,indent=2)+'\n')
    output=ROOT/'figures';output.mkdir(exist_ok=True)
    summary=ROOT/'data/derived/corrected_length_trends.csv'
    with summary.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(fits[0]));writer.writeheader();writer.writerows(fits)
    correlations=correlation_summary(rows)
    correlation_validation=correlation_check(rows,correlations)
    correlation_path=ROOT/'data/derived/length_correlations.csv'
    with correlation_path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(correlations[0]));writer.writeheader();writer.writerows(correlations)
    correlation_validation_path=ROOT/'data/derived/correlation_validation.json'
    correlation_validation_path.write_text(json.dumps(correlation_validation,indent=2)+'\n')
    prep_path=main_figure(rows,correlations,manifest,output);supplement(rows,output,status['generated_at'],manifest)
    provenance=dict(generated_at=f.datetime.now(f.timezone.utc).isoformat(),producer='finished',
            status='three_test_shared_checkpoint_complete',source_generated_at=status['generated_at'],
            campaign=str(campaign.relative_to(ROOT)),debelly_boundary_method=manifest.get('debelly_boundary_method','deprecated_centered'),input_hashes={str(p.relative_to(ROOT)):f.digest(p) for p in [results/'endpoints.csv',results/'validation.json',results/'independent_audit.json',campaign/'manifest.json',prep_path]},
            renderer_sha256=f.digest(Path(__file__)),summary_sha256=f.digest(summary),
            main_summary='data/derived/length_correlations.csv',
            main_summary_sha256=f.digest(correlation_path),main_statistic='Spearman rank correlation',
            correlation_inclusion_rule=CORRELATION_RULE,
            correlation_validation_sha256=f.digest(correlation_validation_path),
            output_hashes={p.name:f.digest(p) for p in output.glob('*') if p.suffix in ['.png','.pdf','.svg'] and p.stem in ['main_polarity','supplement_width_response']})
    (ROOT/'data/derived/figure_generation_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (output/'CAPTIONS.md').write_text('''# Main polarity figure

A. Recorded Wave-pinning preparation at length 15: nominal left cosine patch,
withdrawal and verified basal endpoint. The relative-20% width is 3 here.
All three response tests use the exact nominal relative-20% preparation checkpoint
at each model/length, with no reconditioning or reflection.

B. Two finite right cosine patches have fixed physical width 3 (Spring native
0.2) or relative width 20%. Their peak stays independent of length. Goryachev's
single fixed multiplier is selected by an audited bounded pilot, shared across
both patches and all lengths; every other model retains its nominal peak.
The third challenge is the existing finite g*x countergradient, with constant
physical slope and length-dependent peak. All pulses withdraw to basal input.

C. Spearman rank correlation between nominal length (native resting L0 for
Spring) and shortest passing tested cue duration, separately for each scalar
gate. The same inclusion rule applies to every model and challenge: include
every audited scalar success, irrespective of its separate native-field tail
flag, with at least four distinct successful lengths and nonconstant durations.
Tied durations receive average ranks. No significance or effect-size cutoff is
applied; weak, zero and negative correlations receive colors on the same fixed
-1 to +1 scale. Constant durations have undefined correlation; fewer than four
successful lengths is unavailable. Neither is assigned an artificial zero.
All 30 current cells are available. In-cell n counts the included scalar
successes. Native-field flags remain separate, visible as open markers in
Supplement 1 and as pass/fail/unknown counts in data/derived/length_correlations.csv.
Censored and preparation-ineligible lengths are excluded, never assigned a
duration. These success-conditional summaries can omit difficult lengths and
summarize monotonic association rather than curve shape or a scaling law.
Duration excludes preparation and the 10,000-unit recovery. The 1% strength
gate compares scalar polarity, not full-field convergence. Coefficients are
independently checked using pairwise average ranks and a centered dot product.
The older field-stable-only OLS table is retained in corrected_length_trends.csv
for provenance and is not displayed in panel C.
Statistical definition: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html

# Supplement 1

Five model rows by exactly three response-test columns. Orientation requires all
201 final-500 basal samples to exceed 0.05, with contrast at least 0.05 and
orientation range below 0.005. Strength additionally requires polarity within
1% of reversed initial scalar magnitude. Open success markers fail the separate
0.1% native-field range/projected-drift flag. Triangles mark explicit censor
limits; ineligible preparation is separate. Each endpoint is the shortest tested
passing pulse; every sampled fail-to-pass boundary was refined within 1%.
Duration excludes preparation and 10,000-unit recovery. Exact unchanged trials
can be imported with their immutable native arrays; altered stimuli or starting
checkpoints are integrated afresh. All arrays and imported identities are audited.
''')
    if manifest.get('debelly_boundary_method')=='independent_endpoints':
        with (output/'CAPTIONS.md').open('a') as stream:
            stream.write('\n## Corrected deBelly boundary method\n\nThe deBelly row replaces all centered-model results with fresh relative-20% establishment/release and reversal simulations using independent local Rac/tension endpoint commands (center_boundary_motion=False). Reflection symmetry is retained; inactive endpoints remain prescribed constraints, not stress-free migration. The unchanged four models retain their audited native trials. Old centered-model figures are deprecated and archived.\n')
    import build
    build.index()
    print(json.dumps(dict(output=str(output),fits=check,correlations=correlation_validation,status=provenance['status'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,default=ROOT/'data/runs/debelly_corrected_three_test_20261008')
    render(parser.parse_args().campaign.resolve())
