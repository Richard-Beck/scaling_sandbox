"""S10 figures for the predeclared A6 export, with direct/proxy separation."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring='Spring')
COLORS=dict(zip(NAMES,['#3974a0','#c58428','#76518c','#2b8a77']))
PRIMARY=['multiple_front_fraction','region_count_mean']
RAW_TOTAL_COUNTS={'observed_births','observed_takeovers','observed_pairwise_takeovers','unscorable_takeover_tracks'}
LABELS=dict(multiple_front_fraction='Multiple-front fraction',region_count_mean='Mean region count')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def pretty(name):return LABELS.get(name,name.replace('_',' '))
def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8)


def summary_pages(pdf,points,cross,folds,out,context=''):
    fig,axes=plt.subplots(2,2,figsize=(10.5,10.5))
    for i,statistic in enumerate(PRIMARY):
        for j,outcome in enumerate(['direct_tau','Dtheta_proxy']):
            ax=axes[i,j]
            frame=points[(points.statistic==statistic)&(points.imaging_view=='primary')&(points.outcome==outcome)]
            for model in NAMES:
                part=frame[(frame.model==model)&frame.feature_ratio.notna()&frame.log_persistence_ratio.notna()]
                for forced in [False,True]:
                    p=part[part.forced==forced]
                    ax.scatter(p.feature_ratio,p.log_persistence_ratio,facecolor='none' if forced else COLORS[model],
                        edgecolor=COLORS[model],marker='s' if forced else 'o',s=35,linewidth=.9)
                if 'forcing_status' in part:
                    unknown=part[part.forcing_status=='unavailable_unknown']
                    ax.scatter(unknown.feature_ratio,unknown.log_persistence_ratio,color=COLORS[model],marker='x',s=30,linewidth=.9)
            rank=cross[(cross.statistic==statistic)&(cross.imaging_view=='primary')&(cross.outcome==outcome)&~cross.forced_inclusive].iloc[0]
            rho=f'{rank.spearman:.3f}' if pd.notna(rank.spearman) else 'unavailable'
            ax.set_title(f'{pretty(statistic)} · '+('direct τ ratio' if outcome=='direct_tau' else '$D_{\\theta}$ proxy sensitivity')+
                f'\nNonforced Spearman {rho}; signs {int(rank.successful_folds)}/4 ({int(rank.unavailable_folds)} unavailable)',fontsize=9.5,weight='bold')
            ax.set_xlabel('Within-model large / small feature quotient',fontsize=9)
            ax.set_ylabel('log(large / small persistence)',fontsize=9);clean(ax)
    handles=[Line2D([],[],color=COLORS[m],marker='o',ls='',label=NAMES[m]) for m in NAMES]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.95),ncol=4,frameon=False,fontsize=9)
    primary=cross[cross.statistic.isin(PRIMARY)&(cross.imaging_view=='primary')&(cross.outcome=='direct_tau')&~cross.forced_inclusive]
    supported=primary.A6_primary_predictor_support.map(lambda value:str(value).lower()=='true').any()
    fig.suptitle('S10 · Direct-τ A6 prediction criterion '+('is met' if supported else 'is not met'),fontsize=13,y=.98)
    fig.text(.5,.018,context+' Crosses: unknown forcing, excluded from primary.\n'
        'Filled circles: nonforced primary conditions; open squares: forced context. No absolute native quantity is pooled across models.\n'
        'Direct resolved τ ratios and $D_{\\theta}$-derived proxy sensitivity are separate. Marginal whole-seed CIs remain in the source tables.\n'
        'Primary support requires Spearman ≤ −0.6 and negative training/held-out sign in at least 3 of 4 planned folds.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.08,1,.93));pdf.savefig(fig)
    fig.savefig(out/'A6_primary_scatter.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10.5,8.2))
    palette=ListedColormap(['#d7dce0','#e3a16a','#6cb29f'])
    for i,statistic in enumerate(PRIMARY):
        for j,outcome in enumerate(['direct_tau','Dtheta_proxy']):
            ax=axes[i,j];part=folds[(folds.statistic==statistic)&(folds.imaging_view=='primary')&(folds.outcome==outcome)&~folds.forced_inclusive]
            matrix=np.array([[dict(unavailable=0,failure=1,success=2)[part[part.held_out_model==m].iloc[0].status] for m in NAMES]])
            ax.imshow(matrix,cmap=palette,vmin=-.5,vmax=2.5,aspect='auto')
            for column,model in enumerate(NAMES):
                row=part[part.held_out_model==model].iloc[0]
                train=f'{row.training_spearman:.2f}' if pd.notna(row.training_spearman) else 'NA'
                held=f'{row.held_out_spearman:.2f}' if pd.notna(row.held_out_spearman) else 'NA'
                ax.text(column,0,f'{row.status}\ntrain {train}\nheld {held}\nn={int(row.held_out_n)}',ha='center',va='center',fontsize=8)
            ax.set_xticks(range(4),[NAMES[m] for m in NAMES],rotation=20,ha='right',fontsize=8);ax.set_yticks([])
            ax.set_title(pretty(statistic)+'\n'+('Direct τ' if outcome=='direct_tau' else '$D_{\\theta}$ proxy sensitivity'),fontsize=10)
    fig.suptitle('S10 · Every planned held-out-model fold remains in the denominator',fontsize=13)
    fig.text(.5,.025,context+'\nGreen: negative training and held-out signs. Orange: evaluated sign failure. Grey: insufficient/constant/censored data.\n'
        'Each training and held-out sample needs at least 3 matched input fractions and nonconstant predictor/outcome.\n'
        'Absence of a qualifying Goryachev or other model is explicit; no missing fold is counted as successful.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.10,1,.94));pdf.savefig(fig);fig.savefig(out/'A6_four_fold_audit.png',dpi=180);plt.close(fig)


def within_pages(pdf,points,within,context=''):
    primary=points[points.imaging_view=='primary']
    for statistic in sorted(primary.statistic.unique()):
        fig,axes=plt.subplots(4,3,figsize=(10.5,13.5))
        for i,model in enumerate(NAMES):
            for j,(transform,column,label) in enumerate([('small','small_estimate','Small-cell feature (native units)'),
                ('large','large_estimate','Large-cell feature (native units)'),('large_over_small','feature_ratio','Within-model large / small feature')]):
                ax=axes[i,j]
                for outcome,marker in [('direct_tau','o'),('Dtheta_proxy','^')]:
                    part=primary[(primary.statistic==statistic)&(primary.model==model)&(primary.outcome==outcome)]
                    if 'observational_pair_status' in part:part=part[part.observational_pair_status=='matched']
                    for forced in [False,True]:
                        p=part[(part.forced==forced)&part[column].notna()&part.log_persistence_ratio.notna()]
                        ax.scatter(p[column],p.log_persistence_ratio,marker=marker,s=20,
                            facecolor='none' if forced else COLORS[model],edgecolor=COLORS[model],linewidth=.8)
                    if 'forcing_status' in part:
                        unknown=part[(part.forcing_status=='unavailable_unknown')&part[column].notna()&part.log_persistence_ratio.notna()]
                        ax.scatter(unknown[column],unknown.log_persistence_ratio,color=COLORS[model],marker='x',s=20,linewidth=.8)
                rank=within[(within.statistic==statistic)&(within.imaging_view=='primary')&(within.model==model)&
                    (within.feature_transform==transform)&(within.outcome=='direct_tau')&~within.forced_inclusive].iloc[0]
                value=f'{rank.spearman:.2f}' if pd.notna(rank.spearman) else 'NA'
                ax.set_title(f'{NAMES[model]} · direct Spearman {value} (n={int(rank.n)})',fontsize=8)
                ax.set_xlabel(label,fontsize=8);ax.set_ylabel('log persistence quotient',fontsize=8);clean(ax)
        fig.suptitle('S10 · '+pretty(statistic),fontsize=13)
        pending=(primary[(primary.statistic==statistic)].feature_status=='pending_strict_global_takeover_replay').any()
        fig.text(.5,.018,('TAKEOVER FEATURE HELD PENDING STRICT GLOBAL REPLAY.\n' if pending else '')+
            ('Raw event totals retain actual exposure; size quotients are unavailable for unequal horizons.\n' if statistic in RAW_TOTAL_COUNTS else '')+
            context+' Crosses: unknown forcing, excluded from primary.\n'+
            '○ Direct resolved τ ratio; △ $D_{\\theta}$ proxy sensitivity. Filled: nonforced; open: forced context.\n'
            'All small/large native quantities remain within their own model. Cross-model tests use only size quotients.\n'
            'NA retains insufficient, constant, censored or unscorable measurements; native caps are never substituted.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.075,1,.955));pdf.savefig(fig);plt.close(fig)


def sensitivity_pages(pdf,cross,context=''):
    for inclusive in [False,True]:
        part=cross[cross.statistic.isin(PRIMARY)&(cross.forced_inclusive==inclusive)]
        views=sorted(part.imaging_view.unique(),key=lambda value:(value!='primary',value))
        rows=[]
        for view in views:
            label='Primary imaging view' if view=='primary' else view.replace('threshold','q=').replace('_smooth0',' / raw').replace('_smooth1',' / optical').replace('_cadence',' / Δt=')+'×'
            values=[label]
            for statistic in PRIMARY:
                for outcome in ['direct_tau','Dtheta_proxy']:
                    matching=part[(part.imaging_view==view)&(part.statistic==statistic)&(part.outcome==outcome)]
                    row=matching.iloc[0] if len(matching) else None
                    rho=f'{row.spearman:.3f}' if row is not None and pd.notna(row.spearman) else 'NA'
                    values.extend([rho,f'{int(row.successful_folds)}/4 ({int(row.unavailable_folds)} NA)' if row is not None else '0/4 (4 NA)'])
            rows.append(values)
        fig,ax=plt.subplots(figsize=(10.5,13.5));ax.axis('off')
        labels=['A5 imaging view','Multi-front\nτ rank','τ\nfolds','Multi-front\nproxy rank','Proxy\nfolds','Count\nτ rank','τ\nfolds','Count\nproxy rank','Proxy\nfolds']
        table=ax.table(cellText=rows,colLabels=labels,colWidths=[.32]+[.085]*8,cellLoc='center',loc='upper center',bbox=[0,.10,1,.84])
        table.auto_set_font_size(False);table.set_fontsize(7)
        for (row,column),cell in table.get_celld().items():
            if row==0:cell.set_facecolor('#e7edf1');cell.set_text_props(weight='bold')
            elif row==1:cell.set_facecolor('#e6f1f4')
        fig.suptitle('S10 · Primary predictor sensitivity · '+('forced included' if inclusive else 'forced excluded'),fontsize=12)
        fig.text(.5,.035,context+'\nRank entries are Spearman coefficients; multi-front denotes the fraction of time with at least two qualified regions.\n'
            'Every threshold, raw/physical-kernel and observation-cadence view is retained. Δt is relative to the declared common imaging cadence.\n'
            'Direct resolved τ and $D_{\\theta}$-proxy associations remain separate. All four folds are planned; NA counts remain in the denominator.\n'
            'Only the primary, forced-excluded, direct-τ rows enter the declared H2 A6 gate. Other rows are sensitivity checks.',ha='center',fontsize=8)
        pdf.savefig(fig);plt.close(fig)


def export_page_data(source,out,tables,gate):
    points,cross,folds,within=[tables[name] for name in ['condition_pairs','cross_model','LOMO_folds','within_model']]
    primary=points[points.imaging_view=='primary'];specs=[]
    tagged=lambda frame,kind:frame.assign(record_kind=kind)
    summary=primary[primary.statistic.isin(PRIMARY)].copy()
    summary['actually_drawn']=summary.feature_ratio.notna()&summary.log_persistence_ratio.notna()
    specs.append(('Primary direct/proxy scatter',pd.concat([tagged(summary,'point'),tagged(cross[cross.statistic.isin(PRIMARY)&(cross.imaging_view=='primary')],'rank')],ignore_index=True,sort=False)))
    specs.append(('Four planned held-out-model folds',tagged(folds[folds.statistic.isin(PRIMARY)&(folds.imaging_view=='primary')&~folds.forced_inclusive],'plotted_fold')))
    for statistic in sorted(primary.statistic.unique()):
        part=primary[primary.statistic==statistic];frames=[]
        for transform,column in [('small','small_estimate'),('large','large_estimate'),('large_over_small','feature_ratio')]:
            frame=part.copy();frame['feature_transform']=transform;frame['x_value']=frame[column]
            matched=frame.observational_pair_status=='matched' if 'observational_pair_status' in frame else np.ones(len(frame),dtype=bool)
            frame['actually_drawn']=matched&frame.x_value.notna()&frame.log_persistence_ratio.notna();frames.append(tagged(frame,'point'))
        ranks=within[(within.statistic==statistic)&(within.imaging_view=='primary')&(within.outcome=='direct_tau')&~within.forced_inclusive]
        specs.append((pretty(statistic),pd.concat([*frames,tagged(ranks,'annotated_rank')],ignore_index=True,sort=False)))
    for inclusive in [False,True]:
        part=cross[cross.statistic.isin(PRIMARY)&(cross.forced_inclusive==inclusive)]
        rows=[]
        for view in sorted(part.imaging_view.unique(),key=lambda value:(value!='primary',value)):
            for statistic in PRIMARY:
                for outcome in ['direct_tau','Dtheta_proxy']:
                    matching=part[(part.imaging_view==view)&(part.statistic==statistic)&(part.outcome==outcome)]
                    rows.append(matching.iloc[0].to_dict() if len(matching) else dict(statistic=statistic,imaging_view=view,outcome=outcome,
                        forced_inclusive=inclusive,spearman=np.nan,status='unavailable',reason='observational_view_absent',
                        successful_folds=0,planned_folds=4,unavailable_folds=4))
        specs.append(('Primary predictor sensitivity, '+('forced included' if inclusive else 'forced excluded'),tagged(pd.DataFrame(rows),'plotted_sensitivity_row')))
    directory=out/'page_data';directory.mkdir();records=[]
    for page,(label,frame) in enumerate(specs,1):
        path=directory/f'S10_page_{page:03d}.csv';frame.to_csv(path,index=False,float_format='%.17g')
        roundtrip=pd.read_csv(path,float_precision='round_trip')
        for field in frame.select_dtypes(include='floating').columns:
            np.testing.assert_array_equal(np.asarray(roundtrip[field],float),np.asarray(frame[field],float))
        records.append(dict(page=page,label=label,csv=str(path),csv_sha256=sha(path),rows=len(frame),float64_roundtrip_equal=True))
    receipt=dict(passed=True,pages=len(records),analysis_arm=gate.get('analysis_arm','E0'),sigma=gate.get('sigma',1.),
        source_complete_sha256=sha(source/'complete.json'),PDF=str(out/'figure4_S10_A6.pdf'),PDF_sha256=sha(out/'figure4_S10_A6.pdf'),
        native_horizon_pairs=gate.get('native_horizon_pairs'),source_renderer_sha256=sha(__file__),page_data=records,
        scope='Every plotted point/transform, annotated rank, planned fold and sensitivity table row; unavailable coverage is retained. Numeric floats roundtrip at17 significant digits.')
    (out/'page_data_index.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    return receipt


def exhaustive_pdf_word_bounds(path,output):
    root=ET.fromstring(subprocess.check_output(['pdftotext','-bbox',str(path),'-']))
    ns='{http://www.w3.org/1999/xhtml}';pages=[];bad=[]
    for number,page in enumerate(root.iter(ns+'page'),1):
        width,height=float(page.attrib['width']),float(page.attrib['height']);words=list(page.iter(ns+'word'));outside=[]
        for word in words:
            x0,y0,x1,y1=[float(word.attrib[key]) for key in ['xMin','yMin','xMax','yMax']]
            if x0<0 or y0<0 or x1>width or y1>height:
                outside.append(dict(page=number,text=word.text,bbox=[x0,y0,x1,y1],page_size=[width,height]))
        bad.extend(outside);pages.append(dict(page=number,words=len(words),out_of_bounds_words=outside))
    result=dict(passed=not bad,PDF_sha256=sha(path),pages_checked=len(pages),
        words_checked=sum(p['words'] for p in pages),strict_tolerance_points=0,
        all_pages=pages,out_of_bounds_words=bad,method='Exhaustive pdftotext -bbox XHTML word bounds against each actual PDF page rectangle.')
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    if bad:raise ValueError('S10 PDF word bounds failed; complete diagnostic retained')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--allow-provisional',action='store_true');args=p.parse_args()
    source,out=args.analysis.resolve(),args.output.resolve();gate=read(source/'complete.json')
    assert gate['complete'] and gate['passed']
    if gate['strict_global_takeover_features_pending'] and not args.allow_provisional:
        raise ValueError('Final S10 figures require the accepted strict-global A5 export')
    for name,expected in gate['outputs'].items():assert sha(source/name)==expected
    tables={name:pd.read_csv(source/(name+'.csv'),float_precision='round_trip') for name in ['condition_pairs','cross_model','LOMO_folds','within_model']}
    if out.exists():raise FileExistsError('S10 rendering uses a fresh immutable output directory')
    out.mkdir(parents=True)
    declared_horizons='native_horizon_small' in tables['condition_pairs'] and tables['condition_pairs'].native_horizon_small.notna().any()
    context=f"Arm {gate.get('analysis_arm','E0')} · σ={gate.get('sigma',1.):g} · "+('actual per-size horizons retained in the source CSV.' if declared_horizons else 'native horizons follow the canonical E0 source inventory.')
    with PdfPages(out/'figure4_S10_A6.pdf') as pdf:
        summary_pages(pdf,tables['condition_pairs'],tables['cross_model'],tables['LOMO_folds'],out,context)
        within_pages(pdf,tables['condition_pairs'],tables['within_model'],context)
        sensitivity_pages(pdf,tables['cross_model'],context)
    bounds=exhaustive_pdf_word_bounds(out/'figure4_S10_A6.pdf',out/'pdf_word_bounds_audit.json')
    shutil.copy2(__file__,out/'frozen_renderer.py')
    pages=4+tables['condition_pairs'][tables['condition_pairs'].imaging_view=='primary'].statistic.nunique()
    page_data=export_page_data(source,out,tables,gate);assert page_data['pages']==pages
    assert bounds['pages_checked']==pages
    provenance=dict(generated_utc=datetime.now(timezone.utc).isoformat(),pages=int(pages),source_complete_sha256=sha(source/'complete.json'),
        analysis_arm=gate.get('analysis_arm','E0'),sigma=gate.get('sigma',1.),native_horizon_pairs=gate.get('native_horizon_pairs'),
        renderer_sha256=sha(Path(__file__)),strict_takeover_features_pending=gate['strict_global_takeover_features_pending'],
        outputs={p.name:sha(p) for p in out.iterdir() if p.suffix in ['.png','.pdf']},
        page_data_index_sha256=sha(out/'page_data_index.json'),all_page_CSV_float64_roundtrip_equal=True,
        exhaustive_PDF_word_bounds_sha256=sha(out/'pdf_word_bounds_audit.json'),
        scope='All main-imaging A6 statistics, direct/proxy and forced markings; complete threshold/cadence/smoothing rank tables remain retained in the immutable source export.')
    (out/'generation.json').write_text(json.dumps(provenance,indent=2)+'\n');print(json.dumps(provenance))


if __name__=='__main__':main()
