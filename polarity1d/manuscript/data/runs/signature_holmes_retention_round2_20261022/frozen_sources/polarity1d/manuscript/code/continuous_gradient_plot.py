"""A single PDF for continuous fixed-slope strength arrival across five models."""
from pathlib import Path
import csv,json,sys,textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
CODE=Path(__file__).resolve().parent
sys.path.insert(0,str(CODE))
import continuous_gradient as c
n=c.n
LABELS={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring'}
BLUE='#0072B2';GRAY='#6f7780';INK='#202c38';ORANGE='#D55E00'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})


def header(fig,title,subtitle,stamp):
    fig.suptitle(title,fontsize=17,weight='bold',x=.055,ha='left',y=.97,color=INK)
    fig.text(.055,.91,subtitle,fontsize=9,color=GRAY)
    fig.text(.055,.025,f'Data generated {stamp} | producer finished | independent native-array audit passed',fontsize=7,color=GRAY)


def build(out):
    audit=c.read(out/'results/independent_audit.json')
    if not audit['passed'] or not audit['complete']:raise RuntimeError('Final native-array audit required')
    manifest=c.verify(out);records=c.read(out/'results/records.json');rows=records['cases'];stamp=records['generated_at']
    figure=n.ROOT/'figures/continuous_fixed_slope';figure.mkdir(parents=True,exist_ok=True)
    previous=list(csv.DictReader((c.SOURCE/'results/endpoints.csv').open()))
    final=figure/'continuous_gradient_report.pdf'
    with PdfPages(final) as pdf:
        fig,axs=plt.subplots(2,3,figsize=(11.7,8.3));fig.subplots_adjust(left=.075,right=.97,bottom=.14,top=.80,wspace=.29,hspace=.46)
        header(fig,'Continual gradient: time to regain the original polarity magnitude',
               'Same fixed physical slopes and exact basal checkpoints; the input stays on throughout.',stamp)
        for ax,key in zip(axs.flat,manifest['models']):
            rs=sorted([r for r in rows if r['model']==key],key=lambda r:r['length'])
            success=[r for r in rs if r['status']=='success'];fail=[r for r in rs if r['status']=='not_reached_within_horizon']
            lengths=[r['length'] for r in rs]
            values=[r.get('duration') if r['status']=='success' else np.nan for r in rs]
            ax.plot(lengths,values,'o-',color=BLUE,ms=4,lw=1.5)
            for r in success:ax.vlines(r['length'],r['lower_time'],r['upper_time'],color=BLUE,lw=2)
            if fail:ax.scatter([r['length'] for r in fail],[20000]*len(fail),marker='^',facecolors='none',edgecolors=ORANGE,s=27)
            old=sorted([r for r in previous if r['model']==key and r['condition']=='relative_20pct' and r['definition']=='full_polarity' and r['status']=='success'],key=lambda r:float(r['length']))
            ax.plot([float(r['length']) for r in old],[float(r['endpoint']) for r in old],'--',color=GRAY,lw=1.2)
            ax.set_yscale('log');ax.set_xlabel('Native cell length');ax.set_ylabel('Input duration (native time)')
            ax.set_title(f'{LABELS[key]} | g = {manifest["slope_by_model"][key]:.6g}',fontsize=11,weight='bold')
            ax.grid(axis='y',alpha=.18)
            if fail:ax.set_ylim(top=45000)
            ax.text(.02,.96,f'{len(success)} / {len(success)+len(fail)} eligible reach target',transform=ax.transAxes,fontsize=8,va='top',color=GRAY)
        ax=axs.flat[-1];ax.axis('off')
        lines=[Line2D([0],[0],color=BLUE,marker='o',label='Continual input: first target crossing'),
               Line2D([0],[0],color=GRAY,ls='--',label='Earlier finite pulse: strength gate'),
               Line2D([0],[0],color=ORANGE,marker='^',mfc='none',ls='none',label='Target not reached by 20,000')]
        ax.legend(handles=lines,loc='upper left',frameon=False,fontsize=8.5)
        text='Target: rightward P at least the magnitude of the starting leftward P. Blue times are first passing native samples; bars show crossing brackets.\n\nOne relative-20% preparation per model and length. Unpolarized preparations are excluded before integration.\n\nThe dashed comparison used withdrawal and a later stable, 1%-strength test. It measures a different outcome.'
        ax.text(.03,.66,'\n'.join(textwrap.fill(p,45) for p in text.split('\n')),va='top',fontsize=8.5,linespacing=1.4)
        fig.text(.075,.073,'No withdrawal, confirmation window or amplitude retuning. Native time units differ between models; censoring means only that the target was not reached within this horizon.',fontsize=8)
        pdf.savefig(fig);fig.savefig(figure/'continuous_gradient_duration.png',dpi=180);plt.close(fig)

        fig,axs=plt.subplots(2,3,figsize=(11.7,8.3));fig.subplots_adjust(left=.075,right=.97,bottom=.14,top=.80,wspace=.30,hspace=.43)
        header(fig,'What the polarity score does while the gradient remains on',
               'All eligible lengths; score divided by the original left-polarity magnitude. The target is +1.',stamp)
        for ax,key in zip(axs.flat,manifest['models']):
            rs=sorted([r for r in rows if r['model']==key and r['status']!='preparation_ineligible'],key=lambda r:r['length'])
            norm=plt.Normalize(min(r['length'] for r in rs),max(r['length'] for r in rs));colors=plt.get_cmap('viridis')
            for r in rs:
                with np.load(out/'work/cases'/f'case_{r["case_index"]:03d}'/'trajectory.npz') as z:
                    t=z['time'];o=z['orientation']/r['target_orientation'];ax.plot(t[1:],o[1:],color=colors(norm(r['length'])),lw=.8,alpha=.85)
                    if r['status']=='success':ax.plot(r['duration'],r['orientation_at_crossing']/r['target_orientation'],'o',color=colors(norm(r['length'])),ms=3)
            ax.axhline(1,color=ORANGE,ls='--',lw=1.2);ax.axhline(0,color=GRAY,lw=.5)
            ax.set_xscale('log');ax.set_xlabel('Time under constant gradient');ax.set_ylabel('P / |initial P|');ax.set_title(LABELS[key],fontsize=11,weight='bold');ax.grid(alpha=.12)
            cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=colors),ax=ax,fraction=.046,pad=.03);cb.set_label('Length',fontsize=8)
        ax=axs.flat[-1];ax.axis('off')
        ax.text(.02,.96,'Reading these histories',fontsize=12,weight='bold',va='top',color=INK)
        notes='Each curve begins at -1. Zero is the sign change; +1 is the requested recovered magnitude.\n\nA dot marks the first passing sample. A trajectory is retained through the segment containing that crossing. Cases without a crossing run through 20,000.\n\nA positive plateau below +1 is rightward bias that does not meet this assay. Continued stimulation can change profile shape and polarity magnitude.\n\nThis readout does not establish persistence after withdrawal or a complete reflected steady profile.'
        ax.text(.02,.84,'\n'.join(textwrap.fill(p,46) for p in notes.split('\n')),fontsize=9,va='top',linespacing=1.4)
        pdf.savefig(fig);plt.close(fig)

        fig=plt.figure(figsize=(11.7,8.3));header(fig,'Protocol, results and numerical resolution',
                 'One shared outcome rule for all five models; original native physics and cue coordinates retained.',stamp)
        total_success=sum(r['status']=='success' for r in rows);total_fail=sum(r['status']=='not_reached_within_horizon' for r in rows)
        fig.text(.07,.82,f'97 eligible integrations completed: {total_success} reached the target; {total_fail} did not reach it by 20,000.\n13 preparation-ineligible lengths retained as exclusions; no missing cases or execution errors.',fontsize=11,linespacing=1.6)
        ax=fig.add_axes([.07,.53,.86,.22]);ax.axis('off')
        table_rows=[]
        for key in manifest['models']:
            rs=[r for r in rows if r['model']==key];succ=[r for r in rs if r['status']=='success'];fail=[r for r in rs if r['status']=='not_reached_within_horizon'];bad=[r for r in rs if r['status']=='preparation_ineligible']
            table_rows.append([LABELS[key],len(succ),len(fail),len(bad),f'{min(r["duration"] for r in succ):.6g}–{max(r["duration"] for r in succ):.6g}' if succ else 'Not reached'])
        table=ax.table(cellText=table_rows,colLabels=['Model','Reached target','Not reached by cap','Preparation ineligible','Passing duration range'],loc='center',cellLoc='center',colWidths=[.20,.17,.20,.20,.23]);table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.6)
        for (row,col),cell in table.get_celld().items():
            cell.set_edgecolor('#dce3e8')
            if row==0:cell.set_facecolor('#eaf1f7');cell.set_text_props(weight='bold')
        paragraphs=[
            'Starting state: exact qualifying saved relative-20% basal-release checkpoint, without gradient conditioning, reflection or re-equilibration. The same selected model slope from the previous finite-gradient experiment stays on from time zero.',
            'Outcome: first recorded rightward normalized active/Rac first moment P >= |P initially|. There is no 1% tolerance band or additional contrast/stability gate. Native contrast is saved at the crossing; this scalar measure does not require complete spatial-profile replacement.',
            'Sampling: geometrically increasing integration segments retain closely spaced native profiles. deBelly output times respect its original 0.5-unit staggered step. Every success stores its last failing/first passing bracket and an interpolation estimate separately; only a recorded passing value determines success.',
            'Audit: independently reconstructed first moments from every native activity profile, checked exact initial fields/checkpoint hashes, first crossings, complete coverage, conserved material and unchanged constant-gradient construction. Scientific dependencies and input manifests were frozen before submission. All integrations ran on SLURM QoS small, one CPU each.'
        ]
        for i,text in enumerate(paragraphs):fig.text(.07,.475-i*.095,textwrap.fill(text,143),fontsize=9,va='top',linespacing=1.4)
        pdf.savefig(fig);plt.close(fig)
    n.write(out/'results/figure_provenance.json',dict(generated_at=n.now(),producer='finished',report=str(final.relative_to(n.ROOT)),
            report_sha256=n.sha(final),audit_sha256=n.sha(out/'results/independent_audit.json'),
            endpoints_sha256=n.sha(out/'results/endpoints.csv'),manifest_sha256=n.sha(out/'manifest.json')))
    print(final)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=c.DEFAULT);build(p.parse_args().out.resolve())
