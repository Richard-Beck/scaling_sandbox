"""Assemble one vector PDF explaining the three fixed-slope duration curves."""
from pathlib import Path
import csv
import json
import textwrap
import subprocess
import hashlib
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/runs/fixed_slope_mechanisms_20261006'
FIG=ROOT/'figures/fixed_slope_mechanisms'
CAMPAIGN=ROOT/'data/runs/fixed_slope_reversal_20261006'
MODELS=['wave_pinning','debelly','holmes_model3']
LABELS=['Wave-pinning','deBelly','Holmes 3']
LENGTHS=[[5.,13.,39.],[18.,23.,39.],[10.,28.,50.]]
BLUE='#0072B2';ORANGE='#D55E00';GRAY='#59636d';INK='#202c38'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})


def page(title,subtitle=''):
    fig=plt.figure(figsize=(11,7),facecolor='white')
    fig.text(.055,.952,title,fontsize=18,weight='bold',color=INK,va='top')
    if subtitle:fig.text(.055,.903,subtitle,fontsize=10,color=GRAY,va='top')
    fig.text(.945,.025,'6 October 2026 | completed analysis',fontsize=7,color=GRAY,ha='right')
    return fig


def prose(fig,x,y,text,width=70,size=10,color=INK):
    wrapped='\n'.join(textwrap.fill(p,width=width) if p else '' for p in text.split('\n'))
    fig.text(x,y,wrapped,fontsize=size,color=color,va='top',linespacing=1.45)


def intro():
    rows=list(csv.DictReader((CAMPAIGN/'results/endpoints.csv').open()))
    manifest=json.loads((CAMPAIGN/'manifest.json').read_text())
    pdf=FIG/'report_introduction.pdf'
    with PdfPages(pdf) as output:
        fig=page('Why does the required gradient pulse first fall, then rise?',
                 'Three native models | one basal preparation | exact passing and neighboring failing trials')
        for i,(model,label,lengths) in enumerate(zip(MODELS,LABELS,LENGTHS)):
            ax=fig.add_axes([.065+i*.315,.455,.275,.335])
            for gate,color,marker,labelgate in [('orientation',BLUE,'o','Rightward-bias criterion'),('full_polarity',ORANGE,'s','Initial scalar strength recovered')]:
                rs=[r for r in rows if r['model']==model and r['condition']=='relative_20pct' and r['definition']==gate and r['status']=='success']
                ax.plot([float(r['length']) for r in rs],[float(r['endpoint']) for r in rs],color=color,marker=marker,ms=3,lw=1.3,label=labelgate)
                for r in rs:
                    if r['selected_field_tail_pass']=='False':ax.plot(float(r['length']),float(r['endpoint']),marker,mfc='white',mec=color,ms=4)
            for length in lengths:
                rs=[r for r in rows if r['model']==model and r['condition']=='relative_20pct' and r['definition']=='orientation' and float(r['length'])==length]
                y=float(rs[0]['endpoint']);ax.scatter(length,y,s=70,facecolors='none',edgecolors=INK,lw=1.3,zorder=5)
            ax.set_title(label,fontsize=12,weight='bold');ax.set_xlabel('Cell length (native units)');ax.set_ylabel('Minimum tested pulse duration');ax.grid(axis='y',alpha=.18)
            ax.tick_params(labelsize=8)
            fig.text(.065+i*.315,.365,'Selected lengths: '+', '.join(f'{L:g}' for L in lengths),fontsize=10,weight='bold',color=INK)
            fig.text(.065+i*.315,.335,f'Fixed physical slope: {manifest["slope_by_model"][model]:.6g}',fontsize=9,color=GRAY)
        handles,labels=fig.axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.87),ncol=2,frameon=False,fontsize=9)
        prose(fig,.065,.268,'The curves measure how long the input must stay on to meet a later recovery test. They do not directly measure how quickly a new right front appears. The spatial evidence shows why this distinction matters.',width=126,size=11)
        prose(fig,.065,.156,'Circled points define the snapshot triplets: descending arm, ordinary-criterion minimum, ascending arm. Holmes has a tied ordinary minimum at lengths 28 and 30; its stricter strength minimum is at length 24, so the same selected cell is not the minimum of both criteria.',width=144,size=9)
        fig.text(.065,.055,'Open markers: separate native-field tail flag fails. Native time units differ by model. Existing complete campaign: 74,295 trials; independent audit passed.',fontsize=8,color=GRAY)
        output.savefig(fig);plt.close(fig)

        fig=page('What was held fixed, and what counts as reversal?',
                 'The input rule is shared. The reported time includes neither preparation nor basal recovery.')
        ax=fig.add_axes([.055,.64,.89,.19]);ax.set_axis_off()
        titles=['Saved left-polarized state','Finite rightward gradient','Input off: native recovery','Final 500-unit readout']
        bodies=['One exact relative-20%\npreparation checkpoint','Added input at position x:\nslope multiplied by x','10,000 native time units\nunder the original model','Rightward bias, contrast,\nstability; then strength']
        for i,(title,body) in enumerate(zip(titles,bodies)):
            x=i*.255
            ax.add_patch(FancyBboxPatch((x,.16),.22,.74,boxstyle='round,pad=.012',edgecolor=BLUE,facecolor='#f4f8fb',lw=1.2))
            ax.text(x+.11,.79,title,ha='center',va='top',fontsize=8.3,weight='bold')
            ax.text(x+.11,.47,body,ha='center',va='center',fontsize=8.2,linespacing=1.5)
            if i<3:ax.annotate('',xy=(x+.244,.53),xytext=(x+.222,.53),arrowprops=dict(arrowstyle='->',color=GRAY))
        ax.set_xlim(-.02,1.02)
        prose(fig,.06,.59,'A fixed slope means a longer cell receives a higher right-end input and more total stimulation. That makes local right-front activation easier. It does not guarantee that the old front, shared inactive pools, or cortical memory have been replaced.',width=140,size=10)
        prose(fig,.06,.47,'Orientation P is activity-weighted position from -1 (left) to +1 (right). Ordinary criterion: all 201 final-window samples have P > 0.05, contrast >= 0.05, and P range < 0.005. Strength adds P within 1% of initial left magnitude. The separate field flag requires pointwise temporal range and projected 500-unit drift <= 0.1% of frozen preparation scales.',width=145,size=9.5)
        audit=json.loads((DATA/'preparation_audit/audit.json').read_text())['models']
        ax=fig.add_axes([.06,.22,.88,.15]);ax.set_axis_off()
        cell=[]
        for model,label in zip(MODELS,LABELS):
            r=audit[model];cell.append([label,str(r['eligible_lengths']),f'{r["maximum_scaled_preparation_difference"]:.2g}','Exactly identical'])
        table=ax.table(cellText=cell,colLabels=['Model','Eligible lengths','Largest field difference / scale','Both duration thresholds across 3 starts'],loc='center',cellLoc='center',colWidths=[.18,.15,.28,.39])
        table.auto_set_font_size(False);table.set_fontsize(8.5);table.scale(1,1.65)
        for (row,col),c in table.get_celld().items():
            c.set_edgecolor('#dce3e8')
            if row==0:c.set_facecolor('#eaf1f7');c.set_text_props(weight='bold')
        prose(fig,.06,.157,'The three preparations add no different duration outcomes at these 55 model-length pairs. Wave-pinning and Holmes saved states are numerically indistinguishable; deBelly differences remain below the 0.0001 preparation tolerance. One preparation is sufficient for this investigation.',width=143,size=9.2)
        fig.text(.06,.055,'All replays start from original checkpoints; no reflection or reconditioning. Longer recovery and altered transport/MCA controls are labeled diagnostic experiments.',fontsize=8,color=GRAY)
        output.savefig(fig);plt.close(fig)
    return pdf


BRIEFS={
    'wave_pinning':dict(
        title='Wave-pinning: start a new cap, then resolve competing caps',
        architecture='Active protein reinforces its own activation. Active and inactive forms share one conserved total; activity is localized, while the inactive form supplies distant regions by diffusion.',
        cards=[
            ('Descending arm | L = 5','Pulse 30.75; neighbor 30.5 fails.','The ramp is weakest here. A sufficient pulse leaves a viable right cap; the slightly shorter pulse loses polarity and becomes homogeneous.'),
            ('Minimum | L = 13','Pulse 13.875; neighbor 13.75 fails.','The stronger ramp starts right activity sooner. Two end caps compete after withdrawal; neighboring pulses select opposite surviving ends.'),
            ('Ascending arm | L = 39','Ordinary 39; strength 49.5.','Right activation is fastest, but broad activation leaves an interior competitor and an unequal inactive supply. Eliminating that competitor takes additional recovery.')],
        test='Two interventions separate the arms. At L=5, matching the L=13 ramp peak by raising slope 2.6-fold makes a 13.875-unit pulse pass both gates. At L=39, increasing only inactive diffusion after the exact native cutoff removes the competitor rapidly. The unforced left state remains polarized, with a 0.668% scalar-strength shift.',
        scope='The rising arm includes an observation-window effect. At L=39, the 38.75 ordinary failure becomes right-biased enough by recovery 50,000. The 49.25 strength failure clears its interior competitor during the official final window, so it fails stability despite eventually reaching a clean right cap.',
        guides=[
            'Spatial snapshots: separate the cue and recovery clocks; follow the new right cap and old interior remnant.',
            'Neighboring pulses: homogeneous loss at L=5, old-cap return at L=13, and a slow scalar crossing at L=39.',
            'Mechanism controls: matched short-cell ramp peak and identical-cutoff inactive mobility test the two bottlenecks.',
            'Full-strength cleanup: pulse 49.25 clears its competitor too late. A 401-cell check preserves the sequence, with changed timing.']),
    'debelly':dict(
        title='deBelly: turn Rac quickly, rewrite cortex memory slowly',
        architecture='MCA means membrane-cortex attachment. Low bound MCA promotes Rac; Rho gradients drive cortex flow, which transports MCA. Attachment and tension feed back to Rac/Rho, storing a coupled chemical and mechanical memory.',
        cards=[
            ('Descending arm | L = 18','Pulse 737.35; neighbor 732.30 fails.','The initial orientation is only -0.05010. Both neighbors ultimately recover right; the failure stays just below +0.05 at the official horizon. Weak polarity makes the gate demanding.'),
            ('Minimum | L = 23','Ordinary 172.90; strength 180.47.','The ordinary failing neighbor is already right-polarized, but still changes by 0.00771 in the final window, above the 0.005 limit. The sampled minimum includes a relaxation requirement.'),
            ('Ascending arm | L = 39','Ordinary 376.20; strength 383.78.','Rac turns right quickly under the cue. The old right-biased attachment imprint is larger; neighboring pulses then recover opposite Rac polarities, even over extended recovery.')],
        test='Freezing all MCA changes only during the cue retains a strong right Rac response, yet native recovery returns left at L=23 and 39. Reducing only the initial odd MCA imprint at L=39 lets pulse 300 succeed where the original state fails; pulses 173 and 250 still fail. MCA change is necessary in the tested examples, and the larger imprint contributes, but its first moment alone is insufficient.',
        scope='The native MCA moment decreases by about 0.085 by cue time 170 at both L=23 and 39; L=39 starts at 0.1798 versus 0.1225. The long cell has more old imprint to undo, without evidence of a slower early rewriting rate. The spatial attachment profile, friction, flow and other original fields matter together.',
        guides=[
            'On-cue snapshots: Rac responds within 5–7.5 units; follow the slower attachment and mechanical changes.',
            'Basal snapshots: fields still evolve at successful scalar readouts. Extended recovery is not convergence certification.',
            'Mechanism tests: amplitude gate, range gate and fate boundary differ. Homogeneous-rest growth gets faster at longer length.']),
    'holmes_model3':dict(
        title='Holmes 3: local switching outruns global pool competition',
        architecture='Rac and Rho inhibit one another. Each has a conserved active-plus-inactive pool; slowly diffusing active domains draw on faster inactive supply. Local Rac switching and removal of the old Rac domain are separate processes.',
        cards=[
            ('Descending arm | L = 10','Both criteria: pulse 8.0625.','The right tip reaches the old Rac plateau at 4.76. Adjacent pulses select a single left or right domain after recovery.'),
            ('Ordinary minimum | L = 28','Ordinary 3.14375; strength 3.67685.','Right seeding takes only 1.16, but the ordinary passing state keeps two Rac-rich domains. Slow elimination of the interior domain makes the recovery horizon matter.'),
            ('Ascending arm | L = 50','Ordinary 7.25; strength 10.125.','Right seeding is fastest: 0.54. The displaced incumbent survives inside the cell; the ordinary two-domain pattern persists to recovery 100,000.')],
        test='Increasing both inactive diffusivities from 50 to 5,000, with the same original starts, cue and reactions, changes the tested strength durations to 7.8125, 2.4375 and 1.859375. The long-length rise disappears. Zero-pulse controls stay left with less than 1% scaled state excursion, supporting a role for pool communication.',
        scope='Ordinary reversal at L=28 and 50 means weak right bias, P about 0.0505 and 0.0510 versus initial magnitudes about 0.326 and 0.327. It can accept two Rac domains. The faster-transport diagnostic changes physical dynamics and perturbs equilibrium slightly; it does not isolate Rac versus Rho or prove a universal diffusion-time law.',
        guides=[
            'Spatial snapshots: track Rac/Rho through cue and recovery; identify the new right domain and displaced old domain.',
            'Pool mechanism: faster tip switching, larger inactive gradients, mixing rescue, and 401-cell domain checks.',
            'Withdrawal fates: L=28 has slow eventual reversals; L=50 retains weakly biased two-domain patterns to 100,000.'])
}


def model_brief(model,label):
    b=BRIEFS[model];pdf=FIG/model/'report_interpretation.pdf'
    fig=page(b['title'],'Physical interpretation and a guide to the following spatial evidence')
    prose(fig,.06,.825,b['architecture'],width=119,size=10.5)
    for i,(title,dose,body) in enumerate(b['cards']):
        x=.055+i*.30
        ax=fig.add_axes([x,.455,.28,.245]);ax.set_axis_off()
        ax.add_patch(FancyBboxPatch((.01,.02),.98,.96,boxstyle='round,pad=.005',edgecolor='#bfcdd8',facecolor='#f4f8fb'))
        ax.text(.045,.90,title,fontsize=10,weight='bold',color=BLUE,va='top')
        ax.text(.045,.76,dose,fontsize=8.5,color=GRAY,va='top')
        ax.text(.045,.60,textwrap.fill(body,43),fontsize=9.1,color=INK,va='top',linespacing=1.4)
    fig.text(.06,.425,'Tests of the explanation',fontsize=10,weight='bold',color=BLUE)
    prose(fig,.06,.395,b['test'],width=146,size=9)
    prose(fig,.06,.285,b['scope'],width=150,size=8.7)
    fig.text(.06,.182,'Read the following figure pages in order',fontsize=9.5,weight='bold',color=BLUE)
    for i,guide in enumerate(b['guides']):
        prose(fig,.075,.153-i*.028,f'{i+1}. {guide}',width=159,size=8.0)
    fig.savefig(pdf);plt.close(fig)
    return pdf


def synthesis():
    pdf=FIG/'report_synthesis.pdf'
    with PdfPages(pdf) as output:
        fig=page('The same curve shape does not imply the same bottleneck',
                 'Interpret the spatial state at withdrawal, the subsequent recovery, and the acceptance rule together.')
        statements=[
            ('Wave-pinning','Short cell: weak right-end drive makes a self-sustaining right cap harder to establish. Intermediate cell: right activation is easier without creating a persistent interior competitor. Long cell: rapid right seeding coexists with a displaced old cap and an unequal inactive pool; slow post-cue competition can dominate the readout.',
             'Direct support: matched-peak short-cell rescue, native snapshots, 50,000-unit recovery, and identical-cutoff inactive-diffusion controls.'),
            ('deBelly','The left arm is strongly affected by the assay: the smallest eligible state lies just above the orientation cutoff, and a right-polarized neighbor can still fail amplitude or tail-range tests. At long length, quick Rac redirection is insufficient; the original MCA/cortical imprint must be weakened enough for native feedback to select the opposite state.',
             'Direct support: extended recovery; cue-only MCA freeze; partial rescue by reducing the original attachment imprint, with insufficient shorter pulses.'),
            ('Holmes 3','Higher right-end input makes the new Rac-rich region appear faster. Removing the old Rac-rich region is a separate conserved-pool competition problem. Larger cells leave an interior remnant; weak scalar reversal can accept that multidomain state. Faster inactive transport removes the rising strength-duration arm.',
             'Direct support: four-field snapshots, adjacent failures, long native recovery, and faster-transport tests with unchanged reactions and input.')]
        for i,(model,statement,evidence) in enumerate(statements):
            y=.80-i*.22
            fig.text(.06,y,model,fontsize=12,weight='bold',color=BLUE)
            prose(fig,.22,y+.004,statement,width=112,size=10)
            prose(fig,.22,y-.118,evidence,width=133,size=8.5,color=GRAY)
        prose(fig,.06,.12,'What remains conditional: these are selected native regimes and a fixed finite recovery window. Transport and MCA interventions change model dynamics; they identify involvement, not a universal scaling law. Scalar-strength recovery is not proof of reflected full-profile identity, and a stable tail is not a separately verified final steady state.',width=145,size=9)
        output.savefig(fig);plt.close(fig)
    return pdf


def main():
    FIG.mkdir(parents=True,exist_ok=True)
    pages=[intro()];assets=[]
    for model in MODELS:
        path=FIG/model/'assets.json'
        if not path.exists():path=FIG/model/'ordered_assets.json'
        payload=json.loads(path.read_text())
        items=payload.get('pages',payload.get('ordered_pages')) if isinstance(payload,dict) else payload
        pages.append(model_brief(model,LABELS[MODELS.index(model)]))
        for item in items:
            filename=item['pdf'] if isinstance(item,dict) else item
            p=Path(filename)
            if not p.is_absolute():p=path.parent/p
            if not p.exists():raise FileNotFoundError(p)
            pages.append(p)
            assets.append(dict(model=model,pdf=str(p.relative_to(ROOT)),
                               caption=item.get('caption',''),claim=item.get('main_claim',item.get('claim',''))))
    pages.append(synthesis())
    final=FIG/'length_duration_mechanisms.pdf'
    subprocess.run(['pdfunite',*[str(p) for p in pages],str(final)],check=True)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',models=MODELS,preparation='relative_20pct',
                  selected_lengths=dict(zip(MODELS,LENGTHS)),report_sha256=sha(final),builder_sha256=sha(Path(__file__)),
                  component_sha256={str(p.relative_to(ROOT)):sha(p) for p in pages},source_campaign_manifest_sha256=sha(CAMPAIGN/'manifest.json'),
                  preparation_audit_sha256=sha(DATA/'preparation_audit/audit.json'),native_arrays='See each model subfolder provenance and replay validations')
    manifest['figures']=assets
    manifest['model_provenance_sha256']={}
    for model in MODELS:
        p=DATA/model/'final_provenance.json'
        if not p.exists():p=DATA/model/'provenance.json'
        manifest['model_provenance_sha256'][str(p.relative_to(ROOT))]=sha(p)
    (DATA/'report_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(final)


if __name__=='__main__':main()
