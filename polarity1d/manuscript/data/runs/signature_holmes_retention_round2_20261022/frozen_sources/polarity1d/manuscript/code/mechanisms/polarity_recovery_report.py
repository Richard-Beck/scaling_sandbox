"""Polarity and two-pole trajectories from audited native three-test histories."""
from pathlib import Path
import csv, gzip, hashlib, html, json, textwrap
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from three_test_report import Evidence, OUT as SOURCE, GEOMETRIES, NAMES
from initial_polarity_report import amount, LABELS, ROOT

OUT=ROOT/'figures/polarity_recovery'
DATA=ROOT/'data/runs/polarity_recovery_report_20261008'
BLUE='#0072B2'; ORANGE='#D55E00'; PURPLE='#8E5AA7'
STAMP=datetime.now(timezone.utc).isoformat()
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def pole(u,q,L,contrast,spring):
    """Boundary-connected pole, within its own half, above fixed initial contrast."""
    level=float(u.min()+.5*contrast)
    results=[]
    for values in [u,u[::-1]]:
        if values[0]<level:
            results.append(dict(width=0.,peak=np.nan,mean=np.nan,total=0.));continue
        below=np.flatnonzero(values<level)
        if not len(below):fraction=.5
        else:
            i=int(below[0]);assert i>0
            fraction=min(.5,float(q[i-1]+(level-values[i-1])/(values[i]-values[i-1])*(q[i]-q[i-1])))
        width=fraction*L
        total=amount(values,q,width,L,spring)
        # Peak of the same connected region, including interpolated edge.
        inside=q<=fraction
        peak=max(float(values[0]),float(values[inside].max()) if inside.any() else float(values[0]))
        mean=total/width
        assert 0<width<=.5*L+1e-12 and total>=0 and mean<=peak+1e-9
        results.append(dict(width=width,peak=peak,mean=mean,total=total))
    return results

def traces(evidence,record):
    c=record['job']['case'];model=c['model'];D=record['job']['duration']
    with evidence.history(record) as z:
        assert 'pulse_time' in z
        time=np.r_[z['pulse_time'],D+z['recovery_time'][1:]]
        u=np.r_[z['pulse_active'],z['recovery_active'][1:]]
        orient=np.r_[z['pulse_orientation'],z['recovery_orientation'][1:]]
        q=np.linspace(0,1,u.shape[1]) if model=='spring' else (np.arange(u.shape[1])+.5)/u.shape[1]
        w=np.ones(len(q))/len(q)
        if model=='spring':w[:]=1/(len(q)-1);w[[0,-1]]*=.5
        independently=(u*(2*q-1))@w/(u@w)
        np.testing.assert_allclose(independently,orient,rtol=1e-9,atol=1e-11)
        lengths=np.r_[z['pulse_field_length'],z['recovery_field_length'][1:]] if model=='spring' else np.full(len(time),c['length'])
        with np.load(ROOT/c['checkpoint']) as prep:
            assert digest(ROOT/c['checkpoint'])==c['checkpoint_sha256']
            np.testing.assert_allclose(u[0],prep['active'][-1],rtol=0,atol=1e-12)
            if model=='debelly':
                for name in [k for k in prep.files if k.startswith('state_')]:np.testing.assert_array_equal(z['originalstate_'+name[6:]],prep[name])
            else:np.testing.assert_array_equal(z['originalstate'],prep['state'])
        mca=None
        if model=='debelly':
            b=np.r_[z['pulse_field_mca_bound'],z['recovery_field_mca_bound'][1:]]
            mca=(b*(2*q-1))@w/(b@w)
        assert np.all(np.diff(time)>0) and time[0]==0 and np.isclose(time[-1],D+10000)
        result=dict(time=time,active_polarity=orient,physical_length=lengths,mca_polarity=mca)
        for side in ['old','new']:
            for key in ['width','peak','mean','total']:result[f'{side}_{key}']=np.empty(len(time))
        for i,values in enumerate(u):
            old,new=pole(values,q,float(lengths[i]),float(np.ptp(u[0])),model=='spring')
            for side,entry in [('old',old),('new',new)]:
                for key,value in entry.items():result[f'{side}_{key}'][i]=value
        return result

def main():
    OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    E=Evidence(SOURCE);pages=[];chosen=[];selection={}
    models=list(LABELS)
    for model in models:
        anchors=E.anchors(model,'fixed_relative')
        selection[model]=[anchors[i]['job']['case']['length'] for i in [0,len(anchors)//2,-1]]
    csvpath=DATA/'pole_trajectories.csv.gz'
    with PdfPages(OUT/'polarity_recovery.pdf') as pdf, gzip.open(csvpath,'wt',newline='') as stream:
        writer=None
        def save(fig,title,caption,slug):
            fig.suptitle(title,fontsize=20,weight='bold',y=.98)
            lines=textwrap.wrap(caption,165)
            cy=.04+(len(lines)*13.5/72+.12)/fig.get_figheight()
            fig.text(.065,cy,'\n'.join(lines),fontsize=9,va='top',linespacing=1.5)
            fig.text(.065,.015,'Cue and recovery · exact native histories · five models, three cue classes',fontsize=9,color='#526674')
            fig.text(.95,.015,str(len(pages)+1),fontsize=9,ha='right')
            fig.savefig(OUT/f'{slug}.png',dpi=125);pdf.savefig(fig);plt.close(fig)
            pages.append(dict(title=title,caption=caption,image=f'{slug}.png'))
        fig=plt.figure(figsize=(15,13))
        paragraphs=[('How to read the trajectories','Each of the following fifteen pages is one model × cue class. Its three columns show short, middle and long cells. All plots use logarithmic elapsed time since cue onset. Gray shading is cue on; the vertical dashed line is cue withdrawal. Time continues through 10,000 native recovery units; it does not reset at withdrawal. The t=0 state cannot appear on a log axis; its polarity is printed above each column and is retained in the exported data. All subsequent recorded samples are plotted without a time offset.'),
        ('Which trials are used?','Native orientation-pass replays at the shortest tested successful pulse, or the tested maximum pulse when orientation was censored. These exact trials were independently audited in the completed three-test mechanism campaign. Three lengths per model are held the same across all cue classes: '+ '; '.join(LABELS[m]+' '+', '.join(f'{v:g}' for v in selection[m]) for m in models)+'. Durations vary with model, length and cue class. Fixed support is 3 (Spring 0.2), relative support is 20%, and the gradient has constant native physical slope.'),
        ('Polarity and recovery','Active/Rac polarity is its signed normalized spatial first moment: negative = old left pole, positive = new right pole. deBelly adds a separate bound-MCA polarity row using the same definition; positive MCA polarity means opposite-end attachment enrichment. Column titles report the original scalar reversal outcome and the independent field-tail stability flag. A stable scalar outcome does not establish a newly verified full-field steady state. Censored cases remain visible, including unsuccessful long fixed-patch trajectories.'),
        ('A consistent definition of both poles','At each sample, the pole must connect to its own end and exceed the current minimum active concentration plus half the initial active contrast. The contrast increment is frozen at t=0 so fading noise does not become a pole. Each region is confined to its own half of the cell, preventing overlap during broad activation. Its half-contrast edge is interpolated between native samples. Old = left end; new = right end throughout, regardless of polarity sign.'),
        ('Width, peak, mean and total','Blue is the old pole; orange is the new pole. Width is the extent of the defined region in native domain coordinates. Peak and mean are concentrations within that region; total is its integral, including background, in concentration × length units. A missing pole has zero width and total, with peak and mean undefined (a gap, not zero concentration). Finite-volume amounts use cell averages and fractional edge-cell overlap; Spring uses trapezoidal integration over its instantaneous physical length. deBelly chemical amounts use the native nominal domain; no physical cross-section or molecule conversion is assumed.'),
        ('Sampling and provenance','These are the stored histories, with 20 equal cue-time intervals (more for sufficiently long pulses) and recovery samples every 2.5 units. Very rapid recovery can occur between stored samples, particularly for Spring; plotted connections do not supply additional temporal resolution. Archive hashes, starting-checkpoint hashes, exact restart states and initial active arrays (absolute tolerance 1e-12) and all active first moments are checked again by this renderer. Source producer finished '+E.status['generated_at']+'. Report generated '+STAMP+'. Existing native data and outcomes are unchanged.')]
        y=.92
        for heading,body in paragraphs:
            fig.text(.065,y,heading,fontsize=14,weight='bold',va='top');y-=.03
            lines=textwrap.wrap(body,155);fig.text(.065,y,'\n'.join(lines),fontsize=11,va='top',linespacing=1.45)
            y-=len(lines)*.0185+.024
        assert y>.11,y
        save(fig,'Polarity, cue withdrawal and recovery','All 45 selected trajectories are native replays from the current audited comparison. The plots show recorded dynamics and retain the original outcome classifications.','01_reading_guide')
        for model in models:
            for geometry in GEOMETRIES:
                records=[E.anchor(model,geometry,L) for L in selection[model]]
                nrows=6 if model=='debelly' else 5
                fig,axes=plt.subplots(nrows,3,figsize=(15,3*nrows+2),squeeze=False)
                for col,record in enumerate(records):
                    job=record['job'];c=job['case'];D=job['duration'];t=traces(E,record)
                    chosen.append(dict(slot=job['slot'],model=model,cue_class=geometry,length=c['length'],duration=D,archive_sha256=record['archive_sha256'],checkpoint_sha256=c['checkpoint_sha256'],metrics=record['metrics']))
                    mask=t['time']>0;times=t['time'][mask]
                    specs=[('active_polarity','Active / Rac polarity','Normalized first moment')]
                    if model=='debelly':specs.append(('mca_polarity','Bound-MCA polarity','Normalized first moment'))
                    specs += [('width','Pole width','Native length'),('peak','Peak active / Rac','Concentration'),('mean','Mean active / Rac','Concentration'),('total','Total active / Rac material','Concentration × length')]
                    for row,(key,label,unit) in enumerate(specs):
                        ax=axes[row,col]
                        if key.endswith('polarity'):
                            ax.plot(times,t[key][mask],color=PURPLE if key=='mca_polarity' else BLUE,lw=1.7)
                            ax.axhline(0,color='#999999',lw=.7)
                        else:
                            for side,color in [('old',BLUE),('new',ORANGE)]:ax.plot(times,t[f'{side}_{key}'][mask],color=color,lw=1.7,label=side.capitalize()+' pole')
                            ax.set_ylim(bottom=0)
                            if row==len(specs)-4:ax.legend(frameon=False,fontsize=9,loc='best')
                        ax.set_xscale('log');ax.set_xlim(times[0],times[-1]);ax.axvspan(times[0],D,color='#777777',alpha=.11);ax.axvline(D,color='#555555',ls='--',lw=.85)
                        ax.set_ylabel(label+'\n'+unit if col==0 else unit,fontsize=10);ax.grid(alpha=.15)
                        ax.set_xlabel('Time since cue onset (log axis)',fontsize=9)
                        ax.tick_params(labelsize=9)
                        ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useOffset=False)
                    outcome='retained reversal' if record['metrics']['orientation_success'] else 'no retained reversal (censored)'
                    field='field tail passes' if record['metrics']['field_tail_pass'] else 'field tail still changing'
                    p0=f'Active P0 = {t["active_polarity"][0]:.3g}'
                    if model=='debelly':p0+=f'; MCA P0 = {t["mca_polarity"][0]:.3g}'
                    axes[0,col].set_title(f'{"L0" if model=="spring" else "L"} = {c["length"]:g}; cue duration = {D:.4g}\n{outcome}; {field}\n{p0}',fontsize=10.5,pad=12)
                    for i,time in enumerate(t['time']):
                        entry=dict(model=model,cue_class=geometry,domain_length=c['length'],slot=job['slot'],cue_duration=D,time_since_onset=float(time),phase='cue' if time<=D else 'recovery',active_polarity=float(t['active_polarity'][i]),mca_polarity=float(t['mca_polarity'][i]) if t['mca_polarity'] is not None else '',physical_length=float(t['physical_length'][i]))
                        for side in ['old','new']:
                            for key in ['width','peak','mean','total']:
                                value=t[f'{side}_{key}'][i];entry[f'{side}_{key}']=float(value) if np.isfinite(value) else ''
                        if writer is None:writer=csv.DictWriter(stream,fieldnames=list(entry));writer.writeheader()
                        writer.writerow(entry)
                fig.subplots_adjust(left=.09,right=.97,top=.925,bottom=.115,hspace=.62,wspace=.3)
                caption='Gray shading: cue on. Dashed vertical line: withdrawal; basal recovery follows on the same elapsed-time axis. Blue: old left pole; orange: new right pole. Pole detection uses the current minimum plus half the initial contrast, confined to each end’s half-domain. Absent poles have zero width/material and undefined peak/mean. All times and quantities retain native model units. Full-field tail flags in the titles are distinct from retained scalar reversal.'
                if model=='debelly':caption+=' MCA polarity is bound-attachment polarity; it need not follow the sign of Rac polarity.'
                save(fig,LABELS[model]+' · '+('Fixed patch (support 0.2)' if model=='spring' and geometry=='fixed_physical' else NAMES[geometry])+': cue and recovery',caption,f'{len(pages)+1:02d}_{model}_{geometry}')
                print('Rendered',model,geometry,flush=True)
    assert len(chosen)==45 and len(pages)==16
    pdfpath=OUT/'polarity_recovery.pdf'
    (DATA/'report_provenance.json').write_text(json.dumps(dict(generated_at=STAMP,producer='finished',source_generated_at=E.status['generated_at'],page_count=len(pages),selected_lengths=selection,selected_trajectories=chosen,archive_and_checkpoint_hashes_pass=True,exact_restart_states_pass=True,initial_active_arrays_atol_1e12_pass=True,independent_active_moments_pass=True,pdf_sha256=digest(pdfpath),renderer_sha256=digest(Path(__file__)),trajectory_csv_sha256=digest(csvpath),source_audit_sha256=digest(SOURCE/'results/spatial_audit.json')),indent=2)+'\n')
    sections=''.join(f'<section><h2>{html.escape(p["title"])}</h2><img loading="lazy" src="{p["image"]}" alt="{html.escape(p["title"])}"><p>{html.escape(p["caption"])}</p></section>' for p in pages)
    (OUT/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Polarity and recovery</title><style>body{max-width:1350px;margin:32px auto;padding:0 24px;font:17px/1.6 system-ui;color:#202c35}img{max-width:100%}section{border-top:1px solid #ddd;margin-top:32px}a{color:#176487}</style><h1>Polarity during cue and recovery</h1><p><a href="polarity_recovery.pdf">16-page PDF</a> · <a href="../../data/runs/polarity_recovery_report_20261008/pole_trajectories.csv.gz">All trajectory measurements (compressed CSV)</a></p><p>Generated '+STAMP+'; source simulation producer finished.</p>'+sections+'</html>')
    print(json.dumps(dict(pdf=str(pdfpath),pages=len(pages),trajectories=len(chosen)),indent=2))

if __name__=='__main__':main()
