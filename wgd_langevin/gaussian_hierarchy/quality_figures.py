"""Observed intervals versus frozen analytical and sampled predictions."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from model import OUT,STUDIES,FIT_OBS,CHECK_OBS,SCENARIOS,SCENARIO_LABELS

LABELS={'speed':'Speed (µm/h)','path_length':'Path length (µm)','endpoint':'Endpoint distance (µm)','directness':'Directness','parallel':'Signed migration along gradient','perpendicular':'Signed migration across gradient','msd':'MSD (µm²)'}
COLORS=['#263746','#e17c19','#157dba']

def generate():
    trackpath=OUT/'tracks_primary.npz'
    if trackpath.exists():
        tracks=np.load(trackpath)
        fig,axes=plt.subplots(2,5,figsize=(16,7))
        for i,study in enumerate(STUDIES):
            lim=max(abs(tracks[f'{study}_{g}']).max() for g in (0,1))*1.05
            for g in (0,1):
                ax=axes[g,i]
                for track in tracks[f'{study}_{g}'][:16]:ax.plot(track[:,0],track[:,1],lw=.8,alpha=.8)
                ax.set(xlim=(-lim,lim),ylim=(-lim,lim),aspect='equal',title=study+(' WGD−' if g==0 else ' WGD+'),xlabel='x (µm)',ylabel='y (µm)')
        fig.suptitle('Frozen analytical fit: 10-minute persistence prior, primary pooling')
        fig.tight_layout();fig.savefig(OUT/'tracks_primary.png',dpi=150);plt.close(fig)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    for scenario in SCENARIOS:
        path=OUT/f'validation_{scenario}.json'
        if not path.exists():continue
        data=json.loads(path.read_text())['comparison']
        for study in STUDIES:
            rows=[r for r in data if r['study']==study]
            hosny=study=='Hosny_baseline'
            keys=list(dict.fromkeys((r['metric'], 'quartiles' if hosny else r['role']) for r in rows))
            cols=min(3,len(keys));nr=(len(keys)+cols-1)//cols
            fig,axs=plt.subplots(nr,cols,figsize=(5*cols,4*nr+1.3),squeeze=False)
            for ax,(metric,role) in zip(axs.flat,keys):
                rr=[r for r in rows if r['metric']==metric and (hosny or r['role']==role)]
                for g in (0,1):
                    group=[r for r in rr if r['group']==g]
                    if hosny:
                        group=sorted(group,key=lambda r:r['quantile'])
                        assert [r['quantile'] for r in group]==[0,1,2]
                        values=[np.array([r['obs'] for r in group]),np.array([r['analytic'] if r['analytic'] is not None else np.nan for r in group]),np.array([(r['tracks_seed1']+r['tracks_seed2'])/2 for r in group])]
                        for j,v in enumerate(values):
                            x=g+(j-1)*.19
                            ax.plot(x,v[1],['o','D','s'][j],color=COLORS[j],ms=6,zorder=4)
                            if np.isfinite(v[[0,2]]).all():
                                ax.vlines(x,v[0],v[2],color=COLORS[j],lw=2,linestyle='--' if j==1 else '-')
                                ax.plot([x,x],[v[0],v[2]],'_',color=COLORS[j],ms=10)
                    else:
                        r=group[0]
                        obslist=(FIT_OBS if role=='fit' else CHECK_OBS)[study]
                        o=next(o for o in obslist if o['group']==g and o['metric']==metric)
                        ax.errorbar(g-.19,r['obs'],yerr=[[r['obs']-o['low']],[o['high']-r['obs']]],fmt='o',color=COLORS[0],capsize=5,lw=2)
                        ax.plot(g,r['analytic'],'D',color=COLORS[1],ms=6)
                        # Both independent validation ensembles are shown, not an uncertainty interval.
                        ax.plot([g+.17,g+.21],[r['tracks_seed1'],r['tracks_seed2']],'s',color=COLORS[2],ms=5)
                title=LABELS[metric]
                if study=='Xuan2018' and metric=='msd':title+=' · '+('11.11 h lag' if role=='fit' else '10 min lag')
                status=('Medians fitted; quartiles held out' if metric in ['path_length','directness','parallel'] else 'All quantiles held out') if hosny else ('Fitted measurement' if role=='fit' else 'Held-out check')
                ax.set_title(title+'\n'+status,fontsize=11,pad=12)
                ax.set_xticks([0,1],['WGD− / reference','WGD+']);ax.set_xlim(-.5,1.5)
                ax.grid(axis='y',alpha=.18)
                if metric in ['parallel','perpendicular']:ax.axhline(0,color='#999',lw=.8)
                else:ax.set_ylim(bottom=0)
            for ax in list(axs.flat)[len(keys):]:ax.set_visible(False)
            fig.suptitle(f'{study.replace("_baseline", " · regular baseline")} — {SCENARIO_LABELS[scenario]}',fontsize=16,y=.99)
            fig.legend(handles=[Line2D([],[],color=c,marker=m,linestyle='none',label=l) for c,m,l in zip(COLORS,['o','D','s'],['Observed','Analytical approximation','Generated tracks'])],loc='upper center',bbox_to_anchor=(.5,.945),ncol=3,frameon=False)
            note=('Points: medians. Bars: 25th–75th percentiles, not confidence intervals.\nTrack quantiles averaged across two ensembles; no analytical path-length quartiles available.' if hosny else 'Observed bars: scraped reading intervals, not confidence intervals.\nBlue squares: two independent track ensembles (often overlapping).')
            fig.text(.5,.018,note,ha='center',fontsize=10,color='#445566')
            fig.tight_layout(rect=(0,.10,1,.86 if nr==1 else .9))
            for ext in ['png','svg']:fig.savefig(OUT/f'quality_{scenario}_{study}.{ext}',dpi=160)
            plt.close(fig)

if __name__=='__main__':generate()
