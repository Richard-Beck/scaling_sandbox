"""Static exportable comparisons, with observations and excluded metrics explicit."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from publish import OUT,FAMILIES,publish,status,atomic
COLORS=['#0072b2','#e69f00','#009e73','#cc79a7']
LABELS=['Gaussian','Memory','Coupled','Self-activation']
UNITS=dict(speed='µm/h',msd='µm²',path_length='µm',endpoint='µm',directness='',parallel='',perpendicular='')
def load():
    manifest=json.loads((OUT/'manifest.json').read_text())
    return manifest,{j['key']:json.loads((OUT/(j['key']+'.json')).read_text()) for j in manifest['jobs'] if (OUT/(j['key']+'.json')).exists()}
def main():
    manifest,data=load();studies=list(dict.fromkeys(j['study'] for j in manifest['jobs']))
    for study in studies:
        results=[data.get(f+'_persistence_'+study) for f in FAMILIES]
        available=[r for r in results if r]
        if not available:continue
        first=available[0];metrics=list(dict.fromkeys((role,o['metric']) for role in ('rows','checks') for o in first[role]))
        fig,axes=plt.subplots(1,len(metrics),figsize=(max(9,4*len(metrics)),4.4),squeeze=False)
        for ax,(role,metric) in zip(axes[0],metrics):
            base=[o for o in first[role] if o['metric']==metric]
            for o in base:
                x=o['group']+(o.get('quantile',1)-1)*.24
                if 'low' in o:
                    ax.errorbar(x,o['obs'],yerr=[[o['obs']-o['low']],[o['high']-o['obs']]],fmt='ks',capsize=5,ms=5,zorder=5)
                else:ax.scatter(x,o['obs'],c='black',s=32,marker='s',zorder=5)
            for fi,r in enumerate(results):
                if not r:continue
                for i,o in enumerate(r[role]):
                    if o['metric']!=metric:continue
                    x=o['group']+(o.get('quantile',1)-1)*.24+(fi-1.5)*.035
                    ax.scatter(x,o['predicted'],c=COLORS[fi],marker='o',s=32,zorder=4)
                    if role=='rows':
                        values=[v['rows'][i]['predicted'] for v in r['validations']]
                        ax.plot([x,x],[min(values),max(values)],color=COLORS[fi],lw=3,alpha=.6)
                        ax.scatter([x]*len(values),values,c=COLORS[fi],marker='_',s=45)
            ax.set_xticks([0,1],['WGD−','WGD+']);ax.set_xlim(-.5,1.5);ax.set_ylabel(UNITS[metric]);ax.grid(axis='y',alpha=.2)
            label=metric.replace('_',' ')
            if study=='Xuan2018' and metric=='msd':label+=' at '+('11.11 h' if role=='rows' else '1/6 h')
            ax.set_title(label+'\n'+('fitted' if role=='rows' else 'excluded check'),fontsize=11)
        handles=[Line2D([],[],marker='s',color='black',ls='',label='Observed / reading interval')]+[Line2D([],[],marker='o',color=c,ls='',label=l) for c,l in zip(COLORS,LABELS)]
        fig.legend(handles=handles,loc='lower center',ncol=5,frameon=False,fontsize=9)
        subtitle=' — vertical coloured ranges show independent numerical checks'
        if study.startswith('Hosny'):subtitle+='\nWithin each group: left Q25, middle median, right Q75'
        fig.suptitle(study+subtitle,fontsize=13)
        fig.tight_layout(rect=(0,.09,1,.91));fig.savefig(OUT/(study+'.png'),dpi=160);fig.savefig(OUT/(study+'.pdf'));plt.close(fig)
    scores=np.full((len(studies),4),np.nan);gains=scores.copy()
    for i,s in enumerate(studies):
        for j,f in enumerate(FAMILIES):
            r=data.get(f+'_persistence_'+s);n=data.get(f+'_shared_'+s)
            if r:scores[i,j]=r['score']
            if r and n and n['score']>0:gains[i,j]=100*(1-r['score']/n['score'])
    fig,axes=plt.subplots(1,2,figsize=(13,7),layout='constrained')
    for ax,mat,title,cmap in [(axes[0],np.log10(1+scores),'log10(1 + descriptive score)','viridis_r'),(axes[1],gains,'Improvement over shared null (%)','viridis')]:
        im=ax.imshow(mat,aspect='auto',cmap=cmap);fig.colorbar(im,ax=ax,shrink=.7)
        ax.set_xticks(range(4),LABELS,rotation=20,ha='right');ax.set_yticks(range(len(studies)),studies);ax.set_title(title)
        for i in range(len(studies)):
            for j in range(4):
                if np.isfinite(mat[i,j]):ax.text(j,i,f'{scores[i,j]:.2g}' if ax==axes[0] else f'{gains[i,j]:.1f}',ha='center',va='center',color='black',bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=2))
    fig.suptitle('Persistence-only fits: scores compare variants within each experiment')
    fig.savefig(OUT/'comparison.png',dpi=170);fig.savefig(OUT/'comparison.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,7),layout='constrained')
    for ax,name in zip(axes,['tp','tm']):
        for i,s in enumerate(studies):
            for j,f in enumerate(FAMILIES):
                r=data.get(f+'_persistence_'+s)
                if r and name in r['ratios']:ax.scatter(r['ratios'][name],i+(j-1.5)*.12,c=COLORS[j],s=30)
        ax.axvline(1,c='gray',ls='--');ax.set_xscale('log');ax.set_yticks(range(len(studies)),studies);ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
        ax.set_title(name+': WGD+ / WGD−');ax.set_xlabel('Best-fit ratio (not an identified effect or confidence interval)')
    fig.legend(handles=[Line2D([],[],marker='o',color=c,ls='',label=l) for c,l in zip(COLORS,LABELS)],loc='outside lower center',ncol=4,frameon=False)
    fig.savefig(OUT/'persistence_ratios.png',dpi=170);fig.savefig(OUT/'persistence_ratios.pdf');plt.close(fig)
    metadata={s:sorted(k for k,r in data.items() if r['study']==s and r['indices']) for s in studies}
    metadata['all_keys']=sorted(data)
    atomic('figure_manifest.json',json.dumps(metadata,indent=2)+'\n')
    publish()
if __name__=='__main__':main()
