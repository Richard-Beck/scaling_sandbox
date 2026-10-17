"""Parameter point estimates, geometric group centers, and WGD fold changes."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
from model import OUT,STUDIES,SCENARIOS,SCENARIO_LABELS,decode

NAMES=['Xuan 2018','Mallin 2023','Jemaa 2017','Jemaa 2023','Hosny baseline']
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7','#8a6b00']

def generate():
    fits=[json.loads(p.read_text()) for p in OUT.glob('fit_*.json')]
    for scenario in SCENARIOS:
        candidates=[f for f in fits if f['scenario']==scenario]
        if not candidates:continue
        best=min(candidates,key=lambda f:f['total_score'])
        pars,drift,_=decode(best['z']);values=pars.copy();values[:,:,0]*=60
        means=np.exp(np.log(values).mean(axis=0))
        # Check that plotted geometric means match the profiled hierarchy centers.
        for j,key in enumerate(['tau','sigma']):
            h=best['hyperparameters'][key];factor=60 if j==0 else 1
            assert np.allclose(means[:,j],[factor*h['population_baseline'],factor*h['population_baseline']*h['population_wgd_fold']])
        fig,axs=plt.subplots(1,3,figsize=(15,6))
        for j in (0,1):
            ax=axs[j]
            for i,(name,color) in enumerate(zip(NAMES,COLORS)):
                x=np.array([0,1])+(i-2)*.04
                ax.plot(x,values[i,:,j],'-o',color=color,lw=1,ms=7,alpha=.9)
            ax.plot([0,1],means[:,j],'D--',color='#111111',ms=11,lw=2,zorder=10)
            for g in (0,1):ax.annotate(f'{means[g,j]:.2f}',(g,means[g,j]),xytext=(10,9),textcoords='offset points',fontsize=10,fontweight='bold',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1))
            ax.set_yscale('log');ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.set_title(['Persistence time','Velocity-noise strength'][j]);ax.set_ylabel(['Minutes · logarithmic scale','µm/h^(3/2) · logarithmic scale'][j])
        axs[2].plot([0,1],drift,'-o',color=COLORS[-1],ms=8)
        for g in (0,1):axs[2].annotate(f'{drift[g]:.2f}',(g,drift[g]),xytext=(7,7),textcoords='offset points')
        axs[2].set_ylim(min(0,float(drift.min()))-.15,max(0,float(drift.max()))+.65)
        axs[2].set_title('Hosny signed drift\nOnly one directed assay; unpooled');axs[2].set_ylabel('µm/h');axs[2].axhline(0,color='#aaa',lw=.8)
        for ax in axs:
            ax.set_xticks([0,1],['WGD−','WGD+']);ax.set_xlim(-.3,1.4);ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
        handles=[Line2D([],[],color=c,marker='o',label=n) for c,n in zip(COLORS,NAMES)]
        handles.append(Line2D([],[],color='black',marker='D',linestyle='--',label='Group geometric mean'))
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.9),ncol=3,frameon=False)
        fig.suptitle('Parameter estimates — '+SCENARIO_LABELS[scenario],fontsize=16,y=.99)
        fig.text(.5,.025,'Colored points: experiment-specific estimates; connecting lines pair WGD conditions. Black diamonds: means on the log scale.\nPoint estimates only: no uncertainty intervals. Geometric means are not arithmetic means of a lognormal population.',ha='center',fontsize=10)
        fig.tight_layout(rect=(0,.12,1,.79))
        for ext in ['png','svg']:fig.savefig(OUT/f'parameters_{scenario}.{ext}',dpi=160)
        plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(11,5.5))
        folds=pars[:,1,:]/pars[:,0,:]
        for j,ax in enumerate(axs):
            for i,c in enumerate(COLORS):ax.scatter(folds[i,j],4-i,color=c,s=65,zorder=5)
            pop=best['hyperparameters'][['tau','sigma'][j]]['population_wgd_fold']
            ax.scatter(pop,5.2,color='black',marker='D',s=100,zorder=5)
            ax.annotate(f'{pop:.2f}×',(pop,5.2),xytext=(8,6),textcoords='offset points',fontsize=10)
            ax.axvline(1,color='#999',ls='--');ax.set_xscale('log');ax.grid(axis='x',alpha=.15)
            ax.set_yticks([5.2,4,3,2,1,0],['Group geometric mean']+NAMES);ax.set_ylim(-.5,6)
            ax.set_xlabel('WGD+ / WGD− · logarithmic scale');ax.set_title(['Persistence fold change','Noise-strength fold change'][j]);ax.spines[['top','right']].set_visible(False)
            low,high=ax.get_xlim();ax.set_xlim(min(low,.8),max(high,1.25))
        fig.suptitle('WGD effects — '+SCENARIO_LABELS[scenario],fontsize=14)
        fig.text(.5,.02,'Dashed line: no WGD effect. Each point is a fitted value, not an independent replicate or an uncertainty interval.',ha='center',fontsize=10)
        fig.tight_layout(rect=(0,.07,1,.93))
        for ext in ['png','svg']:fig.savefig(OUT/f'parameter_effects_{scenario}.{ext}',dpi=160)
        plt.close(fig)

if __name__=='__main__':generate()
