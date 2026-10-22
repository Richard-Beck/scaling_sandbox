"""Standalone figure assets and measured summary for the mechanism report."""
from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter

MANUSCRIPT=Path(__file__).resolve().parents[2]
OUT=MANUSCRIPT/'data/runs/fixed_slope_mechanisms_20261006/wave_pinning'
FIG=MANUSCRIPT/'figures/fixed_slope_mechanisms/wave_pinning'
FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'legend.fontsize':7.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
COLORS=['#27627e','#247a51','#a64b35']
def load(k):
    with np.load(OUT/f'trajectory_{k:02d}.npz') as z:return {a:z[a].copy() for a in z.files}
Z={k:load(k) for k in range(8)}
def snap(z,t):return int(np.argmin(abs(z['time']-t)))
def finish(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
def panel(ax,label):ax.text(-.14,1.07,label,transform=ax.transAxes,fontweight='bold',fontsize=11)
def timeaxis(ax):
    ax.set_xscale('symlog',linthresh=20.)
    ax.set_xticks([0,25,100,1000,10000,50000]);ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.tick_params(axis='x',labelsize=7)
    ax.set_xlabel('Time after gradient withdrawal (native units)')

def spatial():
    fig,axes=plt.subplots(3,3,figsize=(11.5,8.2),layout='constrained')
    fig.suptitle('Wave-pinning: stronger long-cell input starts activation faster, but leaves competing domains',fontsize=14)
    labels=['Short / decreasing arm','Minimum duration','Long / increasing arm']
    for col,k in enumerate([0,2,4]):
        z=Z[k];L=float(z['length']);D=float(z['pulse_duration']);x=z['x'];color=COLORS[col]
        ax=axes[0,col];ax.plot(x,z['prep_active'],color=color,lw=2,label='Initial left-active state')
        ax.set_title(f'{labels[col]}\nL = {L:g}; shortest orientation pulse = {D:g}')
        ax.set_ylabel('Active protein');ax.set_ylim(0,1.3);ax.set_xlim(0,L)
        b=ax.twinx();b.spines['right'].set_visible(True);b.plot(x,float(z['slope'])*x,color='#c78a22',lw=1.6,ls='--',label='Added gradient')
        b.set_ylim(0,.28);b.set_ylabel('Added activation input',color='#a36f19');b.tick_params(axis='y',colors='#a36f19');ax.legend(loc='upper left',frameon=False)
        ax.text(.98,.10,f'Ramp peak = {float(z["slope"])*L:.3f}',transform=ax.transAxes,ha='right',va='bottom',fontsize=8,color='#a36f19')
        ax.set_xlabel('Position (native length units)');panel(ax,chr(65+col))
        ax=axes[1,col];times=[0.,D/2,D];cs=['#b7b7b7','#6ca1b5',color]
        for t,c in zip(times,cs):
            i=snap(z,t);ax.plot(x,z['active'][i],color=c,lw=1.8,label=f't = {float(z["time"][i]):g}'+(' (cutoff)' if t==D else ''))
        start=float(z['time'][np.flatnonzero(z['active'][:,-1]>.6)[0]])
        ax.set_ylim(0,1.3);ax.set_xlim(0,L);ax.set_ylabel('Active protein');ax.set_xlabel('Position');ax.set_title(f'During gradient; right edge first > 0.6 at t = {start:g}',fontsize=9);ax.legend(frameon=False,loc='upper left');panel(ax,chr(68+col))
        ax=axes[2,col];times=[0,25,100,1000,10000];cs=['#a4a4a4','#dfb677','#969cc3','#448783','#173e52']
        for t,c in zip(times,cs):
            i=snap(z,D+t);ax.plot(x,z['active'][i],color=c,lw=1.6,label=f'withdrawal + {t:g}')
        ax.set_ylim(0,1.3);ax.set_xlim(0,L);ax.set_ylabel('Active protein');ax.set_xlabel('Position');ax.set_title('Recovery under basal input');ax.legend(frameon=False,loc='upper left');panel(ax,chr(71+col))
    finish(fig,'wave_pinning_spatial')

def decisions():
    fig,axes=plt.subplots(3,3,figsize=(11.5,8.3),layout='constrained')
    fig.suptitle('A slightly shorter pulse fails for three different reasons',fontsize=14)
    for col,(pk,fk) in enumerate([(0,1),(2,3),(4,5)]):
        for row,k in enumerate([pk,fk]):
            z=Z[k];D=float(z['pulse_duration']);L=float(z['length']);mask=z['time']>=D;t=z['time'][mask]-D
            keep=np.unique(np.r_[np.flatnonzero(mask)[0],np.flatnonzero(mask)[::4],len(z['time'])-1]);t=z['time'][keep]-D
            ax=axes[row,col];a=ax.pcolormesh(z['x'],np.log10(1+t),z['active'][keep],cmap='viridis',vmin=.17,vmax=1.15,shading='auto',rasterized=True)
            ticks=[0,25,100,1000,10000] + ([50000] if L==39 else [])
            ax.set_yticks(np.log10(1+np.array(ticks)),[str(t) for t in ticks]);ax.set_ylim(0,np.log10(1+(50000 if L==39 else 10000)))
            ax.set_ylabel('Time after withdrawal');ax.set_xlabel('Position');ax.set_title(f'L = {L:g}; '+('passing' if row==0 else 'failing')+f' pulse {D:g}')
            ax.axhline(np.log10(10001),color='white',ls='--',lw=.8);panel(ax,chr(65+row*3+col))
        ax=axes[2,col]
        for k,ls in [(pk,'-'),(fk,'--')]:
            z=Z[k];D=float(z['pulse_duration']);mask=z['time']>=D
            ax.plot(z['time'][mask]-D,z['orientation'][mask],color=COLORS[col],ls=ls,lw=1.6,label=f'Pulse {D:g}')
        ax.axhline(.05,color='black',ls=':',lw=.9,label='Orientation cutoff 0.05');ax.axhline(0,color='#ccc',lw=.6)
        ax.axvspan(9500,10000,color='#bbbbbb',alpha=.25,label='Official final-500 window')
        ax.set_xlim(0,50000 if col==2 else 10000);timeaxis(ax);ax.set_ylim(-.46,.46) if col<2 else ax.set_ylim(.025,.125)
        ax.set_ylabel('Active-protein orientation');ax.legend(frameon=False,loc='best');panel(ax,chr(71+col))
        texts=['30.5 → homogeneous; 30.75 → right cap','13.75 → left cap; 13.875 → right cap','38.75 crosses 0.05 after the official horizon']
        ax.set_title(texts[col],fontsize=9)
    fig.colorbar(a,ax=axes[:2,:],label='Active protein',shrink=.6,location='right')
    finish(fig,'wave_pinning_decisions')

def mechanisms():
    fig,axes=plt.subplots(2,3,figsize=(11.5,7.6),layout='constrained')
    fig.suptitle('Wave-pinning: finite pool redistribution creates a long-lived interior competitor',fontsize=14)
    ax=axes[0,0]
    rows=list(csv.DictReader((MANUSCRIPT/'data/runs/fixed_slope_reversal_20261006/results/endpoints.csv').open()))
    for gate,ls,label in [('orientation','-','Persistent orientation'),('full_polarity','--','1% scalar strength')]:
        a=[r for r in rows if r['model']=='wave_pinning' and r['condition']=='relative_20pct' and r['definition']==gate and r['status']=='success']
        ax.plot([float(r['length']) for r in a],[float(r['endpoint']) for r in a],ls=ls,color='#374456',label=label)
        for row in a:
            field=row['selected_tail_field_pass_0_1_percent']=='True'
            ax.plot(float(row['length']),float(row['endpoint']),'o',ms=3,markerfacecolor='#374456' if field else 'white',markeredgecolor='#374456')
    for k,c in zip([0,2,4],COLORS):
        z=Z[k];ax.scatter(float(z['length']),float(z['pulse_duration']),c=c,s=60,zorder=5)
    ax.set_xlabel('Length');ax.set_ylabel('Shortest passing tested pulse');ax.legend(frameon=False);ax.set_title('Scalar readouts; open = field-tail flag fails',fontsize=9);panel(ax,'A')
    ax=axes[0,1]
    for k,c in zip([0,2,4],COLORS):
        z=Z[k];D=float(z['pulse_duration']);i=snap(z,D);ax.plot(z['x']/float(z['length']),z['inactive'][i],color=c,lw=2,label=f'L={float(z["length"]):g}, pulse={D:g}')
    ax.set_xlabel('Relative position x/L');ax.set_ylabel('Inactive protein at cutoff');ax.set_title('Long ramp depletes the pool unevenly');ax.legend(frameon=False);panel(ax,'B')
    ax=axes[0,2]
    for k,c,label in [(4,'#27627e','Native slope; peak 0.033'),(0,'#278363','Slope ×2.6; peak 0.087')]:
        with np.load(OUT/f'cue_control_{k:02d}.npz') as q:
            D=float(q['duration']);mask=(q['time']>=D)&(q['time']<=D+1000)
            ax.plot(q['time'][mask]-D,q['orientation'][mask],color=c,lw=1.8,label=label)
    ax.axhline(.05,color='#888',ls=':',lw=.8);ax.set_xlabel('Time after withdrawal');ax.set_ylabel('Orientation');ax.set_title('Short arm: L=5, same pulse 13.875',fontsize=9);ax.legend(frameon=False,loc='lower right');panel(ax,'C')
    ax=axes[1,0]
    for k,ls in [(6,'-'),(7,'--')]:
        z=Z[k];D=float(z['pulse_duration']);mask=z['time']>=D
        ax.plot(z['time'][mask]-D,z['orientation'][mask],color='#a64b35',ls=ls,lw=1.8,label=f'Pulse {D:g}')
    ref=abs(json.loads((OUT/'trajectory_06.json').read_text())['job']['case']['initial_orientation'])
    ax.axhspan(.99*ref,1.01*ref,color='#43885a',alpha=.2,label='1% initial-strength band');ax.axvspan(9500,10000,color='#bbb',alpha=.25)
    ax.set_xlim(0,15000);timeaxis(ax);ax.set_ylim(.14,.46);ax.set_ylabel('Orientation');ax.set_title('Full-strength threshold includes coarsening');ax.legend(frameon=False,loc='upper left');panel(ax,'D')
    ax=axes[1,1]
    z=Z[4];D=float(z['pulse_duration']);mask=(z['time']>=D)&(z['time']<=D+1000)
    ax.plot(z['time'][mask]-D,z['orientation'][mask],color='#a64b35',lw=2,label='Native inactive diffusion = 10')
    for k,c,label in [(4,'#278363','Inactive diffusion = 100'),(0,'#3967aa','Inactive diffusion = 1000')]:
        with np.load(OUT/f'control_{k:02d}.npz') as q:
            m=q['recovery_time']<=1000;ax.plot(q['recovery_time'][m],q['orientation'][m],color=c,lw=1.5,label=label)
    ax.set_xlim(0,1000);ax.set_xlabel('Time after withdrawal');ax.set_ylabel('Orientation');ax.set_title('Same cutoff; diffusion changes only after withdrawal',fontsize=9);ax.legend(frameon=False,loc='lower right');panel(ax,'E')
    ax=axes[1,2]
    for k,c,label in [(4,'#a64b35','Native'),(4,'#278363','Inactive diffusion ×10')]:
        if label=='Native':
            z=Z[k];i=snap(z,float(z['pulse_duration'])+1000);x=z['x'];u=z['active'][i]
        else:
            with np.load(OUT/'control_04.npz') as q:i=np.argmin(abs(q['recovery_time']-1000));x=q['x'];u=q['active'][i]
        ax.plot(x,u,color=c,lw=2,label=label)
    with np.load(OUT/'mesh_control_00.npz') as q:
        i=np.argmin(abs(q['time']-float(q['duration'])-1000));ax.plot(q['x'],q['active'][i],color='#5b5b5b',ls='--',lw=1,label='Native, 401-cell mesh check')
    ax.set_xlabel('Position');ax.set_ylabel('Active protein');ax.set_ylim(0,1.3);ax.set_title('At withdrawal +1000: competitor or one cap');ax.legend(frameon=False);panel(ax,'F')
    finish(fig,'wave_pinning_mechanism')

def full_snapshots():
    fig,axes=plt.subplots(2,2,figsize=(11.5,7.6),layout='constrained')
    fig.suptitle('L=39: the full-strength threshold depends on when the interior domain disappears',fontsize=14)
    for col,k in enumerate([6,7]):
        z=Z[k];D=float(z['pulse_duration']);ax=axes[0,col]
        for target,c in zip([25,1000,9000,9500,10000],['#d7b689','#899eaf','#d28b66','#a8483d','#173e52']):
            i=snap(z,D+target);ax.plot(z['x'],z['active'][i],color=c,lw=1.7,label=f'withdrawal + {target:g}')
        ax.set_xlabel('Position');ax.set_ylabel('Active protein');ax.set_ylim(0,1.3);ax.set_title(f'Pulse {D:g}: '+('passes both scalar gates' if col==0 else 'fails both scalar gates at official horizon'))
        ax.legend(frameon=False,loc='upper left');panel(ax,'A' if col==0 else 'B')
        ax=axes[1,col]
        mask=(z['time']>=D)&(z['time']<=D+15000);idx=np.flatnonzero(mask)[::2]
        a=ax.pcolormesh(z['x'],z['time'][idx]-D,z['active'][idx],cmap='viridis',vmin=.17,vmax=1.15,shading='auto',rasterized=True)
        ax.axhspan(9500,10000,facecolor='white',alpha=.3);ax.axhline(9500,color='white',ls='--',lw=.7);ax.axhline(10000,color='white',ls='--',lw=.7)
        ax.set_xlabel('Position');ax.set_ylabel('Time after withdrawal');ax.set_title('Shaded strip: official final-500 gate window');panel(ax,'C' if col==0 else 'D')
    fig.colorbar(a,ax=axes[1,:],label='Active protein',shrink=.8,location='right')
    finish(fig,'wave_pinning_full_strength')

def summarize():
    rows=[]
    for k in range(8):
        z=Z[k];D=float(z['pulse_duration']);L=float(z['length']);j=json.loads((OUT/f'trajectory_{k:02d}.json').read_text())
        row=dict(slot=k,length=L,duration=D,gates=j['job']['gates'],**j['metrics'])
        row['right_edge_first_above_0_6']=float(z['time'][np.flatnonzero(z['active'][:,-1]>.6)[0]])
        i=snap(z,D);row['inactive_cutoff_range']=float(np.ptp(z['inactive'][i]));row['mean_inactive_cutoff']=float(z['inactive'][i].mean());row['cutoff_orientation']=float(z['orientation'][i])
        row['cutoff_inactive_left']=float(z['inactive'][i,0]);row['cutoff_inactive_right']=float(z['inactive'][i,-1])
        snapshots=[]
        for target in [D,D+25,D+100,D+1000,D+10000,D+50000]:
            i=snap(z,target);u=z['active'][i];idx=np.flatnonzero(u>.6);components=np.split(idx,np.flatnonzero(np.diff(idx)>1)+1)
            snapshots.append(dict(time=float(z['time'][i]),recovery_time=float(z['time'][i]-D),orientation=float(z['orientation'][i]),regions_active_above_0_6=[dict(left=float(z['x'][c[0]]),right=float(z['x'][c[-1]])) for c in components if len(c)]))
        row['snapshots']=snapshots;rows.append(row)
    result=dict(generated_at=json.loads((OUT/'trajectory_07.json').read_text())['generated_at'],native_reconstruction_audit='All8 replayed scalar gates/field flags and final orientations match frozen trial records. Mass conserved to<=1.46e-13 relative.',rows=rows,
                controls=[],cue_controls=[],mesh_controls=[])
    for k in range(8):
        with np.load(OUT/f'control_{k:02d}.npz') as q:
            result['controls'].append(dict(slot=k,du=float(q['du']),dv=float(q['dv']),source_slot=int(q['source_slot']),orientation_100=float(q['orientation'][np.argmin(abs(q['recovery_time']-100))]),orientation_1000=float(q['orientation'][np.argmin(abs(q['recovery_time']-1000))]),orientation_10000=float(q['orientation'][-1])))
    for k in range(5):result['cue_controls'].append(json.loads((OUT/f'cue_control_{k:02d}.json').read_text()))
    for k in range(2):result['mesh_controls'].append(json.loads((OUT/f'mesh_control_{k:02d}.json').read_text()))
    (OUT/'mechanism_summary.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':spatial();decisions();mechanisms();full_snapshots();summarize()
