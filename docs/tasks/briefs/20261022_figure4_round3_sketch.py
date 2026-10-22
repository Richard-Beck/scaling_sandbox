import json,csv,math,collections,sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
D=sys.argv[1]
plt.rcParams.update({'font.size':7,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
 'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
sizes=[20,30,40,50,60]
cmap=plt.get_cmap('Blues'); scol={s:cmap(0.35+0.6*i/4) for i,s in enumerate(sizes)}
rng=np.random.default_rng(3)
fig=plt.figure(figsize=(7.2,9.6))
gs=GridSpec(4,1,figure=fig,height_ratios=[1.15,0.8,1.0,1.0],hspace=0.55)
T=24000; dt=20; nt=T//dt
def sim(C,Ds,lam):
    x=np.empty(nt); x[0]=rng.uniform(0,C); jumps=[]
    for i in range(1,nt):
        if rng.random()<lam*dt: x[i]=rng.uniform(0,C); jumps.append(i)
        else: x[i]=(x[i-1]+math.sqrt(2*Ds*dt)*rng.standard_normal())%C
    return x,jumps
def kymo(ax,x,C):
    grid=np.linspace(0,C,int(C*4))
    d=np.abs(grid[None,:]-x[:,None]); d=np.minimum(d,C-d)
    ax.imshow(np.exp(-0.5*(d/1.5)**2),aspect='auto',cmap='Greys',origin='lower',extent=[0,C,0,T/1000],interpolation='nearest')
    ax.set_xticks([]); ax.tick_params(length=2,pad=1)
Ds=2e-3; lam0=1/7000
nulls=[('Drifting peak\n(same physical $D_s$)',lambda C:(Ds,0)),
       ('Jumping peak\n(same rate)',lambda C:(0,lam0)),
       ('Jumping peak\n(rate ∝ perimeter)',lambda C:(0,lam0*C/40))]
sub=gs[0].subgridspec(1,9,width_ratios=[1,3,0.55,1,3,0.55,1,3,0.0001],wspace=0.15)
col=0
for j,(name,par) in enumerate(nulls):
    for C in (20,60):
        ax=fig.add_subplot(sub[0,col]); col+=1
        x,_=sim(C,*par(C)); kymo(ax,x,C)
        if col==1: ax.set_ylabel('time (×10³)')
        else: ax.set_yticklabels([])
        ax.set_xlabel('C%d'%C,labelpad=1)
        if C==20: ax.text(0,1.04,name,transform=ax.transAxes,va='bottom',ha='left',fontsize=7)
    col+=1
fig.text(0.04,0.955,'A',fontsize=10,fontweight='bold')
fig.text(0.13,0.955,'sketch of A: null simulations, small (C20) and large (C60) rings, width ∝ circumference',fontsize=6.5,color='0.35')
sub=gs[1].subgridspec(1,3,wspace=0.35); ell=np.linspace(0,24000,200)
for j,(name,par) in enumerate(nulls):
    ax=fig.add_subplot(sub[0,j])
    for C,ls in ((20,'-'),(60,':')):
        d,l=par(C); k=(2*math.pi/C)**2*d+l
        ax.plot(ell/1000,np.exp(-k*ell),color='k',ls=ls,lw=1.2,label='C%d'%C)
    ax.set_ylim(0,1.02); ax.set_xlabel('lag ℓ (×10³)',labelpad=1)
    if j==0: ax.set_ylabel('heading autocorr. C(ℓ)'); ax.legend(frameon=False,fontsize=6,loc='lower left')
    k20=(2*math.pi/20)**2*par(20)[0]+par(20)[1]; k60=(2*math.pi/60)**2*par(60)[0]+par(60)[1]
    a=-math.log(k60/k20)/math.log(3)
    ax.text(0.97,0.95,('α = %+.0f'%a).replace('-0','0')+('\n(curves coincide)' if abs(a)<0.1 else ''),transform=ax.transAxes,ha='right',va='top')
f=np.logspace(-3,-0.5,200); Lref=2000; a=8.8; b=135.
def k_of(C,f,mode):
    kd=(2*math.pi/C)**2*a*f**2
    if mode=='drift': return kd
    if mode=='drift+const': return kd+b*f**4
    if mode=='drift+perim': return kd+b*(C/40)*f**4
modes=[('drift only','drift'),('drift + same-rate jumps','drift+const'),('drift + per-perimeter jumps','drift+perim')]
sub=gs[2].subgridspec(1,4,wspace=0.45,width_ratios=[1,1,1,1.1])
for j,(lab,mode) in enumerate(modes):
    ax=fig.add_subplot(sub[0,j])
    for s in sizes:
        h=1.02
        al=-(np.log(k_of(s*h,f,mode))-np.log(k_of(s/h,f,mode)))/np.log(h*h)
        ax.plot(f,al,color=scol[s],lw=1.2)
    for g in (2,0,-1): ax.axhline(g,color='0.75',lw=0.6,zorder=0)
    ax.set_xscale('log'); ax.set_ylim(-1.3,2.5); ax.set_yticks([-1,0,1,2])
    ax.set_xlabel('input f',labelpad=1); ax.set_title(lab,fontsize=7,loc='left')
    if j==0:
        ax.set_ylabel('local size exponent α');
        for s in (20,60): ax.text(1.2e-3,2.12 if s==20 else 1.85,'',)
ax=fig.add_subplot(sub[0,3])
Cg=np.linspace(20,60,41)
for lab,mode in modes:
    fh=[]
    for C in Cg:
        ff=np.logspace(-3,0,4000); k=k_of(C,ff,mode); fh.append(ff[np.argmax(np.exp(-k*Lref)<0.5)])
    ax.plot(Cg,fh,color='k',lw=1.2,ls={'drift':'-','drift+const':'--','drift+perim':':'}[mode],label=lab)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks([20,30,40,60]); ax.set_xticklabels(['20','30','40','60'])
ax.set_xlabel('circumference C',labelpad=1); ax.set_ylabel('noise tolerance f½'); ax.legend(frameon=False,fontsize=5.5,loc='upper left')
yB=gs[2].get_position(fig).y1+0.035; fig.text(0.04,yB,'B',fontsize=10,fontweight='bold')
fig.text(0.13,yB,'what C (left three) and E (right) would look like under the nulls; assumes D_s ∝ f², jump rate ∝ f⁴',fontsize=6.5,color='0.35')
ad=json.load(open(sys.argv[2] if len(sys.argv)>2 else 'alpha_data.json'))
models=[('wave_pinning','Wave-pinning',1),('holmes_model3','Holmes',1),('spring_ring','Spring',20),('goryachev','Goryachev',1)]
sub=gs[3].subgridspec(1,5,wspace=0.45,width_ratios=[1,1,1,1,1.1])
for j,(m,lab,sc) in enumerate(models):
    ax=fig.add_subplot(sub[0,j])
    pts=[o for o in ad if o['model']==m]
    for s in sizes:
        P=sorted([o for o in pts if abs(o['size']*sc-s)<1e-6],key=lambda o:o['f'])
        if not P: continue
        ax.plot([o['f'] for o in P],[o['alpha'] for o in P],'-',color=scol[s],lw=0.8,alpha=0.9)
        for o in P:
            ax.plot(o['f'],o['alpha'],'o',ms=3.2,mfc=scol[s] if not o['forced'] else 'white',mec=scol[s],mew=0.8)
    for g in (2,0,-1): ax.axhline(g,color='0.75',lw=0.6,zorder=0)
    ax.set_xscale('log'); ax.set_ylim(-2,5.2); ax.set_xlim(2e-4,1.5)
    ax.set_title(lab,fontsize=7,loc='left'); ax.set_xlabel('input f',labelpad=1)
    if j==0: ax.set_ylabel('local size exponent α\n(rough, existing data)')
ax=fig.add_subplot(sub[0,4])
fh=collections.defaultdict(dict)
for r in csv.DictReader(open(D+'/main_noise_tolerance.csv')):
    if r['arm']=='common': fh[r['model']][float(r['C_label'])]=float(r['fhalf'])
cl=collections.defaultdict(dict)
for r in csv.DictReader(open(D+'/condition_lag_metrics.csv')):
    if r['arm']=='common' and r['view']=='main' and r['model']=='goryachev' and float(r['lag'])==2000:
        try: cl[float(r['C_label'])][float(r['f'])]=float(r['C'])
        except: pass
for s,dd in cl.items():
    fs=sorted(dd)
    for f0,f1 in zip(fs,fs[1:]):
        if dd[f0]>=0.5>dd[f1]:
            t=(dd[f0]-0.5)/(dd[f0]-dd[f1]); fh['goryachev'][s]=math.exp(math.log(f0)+t*math.log(f1/f0)); break
mk={'wave_pinning':'o','holmes_model3':'s','spring_ring':'^','goryachev':'D'}
for m,lab,sc in models:
    xs=sorted(fh[m]); y=[fh[m][x] for x in xs]; xs=[x*sc for x in xs]
    y=np.array(y)/y[0]
    ax.plot(xs,y,'-',marker=mk[m],ms=3.2,color='k',lw=0.8,mfc='white' if m=='goryachev' else 'k',label=lab)
ax.plot([20,60],[1,3],color='0.6',lw=0.8,ls='--'); ax.text(58,1.9,'f½ ∝ C\n(dashed)',fontsize=5.5,color='0.4',ha='right')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks([20,30,40,60]); ax.set_xticklabels(['20','30','40','60'])
ax.set_yticks([1,2,3,4]); ax.set_yticklabels(['1','2','3','4'])
ax.set_xlabel('circumference C (Spring ×20)',labelpad=1); ax.set_ylabel('f½ / f½(C20)'); ax.legend(frameon=False,fontsize=5.5,loc='upper left')
yC=gs[3].get_position(fig).y1+0.035; fig.text(0.04,yC,'C',fontsize=10,fontweight='bold')
fig.text(0.13,yC,'existing round-2 data: rough local α (3-size fit, lag chosen per condition; open = forced); f½ fold vs size (Goryachev open, forced included)',fontsize=6.5,color='0.35')
kx=fig.add_axes([0.80,0.957,0.16,0.012]);
for i,s in enumerate(sizes): kx.add_patch(plt.Rectangle((i,0),0.9,1,color=scol[s])); kx.text(i+0.45,1.3,'C%d'%s,ha='center',fontsize=5.5)
kx.set_xlim(0,5); kx.set_ylim(0,1); kx.axis('off')
fig.savefig('figure4_nulls_sketch.png',dpi=220,bbox_inches='tight')
