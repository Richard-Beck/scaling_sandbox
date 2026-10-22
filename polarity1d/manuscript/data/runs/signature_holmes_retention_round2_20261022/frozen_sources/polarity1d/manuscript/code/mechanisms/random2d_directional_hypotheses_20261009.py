"""Event-first assay design with one stochastic synthetic event/track dataset.

No model fits or biological observations. Real HT1080 geometry is a template
only. No mechanism-specific size effect is imposed. Current and historical
predictor windows are separated from every future outcome window.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, FancyArrowPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
FIG = ROOT/'figures/random2d_assay_mockup_20261009/directional_v3'
DATA = ROOT/'data/random2d_assay_mockup_20261009/directional_v3'
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
INK, MUTED, TEAL, ORANGE, PURPLE = '#203342', '#63727c', '#008c91', '#d46a35', '#8054a1'
RULE, PALE = '#dce4e9', '#f6f8fa'
COLORS = {'small': TEAL, 'large': ORANGE}
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11,
    'text.color':INK, 'axes.labelcolor':INK, 'axes.edgecolor':'#a4b4be',
    'xtick.color':MUTED, 'ytick.color':MUTED, 'axes.spines.top':False,
    'axes.spines.right':False, 'xtick.labelsize':9, 'ytick.labelsize':9,
    'axes.labelsize':10, 'pdf.fonttype':42, 'svg.fonttype':'none',
    'savefig.facecolor':'white'})
SEED = 912713
rng = np.random.default_rng(SEED)
dt, horizon = .20, 60.
times = np.arange(0, horizon+dt/2, dt)
events, cells, snapshots, predictions, tracks = [], [], [], [], {}


def wrapped(a):
    return np.angle(np.exp(1j*a))


for group in ['small', 'large']:
    for ci in range(24):
        ident = f'{group[0].upper()}{ci+1:02d}'
        perimeter = rng.lognormal(np.log(88 if group=='small' else 145), .15)
        cell_area = perimeter**2 / (4*np.pi) * rng.uniform(.57, .78)
        # Rate per perimeter is identically distributed in the two size bins.
        rate = .0065 * perimeter * rng.lognormal(0, .32)
        n = rng.poisson(rate*(horizon+5))
        starts = np.sort(rng.uniform(-5, horizon, n))
        durations = rng.lognormal(np.log(1.15)+rng.normal(0,.24), .67, n)
        peaks = rng.lognormal(np.log(9)+rng.normal(0,.25), .78, n)
        widths = rng.lognormal(np.log(3.6), .48, n)
        baseline_heading = rng.uniform(-np.pi, np.pi)
        # Directions are mixed: front-biased and freely oriented events coexist.
        dirs = np.where(rng.random(n)<.58,
                        rng.vonmises(baseline_heading, 1.5, n),
                        rng.uniform(-np.pi, np.pi, n))
        active_area = np.zeros((len(times), 2))
        total_area = np.zeros(len(times))
        local_events = []
        for j,(start,duration,peak,width,angle) in enumerate(zip(starts,durations,peaks,widths,dirs)):
            end = start+duration
            if end <= 0:
                continue
            phase = (times-start)/duration
            pulse = np.where((phase>0)&(phase<1),
                peak*np.sin(np.pi*np.clip(phase,0,1))**1.1, 0.)
            pulse *= np.clip(1+rng.normal(0,.14,len(times)), .35, 1.65)
            active_area += pulse[:,None]*np.array([np.cos(angle),np.sin(angle)])
            total_area += pulse
            row = dict(cell=ident, size_group=group, event=f'{ident}_{j:03d}',
                start_min=float(start), end_min=float(end), observed_start_min=float(max(0,start)),
                observed_end_min=float(min(horizon,end)),
                observed_duration_min=float(min(horizon,end)-max(0,start)),
                left_censored=bool(start<0), right_censored=bool(end>horizon),
                synthetic_duration_min=float(duration), footprint_width_um=float(width),
                observed_peak_area_um2=float(pulse.max()), direction_rad=float(angle))
            events.append(row)
            local_events.append((row,pulse))

        # Stochastic motility with history dependence by construction, identical
        # memory-time distributions across size bins. It is an illustration of
        # data structure, not a biochemical or mechanical model.
        tau = rng.uniform(3,7)
        memory, vel = np.zeros((len(times),2)), np.zeros((len(times),2))
        vel[0] = .8*np.array([np.cos(baseline_heading),np.sin(baseline_heading)])
        for k in range(1,len(times)):
            drive = active_area[k]/max(rate*9,1)
            memory[k] = memory[k-1] + dt*(drive-memory[k-1])/tau
            vel[k] = vel[k-1] + dt*(-vel[k-1]/.85 + .80*memory[k]) + rng.normal(0,.48,2)*np.sqrt(dt)
        pos = np.cumsum(vel*dt,axis=0)
        speed = np.linalg.norm(vel,axis=1)
        noise = np.zeros(len(times))
        for k in range(1,len(times)):
            noise[k] = .90*noise[k-1]+rng.normal(0,.07)
        anisotropy = np.clip(.25+.24*speed+noise+rng.normal(0,.07), .05,.9)
        axis_alignment = np.clip(rng.normal(0,.32,len(times))+noise, -np.pi/2,np.pi/2)
        tracks[ident] = dict(time=times,velocity=vel,position=pos,anisotropy=anisotropy,
            axis_alignment=axis_alignment,total_area=total_area)
        complete = [r for r,p in local_events if not r['left_censored'] and not r['right_censored']]
        onsets = sum(0<=r['start_min']<horizon for r,p in local_events)
        cells.append(dict(cell=ident,size_group=group,perimeter_um=float(perimeter),
            area_um2=float(cell_area),observation_min=horizon,onsets=onsets,
            events_per_10min=onsets/horizon*10,
            events_per_10min_per_100um=onsets/horizon*10/perimeter*100,
            median_complete_duration_min=float(np.median([r['observed_duration_min'] for r in complete])),
            median_peak_area_um2=float(np.median([r['observed_peak_area_um2'] for r in complete])),
            median_width_um=float(np.median([r['footprint_width_um'] for r in complete])),
            censored_events=sum(r['left_censored'] or r['right_censored'] for r,p in local_events)))
        for cutoff in [15.,25.,35.,45.]:
            k = int(round(cutoff/dt))
            v0 = vel[k-3:k+1].mean(axis=0)  # trailing 36 s, ending at cutoff
            heading = np.arctan2(v0[1],v0[0])
            forward = np.array([np.cos(heading),np.sin(heading)])
            lateral = np.array([-np.sin(heading),np.cos(heading)])
            vprev = vel[k-6:k-3].mean(axis=0)
            turn = wrapped(heading-np.arctan2(vprev[1],vprev[0]))/.6
            hist_start,hist_end = cutoff-10,cutoff-1
            mask = (times>=hist_start)&(times<hist_end)
            hist_count,hist_duration,hist_q,hist_signed = 0,0.,0.,0.
            directional = np.zeros(4)
            for r,pulse in local_events:
                if hist_start<=r['start_min']<hist_end:
                    hist_count += 1
                hist_duration += max(0,min(r['end_min'],hist_end)-max(r['start_min'],hist_start))
                q = float(pulse[mask].sum()*dt)
                rel = wrapped(r['direction_rad']-heading)
                hist_q += q
                hist_signed += q*np.sin(rel)
                directional[int(np.floor((rel+np.pi/4)%(2*np.pi)/(np.pi/2)))] += q
            imbalance = hist_signed/max(hist_q,1e-9)
            future_start,future_end = cutoff+1,cutoff+3
            fa,fb = [int(round(z/dt)) for z in [future_start,future_end]]
            delta = pos[fb]-pos[fa]
            current = dict(cell=ident,size_group=group,cutoff_min=cutoff,
                current_speed_um_min=float(np.linalg.norm(v0)),
                current_heading_rad=float(heading), current_anisotropy=float(anisotropy[k]),
                current_axis_alignment_rad=float(axis_alignment[k]),
                current_turn_rate_rad_min=float(turn), current_area_um2=float(cell_area),
                current_perimeter_um=float(perimeter), current_edge_area_um2=float(total_area[k]),
                history_start_min=hist_start,history_end_min=hist_end,
                past_onsets=hist_count,past_event_time_min=float(hist_duration),
                past_area_time_um2_min=float(hist_q),past_lateral_imbalance=float(imbalance),
                future_start_min=future_start,future_end_min=future_end,
                future_forward_displacement_um=float(delta@forward),
                future_lateral_displacement_um=float(delta@lateral),
                future_anisotropy=float(anisotropy[fb]))
            predictions.append(current)
            if cutoff==25:
                snapshots.append(current)
            if ident=='S07' and cutoff==25:
                history_example = (current,directional)


def write_csv(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


write_csv(DATA/'synthetic_events.csv',events)
write_csv(DATA/'synthetic_cell_summaries.csv',cells)
write_csv(DATA/'synthetic_state_history_future.csv',predictions)
np.savez_compressed(DATA/'synthetic_tracks.npz',**{
    f'{ident}_{k}':v for ident,track in tracks.items() for k,v in track.items()})

fig=plt.figure(figsize=(18.4,16.2),facecolor='white')
page_text=[]


def txt(x,y,s,**kw):
    a=fig.text(x,y,s,**kw);page_text.append(a);return a


def box(x,y,w,h,face='white'):
    fig.add_artist(FancyBboxPatch((x,y),w,h,transform=fig.transFigure,
        boxstyle='round,pad=0.007,rounding_size=0.008',facecolor=face,
        edgecolor=RULE,lw=.8,zorder=-4))


def header(x,y,letter,title,subtitle=None):
    txt(x,y,letter,fontsize=17,fontweight='bold',va='top')
    txt(x+.024,y,title,fontsize=14,fontweight='bold',va='top')
    if subtitle:
        txt(x+.024,y-.022,subtitle,fontsize=10.2,color=MUTED,va='top')


def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(.7)
    ax.tick_params(length=3,width=.7)
    ax.grid(axis='y',color=RULE,lw=.65,alpha=.7)
    ax.set_axisbelow(True)


callouts=[]


def badge(x,y,code,fontsize=8.0):
    return txt(x,y,code,ha='center',va='center',fontsize=fontsize,
        fontweight='bold',color=PURPLE,
        bbox=dict(boxstyle='circle,pad=.20',facecolor='#f0eaf5',
                  edgecolor='#cfc0dc',linewidth=.55))


def callout(x,y,direction,label,codes=(),fontsize=9.4,panel=''):
    lo,hi=(y-.005,y+.005) if direction=='up' else (y+.005,y-.005)
    fig.add_artist(FancyArrowPatch((x,lo),(x,hi),transform=fig.transFigure,
        arrowstyle='-|>',mutation_scale=8,linewidth=1.1,color=PURPLE,
        shrinkA=0,shrinkB=0))
    a=txt(x+.009,y,label,fontsize=fontsize,va='center')
    width=a.get_window_extent(fig.canvas.get_renderer()).width/fig.bbox.width
    for i,code in enumerate(codes):badge(x+.009+width+.010+i*.014,y,code)
    callouts.append(dict(panel=panel,direction_in_large_cells=direction,
        measurement=label,mechanism_codes=list(codes),status='conditional hypothesis'))


txt(.038,.972,'Protrusion statistics, cell state and behavioral memory',fontsize=24,fontweight='bold',va='top')
txt(.038,.939,'Raw events → current observable state → past activity that might predict the future',fontsize=14,color=MUTED,va='top')
txt(.038,.912,'NOISY SYNTHETIC DATA  /  ARROW INSETS: CANDIDATE CHANGES IN LARGE CELLS',
    fontsize=10.2,color=PURPLE,fontweight='bold',va='top')
for x,col,label in [(.793,TEAL,'Small'),(.875,ORANGE,'Large')]:
    fig.add_artist(Line2D([x,x+.022],[.943,.943],transform=fig.transFigure,color=col,lw=3))
    txt(x+.028,.943,label,fontsize=11,va='center')

# A: actual observable record, with event duration drawn as intervals.
box(.038,.755,.925,.135,PALE)
header(.049,.879,'A','Build an event record before assigning a mechanism',
    'Record onset/end, width, gained area and direction; flag incomplete events')
source=ROOT/'data/random2d_assay_mockup_20261009/ht1080_outline_examples.json'
raw=np.array(json.loads(source.read_text())['records'][0]['outline_xy_um'][2],dtype=float)
raw[:,1]*=-1;raw-=raw.mean(axis=0)
_,e=np.linalg.eigh(np.cov(raw.T));raw=raw@e[:,::-1]
if np.std(raw[raw[:,0]>0,1])<np.std(raw[raw[:,0]<0,1]):raw[:,0]*=-1
raw/=np.ptp(raw[:,0]);raw*=2
a=np.arctan2(raw[:,1],raw[:,0]);d=wrapped(a-1.25)
current=raw.copy();current[:,1]+=.35*np.exp(-.5*(d/.3)**2)
ax=fig.add_axes([.048,.769,.190,.077]);ax.set_aspect('equal');ax.axis('off')
ax.plot(*np.vstack([raw,raw[0]]).T,color='#91a1ab',ls=(0,(3,2)),lw=1.3)
ax.fill(*current.T,color='#e2e9ee',zorder=1)
ax.plot(*np.vstack([current,current[0]]).T,color=INK,lw=1.3)
sel=np.where(abs(d)<.46)[0]
ax.fill(*np.vstack([current[sel],raw[sel][::-1]]).T,color=PURPLE,alpha=.35,zorder=3)
ax.plot(*current[sel].T,color=PURPLE,lw=3)
p=current[np.argmin(abs(d))]
ax.annotate('width + area',xy=p,xytext=(-1.4,1.1),fontsize=9,color=PURPLE,
    arrowprops=dict(arrowstyle='-',color=PURPLE,lw=.8))
ax.annotate('',xy=(.80,-.73),xytext=(-.05,-.73),
    arrowprops=dict(arrowstyle='-|>',color=INK,lw=1.25))
ax.text(-.12,-.73,'heading',fontsize=8.5,color=MUTED,va='center',ha='right')
ax.set(xlim=(-1.5,1.3),ylim=(-.9,1.3))
txt(.246,.835,'One row per event',fontsize=11,fontweight='bold',va='top')
txt(.246,.815,'Cell ID · start/end · size\nDirection · split/merge flags\nCensoring · detection quality',fontsize=10,color=MUTED,va='top',linespacing=1.55)
ax=fig.add_axes([.463,.782,.470,.060])
for yi,ident in enumerate(['S07','S14','L07','L14']):
    group='small' if ident[0]=='S' else 'large'
    visible=[]
    lane_ends=[]
    for r in events:
        if r['cell']==ident and r['observed_end_min']>5 and r['observed_start_min']<20:
            start,end=max(5,r['observed_start_min']),min(20,r['observed_end_min'])
            free=[i for i,e in enumerate(lane_ends) if e<=start]
            lane=free[0] if free else len(lane_ends)
            if lane==len(lane_ends):lane_ends.append(end)
            else:lane_ends[lane]=end
            visible.append((r,start,end,lane))
    for r,start,end,lane in visible:
        ypos=yi+(lane-(len(lane_ends)-1)/2)*(.55/max(1,len(lane_ends)))
        ax.plot([start,end],[ypos,ypos],color=COLORS[group],alpha=.82,
            lw=1.0+min(3,r['observed_peak_area_um2']/10),solid_capstyle='butt')
        if r['start_min']>=5:
            ax.plot([start,start],[ypos-.035,ypos+.035],color=COLORS[group],lw=.7)
ax.set(xlim=(5,20),ylim=(-.5,3.5),yticks=range(4),yticklabels=['S07','S14','L07','L14'],
    xticks=[5,10,15,20],xlabel='Time (min)')
ax.invert_yaxis();ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
txt(.934,.850,'Bar length: duration · Thickness: gained area',fontsize=9,color=MUTED,va='bottom',ha='right')

# B: broad, per-cell distributions, no constructed mechanistic sign.
box(.038,.492,.925,.242)
header(.049,.722,'B','Start with number, duration and size',
    'Solid dots: per-cell summaries. Faint dots: individual events. Keep all attempts, including those that disappear.')
metrics=[('events_per_10min','Number','Events / 10 min'),
         ('events_per_10min_per_100um','Number / perimeter','Events / 10 min / 100 µm'),
         ('median_complete_duration_min','Duration','Event duration (min)'),
         ('median_peak_area_um2','Size','Peak gained area (µm²)')]
jitter_rng=np.random.default_rng(47)
for j,(key,title,ylabel) in enumerate(metrics):
    ax=fig.add_axes([.080+j*.228,.554,.170,.119]);clean(ax)
    for gi,group in enumerate(['small','large']):
        if j in [2,3]:
            event_key='observed_duration_min' if j==2 else 'observed_peak_area_um2'
            ev=np.array([r[event_key] for r in events if r['size_group']==group
                and not r['left_censored'] and not r['right_censored']])
            ax.scatter(gi+jitter_rng.uniform(-.23,.23,len(ev)),ev,s=5,
                color=COLORS[group],alpha=.09,linewidths=0,zorder=2)
        vals=np.array([r[key] for r in cells if r['size_group']==group])
        x=gi+jitter_rng.uniform(-.17,.17,len(vals))
        ax.scatter(x,vals,s=24,color=COLORS[group],alpha=.65,edgecolors='white',linewidths=.3,zorder=3)
        ax.plot([gi-.20,gi+.20],[np.median(vals)]*2,color=INK,lw=1.8,zorder=5)
    ax.set(xlim=(-.45,1.45),xticks=[0,1],xticklabels=['Small','Large'],ylabel=ylabel)
    ax.set_title(title,fontsize=12,fontweight='bold',loc='left',pad=5)
    if j in [2,3]:
        ax.set_yscale('log')
        ax.tick_params(axis='y',which='minor',length=0)
        ax.set_yticks([.2,1,5,20] if j==2 else [1,10,100])
        ax.set_yticklabels(['0.2','1','5','20'] if j==2 else ['1','10','100'])
callout(.062,.522,'up','Perimeter-only count',panel='B-number')
txt(.071,.506,'Geometry / exposure alone',fontsize=8.7,color=MUTED,va='center')
callout(.290,.522,'down','Detectable starts / perimeter',['P'],panel='B-rate')
callout(.290,.506,'down','Persistent off-axis events',['P','R'],panel='B-rate')
callout(.518,.522,'down','Aborted-event lifetime',['R'],panel='B-duration')
callout(.518,.506,'up','Activity before takeover',['W','T','M'],panel='B-duration')
callout(.746,.522,'down','Early gained area',['P'],panel='B-size')
callout(.746,.506,'up','Spread before retained turn',['W'],panel='B-size')

# C: deliberately continuous, observable current state, not inferred polarity.
for x in [.038,.354,.670]:box(x,.155,.293,.322)
header(.049,.466,'C','Define state and its occupancy',
    'Motion, shape and angular edge activity')
ax=fig.add_axes([.087,.345,.221,.078]);clean(ax)
for group in ['small','large']:
    rows=[r for r in snapshots if r['size_group']==group]
    ax.scatter([r['current_anisotropy'] for r in rows],[r['current_speed_um_min'] for r in rows],
        s=22,color=COLORS[group],alpha=.65,edgecolors='white',linewidths=.25)
ax.set(xlabel='Contour anisotropy',ylabel='Current speed (µm/min)',xlim=(0,1))
txt(.052,.310,'Also retain heading, turning rate and axis alignment.',fontsize=9.6,va='top')
txt(.052,.292,'Keep state continuous; avoid an assumed “polarity” class.',fontsize=9.3,color=MUTED,va='top')

# D: history summaries have distinct units; never one arbitrary combined score.
header(.365,.466,'D','Summarize the preceding history',
    'Count, event-time and size-time remain separate')
example,directional=history_example
ax=fig.add_axes([.375,.380,.245,.040]);ax.set(xlim=(-11,4),ylim=(0,1));ax.axis('off')
for lo,hi,y,label,color in [(-10,-1,.65,'past events\n−10 to −1 min',PURPLE),
                          (-.6,0,.65,'state\n−0.6 to 0',INK),
                          (1,3,.65,'future\n+1 to +3',TEAL)]:
    ax.plot([lo,hi],[y,y],color=color,lw=7,solid_capstyle='butt')
    ax.text((lo+hi)/2,.09,label,ha='center',va='center',fontsize=7.8,color=color,linespacing=1.3)
ax.axvline(0,color=INK,lw=.8,ls=':')
ax.text(0,1.05,'t = 0',ha='center',fontsize=8.4,color=INK)
txt(.371,.363,f"Onsets: {example['past_onsets']}\nSummed active time: {example['past_event_time_min']:.1f} event·min\nArea-time: {example['past_area_time_um2_min']:.0f} µm²·min",
    fontsize=10.0,va='top',linespacing=1.5)
txt(.371,.301,'Keep direction + recency; compare several history windows.',fontsize=9.3,color=MUTED,va='top')
txt(.371,.285,'Opposing protrusions need separate channels.',fontsize=9.3,color=MUTED,va='top')

# E: noisy future outcomes from exactly the same event/track dataset.
header(.681,.466,'E','Ask what happens next',
    'Does history add information beyond current state?')
ps=np.array([r['current_speed_um_min'] for r in predictions])
pa=np.array([r['current_anisotropy'] for r in predictions])
pt=np.abs([r['current_turn_rate_rad_min'] for r in predictions])
selected=(ps>=np.quantile(ps,.15))&(ps<=np.quantile(ps,.85))&\
         (pa>=np.quantile(pa,.10))&(pa<=np.quantile(pa,.90))&(pt<=np.quantile(pt,.80))
ax=fig.add_axes([.721,.345,.218,.078]);clean(ax)
for group in ['small','large']:
    rows=[r for i,r in enumerate(predictions) if selected[i] and r['size_group']==group]
    ax.scatter([r['past_lateral_imbalance'] for r in rows],
        [r['future_lateral_displacement_um'] for r in rows],
        s=17,color=COLORS[group],alpha=.54,edgecolors='white',linewidths=.25)
ax.axhline(0,color=RULE,lw=1)
ax.set(xlabel='Past lateral area-time bias',ylabel='Future lateral displacement (µm)',xlim=(-1.05,1.05))
txt(.684,.310,'Example state stratum; each dot is a forecast window.',fontsize=9.5,va='top')
txt(.684,.292,'Compare state-only vs state + history on held-out cells.',fontsize=9.3,color=MUTED,va='top')

# Directional hypothesis insets sit outside the synthetic data axes. They can
# refer to additional observables, without assigning signs to speed/anisotropy
# or conflating heading return with signed lateral displacement.
for x in [.049,.365,.681]:box(x,.174,.271,.094,'#faf7fc')
txt(.057,.261,'Large-cell state discrepancies to test',fontsize=9.2,fontweight='bold',color=PURPLE,va='center')
callout(.060,.241,'up','Old-heading dwell time',['R','M'],panel='C')
callout(.060,.216,'up','Competing-front occupancy',['W','T'],panel='C')
callout(.060,.191,'up','Unretained deflections',['M'],panel='C')
txt(.373,.261,'For the same retained-turn probability',fontsize=9.2,fontweight='bold',color=PURPLE,va='center')
callout(.376,.241,'up','Directed area-time needed',['R','W','M'],panel='D')
callout(.376,.216,'up','Useful history horizon?',['M'],panel='D')
callout(.376,.191,'down','Load if mechanics facilitates',['X'],panel='D')
txt(.689,.261,'At matched state and event exposure',fontsize=9.2,fontweight='bold',color=PURPLE,va='center')
callout(.692,.245,'down','Retained turns / attempt',['R','W','M'],panel='E')
callout(.692,.225,'up','Return to old heading',['R','M'],panel='E')
callout(.692,.205,'up','Redirection delay',['T','M'],panel='E')
callout(.692,.185,'down','Delay if mechanics facilitates',['X'],panel='E')

# Compact lookup, with no comparison table or data rows.
txt(.038,.134,'Mechanism key',fontsize=11.5,fontweight='bold',va='center')
txt(.161,.134,'Arrows are conditional hypotheses; competing mechanisms can push the same measurement in opposite directions.',
    fontsize=9.8,color=MUTED,va='center')
mechanisms=[
    ('P','Depleted available\npool','Goryachev · Fig. 3'),
    ('R','Stronger incumbent\nrepair','Goryachev · Fig. 4'),
    ('W','Wider polarized\ndomain','Wave-pinning · Figs. 2, 6'),
    ('T','Redistribution /\nseparation burden','Figs. 1, 2, 5'),
    ('M','Slow attachment\nmemory','deBelly · Figs. 8–11'),
    ('X','Facilitating flow /\nextension feedback','deBelly / Spring · Figs. 1, 9–11'),
]
for i,(code,label,refs) in enumerate(mechanisms):
    x=.048+i*.156
    badge(x,.098,code,fontsize=10)
    txt(x+.014,.111,label,fontsize=10,fontweight='bold',va='top',linespacing=1.45)
    txt(x+.014,.069,refs,fontsize=8.0,color=MUTED,va='center')

txt(.038,.031,'All points and event records are synthetic. History dependence is built into the illustrative generator; the figure provides no evidence for biological memory.',
    fontsize=9.7,color=MUTED,va='center')
txt(.038,.014,'HT1080 outline template only · Bosgraaf & Van Haastert (2009), doi:10.1371/journal.pone.0005253 · Mechanism figure numbers refer to length_mechanisms.pdf.',
    fontsize=9.0,color=MUTED,va='center')

CAPTION='''Event-first design for a reporter-free 2D random-migration assay, with directional large-cell hypotheses. Numeric data are the same stochastic illustration as event-first v2: no new mechanistic size shifts have been imposed on the dot clouds. Purple arrow insets propose additional size discrepancies to test and identify contributing mechanisms through six letter badges. Different arrows are conditional, sometimes opposing hypotheses, not a single joint predicted phenotype. The previous v1 and v2 figures are retained. No comparison table is included.

A. Detect connected regions of protrusive boundary motion over time. Store onset/end, footprint width, gained area and extension rate, event location and direction, cell identity, uncertainty and split/merge flags. A practical event definition must be calibrated against contour noise and manually reviewed movies, with the same physical detection rules across sizes. Store laboratory-frame boundary displacement and body-registered deformation separately: steady front advance and shape change are related but different observations. Registration does not identify permanent membrane material. Event merging and splitting require a declared lineage rule so segmentation changes do not become new biological onsets. Raster bars depict active intervals; thickness depicts observed event peak gained area. The shown raster is an excerpt from the same dataset as all other plots, not a model-specific example.

B. Begin with the raw count, duration and size distributions. Count onsets per observation time, with a separate normalization by measured perimeter-time exposure; do not assume perimeter normalization is the biological answer. Retain the number of simultaneously active regions and inter-event waiting times. Duration means all detected lifetimes, not only protrusions that lead to turns. Track left/right censoring; use survival analysis for full duration distributions and report complete-event summaries as descriptive quantities. For duration and size, faint dots show individual complete events, solid dots show per-cell medians (24 cells in each display bin), and dark horizontal bars show the median of those cell summaries. Duration/area axes use log spacing to retain the broad event tails. Censored events remain explicitly recorded in the event table. Size includes physical width, gained area, extension depth, peak size and time course; also retain fractional perimeter/area measurements without replacing the physical ones. For the two count plots, each solid dot is one cell. These broad distributions intentionally carry no prescribed mechanism-specific length shifts. Synthetic event rate per perimeter, event-duration, size and memory-time distributions are identical across size bins; a larger perimeter can yield more events per cell without altered local event dynamics. Use cell-level uncertainty or hierarchical models rather than treating every event as an independent biological replicate.

C. Start with a continuous observable state S(t) rather than an assumed biochemical-polarity label. The minimal practical state includes short-window velocity/speed/heading, current turning rate, contour area/perimeter, shape anisotropy, the shape-axis angle relative to heading, and current total positive edge activity plus its angular distribution. The scatter shows two coordinates only: anisotropy and speed. Anisotropy and an axis cannot supply a front/back sign by themselves. A weak or undefined heading at low speed is retained as a state/uncertainty flag; do not silently discard paused cells or orient events to noisy velocities. Candidate state descriptions should be checked for reproducibility and short-term predictive utility. Expand this state if a simpler one leaves obvious ongoing turning or current edge competition unaccounted for. Size is a continuous explanatory covariate, not matched away between display bins. Dimensionless shape descriptors can be matched across cell sizes. The synthetic anisotropy is a generated track descriptor, not obtained from a simulation of the illustrated HT1080 contour.

D. Keep past event count, summed active duration and size exposure distinct. Summed active duration is the time integral of the number of simultaneously active events; overlap can exceed elapsed clock time. Area-time is the time integral of event-associated gained area. Do not multiply by duration again after integrating area-time. Also retain physical width-time, early growth, event order and inter-event intervals when useful. Compute these summaries by direction (for example front/left/back/right relative to current reliable heading), along with total load, so strong opposing events are not confused with no activity. Compare rectangular history windows and recency-weighted summaries; select window lengths or decay times using training cells only. The displayed values describe one actual synthetic window (cell S07, cutoff 25 min). The example history window is t−10 to t−1 min; current velocity uses t−0.6 to t; future displacement uses t+1 to t+3 min. These windows are illustrative choices, not known biological timescales. Clip history contributions at the history cutoff: ongoing events contribute only their observed activity and elapsed duration, never their eventual lifetime or future maximum size. Completed-event duration statistics in B and accrued predictor duration in D are different quantities.

E. Ask whether these historical summaries improve prediction of future motion/shape after conditioning on current state and continuous cell size. Outcomes can be future forward/lateral displacement relative to the current heading, speed, shape anisotropy, or subsequent event count/duration/size. Future displacement is meaningful when a future cell stops and avoids defining all outcomes only for successful turns. The noisy scatter illustrates one relationship between preceding lateral area-time imbalance and future lateral displacement, in a broad common current-state stratum (central speed and anisotropy ranges, excluding the strongest current turn rates). This stratum is illustrative matching, not a complete adjustment for state. The imbalance is summed area-time times sine of event angle relative to current heading, divided by total area-time; the full predictor set must retain total area-time and opposite-direction channels, not only this signed scalar. Cases with unreliable current heading need nondirectional outcomes or explicit uncertainty handling. No regression line, fitted accuracy, significance claim or biological inference is plotted. The synthetic motility generator contains history dependence by construction and uses identical memory-time distributions across size groups. The intended empirical test compares current-state-only forecasts with current-state-plus-history forecasts, holding out entire cells and contiguous trajectory blocks; tune history windows only in training data. A positive result would show useful information beyond the chosen observable state, without uniquely identifying an intracellular species or physical memory mechanism.

Directional insets and mechanism key. P is reduced available inactive concentration in the characterized Goryachev regime (Fig. 3); R is stronger immediate incumbent restoration per amount lost (Goryachev, Fig. 4); W is growing polarized-domain width and retained replacement requiring domain expansion (Wave-pinning, Figs. 2, 6); T is greater redistributed material/target distance and the measured extra-separation and fixed-gap penalties (Figs. 1, 2, 5); M is slower membrane–cortex attachment/mechanical rewriting compared with early Rac response (corrected deBelly, Figs. 8–11); X denotes regimes where increased flow or extension feedback facilitates reversal/redistribution (deBelly/Spring, Figs. 1, 9–11). X is conditional: extension or force transmission is not universally facilitating, and the report does not establish a unique sign for every turnover.

B inset interpretation. A larger boundary can increase onsets per cell even when the local event rate is unchanged; this is an exposure effect, not one of the intracellular mechanism badges. P may reduce detectable onsets per perimeter and early gained area at matched initial footprint/age. P/R may suppress persistent off-axis competitors at comparable present state; R here refers to competitor rejection/maintenance, not a directly established molecular nucleation-rate change. R could shorten aborted-competitor lifetimes. W/T/M instead motivate more sustained activity before takeover, a separate outcome-conditioned clock that must not be equated with the median duration of all events. W motivates greater spatial spread before retained redirection, not a universal increase in every protrusion's maximum size. The arrows sit beneath raw-statistic displays but refer to these explicitly distinguished event classes and conditional readouts.

C inset interpretation. Current speed and contour anisotropy receive no imposed directional predictions. Instead, proposed size differences concern occupancies and transitions defined from the measured state: old-heading dwell time may increase under incumbent restoration or slow attachment memory (R/M); the fraction of time with spatially separated concurrently active fronts may increase when replacement/cleanup takes longer (W/T); brief heading deflections may fail to become retained turns more often when visible initiation runs ahead of mechanical rewriting (M). Define multi-front occupancy from at least two active contour sectors separated by a frozen angular/distance criterion, and define old-heading dwell and unretained deflections with a fixed speed floor, angle bands and hold time. These are additional empirically testable state/transition summaries, not direct implications of the two-coordinate scatter alone.

D inset interpretation. At a comparable retained-turn probability, larger cells may require more directed protrusive area-time (R/W/M), reflecting stronger incumbent resistance, more domain spread or a mechanical imprint to rewrite. This is a conditional outcome relation, not an established universal scalar threshold; counts, active duration, area-time and direction remain separate candidate predictors. M motivates testing a longer useful look-back window, marked with a question mark because the report's physical redistribution clock does not establish an empirically optimal statistical forecast window. Under facilitating flow/extension regimes, X could oppose increased exposure requirements. Tune every history window on training cells, and test size interactions on held-out cells.

E inset interpretation. At matched measured state and fixed physical event exposure, R/W/M could reduce retained turns per attempt; R/M could increase returns toward the old heading; T/M could increase stable-redirection lag. Facilitating mechanics X could reduce that lag, providing an explicit opposing prediction. Heading return, turn retention and response lag are proposed additional future readouts; they are not interchangeable with a positive/negative shift on the plotted signed lateral-displacement axis. Use equal-angle control heading bands to distinguish old-heading memory from generally increased wandering. The direction and magnitude of every large-small difference must be measured; size trends can be nonmonotonic.

Reference: Bosgraaf L, Van Haastert PJM (2009), The Ordered Extension of Pseudopodia by Amoeboid Cells in the Absence of External Cues, PLOS ONE 4(4):e5253, https://doi.org/10.1371/journal.pone.0005253.
'''
(FIG/'CAPTION.txt').write_text(CAPTION)
(DATA/'directional_hypotheses.json').write_text(json.dumps(callouts,indent=2)+'\n')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for artist in page_text:
    b=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    assert b.x0>=0 and b.y0>=0 and b.x1<=1 and b.y1<=1,artist.get_text()
outputs={}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for ext in ['pdf','svg','png']:
    p=FIG/f'random2d_measurements.{ext}';fig.savefig(p,dpi=210);outputs[p.name]=sha(p)
plt.close(fig)
manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
    revision='directional_v3',seed=SEED,synthetic=True,mechanism_specific_size_shifts=False,
    hypothetical_directional_callouts=len(callouts),
    history_dependence_imposed_for_illustration=True,cells=len(cells),events=len(events),
    forecast_windows=len(predictions),forecast_scatter_windows=int(selected.sum()),
    output_sha256=outputs,data_sha256={p.name:sha(p) for p in DATA.glob('*') if p.is_file()},
    source_sha256={str(Path(__file__).relative_to(REPO)):sha(__file__),
        str(source.relative_to(REPO)):sha(source),
        'polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf':
            sha(ROOT/'figures/length_mechanisms_20261008/length_mechanisms.pdf')},
    caption_sha256=sha(FIG/'CAPTION.txt'),
    reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/random2d_directional_hypotheses_20261009.py')
(FIG/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(FIG/'README.txt').write_text('Directional v3 is the current assay-design figure.\n\nOpen random2d_measurements.pdf (vector), .svg (editable), or .png (preview).\nPurple arrow insets propose size discrepancies in panels B–E, linked to a six-item mechanism key.\nAll numerical values remain the same synthetic illustration as v2. The arrows are hypotheses, not fits.\nCAPTION.txt defines every arrow and observable. Data and machine-readable hypotheses are in ../../../data/random2d_assay_mockup_20261009/directional_v3/.\n\nReproduce from repository root:\n'+manifest['reproduction']+'\n')
print(json.dumps({k:manifest[k] for k in ['revision','cells','events','forecast_windows','forecast_scatter_windows']}))
