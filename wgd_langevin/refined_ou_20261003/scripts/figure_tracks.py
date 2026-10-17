"""Shared representative-track selection and plotting for the refined OU fit."""
from pathlib import Path
import csv, hashlib, json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[3]
FIT_DIR = ROOT / 'wgd_langevin/refined_ou_20261003/results'
OUT = ROOT / 'wgd_langevin/refined_ou_20261003/figures'
OUT.mkdir(parents=True, exist_ok=True)
NAVY = '#19394d'
MINUS = '#168b87'
PLUS = '#d87529'
GRAY = '#6d7c86'
CAD, DURATION, LAG = 1/6, 12., 11.1111111111
STUDY_INDEX = 0  # Xuan2018, whose fitted measures include sampled speed and full population MSD curves.

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42})

# Use the same exact integrated-OU displacement function as the fit pipeline.
sys.path.insert(0, str(ROOT / 'wgd_langevin/persistence_only'))
from campaign import displacement_variance

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def track_metrics(track):
    step=np.diff(track,axis=0)
    path=float(np.linalg.norm(step,axis=1).sum())
    endpoint=float(np.linalg.norm(track[-1]-track[0]))
    return {'path_length_um':path,'mean_sampled_speed_um_h':path/DURATION,
            'endpoint_displacement_um':endpoint,'directness':endpoint/path}

def choose_representative(tracks):
    """Closest robustly standardized candidate to group medians.

    Score speed, directness, and endpoint displacement; path length is omitted
    because it is speed multiplied by the fixed 12 h duration.
    """
    metrics=[track_metrics(t) for t in tracks]
    keys=['mean_sampled_speed_um_h','directness','endpoint_displacement_um']
    a=np.array([[m[k] for k in keys] for m in metrics])
    med=np.median(a,axis=0)
    q25,q75=np.quantile(a,[.25,.75],axis=0)
    scale=np.maximum((q75-q25)/1.349,np.maximum(np.abs(med)*1e-8,1e-10))
    score=np.sum(((a-med)/scale)**2,axis=1)
    i=int(np.argmin(score))
    return tracks[i],metrics[i],{'index':i,'score':float(score[i]),
        'group_medians':dict(zip(keys,map(float,med))),
        'robust_scales_iqr_over_1.349':dict(zip(keys,map(float,scale))),
        'scored_metrics':keys,'path_length_omitted_as_redundant_with_speed':True}

def msd_expected(t,tau,eta):
    # Exact two-dimensional expectation; eta is the Wiener coefficient.
    return 2*eta**2*displacement_variance(np.asarray(t),np.array([1/tau]),np.array([tau/2]))

def arrow(ax,a,b,color=GRAY,lw=1.5,style='-|>'):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,
                                 color=color,lw=lw,shrinkA=0,shrinkB=0))

def axes_track(ax,track,color,linewidth=1.55,markers=True):
    ax.plot(track[:,0],track[:,1],color=color,lw=linewidth,zorder=2)
    ax.scatter(*track[0],s=24,marker='o',facecolor='white',edgecolor=color,lw=1.2,zorder=4)
    ax.scatter(*track[-1],s=38,marker='*',facecolor=color,edgecolor='white',lw=.4,zorder=5)

def fmt_metrics(m):
    return (f"Path length: {m['path_length_um']:.0f} µm\n"
            f"Mean sampled speed: {m['mean_sampled_speed_um_h']:.1f} µm/h\n"
            f"End displacement: {m['endpoint_displacement_um']:.0f} µm\n"
            f"Directness: {m['directness']:.2f}")

def style_axis(ax,lim):
    ax.set_aspect('equal',adjustable='box');ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
    ax.set_xlabel('Position x (µm)');ax.set_ylabel('Position y (µm)')
    ax.grid(alpha=.15,lw=.5)

