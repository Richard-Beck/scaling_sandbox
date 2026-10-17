#!/usr/bin/env python3
"""Moving laboratory-coordinate kymographs of the fixed-substrate clutch trials."""
import json
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from polarity1d.contract import Trajectory
from build_debelly_figure2_kymograph_supplement import _edges

CAMPAIGN = ROOT/'results/campaigns/debelly_fixed_substrate_L21'
OUTPUT = ROOT/'results/figures'
FIELDS = (
    ('rac', 'Rac', 'viridis', False),
    ('rho', 'Rho', 'magma', False),
    ('tension', 'Membrane tension', 'coolwarm', True),
    ('mca_bound', 'Bound MCA', 'cividis', False),
    ('local_strain', 'Local length change (%)', 'coolwarm', True),
    ('velocity', 'Cortical velocity: cell frame (µm/s)', 'RdBu_r', True),
    ('cortical_velocity_lab', 'Cortical velocity: lab frame (µm/s)', 'RdBu_r', True),
    ('clutch_bound', 'Bound clutch fraction', 'viridis', False),
    ('clutch_extension', 'Mean bound extension (µm)', 'RdBu_r', True),
    ('clutch_traction_lab_density', 'Clutch traction on substrate (pN/µm)', 'RdBu_r', True),
)


def prepare(density):
    tr = Trajectory.load(CAMPAIGN/f'density_{density}/trajectory.npz')
    f = tr.fields
    edge_u = np.column_stack((f['left_boundary_displacement'][:,0],
        .5*(f['displacement'][:,:-1]+f['displacement'][:,1:]),
        f['right_boundary_displacement'][:,0]))
    reference_edges = _edges(tr.x)
    edges = reference_edges[None,:]+edge_u+f['cell_translation'][:,:1]
    f['local_strain'] = 100*(np.diff(edges,axis=1)/np.diff(reference_edges)-1)
    length = edges[:,-1]-edges[:,0]
    np.testing.assert_allclose(length,21+f['boundary_separation'][:,0],atol=1e-10)
    if np.any(np.diff(edges,axis=1)<=0):
        raise ValueError('Folded display geometry')
    return tr, edges, length


def panel(ax, tr, edges, field, cmap, norm):
    t = tr.time-200
    indices = np.unique(np.r_[np.flatnonzero(t<=800),np.arange(0,len(t),4),len(t)-1])
    te = _edges(t[indices])
    te[0],te[-1] = t[0],t[-1]
    mapped = np.column_stack([np.interp(te,t,edges[:,j]) for j in range(edges.shape[1])])
    mesh = ax.pcolormesh(np.broadcast_to(te,mapped.T.shape),mapped.T,
                        tr.fields[field][indices].T,cmap=cmap,norm=norm,
                        shading='flat',rasterized=True)
    ax.plot(t,edges[:,0],color='black',lw=.65)
    ax.plot(t,edges[:,-1],color='black',lw=.65)
    ax.axvline(0,color='gray',ls='--',lw=.7)
    ax.axvline(160,color='gray',ls='--',lw=.7)
    ax.set_xscale('symlog',linthresh=200)
    ax.set_xlim(t[0],t[-1])
    ax.set_ylim(edges.min()-.5,edges.max()+.5)
    ax.set_xlabel('Time from cue onset (s; symlog)')
    return mesh


def save(fig, stem):
    for ext in ('png','pdf'):
        fig.savefig(OUTPUT/f'{stem}.{ext}',dpi=220)
    plt.close(fig)


def main():
    global CAMPAIGN
    parser=argparse.ArgumentParser()
    parser.add_argument('--local-rac', action='store_true')
    parser.add_argument('--uncentered', action='store_true')
    parser.add_argument('--passive-rear', action='store_true')
    args=parser.parse_args()
    if args.passive_rear or args.local_rac:
        args.uncentered=True
    suffix='_uncentered' if args.uncentered else ''
    annotation=' · NO CENTERING' if args.uncentered else ''
    if args.uncentered:
        CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21_uncentered'
    if args.passive_rear:
        CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21_passive_rear'
        suffix='_passive_rear'
        annotation=' · PASSIVE REAR (OPTION 2)'
    if args.local_rac:
        CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21_local_rac'
        suffix='_local_rac'
        annotation=' · LOCAL RAC BOUNDARIES · SUBSTRATE FRAME'
    plt.rcParams.update({'font.size':8,'axes.titlesize':10,'pdf.fonttype':42})
    OUTPUT.mkdir(parents=True,exist_ok=True)
    runs = {d:prepare(d) for d in (1,10,50)}
    if args.local_rac and any('early_termination' in tr.metadata for tr,_,_ in runs.values()):
        annotation += ' · TRUNCATED BEFORE GRID CROSSING'
    norms = {}
    for name,_,_,signed in FIELDS:
        lo = min(tr.fields[name].min() for tr,_,_ in runs.values())
        hi = max(tr.fields[name].max() for tr,_,_ in runs.values())
        if signed:
            bound = max(abs(lo),abs(hi),1e-12)
            norms[name] = TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound)
        else:
            norms[name] = Normalize(0,1) if name=='clutch_bound' else Normalize(lo,hi)
    for density,(tr,edges,length) in runs.items():
        fig = plt.figure(figsize=(21,11),layout='constrained')
        grid = fig.add_gridspec(3,10,height_ratios=[1,1,.5])
        for index,(name,label,cmap,_) in enumerate(FIELDS):
            row,col = divmod(index,5)
            ax = fig.add_subplot(grid[row,2*col:2*col+2])
            m = panel(ax,tr,edges,name,cmap,norms[name])
            ax.set_title(label.replace('cell frame','substrate frame' if args.local_rac else 'reference frame') if args.uncentered else label)
            if col==0:
                ax.set_ylabel('Laboratory position (µm)')
            fig.colorbar(m,ax=ax,orientation='horizontal',fraction=.05,pad=.035)
        for columns,values,label,color in (
            (slice(0,5),.5*(edges[:,0]+edges[:,-1])-10.5,'Cell-midpoint displacement (µm)','#8c2d04'),
            (slice(5,10),length,'Total length (µm)','#245e89'),
        ):
            ax=fig.add_subplot(grid[2,columns])
            ax.plot(tr.time-200,values,color=color,lw=1.5)
            ax.axvspan(0,160,color='#f2b134',alpha=.2)
            ax.set_xscale('symlog',linthresh=200)
            ax.set_xlim(-200,tr.time[-1]-200)
            ax.set_xlabel('Time from cue onset (s; symlog)')
            ax.set_ylabel(label)
        fig.suptitle(f'L₀ = 21 µm · fixed substrate clutches · density = {density} sites/µm{annotation}\n'
                     + ('Displacements solved in substrate frame; illustrative clutch parameters' if args.local_rac else 'Coupled dynamics; zero net substrate traction; illustrative, uncalibrated clutch parameters'),fontsize=14)
        save(fig,f'Figure2_L21_fixed_substrate_clutch_density_{density}{suffix}')

    fig,axes=plt.subplots(3,3,figsize=(14,13),layout='constrained')
    comparison = [FIELDS[i] for i in (0,5,9)]
    common_ylim = (min(edges.min() for _,edges,_ in runs.values())-.5,
                   max(edges.max() for _,edges,_ in runs.values())+.5)
    for row,(density,(tr,edges,length)) in enumerate(runs.items()):
        for col,(name,label,cmap,_) in enumerate(comparison):
            m=panel(axes[row,col],tr,edges,name,cmap,norms[name])
            axes[row,col].set_ylim(*common_ylim)
            axes[row,col].set_title(label.replace('cell frame','substrate frame' if args.local_rac else 'reference frame') if args.uncentered else label)
            if col==0:
                axes[row,col].set_ylabel(f'{density} sites/µm\nLaboratory position (µm)')
    for col,(name,label,cmap,_) in enumerate(comparison):
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norms[name],cmap=cmap),
                     ax=axes[:,col].tolist(),orientation='horizontal',fraction=.025,pad=.025,
                     label=label.replace('cell frame','substrate frame' if args.local_rac else 'reference frame') if args.uncentered else label)
    fig.suptitle(f'L₀ = 21 µm · moving kymographs with fixed substrate clutches{annotation}\n'
                 'Shared color and position scales; black curves mark the cell boundaries',fontsize=14)
    save(fig,'Figure2_L21_fixed_substrate_clutch_comparison'+suffix)

    # Raw fixed-lattice histories independently display anchoring to substrate.
    fig,axes=plt.subplots(3,2,figsize=(11,12),layout='constrained')
    for row,(density,(tr,edges,_)) in enumerate(runs.items()):
        with np.load(CAMPAIGN/f'density_{density}/fixed_substrate_history.npz') as z:
            t=z['time']-200
            ids=np.unique(np.r_[np.flatnonzero(t<=800),np.arange(0,len(t),4),len(t)-1])
            arrays=[z['bound_fraction'],z['extension_moment']]
            for col,(a,label) in enumerate(zip(arrays,('Bound fraction','Bound extension moment n × e (µm)'))):
                bound=max(float(np.nanmax(np.abs(a))),1e-12)
                norm=Normalize(0,1) if col==0 else TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound)
                ax=axes[row,col]
                m=ax.pcolormesh(_edges(t[ids]),_edges(z['substrate_x']),a[ids].T,
                                cmap='viridis' if col==0 else 'RdBu_r',norm=norm,rasterized=True)
                for j in (0,-1): ax.plot(tr.time-200,edges[:,j],color='black',lw=.7)
                ax.set_xscale('symlog',linthresh=200)
                ax.set_xlim(-200,t[-1]);ax.set_ylim(edges.min()-.5,edges.max()+.5)
                ax.set_title(label);ax.set_xlabel('Time from cue onset (s; symlog)')
                ax.set_ylabel(f'{density} sites/µm\nFixed substrate position (µm)')
                fig.colorbar(m,ax=ax,orientation='horizontal',fraction=.05)
    fig.suptitle(f'Adhesion states on the stationary substrate grid{annotation}\nWhite regions are outside the recorded contact footprint',fontsize=13)
    save(fig,'Figure2_L21_fixed_substrate_adhesion_history'+suffix)
    if args.passive_rear:
        fig,axes=plt.subplots(3,3,figsize=(16,11),layout='constrained')
        for row,(density,(tr,edges,length)) in enumerate(runs.items()):
            t=tr.time-200
            f=tr.fields
            for side in ('left','right'):
                axes[row,0].plot(t,f[f'{side}_boundary_displacement'][:,0],label=side)
            for field,label in (('passive_membrane_stress','Membrane'),('passive_cortical_stress','Cortex'),('right_total_boundary_stress','Total')):
                axes[row,1].plot(t,f[field][:,0],label=label)
            axes[row,1].axhline(tr.metadata['resting_total_stress'],color='gray',ls=':',label='Resting total')
            axes[row,2].plot(t,f['passive_membrane_velocity'][:,0],label='Membrane',lw=2)
            axes[row,2].plot(t,f['passive_cortical_velocity'][:,0],label='Cortex',ls='--',lw=1)
            for col,title in enumerate(('Endpoint displacement (µm)','Passive-edge stress','Passive-edge velocity (µm/s)')):
                ax=axes[row,col]
                ax.set_title(title)
                ax.set_xscale('symlog',linthresh=200)
                ax.set_xlim(-200,t[-1])
                ax.axvspan(0,160,color='#f2b134',alpha=.15)
                ax.set_xlabel('Time from cue onset (s; symlog)')
                ax.legend(frameon=False,fontsize=7)
            axes[row,0].set_ylabel(f'{density} sites/µm')
        fig.suptitle('Option 2 boundary audit: passive edge moves with cortex at fixed resting total stress',fontsize=13)
        save(fig,'Figure2_L21_passive_rear_boundary_audit')
    summary={str(d):json.loads((CAMPAIGN/f'density_{d}/summary.json').read_text()) for d in runs}
    (CAMPAIGN/'figure_summary.json').write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'runs':summary},indent=2)+'\n')
    print(OUTPUT/f'Figure2_L21_fixed_substrate_clutch_comparison{suffix}.png')


if __name__=='__main__':
    main()
