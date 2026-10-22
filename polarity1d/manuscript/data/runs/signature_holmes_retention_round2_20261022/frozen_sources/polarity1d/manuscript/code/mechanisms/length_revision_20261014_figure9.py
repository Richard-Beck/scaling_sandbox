"""Restore a mechanics row beneath the audited single-page Figure 9 budget.

All chemical, mechanical-flow and bound-budget curves reuse the hash-verified
13 October exports without integration. Drawings resample the unchanged native
relative-20% L31 history. No old constant-slope illustration data are reused.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

from length_revision_20261013_figure9 import COLORS as PRIOR_COLORS, LABELS, TERMS
from length_dev_figure_y_mass_balance import draw_dashpot, draw_spring

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'figures/length_mechanisms_20261008'
PRIOR=REPORT/'revision_20261013/debelly_figure9'
OUT=REPORT/'revision_20261014/debelly_figure9'
STEM='Figure_9_bound_MCA_balance'
LENGTH=31.
OLD_HEIGHT=13.7
NEW_HEIGHT=16.7
EXAGGERATION=2.5
SLIP_WINDOW=10.
NODES=9
SOURCES={}
COLORS={k:PRIOR_COLORS[k] for k in TERMS}
FLOW_COLORS=dict(endpoint='#6A3D9A',rho='#17BECF',membrane='#8C564B',total='#202020')
CORTEX_COLOR='#334E68'
GLOBAL_COLORS=dict(left='#BDAA22',right='#B2182B',split='#A6761D',total='#202020',uniform='#CC79A7')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(path):
    path=Path(path)
    SOURCES[str(path.relative_to(ROOT))]=sha(path)
    return path


def read(path):
    return json.loads(register(path).read_text())


def csvread(path):
    with register(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def csvwrite(name,records):
    with (OUT/name).open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(records[0]));w.writeheader();w.writerows(records)


def commands(fields,p,threshold):
    rac=fields['rac'][..., [0,-1]];tension=fields['tension'][..., [0,-1]]
    switch=.5*(1+np.tanh((rac-threshold)/p['protrusion_switch_scale']))
    result=p['rac_to_rho_coupling']*p['protrusion_speed']*switch/(1+tension**2/p['protrusion_tension_scale'])
    result[...,1]*=-1
    return result


def load():
    old=read(PRIOR/(STEM+'.manifest.json'))
    for name,digest in old['artifacts'].items():
        assert sha(register(PRIOR/name))==digest,name
    for name,digest in old['source_sha256'].items():
        assert sha(register(ROOT/name))==digest,name
    audit=read(PRIOR/(STEM+'_audit.json'))
    p=read(PRIOR/(STEM+'_parameters.json'))
    chapter=read(PRIOR/'chapter.json')
    profiles=csvread(PRIOR/(STEM+'_profiles_and_rates.csv'))
    trace=csvread(PRIOR/(STEM+'_global_pool_rates.csv'))
    threshold_path=ROOT/'figures/length_mechanisms_20261008/dev/deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv'
    threshold=float(next(r for r in csvread(threshold_path) if int(r['length_um'])==31)['rac_switch_threshold'])
    archive=register(ROOT/audit['archive'])
    assert sha(archive)==audit['archive_sha256']
    with np.load(archive) as z:
        x=z['x'].copy();t=z['pulse_time'].copy()
        ff={k:z['pulse_field_'+k].copy() for k in audit['exact_initial_checkpoint_and_phase_continuity_fields']}
    assert t[-1]==audit['cue_duration_s'] and p['center_boundary_motion'] is False
    nodes=np.linspace(0,LENGTH,NODES)
    snapshots=[];drawrows=[];summary=[]
    for key,target,label in [('basal',0.,'Basal state'),('early',30.,'30 s after cue onset'),('later',180.,'180 s after cue onset')]:
        rr=[r for r in profiles if r['snapshot']==key]
        assert len(rr)==len(x)==201
        i=int(np.argmin(abs(t-target)));assert t[i]==target
        state={k:v[i] for k,v in ff.items()}
        previous=state if i==0 else {k:v[i-1] for k,v in ff.items()}
        completed_edge=commands(previous,p,threshold)
        next_edge=commands(state,p,threshold)
        actual=lambda name:np.array([float(r[name]) for r in rr])
        np.testing.assert_array_equal(actual('bound_MCA'),state['mca_bound'])
        np.testing.assert_array_equal(actual('cortical_velocity_um_s'),state['velocity'])
        np.testing.assert_array_equal(actual('reference_position_um'),x)
        b={name:actual(name) for name in TERMS}
        comp={name:actual(name+'_velocity_um_s') for name in ['endpoint','rho','membrane']}
        endpoints=np.array([state[k][0] for k in ['left_boundary_displacement','right_boundary_displacement']])
        vm=np.zeros_like(x) if i==0 else (state['displacement']-previous['displacement'])/(t[i]-t[i-1])
        membrane_edge=np.zeros(2) if i==0 else np.array([(state[k][0]-previous[k][0])/(t[i]-t[i-1]) for k in ['left_boundary_displacement','right_boundary_displacement']])
        if i==0:
            slip=SLIP_WINDOW*state['velocity']
            endpoint_slip=SLIP_WINDOW*completed_edge
            start=0.
        else:
            start=target-SLIP_WINDOW
            first=int(np.argmin(abs(t-start)));assert t[first]==start
            integral=np.trapezoid(ff['velocity'][first:i+1],t[first:i+1],axis=0)
            slip=integral-(state['displacement']-ff['displacement'][first])
            # Each completed-step command was computed from its preceding state.
            edges=commands({k:v[first-1:i] for k,v in ff.items()},p,threshold)
            edge_integral=np.trapezoid(edges,t[first:i+1],axis=0)
            old_ends=np.array([ff[k][first,0] for k in ['left_boundary_displacement','right_boundary_displacement']])
            endpoint_slip=edge_integral-(endpoints-old_ends)
        node_u=np.interp(nodes,np.r_[0,x,LENGTH],np.r_[endpoints[0],state['displacement'],endpoints[1]])
        node_slip=np.interp(nodes,np.r_[0,x,LENGTH],np.r_[endpoint_slip[0],slip,endpoint_slip[1]])
        node_v=np.interp(nodes,np.r_[0,x,LENGTH],np.r_[completed_edge[0],state['velocity'],completed_edge[1]])
        node_b=np.interp(nodes,x,state['mca_bound'])
        node_t=np.interp(nodes,x,state['tension'])
        segment_t=np.interp(.5*(nodes[:-1]+nodes[1:]),x,state['tension'])
        xm=nodes+EXAGGERATION*node_u
        ghost=xm+EXAGGERATION*node_slip
        assert np.all(np.diff(xm)>0) and np.all(np.diff(ghost)>0)
        frame=dict(key=key,label=label,time=target,state=state,budget=b,components=comp,
                   nodes=nodes,membrane_x=xm,ghost_x=ghost,node_slip=node_slip,node_v=node_v,
                   node_b=node_b,node_t=node_t,segment_t=segment_t,
                   completed_edge=completed_edge,next_edge=next_edge,membrane_edge=membrane_edge,
                   node_membrane_velocity=np.interp(nodes,np.r_[0,x,LENGTH],np.r_[membrane_edge[0],vm,membrane_edge[1]]))
        snapshots.append(frame)
        for j,xx in enumerate(nodes):
            drawrows.append(dict(snapshot=key,time_after_onset_s=target,reference_position_um=float(xx),
                native_membrane_displacement_um=float(node_u[j]),plotted_membrane_position_um=float(xm[j]),
                plotted_schematic_cortex_position_um=float(ghost[j]),
                plotted_cortex_marker_position_um=float(ghost[j]),
                fixed_reference_10s_relative_slip_um=float(node_slip[j]),
                plotted_10s_slip_marker_position_um=float(ghost[j]),
                cortical_velocity_um_s=float(node_v[j]),membrane_velocity_um_s=float(frame['node_membrane_velocity'][j]),
                sampled_bound_MCA=float(node_b[j]),sampled_linker_drag=float(p['linker_friction']*node_b[j]),
                sampled_native_tension_pN=float(node_t[j]),geometry_exaggeration=EXAGGERATION,
                flow_arrow_um_per_um_s=40.))
        summary.append(dict(snapshot=key,time_after_onset_s=target,slip_window_start_s=start,
            slip_window_end_s=target,basal_stationary_continuation=i==0,
            completed_left_cortical_command_um_s=float(completed_edge[0]),completed_right_cortical_command_um_s=float(completed_edge[1]),
            next_left_cortical_command_um_s=float(next_edge[0]),next_right_cortical_command_um_s=float(next_edge[1]),
            left_membrane_velocity_um_s=float(membrane_edge[0]),right_membrane_velocity_um_s=float(membrane_edge[1]),
            native_left_endpoint_displacement_um=float(endpoints[0]),native_right_endpoint_displacement_um=float(endpoints[1]),
            bound_uniform_replacement_conc_s=float(actual('uniform_replacement')[0]),
            bound_correction_fraction=p['mca_on']/(p['mca_on']+p['mca_off']),drawing_nodes=NODES,native_cells=len(x)))
    trace=[{k:float(v) for k,v in r.items()} for r in trace]
    return x,p,snapshots,trace,audit,chapter,drawrows,summary


def panel(ax,letter):
    off=-.045 if ax.get_position().width>.5 else -.14
    ax.text(off,1.04,letter,transform=ax.transAxes,fontsize=13,weight='bold',va='bottom')
    ax.grid(axis='y',alpha=.15,lw=.65);ax.tick_params(labelsize=10)


def render(x,p,snapshots,trace):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                        'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    original=plt.figure(figsize=(13.1,OLD_HEIGHT))
    gs=original.add_gridspec(4,3,left=.082,right=.985,bottom=.062,top=.935,
                            height_ratios=[1.15,.95,1.4,.95],hspace=.58,wspace=.30)
    boxes={(r,c):gs[r,c].get_position(original) for r in range(4) for c in range(3)}
    global_box=gs[3,:].get_position(original)
    plt.close(original)
    fig=plt.figure(figsize=(13.1,NEW_HEIGHT))
    extra=NEW_HEIGHT-OLD_HEIGHT
    def axes(row,col):
        box=boxes[row,col]
        shift=extra if row<3 else 0.
        return fig.add_axes([box.x0,(box.y0*OLD_HEIGHT+shift)/NEW_HEIGHT,box.width,box.height*OLD_HEIGHT/NEW_HEIGHT])
    layouts=[];line_mapping=[];budget_handles=[];flow_handles=[]
    bound_max=max(float(s['state']['mca_bound'].max()) for s in snapshots)*1.08
    for j,s in enumerate(snapshots):
        ax=axes(0,j);panel(ax,'ABC'[j])
        ax.plot(x,s['state']['mca_bound'],color='#222222',lw=2.4)
        ax.set(xlim=(0,LENGTH),ylim=(0,bound_max),xticks=[0,15.5,31],title=s['label']);ax.set_xlabel('Reference position (µm)')
        if j==0:ax.set_ylabel('Bound MCA\n(native conc.)')
        ax.axvspan(.8*LENGTH,LENGTH,color='#EEEEEE',zorder=-2)
        layouts.append(layout_record(ax,'ABC'[j],fig))
        ax=axes(1,j);panel(ax,'DEF'[j]);ax.axhline(0,color='#999999',lw=.7)
        for k,label in [('endpoint','Endpoint forcing'),('rho','Rho stress gradient'),('membrane','Membrane motion')]:
            line=ax.plot(x,1000*s['components'][k],color=FLOW_COLORS[k],lw=1.9,label=label)[0]
            if j==1:flow_handles.append(line)
        line=ax.plot(x,1000*s['state']['velocity'],color=FLOW_COLORS['total'],ls='--',lw=2.2,label='Cortical flow')[0]
        if j==1:flow_handles.append(line)
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)')
        if j==0:ax.set_ylabel('Cortical velocity\n(nm/s; rightward +)')
        layouts.append(layout_record(ax,'DEF'[j],fig))
        ax=axes(2,j);panel(ax,'GHI'[j]);ax.axhline(0,color='#999999',lw=.7)
        for k in TERMS:
            expected=1e4*s['budget'][k]
            line=ax.plot(x,expected,color=COLORS[k],label=LABELS[k],ls='--' if k=='net' else '-',lw=2.3 if k=='net' else 1.65,zorder=5 if k=='net' else 3)[0]
            np.testing.assert_array_equal(line.get_ydata(),expected)
            if j==1:budget_handles.append(line)
            line_mapping.append(dict(panel='GHI'[j],snapshot=s['key'],CSV_column=k,label=line.get_label(),
                                     color=line.get_color(),linestyle=line.get_linestyle(),
                                     plotted_rate_scale=1e4,plotted_y_sha256=hashlib.sha256(expected.tobytes()).hexdigest()))
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)')
        if j==0:ax.set_ylabel('Bound-MCA rate\n(10⁻⁴ native conc./s)\nsymlog outside ±1',fontsize=10)
        ax.axvspan(.8*LENGTH,LENGTH,color='#EEEEEE',zorder=-2)
        limit=max(float(abs(s['budget'][k]).max()) for k in TERMS)*1e4*1.10
        ax.set_yscale('symlog',linthresh=1.,linscale=1.);ax.set_ylim(-limit,limit)
        ticks=[v for v in [-100,-10,-1,0,1,10,100] if abs(v)<limit]
        ax.set_yticks(ticks);ax.set_yticklabels([str(v) for v in ticks])
        layouts.append(layout_record(ax,'GHI'[j],fig))
    fig.legend(handles=flow_handles,
               loc='upper center',bbox_to_anchor=(.535,((boxes[1,0].y1+.050)*OLD_HEIGHT+extra)/NEW_HEIGHT),
               frameon=False,ncol=4,fontsize=10,columnspacing=2.)
    fig.legend(handles=budget_handles,loc='upper center',
               bbox_to_anchor=(.535,(boxes[2,0].y1*OLD_HEIGHT+extra+.85)/NEW_HEIGHT),
               frameon=False,ncol=3,fontsize=11,columnspacing=2.2)
    # The original G-I axes stay exactly the same physical size and position
    # relative to all A-I axes. Only the new row and the global panel intervene.
    scene_top=boxes[2,0].y0*OLD_HEIGHT+extra-.70
    scene_height=1.70
    scene_bottom=scene_top-scene_height
    positions=[]
    for s in snapshots:
        positions.extend([s['membrane_x'],s['ghost_x'],s['ghost_x']+40*s['node_v'],
                          s['ghost_x'][[0,-1]]+40*s['completed_edge']])
    lo=float(np.concatenate(positions).min())-.7;hi=float(np.concatenate(positions).max())+.7
    max_t=max(5.,max(float(abs(s['state']['tension']).max()) for s in snapshots))
    tn=TwoSlopeNorm(vmin=-max_t,vcenter=0,vmax=max_t);bn=Normalize(0,12)
    tc=plt.get_cmap('PuOr_r');bc=plt.get_cmap('cividis')
    scene_axes=[]
    for j,s in enumerate(snapshots):
        bb=boxes[2,j]
        ax=fig.add_axes([bb.x0,scene_bottom/NEW_HEIGHT,bb.width,scene_height/NEW_HEIGHT]);scene_axes.append(ax)
        ax.set(xlim=(lo,hi),ylim=(-.65,1.75),xticks=[],yticks=[1.5,1.,.12])
        ax.set_yticklabels(['Reference','Membrane','Cortex'] if j==0 else ['','',''])
        ax.tick_params(axis='y',length=0,labelsize=10,pad=3);ax.spines[:].set_visible(False)
        ax.text(-.14,1.02,'JKL'[j],transform=ax.transAxes,fontsize=13,weight='bold')
        ax.set_title(['Basal','30 s','180 s'][j],fontsize=11,pad=5)
        layouts.append(layout_record(ax,'JKL'[j],fig))
    fig.canvas.draw()
    for ax,s in zip(scene_axes,snapshots):
        xx=s['membrane_x'];xc=s['ghost_x'];nodes=s['nodes']
        ax.plot(nodes,np.full(NODES,1.5),color='#CCCCCC',lw=.7,zorder=0)
        ax.scatter(nodes,np.full(NODES,1.5),s=11,color='#BBBBBB',zorder=1)
        for j in [0,4,8]:ax.plot([nodes[j],xx[j]],[1.46,1.12],color='#BBBBBB',ls=':',lw=.7,zorder=0)
        for positions,y,col in [(xx,1.,'#D99C8A'),(xc,.12,CORTEX_COLOR)]:
            ax.add_patch(FancyBboxPatch((positions[0]-.1,y-.10),positions[-1]-positions[0]+.2,.2,
                         boxstyle='round,pad=0.02,rounding_size=0.08',facecolor=col,edgecolor=col,alpha=.10,lw=.7,zorder=0))
        for j in range(NODES-1):
            draw_spring(ax,xx[j],xx[j+1],1.,tc(tn(s['segment_t'][j])))
            draw_dashpot(ax,[xc[j],.12],[xc[j+1],.12],CORTEX_COLOR,width=1.1,halfwidth_pt=2.)
        for j,value in enumerate(s['node_b']):
            draw_dashpot(ax,[xc[j],.22],[xx[j],.90],bc(bn(p['linker_friction']*value)),
                         width=.7+2.2*value/.1,halfwidth_pt=2.5+3.0*value/.1)
        ax.scatter(xx,np.full(NODES,1.),s=18,facecolor='white',edgecolor='#666666',lw=.7,zorder=5)
        ax.scatter(xc,np.full(NODES,.12),s=15,color=CORTEX_COLOR,marker='s',zorder=5)
        ax.quiver(xc,np.full(NODES,-.12),40*s['node_v'],np.zeros(NODES),angles='xy',scale_units='xy',scale=1,
                  color=FLOW_COLORS['total'],width=.0047,headwidth=3.5,zorder=4)
        ax.quiver(xc[[0,-1]],np.full(2,.12),40*s['completed_edge'],np.zeros(2),angles='xy',scale_units='xy',scale=1,
                  color=FLOW_COLORS['endpoint'],width=.0065,headwidth=3.5,zorder=6)
        for j in [1,4,7]:
            ax.text(xx[j],1.26,'+',color='#CC79A7',fontsize=12,ha='center',va='center')
            ax.text(.5*(xx[j]+xc[j])+.8,.53,'+',color='#CC79A7',fontsize=12,ha='center',va='center')
    key_y=scene_bottom-.08
    fig.legend(handles=[Line2D([],[],color=CORTEX_COLOR,lw=2,label='Viscous cortex · 10-s relative motion'),
                        Line2D([],[],color=FLOW_COLORS['total'],marker='>',ls='',label='Recorded cortical flow'),
                        Line2D([],[],color=FLOW_COLORS['endpoint'],marker='>',ls='',label='Endpoint command'),
                        Line2D([],[],color='#CC79A7',marker='+',ls='',label='Uniform correction: both pools')],
               loc='upper center',bbox_to_anchor=(.535,key_y/NEW_HEIGHT),frameon=False,ncol=2,fontsize=10,columnspacing=2.,
               labelspacing=.30,handlelength=1.8,handletextpad=.5,borderaxespad=0.)
    cb_y=scene_bottom-.67
    for rect,norm,cmap,label,ticks in [([.15,cb_y/NEW_HEIGHT,.25,.055/NEW_HEIGHT],tn,tc,'Spring color: tension (pN)',[-max_t,0,max_t]),
        ([.61,cb_y/NEW_HEIGHT,.25,.055/NEW_HEIGHT],bn,bc,'MCA drag (pN s / µm²)',[0,6,12])]:
        cax=fig.add_axes(rect);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=ticks)
        cb.ax.tick_params(labelsize=10,length=2);cb.set_label(label,fontsize=10,labelpad=1)
    fig.text(.535,(cb_y-.55)/NEW_HEIGHT,
             'Membrane displacement / local 10-s slip ×2.5 · flow arrows: 4 µm = 0.1 µm/s',
             fontsize=10,ha='center',va='bottom',color='#334E68')
    ax=fig.add_axes([global_box.x0,global_box.y0*OLD_HEIGHT/NEW_HEIGHT,global_box.width,global_box.height*OLD_HEIGHT/NEW_HEIGHT]);panel(ax,'M')
    tt=np.array([r['time_start_s'] for r in trace])
    for k,label,color,ls in [('total_left_outflow','Left signed outward flux',GLOBAL_COLORS['left'],'-'),
        ('total_right_outflow','Right signed outward flux',GLOBAL_COLORS['right'],'-'),('total_boundary_loss','Sum of both boundaries',GLOBAL_COLORS['total'],'--'),
        ('total_uniform_replacement','Uniform replacement',GLOBAL_COLORS['uniform'],'-'),('total_split_residual','Sequential split residual',GLOBAL_COLORS['split'],'--')]:
        ax.plot(tt,1e3*np.array([r[k] for r in trace]),color=color,lw=2.,ls=ls,label=label)
    ax.set(xlim=(0,trace[-1]['time_end_s']),xlabel='Time after cue onset (s)',ylabel='Total-pool rate\n(10⁻³ native conc.·µm/s)')
    ax.set_xticks([0,60,120,180,240,trace[-1]['time_end_s']],['0','60','120','180','240','285']);ax.get_xticklabels()[-1].set_ha('right')
    ax.legend(frameon=False,fontsize=9.7,ncol=3,loc='lower center',bbox_to_anchor=(.5,1.02))
    layouts.append(layout_record(ax,'M',fig))
    fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_visible() and obj.get_text():
            bb=obj.get_window_extent(renderer)
            if bb.x0<-1 or bb.y0<-1 or bb.x1>fig.bbox.width+1 or bb.y1>fig.bbox.height+1:outside.append(obj.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:fig.savefig(OUT/(STEM+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    return dict(figure_width_in=13.1,figure_height_in=NEW_HEIGHT,previous_figure_height_in=OLD_HEIGHT,
                drawing_node_count=NODES,native_cells=201,illustration_xrange_um=[lo,hi],
                mechanical_scene_height_in=scene_height,axes=layouts,plot_text_within_page=True,
                signed_spring_tension_color_range_pN=[-max_t,max_t],linker_drag_color_range=[0,12],
                palette=dict(bound_process=COLORS,mechanical_source=FLOW_COLORS,global_pool=GLOBAL_COLORS,cortex_structure=CORTEX_COLOR,
                             spring_tension='PuOr_r',MCA_drag='cividis'),
                bound_process_actual_line_mapping=line_mapping,process_legend_from_actual_H_lines=True)


def layout_record(ax,letter,fig):
    bb=ax.get_position();ww,hh=fig.get_size_inches()
    return dict(panel=letter,left_in=float(bb.x0*ww),bottom_in=float(bb.y0*hh),
                width_in=float(bb.width*ww),height_in=float(bb.height*hh),
                xlabel_fontsize=ax.xaxis.label.get_fontsize(),ylabel_fontsize=ax.yaxis.label.get_fontsize(),
                title_fontsize=ax.title.get_fontsize(),tick_fontsize=[t.get_fontsize() for t in ax.get_xticklabels()])


def render_companion(x,trace):
    """Palette-matched retained withdrawal companion, without native replay."""
    profiles=csvread(PRIOR/(STEM+'_profiles_and_rates.csv'))
    faces=csvread(PRIOR/(STEM+'_face_fluxes.csv'))
    stages=[]
    for key in ['early','withdrawal','recovery']:
        rr=[r for r in profiles if r['snapshot']==key]
        ff=[r for r in faces if r['snapshot']==key]
        assert len(rr)==201 and len(ff)==202
        stages.append(dict(key=key,label=rr[0]['label'],budget={k:np.array([float(r[k]) for r in rr]) for k in TERMS},
            advective_flux=np.array([float(r['bound_advective_flux']) for r in ff]),
            diffusive_flux=np.array([float(r['bound_diffusive_flux']) for r in ff]),
            face_x=np.array([float(r['reference_face_position_um']) for r in ff])))
    fig=plt.figure(figsize=(13.1,14.5))
    gs=fig.add_gridspec(4,3,left=.085,right=.98,bottom=.07,top=.935,hspace=.66,wspace=.36)
    flux_handles=[];budget_handles=[];mapping=[]
    for j,s in enumerate(stages):
        ax=fig.add_subplot(gs[0,j]);panel(ax,'ABC'[j])
        for k,label,color in [('advective_flux','Bound advection','#E69F00'),('diffusive_flux','Bound diffusion',COLORS['bound_diffusion'])]:
            line=ax.plot(s['face_x'],1e3*s[k],color=color,lw=2. if k=='advective_flux' else 1.6,label=label)[0]
            if j==1:flux_handles.append(line)
        ax.axhline(0,color='#AAAAAA',lw=.7)
        ax.set(title=s['label'],xlabel='Reference face position (µm)',xlim=(0,LENGTH),xticks=[0,15.5,31],
               ylabel='Bound-MCA flux\n(10⁻³ native conc.·µm/s)')
        ax=fig.add_subplot(gs[1,j]);panel(ax,'DEF'[j])
        for k in TERMS:
            values=1e4*s['budget'][k]
            line=ax.plot(x,values,color=COLORS[k],label=LABELS[k],ls='--' if k=='net' else '-',lw=1.7)[0]
            np.testing.assert_array_equal(line.get_ydata(),values)
            if j==1:budget_handles.append(line)
            mapping.append(dict(panel='DEF'[j],snapshot=s['key'],CSV_column=k,label=line.get_label(),color=line.get_color(),
                                plotted_y_sha256=hashlib.sha256(values.tobytes()).hexdigest()))
        ax.set(xlim=(0,LENGTH),xticks=[0,15.5,31],xlabel='Reference position (µm)',ylabel='Bound-MCA rate\n(10⁻⁴ native conc./s)')
        ax.axhline(0,color='#AAAAAA',lw=.7)
        if s['key']=='withdrawal':
            inset=ax.inset_axes([.08,.10,.69,.50])
            for k in TERMS:inset.plot(x,1e4*s['budget'][k],color=COLORS[k],ls='--' if k=='net' else '-',lw=1.)
            selected=x<=30.;limit=1.1e4*max(float(abs(s['budget'][k][selected]).max()) for k in TERMS)
            inset.set(xlim=(0,30),ylim=(-limit,limit),xticks=[0,15,30])
            inset.set_title('Bulk view: x ≤ 30 µm',fontsize=9,pad=2);inset.axhline(0,color='#AAAAAA',lw=.5);inset.tick_params(labelsize=8,pad=1)
        ax=fig.add_subplot(gs[2,j]);panel(ax,'GHI'[j]);positions=[0,1,20,180,199,200];xx=np.arange(6)
        incoming=np.array([s['advective_flux'][k] for k in positions])/ (LENGTH/201.)*1e4
        outgoing=-np.array([s['advective_flux'][k+1] for k in positions])/(LENGTH/201.)*1e4
        ax.bar(xx-.18,incoming,.34,color=GLOBAL_COLORS['left'],label='Left-face signed influx')
        ax.bar(xx+.18,outgoing,.34,color=GLOBAL_COLORS['right'],label='Right-face signed influx')
        ax.plot(xx,incoming+outgoing,'ko',ms=4,label='Net advection');ax.axhline(0,color='#AAAAAA',lw=.7)
        ax.set(xticks=xx,xticklabels=[f'{x[k]:.2f}' for k in positions],xlabel='Cell-center position (µm)',ylabel='Cell advective rate\n(10⁻⁴ native conc./s)')
        ax.tick_params(axis='x',labelrotation=30,labelsize=9)
    fig.legend(handles=flux_handles,loc='upper center',bbox_to_anchor=(.535,(gs[0,0].get_position(fig).y1*14.5+.55)/14.5),
               frameon=False,ncol=2,fontsize=11)
    fig.legend(handles=budget_handles,loc='upper center',bbox_to_anchor=(.535,(gs[1,0].get_position(fig).y1*14.5+.84)/14.5),
               frameon=False,ncol=3,fontsize=11,columnspacing=2.)
    fig.legend(handles=[Line2D([],[],color=GLOBAL_COLORS['left'],lw=6,label='Left-face signed influx'),
                        Line2D([],[],color=GLOBAL_COLORS['right'],lw=6,label='Right-face signed influx'),
                        Line2D([],[],color='black',marker='o',ls='',label='Net advection')],
               loc='upper center',bbox_to_anchor=(.535,(gs[2,0].get_position(fig).y1*14.5+.55)/14.5),
               frameon=False,ncol=3,fontsize=10,columnspacing=2.)
    ax=fig.add_subplot(gs[3,:]);panel(ax,'J');tt=np.array([r['time_start_s'] for r in trace])
    for k,label,color,ls in [('total_left_outflow','Left signed outward flux',GLOBAL_COLORS['left'],'-'),
                           ('total_right_outflow','Right signed outward flux',GLOBAL_COLORS['right'],'-'),
                           ('total_uniform_replacement','Uniform replacement',GLOBAL_COLORS['uniform'],'-'),
                           ('total_split_residual','Sequential split residual',GLOBAL_COLORS['split'],'--')]:
        ax.plot(tt,1e3*np.array([r[k] for r in trace]),color=color,lw=2.,ls=ls,label=label)
    ax.set(xlim=(0,trace[-1]['time_end_s']),xlabel='Time after cue onset (s)',ylabel='Total-pool rate\n(10⁻³ native conc.·µm/s)')
    ax.set_xticks([0,60,120,180,240,trace[-1]['time_end_s']],['0','60','120','180','240','285']);ax.get_xticklabels()[-1].set_ha('right')
    ax.legend(frameon=False,fontsize=10,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.02))
    fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_visible() and obj.get_text():
            bb=obj.get_window_extent(renderer)
            if bb.x0<-1 or bb.y0<-1 or bb.x1>fig.bbox.width+1 or bb.y1>fig.bbox.height+1:outside.append(obj.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:fig.savefig(OUT/(STEM+'_companion.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    caption=(PRIOR/(STEM+'_companion_caption.txt')).read_text().strip()
    caption+=' The process legend now uses the same canonical colors as the main bound-MCA budget. Full bound-advection flux is gold, matching Figure 11; bound diffusion is gray. Left/right face fluxes and global boundary fluxes use olive/crimson, numerical split residual uses ochre, and uniform correction uses pink. Each key is beside its panel group. All curves and bar values retain their original audited data.'
    (OUT/(STEM+'_companion_caption.txt')).write_text(caption+'\n')
    return dict(passed=True,new_integrations=0,retained_snapshots=['early','withdrawal','recovery'],native_cells=201,
                actual_process_line_mapping=mapping,plot_text_within_page=True,
                full_bound_advection_color='#E69F00',bound_diffusion_color=COLORS['bound_diffusion'],
                signed_face_colors=dict(left=GLOBAL_COLORS['left'],right=GLOBAL_COLORS['right']))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for name in ['length_revision_20261014_figure9.py','length_revision_20261013_figure9.py','length_dev_figure_y_mass_balance.py',
                 'length_dev_mca_flow_budget.py','length_dev_basal_front_tension_decomposition.py']:
        register(Path(__file__).with_name(name))
    x,p,snapshots,trace,native_audit,chapter,drawrows,summary=load()
    for suffix in ['_profiles_and_rates.csv','_face_fluxes.csv','_regional_rates.csv','_global_pool_rates.csv','_snapshot_diagnostics.csv','_parameters.json']:
        shutil.copy2(PRIOR/(STEM+suffix),OUT/(STEM+suffix))
        assert sha(OUT/(STEM+suffix))==sha(PRIOR/(STEM+suffix))
    csvwrite(STEM+'_drawing_nodes.csv',drawrows);csvwrite(STEM+'_drawing_snapshots.csv',summary)
    layout=render(x,p,snapshots,trace)
    companion_audit=render_companion(x,trace)
    caption=chapter['caption'].replace('J) Total-pool','M) Total-pool')
    caption=caption.replace('M) Total-pool',
        'J–L) Spring/dashpot illustrations at the same three times. Membrane nodes use native displacement magnified ×2.5; '
        'spring color shows signed tension. Navy axial dashpots depict cortical viscosity, and transverse dashpots '
        'encode local bound-MCA drag. The cortex marker row is displaced relative to the membrane by the 10-s fixed-reference '
        'relative-motion integral, also ×2.5; inclined viscous couplers connect the offset rows. These local-window positions '
        'do not represent stored elastic MCA strain or permanently tracked particles. '
        'Black arrows show cortical flow (4 µm per 0.1 µm/s), purple arrows the prescribed endpoint commands, and pink '
        'plus signs the uniform correction to both MCA pools. Nine drawing nodes resample the unchanged 201-cell fields. M) Total-pool')
    chapter['caption']=caption;chapter['image']=str((OUT/(STEM+'.png')).relative_to(REPORT))
    chapter['manifest']=str((OUT/(STEM+'.manifest.json')).relative_to(REPORT))
    chapter['paragraphs']=[p.replace('Figure 9J','Figure 9M') for p in chapter['paragraphs']]
    chapter['paragraphs'][0]+=' The spring/dashpot views retain membrane deformation, local cortical motion and MCA-dependent drag alongside these quantitative balances (Figure 9J–L).'
    chapter['sources']=[[str((OUT/(STEM+'_methods.txt')).relative_to(REPORT)),'Exact bound-MCA and drawing methods'],
        [str((OUT/(STEM+'_companion.pdf')).relative_to(REPORT)),'Flux, endpoint cells and recovery companion'],
        [str((OUT/(STEM+'_audit.json')).relative_to(REPORT)),'Native-data and layout audit']]
    methods=(PRIOR/(STEM+'_methods.txt')).read_text().replace('in J.','in M.').replace('in J\n','in M\n')
    methods=methods.replace('Verification. All 571 retained cue intervals',
        'Retained 13 October verification. All 571 retained cue intervals')
    methods=methods.replace('No new long integration, altered parameter, source input, or original artifact is used.',
        'No new long integration is executed; original parameters, source inputs and artifacts are preserved.')
    methods+='''\nRestored mechanics row, 2026-10-14. The original audited bound-MCA and global-rate CSVs are copied byte-for-byte. No integration is executed. Main A-I keep their original 13 October axes widths/heights, line inputs, limits, ticks and font sizes. The source canvas grows from 13.1 x 13.7 to 13.1 x 15.7 inches, with one row of scenes J-L between G-I and the global panel, now M. Report scaling remains by width only.\n\nDrawing geometry. Nine nominal nodes from 0 to 31 µm resample the original 201-cell fields. The gray row marks those undeformed reference nodes. Membrane nodes are drawn at x+2.5*u(x), including the actual independent endpoint displacements; no endpoint centering is applied. Springs connect those nodes, with a common signed native-tension color scale, sampled at nominal segment midpoints. The purple axial-dashpot row is a topology schematic registered to the membrane nodes for legibility, not an accumulated cortical-displacement field. Green MCA dashpots join the rows vertically. Their color and barrel/line widths encode linker_friction*bound_MCA sampled at the nominal nodes. The width rule includes a visible baseline and is not literally proportional to drag. No coupler extension or elastic MCA force is inferred. Native field sampling occurs at reference positions before any drawing-coordinate shift.\n\nKinematics. Purple arrows encode the completed stored cortical flow on one common scale, 4 drawn µm per 0.1 µm/s. Orange endpoint arrows encode the commands used in that completed mechanical solve, from its preceding Rac/tension state. Membrane endpoint velocity is recorded separately. Hollow diamonds below the schematic cortex mark x+2.5*u+2.5*s, where s is the time integral of recorded cortical velocity at the fixed nominal position over the previous 10 s, minus the corresponding membrane displacement change. Time quadrature is trapezoidal on the original native half-second outputs; endpoint commands use completed-step times. The basal window is its stationary continuation, with zero membrane velocity and s=10*stored cortical velocity. These local fixed-reference markers are not advected material particles or MCA bond extension. All three scenes share one coordinate extent that includes deformed endpoints, velocity arrows and slip markers. The global gray reference row supplies the displacement gauge.\n\nCorrection glyphs. Pink plus signs mark spatially uniform concentration correction in both pools, rather than a directed mechanical force or local binding reaction. Their locations and sizes are schematic; correction magnitudes remain in G-I and M. The native resting partition assigns 78.125% of the total correction to bound MCA and 21.875% to unbound MCA. The main mechanical drawings use only the same relative-20% L31 trajectory as all other main panels; no constant-slope drawing input is reused. Existing withdrawal/recovery evidence and its native endpoint-timing explanation remain linked unchanged.\n'''
    methods=methods.replace('13.1 x 15.7','13.1 x 16.7').replace('The purple axial-dashpot row','The navy axial-dashpot row')
    methods=methods.replace('Green MCA dashpots join the rows vertically.','MCA dashpots, colored on the distinct cividis drag scale, join the rows vertically.')
    methods=methods.replace('Purple arrows encode the completed stored cortical flow','Black arrows encode the completed stored cortical flow')
    methods=methods.replace('Orange endpoint arrows encode','Purple endpoint arrows encode')
    methods=methods.replace('The navy axial-dashpot row is a topology schematic registered to the membrane nodes for legibility, not an accumulated cortical-displacement field.',
        'The navy axial-dashpot row occupies data-derived local 10-s cortex illustration positions x+2.5*u+2.5*s. Its shift relative to the membrane is the retained local-window motion, not a stored cumulative cortical-displacement field.')
    methods=methods.replace('MCA dashpots, colored on the distinct cividis drag scale, join the rows vertically.',
        'Inclined MCA viscous dashpots, colored on the distinct cividis drag scale, connect the membrane positions x+2.5*u to the displaced cortex positions x+2.5*u+2.5*s.')
    methods=methods.replace('Hollow diamonds below the schematic cortex mark x+2.5*u+2.5*s,',
        'Cortex squares on the displaced axial-dashpot row occupy x+2.5*u+2.5*s,')
    methods=methods.replace('The global gray reference row supplies the displacement gauge.',
        'The undeformed gray reference row supplies the displacement gauge; membrane spring spans therefore change with the native displacement field, and the cortex row itself shows the retained relative motion.')
    methods=methods.replace('No coupler extension or elastic MCA force is inferred.',
        'Connector endpoints follow these two data-derived rows; connector tilt shows the local 10-s relative offset. The coupling is viscous, without a stored elastic MCA strain field.')
    methods=methods.replace('Existing withdrawal/recovery evidence and its native endpoint-timing explanation remain linked unchanged.',
        'The withdrawal/recovery companion retains its audited values and native endpoint-timing explanation, with palettes and local keys matched to the main figure.')
    methods+='\nLegend mapping. The budget key is directly above G-I and uses the actual six Line2D objects from H, with each CSV column, label, color and plotted-value hash exported in the audit. Mechanical forces use purple endpoint, cyan Rho stress and brown membrane motion; those hues do not duplicate the blue/orange/gray process-rate curves. The global key uses olive left flux, crimson right flux, ochre numerical split residual, pink uniform correction and black summed flux. The spring-tension palette is PuOr_r, matching the same field in Figure 8; drag uses cividis. Keys are adjacent to the panels they explain.\n'
    (OUT/(STEM+'_caption.txt')).write_text(caption+'\n');(OUT/(STEM+'_methods.txt')).write_text(methods)
    audit=dict(passed=True,new_integrations=0,unchanged_native_audit=native_audit,layout=layout,
        matched_snapshot_bound_MCA_and_flow_exactly=True,prior_budget_CSVs_byte_identical=True,
        prior_artifact_and_source_hashes_verified=True,drawing_snapshots=summary,
        geometry_definition='membrane=x+2.5*u; cortex=membrane+2.5*local_10s_relative_slip',
        first_variant_preserved='archive_20261014_registered_geometry',companion=companion_audit)
    (OUT/(STEM+'_audit.json')).write_text(json.dumps(audit,indent=2)+'\n')
    (OUT/'chapter.json').write_text(json.dumps(chapter,indent=2,ensure_ascii=False)+'\n')
    stamp=datetime.now(timezone.utc).isoformat()
    body='<!doctype html><meta charset="utf-8"><title>Restored Figure 9 mechanics</title><style>body{font:18px/1.5 sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem}img{width:100%}a{color:#17609a}</style>'
    body+='<h1>Bound-MCA balance with membrane and cortex mechanics</h1><p>Generated '+stamp+'; producer finished.</p>'
    body+=''.join('<p>'+html.escape(p)+'</p>' for p in chapter['paragraphs'])
    body+='<img src="'+STEM+'.png" alt="Bound-MCA profiles, native process budgets and restored spring/dashpot mechanics"><p>'+html.escape(caption).replace('*Figure 9*','<em>Figure 9</em>')+'</p>'
    body+='<p><a href="'+STEM+'_methods.txt">Methods</a> · <a href="'+STEM+'.pdf">Single-page vector PDF</a> · <a href="'+STEM+'_audit.json">Audit</a></p>'
    (OUT/'index.html').write_text(body+'\n')
    artifacts={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!=STEM+'.manifest.json'}
    manifest=dict(generated_at=stamp,producer='finished',passed=True,renderer=str(Path(__file__).relative_to(ROOT)),
        renderer_sha256=sha(__file__),artifacts=artifacts,sha256=artifacts,sources=SOURCES,source_sha256=SOURCES,
        caption=caption,audit=audit,new_integrations=0,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_revision_20261014_figure9.py')
    (OUT/(STEM+'.manifest.json')).write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    print(OUT/(STEM+'.png'));print(json.dumps(layout,indent=2))


if __name__=='__main__':main()
