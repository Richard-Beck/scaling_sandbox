"""Render one reusable polarity-vector component from the eight retained traces."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/global_polarity_vector_20261017'
COLORS={'wave_pinning':'#3976a8','goryachev':'#c96926','holmes_model3':'#8560a4'}
NAMES={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','holmes_model3':'Holmes 3'}
SOURCES={}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(2**20),b''):h.update(b)
    return h.hexdigest()

def source(path):
    path=Path(path).resolve();SOURCES[str(path.relative_to(REPO))]=sha(path);return path

def read(name):
    with source(OUT/'supporting/data'/name).open() as stream:return list(csv.DictReader(stream))

def choose(rows,model,C,d):
    rr=[r for r in rows if r['model']==model and float(r['circumference'])==C and np.isclose(float(r['fraction']),d/C)]
    return {key:np.array([float(r[key]) for r in rr]) for key in ('time','vector_x','vector_y','relative_magnitude','total_active')}

def cartesian(ax,title):
    theta=np.linspace(0,2*np.pi,361)
    ax.plot(np.cos(theta),np.sin(theta),color='#d1d1d1',lw=.8,zorder=0)
    ax.axhline(0,color='#dddddd',lw=.6);ax.axvline(0,color='#dddddd',lw=.6)
    ax.set(xlim=(-1.15,1.15),ylim=(-1.15,1.15),xticks=[-1,0,1],yticks=[-1,0,1],
        xlabel='Activity component Re P',ylabel='Activity component Im P')
    ax.set_aspect('equal',adjustable='box');ax.set_title(title,loc='left',fontsize=9,pad=8)

def vector_path(ax,z,color,label,style='-'):
    keep=z['time']<=4000
    ax.plot(z['vector_x'][keep],z['vector_y'][keep],color=color,lw=1.35,ls=style,label=label)
    for time,marker in ((1000,'s'),(4000,'^')):
        i=np.flatnonzero(z['time']==time);assert len(i)==1;i=i[0]
        ax.plot(z['vector_x'][i],z['vector_y'][i],marker=marker,ms=5,color=color,mec='white',mew=.35,zorder=6)
    ax.plot(1,0,'o',mfc='white',mec='black',ms=5,mew=.7,zorder=8)

def draw(rows):
    fig,axes=plt.subplots(2,2,figsize=(8.8,7.2))
    fig.subplots_adjust(left=.09,right=.98,bottom=.13,top=.86,hspace=.60,wspace=.40)
    a,b,c,d=axes.flat
    cartesian(a,'A  Common geometry: C40, d16')
    for model in COLORS:
        z=choose(rows,model,40,16)
        vector_path(a,z,COLORS[model],NAMES[model])
        keep=z['time']<=4000
        b.plot(z['time'][keep],z['relative_magnitude'][keep],color=COLORS[model],lw=1.35,label=NAMES[model])
    a.legend(loc='lower left',frameon=False,fontsize=7,handlelength=1.7,labelspacing=.4)
    cancel=choose(rows,'holmes_model3',40,20);keep=cancel['time']<=4000
    b.plot(cancel['time'][keep],cancel['relative_magnitude'][keep],color='#707070',ls='--',lw=1.3,label='Holmes C40 d20: cancellation')
    b.axhline(1,color='#bdbdbd',lw=.7,ls=':')
    b.axvline(1000,color='#888888',lw=.75,ls='--')
    b.set(xlim=(0,4000),ylim=(0,1.15),xticks=[0,1000,2000,4000],yticks=[0,.5,1],
        xlabel='Time from cue onset (native units)',ylabel='Relative magnitude |P|')
    b.set_title('B  Strength can fall without losing activity',loc='left',fontsize=9,pad=8)
    b.legend(loc='upper right',bbox_to_anchor=(1,.80),frameon=False,fontsize=6.7,labelspacing=.3,handlelength=1.8)
    b.text(1460,.15,'Antipodal Holmes: two sites\ntotal active Rac 33.13→33.05',fontsize=6.6,color='#555555',va='center')
    for C,style in ((40,'-'),(80,'--')):
        z=choose(rows,'goryachev',C,C/2);keep=z['time']<=4000
        assert np.max(np.abs(z['vector_y']))<1e-12
        c.plot(z['time'][keep],z['vector_x'][keep],color=COLORS['goryachev'],lw=1.5,ls=style,
            label=f'C{C}, d{C//2}: '+('replacement' if C==40 else 'incumbent return'))
    c.axvline(1000,color='#888888',lw=.75,ls='--');c.axhline(0,color='#bdbdbd',lw=.8)
    c.set(xlim=(0,4000),ylim=(-1.15,1.15),xticks=[0,1000,2000,4000],yticks=[-1,0,1],
        xlabel='Time from cue onset (native units)',ylabel='Signed component Re P (Im P ≈ 0)')
    c.set_title('C  Goryachev: replacement versus return',loc='left',fontsize=9,pad=8)
    c.legend(loc='lower right',bbox_to_anchor=(1,.2),frameon=False,fontsize=7,handlelength=2.4,labelspacing=.5)
    cartesian(d,'D  Holmes: distinct local outcomes')
    for C,distance,color,style,label in ((40,12,'#2a9288','-','C40 d12: merged region'),
            (40,16,COLORS['holmes_model3'],'-','C40 d16: two sites'),
            (30,12,'#65747a','--','C30 d12: cue single (size changes)')):
        vector_path(d,choose(rows,'holmes_model3',C,distance),color,label,style)
    d.legend(loc='lower left',frameon=False,fontsize=6.5,handlelength=1.7,labelspacing=.4)
    fig.text(.54,.975,'Global activity polarity through site competition',ha='center',va='top',fontsize=10,weight='bold')
    fig.text(.54,.932,r'$P(t)=M_1(t)/M_1(0)$: initial vector (1,0); magnitude may exceed 1.',ha='center',fontsize=8)
    fig.text(.54,.902,'A/D: paths in activity-vector space. All panels show 0–4000 native units; full records reach 14,000.',ha='center',fontsize=7.2)
    marker_handles=[Line2D([],[],marker='o',mfc='white',mec='black',linestyle='none',ms=5,label='Initial (1,0)'),
        Line2D([],[],marker='s',color='#555555',linestyle='none',ms=5,label='Cue withdrawal: 1000'),
        Line2D([],[],marker='^',color='#555555',linestyle='none',ms=5,label='Displayed endpoint: 4000')]
    fig.legend(handles=marker_handles,loc='lower center',bbox_to_anchor=(.54,.055),ncol=3,frameon=False,fontsize=7,columnspacing=1.4)
    fig.text(.54,.023,'Width-3 cues: WP peak 0.1; Goryachev 6.4; Holmes 1. Native input/time scales are model-specific.\nC40 Holmes d12/d16 share the exact full basal state; other size comparisons change material, architecture and routes.',
        ha='center',va='center',fontsize=7,linespacing=1.4)
    artifacts=[]
    for ext in ('.png','.pdf','.svg'):
        p=layout.dev_output_path(OUT,'global_polarity_component'+ext)
        fig.savefig(p,dpi=200 if ext=='.png' else None,bbox_inches='tight');artifacts.append(p)
    plt.close(fig);return artifacts

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,
        'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    source(__file__);source(OUT/'DESIGN.md')
    audit_path=source(OUT/'supporting/metadata/global_polarity_audit.json')
    audit=json.loads(audit_path.read_text());assert audit['cases']==8 and audit['sampled_frames']==20108
    original=source(ROOT/'code/mechanisms/global_polarity_vector_20261017.py')
    assert sha(original)==audit['generator_sha256']
    for rel,digest in audit['outputs'].items():assert sha(OUT/rel)==digest
    rows=read('global_polarity_traces.csv');summary=read('global_polarity_summary.csv');sources=read('global_polarity_sources.csv')
    assert len(rows)==20108 and len(summary)==8
    for r in rows:np.testing.assert_allclose(np.hypot(float(r['vector_x']),float(r['vector_y'])),float(r['relative_magnitude']),rtol=2e-15,atol=1e-16)
    for r in sources:assert sha(source(REPO/r['path']))==r['sha256']
    artifacts=draw(rows)
    checks=dict(no_new_integration=True,retained_histories=8,retained_samples=len(rows),displayed_time_window=[0,4000],
        full_source_endpoint=14000,cartesian_limits=[-1.15,1.15],magnitude_limits=[0,1.15],
        definition='M1(t)/M1(0), complex basal normalization; no per-frame normalization',
        no_angles_or_assigned_velocities=True,antipodal_Holmes_cancellation_shown=True,
        source_code_design_and_CSVs_preserved=True,source_archive_hashes_verified=len(sources))
    path=layout.dev_output_path(OUT,'global_polarity_component_verification.json');path.write_text(json.dumps(checks,indent=2)+'\n');artifacts.append(path)
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root='../..',source_sha256=SOURCES,
        sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},checks=checks)
    path=layout.dev_output_path(OUT,'global_polarity_component.manifest.json');path.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(artifacts=len(artifacts),checks=checks),indent=2))

if __name__=='__main__':main()
