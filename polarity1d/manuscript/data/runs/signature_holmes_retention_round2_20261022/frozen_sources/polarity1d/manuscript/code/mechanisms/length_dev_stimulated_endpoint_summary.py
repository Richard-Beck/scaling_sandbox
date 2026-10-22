"""Cue-averaged protrusion/tension and the native tension-inhibition law."""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import length_dev_output_layout as layout

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/stimulated_endpoint'
RUN=ROOT/'data/runs/length_narrative_memory_20261009'
STEM='stimulated_endpoint_length_summary'
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
GEOS=['fixed_physical','fixed_relative','fixed_slope']
LABELS=['Fixed patch: 3 µm','Relative patch: 20%','Constant slope: 0.24']
COLORS=['#2978b8','#d66a24','#45905c']
SOURCES={};ARTIFACTS={}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def register(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return p
def read(p):return json.loads(register(p).read_text())
def output(name):return layout.dev_output_path(OUT,name)
def keep(p):ARTIFACTS[str(p.relative_to(OUT))]=sha(p)
def write(name,value):
    p=output(name);p.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n');keep(p)
def csvwrite(name,rows):
    p=output(name)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    keep(p)

def main():
    register(Path(__file__));register(Path(layout.__file__))
    parent=read(output('stimulated_endpoint.manifest.json'))
    assert parent['independent_audit_passed'] and parent['producer']=='finished'
    for name,digest in parent['sha256'].items():assert sha(OUT/name)==digest,name
    audit=read(output('stimulated_endpoint_audit.json'));assert audit['passed']
    path=register(output('stimulated_endpoint_traces.npz'))
    assert sha(path)==parent['sha256'][str(path.relative_to(OUT))]
    manifest=read(RUN/'manifest.json');assert manifest['center_boundary_motion'] is False
    summary=[];checks=[]
    with np.load(path) as z:
        for geo in GEOS:
            for L in LENGTHS:
                pre=f'{geo}_L{L}_pulse_';t=z[pre+'time'];T=z[pre+'tension'];drive=z[pre+'drive'];s=z[pre+'switch']
                D=float(z[f'{geo}_L{L}_cue_duration_s'])
                assert t[0]==0 and t[-1]==D and np.all(np.diff(t)>0)
                assert np.isfinite(T).all() and np.isfinite(drive).all()
                expected=.15*s/(1+T*T/10)
                np.testing.assert_allclose(drive,expected,rtol=1e-14,atol=1e-16)
                def average(y):
                    value=float(np.trapezoid(y,t)/D)
                    independent=float(np.sum(np.diff(t)*(y[:-1]+y[1:])/2)/D)
                    assert abs(value-independent)<1e-12
                    assert y.min()-1e-12<=value<=y.max()+1e-12
                    return value
                row=dict(condition=geo,length_um=L,cue_duration_s=D,
                    mean_protrusion_drive_um_s=average(drive),mean_signed_endpoint_tension_pN=average(T),
                    minimum_endpoint_tension_pN=float(T.min()),maximum_endpoint_tension_pN=float(T.max()),
                    rms_endpoint_tension_pN=float(np.sqrt(average(T*T))),mean_rac_switch=average(s),
                    mean_saturated_drive_um_s=average(.15/(1+T*T/10)),
                    drive_at_mean_tension_um_s=float(.15/(1+average(T)**2/10)),
                    native_field_tail_pass=bool(z[f'{geo}_L{L}_field_tail_pass']))
                summary.append(row)
                checks.append(dict(condition=geo,length_um=L,passed=True,exact_cue_interval=True,
                    time_weighted_average_independently_checked=True,pointwise_drive_law_checked=True))
    # Domain of panel C: actual left AND right endpoint-adjacent tension samples
    # during the cue AND complete recovery, not an arbitrary plotting interval.
    endpoints=[]
    for j in manifest['jobs']:
        c=j['case'];geo=c['condition'];L=int(c['length'])
        path=ROOT/j['source'] if j['source'] else RUN/f"workers/{manifest['missing_jobs'].index(j):03d}.npz"
        rel=str(path.relative_to(ROOT));assert sha(path)==parent['source_sha256'][rel]
        register(path)
        with np.load(path) as z:
            for phase in ['pulse','recovery']:
                a=z[phase+'_field_tension']
                for side,col in [('left',0),('right',-1)]:
                    v=a[:,col];assert np.isfinite(v).all()
                    endpoints.append(dict(condition=geo,length_um=L,phase=phase,endpoint=side,
                        minimum_tension_pN=float(v.min()),maximum_tension_pN=float(v.max())))
        print(geo,L,'verified',flush=True)
    low=min(r['minimum_tension_pN'] for r in endpoints)
    high=max(r['maximum_tension_pN'] for r in endpoints)
    assert low<0<high
    tension=np.unique(np.r_[np.linspace(low,high,1501),-np.sqrt(10),0,np.sqrt(10)])
    curve=.15/(1+tension*tension/10)
    np.testing.assert_allclose(.15/(1+np.array([0.,-np.sqrt(10),np.sqrt(10)])**2/10),[.15,.075,.075])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,3,figsize=(13.5,4.8),layout='constrained')
    for ax,letter in zip(axs,'ABC'):
        ax.text(-.13,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
        ax.grid(alpha=.13)
    for geo,label,color in zip(GEOS,LABELS,COLORS):
        rows=[r for r in summary if r['condition']==geo]
        xs=[r['length_um'] for r in rows]
        axs[0].plot(xs,[r['mean_protrusion_drive_um_s'] for r in rows],'-o',color=color,lw=1.7,ms=4,label=label)
        axs[1].plot(xs,[r['mean_signed_endpoint_tension_pN'] for r in rows],'-o',color=color,lw=1.7,ms=4)
        axs[1].fill_between(xs,[r['minimum_endpoint_tension_pN'] for r in rows],
            [r['maximum_endpoint_tension_pN'] for r in rows],color=color,alpha=.11,lw=0)
    axs[0].set(xlabel='Nominal cell length (µm)',ylabel='Mean protrusion drive during cue (µm/s)',ylim=(0,.15),xlim=(17,40))
    axs[1].set(xlabel='Nominal cell length (µm)',ylabel='Stimulated-end tension during cue (pN)',xlim=(17,40))
    axs[1].axhline(0,color='#888',lw=.7)
    for ax in axs[:2]:ax.set_xticks([18,23,27,31,35,39])
    axs[2].plot(tension,curve,color='#343d50',lw=2)
    axs[2].axvline(0,color='#999',lw=.7)
    axs[2].plot([-np.sqrt(10),np.sqrt(10)],[.075,.075],'o',ms=4,color='#343d50')
    axs[2].set(xlabel='Endpoint tension, T (pN)',ylabel='Protrusion drive (µm/s)',xlim=(low,high),ylim=(0,.184))
    axs[2].text(.5,.97,'drive(T) = 0.15 / (1 + T²/10)\nRac switch saturated',ha='center',va='top',
        transform=axs[2].transAxes,fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.95,pad=3))
    handles=[Line2D([],[],color=c,lw=1.7,marker='o',ms=4,label=l) for c,l in zip(COLORS,LABELS)]
    handles += [Line2D([],[],color='#555',lw=1.7,label='Time-weighted cue mean'),
                Patch(facecolor='#888',alpha=.18,label='Tension range during cue')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=9,frameon=False)
    for ext in ['png','pdf','svg']:
        path=output(STEM+'.'+ext);fig.savefig(path,dpi=220,facecolor='white');keep(path)
    plt.close(fig)
    csvwrite(STEM+'_measurements.csv',summary)
    csvwrite(STEM+'_endpoint_tension_ranges.csv',endpoints)
    csvwrite(STEM+'_drive_tension_curve.csv',[dict(tension_pN=float(T),saturated_protrusion_drive_um_s=float(v)) for T,v in zip(tension,curve)])
    caption=('A: time-weighted mean stimulated-right-end protrusion drive over the complete finite cue, from onset to withdrawal, at all twelve eligible lengths and all three native cue geometries. Each trial uses its own canonical shortest tested 1%-strength-passing duration; these are not averages over a common exposure duration. The actual Rac switch is retained, including startup. B: time-weighted mean signed tension at the same stimulated-right-end native sample and over the same cue interval. Shading spans the recorded temporal minimum and maximum at each length; it is not an uncertainty interval. Tension is sampled in the last native finite-volume cell, x=L−L/(2×201), the sample used by the endpoint protrusion law. Averages use trapezoidal integration of all retained samples, including the exact fractional cutoff, divided by cue duration. The averaged drive is evaluated pointwise before integration; it is not the drive at mean tension. RMS tension and saturated-drive means are exported as additional diagnostics, since sign cancellation and the nonlinear law prevent deriving A from B alone. C: the native rightward protrusion law with saturated Rac switch, drive(T)=0.15/(1+T²/10), over the full recorded range of left and right endpoint-adjacent tensions during both cue and complete 10000-s recovery across these same 36 trials. Both tension signs inhibit equally. Dots mark half-maximal drive at T=±sqrt(10) pN. The complete law multiplies this curve by the Rac switch, 0.5×(1+tanh((Rac−1.3×resting Rac)/0.01)); every cue-on trace is saturated from 10 s onward. The drive is evaluated from saved states for the next native mechanical step, preserving the preceding figure’s up-to-0.5-s stagger convention. Only cue means and the tension-response law appear in this revised figure. Existing trajectories are preserved. No new integration.')
    path=output(STEM+'_caption.txt');path.write_text(caption+'\n');keep(path)
    write(STEM+'_audit.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
        trajectories=len(checks),checks=checks,new_integrations=0,observed_endpoint_tension_range_pN=[low,high],
        endpoint_range_scope='Both ends, cue and full recovery, all 36 native reversal trials',
        endpoint_range_records=len(endpoints),saturated_law_zero_and_half_max_checks_passed=True))
    for name,digest in ARTIFACTS.items():assert sha(OUT/name)==digest
    write(STEM+'.manifest.json',dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        independent_audit_passed=True,new_integrations=0,trajectories=36,legend=caption,
        sha256=dict(ARTIFACTS),source_sha256=SOURCES,reproduction=[
            'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))]))
    print(output(STEM+'.png'))

if __name__=='__main__':main()
