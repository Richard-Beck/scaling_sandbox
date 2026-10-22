"""Six-panel Holmes spatial-outcome candidate from retained matched-grid data."""
from pathlib import Path
from datetime import datetime,timezone
import os,sys,json,csv,hashlib,shutil,re
os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.signal import find_peaks,peak_widths
from length_dev_output_layout import dev_output_path,dev_artifact_key
ROOT=Path(__file__).resolve().parents[2]
WORKSPACE=ROOT.parents[1]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/main_candidate_2_v2_20261017'
SIZE=ROOT/'data/runs/holmes_size_bridge_20261017'
FINE=ROOT/'data/runs/holmes_distance_robustness_20261017'
SOURCES={};FILES=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(p,expected=None):
    p=Path(p);digest=sha(p)
    if expected is not None:assert digest==expected,p
    SOURCES[str(p.relative_to(ROOT))]=digest;return p

def read(p):return json.loads(source(p).read_text())
def savejson(name,value):
    p=dev_output_path(OUT,name);p.write_text(json.dumps(value,indent=2)+'\n');FILES.append(p);return p

def savecsv(name,rows):
    p=dev_output_path(OUT,name)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    FILES.append(p);return p

def load_case(C,d,folder,basal_path,label):
    run=SIZE if C==30 else FINE;status=read(folder/'status.json');assert status['complete'] and status['status']=='complete'
    expected=status.get('checkpoint_sha256',status.get('input_sha256'))
    with np.load(source(basal_path,expected)) as z:basal=z['state'].copy();x=z['x'].copy()
    n=basal.shape[1];dx=C/n;assert dx==.025;np.testing.assert_allclose(x,(np.arange(n)+.5)*dx,rtol=0,atol=1e-14)
    assert int(np.argmax(basal[0]))==0
    pools=np.array([basal[:2].sum(),basal[2:].sum()])*dx
    np.testing.assert_allclose(pools,[7.5*C,3.1*C],rtol=2e-14)
    phases={};prev=basal;error=0
    for record in status['phases']:
        with np.load(source(run/record['path'],record['sha256'])) as z:
            t=z['time'].copy();fields=z['states'].copy();np.testing.assert_array_equal(fields[0],prev)
            current=np.column_stack([fields[:,:2].sum(axis=(1,2)),fields[:,2:].sum(axis=(1,2))])*dx
            error=max(error,float(np.max(np.abs(current/pools-1))));assert error<1e-7
            phase=record['phase'];assert t[-1]==(1000 if phase=='cue' else 13000)
            phases[phase]=dict(time=t,primary=fields[:,0].copy());prev=fields[-1].copy()
    cue,rec=phases['cue'],phases['recovery'];t=np.r_[cue['time'],1000+rec['time'][1:]];u=np.r_[cue['primary'],rec['primary'][1:]];keep=t<=4000
    display=(x-x[0]+10)%C-10;order=np.argsort(display)
    archive=dev_output_path(OUT,f'C{C}_d{d}_plotted_native_series.npz')
    np.savez_compressed(archive,x=x,time_after_cue_onset=t[keep],active_rac=u[keep],basal=basal[0],withdrawal=cue['primary'][-1],final_13000_cue_free=rec['primary'][-1],circumference=C,distance=d,dx=dx)
    FILES.append(archive)
    peaks,_=find_peaks(np.tile(rec['primary'][-1],3),prominence=.05);peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
    return dict(C=C,d=d,label=label,x=x,basal=basal,phases=phases,time=t[keep],primary=u[keep],display=display,order=order,dx=dx,
        final_peaks=peaks,pool_error=error,pools=pools,basal_path=basal_path,folder=folder)

def main():
    source(Path(__file__));source(Path(__file__).with_name('length_dev_output_layout.py'))
    manifests=[(SIZE,read(SIZE/'manifest.json')),(FINE,read(FINE/'manifest.json'))]
    relevant=['code/core/models/holmes_model3.py','code/core/src/polarity1d/holmes.py','code/mechanisms/ring_distance_prototype_20261017.py','code/mechanisms/chemical_ring_distance_pilot_20261017.py']
    for run,m in manifests:
        for key in relevant:source(run/'frozen_sources'/key,m['source_sha256'][key])
    for key in relevant:assert manifests[0][1]['source_sha256'][key]==manifests[1][1]['source_sha256'][key]
    cases=[load_case(30,12,SIZE/'workers/00',SIZE/'workers/00/qualified_basal.npz','Replacement'),
           load_case(40,12,FINE/'workers/fine_0',FINE/'inputs/fine_basal.npz','Merger'),
           load_case(40,16,FINE/'workers/fine_1',FINE/'inputs/fine_basal.npz','Two retained sites')]
    np.testing.assert_array_equal(cases[1]['basal'],cases[2]['basal'])
    # Archive the previous discussion once, leaving its figure and generator unchanged.
    discussion=WORKSPACE/'current manuscript figs/main_candidate_2_20261017.md';v1=discussion.with_name('main_candidate_2_v1_20261017.md')
    if not v1.exists():shutil.copy2(discussion,v1)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,3,figsize=(8,5.5),sharex=True,sharey='row')
    fig.subplots_adjust(left=.08,right=.91,bottom=.10,top=.87,hspace=.29,wspace=.17)
    vmax=np.ceil(max(c['primary'].max() for c in cases)*10)/10;norm=Normalize(0,vmax)
    rows=[];summary=[];peak_rows=[]
    for col,c in enumerate(cases):
        C,d=c['C'],c['d'];x,o=c['display'],c['order'];ax=axes[0,col]
        ax.pcolormesh(x[o],c['time'],c['primary'][:,o],cmap='viridis',norm=norm,shading='auto',rasterized=True)
        ax.axhline(1000,ls='--',lw=.7,color='white');ax.text(-9.2,1080,'Cue off',color='white',fontsize=7,va='bottom');ax.set_title(f'{chr(65+col)}   C = {C}, d = {d} µm',loc='left',fontsize=9.5,pad=11)
        ax.set(ylim=(0,4000),yticks=[0,1000,2000,3000,4000]);ax.text(0,4050,'old',ha='center',va='bottom',fontsize=7,color='#777777');ax.text(d,4050,'cue',ha='center',va='bottom',fontsize=7,color='#9b6a00')
        ax=axes[1,col]
        snapshots=[('Basal (0 s)',c['basal'][0],'#6e6e6e','--'),('Cue off (1000 s)',c['phases']['cue']['primary'][-1],'#bf7908','-'),('14,000 s (13,000 s cue-free)',c['phases']['recovery']['primary'][-1],'#087f8c','-')]
        for label,values,color,style in snapshots:
            tiled=np.tile(values,3);allpeaks,props=find_peaks(tiled,prominence=.05);widths=peak_widths(tiled,allpeaks,rel_height=.5)[0]*c['dx']
            for k,j in enumerate(allpeaks):
                if len(values)<=j<2*len(values):
                    cell=int(j-len(values));peak_rows.append(dict(circumference=C,distance=d,phase=label,peak_cell=cell,arc_from_original_incumbent=float(x[cell]),peak_height=float(values[cell]),local_prominence=float(props['prominences'][k]),local_half_prominence_width=float(widths[k])))
            ax.plot(x[o],values[o],color=color,ls=style,lw=1.3,label=label)
            for j in o:rows.append(dict(circumference=C,distance=d,phase=label,native_x=float(c['x'][j]),arc_from_original_incumbent=float(x[j]),active_rac=float(values[j])))
        ax.set_title(f'{chr(68+col)}   '+c['label'],loc='left',fontsize=9.5,pad=9);ax.set(ylim=(0,2.6),yticks=[0,1,2]);ax.set_xlabel('Arc from original incumbent (µm)')
        for row in range(2):
            ax=axes[row,col];ax.axvline(0,ls=':',lw=.7,color='white' if row==0 else '#777777');ax.axvline(d,ls=':',lw=.7,color='#F4C66A' if row==0 else '#b28015');ax.set(xlim=(-10,30),xticks=[-10,0,10,20,30])
            if C==30:
                ax.axvspan(20,30,color='#eeeeee',zorder=0);ax.axvline(20,color='#aaaaaa',lw=.5)
        axes[0,col].set_facecolor('#eeeeee')
        if C==30:axes[0,col].text(25,2100,'Outside\nC = 30',ha='center',va='center',fontsize=7,color='#777777')
        final=c['phases']['recovery']['primary'][-1];summary.append(dict(circumference=C,distance=d,cells=len(x),dx=c['dx'],outcome=c['label'],rac_pool=float(c['pools'][0]),rho_pool=float(c['pools'][1]),final_peak_count=len(c['final_peaks']),final_peak_arcs=';'.join(f'{x[j]:.5f}' for j in c['final_peaks']),max_pool_relative_error=c['pool_error']))
    axes[0,0].set_ylabel('Time after cue onset (0–4000 s)');axes[1,0].set_ylabel('Active Rac (µM)')
    cax=fig.add_axes([.93,.537,.014,.333]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax);cb.set_label('Active Rac (µM)',fontsize=8);cb.ax.tick_params(labelsize=7)
    h,l=axes[1,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',bbox_to_anchor=(.5,.995),ncol=3,frameon=False,fontsize=8)
    for ext in ('png','pdf','svg'):
        path=dev_output_path(OUT,'main_candidate_2_v2.'+ext);fig.savefig(path,dpi=260,facecolor='white');FILES.append(path)
    plt.close(fig)
    savecsv('native_profile_snapshots.csv',rows);savecsv('case_summary.csv',summary);savecsv('native_peak_widths.csv',peak_rows)
    caption='Holmes spatial outcomes after challenger activation. Columns show C = 30 µm, d = 12 µm; C = 40 µm, d = 12 µm; and C = 40 µm, d = 16 µm, all at h = 0.025 µm. The C40 columns restart one identical fresh qualified basal array. A–C retain raw active Rac through 1000 s cue-on plus 3000 s cue-free; dashed horizontal lines mark withdrawal, dotted vertical lines the original incumbent and cue. D–F show native basal, withdrawal and endpoint after 13,000 cue-free seconds. All panels keep the original incumbent at zero and common physical axes; periodic edges connect at −10/20 µm for C30 and −10/30 µm for C40. Gray area lies outside the C30 displayed ring interval. The equal local cue has total support 3 µm and peak 1 µM/s. Replacement, merger and two retained sites follow substantial local challenger activation. Changing C also changes total material, incumbent organization, the complementary route and local basal state; it is not a pure incumbent-width intervention. Active-plus-inactive Rac/Rho pools are conserved, while active Rac alone is not. Two-site retention is a finite-window observation, not branch stability or migration persistence. At cue withdrawal, original-side local half-prominence widths are 3.94, 8.16 and 8.23 µm (basal 12.38, 16.54 and 16.54). This motivates a residual-incumbent capacity hypothesis, without identifying a width threshold or separating material from geometry. The similar C40 residual widths do not themselves explain their different fates. The predefined C50,d20 outcome remains supporting.'
    p=dev_output_path(OUT,'caption.txt');p.write_text(caption+'\n');FILES.append(p)
    # A concise standalone page accompanies the actual figure immediately.
    rel=lambda path:os.path.relpath(path,discussion.parent)
    figure=OUT/'main_candidate_2_v2.png';support=ROOT/'figures/length_mechanisms_20261008/dev'
    md=f'''# Candidate 2, version 2: activation does not determine the spatial outcome

This provisional six-panel candidate extends the [archived four-panel version 1]({v1.name}). It is a result draft for fresh-context critique, not a selected final main figure.

![Replacement, merger and retained sites]({rel(figure)})

## Proposed message

All three challenged regions activate, but later spatial outcomes differ: replacement at C = 30 µm, d = 12 µm; merger at C = 40 µm, d = 12 µm; and two sites retained through the observed cue-free window at C = 40 µm, d = 16 µm. These outcomes cannot be reduced to challenger activation or a single clock measuring loss in the original core. Raw spatial histories are needed to distinguish loss at the old position, bridging between sites and coexistence.

## Why these panels convey the message

The top row shows full raw active-Rac histories through 4000 s, preserving 1000 s with cue and the first 3000 s after withdrawal. The bottom row compares basal, withdrawal and 13,000-s cue-free endpoint profiles. It exposes the distinction between a final cap near the cue, a cap between the two sites and two narrower caps. All columns use one physical frame with the original incumbent at zero, native µM/µm/s units and one color scale. Gray space in the first column is outside its 30-µm periodic interval; no 40-µm geometry is implied there.

The first comparison changes circumference from 30 to 40 µm while keeping physical cue separation at 12 µm. The second keeps the same 40-µm circumference and **exact same N = 1600 full basal state** while moving the cue from 12 to 16 µm. The two periodic routes are respectively 12/18, 12/28 and 16/24 µm. All three runs use matched physical cell width h = 0.025 µm, unchanged native kinetics/diffusion, fresh one-site basal preparation and the same support-3-µm, peak-1-µM/s cue for 1000 s, followed by 13,000 s without cue.

The [native peak-width export]({rel(OUT/'supporting/data/native_peak_widths.csv')}) describes the residual original-side regions: local half-prominence widths change from 12.38 to 3.94 µm at C30,d12; from 16.54 to 8.16 µm at C40,d12; and from 16.54 to 8.23 µm at C40,d16 between basal and withdrawal. Basal peak concentrations are near 1.484 µM. These observations motivate a residual-incumbent capacity or architecture hypothesis. The similar C40 residual widths do not themselves explain their different fates. They do not identify a width threshold or separate material effects from geometry. Width is measured at half each contemporaneous local peak's prominence using its valleys in a periodically tiled native profile, not a frozen basal gate.

## Support, selection and scope

This is a combined size-and-separation comparison, not isolated incumbent-width causality. C30 and C40 differ in total Rac/Rho material (225/93 versus 300/124 µM·µm), basal incumbent width, the complementary route, local basal chemistry and cue fraction of the circumference. The native pool densities remain 7.5/3.1 µM. Active Rac alone is not conserved. The figure is a within-model chemical pattern result and does not demonstrate migration persistence or a biological WGD calibration.

All predefined size-path cases are disclosed in the [size report]({rel(support/'holmes_size_bridge_20261017/supporting/metadata/RESULTS.md')}): C30,d12 gives replacement; C40,d16 and C50,d20 retain two sites over the recorded window. C50 is supporting rather than omitted as a contrary outcome. The separate fixed-C40 d12 result gives merger. The [original three-model distance pilot]({rel(support/'chemical_ring_distance_20261017/supporting/metadata/RESULTS.md')}) retains every original distance and numerical control; no failed or unresolved case is relabelled as a biological outcome.

The **figure uses fine-grid h = 0.025 µm trajectories**. The [challenge mesh report]({rel(support/'holmes_distance_robustness_20261017/supporting/metadata/MESH_RESULTS.md')}) finds the same d12/d16 endpoint modes at h = 0.05/0.025 µm, with closely matching final profiles but larger transient differences around d12 merger. Exact merger time is not claimed. Separately, the [conserved-Rac imbalance continuations]({rel(support/'holmes_distance_robustness_20261017/supporting/metadata/IMBALANCE_RESULTS.md')}) use the **N = 800, h = 0.05 µm d16 endpoint**, not the fine-grid states plotted here: ±1% local active-Rac redistribution and a matched unperturbed control retain two sites for a further 13,000 s. This supports finite-window robustness under that bounded perturbation, not asymptotic or global branch stability. Slow competition, exact transition locations and causal attribution remain open.

## Raw evidence and reproducibility

- [Fine figure source summary]({rel(OUT/'supporting/data/case_summary.csv')}) and [native profile snapshots]({rel(OUT/'supporting/data/native_profile_snapshots.csv')}); the three compact NPZ files beside these retain every plotted heatmap sample in native coordinates and the three raw profiles.
- [C30 full native fields]({rel(SIZE/'workers/00')}), [C40,d12 fields]({rel(FINE/'workers/fine_0')}) and [C40,d16 fields]({rel(FINE/'workers/fine_1')}). Basal/phase identities, conserved pools, matched h and relevant frozen source hashes were checked during rendering.
- [Figure caption]({rel(OUT/'supporting/metadata/caption.txt')}), [vector PDF]({rel(OUT/'supporting/vectors/main_candidate_2_v2.pdf')}), [verification]({rel(OUT/'supporting/metadata/verification.json')}), [source/output manifest]({rel(OUT/'supporting/metadata/manifest.json')}) and [generator]({rel(Path(__file__))}). No new integration occurs in the generator.
'''
    discussion.write_text(md);p=dev_output_path(OUT,'RESULTS.md');p.write_text(md.replace(']('+str(v1.name)+')',']('+os.path.relpath(v1,p.parent)+')'))
    # Use an exact local copy of the discussion, adjusting every relative link.
    for url in re.findall(r'\]\(([^)]+)\)',md):
        target=discussion.parent/url
        md=md.replace(']('+url+')',']('+os.path.relpath(target,p.parent)+')')
    p.write_text(md);FILES.append(p)
    savejson('verification.json',dict(passed=True,matched_dx=.025,C40_basals_exactly_identical=True,exact_basal_and_withdrawal_continuity=True,relevant_native_source_hashes_match_between_runs=True,
        max_conserved_pool_relative_error=max(c['pool_error'] for c in cases),cases=summary,source_hash_count=len(SOURCES),peak_width_definition="SciPy peak_widths on threefold tiled native primary field, rel_height=0.5, local peaks prominence>=0.05 µM; descriptive widths, no front gate.",no_new_integration=True))
    savejson('manifest.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',artifact_root='../..',status='Provisional version 2 awaiting fresh critique',source_sha256=SOURCES,
        sha256={dev_artifact_key(OUT,p):sha(p) for p in FILES},reproduction='/home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__)),discussion=str(discussion),archived_v1=str(v1)))
    print(json.dumps(dict(figure=str(figure),discussion=str(discussion),cases=summary),indent=2))
if __name__=='__main__':main()
