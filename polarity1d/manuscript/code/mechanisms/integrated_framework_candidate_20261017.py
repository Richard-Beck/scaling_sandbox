"""Integrated enlargement/competition framework using retained audited evidence only."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,os,sys,re,shutil
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from length_dev_output_layout import dev_output_path,dev_artifact_key
import initial_polarity_report as basal_report
import main_candidate_1_20261017 as chemical
ROOT=Path(__file__).resolve().parents[2]
WORKSPACE=ROOT.parents[1]
BASE=ROOT/'figures/length_mechanisms_20261008'
OUT=BASE/'dev/integrated_framework_candidate_20261017'
INITIAL=ROOT/'data/runs/initial_polarity_length_corrected_20261008'
CHEM=BASE/'dev/main_candidate_1_20261017'
RING=BASE/'dev/main_candidate_2_v2_20261017'
DISCUSSION=WORKSPACE/'current manuscript figs/integrated_framework_candidate_20261017.md'
MODELS={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring'}
PAIRS={'wave_pinning':(19,39),'goryachev':(19,39),'debelly':(19,39),'holmes_model3':(20,40),'spring':(1,2)}
COLORS={'goryachev':'#893f82','wave_pinning':'#007d73'}
SOURCES={};FILES=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(p,expected=None):
    p=Path(p);value=sha(p)
    if expected is not None:assert value==expected,p
    SOURCES[os.path.relpath(p,ROOT)]=value;return p

def read(p):return json.loads(source(p).read_text())
def csvread(p,expected=None):return list(csv.DictReader(source(p,expected).open()))
def load(p,expected=None):
    with np.load(source(p,expected)) as z:return {key:z[key].copy() for key in z.files}
def csvwrite(name,rows):
    p=dev_output_path(OUT,name)
    with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    FILES.append(p);return p

def jsonwrite(name,value):
    p=dev_output_path(OUT,name);p.write_text(json.dumps(value,indent=2)+'\n');FILES.append(p);return p

def polish_prose(text,mapping):
    """Format visible prose while preserving every Markdown URL identifier."""
    pieces=re.split(r'(\]\([^)]+\))',text)
    for index in range(0,len(pieces),2):
        for old,new in mapping.items():pieces[index]=pieces[index].replace(old,new)
    return ''.join(pieces)


def collect_basal():
    prov=read(INITIAL/'report_provenance.json');metrics=csvread(INITIAL/'steady_state_metrics.csv',prov['summary_csv_sha256'])
    source(Path(basal_report.__file__))
    snapshots=[];folds=[];checks=[]
    for model,(lo,hi) in PAIRS.items():
        selected=[next(r for r in metrics if r['model']==model and float(r['domain_length'])==L) for L in (lo,hi)]
        for row in selected:
            z=load(ROOT/row['checkpoint'],row['checkpoint_sha256']);fields={k[6:]:z[k][-1].copy() for k in z if k.startswith('field_')}
            q=np.linspace(0,1,len(z['x'])) if model=='spring' else z['x']/float(row['domain_length'])
            case=dict(model=model,length=float(row['domain_length']),initial_orientation=float(row['active_polarity']),checkpoint=row['checkpoint'],checkpoint_sha256=row['checkpoint_sha256'])
            measured=basal_report.measure(case,fields,q)
            for key in ('physical_length','old_pole_excess_active','old_pole_width','old_pole_peak_active','background_active'):
                np.testing.assert_allclose(measured[key],float(row[key]),rtol=2e-12,atol=2e-13)
            name='rac' if model=='debelly' else 'rac_active' if model=='holmes_model3' else 'active'
            for j,value in enumerate(fields[name]):snapshots.append(dict(model=model,nominal_size=float(row['domain_length']),physical_length=float(row['physical_length']),position_fraction=float(q[j]),physical_position=float(q[j]*float(row['physical_length'])),active_concentration=float(value)))
            checks.append(dict(model=model,nominal_size=float(row['domain_length']),checkpoint_sha256=row['checkpoint_sha256'],native_metrics_reconstructed=True))
        a,b=selected
        folds.append(dict(model=model,model_label=MODELS[model],nominal_size_small=lo,nominal_size_large=hi,physical_size_small=float(a['physical_length']),physical_size_large=float(b['physical_length']),
            physical_size_fold=float(b['physical_length'])/float(a['physical_length']),
            excess_active_small=float(a['old_pole_excess_active']),excess_active_large=float(b['old_pole_excess_active']),excess_active_fold=float(b['old_pole_excess_active'])/float(a['old_pole_excess_active']),
            width_small=float(a['old_pole_width']),width_large=float(b['old_pole_width']),width_fold=float(b['old_pole_width'])/float(a['old_pole_width']),
            peak_concentration_small=float(a['old_pole_peak_active']),peak_concentration_large=float(b['old_pole_peak_active']),peak_concentration_fold=float(b['old_pole_peak_active'])/float(a['old_pole_peak_active']),
            checkpoint_small=a['checkpoint'],checkpoint_large=b['checkpoint']))
    return metrics,folds,snapshots,checks

def collect_functions(metrics):
    manifest=read(CHEM/'supporting/metadata/manifest.json');path=CHEM/'supporting/data/length_measurements.csv'
    rows=csvread(path,manifest['artifacts'][dev_artifact_key(CHEM,path)])
    for module in (chemical,chemical.retained):source(Path(module.__file__),manifest['source_sha256'][str(Path(module.__file__).relative_to(ROOT))])
    for model in ('goryachev','wave_pinning'):source(ROOT/f'code/core/models/{model}.py',manifest['source_sha256'][f'code/core/models/{model}.py'])
    checks=[];max_net_error=0.;output=[]
    for row in rows:
        model=row['model'];L=float(row['length_native']);metric=next(r for r in metrics if r['model']==model and float(r['domain_length'])==L)
        cp=load(ROOT/row['checkpoint'],row['checkpoint_sha256']);assert row['checkpoint_sha256']==metric['checkpoint_sha256']
        p=ROOT/row['perturbation_archive'];sham_path=p.with_name('epsilon_0.npz')
        z=load(p,manifest['source_sha256'][str(p.relative_to(ROOT))]);s=load(sham_path,manifest['source_sha256'][str(sham_path.relative_to(ROOT))])
        np.testing.assert_array_equal(z['originalstate'],cp['state']);np.testing.assert_array_equal(s['originalstate'],cp['state'])
        np.testing.assert_allclose(z['removed_inactive'],.005*cp['state'][1],rtol=0,atol=1e-13)
        np.testing.assert_array_equal(z['returnedstate'],np.stack([z['pulse_active'][-1],z['pulse_inactive'][-1]+z['removed_inactive']]))
        assert z['pulse_time'][-1]==.1 and z['recovery_time'][0]==0
        np.testing.assert_allclose(z['recovery_active'][0],z['returnedstate'][0],rtol=0,atol=1e-10)
        width=float(row['fixed_core_width_native']);dx=L/len(z['x']);w=np.clip(width-np.arange(len(z['x']))*dx,0,dx)
        np.testing.assert_allclose(width,float(metric['old_pole_width']),rtol=2e-13)
        up,vp=z['returnedstate'];us,vs=s['returnedstate'];damage=float((us-up)@w);assert damage>0
        rp,dp=chemical.local_rates(model,up,vp,L);rs,ds=chemical.local_rates(model,us,vs,L)
        net=float((rp-rs+dp-ds)@w);gain=net/damage
        np.testing.assert_allclose(gain,float(row['net_gain_per_damage_native']),rtol=3e-7,atol=3e-8)
        max_net_error=max(max_net_error,abs(gain-float(row['net_gain_per_damage_native'])))
        activation=float(cp['state'][1,-1]);np.testing.assert_allclose(activation,float(metric['opposite_inactive']),rtol=2e-14)
        native_model=chemical.retained.MODELS[model][1]
        u,v=cp['state'];delta=chemical.retained.reaction(native_model,u,v,.1)-chemical.retained.reaction(native_model,u,v,0)
        np.testing.assert_allclose(delta,.1*v,rtol=2e-12,atol=2e-14)
        substrate=float((chemical.activation_factor(model,up)*(vp-vs))@w)/damage
        core_surplus=float((vp-vs)@w);exposure=core_surplus/damage;conversion=substrate/exposure
        c=.01733 if model=='goryachev' else 1
        active=float(((chemical.activation_factor(model,up)-chemical.activation_factor(model,us))*vs-c*(up-us))@w)/damage
        diffusion=float((dp-ds)@w)/damage
        np.testing.assert_allclose(active+substrate+diffusion,gain,rtol=3e-7,atol=3e-8)
        if model=='goryachev':assert active<0 and substrate>0
        output.append(dict(model=model,length_native=L,initial_local_cue_gain_native=activation,
            immediate_net_core_repair_per_damage_native=gain,active_core_damage_native=damage,inactive_substrate_contribution_per_damage_native=substrate,
            active_field_contribution_per_damage_native=active,native_diffusion_contribution_per_damage_native=diffusion,
            weighted_conversion_coefficient_native=conversion,core_inactive_surplus_native=core_surplus,core_surplus_exposure_per_damage=exposure,
            inactive_removal_fraction=.005,withholding_duration_native=.1,checkpoint=row['checkpoint'],checkpoint_sha256=row['checkpoint_sha256'],perturbation_archive=row['perturbation_archive']))
        checks.append(dict(model=model,length_native=L,exact_checkpoint_and_return_verified=True,native_cue_identity_verified=True,net_core_repair_reconstructed=True))
    assert len(output)==36
    for model in COLORS:
        group=sorted([r for r in output if r['model']==model],key=lambda r:r['length_native']);assert len(group)==18
        reference=next(r for r in group if r['length_native']==9)
        for row in group:
            row['cue_gain_L9_reference_native']=reference['initial_local_cue_gain_native'];row['net_repair_L9_reference_native']=reference['immediate_net_core_repair_per_damage_native']
            row['cue_gain_relative_to_L9']=row['initial_local_cue_gain_native']/row['cue_gain_L9_reference_native']
            row['net_repair_relative_to_L9']=row['immediate_net_core_repair_per_damage_native']/row['net_repair_L9_reference_native']
        matched_reference=next(r for r in group if r['length_native']==19)
        for row in group:
            row['cue_gain_L19_reference_native']=matched_reference['initial_local_cue_gain_native'];row['net_repair_L19_reference_native']=matched_reference['immediate_net_core_repair_per_damage_native']
            row['cue_gain_relative_to_L19']=row['initial_local_cue_gain_native']/row['cue_gain_L19_reference_native']
            row['net_repair_relative_to_L19']=row['immediate_net_core_repair_per_damage_native']/row['net_repair_L19_reference_native']
    return output,checks,max_net_error

def matched_factorization(functions,folds):
    rows=[];summaries=[]
    for model in COLORS:
        pair=[next(r for r in functions if r['model']==model and r['length_native']==L) for L in (19,39)]
        architecture=next(r for r in folds if r['model']==model)
        for r in pair:
            rows.append(dict(model=model,length_native=r['length_native'],weighted_conversion_coefficient_native=r['weighted_conversion_coefficient_native'],core_surplus_exposure_per_damage=r['core_surplus_exposure_per_damage'],
                inactive_substrate_contribution_per_damage_native=r['inactive_substrate_contribution_per_damage_native'],net_repair_per_damage_native=r['immediate_net_core_repair_per_damage_native'],
                initial_local_cue_gain_native=r['initial_local_cue_gain_native'],core_inactive_surplus_native=r['core_inactive_surplus_native'],core_active_damage_native=r['active_core_damage_native'],perturbation_archive=r['perturbation_archive']))
        a,b=pair
        coefficient=b['weighted_conversion_coefficient_native']/a['weighted_conversion_coefficient_native'];exposure=b['core_surplus_exposure_per_damage']/a['core_surplus_exposure_per_damage'];substrate=b['inactive_substrate_contribution_per_damage_native']/a['inactive_substrate_contribution_per_damage_native']
        np.testing.assert_allclose(coefficient*exposure,substrate,rtol=3e-15)
        summaries.append(dict(model=model,length_small=19,length_large=39,peak_concentration_fold=architecture['peak_concentration_fold'],weighted_conversion_coefficient_fold=coefficient,core_surplus_exposure_per_damage_fold=exposure,
            inactive_substrate_term_fold=substrate,net_repair_per_damage_fold=b['immediate_net_core_repair_per_damage_native']/a['immediate_net_core_repair_per_damage_native'],initial_local_cue_gain_fold=b['initial_local_cue_gain_native']/a['initial_local_cue_gain_native']))
    return rows,summaries


def collect_rings():
    m=read(RING/'supporting/metadata/manifest.json');cases=[]
    for C,d,outcome in ((30,12,'Replacement'),(40,12,'Merger'),(40,16,'Retained sites')):
        name=f'C{C}_d{d}_plotted_native_series.npz';path=RING/'supporting/data'/name;z=load(path,m['sha256'][dev_artifact_key(RING,path)])
        assert float(z['dx'])==.025 and z['time_after_cue_onset'][-1]==4000
        n=len(z['x']);np.testing.assert_allclose(z['x'],(np.arange(n)+.5)*.025,rtol=0,atol=1e-14)
        np.testing.assert_array_equal(z['active_rac'][0],z['basal'])
        display=(z['x']-z['x'][0]+10)%C-10
        cases.append(dict(C=C,d=d,outcome=outcome,**z,display=display,order=np.argsort(display)))
        target=dev_output_path(OUT,name);shutil.copy2(path,target);FILES.append(target)
    np.testing.assert_array_equal(cases[1]['basal'],cases[2]['basal'])
    # Verify the native archival sources retained behind the compact exact exports.
    for rel,digest in m['source_sha256'].items():
        if rel.endswith(('/cue.npz','/recovery.npz','qualified_basal.npz','fine_basal.npz')):source(ROOT/rel,digest)
    return cases

def draw(folds,functions,rings):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.titleweight':'normal','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(8.6,6.7));fig.subplots_adjust(left=.13,right=.91,bottom=.08,top=.85,hspace=.67,wspace=.47)
    fig.text(.52,.98,'Incumbent enlargement: functional consequences and spatial modifiers',ha='center',va='top',fontsize=10,weight='bold')
    ax=axes[0,0];metrics=[('excess_active_fold','o','Excess active amount',-.19),('width_fold','s','Width',0),('peak_concentration_fold','^','Peak concentration',.19)]
    for key,marker,label,offset in metrics:
        ax.scatter([r[key] for r in folds],np.arange(5)+offset,marker=marker,s=22,facecolors='#333333' if marker=='o' else 'white',edgecolors='#333333',linewidths=.8,zorder=4,label=label)
    names=[MODELS[r['model']]+'\n'+('L₀ 1→2' if r['model']=='spring' else f"{r['nominal_size_small']}→{r['nominal_size_large']}") for r in folds]
    ax.set(yticks=np.arange(5),yticklabels=names,xlim=(.8,3.2),xticks=[1,2,3],ylim=(4.5,-.6));ax.tick_params(axis='y',labelsize=7.2,length=0)
    ax.axvline(1,color='#bbbbbb',lw=.7,ls='--');ax.set_xlabel('Larger / smaller size value')
    ax.set_title('A   Incumbent organization',loc='left',fontsize=9,pad=27)
    handles,labels=ax.get_legend_handles_labels();labels[0]='A markers: excess active amount'
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.95),ncol=3,frameon=False,fontsize=7.2,handletextpad=.5,columnspacing=1.5)
    for index,(key,title,ylabel,lim) in enumerate((('cue_gain_relative_to_L19','B   Remote cue activation','Initial cue gain / L = 19 value',(.5,2.08)),('net_repair_relative_to_L19','C   Conditional core recovery','Immediate net repair / L = 19 value',(0,4.2))),start=1):
        ax=axes[0,index]
        for model,color in COLORS.items():
            group=sorted([r for r in functions if r['model']==model],key=lambda r:r['length_native']);x=[r['length_native'] for r in group];y=[r[key] for r in group]
            ax.plot(x,y,color=color,lw=1.4,marker='o',ms=2.2,label=MODELS[model])
            chosen=[r for r in group if r['length_native'] in (19,39)];ax.scatter([r['length_native'] for r in chosen],[r[key] for r in chosen],s=18,facecolors='white',edgecolors=color,lw=.9,zorder=5)
            digits=3 if index==1 else 2
            ax.annotate(f"{group[-1][key]:.{digits}f}×",(39,group[-1][key]),xytext=(-2,-12 if index==1 and model=='goryachev' else 5),textcoords='offset points',ha='right',color=color,fontsize=7)
        ax.axhline(1,color='#bbbbbb',lw=.7,ls='--');ax.set(xlim=(8,41),ylim=lim,xticks=[9,19,29,39],xlabel='Interval length, L (native)',ylabel=ylabel)
        ax.set_title(title,loc='left',fontsize=9,pad=27)
        subtitle='Equal local cue; initial rate gain' if index==1 else '0.5% inactive withheld for 0.1\nthen exact return'
        ax.text(.5,1.02,subtitle,transform=ax.transAxes,ha='center',va='bottom',fontsize=6.3)
        if index==1:ax.legend(frameon=False,fontsize=7,loc='lower left',handlelength=1.7)
        else:
            gp=[next(r for r in functions if r['model']=='goryachev' and r['length_native']==L) for L in (19,39)]
            cf=gp[1]['weighted_conversion_coefficient_native']/gp[0]['weighted_conversion_coefficient_native']
            ax.text(.04,.94,f'Goryachev conversion: {cf:.2f}× (19→39)',transform=ax.transAxes,ha='left',va='top',fontsize=6.2,color=COLORS['goryachev'])
    fig.text(.52,.485,'Spatial interaction as a modifier · separate Holmes ring tests',ha='center',fontsize=9,weight='bold')
    fig.text(.325,.445,'D→E: change size; d = 12 fixed',ha='center',fontsize=7)
    fig.text(.70,.445,'E→F: same C/state; move cue',ha='center',fontsize=7)
    norm=Normalize(0,np.ceil(max(r['active_rac'].max() for r in rings)*10)/10)
    for col,c in enumerate(rings):
        ax=axes[1,col];x,o=c['display'],c['order'];C,d=c['C'],c['d']
        ax.pcolormesh(x[o],c['time_after_cue_onset'],c['active_rac'][:,o],norm=norm,cmap='viridis',shading='auto',rasterized=True)
        ax.axvline(0,color='white',ls=':',lw=.7);ax.axvline(d,color='#f4c66a',ls=':',lw=.7);ax.axhline(1000,color='white',ls='--',lw=.7)
        ax.text(-9.3,1070,'Cue off',color='white',fontsize=6.5,va='bottom')
        ax.set_title(f"{chr(68+col)}   C={C}, d={d} µm: {c['outcome']}",loc='left',fontsize=8,pad=7)
        ax.set(xlim=(-10,30),ylim=(0,4000),xticks=[-10,0,10,20,30],yticks=[0,1000,2000,3000,4000],xlabel='Arc from original incumbent (µm)')
        if col==0:ax.set_ylabel('Time after cue onset (s)')
        else:ax.set_yticklabels([])
        ax.set_facecolor('#eeeeee')
        if C==30:ax.axvspan(20,30,color='#eeeeee',zorder=0);ax.text(25,2300,'Outside\nC=30',ha='center',va='center',fontsize=6.5,color='#666666')
        for spine in ax.spines.values():spine.set_visible(True)
    box=axes[1,2].get_position();cax=fig.add_axes([.935,box.y0,.012,box.height]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax);cb.ax.tick_params(labelsize=7);cb.set_label('Active Rac (µM)',fontsize=8)
    for ext in ('png','pdf','svg'):
        p=dev_output_path(OUT,'integrated_framework_candidate.'+ext);fig.savefig(p,dpi=270,facecolor='white');FILES.append(p)
    plt.close(fig)

def main():
    source(Path(__file__));source(Path(__file__).with_name('length_dev_output_layout.py'))
    metrics,folds,snapshots,basalchecks=collect_basal();functions,checks,max_error=collect_functions(metrics);factor_rows,factor_folds=matched_factorization(functions,folds);rings=collect_rings()
    review=WORKSPACE/'current manuscript figs/reviews/integrated_cohesion_r1.md';source(review)
    final=read(BASE/'revision_20261014_draft_alignment/equation_typesetting/final_audit.json');source(BASE/'length_mechanisms.pdf',final['pdf_sha256'])
    csvwrite('all_five_incumbent_size_folds.csv',folds);csvwrite('all_five_native_basal_profiles.csv',snapshots);csvwrite('conditional_functional_readouts.csv',functions)
    csvwrite('matched_L19_L39_factorization.csv',factor_rows);csvwrite('matched_L19_L39_folds.csv',factor_folds)
    jsonwrite('matched_L19_L39_factorization.json',dict(definition='At exact return, weighted conversion=integral_core A(up)*delta_v / integral_core delta_v; exposure=integral_core delta_v / core active deficit; substrate term=conversion*exposure. Net repair also includes the active-field and native diffusion terms.',native_state_rows=factor_rows,folds=factor_folds,source_scope='Four existing pulse/sham returned-state pairs; no integration or independently evolved ablation.'))
    draw(folds,functions,rings)
    caption='Integrated provisional framework. A: paired-size basal incumbents across all five models; separate symbols show excess-active amount, boundary-connected half-contrast width and absolute peak concentration, each divided by its smaller-size value. These are separate observables, not a width-times-peak decomposition. Pairs are19→39 for Wave-pinning/Goryachev/corrected deBelly,20→40 for Holmes3 and resting L₀1→2 for Spring (physical domain1.048990→2.049367). They are roughly twofold nominal enlargements, without WGD calibration. More active material is the common observation; organization differs, and Holmes peak concentration decreases slightly. B–C: model-specific consequences in separate retained Goryachev/Wave-pinning interval states at all18 lengths9–39. Both readouts are divided by their own model’s L19 value, matching A’s19→39 Goryachev/Wave-pinning comparison; all36 L9–39 cases and prior L9 ratios remain in the raw CSV. B: initial added opposite-end activation per unit equal local cue equals native remote inactive concentration. C: actual immediate net fixed-core repair rate per active damage, including native diffusion, after0.5% of the whole inactive profile is withheld for0.1 native time unit and returned exactly pointwise. Positive Goryachev chemical repair is mediated by the endogenous inactive surplus; its active-field contribution at sham substrate is slightly negative. For matched19→39, Goryachev peak concentration grows1.9033×, surplus-weighted conversion coefficient3.5644× and core surplus exposure per damage0.9548×, giving3.4032× substrate contribution and3.4144× actual net repair. Wave-pinning’s corresponding conversion/exposure/product/net folds are1.0246/1.0175/1.0425/1.0630. These are properties of coupled returned states, not an independently evolved causal partition. This is conditional recovery, not global or intrinsic stability. D–F: separate matched-h=0.025-µm Holmes ring tests show raw active Rac through1000 s cue-on plus3000 s cue-free; all challengers activate, followed by replacement, merger or two retained sites. D→E changes circumference while keeping physical d12; only E→F holds circumference and exact basal state fixed while moving the cue. Changing size also changes total material, incumbent organization, complementary route and local state. The cue has total support3µm,peak1µM/s. White/gold dotted lines mark original/cue sites; labelled dashed lines mark withdrawal. Periodic edges connect at−10/20 for C30 and−10/30 for C40; gray is outside C30’s displayed interval. Full14000-s histories and endpoints, C50 and other models remain supporting. The bottom experiments are not the preparations used in A–C. No universal activation suppression/recovery effect, branch stability, causal partition or migration-persistence mechanism is established.'
    p=dev_output_path(OUT,'caption.txt');p.write_text(caption+'\n');FILES.append(p)
    verification=dict(passed=True,new_integrations=0,all_five_checkpoint_checks=basalchecks,conditional_cases=checks,conditional_case_count=len(functions),max_net_repair_reconstruction_error=max_error,
        normalization='A each larger/smaller pair; B/C within model relative to exact L19. All36 native values and prior L9 reference/ratio columns retained in CSV.',ring_matched_dx=.025,C40_ring_basals_identical=True,pdf_sha256=final['pdf_sha256'],source_hash_count=len(SOURCES))
    jsonwrite('verification.json',verification)
    mapping=[dict(panel='A',existing_PDF='Figure2A–C',change='Five paired-size folds of separate amount/width/absolute peak observables; raw whole-length evidence remains linked.'),dict(panel='B',existing_PDF='Figure3A/D chemical substrate evidence',change='Explicit native initial local-cue derivative at remote endpoint, all36 cases, matched within-model L19 normalization; prior L9 ratios retained.'),dict(panel='C',existing_PDF='Figure4B/D reaction-only rebound',change='Actual net reaction-plus-diffusion repair per core damage and exact substrate-surplus qualification; conditional return protocol.'),dict(panel='D–F',existing_PDF='No direct counterpart; new reviewed periodic-domain pilot',change='Spatial interaction modifier, not a replacement headline or extension of interval boundary geometry.')]
    jsonwrite('pdf_to_integrated_mapping.json',mapping)
    gaps=['A amounts are one-dimensional concentration integrals and are not normalized by domain length; larger amount alone does not isolate a size-specific feedback.','A is descriptive across five models; B/C are model-specific exemplars, not demonstrated for all five.','Excess amount, half-contrast width and absolute peak concentration are not an exact multiplicative decomposition.','Interval basal evidence and periodic challenge tests use different geometry/preparations; no paired causal path is inferred between them.','Only selected C40 ring pair isolates cue position; circumference changes also alter material, state and both routes.','D–F show finite histories, not stable coexistence, migration persistence or biological WGD calibration.','Full ring endpoints/late histories and predefined C50 outcome remain linked; precise merger time is mesh-sensitive.']
    jsonwrite('evidence_gaps.json',gaps)
    rel=lambda p:os.path.relpath(p,DISCUSSION.parent)
    table='\n'.join(f"| {r['model_label']} | {'L₀ ' if r['model']=='spring' else ''}{r['nominal_size_small']}→{r['nominal_size_large']} | {r['excess_active_fold']:.4f} | {r['width_fold']:.4f} | {r['peak_concentration_fold']:.4f} |" for r in folds)
    fends={m:next(r for r in functions if r['model']==m and r['length_native']==39) for m in COLORS}
    md=f'''# Integrated framework: enlargement changes the incumbent; separation modifies interaction

This is a provisional cohesive figure draft using existing evidence only. The previous [chemical candidate](main_candidate_1_20261017.md) and [spatial candidate](main_candidate_2_20261017.md) remain unchanged as supporting ingredients. No new simulation, fitting or duration optimization was performed.

![Integrated enlargement and interaction framework]({rel(OUT/'integrated_framework_candidate.png')})

## Main message and why both rows are needed

Across the five retained model preparations, roughly doubling nominal size increases the excess active material in the incumbent, with model-specific changes in width and peak concentration. A larger incumbent can therefore present a different competition problem. The functional consequences are demonstrated conditionally in Goryachev and Wave-pinning: changed remote substrate availability affects initial activation, and the disturbed incumbent can recapture redistributed substrate after exact return. The later interaction also depends on spatial arrangement, as the separate Holmes tests show. Separation is a modifier of this framework, not an alternative headline. All-five incumbent growth is established in these retained preparations; generalized greater stability is not. The spatial results cannot fill that gap. Organizing these ingredients into one figure does not prove that the causal mechanism is complete.

The top row makes the hierarchy explicit. **A is the common descriptive observation across all five models. B and C are model-specific functional exemplars**, not claims proved across the five. Marker shapes in A distinguish amount, width and peak concentration; Goryachev/Wave-pinning colors in B/C identify models. All folds retain raw values and exact references in the linked data. The bottom row preserves native spatial histories because activation or an original-core-loss clock cannot distinguish replacement, merger and retention of separate sites.

## A: increased active material, different organization

| Model | Nominal size pair | Excess-active amount fold | Half-contrast width fold | Absolute peak concentration fold |
|---|---|---:|---:|---:|
{table}

Spring uses resting length L₀1→2; its measured physical domain changes1.048990→2.049367, so amounts and widths use that physical domain. The other pairs are19→39 or20→40 in their native length units. These are convenient near-doubling comparisons, not a biological WGD fold calibration. Holmes peak concentration decreases slightly; “larger incumbent” means more active material or increased spatial scale, not universally higher peak amplitude.

Excess-active amount is the concentration integral above the field minimum inside the boundary-connected half-contrast core. Width and absolute peak concentration are separate observables, **not an exact width×peak factorization**; background changes matter, especially in deBelly. The amount is not divided by domain length, so its increase alone does not isolate a size-specific feedback mechanism. Full retained length series and native profiles remain in the existing PDF Figure2 and linked source tables.

## B–C: conditional effects on competition

All36 retained interval cases (18 per model, L9–39) appear. Both panels normalize each model to its own exact L19 readout, matching A’s19→39 chemical-model pair; native values and references remain in [the functional CSV]({rel(OUT/'supporting/data/conditional_functional_readouts.csv')}). Open points mark L19/L39 to connect with A.

B measures initial extra activation per unit equal local opposing cue. Native kinetics give f(u,v;s)−f(u,v;0)=s·v, so this readout is the remote inactive concentration. By L39, relative gain is {fends['goryachev']['cue_gain_relative_to_L19']:.4f} in Goryachev and {fends['wave_pinning']['cue_gain_relative_to_L19']:.4f} in Wave-pinning. It is an initial chemical response, not a front-establishment threshold or an evolved takeover probability.

C measures the **actual immediate net repair rate per active material lost**, including native reaction and diffusion in the original fixed core, relative to a matched sham. The retained experiment removes0.5% of the entire original inactive profile for0.1 native time unit and returns exactly that profile pointwise. By L39, this conditional gain is {fends['goryachev']['net_repair_relative_to_L19']:.2f}× its L19 value for Goryachev and {fends['wave_pinning']['net_repair_relative_to_L19']:.2f}× for Wave-pinning. Goryachev’s positive chemical rebound is mediated by endogenous inactive surplus; its active-field contribution at sham inactive concentration is slightly negative. This is conditional recapture after the specified disturbance, **not global or generic intrinsic stability**. The earlier chemical candidate retains actual damage/recovery histories and the detailed returned-state allocation. All raw native values and the original L9 reference/ratio columns remain unchanged in support; this revision changes normalization, not the experiment.

The visible conversion annotation links organization to chemistry over the **same19→39 pair**. Goryachev’s peak concentration increases1.9033-fold, and the surplus-weighted native conversion coefficient increases3.5644-fold. Core inactive surplus per core damage changes0.9548-fold, giving3.4032-fold substrate contribution; actual net repair changes3.4144-fold. Wave-pinning’s nearly unchanged peak concentration (1.00024-fold) accompanies1.0246-fold weighted conversion,1.0175-fold exposure,1.0425-fold substrate contribution and1.0630-fold net repair. The [matched factorization rows]({rel(OUT/'supporting/data/matched_L19_L39_factorization.csv')}), [folds]({rel(OUT/'supporting/data/matched_L19_L39_folds.csv')}) and [definitions/source JSON]({rel(OUT/'supporting/metadata/matched_L19_L39_factorization.json')}) retain full precision.

“Conversion” is the native coefficient A(u)=∂f/∂v, averaged with the actual positive returned inactive surplus as weight. This is not a new binding process or independently evolved ablation. The coefficient×exposure identity describes the substrate term; the net repair also includes active-field and diffusion changes. A broader incumbent can remain functionally important even if its local coefficient changes little. More material by itself is the common state observation, not the feedback mechanism.

## D–F: spatial arrangement modifies later interaction

These are separate native Holmes periodic-ring preparations at matched h=0.025µm, not the interval checkpoints in A–C. All three challenged regions activate. Through the displayed1000-s cue-on plus3000-s cue-free interval, C30,d12 loses the original site while the challenger remains; C40,d12 fills the intervening valley into a merged region; C40,d16 retains two separated sites. White/gold dotted lines identify original/cue positions and “Cue off” labels mark withdrawal at1000s. Concentration is raw µM on one shared color scale. Common old-at-zero coordinates retain the30-µm versus40-µm domain difference; gray space is outside C30’s displayed periodic interval.

D→E changes circumference while holding physical cue separation12µm. It also changes total material, incumbent width, local state, cue fraction and the complementary route (18→28µm). E→F alone changes cue position at fixed C40 and the **exact same full basal state**, with routes12/28→16/24µm. These interventions do not isolate diffusion travel time or establish incumbent width as the sole cause. The common cue has support3µm, peak1µM/s and duration1000s, followed by13000s without cue.

The [full spatial candidate and endpoints](main_candidate_2_20261017.md), [size-path report]({rel(BASE/'dev/holmes_size_bridge_20261017/supporting/metadata/RESULTS.md')}), [mesh checks]({rel(BASE/'dev/holmes_distance_robustness_20261017/supporting/metadata/MESH_RESULTS.md')}) and [N800 conserved-Rac imbalance continuation]({rel(BASE/'dev/holmes_distance_robustness_20261017/supporting/metadata/IMBALANCE_RESULTS.md')}) remain support. The latter uses h=0.05µm endpoint states, not the fine-grid states plotted here. The predefined C50,d20 case also retains two sites through the observed window; it is disclosed in the size report. No stable branch, universal distance law or precise merger timing is claimed.

## Relation to the current PDF and remaining gaps

A condenses the common material/organization evidence from current PDF Figure2. B makes the initial cue-gain implication of Figure3’s inactive-depletion evidence explicit. C uses the same removal/return evidence as Figure4 but plots actual net repair including diffusion, while retaining its substrate-surplus interpretation. D–F add the separately reviewed circular-domain tests; they have no direct PDF counterpart. [Exact panel mapping]({rel(OUT/'supporting/metadata/pdf_to_integrated_mapping.json')}) records these distinctions.

The figure offers one framework with different evidence levels: common incumbent enlargement; demonstrated model-specific chemical consequences; spatial modification of later interactions. It does not prove all-five suppression of competitor activation, all-five repair enhancement, a causal partition, migration persistence or a biological WGD calibration. Whether these conditional mechanisms explain the complete length response remains open. [Evidence gaps]({rel(OUT/'supporting/metadata/evidence_gaps.json')}) and the original ingredient documents retain the limitations.

## Sources and verification

- [All-five raw pairs/folds]({rel(OUT/'supporting/data/all_five_incumbent_size_folds.csv')}), [native basal profiles]({rel(OUT/'supporting/data/all_five_native_basal_profiles.csv')}) and [36 raw functional readouts]({rel(OUT/'supporting/data/conditional_functional_readouts.csv')}). Compact native ring NPZs beside these preserve every plotted heatmap sample.
- [Verification]({rel(OUT/'supporting/metadata/verification.json')}) reconstructs all ten A checkpoints and36 exact returned-state budgets, checks native cue identities and normalization, and confirms identical selected C40 basal arrays. [Source/output manifest]({rel(OUT/'supporting/metadata/manifest.json')}) identifies every relevant immutable input; source PDF identity is checked.
- [Caption]({rel(OUT/'supporting/metadata/caption.txt')}), [vector PDF]({rel(OUT/'supporting/vectors/integrated_framework_candidate.pdf')}) and [generator]({rel(Path(__file__))}). Existing candidates and retained sources are unchanged.
'''
    md+='\n## Fresh cohesion review and response\n\nThe [independent cohesion review](reviews/integrated_cohesion_r1.md) accepts this hierarchy for discussion while identifying two essential clarifications. This revision aligns B/C with A’s L19→39 chemical-model comparison, adds the visible conversion-coefficient link, and retains the common observation versus generalized-stability gap. The figure remains a provisional framework; the review does not establish a complete causal explanation. D–F and their separate experimental scope are unchanged.\n'
    # Human prose spacing, without changing formulas or data identifiers.
    substitutions={'Pairs are19':'Pairs are 19','pairs are19':'pairs are 19','models. B–C':'models. B–C','Holmes3':'Holmes 3','L₀1':'L₀ 1','domain1.048990':'domain 1.048990','Pairs are19':'Pairs are 19','pairs are19':'pairs are 19','or20':'or 20','PDF Figure2':'PDF Figure 2','Figure3':'Figure 3','Figure4':'Figure 4','All36':'All 36','all ten A checkpoints and36':'all ten A checkpoints and 36','removes0.5':'removes 0.5','for0.1':'for 0.1','displayed1000':'displayed 1000','plus3000':'plus 3000','at1000s':'at 1000 s','the30':'the 30','versus40':'versus 40','separation12µm':'separation 12 µm','routes12/28':'routes 12/28','support3µm':'support 3 µm','peak1µM/s':'peak 1 µM/s','duration1000s':'duration 1000 s','by13000s':'by 13,000 s','state, with routes12':'state, with routes 12','h=0.025µm':'h = 0.025 µm','h=0.05µm':'h = 0.05 µm','physical domain changes1':'physical domain changes 1','time unit and returns':'time unit and returns','A’s19':'A’s 19','matched19':'matched 19','same19':'same 19','grows1.9033':'grows 1.9033','coefficient3.5644':'coefficient 3.5644','damage0.9548':'damage 0.9548','giving3.4032':'giving 3.4032','and3.4144':'and 3.4144','are1.0246':'are 1.0246','increases1.9033':'increases 1.9033','increases3.5644':'increases 3.5644','changes0.9548':'changes 0.9548','giving3.4032':'giving 3.4032','changes3.4144':'changes 3.4144','accompanies1.0246':'accompanies 1.0246',',1.0175':', 1.0175',',1.0425':', 1.0425','and1.0630':'and 1.0630','and prior L9':'and prior L9','deBelly,20':'deBelly, 20','all18 lengths9–39':'all 18 lengths 9–39','after0.5%':'after 0.5%','through1000':'through 1000','µm,peak':'µm, peak','at−10':'at −10','and−10':'and −10','Full14000-s':'Full 14,000-s','all36':'all 36','L9':'L = 9','L19':'L = 19','L39':'L = 39','L9–39':'L = 9–39'}
    md=polish_prose(md,substitutions)
    for a,b in substitutions.items():caption=caption.replace(a,b)
    p=dev_output_path(OUT,'caption.txt');p.write_text(caption+'\n')
    final_prose_cleanups={'18→28µm': '18→28 µm', '16/24µm': '16/24 µm', 'C30,d12': 'C = 30, d = 12', 'C40,d12': 'C = 40, d = 12', 'C40,d16': 'C = 40, d = 16', 'C50,d20': 'C = 50, d = 20', 'fixed C40': 'fixed C = 40', 'for C30': 'for C = 30', 'for C40': 'for C = 40', 'C30’s': 'C = 30’s', 'C40 and': 'C = 40 and', 'matched-h=0.025-µm': 'h = 0.025 µm', 'physical d12': 'physical d = 12', 'all36': 'all 36', 'all18': 'all 18', 'lengths9–39': 'lengths 9–39', 'after0.5%': 'after 0.5%', 'through1000': 'through 1000', 'Full14000-s': 'Full 14,000-s'}
    md=polish_prose(md,final_prose_cleanups)
    for old,new in final_prose_cleanups.items():caption=caption.replace(old,new)
    md=md.replace('f(u,v;s)−f(u,v;0)=s·v',r'$f(u,v;s)-f(u,v;0)=s\,v$')
    p=dev_output_path(OUT,'caption.txt');p.write_text(caption+'\n')
    DISCUSSION.write_text(md)
    support_md=dev_output_path(OUT,'RESULTS.md');local=md
    for url in re.findall(r'\]\(([^)]+)\)',md):
        target=os.path.relpath(DISCUSSION.parent/url.strip('<>'),support_md.parent);target='<'+target+'>' if ' ' in target else target
        local=local.replace(']('+url+')',']('+target+')')
    support_md.write_text(local);FILES.append(support_md)
    jsonwrite('manifest.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',artifact_root='../..',status='Matched-pair cohesion revision; pending bounded final review',new_integrations=0,source_sha256=SOURCES,sha256={dev_artifact_key(OUT,p):sha(p) for p in FILES},discussion=str(DISCUSSION),reproduction='/home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))))
    print(json.dumps(dict(figure=str(OUT/'integrated_framework_candidate.png'),discussion=str(DISCUSSION),folds=folds,functional_cases=len(functions),source_hashes=len(SOURCES)),indent=2))
if __name__=='__main__':main()
