"""Audit spatial replays, measure competing clocks, and report all 15 cells.

Front/mass measurements and fixed-dose interventions are diagnostic; they do
not replace the independently audited shortest tested pulse durations.
"""
from pathlib import Path
import argparse, csv, html, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parents[2]
CODE=ROOT/'code';sys.path.insert(0,str(CODE))
import figure_builder as f
from three_test_figures import TESTS,LABELS
from datetime import datetime,timezone

COLORS=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']
OUT=ROOT/'data/runs/three_test_mechanisms_20261006'
FIG=ROOT/'figures/three_test_mechanisms'


def read(path):return json.loads(path.read_text())
def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def finite(value):return float(value) if value is not None and np.isfinite(value) else None
def fmt(value):return 'not observed' if value is None else f'{value:.4g}'
def savecsv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def polarity(u,model):
    q=np.linspace(0,1,u.shape[1]) if model=='spring' else (np.arange(u.shape[1])+.5)/u.shape[1]
    if model=='spring':return np.trapezoid(u*(2*q-1),q,axis=1)/np.trapezoid(u,q,axis=1)
    return (u@(2*q-1))/u.sum(axis=1)


def first(time,good):
    indices=np.flatnonzero(good)
    return float(time[indices[0]]) if len(indices) else None


def extent(u,q,threshold,length):
    """Boundary-connected high domain, with linear edge interpolation."""
    high=u>=threshold;counts=np.argmin(high,axis=1);counts[high.all(axis=1)]=u.shape[1]
    values=np.zeros(len(u));values[counts==u.shape[1]]=1.
    for i in np.flatnonzero((counts>0)&(counts<u.shape[1])):
        j=counts[i];den=u[i,j-1]-u[i,j]
        level=float(np.asarray(threshold[i]).reshape(-1)[0]) if np.ndim(threshold)>0 else threshold
        fraction=(u[i,j-1]-level)/den if den>0 else .5
        values[i]=q[j-1]+fraction*(q[j]-q[j-1])
    return values*length


def full_history(z):
    if 'pulse_time' not in z:
        return z['recovery_time'],z['recovery_active'],z['recovery_orientation'],0
    D=z['pulse_time'][-1]
    return np.r_[z['pulse_time'],D+z['recovery_time'][1:]],np.r_[z['pulse_active'],z['recovery_active'][1:]],np.r_[z['pulse_orientation'],z['recovery_orientation'][1:]],len(z['pulse_time'])-1


def measurements(z,job):
    c=job['case'];model=c['model'];L=c['length'];D=job['duration'];time,u,o,cut=full_history(z)
    q=np.linspace(0,1,u.shape[1]) if model=='spring' else (np.arange(u.shape[1])+.5)/u.shape[1]
    length=np.r_[z['pulse_field_length'],z['recovery_field_length'][1:]] if model=='spring' and 'pulse_time' in z else z['recovery_field_length'] if model=='spring' else np.full(len(u),L)
    weights=np.ones(len(q))/len(q)
    if model=='spring':weights[:]=1/(len(q)-1);weights[[0,-1]]*=.5
    mid=(q>=.475)&(q<=.525)
    background=np.median(u[:,mid],axis=1)
    excess=np.maximum(u-background[:,None],0);left=q<.5;right=q>.5
    leftmass=(excess[:,left]@weights[left])*length;rightmass=(excess[:,right]@weights[right])*length
    threshold=background+.5*(u[0].max()-background[0])
    width=extent(u,q,threshold[:,None],length);newwidth=extent(u[:,::-1],q,threshold[:,None],length)
    initialmass=float(leftmass[0]);initialwidth=float(width[0]);ratio=leftmass/initialmass if initialmass>0 else np.zeros(len(u))
    t50=first(time,ratio<=.5);t90=first(time,ratio<=.1);nuc=first(time,u[:,right].max(axis=1)>=threshold)
    frontmask=(width<=.8*initialwidth)&(width>=.2*initialwidth)&(time<=(first(time,width<=.2*initialwidth) or time[-1]))
    speed=None;rsq=None
    if frontmask.sum()>=4 and np.ptp(time[frontmask])>0:
        x=time[frontmask];y=width[frontmask];slope,intercept=np.linalg.lstsq(np.c_[x,np.ones(len(x))],y,rcond=None)[0]
        speed=-float(slope);rsq=1-float(np.sum((y-slope*x-intercept)**2))/float(np.sum((y-y.mean())**2)) if np.ptp(y)>0 else None
    pools=[k for k in z.files if k.startswith('pulse_field_') and ('inactive' in k or k=='pulse_field_mca_unbound')]
    poolgradient=max([float(np.max(np.ptp(z[k],axis=1)/(np.mean(abs(z[k]),axis=1)+1e-15))) for k in pools],default=None)
    tzero=first(time,o>=0)
    m=dict(slot=job['slot'],model=model,geometry=c['condition'],length=L,kind=job['kind'],duration=D,
           initial_old_pole_excess_amount=initialmass,initial_old_domain_width=initialwidth,
           initial_old_domain_fraction=initialwidth/L,initial_active_peak=float(u[0].max()),
           new_pole_activation_time=nuc,old_pole_half_removal_time=t50,old_pole_90pct_removal_time=t90,
           first_rightward_time=tzero,old_removal_front_speed=finite(speed),front_linear_r_squared=finite(rsq),
           reliable_front_speed=finite(speed) if rsq is not None and rsq>=.8 and speed>0 else None,
           old_amount_fraction_at_cutoff=float(ratio[cut]),new_amount_over_initial_old_at_cutoff=float(rightmass[cut]/initialmass),
           old_domain_fraction_at_cutoff=float(width[cut]/L),new_domain_fraction_at_cutoff=float(newwidth[cut]/L),
           maximum_relative_pool_gradient_during_pulse=poolgradient,
           right_absolute_half_height_time=first(time,u[:,right].max(axis=1)>=u[0].min()+.5*np.ptp(u[0])),
           orientation_at_cutoff=float(o[cut]),orientation_after_recovery=float(o[-1]),
           pulse_ends_before_old_half_removed=bool(t50 is None or D<t50))
    phase='pulse' if 'pulse_time' in z else 'recovery'
    inactive_name='rac_inactive' if model=='holmes_model3' else 'inactive'
    inactive=z[phase+'_field_'+inactive_name][0] if model!='debelly' else None
    m['initial_inactive_right']=float(inactive[-1]) if inactive is not None else None
    m['initial_inactive_mean']=float(inactive.mean()) if inactive is not None else None
    domain=float(length[0]);position=domain if model=='spring' else float(z[phase+'_x'][-1])
    if c['condition']=='fixed_slope':
        slope=c['gradient_slope']*((1. if model=='spring' else 15.)/L if job['kind']=='matched_gradient_peak' else 1.)
        signal=slope*position
    else:
        support=.2*domain if c['condition']=='fixed_relative' else c['width']
        signal=c['response_amplitude']*.5*(1+np.cos(np.pi*(domain-position)/support))
    if D==0:signal=0.
    m['initial_right_added_stimulus']=float(signal)
    source=signal*inactive[-1]/(7.5 if model=='holmes_model3' else 1.) if inactive is not None else .04*signal/(1+z[phase+'_field_rho'][0,-1]**2)
    m['initial_right_cue_activation_source']=float(source)
    # Average loss of the declared contrast amount, not a universal molecular flux.
    m['mean_old_excess_loss_rate_to_half']=initialmass/(2*t50) if t50 is not None and t50>0 else None
    if model=='debelly':
        b=np.r_[z['pulse_field_mca_bound'],z['recovery_field_mca_bound'][1:]] if 'pulse_time' in z else z['recovery_field_mca_bound']
        moment=(b@(2*q-1))/b.sum(axis=1)
        m.update(initial_mca_bound_orientation=float(moment[0]),mca_bound_orientation_at_cutoff=float(moment[cut]),mca_bound_orientation_after_recovery=float(moment[-1]))
    else:m.update(initial_mca_bound_orientation=None,mca_bound_orientation_at_cutoff=None,mca_bound_orientation_after_recovery=None)
    if model=='spring':
        m.update(initial_extension=float(length[0]-L),maximum_extension_during_pulse=float(np.max(length[:cut+1]-L)),extension_at_cutoff=float(length[cut]-L))
    else:m.update(initial_extension=None,maximum_extension_during_pulse=None,extension_at_cutoff=None)
    traces=dict(time=time,old_fraction=ratio,new_over_initial=rightmass/initialmass,old_width=width,new_width=newwidth,orientation=o,cut=cut)
    return m,traces


def audit(z,record):
    job=record['job'];c=job['case'];model=c['model'];metrics=record['metrics']
    original=ROOT/c['checkpoint']
    assert f.digest(original)==c['checkpoint_sha256']==record['checkpoint_sha256']
    with np.load(original) as prep:
        if model=='debelly':
            for name in [k for k in prep.files if k.startswith('state_')]:np.testing.assert_array_equal(z['originalstate_'+name[6:]],prep[name])
        else:np.testing.assert_array_equal(z['originalstate'],prep['state'])
        if job['kind']=='native':
            phase='pulse' if 'pulse_time' in z else 'recovery'
            for name in [k for k in prep.files if k.startswith('field_')]:
                if model=='debelly' and name[6:] not in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement','boundary_separation']:
                    continue  # Derived rates include the newly applied cue.
                key=phase+'_'+name
                if key in z:np.testing.assert_allclose(z[key][0],prep[name][-1],rtol=0,atol=1e-12)
    for phase in ['pulse','recovery']:
        if phase+'_time' not in z:continue
        u=z[phase+'_active'];assert np.isfinite(u).all()
        np.testing.assert_allclose(polarity(u,model),z[phase+'_orientation'],rtol=1e-9,atol=1e-11)
        np.testing.assert_allclose(np.ptp(u,axis=1),z[phase+'_contrast'],rtol=1e-10,atol=1e-12)
    time=z['recovery_time'];assert np.array_equal(time,np.arange(4001)*2.5)
    o=z['recovery_orientation'][-201:];q=z['recovery_contrast'][-201:]
    passed=bool(o.min()>.05 and q.min()>=.05 and np.ptp(o)<.005);error=float(np.max(abs(o/abs(c['initial_orientation'])-1)))
    assert passed==metrics['orientation_success'] and bool(passed and error<=.01)==metrics['full_polarity_success']
    if model=='debelly':
        names=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension']
        fields={name:z['recovery_field_'+name].copy() for name in names}
        fields['displacement']-=fields['displacement'].mean(axis=1)[:,None]
        fields['physical_length']=c['length']+z['recovery_field_boundary_separation']
    else:fields={key[len('recovery_field_'):]:z[key] for key in z.files if key.startswith('recovery_field_')}
    maximum=0.
    for name,values in fields.items():
        tail=values[-201:];scale=c['fixed_scales'][name]
        maximum=max(maximum,float(np.ptp(tail,axis=0).max())/scale,float(abs(tail[-1]-tail[-2]).max())/2.5*500/scale)
    assert bool(maximum<=.001)==metrics['field_tail_pass']
    np.testing.assert_allclose(maximum,metrics['field_tail_max_scaled'],rtol=1e-10,atol=1e-12)
    if job['kind']=='native':
        for key in ['orientation_success','full_polarity_success','field_tail_pass']:assert metrics[key]==job['expected'][key]
        np.testing.assert_allclose(float(o[-1]),job['expected']['final_orientation'],rtol=0,atol=1e-9)
    # Independent conservation audit follows each model's native convention.
    groups=[('rac_active','rac_inactive'),('rho_active','rho_inactive')] if model=='holmes_model3' else [('mca_bound','mca_unbound')] if model=='debelly' else [('active','inactive')]
    errors=[]
    for names in groups:
        for phase in ['pulse','recovery']:
            if phase+'_time' not in z:continue
            pool=sum(z[phase+'_field_'+name] for name in names)
            if model=='spring':mass=np.trapezoid(pool,np.linspace(0,1,pool.shape[1]),axis=1)*z[phase+'_field_length']
            else:mass=pool.sum(axis=1)
            errors.append(float(np.max(abs(mass/mass[0]-1))))
    assert max(errors)<1e-5,'Conservation error'
    return max(errors)


def ratios(rows,key):
    vals=[r[key] for r in sorted(rows,key=lambda r:r['length']) if r[key] is not None and r[key]>0]
    return vals[-1]/vals[0] if len(vals)>1 else None


def cell_summary(model,geometry,records,measurements,endpoints):
    native=[x for x in measurements if x['model']==model and x['geometry']==geometry and x['kind']=='native']
    passed=[x for x in native if 'orientation_pass' in records[x['slot']]['job']['labels']]
    if not passed:passed=native
    scalar=[x for x in endpoints if x['model']==model and x['condition']==geometry and x['definition']=='orientation']
    strength=[x for x in endpoints if x['model']==model and x['condition']==geometry and x['definition']=='full_polarity']
    controls=[x for x in measurements if x['model']==model and x['geometry']==geometry and x['kind']!='native']
    successful=[x for x in scalar if x['status']=='success']
    result=dict(model=model,geometry=geometry,n_success=len(successful),n_censored=sum(x['status']=='censored' for x in scalar),n_ineligible=sum(x['status']=='ineligible' for x in scalar),
                sampled_lengths=sorted(set(x['length'] for x in native)),initial_old_width_long_short_ratio=ratios(passed,'initial_old_domain_width'),
                initial_old_amount_long_short_ratio=ratios(passed,'initial_old_pole_excess_amount'),front_speed_long_short_ratio=ratios(passed,'reliable_front_speed'),
                old_half_time_long_short_ratio=ratios(passed,'old_pole_half_removal_time'),activation_time_long_short_ratio=ratios(passed,'right_absolute_half_height_time'),
                initial_inactive_right_long_short_ratio=ratios(passed,'initial_inactive_right'),initial_activation_source_long_short_ratio=ratios(passed,'initial_right_cue_activation_source'),
                old_excess_loss_rate_long_short_ratio=ratios(passed,'mean_old_excess_loss_rate_to_half'),
                pulse_duration_long_short_ratio=ratios(passed,'duration'),successful_pulses_ending_before_old_half_removed=sum(x['pulse_ends_before_old_half_removed'] for x in passed),
                number_of_native_selected_orientation_replays=len(passed))
    if successful:
        low=min(successful,key=lambda x:x['endpoint']);high=max(successful,key=lambda x:x['endpoint'])
        result.update(shortest_duration_length=low['length'],shortest_duration=low['endpoint'],longest_duration_length=high['length'],longest_duration=high['endpoint'])
    else:result.update(shortest_duration_length=None,shortest_duration=None,longest_duration_length=None,longest_duration=None)
    paired={(x['length'],x['definition']):x for x in scalar+strength}
    gate_ratios=[paired[x['length'],'full_polarity']['endpoint']/x['endpoint'] for x in successful if paired[x['length'],'full_polarity']['status']=='success']
    result['maximum_strength_orientation_duration_ratio']=max(gate_ratios,default=None)
    result['controls']={kind:[dict(length=x['length'],duration=x['duration'],orientation_success=records[x['slot']]['metrics']['orientation_success'],
                    final_orientation=x['orientation_after_recovery'],old_half_removal=x['old_pole_half_removal_time'],activation_time=x['new_pole_activation_time'],
                    old_fraction_at_cutoff=x['old_amount_fraction_at_cutoff']) for x in controls if x['kind']==kind] for kind in sorted(set(x['kind'] for x in controls))}
    return result


def interpret(cell):
    m,g=cell['model'],cell['geometry']
    common={
      'wave_pinning':('Expanding domains and competition for a conserved pool','At larger lengths the old active domain contains more material, while the initial inactive concentration stays nearly unchanged. The cue first creates a competitor; transport and competition after withdrawal determine whether it persists. Mobility can change commitment without changing the original cue.',
                       'Measure right activation, old-domain retreat and inactive gradients separately. Increasing inactive mobility should change long-cell cleanup or commitment more than the initial local cue source.'),
      'goryachev':('Mass sequestration in a concentrated self-amplifying pole','The incumbent pole stays about the same physical width but becomes much taller with length. Its inactive substrate at the opposite end falls sharply. The native added activation term is cue times inactive concentration, so equal patch peaks do not produce equal biochemical forcing. Slow competition and reconcentration of the conserved pool remain after cue initiation.',
                    'Measure old-pole height and opposite-end inactive abundance, rather than assuming a widening domain. Restoring substrate availability or using a separate fixed-total-material preparation should reduce the length dependence of initial patch activation; this proposed test is not part of the native comparison.'),
      'holmes_model3':('Spatial coverage plus coupled Rac/Rho domain competition','Initial Rac plateau height and inactive availability are almost length independent, but the Rac domain expands. Mutual inhibition with Rho and the conserved pools supports competing spatial domains. A localized right seed need not overwrite the incumbent, even when its right activity rises quickly.',
                         'Track both Rac and Rho boundaries and inactive-pool gradients. Broader support at the same peak should rescue spatial-coverage failures; faster inactive transport can improve gradient recovery without removing every localized-patch failure.'),
      'debelly':('Rewriting a spatial adhesion/mechanical imprint','The fast Rac profile can respond while MCA still carries the old state. Freezing MCA redistribution during the pulse makes every sampled native-successful dose fail in this geometry when native recovery resumes. Changing unbound-MCA diffusion often removes retained polarization rather than simply speeding reversal, which also shows that diffusion helps determine memory strength.',
                 'At matched Rac displacement, measure the bound/unbound MCA distribution and mechanical fields at cutoff. A Rac-only transient should relapse when the adhesion imprint is not rewritten; interventions on adhesion transport can change both the reversal threshold and the existence of persistent polarity.'),
      'spring':('A global extension-dependent gate coupled to polarity competition','The domain grows with length while the model links its extension to global deactivation and dilution. This is a global feedback state, distinct from deBelly adhesion memory. Cue-only removal of that feedback tests commitment at the native pulse and restores the original physical equations during recovery.',
                 'Measure physical length and extension during the cue together with Rac polarity. Changing mechanical restoring/extension feedback should shift the dose required for commitment, with effects depending on spatial cue geometry.')}
    title,text,prediction=common[m]
    if g=='fixed_physical':
        text+=' The fixed support covers a decreasing fraction of a longer cell. '
        if m=='wave_pinning':text+='The L=39 case is censored natively but passes after fourfold inactive mobility at the same maximum tested pulse; a single remeshing does not rescue it.'
        elif m=='goryachev':text+='The fixed K=64 amplitude rescues all eligible lengths, while the diminishing substrate still reduces initial added activation with length.'
        elif m=='holmes_model3':text+='At L=28 and 50, fourfold inactive mobility does not rescue the maximum native tested pulse; relative support is more effective than this mobility change.'
        elif m=='debelly':text+='The longest native passing dose changes outcome on a 401-point mesh, so the exact longest-length threshold remains resolution sensitive.'
        else:text+='Removing extension feedback defeats the shortest sampled cell but allows the two longer sampled cells to reverse at their native doses.'
    elif g=='fixed_relative':
        text+=' The cue keeps the same footprint fraction and peak as length grows. '
        if m=='wave_pinning':text+='All eligible lengths pass orientation, but the longer-cell full-strength threshold and field stability are distinct. Mobility improves the longest sampled orientation endpoint.'
        elif m=='goryachev':text+='All eligible lengths reverse at the fixed larger peak. Greater support helps without changing the strongly concentrated initial-pole regime.'
        elif m=='holmes_model3':text+='Orientation passes at all lengths, while many long cells remain outside the original-strength gate. This is a domain/shape and readout distinction, not absence of rightward response.'
        elif m=='debelly':text+='Cue-only MCA freezing fails at every sampled length, and faster unbound-MCA diffusion removes polarization at the shorter doses.'
        else:text+='The duration still increases despite constant footprint fraction. Extension feedback matters at the short sampled dose but is not indispensable at the two longer sampled doses.'
    else:
        text+=' The ramp peak increases with length because its physical slope is fixed. '
        if m in ['wave_pinning','holmes_model3']:
            text+='A matched-reference-peak ramp fails at the native long-cell doses, while fourfold inactive mobility preserves reversal and yields stronger recovered polarity. This supports opposing input-strength and pool/domain-cleanup contributions to the nonmonotonic curve.'
        elif m=='goryachev':text+='The growing ramp partly offsets the decline in initial inactive substrate, but removal/competition still slows. Matched-reference-peak ramps fail at the middle and longest native doses. Remeshing preserves the long-cell orientation classification at the selected orientation dose, while its numerical first moment changes.'
        elif m=='debelly':text+='All sampled matched-reference-peak doses still pass orientation, so this control does not identify a critical ramp-peak bottleneck here. MCA freezing fails throughout. The short arm is also conditioned by the scalar gate and relaxation near the preparation-eligibility boundary; native field flags remain explicit.'
        else:text+='Both matching the reference peak and removing extension feedback defeat every sampled dose above the shortest length, and feedback removal also defeats the shortest. Growing input and the global mechanical gate jointly support the flat/slightly decreasing native duration curve.'
    failures={key:sum(not x['orientation_success'] for x in values) for key,values in cell['controls'].items()}
    evidence=f'Old-width long/short={fmt(cell["initial_old_width_long_short_ratio"])}; old-excess amount={fmt(cell["initial_old_amount_long_short_ratio"])}; opposite-end inactive ratio={fmt(cell["initial_inactive_right_long_short_ratio"])}; initial added activation-source ratio={fmt(cell["initial_activation_source_long_short_ratio"])}. Control failures at native doses: '+', '.join(f'{key} {failures[key]}/{len(values)}' for key,values in cell['controls'].items())+'.'
    return dict(mechanism=title,interpretation=text,prediction=prediction,evidence=evidence,
                qualification='Trajectory and fixed-dose evidence in this native regime. Counterfactual thresholds were not fully searched; frozen-scale field flags and mesh sensitivity remain separate from scalar success.')


def text_page(pdf,title,paragraphs):
    import textwrap
    fig=plt.figure(figsize=(11,8.5));fig.text(.07,.94,title,fontsize=17,weight='bold')
    y=.865
    for heading,text in paragraphs:
        fig.text(.07,y,heading,fontsize=11,weight='bold');y-=.035
        lines=textwrap.wrap(text,width=122)
        fig.text(.07,y,'\n'.join(lines),fontsize=9,va='top',linespacing=1.5);y-=.025*len(lines)+.045
    pdf.savefig(fig);plt.close(fig)


def mechanical_plot(model,geometry,records,pdf):
    chosen=sorted([rec for rec in records.values() if rec['job']['case']['model']==model and rec['job']['case']['condition']==geometry and rec['job']['kind']=='native' and 'orientation_pass' in rec['job']['labels']],key=lambda rec:rec['job']['case']['length'])
    fig,axes=plt.subplots(1,3,figsize=(12.5,4.7));fig.suptitle(f'{f.LABELS[model]} · {LABELS[TESTS.index(geometry)]} · mechanical state during the pulse',fontsize=13)
    for i,rec in enumerate(chosen):
        job=rec['job'];L=job['case']['length'];D=job['duration'];color=COLORS[i%len(COLORS)]
        with np.load(OUT/'trajectories'/f'trajectory_{job["slot"]:04d}.npz') as z:
            time=z['pulse_time']/D;q=z['pulse_x']/L
            if model=='debelly':
                b=z['pulse_field_mca_bound'];moment=(b@(2*q-1))/b.sum(axis=1)
                axes[0].plot(time,z['pulse_orientation'],color=color,label=f'L={L:g} Rac');axes[0].plot(time,moment,'--',color=color,label=f'L={L:g} bound MCA')
                axes[1].plot(q,b[0],':',color=color);axes[1].plot(q,b[-1],color=color,label=f'L={L:g}')
                tension=z['pulse_field_tension'];axes[2].plot(time,np.ptp(tension,axis=1),color=color,label=f'L={L:g}')
            else:
                extension=z['pulse_field_length']-L
                axes[0].plot(time,extension,color=color,label=f'L0={L:g}')
                axes[1].plot(time,80*extension,color=color,label=f'L0={L:g}')
                axes[2].plot(80*extension,z['pulse_orientation'],color=color,label=f'L0={L:g}');axes[2].scatter(80*extension[-1],z['pulse_orientation'][-1],color=color,s=20)
    if model=='debelly':
        for ax,title,xlabel,ylabel in zip(axes,['A  Fast Rac and adhesion imprint','B  Bound MCA before / after cue','C  Mechanical tension range'],['Time / applied duration','Position x / L','Time / applied duration'],['Native first moment','Bound MCA','Spatial tension range']):ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        axes[0].legend(frameon=False,fontsize=6,ncol=2)
        footer='Solid: Rac or cutoff profile; dashed: bound-MCA moment; dotted: initial MCA profile. MCA and Rac can have opposite spatial moments. Frozen-MCA controls preserve the original adhesion fields during cue only.'
    else:
        for ax,title,xlabel,ylabel in zip(axes,['A  Domain extension','B  Global feedback contribution','C  Coupled chemical/mechanical state'],['Time / applied duration','Time / applied duration','80 × extension'],['Physical length − L0','Added global deactivation','Rac orientation']):ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        axes[0].legend(frameon=False,fontsize=7)
        footer='The native deactivation rate is 3 plus 80 times extension. This global feedback is distinct from local adhesion memory. The diagnostic removes only the extension-dependent term during cue, then restores native recovery.'
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=7)
    fig.text(.05,.025,footer,fontsize=7);fig.tight_layout(rect=[0,.10,1,.90]);name=f'{model}_{geometry}_mechanics'
    fig.savefig(FIG/(name+'.png'),dpi=180);pdf.savefig(fig);plt.close(fig)
    return name


def cell_plot(model,geometry,records,metrics,ends,pdf):
    native=[rec for rec in records.values() if rec['job']['case']['model']==model and rec['job']['case']['condition']==geometry and rec['job']['kind']=='native']
    passes=[rec for rec in native if 'orientation_pass' in rec['job']['labels'] or 'orientation_censored' in rec['job']['labels']]
    passes=sorted(passes,key=lambda x:x['job']['case']['length'])
    fig,axes=plt.subplots(2,3,figsize=(13,7.8));fig.suptitle(f'{f.LABELS[model]} · {LABELS[TESTS.index(geometry)]}',fontsize=15)
    ax=axes[0,0]
    for gate,color,marker,label in f.GATES:
        rs=sorted([x for x in ends if x['model']==model and x['condition']==geometry and x['definition']==gate],key=lambda x:x['length'])
        ax.plot([x['length'] for x in rs],[x['endpoint'] if x['status']=='success' else np.nan for x in rs],marker+'-',color=color,label=label,ms=3)
        for x in rs:
            if x['status']=='censored':ax.plot(x['length'],x['maximum_tested_duration'],'^',mfc='white',color=color)
    ax.set_yscale('log');ax.set_title('A  Audited finite-pulse durations');ax.set_xlabel('Native length');ax.set_ylabel('Application duration');ax.legend(fontsize=6,frameon=False)
    for j,rec in enumerate(passes):
        job=rec['job'];c=job['case'];L=c['length'];metric=metrics[job['slot']];color=COLORS[j%len(COLORS)]
        with np.load(OUT/'trajectories'/f'trajectory_{job["slot"]:04d}.npz') as z:
            m,tr=measurements(z,job);time=tr['time'];D=job['duration']
            axes[0,1].plot(z['pulse_x'] if 'pulse_x' in z else z['recovery_x'],z['pulse_active'][0] if 'pulse_active' in z else z['recovery_active'][0],color=color,label=f'L={L:g}')
            # Compare pulse-progress clocks without folding in recovery.
            take=time<=D
            axes[0,2].plot(time[take]/D,tr['old_fraction'][take],color=color,label=f'L={L:g}')
            axes[1,0].plot(time[take]/D,tr['new_over_initial'][take],color=color)
            rt=z['recovery_time'];ro=z['recovery_orientation'];axes[1,1].plot(rt,ro,color=color,label=f'L={L:g}')
        for kind,marker in [('native','o'),('mobility_4x','s'),('matched_gradient_peak','^'),('mca_frozen_during_pulse','D'),('feedback_off_during_pulse','v'),('mesh_401','x')]:
            candidates=[x for x in metrics.values() if x['model']==model and x['geometry']==geometry and x['length']==L and x['kind']==kind and (kind!='native' or x['slot']==job['slot'])]
            for x in candidates:axes[1,2].scatter(L,x['orientation_after_recovery'],marker=marker,color=color,s=30)
    axes[0,1].set_title('B  Initial coherent domain / peak');axes[0,1].set_xlabel('Native position');axes[0,1].set_ylabel('Initial active species');axes[0,1].legend(frameon=False,fontsize=7)
    axes[0,2].set_title('C  Old-pole excess during selected pulse');axes[0,2].set_xlabel('Time / applied duration');axes[0,2].set_ylabel('Old amount / original old amount');axes[0,2].axhline(.5,color='#aaaaaa',ls=':',lw=1)
    axes[1,0].set_title('D  New-pole excess during selected pulse');axes[1,0].set_xlabel('Time / applied duration');axes[1,0].set_ylabel('New amount / original old amount')
    axes[1,1].set_title('E  Native polarity after withdrawal');axes[1,1].set_xlabel('Basal recovery time');axes[1,1].set_ylabel('Orientation');axes[1,1].set_xscale('symlog',linthresh=1);axes[1,1].axhline(.05,color='#777777',ls=':',lw=1)
    axes[1,2].set_title('F  Fixed-dose causal diagnostics');axes[1,2].set_xlabel('Native length');axes[1,2].set_ylabel('Orientation after recovery');axes[1,2].axhline(.05,color='#777777',ls=':',lw=1)
    handles=[plt.Line2D([0],[0],marker=m,color='#555555',ls='none',label=l) for m,l in [('o','Native'),('s','Mobility ×4'),('^','Matched ramp peak'),('D','MCA frozen during cue'),('v','Spring feedback off during cue'),('x','401-point mesh')]]
    axes[1,2].legend(handles=handles,fontsize=6,frameon=False,loc='best')
    for ax in axes.flat:ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=7)
    fig.text(.07,.022,'Amounts use excess above the current mid-cell median in each half; domain edges use that background plus half the initial pole contrast. Controls are fixed-dose diagnostics.',fontsize=7)
    fig.tight_layout(rect=[0,.06,1,.94]);name=f'{model}_{geometry}_evidence';fig.savefig(FIG/(name+'.png'),dpi=180);pdf.savefig(fig);plt.close(fig)
    # Full spatial evolution, explicitly showing withdrawal and a separate recovery scale.
    fig,axes=plt.subplots(len(passes),3,figsize=(12,2.2*len(passes)+1.3),squeeze=False)
    for i,rec in enumerate(passes):
        job=rec['job'];L=job['case']['length'];D=job['duration'];color=COLORS[i%len(COLORS)]
        with np.load(OUT/'trajectories'/f'trajectory_{job["slot"]:04d}.npz') as z:
            if 'pulse_time' not in z:continue
            q=z['pulse_x']/L
            for phase,col in [('pulse',0),('recovery',1)]:
                u=z[phase+'_active'];time=z[phase+'_time'];indices=np.unique(np.linspace(0,len(time)-1,min(500,len(time))).astype(int))
                ax=axes[i,col]
                im=ax.pcolormesh(q,time[indices]/D if phase=='pulse' else time[indices],u[indices],shading='auto',cmap='viridis',rasterized=True)
                if phase=='recovery':ax.set_yscale('symlog',linthresh=2.5)
                ax.set_ylabel(f'L={L:g}\n'+('Time / pulse duration' if phase=='pulse' else 'Recovery time'));ax.set_xlabel('Position x / L');fig.colorbar(im,ax=ax,pad=.015).ax.tick_params(labelsize=6)
            ax=axes[i,2]
            for frac,ls in [(0,':'),(.5,'--'),(1,'-')]:
                idx=int(np.argmin(abs(z['pulse_time']-D*frac)));ax.plot(q,z['pulse_active'][idx],ls,color=color,label=f'{frac:g} × pulse')
            ax.plot(q,z['recovery_active'][-1],color='#333333',lw=1.2,label='Recovery endpoint');ax.set_xlabel('Position x / L');ax.set_ylabel('Active species');ax.legend(fontsize=6,frameon=False)
            if i==0:
                for ax,title in zip(axes[i],['Native pulse','Native basal recovery','Actual profile snapshots']):ax.set_title(title,fontsize=10)
    fig.suptitle(f'{f.LABELS[model]} · {LABELS[TESTS.index(geometry)]} · selected native spatial histories',fontsize=13)
    fig.text(.05,.012,'Native scalar-selected pulses. Separate panel color scales preserve concentration changes; these are not independently verified steady states.',fontsize=7)
    fig.tight_layout(rect=[0,.035,1,.96]);name2=f'{model}_{geometry}_fields';fig.savefig(FIG/(name2+'.png'),dpi=180);pdf.savefig(fig);plt.close(fig)
    return name,name2


def run(out):
    global OUT
    OUT=out;FIG.mkdir(parents=True,exist_ok=True)
    manifest=read(out/'manifest.json');campaign=ROOT/manifest['campaign'];status=read(campaign/'results/status.json');assert status['complete']
    records={};rows=[];audits=[]
    for job in read(out/'jobs.json'):
        path=out/'trajectories'/f'trajectory_{job["slot"]:04d}.npz';rec=read(path.with_suffix('.json'))
        assert rec['job']==job and f.digest(path)==rec['archive_sha256']
        with np.load(path) as z:
            conservation=audit(z,rec);m,_=measurements(z,job)
        rows.append(m);records[job['slot']]=rec;audits.append(dict(slot=job['slot'],kind=job['kind'],native_replay=job['kind']=='native',maximum_conservation_error=conservation,passed=True))
    ends=f.read_endpoint_rows(campaign/'results/endpoints.csv');metrics={x['slot']:x for x in rows}
    cells=[cell_summary(m,g,records,rows,ends) for m in f.MODELS for g in TESTS]
    for cell in cells:cell.update(interpret(cell))
    mesh_checks=[]
    for record in records.values():
        job=record['job']
        if job['kind']!='mesh_401':continue
        native=next(x for x in records.values() if x['job']['kind']=='native' and x['job']['case']['case_index']==job['case']['case_index'] and x['job']['duration']==job['duration'])
        mesh_checks.append(dict(model=job['case']['model'],geometry=job['case']['condition'],length=job['case']['length'],duration=job['duration'],
             orientation_classification_preserved=record['metrics']['orientation_success']==native['metrics']['orientation_success'],
             strength_classification_preserved=record['metrics']['full_polarity_success']==native['metrics']['full_polarity_success'],
             native_orientation=native['metrics']['final_orientation'],mesh_orientation=record['metrics']['final_orientation']))
    write(out/'results/mesh_diagnostic.json',dict(orientation_preserved=sum(x['orientation_classification_preserved'] for x in mesh_checks),
          strength_preserved=sum(x['strength_classification_preserved'] for x in mesh_checks),total=len(mesh_checks),
          caution='401-point interpolated initial states at native selected orientation doses (or native censor cap). These are not converged thresholds or native original-strength endpoint mesh tests.',checks=mesh_checks))
    savecsv(out/'results/trajectory_metrics.csv',rows);write(out/'results/cell_evidence.json',cells)
    write(out/'results/spatial_audit.json',dict(passed=True,complete=True,trajectories=len(rows),coverage_cells=15,
          scope='Every archive hash, original full checkpoint identity, native initial fields, profile-derived orientation/contrast, scalar gates, native replay outcomes, and conserved pools.',checks=audits))
    anchors=[x for x in rows if x['kind']=='native' and 'orientation_pass' in records[x['slot']]['job']['labels']]
    before=sum(x['pulse_ends_before_old_half_removed'] for x in anchors)
    # Keep exhaustive diagnostic plots separate from the figure-led reader PDF.
    pdfpath=FIG/'three_test_diagnostics.pdf'
    paragraphs=[('Coverage and assay',f'All five models and all three tests are covered. {len(rows)} native replays and diagnostic controls passed spatial audit. The comparison measures the shortest tested pulse that persists after 10,000 units of recovery. {before} of {len(anchors)} selected successful replay pulses end before half the old-pole contrast excess disappears. Pulse duration is a commitment dose and need not equal complete reorientation time.'),
       ('1. Expanding domains versus spatial cue coverage','Wave-pinning and Holmes keep nearly constant plateau height while the initial coherent domain expands with L. A fixed support covers a decreasing fraction; relative support rescues many long-cell failures without changing peak amplitude. Prediction: varying support at fixed peak should preferentially change large-cell commitment. This geometric burden does not by itself establish constant edge speed or a linear time law.'),
       ('2. Substrate starvation around a concentrated incumbent pole','Goryachev instead retains a narrow pole while its peak rises from about 6 at L=3 to 117 at L=39. Opposite-end inactive abundance falls from about 0.879 to 0.088. The cue acts through that substrate, so a constant patch peak supplies about tenfold less initial activation in the long cell. Prediction: changing substrate availability or the total-material preparation should alter this component; the present assay retains its exact original states.'),
       ('3. Stronger gradient initiation versus redistribution/competition','At fixed physical slope the peak grows with L. In Wave-pinning and Holmes, matched-reference-peak ramps fail at native long-cell doses, while increased inactive mobility yields stronger recovered polarity. The longest fixed-width Wave-pinning case is rescued by mobility at its unchanged maximum pulse. These controls support distinct input-strength and transport/competition contributions, rather than a single erosion clock.'),
       ('4. Spatial adhesion/mechanical memory','deBelly requires rewriting an additional state: cue-only MCA freezing defeats all twelve sampled native-successful orientation doses across the three tests. Fast Rac response alone therefore does not commit those reversals. Unbound-MCA mobility changes can erase polarization itself. Prediction: measure MCA and mechanics at cutoff alongside Rac, and test retention under adhesion perturbation.'),
       ('5. Global extension feedback gates commitment','Spring uses extension-dependent global deactivation and dilution, rather than local MCA memory. Removing extension feedback during the cue defeats all four sampled gradient doses, while some longer patch doses remain successful. Matching the gradient peak also defeats the longer native doses. Prediction: perturbing extension/restoring feedback should affect commitment differently for localized patches and whole-cell ramps.')]
    with PdfPages(pdfpath) as pdf:
        text_page(pdf,'Distinct, testable explanations across the three tests',paragraphs[:4])
        text_page(pdf,'Two distinct mechanical mechanisms',paragraphs[4:]+[
            ('Evidence is model- and test-specific','The following fifteen sections each give an interpretation, evidence and testable prediction. A mechanism supported in one geometry is not automatically assigned to another. Coupled Rac/Rho domains, local adhesion memory and a global extension gate are different state variables and require different measurements.'),
            ('Boundary sensitivity is visible',f'Longest-length remeshing preserves {sum(x["orientation_classification_preserved"] for x in mesh_checks)} of {len(mesh_checks)} orientation classifications and {sum(x["strength_classification_preserved"] for x in mesh_checks)} strength classifications at the native selected orientation doses or censor caps. deBelly fixed-width orientation changes. These are fixed-dose checks, not converged thresholds or mesh tests at every original-strength endpoint. Native scalar gates and full-field steady-state claims remain separate.')])
        text_page(pdf,'Tests of length mechanisms and practical limits',[
            ('The suggested constant-speed erosion mechanism','It predicts an old-domain width that grows with L, similar physical edge speeds at different lengths, and a removal time explained by width/speed. The report measures all three independently. It must not be inferred solely from a positive duration slope. Narrow poles that increase peak height, pulses ending before removal, and major mobility effects require a different account.'),
            ('Separate the stimulus from the internal mechanism','A constant-slope ramp has a peak proportional to length. The matched-peak control deliberately changes the slope only as a diagnostic, outside the native three-test comparison. A reversal of the activation trend under that control supports input-strength acceleration. Fixed and relative patch amplitudes never scale with length.'),
            ('Shape/readout effects can condition apparent length trends','The strength gate uses a normalized first moment, whereas pole width, peak height and total activity measure shape and material. Adjacent passing/failing native pulses and both scalar thresholds are retained. A minimum in the curve can combine faster initiation at larger L with slower cleanup, mechanical rewriting or selective scalar stability. Recovery time is not added to the plotted pulse duration.'),
            ('What fixed-dose controls establish','Interventions test trajectories at the original selected pulse, not a fully refined counterfactual threshold. Mobility changes can alter an equilibrium profile and the 1% strength readout as well as transport. A 401-point remeshing diagnostic checks the longest selected cell; it interpolates rather than independently re-establishing the start. The report preserves outcome changes instead of treating diagnostic completion as robustness.'),
            ('Measurements and evidence','Pole amount integrates positive excess above the contemporaneous median in the central 5% of the cell, separately in each half. This removes broadly elevated activity under gradients. The boundary-connected domain threshold is that background plus half the initial pole contrast. Front fits use the first 80%-to-20% erosion segment, need four samples, and enter speed ratios only with R-squared at least 0.8. First crossings are sampled; absent crossings remain explicit. Raw fits and all quantities are in the CSV.')])
        images=[]
        for cell in cells:
            text_page(pdf,f'{f.LABELS[cell["model"]]} · {LABELS[TESTS.index(cell["geometry"])]}',[
                (cell['mechanism'],cell['interpretation']),('Measured evidence',cell['evidence']),('Testable prediction',cell['prediction']),('Scope of the inference',cell['qualification'])])
            names=cell_plot(cell['model'],cell['geometry'],records,metrics,ends,pdf)
            mechanics=mechanical_plot(cell['model'],cell['geometry'],records,pdf) if cell['model'] in ['debelly','spring'] else None
            images.append((cell,names,mechanics))
    stamp=datetime.now(timezone.utc).isoformat()
    # A table across every cell prevents the mechanism report narrowing to one model/test.
    table=[];sections=[]
    for cell,(name,name2),mechanics in images:
        controls='; '.join(f'{key}: '+', '.join(f'L={x["length"]:g}, P={x["final_orientation"]:.3g}' for x in values) for key,values in cell['controls'].items())
        table.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [f.LABELS[cell['model']],LABELS[TESTS.index(cell['geometry'])],f'{cell["n_success"]} / {cell["n_success"]+cell["n_censored"]}',fmt(cell['initial_old_width_long_short_ratio']),fmt(cell['initial_old_amount_long_short_ratio']),fmt(cell['front_speed_long_short_ratio']),fmt(cell['old_half_time_long_short_ratio']),fmt(cell['activation_time_long_short_ratio'])])+'</tr>')
        mechanical_html=f'<details><summary>Mechanical state during the pulse</summary><img loading="lazy" src="{mechanics}.png"></details>' if mechanics else ''
        sections.append(f'<section id="{name}"><h2>{f.LABELS[cell["model"]]} · {LABELS[TESTS.index(cell["geometry"])]}</h2><h3>{html.escape(cell["mechanism"])}</h3><p>{html.escape(cell["interpretation"])}</p><p><b>Prediction:</b> {html.escape(cell["prediction"])}</p><p>{html.escape(cell["evidence"])}</p><p>{cell["n_success"]} successful eligible lengths; {cell["n_censored"]} censored; {cell["n_ineligible"]} ineligible. Native replay lengths: {cell["sampled_lengths"]}.</p><p>Old width long/short: {fmt(cell["initial_old_width_long_short_ratio"])}; old amount: {fmt(cell["initial_old_amount_long_short_ratio"])}; fitted removal-speed ratio: {fmt(cell["front_speed_long_short_ratio"])}. {cell["successful_pulses_ending_before_old_half_removed"]}/{cell["number_of_native_selected_orientation_replays"]} selected pulses end before half-removal.</p><p>{html.escape(controls)}</p><p>{html.escape(cell["qualification"])}</p><img loading="lazy" src="{name}.png"><details><summary>Full native spatial histories</summary><img loading="lazy" src="{name2}.png"></details>{mechanical_html}</section>')
    links=' · '.join(f'<a href="#{cell["model"]}_{cell["geometry"]}_evidence">{f.LABELS[cell["model"]]} / {LABELS[TESTS.index(cell["geometry"])]}</a>' for cell in cells)
    body=''.join(f'<h2>{heading}</h2><p>{text}</p>' for heading,text in paragraphs)
    (FIG/'index.html').write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Five-model three-test mechanisms</title>
<style>body{{max-width:1320px;margin:30px auto;padding:0 25px;font:16px/1.6 system-ui;color:#20252b}}img{{max-width:100%}}table{{border-collapse:collapse;font-size:13px}}td,th{{padding:8px;border:1px solid #ddd;text-align:left}}section{{border-top:1px solid #ddd;margin-top:40px}}a{{color:#12628c}}</style>
<h1>Mechanisms across all five models × three tests</h1><p>Data generated {stamp}; producer finished. <a href="three_test_mechanisms.pdf">Complete PDF</a> · <a href="../../data/runs/three_test_mechanisms_20261006/results/trajectory_metrics.csv">All measurements</a> · <a href="../../data/runs/three_test_mechanisms_20261006/results/cell_evidence.json">Cell evidence and controls</a> · <a href="../../data/runs/three_test_mechanisms_20261006/results/spatial_audit.json">Independent spatial audit</a></p>{body}
<p>{links}</p><h2>Every cell, with independently measured quantities</h2><p>Ratios compare the longest and shortest available native selected orientation replays. Absent or undersampled front crossings remain unavailable. A positive pulse-duration slope alone does not establish constant-speed removal.</p>
<table><tr><th>Model</th><th>Test</th><th>Success / eligible</th><th>Old width ratio</th><th>Old amount ratio</th><th>Front-speed ratio</th><th>Half-removal time ratio</th><th>Right absolute half-height ratio</th></tr>{''.join(table)}</table>{''.join(sections)}</html>''')
    write(out/'results/status.json',dict(generated_at=stamp,complete=True,producer='finished',trajectories=len(rows),coverage_cells=15,native_replays=sum(x['kind']=='native' for x in rows),spatial_audit_passed=True))
    write(out/'results/report_provenance.json',dict(generated_at=stamp,producer='finished',campaign_endpoint_sha256=f.digest(campaign/'results/endpoints.csv'),renderer_sha256=f.digest(Path(__file__)),spatial_audit_sha256=f.digest(out/'results/spatial_audit.json'),pdf_sha256=f.digest(pdfpath),html_sha256=f.digest(FIG/'index.html')))
    from three_test_report import render
    render(out)
    import build
    build.index()
    print(json.dumps(dict(complete=True,trajectories=len(rows),coverage_cells=15,
                         pdf=str(FIG/'three_test_mechanisms.pdf'),diagnostic_pdf=str(pdfpath)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--report-only',action='store_true',help='Rebuild the reader report from completed audited evidence, without rerunning analysis.')
    args=p.parse_args()
    if args.report_only:
        from three_test_report import render
        render(args.out.resolve())
    else:
        run(args.out.resolve())
