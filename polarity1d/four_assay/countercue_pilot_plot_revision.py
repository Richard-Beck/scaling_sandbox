"""Cue-design illustration and mean/slope length curves from completed pilot data."""
from pathlib import Path
import gzip,json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import countercue_pilot as p
import runner as r

OUT=p.RESULTS
MODELS=r.MODELS
LABELS=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
MEANS=[.5,1.,2.]
SLOPES=[.25,.5,1.]
COLORS={.5:'#0072B2',1.:'#D55E00',2.:'#009E73'}
STYLES={.25:':',.5:'--',1.:'-'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=180)
 plt.close(fig)

def cue_design(cases,stamp):
 fig,axes=plt.subplots(3,2,figsize=(13,10),sharex=True,sharey=True)
 model='wave_pinning';fraction=.2
 reference=next(c['reference_domain_length'] for c in cases if c['model']==model)
 long=max(p.LENGTHS[model]);xmax=1.13*fraction*long/reference
 checks=[]
 for i,L in enumerate(p.LENGTHS[model]):
  fraction_case=next(c for c in cases if c['model']==model and c['length']==L and c['family']=='width' and c['fraction']==fraction and c['width_mode']=='fractional')
  for j,dose in enumerate(['peak_matched','integral_matched']):
   ax=axes[i,j]
   physical_case=next(c for c in cases if c['model']==model and c['length']==L and c['family']=='width' and c['fraction']==fraction and c['width_mode']=='physical' and c['dose_mode']==dose)
   distances=np.linspace(0,min(xmax*reference,L),20001)
   profiles=[]
   for c,color,ls,label in [(fraction_case,'#176d9c','--','Fractional width: 0.2*L'),(physical_case,'#db7b24','-','Physical width: 0.2*Lref')]:
    profile=p.Cue(c,1)(0,L-distances,L)/c['A']
    ax.fill_between(distances/reference,0,profile,color=color,alpha=.13)
    ax.plot(distances/reference,profile,color=color,ls=ls,lw=2,label=label)
    width=fraction*(L if c['width_mode']=='fractional' else reference)
    peak=1 if c['dose_mode']=='peak_matched' else L/reference
    expected_area=.5*peak*width/reference
    area=float(np.trapezoid(profile,distances/reference))
    np.testing.assert_allclose(area,expected_area,rtol=1e-7,atol=1e-10)
    profiles.append((peak,area,width/reference))
   fp,fa,fw=profiles[0];pp,pa,pw=profiles[1]
   if dose=='integral_matched':np.testing.assert_allclose(fa,pa,rtol=1e-7,atol=1e-10)
   checks.append(dict(length=L,dose_mode=dose,fractional_peak_over_A=fp,physical_peak_over_A=pp,fractional_area_over_A_Lref=fa,physical_area_over_A_Lref=pa))
   ax.text(.97,.94,f'Fractional: peak {fp:.3g} A; area {fa:.3g}\nPhysical: peak {pp:.3g} A; area {pa:.3g}',transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
   ax.axvline(fw,color='#176d9c',ls='--',lw=.8,alpha=.5)
   ax.axvline(pw,color='#db7b24',lw=.8,alpha=.5)
   ax.set_ylim(0,2.2);ax.set_xlim(0,xmax);ax.set_ylabel(f'L = {L:g} ({L/reference:.3g}*Lref)\nCue / A')
   if i==0:ax.set_title('Fixed peak A' if j==0 else 'Equal spatial integral at each length',fontsize=12)
   if i==2:ax.set_xlabel('Distance inward from the stimulated end / Lref')
   if i==1:ax.text(.5,.48,'Cues coincide at the reference length',ha='center',transform=ax.transAxes,fontsize=9,color='#444444')
 axes[0,0].legend(loc='center right',fontsize=9)
 fig.suptitle('What is held fixed in the patch-width experiment?',fontsize=15,y=.985)
 fig.text(.5,.945,'Actual wave-pinning pilot examples: Lref = 21, fractional width = 20% of L. Shading shows area; area labels are divided by A*Lref.',ha='center',fontsize=10)
 fig.text(.5,.027,'Left: both peaks stay A; areas differ. Right: the physical-patch peak changes to match the fractional-patch area.\nEqual area means equal between the two cues at a given length; that common area increases with cell length.\nThis is the imposed input integral, not necessarily equal biochemical activation.',ha='center',fontsize=10)
 fig.subplots_adjust(left=.09,right=.98,top=.90,bottom=.14,wspace=.15,hspace=.22)
 save(fig,'cue_width_design')
 p.write(OUT/'cue_width_design_checks.json',dict(generated_at=stamp,passed=True,reference_length=reference,profile_checks=checks))


def mean_cue_design(cases,stamp):
 fig,axes=plt.subplots(2,3,figsize=(15,9),sharey=True)
 model='wave_pinning';checks=[]
 for column,L in enumerate(p.LENGTHS[model]):
  selected=[c for c in cases if c['model']==model and c['length']==L and c['family']=='mean_slope']
  bound=selected[0]['length_bound'];x=np.linspace(0,L,2001)
  half=L/(2*bound)
  top,bottom=axes[:,column]
  top.axvspan(-.55,-half,color='#f0f0f0',zorder=0)
  top.axvspan(half,.55,color='#f0f0f0',zorder=0)
  top.axvline(-half,color='#999999',lw=.8);top.axvline(half,color='#999999',lw=.8)
  for c in selected:
   mean,slope=c['mean_factor'],c['slope_factor']
   cue=p.Cue(c,1)(0,x,L)/c['A']
   np.testing.assert_allclose(np.trapezoid(cue,x)/L,mean,rtol=1e-12,atol=1e-12)
   np.testing.assert_allclose((cue[-1]-cue[0])/L,slope/bound,rtol=1e-12,atol=1e-12)
   np.testing.assert_allclose(cue[len(x)//2],mean,rtol=1e-12,atol=1e-12)
   checks.append(dict(length=L,mean_over_A=mean,slope_over_g0=slope,left_input_over_A=float(cue[0]),right_input_over_A=float(cue[-1]),end_to_end_difference_over_A=float(cue[-1]-cue[0])))
   for ax,position in [(top,(x-L/2)/bound),(bottom,x/L)]:
    ax.plot(position,cue,color=COLORS[mean],ls=STYLES[slope],lw=1.8)
  for mean in MEANS:
   top.plot(0,mean,'o',color=COLORS[mean],ms=4)
   bottom.plot(.5,mean,'o',color=COLORS[mean],ms=4)
  top.axvline(0,color='#777777',lw=.7,ls=':',zorder=0)
  bottom.axvline(.5,color='#777777',lw=.7,ls=':',zorder=0)
  top.set_xlim(-.55,.55);bottom.set_xlim(0,1)
  top.set_title(f'L = {L:g}  ({L/bound:.3g}*length bound)',fontsize=12)
  top.set_xlabel('Cell-centred physical position / length bound')
  bottom.set_xlabel('Fractional position x/L')
  for ax in [top,bottom]:
   ax.set_ylim(0,2.65);ax.grid(axis='y',color='#eeeeee',lw=.6,zorder=0)
 axes[0,0].set_ylabel('Cue / A\nPhysical-coordinate view',fontsize=11)
 axes[1,0].set_ylabel('Cue / A\nSame cues in fractional coordinates',fontsize=11)
 handles=[Line2D([0],[0],color=COLORS[m],lw=2,label=f'Mean = {m:g} A') for m in MEANS]
 handles += [Line2D([0],[0],color='#333333',ls=STYLES[s],lw=2,label=f'Slope = {s:g} g0') for s in SLOPES]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.925),ncol=6,frameon=False,fontsize=10)
 fig.suptitle('Cue design: all nine mean / physical-slope combinations',fontsize=15,y=.985)
 fig.text(.5,.945,'Actual wave-pinning pilot: length bound = 39; g0 = A/39. Colours and line styles match the length-response curves.',ha='center',fontsize=10)
 fig.text(.5,.045,'Within each colour / line-style pair, mean and physical slope stay fixed as length changes.\nLonger cells span a larger end-to-end input difference. The lower row redraws the same cues on x/L.\nDots mark the common mean at the cell centre; grey regions lie outside the cell.',ha='center',fontsize=10)
 fig.subplots_adjust(left=.075,right=.985,top=.84,bottom=.17,wspace=.20,hspace=.30)
 save(fig,'cue_mean_slope_design')
 p.write(OUT/'cue_mean_slope_design_checks.json',dict(generated_at=stamp,passed=True,length_bound=bound,combinations_per_length=9,profile_checks=checks))


def mean_lines(cases,stamp):
 fig=plt.figure(figsize=(15,18))
 outer=fig.add_gridspec(4,2,left=.09,right=.98,bottom=.08,top=.91,hspace=.42,wspace=.32)
 source_rows=[]
 for model_index,(model,label) in enumerate(zip(MODELS,LABELS)):
  spec=outer[model_index//2,model_index%2]
  selected=[c for c in cases if c['model']==model and c['family']=='mean_slope']
  conditions={(mean,slope):sorted([c for c in selected if c['mean_factor']==mean and c['slope_factor']==slope],key=lambda c:c['length']) for mean in MEANS for slope in SLOPES}
  censored_conditions=[key for key,items in conditions.items() if any(c['result']['status']=='censored' for c in items)]
  height=.22*len(censored_conditions)+.3 if censored_conditions else .28
  sub=spec.subgridspec(2,1,height_ratios=[2.8,height],hspace=.16)
  ax=fig.add_subplot(sub[0]);strip=fig.add_subplot(sub[1],sharex=ax)
  lengths=p.LENGTHS[model];dx=lengths[-1]-lengths[0];bound=selected[0]['length_bound']
  successes=[]
  for (mean,slope),items in conditions.items():
   x=np.asarray([c['length'] for c in items]);y=np.asarray([c['result']['endpoint'] if c['result']['status']=='success' else np.nan for c in items])
   ax.plot(x,y,color=COLORS[mean],ls=STYLES[slope],lw=1.5,zorder=2)
   for c in items:
    q=c['result']
    source_rows.append(dict(case_index=c['case_index'],model=model,length=c['length'],mean_factor=mean,slope_factor=slope,physical_slope=slope*c['A']/bound,mean_input=mean*c['A'],status=q['status'],endpoint=q['endpoint'],maximum_tested_duration=q['maximum_tested_duration'],nonmonotone=q['selected_morphology']=='nonmonotone',field_tail_pass=q['selected_tail_field_pass_0_1_percent']))
    if q['status']!='success':continue
    successes.append(q['endpoint'])
    marker='s' if q['selected_morphology']=='nonmonotone' else 'o'
    if q['lower_grid_limited']:marker='v'
    face=COLORS[mean] if q['selected_tail_field_pass_0_1_percent'] else 'white'
    ax.plot(c['length'],q['endpoint'],marker,color=COLORS[mean],mfc=face,mec=COLORS[mean],ms=5,mew=1.15,zorder=3)
  ax.set_xlim(lengths[0]-.08*dx,lengths[-1]+.08*dx)
  ax.set_xticks(lengths);ax.set_yscale('log')
  if successes:
   ax.set_ylim(min(successes)/1.22,max(successes)*1.3)
  else:
   maximum=max(c['result']['maximum_tested_duration'] for c in selected)
   ax.set_ylim(maximum/1000,maximum)
   ax.text(.5,.5,f'No tested pulse passed\nMaximum tested duration: {maximum:g}',transform=ax.transAxes,ha='center',va='center',fontsize=11,color='#555555')
  ax.set_title(f'{label}  |  g0 = A / {bound:g}',fontsize=12)
  ax.set_ylabel('Passing pulse duration\n(native time units)')
  ax.tick_params(axis='x',labelbottom=False)
  ax.grid(axis='y',which='major',color='#eeeeee',lw=.7,zorder=0)
  strip.spines[['top','right','left','bottom']].set_visible(False)
  strip.set_xticks(lengths);strip.set_xlabel('Resting L0' if model=='spring' else 'Nominal cell length',labelpad=2)
  strip.tick_params(axis='y',length=0,labelsize=7)
  if censored_conditions:
   strip.set_ylim(len(censored_conditions)-.5,-.5)
   strip.set_yticks(range(len(censored_conditions)))
   strip.set_yticklabels([f'{mean:g} / {slope:g}' for mean,slope in censored_conditions])
   strip.set_ylabel('No passing pulse\nmean/A / slope/g0',fontsize=7,labelpad=5)
   for row,key in enumerate(censored_conditions):
    mean,slope=key
    strip.axhline(row,color='#f0f0f0',lw=.6,zorder=0)
    for c in conditions[key]:
     if c['result']['status']=='censored':strip.plot(c['length'],row,'x',color=COLORS[mean],ms=5,mew=1.2)
  else:
   strip.set_yticks([]);strip.set_ylim(0,1)
   strip.text(.5,.65,'All nine conditions passed at all three lengths',ha='center',va='center',transform=strip.transAxes,fontsize=8,color='#555555')
 legendax=fig.add_subplot(outer[3,1]);legendax.axis('off')
 mean_handles=[Line2D([0],[0],color=COLORS[m],lw=2,label=f'Mean input = {m:g} A') for m in MEANS]
 slope_handles=[Line2D([0],[0],color='#333333',ls=STYLES[s],lw=2,label=f'Physical slope = {s:g} g0') for s in SLOPES]
 marker_handles=[Line2D([0],[0],color='#555555',marker='o',ls='none',label='Approximately monotone end profile'),Line2D([0],[0],color='#555555',marker='s',ls='none',label='Nonmonotone retained profile'),Line2D([0],[0],color='#555555',marker='o',mfc='white',ls='none',label='Fields still vary >0.1% fixed scale in tail'),Line2D([0],[0],color='#555555',marker='x',ls='none',label='No passing pulse (shown below each panel)')]
 legendax.legend(handles=mean_handles+slope_handles+marker_handles,loc='upper left',frameon=False,fontsize=10,labelspacing=1.0)
 fig.suptitle('Finite countercue: length response across all nine mean / slope combinations',fontsize=15,y=.975)
 fig.text(.5,.943,'Colour = mean input. Line style = physical slope. Each point is the shortest passing tested pulse; recovery time is excluded.',ha='center',fontsize=10)
 fig.text(.5,.023,'The polarity gate can pass with additional active regions or continued spatial relaxation; marker shape and fill expose those outcomes.\nNative time units and y-axis ranges differ between models. Generated '+stamp,ha='center',fontsize=9)
 save(fig,'mean_slope_length_response')
 assert len(source_rows)==189
 r.csvwrite(OUT/'mean_slope_length_response_data.csv',source_rows)
 return source_rows


def dashboard(stamp):
 status=p.read(OUT/'status.json');data_stamp=status['generated_at']
 text=f'''<!doctype html><html><head><meta charset="utf-8"><title>Finite-countercue pilot plots</title><style>body{{font-family:sans-serif;max-width:1500px;margin:30px auto;line-height:1.5}}img{{width:100%}}summary{{cursor:pointer}}</style></head><body><h1>Finite-countercue width and mean/slope pilot</h1><p>Plots generated {stamp}; producer finished. Simulation data generated {data_stamp}; simulation producer finished.</p><p>Exact original cue-free release states. Thresholds retain the original orientation/contrast gate. Field-tail variation and retained morphology are reported separately. Three lengths per model; no length regressions.</p><p><a href="endpoints.csv">All endpoints</a> · <a href="validation.json">Simulation validation</a> · <a href="../manifest.json">Design</a></p>'''
 for name,title,caption in [
  ('cue_width_design','What fixed peak and equal spatial integral mean','Shading shows the input area. Integral matching is between physical and fractional patches at the same cell length; it does not hold the common area constant across lengths.'),
  ('width_thresholds','Length response: physical versus fractional width','Peak-matched and integral-matched patch controls. Hollow upward triangles mark no passing pulse; squares mark nonmonotone retained fields.'),
  ('../../full_reorientation_width_20261002/results/width_full_reorientation','Width companion: full polarity reorientation','Recovered polarity must remain within 1% of the reversed initial magnitude throughout the final 500 recovery units, with the original gates retained. Original width thresholds above are unchanged. <a href="../../full_reorientation_width_20261002/results/index.html">Companion endpoints, definition and validation</a>.'),
  ('../../width_full_lengths_20261002/results/length_response_two_definitions','All lengths: three cues and both reorientation definitions','The full canonical length grid, with 20% fractional patches and the two fixed physical-patch controls. Fixed widths are restricted to patches fitting within half the cell. Both definitions share the same tested durations. <a href="../../width_full_lengths_20261002/results/index.html">Full rollout report, endpoints and validation</a>.'),
  ('cue_mean_slope_design','Cue design: all nine mean and physical slope combinations','The upper row uses a common physical-coordinate scale; the lower row redraws the same cues on x/L. Mean and physical slope stay fixed across lengths for each colour and line style, while the end-to-end difference grows.'),
  ('mean_slope_length_response','Length response: all nine mean and slope combinations','Colour identifies mean input; line style identifies physical slope. Missing durations break the lines; censored conditions are listed below each panel. Hollow points flag field-tail variation.'),
  ('spatial_outcomes','Spatial outcomes at the longest pilot length','Initial, cue-withdrawal and final recovery profiles for three representative cues.'),
  ('recovery_behaviour','Countercue and recovery across lengths','Seven models × three width conditions, with all three lengths overlaid. Polarity is divided by the magnitude of initial settled polarity; time starts at countercue onset and vertical ticks mark withdrawal. Solid traces use the shortest passing pulse; dashed traces use the maximum tested pulse where none passed.')]:
  text+=f'<h2>{title}</h2><p>{caption}</p><a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a><img src="{name}.png">'
 text+='<details><summary>Earlier heatmap view of the same mean/slope data</summary><p>Numbers are shortest passing tested pulse durations. The colour scale is shared across lengths within a model.</p><a href="mean_slope_thresholds.pdf">PDF</a><img src="mean_slope_thresholds.png"></details></body></html>'
 (OUT/'index.html').write_text(text)


def main():
 assert p.read(OUT/'status.json')['complete']
 original_records_sha=r.sha(OUT/'records.json.gz')
 cases=json.load(gzip.open(OUT/'records.json.gz','rt'))['cases'];stamp=r.now()
 cue_design(cases,stamp);mean_cue_design(cases,stamp);data=mean_lines(cases,stamp);dashboard(stamp)
 assert r.sha(OUT/'records.json.gz')==original_records_sha
 p.write(OUT/'plot_revision_provenance.json',dict(generated_at=stamp,source_sha256=r.sha(Path(__file__).resolve()),records_sha256=original_records_sha,simulation_data_unchanged=True,mean_slope_points=len(data),censored_points=sum(x['status']=='censored' for x in data),new_figures=['cue_width_design','cue_mean_slope_design','mean_slope_length_response'],previous_heatmaps_retained=True))
 p.write(OUT/'artifact_hashes.json',{f.name:r.sha(f) for f in sorted(OUT.iterdir()) if f.suffix in ['.png','.pdf','.svg','.csv']})
 print(json.dumps(p.read(OUT/'plot_revision_provenance.json')),flush=True)

if __name__=='__main__':main()
