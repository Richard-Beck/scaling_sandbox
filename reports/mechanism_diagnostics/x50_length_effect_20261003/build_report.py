"""Vector PDF report of previously audited historical reversal diagnostics."""
from pathlib import Path
import csv,json,textwrap,hashlib
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path('/home/4473331/projects/scaling_sandbox')
OUT=Path(__file__).resolve().parent
DATA=ROOT/'tmp/x50_mechanisms_20261003'
FREEZE=ROOT/'tmp/spring_tension_freeze_20261003'
VIEWER=ROOT/'polarity1d/four_assay/results/reversal_field_report'
extra_sources=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'mathtext.fontset':'stix','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
BLUE='#2474b5';ORANGE='#dd7c25';INK='#193449';MUTED='#516575';GREEN='#32815a'
STAMP=datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
def csvrows(path):return list(csv.DictReader(path.open()))
def numeric(rows,key):return np.array([float(r[key]) for r in rows])
def jload(path):return json.loads(path.read_text())
def body(fig,y,text,size=10,width=94,color=INK):
 lines=[]
 for paragraph in text.split('\n'):
  lines.extend(textwrap.wrap(paragraph,width=width,break_long_words=False,break_on_hyphens=False) or [''])
 fig.text(.09,y,'\n'.join(lines),fontsize=size,color=color,va='top',linespacing=1.5)
 return y-len(lines)*.019*(size/10)-.018
def eq(fig,y,text,size=15):fig.text(.10,y,'$'+text+'$',fontsize=size,color=INK,va='top')
def page(title,subtitle=''):
 fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
 fig.add_artist(plt.Line2D([.09,.91],[.942,.942],transform=fig.transFigure,color=BLUE,lw=2))
 title_artist=fig.text(.09,.921,title,fontsize=19,fontweight='bold',color=INK,va='top')
 fig.canvas.draw()
 while title_artist.get_window_extent(fig.canvas.get_renderer()).x1>.92*fig.bbox.width:
  title_artist.set_fontsize(title_artist.get_fontsize()-.5)
  fig.canvas.draw()
 if subtitle:fig.text(.09,.883,subtitle,fontsize=9,color=MUTED,va='top')
 return fig
def finish(pdf,fig,n):
 fig.text(.09,.033,'Reversal mechanisms | historical diagnostic analysis | '+STAMP,fontsize=7,color=MUTED)
 fig.text(.91,.033,str(n),fontsize=8,ha='right',color=MUTED)
 pdf.savefig(fig);plt.close(fig)
def style(ax):ax.grid(alpha=.17);ax.tick_params(labelsize=8)
def table(fig,rect,headers,rows,widths=None):
 ax=fig.add_axes(rect);ax.axis('off')
 t=ax.table(cellText=rows,colLabels=headers,loc='center',cellLoc='left',colLoc='left',colWidths=widths)
 t.auto_set_font_size(False);t.set_fontsize(8);t.scale(1,1.8)
 for (r,c),cell in t.get_celld().items():
  cell.set_edgecolor('#dce4e9');cell.set_linewidth(.6)
  if r==0:cell.set_facecolor(INK);cell.set_text_props(color='white',weight='bold')
  elif r%2:cell.set_facecolor('#eff4f7')
 return t
def axes_pair(fig,y=.22,height=.29):return [fig.add_axes([.10,y,.35,height]),fig.add_axes([.57,y,.35,height])]

PDF=OUT/'reversal_length_effect_mechanisms.pdf'
with PdfPages(PDF,metadata={'Title':'Why reversal depends on length: half-mass dynamics and mechanistic evidence','Author':'Scaling sandbox diagnostic analysis','Subject':'Wave-pinning, Holmes, Spring, and deBelly; audited historical patch trajectories'}) as pdf:
 # 1: interpretation and evidence hierarchy.
 fig=page('Why reversal depends on length','Half-mass dynamics, competing reactions, MCA transport, and a tension-feedback intervention')
 y=body(fig,.837,'Main finding: longer cells reverse more slowly because the internal redistribution needed to commit to the new pole becomes harder. For Wave-pinning, Holmes and Spring, the measured slowdown appears as nearly balanced opposing reactions. For deBelly, cortical transport must rearrange MCA over a greater distance.')
 table(fig,[.09,.57,.82,.18],['Model / marker','Measured rate control','Interpretation'],[
  ['Wave-pinning / active','Cue versus intrinsic conversion','Consistent with mesa competition'],
  ['Holmes / active Rho','Nearly balanced Rho production','Rac competition seen through Rho'],
  ['Spring / active','Opposing conversion; bounded mechanics','Tension helps, but does not remove delay'],
  ['deBelly / bound MCA','Cortical advection and exchange','More transport distance; limited speed']], [.25,.36,.39])
 y=body(fig,.54,'Saturation is a supported interpretation, not an independently isolated cause in the chemical models. The rate identities locate the actual slow processes. The Spring tension-clamp experiment provides a separate causal result: dynamic tension lowers the required pulse duration, particularly at the longer length.')
 y=body(fig,y,'Scope. This report assembles the historical localized-patch experiments inspected in the field explorer, plus their isolated Spring ablation. Their starting states came from the earlier linear-conditioned/released campaign. They are not the corrected finite-patch preparations or length curves of the newer active manuscript campaign.',size=9,color=MUTED)
 body(fig,y,'Reading order: page 2 derives the shared identity; page 3 defines experiments and response criteria; pages 4-7 quantify the four models; page 8 gives the causal tension test; page 9 checks midpoint commitment; page 10 records validation; pages 11-13 give additional equations; page 14 shows the starting profiles.',size=9)
 finish(pdf,fig,1)
 # 2: exact implicit derivation.
 fig=page('The half-mass position and its velocity','Exact identity; the source terms need not conserve the selected field')
 body(fig,.835,'Let c(x,t) be a nonnegative field on a static interval [0,L], with flux J and source R. Define m(t) as the coordinate with half of the field integral on either side. This is a mass median, not a peak position, half-height, or front.')
 eq(fig,.755,r'\partial_t c=-\partial_x J+R,\qquad M(t)=\int_0^L c\,dx')
 eq(fig,.692,r'\int_0^{m(t)}c\,dx=\frac{M(t)}{2}')
 eq(fig,.628,r'c(m)\dot m+\int_0^m\partial_t c\,dx=\frac{1}{2}\int_0^L\partial_t c\,dx')
 eq(fig,.561,r'\mathcal{B}_m[g]\equiv\frac{1}{2}\left(\int_m^L g\,dx-\int_0^m g\,dx\right)')
 eq(fig,.493,r'\dot m=\frac{J(m)-\frac{1}{2}[J(0)+J(L)]+\mathcal{B}_m[R]}{c(m)}',size=16)
 body(fig,.417,'With diffusion only and zero boundary flux, J=-D*c_x and the boundary term vanishes. The source imbalance drives motion toward the side gaining relatively more field mass. The changing total M is essential; dropping it gives an incorrect velocity.')
 eq(fig,.331,r'R=-k(t)c\quad\Longrightarrow\quad\mathcal{B}_m[R]=0')
 body(fig,.273,'A spatially uniform proportional loss cannot move the median instantaneously. It can still change future profiles and therefore future movement. A spatially uniform additive source does not cancel: its contribution is (L/2-m)*R. Small c(m) can make a median move quickly without a physical front travelling at that speed.')
 body(fig,.144,'For normalized position q50=m/L on a static domain, q50_dot=m_dot/L. Spring requires a moving-domain treatment; deBelly retains advective boundary flux and a numerical mass correction. Those cases are derived separately.',size=9,color=MUTED)
 finish(pdf,fig,2)
 # 3: design and definitions.
 fig=page('Which experiments are being explained?','Localized raised-cosine countercues; fixed peak, spatial support, and input-area controls')
 body(fig,.837,'Each experiment begins from its recorded left-polarized full checkpoint, applies a finite right patch for duration T, and then recovers without added cue for 10,000 native units. The reported endpoint is the shortest tested passing pulse; recovery time is excluded.')
 eq(fig,.747,r'S(x,t)=\frac{A}{2}\left[1+\cos\!\left(\pi\frac{L-x}{w}\right)\right],\quad 0\leq L-x\leq w,\quad t<T',size=13)
 body(fig,.693,'Outside the support, or after cutoff, the extra cue is zero. Native basal reactions continue. Spring evaluates the patch on its instantaneous physical domain; deBelly uses nominal solver coordinates.',size=9)
 table(fig,[.09,.46,.82,.17],['Condition','Width','Peak'],[
  ['Fixed width / peak','0.2 times reference domain','A'],
  ['Fractional width / peak','0.2 times current domain','A'],
  ['Fixed width / equal integral','Same fixed support','A times current/reference domain']], [.34,.36,.30])
 body(fig,.433,'Historical reference lengths: WP 21, Holmes 30, deBelly 27, Spring resting length 1.6 (its released physical domain sets the fixed support). Native peak amplitudes are 0.1, 1, 3.6, and 4 respectively. These differ from the current manuscript width conventions.',size=9)
 table(fig,[.09,.235,.82,.13],['Model','Paired native lengths','Tracked field'],[['WP','9 and 25','Active protein'],['Holmes','14 and 24','Active Rho'],['Spring','0.8 and 2.4','Active protein'],['deBelly','19 and 35','Bound MCA']], [.27,.35,.38])
 body(fig,.192,'Orientation success requires stable right orientation and contrast in the final 500 recovery units. Scalar strength replacement additionally requires the reversed orientation magnitude to match its starting magnitude within 1%. Complete-field tail stability is a separate flag; scalar success does not establish full-profile convergence.',size=9)
 finish(pdf,fig,3)
 # 4: wave pinning.
 fig=page('Wave-pinning: opposing conversion slows redistribution','Active u; inactive v; no-flux static domain')
 eq(fig,.834,r'\partial_t u=0.1u_{xx}+\left[0.067+\frac{u^2}{1+u^2}+S\right]v-u')
 eq(fig,.766,r'\dot m=\frac{-0.1u_x(m)+\mathcal{B}_m[0.067v]+\mathcal{B}_m[\frac{u^2}{1+u^2}v]+\mathcal{B}_m[Sv]}{u(m)}',size=13)
 body(fig,.700,'Uniform active decay cancels directly. Cue conversion moves active mass rightward; the established pole and intrinsic conversion oppose it. At the same q50=0.35, normalized speed falls from 0.03345 to 0.00607 for fixed width, and from 0.01710 to 0.00833 for fractional width.',size=10)
 body(fig,.603,'With fractional width, the direct cue contribution is larger in the longer cell (0.04850 versus 0.02541), yet opposing intrinsic reactions leave slower net movement. Direct active diffusion also contributes less. This is evidence for redistribution/competition rather than a smaller applied input fraction.',size=9)
 axs=axes_pair(fig)
 for L,color in [(9,BLUE),(25,ORANGE)]:
  rr=csvrows(DATA/f'wave_pinning/wave_pinning_physical_peak_matched_L{L}_decomposition.csv')
  axs[0].plot(numeric(rr,'time'),numeric(rr,'q50'),color=color,label=f'L={L}')
 axs[0].axhline(.5,color=MUTED,ls=':',lw=1);axs[0].set(xlabel='Cue-on native time',ylabel='Active half-mass position / L',title='Fixed-width trajectories');axs[0].legend(fontsize=8)
 stages=jload(DATA/'wave_pinning/common_stage.json')
 names=['cue','basal','feedback','diffusion','total'];positions=np.arange(5)
 for j,(L,color) in enumerate([(9,BLUE),(25,ORANGE)]):
  z=next(z for z in stages if z['length']==L and z['q50']==.35 and 'fractional_peak' in z['file'])
  axs[1].bar(positions+(j-.5)*.35,[z['q_velocity'][n] for n in names],width=.35,color=color,label=f'L={L}')
 axs[1].axhline(0,color=MUTED,lw=.7);axs[1].set_xticks(positions,names,rotation=40,ha='right');axs[1].set(ylabel='Contribution to dq50/dt',title='Fractional cue at q50=0.35');axs[1].legend(fontsize=8)
 for ax in axs:style(ax)
 body(fig,.15,'The plateaus broaden with little change in height, consistent with saturation-dependent mesa competition. The exact decomposition establishes opposing balances, not a causal isolation of saturation. Substrate diffusion remains embedded in v even where it is not an explicit median-flux term.',size=9)
 finish(pdf,fig,4)
 # 5: Holmes.
 fig=page('Holmes: a late bottleneck in Rho production','The countercue affects the Rho median indirectly through Rac')
 eq(fig,.837,r'p_h=\frac{3.3h_i}{3.1(1+r^3)},\quad \partial_t h=0.1h_{xx}+p_h-h')
 eq(fig,.771,r'\dot m=\frac{-0.1h_x(m)+\mathcal{B}_m[p_h]}{h(m)}')
 body(fig,.708,'Here r is active Rac, h active Rho, and h_i inactive Rho. The cue adds to Rac activation and bypasses inhibition of its basal activation. Rac in turn suppresses Rho production. Uniform Rho decay cancels from the median velocity.',size=10)
 body(fig,.613,'Initial movement is similar at both lengths, but the long-cell production balance stalls late. At fixed-width cutoff L=24, left/right production integrals are 15.42076 and 15.41437: a difference of only 0.0415%. Direct active-Rho diffusion contributes below 1% of integrated movement.',size=9)
 axs=axes_pair(fig)
 rr=csvrows(DATA/'holmes/cue_on_rates.csv')
 for L,color in [(14,BLUE),(24,ORANGE)]:
  chosen=[r for r in rr if r['condition']=='physical_peak_matched' and float(r['length'])==L]
  t=numeric(chosen,'time');axs[0].plot(t,numeric(chosen,'x50_over_L'),color=color,label=f'L={L}')
  speed=-numeric(chosen,'production_rate');valid=(t>0)&(speed>0)
  axs[1].plot(t[valid],speed[valid],color=color,label=f'L={L}')
 axs[0].axhline(.5,color=MUTED,ls=':',lw=1);axs[0].set(xlabel='Cue-on native time',ylabel='Rho half-mass position / L',title='Fixed-width trajectories')
 axs[1].set(xscale='log',yscale='log',xlabel='Cue-on native time',ylabel='Leftward production-driven speed',title='Late near-cancellation')
 for ax in axs:ax.legend(fontsize=8);style(ax)
 body(fig,.15,'Inactive-Rho variation is only about 0.2%, yet removes 49-74% of the small driving imbalance relative to a uniform-pool approximation. The evidence is consistent with broad Rac-pole competition; it does not show that the Rho-inhibition Hill term is itself fully saturated.',size=9)
 finish(pdf,fig,5)
 # 6: Spring.
 fig=page('Spring: reaction competition with bounded mechanics','Material coordinate q; m=q50; physical length L changes with time')
 eq(fig,.839,r'u_t=\frac{0.01}{L^2}u_{qq}+(4+S+H(u))v-\delta u-\frac{\dot L}{L}u',size=14)
 eq(fig,.783,r'H(u)=\frac{5u^6}{1+u^6},\qquad \delta=3+80(L-L_0)')
 eq(fig,.726,r'\dot m=\frac{-\frac{0.01}{L^2}u_q(m)+\mathcal{B}_m[(4+S+H(u))v]}{u(m)}',size=15)
 body(fig,.663,'Here B integrates over material q in [0,1]. Uniform deactivation and dilution cancel directly, but change future profiles. Absolute physical position is X50=a+Lm, with velocity a_dot+L_dot*m+L*m_dot. The plot uses material position, avoiding rigid translation.',size=9)
 eq(fig,.578,r'P(u)=\frac{0.001}{1+e^{-20(u-1)}},\quad\dot L=P(u_l)+P(u_r)-0.02(L-L_0)',size=13)
 axs=axes_pair(fig,y=.205,height=.29)
 sr=csvrows(DATA/'spring/summary.csv')
 names=['cue','basal','autocatalysis','diffusion'];xx=np.arange(4)
 for j,(L,color) in enumerate([(.8,BLUE),(2.4,ORANGE)]):
  z=next(r for r in sr if r['condition']=='fractional_peak_matched' and float(r['L0'])==L)
  axs[0].bar(xx+(j-.5)*.35,[float(z[n+'_integrated_shift']) for n in names],width=.35,color=color,label=f'L0={L}')
 axs[0].set_xticks(xx,names,rotation=35,ha='right');axs[0].axhline(0,color=MUTED,lw=.7);axs[0].set(ylabel='Integrated contribution to q50 shift',title='Fractional cue: opposing conversion');axs[0].legend(fontsize=8)
 us=np.linspace(.5,2.4,200);P=.001/(1+np.exp(-20*(us-1)))
 axs[1].plot(us,P/.001,color=GREEN)
 for L,color in [(.8,BLUE),(2.4,ORANGE)]:
  z=next(r for r in sr if r['condition']=='fractional_peak_matched' and float(r['L0'])==L)
  axs[1].plot(float(z['right_u_cutoff']),float(z['right_P_cutoff'])/.001,'o',color=color,label=f'L0={L} cutoff')
 axs[1].set(xlabel='Boundary active concentration',ylabel='Protrusion / maximum',title='Endpoint sensing saturates');axs[1].legend(fontsize=7)
 for ax in axs:style(ax)
 body(fig,.135,'Intrinsic conversion opposes 7.6% of direct cue displacement in the short fractional-cue example, versus 55.9% in the long one. At long-cell cutoff, challenged-end protrusion exceeds 99.999997% of its maximum. Peak width adds no protrusion; extra height adds almost none.',size=9)
 finish(pdf,fig,6)
 # 7: MCA.
 fig=page('deBelly: MCA rearrangement is transport dominated','Bound-MCA median in the conserved nominal/reference-coordinate measure')
 eq(fig,.838,r'b_t=-\partial_x F+R+C_b,\quad F=vb-0.001b_x,\quad R=5n-1.4b',size=14)
 eq(fig,.771,r'\dot m=\frac{F(m)-\frac{1}{2}[F(0)+F(L)]+\mathcal{B}_m[R]+(L/2-m)C_b}{b(m)}',size=13)
 body(fig,.707,'b is bound MCA, n unbound MCA, and v cortical velocity. Attachment/detachment rates are constant. Rac acts through protrusion and mechanics, generating flow; MCA redistribution then alters Rac activation. Boundary advection is not zero. C_b is the retained numerical uniform mass-restoration rate.',size=9)
 body(fig,.606,'Fixed width: median travel increases 1.285 to 3.868 nominal microns, while average speed decreases 0.00691 to 0.00566 microns/s. The starting normalized median is farther from its release value in the larger cell. Advection dominates; binding exchange opposes part of the displacement.',size=9)
 axs=axes_pair(fig)
 for L,color in [(19,BLUE),(35,ORANGE)]:
  rr=csvrows(DATA/f'debelly/debelly_physical_peak_matched_L{L}_rates.csv');t=numeric(rr,'t')
  axs[0].plot(t,numeric(rr,'q50'),color=color,label=f'L={L}')
 ds=jload(DATA/'debelly/summary.json');names=['advection','exchange','uniform_closure','diffusion'];xx=np.arange(4)
 for j,(L,color) in enumerate([(19,BLUE),(35,ORANGE)]):
  z=next(z for z in ds if z['condition']=='physical_peak_matched' and z['length']==L)
  axs[1].bar(xx+(j-.5)*.35,[z['median_integrated_displacements'][n] for n in names],width=.35,color=color,label=f'L={L}')
 axs[0].axhline(.5,color=MUTED,ls=':',lw=1);axs[0].set(xlabel='Cue-on native time',ylabel='Bound-MCA half-mass position / L',title='Fixed-width trajectories');axs[0].legend(fontsize=8)
 axs[1].axhline(0,color=MUTED,lw=.7);axs[1].set_xticks(xx,['advection','exchange','closure','diffusion'],rotation=35,ha='right');axs[1].set(ylabel='Median displacement (nominal microns)',title='Integrated motion contributions');axs[1].legend(fontsize=8)
 for ax in axs:style(ax)
 body(fig,.15,'The uniform numerical correction contributes materially and must remain visible; it is not an independent biological process. With reference density, physical median X50=m+d(m), giving X50_dot=d_t(m)+(1+d_x(m))*m_dot. These examples fail full-field tail stability despite scalar success.',size=9)
 finish(pdf,fig,7)
 # 8: causal freeze test.
 fig=page('Causal test: freeze the biochemical tension response','Same checkpoint and cue; boundary motion, diffusion and dilution remain dynamic')
 eq(fig,.837,r'\delta_{\rm normal}(t)=3+80[L(t)-L_0]')
 eq(fig,.778,r'\delta_{\rm frozen}=3+80[L(0)-L_0]')
 body(fig,.716,'Freezing at the initial value preserves the starting reaction rate. At the long length, each original shortest passing pulse fails under frozen feedback; longer pulses rescue reversal. Dynamic tension therefore helps competition resolve, but is not indispensable.',size=10)
 table(fig,[.09,.50,.82,.15],['Cue','L0=0.8: normal to frozen','L0=2.4: normal to frozen'],[['Fixed width / peak','0.135 to 0.136 (+0.8%)','2.240 to 3.063 (+37%)'],['Fractional width / peak','0.237 to 0.241 (+1.6%)','1.110 to 1.292 (+16%)'],['Fixed width / equal integral','0.263 to 0.263 (+0%)','1.614 to 2.005 (+24%)']],[.34,.33,.33])
 axs=axes_pair(fig,y=.19,height=.25)
 for k,slot in enumerate([1,3]):
  z=jload(FREEZE/f'result_{slot:02d}.json');D=z['normal']['duration']
  for freeze,color in [(False,BLUE),(True,ORANGE)]:
   name=f'{"frozen" if freeze else "normal"}_D{D:.12g}'
   with np.load(FREEZE/f'case_{slot:02d}'/(name+'.npz')) as tr:axs[k].plot(tr['time'],tr['orientation'],color=color,label='Frozen' if freeze else 'Normal')
  axs[k].axvline(D,color=MUTED,ls=':',label='Cue off');axs[k].axhline(0,color=MUTED,lw=.6);axs[k].set(xscale='symlog',xlabel='Native time',ylabel='Orientation',title=['Fixed width, L0=2.4','Fractional width, L0=2.4'][k]);axs[k].legend(fontsize=7);style(axs[k])
 body(fig,.126,'95 simulations completed. Normal fields reproduce original trajectories exactly; frozen no-cue controls retain their original polarity. Selected frozen thresholds pass both scalar strength and field-tail stability; brackets are within 1%. Short-cell differences are near threshold resolution.',size=9)
 finish(pdf,fig,8)
 # 9: release criterion.
 fig=page('Does midpoint crossing determine cue release?','It tracks redistribution, but is not a universal commitment threshold')
 body(fig,.835,'The median marks equal integrated field mass, not equal peak amplitudes or a basin boundary. Successful and immediately shorter failing pulses can have very similar medians. Release values also differ strongly between lengths, making a universal q50=0.5 rule untenable.')
 table(fig,[.09,.59,.82,.16],['Model','Short-cell successful q50','Long-cell successful q50'],[['Wave-pinning','0.483-0.503','0.748-0.756'],['Holmes active Rho','0.468-0.472','0.411-0.413'],['Spring active','0.487-0.493','0.528-0.537'],['deBelly bound MCA','0.504-0.511','0.511-0.516']],[.35,.325,.325])
 body(fig,.553,'For Holmes, adjacent failing pulses have already crossed the center. For Spring, successful short-cell pulses release before the active median crosses the center. For deBelly, successful pulses release while the MCA median remains on its original side. The complete profile and post-cue dynamics still matter.',size=10)
 body(fig,.436,'A useful separation is required displacement versus average movement rate. For Spring fractional stimulation, the long/short duration ratio is 4.68: 2.33 times greater normalized median displacement multiplied by a 2.01-fold reduction in average normalized speed. The initial rates are nearly equal; the later reaction balance differs.',size=10)
 eq(fig,.313,r'T=\frac{|\Delta q_{50}|}{|\overline{\dot q}_{50}|},\qquad\frac{T_l}{T_s}=\frac{|\Delta q_l|}{|\Delta q_s|}\frac{|\overline{\dot q}_s|}{|\overline{\dot q}_l|}',size=15)
 body(fig,.242,'This is an exact bookkeeping identity for measured net displacement and signed mean velocity, not a predictive scaling law. It separates how far the marker moves from how quickly it moves. Predicting commitment still requires the full starting profile and a release condition tied to subsequent dynamics.',size=9)
 body(fig,.128,'Interpretation: the derivations explain where the measured length-dependent slowing arises. Causal attribution of chemical saturation would require a targeted perturbation that changes saturation while controlling preparation and the other rate terms.',size=9,color=MUTED)
 finish(pdf,fig,9)
 # 10: equations native details and checks.
 fig=page('Validation, model conventions and provenance','Reproducible derivations and retained numerical evidence')
 table(fig,[.09,.68,.82,.14],['Model','Independent identity check','Trajectory integration'],[['WP','Finite-difference rate RMS: 0.6-1.6%','Displacement error: 0.13-0.37%'],['Holmes','Directional derivative error < 1.95e-8','Displacement error about 0.9-5%'],['Spring','Directional derivative error < 4.2e-10','Displacement error about 0.03-0.99%'],['deBelly','Native update residual < 2.48e-16','Median decomposition error < 0.00143%']],[.16,.43,.41])
 y=body(fig,.644,'The four analyses include both lengths and all three historical cue conditions. The figures emphasize fixed-width or common-stage comparisons. Each length uses its own shortest tested scalar-strength passing pulse; these are not matched-duration controls. The Spring ablation is matched within each case.',size=9)
 y=body(fig,y,'Discretization. Holmes and deBelly use native complete-cell mass; WP and Spring use piecewise-linear profiles. Small differences from viewer quadrature are documented in the analyses. Direct active-field diffusion is distinct from indirect diffusion effects on inactive pools. The fixed 201-point grids imply larger physical spacing at larger lengths; systematic mesh sensitivity is outside this report.',size=9)
 y=body(fig,y,'Field convergence. All selected deBelly trajectories and four of six Holmes trajectories fail the separate native-field tail flag despite scalar success. Interpret their rates as measured transients, not certified attractor dynamics. The source and decay cancellations are exact identities independent of that flag.',size=9)
 y=body(fig,y,'Source boundary. The canonical polarity1d/manuscript workflow now uses corrected finite-patch/basal preparation and has completed its separate fresh reversal campaign. This report does not import these historical curves into that comparison. Goryachev, Otsuji and LEGI have no comparable mechanistic analysis here; failure under a tested stimulus is not a universal model exclusion.',size=9)
 y=body(fig,y,'Primary theoretical reference: Chiou J-G et al. (2018), Principles that govern competition or co-existence in Rho-GTPase driven polarization. PLOS Computational Biology 14:e1006095. doi:10.1371/journal.pcbi.1006095. Its saturation/competition framework motivates the interpretation; the numerical balances in this report are derived from the retained implementations.',size=9)
 y=body(fig,y,'Local evidence: tmp/x50_mechanisms_20261003/{wave_pinning,holmes,spring,debelly}; tmp/spring_tension_freeze_20261003. Model sources: frozen width_full_lengths_20261002 campaign. build_report.py regenerates this PDF with the existing scaling_sandbox environment; provenance.json hashes the inputs and output.',size=8,width=110,color=MUTED)
 finish(pdf,fig,10)
 # 11: native chemical equations and moving-boundary mechanics.
 fig=page('Additional governing equations','Parameter values are those retained in the diagnostic implementations')
 body(fig,.837,'Wave-pinning: the inactive field receives the opposite conversion term and diffuses much faster. Total active plus inactive protein is conserved on the static domain.',size=9)
 eq(fig,.761,r'f_u=\left(0.067+\frac{u^2}{1+u^2}+S\right)v-u,\qquad v_t=10v_{xx}-f_u',size=14)
 body(fig,.687,'Holmes: r and r_i denote active/inactive Rac; h and h_i denote active/inactive Rho. The retained no-PI regime uses Rac basal activation 3.5. Added cue S bypasses the Rho-inhibited basal Rac term.',size=9)
 eq(fig,.607,r'f_r=\left[\frac{3.5}{2(1+(h/1.25)^3)}+S\right]\frac{r_i}{7.5}-r',size=15)
 eq(fig,.540,r'f_h=\frac{6.6}{2(1+r^3)}\frac{h_i}{3.1}-h',size=15)
 eq(fig,.474,r'r_t=0.1r_{xx}+f_r,\quad (r_i)_t=50(r_i)_{xx}-f_r',size=14)
 eq(fig,.418,r'h_t=0.1h_{xx}+f_h,\quad (h_i)_t=50(h_i)_{xx}-f_h',size=14)
 body(fig,.349,'Spring: conversion has equal and opposite signs in u and v; both have dilution from domain extension. Boundary positions a and b follow spring restoring motion plus activity-dependent outward protrusion P.',size=9)
 eq(fig,.270,r'v_t=\frac{10}{L^2}v_{qq}-f-\frac{\dot L}{L}v,\quad f=(4+S+H(u))v-\delta u',size=14)
 eq(fig,.203,r'\dot a=0.01(L-L_0)-P(u_l),\quad \dot b=-0.01(L-L_0)+P(u_r)',size=14)
 body(fig,.139,'With no-flux material endpoints, the conserved physical total is L times the integral of (u+v) over q. Native Spring length and time are uncalibrated; the displayed tension feedback is a deactivation-rate increment, not a calibrated force.',size=9)
 finish(pdf,fig,11)
 # 12: exact deBelly discrete closure and coordinate qualification.
 fig=page('deBelly: numerical closure and spatial measure','The quantified components follow the retained staggered finite-volume update')
 body(fig,.837,'A continuum-form transport identity is useful, but the exact numerical attribution uses lagged exchange and a uniform mass correction after transport. Let star denote the uncorrected state, p=5/6.4, N the cell count, and dt=0.5 seconds.',size=9)
 eq(fig,.752,r'\frac{b^*-b^{\rm old}}{\Delta t}=-\operatorname{div}_{\rm FV}F^*+5n^{\rm old}-1.4b^*',size=15)
 eq(fig,.684,r'c=\frac{1}{N}\sum_j\left[(b_j^{\rm old}+n_j^{\rm old})-(b_j^*+n_j^*)\right]',size=14)
 eq(fig,.616,r'b^{\rm new}=b^*+pc,\qquad n^{\rm new}=n^*+(1-p)c',size=15)
 body(fig,.552,'This correction restores total nominal MCA mass after boundary transport and staggered exchange. It is part of the retained numerical model. Its displacement contribution must not be interpreted as a separately established biochemical pathway.',size=9)
 eq(fig,.466,r'\dot m=\frac{\sum_j[\Delta x/2-\ell_j(m)]\,\dot b_j}{b_i}',size=15)
 body(fig,.401,'Here ell_j(m) is the overlap length of cell j with [0,m], and i is the crossing cell. Applying this functional separately to each update component gives the exact discrete decomposition; integrating the rates reproduces measured movement.',size=9)
 eq(fig,.314,r'X_{50}=m+d(m,t),\qquad\dot X_{50}=d_t(m,t)+[1+d_x(m,t)]\dot m',size=14)
 body(fig,.252,'This maps the conserved nominal/reference-density median to physical position using displacement d. Multiplying stored MCA by a deformation Jacobian would define a different mass diagnostic; it is not the conserved measure used in this analysis.',size=9)
 body(fig,.158,'Bound, unbound and total-MCA cutoff medians agree within 0.000032 times length, consistent with rapid exchange. The intrinsic exchange time is 1/(5+1.4)=0.15625 seconds, shorter than the retained timestep. Thus the exact discrete contribution is preferable to an assumed continuous rapid-equilibrium reduction.',size=9)
 finish(pdf,fig,12)
 # 13: mechanical control of MCA transport.
 fig=page('deBelly: what controls the transport velocity?','Mechanical screening, bounded protrusion, tension suppression, and staggered motion')
 body(fig,.837,'The cortical velocity v and membrane velocity w are different fields. Bound MCA sets friction density g; Rho gradients provide active stress. The mechanical solves use nominal coordinates, without a deformation Jacobian.',size=9)
 eq(fig,.755,r'g(v-w)=\eta v_{xx}+\sigma\rho_x,\quad g=\xi\max(b,10^{-10})',size=14)
 eq(fig,.693,r'\xi=100,\quad\eta=10{,}000,\quad\sigma=100,\quad\ell=\sqrt{\frac{\eta}{g}}',size=14)
 body(fig,.627,'At resting bound density about 0.078, the screening length ell is about 36 nominal microns, comparable to the two tested lengths (19 and 35). The ratio L/ell changes the velocity profile; length does more than extend an otherwise identical transport route.',size=9)
 eq(fig,.541,r'g^{\rm old}\frac{d^{\rm new}-d^{\rm old}}{\Delta t}-K d^{\rm new}_{xx}=g^{\rm old}v^{\rm old},\quad K=100',size=13)
 eq(fig,.474,r'w=\frac{d^{\rm new}-d^{\rm old}}{\Delta t},\qquad T=K d_x',size=14)
 body(fig,.411,'Rac activates outward protrusion through a bounded switch; endpoint tension suppresses it. The symmetric diagnostic implementation centers the cortical endpoint velocities. Membrane endpoint velocities follow the separate displacement update.',size=9)
 eq(fig,.328,r'Q(r)=\frac{1+\tanh[(r-1.3r_{\rm rest})/0.01]}{2}',size=14)
 eq(fig,.266,r'P_l=-\frac{0.15Q(r_l)}{1+T_l^2/10},\quad P_r=\frac{0.15Q(r_r)}{1+T_r^2/10}',size=14)
 eq(fig,.203,r'v_l=P_l-\overline{P},\quad v_r=P_r-\overline{P},\quad\overline{P}=(P_l+P_r)/2',size=13)
 body(fig,.140,'These laws identify bounded protrusion, tension inhibition, MCA-dependent friction and Rho stress gradients as transport controls. The separate raw edge update is d_edge,new=0.98*d_edge,old+dt*P_edge, followed by centering. This retained numerical relaxation is not a continuous elastic constitutive law.',size=9)
 finish(pdf,fig,13)
 # 14: the plateau evidence itself, rather than only its interpretation.
 fig=page('Starting profiles: broader chemical poles','Native concentrations are retained; length colors match the decomposition figures')
 body(fig,.837,'For WP, Holmes and Spring, the starting active/Rac peak changes little in height while its physical width increases. The normalized horizontal coordinate helps compare shapes, but physical width is listed explicitly below. This is evidence consistent with mesa saturation; it is not itself a saturation perturbation.',size=9)
 specs=[('wave_pinning',[9,25],'active','Wave-pinning: active'),('holmes_model3',[14,24],'rac_active','Holmes: active Rac'),('spring',[.8,2.4],'active','Spring: active'),('debelly',[19,35],'mca_bound','deBelly: bound MCA')]
 data_rows=[]
 for i,(model,lengths,field,title) in enumerate(specs):
  ax=fig.add_axes([.10+(i%2)*.47,.54-(i//2)*.255,.35,.15])
  vals=[]
  for L,color in zip(lengths,[BLUE,ORANGE]):
   path=VIEWER/f'{model}_physical_peak_matched_L{L:g}_full.npz';extra_sources.append(path)
   with np.load(path) as z:
    q=z['x']/L;u=z['field_'+field][0];height=float(u.max());base=float(u.min());width=float(np.trapezoid((u>=base+.5*(height-base)).astype(float),q)*L)
    ax.plot(q,u,color=color,label=f'L={L:g}');vals.append((height,width))
  ax.set(title=title,xlabel='Nominal x/L (Spring: material q)',ylabel='Native field value');ax.legend(fontsize=7);style(ax)
  if model!='debelly':data_rows.append([title.split(':')[0],f'{vals[0][0]:.3f} to {vals[1][0]:.3f}',f'{vals[0][1]:.3f} to {vals[1][1]:.3f}'])
 table(fig,[.09,.097,.82,.105],['Model','Starting peak: short to long','Half-height width: short to long'],data_rows,[.22,.36,.42])
 body(fig,.071,'Widths are nominal native lengths; Spring is uncalibrated. deBelly is shown as an MCA transport marker, not as evidence for a chemical saturation mechanism.',size=8,width=110,color=MUTED)
 finish(pdf,fig,14)

sources=[p for folder in [DATA,FREEZE] for p in folder.rglob('*') if p.suffix in ['.csv','.json','.txt','.py'] and '__pycache__' not in str(p)]
sources.extend([Path(__file__),ROOT/'polarity1d/manuscript/README.md',ROOT/'polarity1d/manuscript/PROTOCOL.md',ROOT/'polarity1d/manuscript/ASSUMPTIONS.md'])
sources.extend((ROOT/'polarity1d/four_assay/pilots/width_full_lengths_20261002/frozen_sources').rglob('*.py'))
sources.extend(extra_sources)
manifest=dict(generated_at=STAMP,pdf=str(PDF),pages=14,scope='Historical mechanism diagnostic, independent of current manuscript results',inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(sources))},pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(),equations='Vector matplotlib MathText with STIX fonts; body and all figures vector PDF',new_simulations=False)
(OUT/'provenance.json').write_text(json.dumps(manifest,indent=2))
print(PDF)
