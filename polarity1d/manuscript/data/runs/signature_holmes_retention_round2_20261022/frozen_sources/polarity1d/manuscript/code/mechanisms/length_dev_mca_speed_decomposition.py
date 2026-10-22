"""Length dependence of native bound-MCA centroid-speed contributions.

Independent transport and moment identities extend the verified MCA flow budget.
Three coarsely saved cue histories are replayed at the unchanged native step.
Workers run frozen sources under SLURM; original histories remain immutable.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, gzip, hashlib, json, os, shutil, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_mca_flow_budget as base
from length_dev_output_layout import dev_output_path

ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='manuscript')
REPORT=ROOT/'figures/length_mechanisms_20261008'
PARENT=ROOT/'data/runs/length_narrative_memory_20261009'
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
RUN=ROOT/'data/runs/debelly_mca_speed_decomposition_20261010'
OUT=REPORT/'dev/mca_speed_decomposition'
STEM='MCA_speed_decomposition_across_length'
ALL_LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
LENGTHS=ALL_LENGTHS.copy()
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
SOURCES={}
transport,flux_parts,commands,solve=base.transport,base.flux_parts,base.commands,base.solve

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(p):
 p=Path(p);SOURCES[str(p.relative_to(ROOT))]=sha(p);return p
def read(p):return json.loads(source(p).read_text())
def write(p,value):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def csvwrite(p,rows):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def now():return datetime.now(timezone.utc).isoformat()
def output(name):return dev_output_path(OUT,name)

def collect():
    manifest = read(PARENT/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(PARENT/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    parameter_meta = read(REPORT/'dev/resting_state/supporting/metadata/basal_endpoint_metrics_vs_length.manifest.json')
    p = parameter_meta['parameters']
    assert p['center_boundary_motion'] is False
    endpoint = source(REPORT/'dev/deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv')
    with endpoint.open() as stream:
        thresholds = {int(r['length_um']):float(r['rac_switch_threshold']) for r in csv.DictReader(stream)}
    with gzip.open(source(BASE/'results/records.json.gz'), 'rt') as stream:
        records = {r['case_index']:r for r in json.load(stream)['cases']}
    rows, audits, snapshots, trajectories = [], [], {}, {}
    fb = p['mca_on']/(p['mca_on']+p['mca_off'])
    for L in LENGTHS:
        job = next(j for j in manifest['jobs'] if j['case']['length']==L and j['case']['condition']=='fixed_relative')
        case, D = job['case'], float(job['duration'])
        record = records[case['case_index']]
        assert D == record['result']['full_polarity']['endpoint']
        assert D == min(t['duration'] for t in record['tests'] if t['full_polarity_success'])
        assert case['response_amplitude']==3.6 and case['width_mode']=='relative'
        path = ROOT/job['source'] if job['source'] else PARENT/f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        meta = read(path.with_suffix('.json'))
        assert sha(source(path)) == meta['archive_sha256']
        checkpoint = source(ROOT/case['checkpoint'])
        assert sha(checkpoint)==case['checkpoint_sha256']
        with np.load(path) as z:
            x, t = z['x'].copy(), z['pulse_time'].copy()
            fields = {k[len('pulse_field_'):]:z[k].copy() for k in z.files if k.startswith('pulse_field_')}
            with np.load(checkpoint) as cp:
                keys = [k[6:] for k in cp.files if k.startswith('state_')]
                assert len(keys)==11
                for k in keys:
                    np.testing.assert_array_equal(fields[k][0], cp['state_'+k])
                    np.testing.assert_array_equal(fields[k][-1], z['recovery_field_'+k][0])
            tail = z['recovery_time']>=9500
            rac = z['recovery_field_rac'][tail]
            orientation = rac@(2*x/L-1)/rac.sum(axis=1)
            strength_error = float(np.max(abs(orientation/abs(case['initial_orientation'])-1)))
            assert orientation.min()>.05 and np.ptp(rac,axis=1).min()>=.05 and np.ptp(orientation)<.005
            assert strength_error<=.01 and meta['metrics']['full_polarity_success']
        original_t, original_c = t.copy(), fields['mca_bound']@x/fields['mca_bound'].sum(axis=1)
        dense_replay_error = 0.
        dense_replayed = bool(np.diff(t).max()>.5+1e-9)
        if dense_replayed:
            # Recover original half-second steps, checking every retained sample.
            from dataclasses import replace
            import sys
            sys.path[:0] = [str(PARENT/'frozen_sources/code'),str(PARENT/'frozen_sources/code/core')]
            import native as n
            import polarity_establishment as establishment
            native_p = n.ds.DeBellyParameters(**p)
            assert not native_p.center_boundary_motion
            n.model_for = lambda key:n.ds.DeBelly(native_p)
            initial = {key:fields[key][0].copy() for key in n.DB_NAMES}
            cue = establishment.Patch(case['response_amplitude'],case['width'],side='right',width_mode=case['width_mode'],fraction=.2)
            tr = n.run('debelly',L,initial,cue,D,interval=.5,n=case['cells'])
            ids = np.searchsorted(tr.time,t)
            assert np.max(abs(tr.time[ids]-t))<1e-8
            for key in fields:
                error = float(abs(tr.fields[key][ids]-fields[key]).max())
                dense_replay_error = max(dense_replay_error,error)
                assert error<1e-10,(L,key,error)
            t,fields = tr.time,{key:tr.fields[key] for key in fields}
            for module in list(sys.modules.values()):
                name=getattr(module,'__file__',None)
                if name and Path(name).suffix=='.py' and PARENT/'frozen_sources' in Path(name).parents:
                    path_native=source(Path(name))
                    relative=str(path_native.relative_to(PARENT/'frozen_sources'))
                    if relative in manifest['source_sha256']:
                        assert sha(path_native)==manifest['source_sha256'][relative]
        assert t[0]==0 and abs(t[-1]-D)<1e-9 and np.all((np.diff(t)>0)&(np.diff(t)<=.5+1e-9))
        np.testing.assert_allclose(x, (np.arange(len(x))+.5)*L/len(x), atol=1e-12)
        dx = L/len(x)
        c = fields['mca_bound']@x/fields['mca_bound'].sum(axis=1)
        front = x>=.8*L
        total_mass = float(np.sum(fields['mca_bound'][0]+fields['mca_unbound'][0]))
        errors = dict(pool_reconstruction=0.,local_bound_budget=0.,local_total_budget=0.,centroid_budget=0.,
                      total_mass=0.,flow_superposition=0.,uniform_source=0.,flux_telescoping=0.,mechanical_transport_budget=0.,
                      independent_moment_formula=0.)
        first = {k:v[0] for k,v in fields.items()}
        initial_edge = commands(first,p,thresholds[L])
        initial_friction = p['linker_friction']*np.maximum(first['mca_bound'],1e-10)
        old_component = dict(left_drive=solve(initial_friction,p['cortical_viscosity'],np.zeros_like(x),dx,initial_edge[0],0.),
                             right_drive=solve(initial_friction,p['cortical_viscosity'],np.zeros_like(x),dx,0.,initial_edge[1]),
                             rho=solve(initial_friction,p['cortical_viscosity'],p['active_stress']*p['rho_to_rac_coupling']*np.gradient(first['rho'],dx),dx),
                             membrane=np.zeros_like(x))
        initial_remainder = first['velocity']-sum(old_component.values())
        assert abs(initial_remainder).max()<1e-5
        old_component['initial_remainder'] = initial_remainder
        for i in range(1, len(t)):
            dt = float(t[i]-t[i-1])
            old = {k:v[i-1] for k,v in fields.items()}
            new = {k:v[i] for k,v in fields.items()}
            edge = commands(old,p,thresholds[L])
            vm = (new['displacement']-old['displacement'])/dt
            me = [(new[k][0]-old[k][0])/dt for k in ['left_boundary_displacement','right_boundary_displacement']]
            bpre = transport(old['mca_bound'],old['velocity'],edge,p['mca_bound_diffusion'],p['mca_off'],p['mca_on']*old['mca_unbound'],dt,dx)
            upre = transport(old['mca_unbound'],vm,me,p['mca_unbound_diffusion'],p['mca_on'],p['mca_off']*bpre,dt,dx)
            correction = (total_mass-(bpre+upre).sum())/len(x)
            bfinal, ufinal = bpre+fb*correction, upre+(1-fb)*correction
            errors['pool_reconstruction'] = max(errors['pool_reconstruction'],float(abs(bfinal-new['mca_bound']).max()),float(abs(ufinal-new['mca_unbound']).max()))
            flux = flux_parts(bpre,old['velocity'],edge,p['mca_bound_diffusion'],dx)
            rate = {k:-np.diff(v)/dx for k,v in flux.items()}
            rate['exchange'] = p['mca_on']*old['mca_unbound']-p['mca_off']*bpre
            rate['uniform_replenishment'] = np.full(len(x),fb*correction/dt)
            actual = (new['mca_bound']-old['mca_bound'])/dt
            errors['local_bound_budget'] = max(errors['local_bound_budget'],float(abs(sum(rate.values())-actual).max()))
            fu = flux_parts(upre,vm,me,p['mca_unbound_diffusion'],dx)
            unbound_rate = -np.diff(sum(fu.values()))/dx+p['mca_off']*bpre-p['mca_on']*upre+(1-fb)*correction/dt
            observed_total = (new['mca_bound']+new['mca_unbound']-old['mca_bound']-old['mca_unbound'])/dt
            errors['local_total_budget'] = max(errors['local_total_budget'],float(abs(sum(rate.values())+unbound_rate-observed_total).max()))
            errors['total_mass'] = max(errors['total_mass'],float(abs((new['mca_bound']+new['mca_unbound']).sum()-total_mass)*dx))
            errors['uniform_source'] = max(errors['uniform_source'],float(np.ptp((new['mca_bound']-bpre)/dt)))
            total_flux = sum(flux.values())+sum(fu.values())
            errors['flux_telescoping'] = max(errors['flux_telescoping'],float(abs((-np.diff(total_flux)).sum()-(total_flux[0]-total_flux[-1]))))
            # Exact finite-step changing-mass centroid identity. Positive speed
            # means the bound-MCA centroid moves left, toward the new rear.
            weights = -(x-c[i-1])/new['mca_bound'].sum()
            speeds = {k:float(v@weights) for k,v in rate.items()}
            alternate = dict(interior_advection=-float(flux['interior_advection'][1:-1].sum()/new['mca_bound'].sum()),
                             old_end_exit=-float((x[0]-c[i-1])*flux['old_end_exit'][0]/(dx*new['mca_bound'].sum())),
                             new_end_exit=float((x[-1]-c[i-1])*flux['new_end_exit'][-1]/(dx*new['mca_bound'].sum())))
            errors['independent_moment_formula'] = max(errors['independent_moment_formula'],max(abs(alternate[k]-speeds[k]) for k in alternate))
            # Freeze donor selection on the actual full transport velocity.
            # Interior velocity is staggered (previous cortical solve); endpoint
            # transport commands use the current step's Rac/tension inputs.
            actual_vf = .5*(old['velocity'][:-1]+old['velocity'][1:])
            donors = np.where(actual_vf>=0,bpre[:-1],bpre[1:])
            mechanical_rate = {}
            for key,v in old_component.items():
                vf = .5*(v[:-1]+v[1:])
                ff = np.r_[edge[0]*bpre[0] if key=='left_drive' else 0.,vf*donors,
                           edge[1]*bpre[-1] if key=='right_drive' else 0.]
                mechanical_rate[key] = -np.diff(ff)/dx
            errors['mechanical_transport_budget'] = max(errors['mechanical_transport_budget'],
                float(abs(sum(mechanical_rate.values())-sum(rate[k] for k in ['interior_advection','old_end_exit','new_end_exit'])).max()))
            observed_speed = float((c[i-1]-c[i])/dt)
            errors['centroid_budget'] = max(errors['centroid_budget'],abs(sum(speeds.values())-observed_speed))
            friction = p['linker_friction']*np.maximum(new['mca_bound'],1e-10)
            current_split = dict(left_drive=solve(friction,p['cortical_viscosity'],np.zeros_like(x),dx,edge[0],0.),
                                 right_drive=solve(friction,p['cortical_viscosity'],np.zeros_like(x),dx,0.,edge[1]))
            component = dict(endpoint=current_split['left_drive']+current_split['right_drive'],
                             rho=solve(friction,p['cortical_viscosity'],p['active_stress']*p['rho_to_rac_coupling']*np.gradient(old['rho'],dx),dx),
                             membrane=solve(friction,p['cortical_viscosity'],friction*vm,dx))
            errors['flow_superposition'] = max(errors['flow_superposition'],float(abs(sum(component.values())-new['velocity']).max()))
            common = dict(length_um=L,time_start_s=float(t[i-1]),time_end_s=float(t[i]),step_s=dt,
                          cue_duration_s=D,centroid_um=float(c[i]),centroid_speed_um_s=observed_speed,
                          new_front_bound_MCA=float(new['mca_bound'][front].mean()),
                          new_endpoint_bound_MCA=float(new['mca_bound'][-1]),
                          mean_cortical_velocity_um_s=float(new['velocity'].mean()),
                          mean_transport_cortical_velocity_um_s=float(old['velocity'].mean()),
                          left_command_um_s=float(edge[0]),right_command_um_s=float(edge[1]),
                          right_bound_outflux=float(flux['new_end_exit'][-1]),
                          left_bound_outflux=float(-flux['old_end_exit'][0]),
                          bound_mass=float(new['mca_bound'].sum()*dx),
                          bound_uniform_source_conc_s=float(fb*correction/dt),
                          native_field_tail_pass=bool(meta['metrics']['field_tail_pass']))
            common.update({'flow_'+k+'_um_s':float(v.mean()) for k,v in component.items()})
            common.update({'speed_'+k+'_um_s':v for k,v in speeds.items()})
            common.update({'front_depletion_'+k+'_conc_s':-float(v[front].mean()) for k,v in rate.items()})
            common.update({'mechanical_speed_'+k+'_um_s':float(v@weights) for k,v in mechanical_rate.items()})
            common.update({'mechanical_front_depletion_'+k+'_conc_s':-float(v[front].mean()) for k,v in mechanical_rate.items()})
            common.update({'flow_'+k+'_um_s':float(v.mean()) for k,v in current_split.items()})
            rows.append(common)
            old_component = dict(**current_split,rho=component['rho'],membrane=component['membrane'],initial_remainder=np.zeros_like(x))
            if abs(t[i]-60)<1e-9 or abs(t[i]-180)<1e-9:
                label = f'L{L}_t{int(t[i])}'
                snapshots.update({label+'_x':x,label+'_bound':new['mca_bound'],label+'_bound_before':old['mca_bound'],
                                  label+'_velocity':new['velocity'],label+'_transport_velocity':old['velocity'],
                                  label+'_edge':edge,label+'_right_tension':np.array(old['tension'][-1])})
                snapshots.update({label+'_flow_'+k:v for k,v in component.items()})
                snapshots.update({label+'_rate_'+k:v for k,v in rate.items()})
        assert errors['pool_reconstruction']<1e-12,errors
        assert errors['flow_superposition']<1e-11,errors
        assert max(errors[k] for k in errors if k!='flow_superposition')<1e-10,errors
        trajectories.update({f'L{L}_original_time':original_t,f'L{L}_original_centroid':original_c,f'L{L}_time':t,f'L{L}_centroid':c,f'L{L}_bound':fields['mca_bound'],
                             f'L{L}_velocity':fields['velocity'],f'L{L}_tension':fields['tension'],
                             f'L{L}_x':x,f'L{L}_rac':fields['rac'],f'L{L}_rho':fields['rho']})
        audits.append(dict(length_um=L,archive=str(path.relative_to(ROOT)),checkpoint=str(checkpoint.relative_to(ROOT)),
                           cue_duration_s=D,strength_error=strength_error,native_field_tail_pass=bool(meta['metrics']['field_tail_pass']),
                           intervals=len(t)-1,dense_replayed=dense_replayed,dense_replay_max_error=dense_replay_error,errors=errors))
        print('Verified',L,errors,flush=True)
    return rows,audits,snapshots,trajectories,p


def launch():
 assert not RUN.exists(),'Preserve existing evidence; choose a fresh run.'
 for folder in ['workers','logs','frozen_sources/code/mechanisms']:(RUN/folder).mkdir(parents=True,exist_ok=True)
 files=[Path(__file__),Path(base.__file__),Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py'),Path(__file__).with_name('length_dev_output_layout.py')]
 fingerprints={}
 for p in files:
  dest=RUN/'frozen_sources/code/mechanisms'/p.name;shutil.copy2(p,dest);fingerprints[str(dest.relative_to(ROOT))]=sha(dest)
 write(RUN/'manifest.json',dict(generated_at=now(),lengths_um=ALL_LENGTHS,frozen_source_sha256=fingerprints,geometry='fixed_relative',native_step_s=.5))
 write(RUN/'status.json',dict(generated_at=now(),producer='running',complete=False))
 import shlex
 cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([PYTHON,str(RUN/'frozen_sources/code/mechanisms'/Path(__file__).name),'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"'
 job=subprocess.check_output(['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G','--time=00:10:00','--array=0-11','--job-name=mca-speed',
  '--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err'),'--wrap',cmd],text=True).strip()
 write(RUN/'launch.json',dict(generated_at=now(),array=job));print(job,flush=True)

def worker(slot):
 global LENGTHS
 manifest=json.loads((RUN/'manifest.json').read_text())
 for name,h in manifest['frozen_source_sha256'].items():assert sha(ROOT/name)==h
 LENGTHS=[ALL_LENGTHS[slot]]
 for p in [Path(__file__),Path(base.__file__),Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py'),Path(__file__).with_name('length_dev_output_layout.py')]:source(p)
 rows,audits,snapshots,trajectories,p=collect()
 target=RUN/'workers'/f'{slot:03d}'
 csvwrite(target.with_suffix('.csv'),rows)
 np.savez_compressed(target.with_suffix('.npz'),**trajectories)
 write(target.with_suffix('.json'),dict(generated_at=now(),passed=True,checks=audits,parameters=p,source_sha256=SOURCES,
  csv_sha256=sha(target.with_suffix('.csv')),npz_sha256=sha(target.with_suffix('.npz'))))
 print('Finished',LENGTHS[0],len(rows),flush=True)

def summarize(rows,trajectories):
 result,checks=[],[]
 zpath=source(REPORT/'dev/deBelly_redraft/supporting/data/Figure_Z_MCA_distance_speed.csv')
 zmeta=read(zpath.parent.parent/'metadata/Figure_Z_MCA_distance_speed.manifest.json')
 assert sha(zpath)==zmeta['sha256']['supporting/data/'+zpath.name]
 with zpath.open() as f:zrows={int(r['length_um']):r for r in csv.DictReader(f) if r['condition']=='fixed_relative'}
 for L in ALL_LENGTHS:
  rr=[r for r in rows if int(r['length_um'])==L]
  t=trajectories[f'L{L}_time'];c=trajectories[f'L{L}_centroid'];original_t=trajectories[f'L{L}_original_time'];original_c=trajectories[f'L{L}_original_centroid']
  starts=np.array([r['time_start_s'] for r in rr]);ends=np.array([r['time_end_s'] for r in rr]);dt=ends-starts;D=ends[-1]
  # Partial saved intervals use their own time-averaged native-step budget,
  # exactly matching Figure Z's interpolation of the original saved centroid.
  indices=np.searchsorted(original_t,(starts+ends)/2,side='right')-1
  for name,start,end,zkey in [('early',0.,60.,'early_mean_speed_um_s'),('late',D-60.,D,'late_mean_speed_um_s'),('whole',0.,D,'cue_average_speed_um_s')]:
   fraction=np.maximum(0.,np.minimum(original_t[indices+1],end)-np.maximum(original_t[indices],start))/(original_t[indices+1]-original_t[indices])
   w=dt*fraction
   fine=np.maximum(0.,np.minimum(ends,end)-np.maximum(starts,start))
   assert abs(w.sum()-(end-start))<1e-8
   row=dict(length_um=L,window=name,start_s=float(start),end_s=float(end),cue_duration_s=float(D),native_field_tail_pass=bool(rr[0]['native_field_tail_pass']))
   fields=[k for k in rr[0] if k.startswith(('speed_','mechanical_speed_')) or k=='centroid_speed_um_s']
   for key in fields:row[key]=float(w@np.array([r[key] for r in rr])/(end-start))
   row['mechanical_speed_endpoint_um_s']=row['mechanical_speed_left_drive_um_s']+row['mechanical_speed_right_drive_um_s']
   row['speed_advection_um_s']=sum(row['speed_'+k+'_um_s'] for k in ['interior_advection','old_end_exit','new_end_exit'])
   observed=(np.interp(start,original_t,original_c)-np.interp(end,original_t,original_c))/(end-start)
   transport=sum(row['speed_'+k+'_um_s'] for k in ['interior_advection','old_end_exit','new_end_exit','diffusion','exchange','uniform_replenishment'])
   mechanical=sum(row['mechanical_speed_'+k+'_um_s'] for k in ['left_drive','right_drive','rho','membrane','initial_remainder'])+sum(row['speed_'+k+'_um_s'] for k in ['diffusion','exchange','uniform_replenishment'])
   fine_observed=(np.interp(start,t,c)-np.interp(end,t,c))/(end-start)
   error=max(abs(observed-row['centroid_speed_um_s']),abs(observed-transport),abs(observed-mechanical),abs(observed-float(zrows[L][zkey])))
   assert error<1e-10,(L,name,error)
   row['observed_speed_per_length_s']=float(observed/L)
   row['native_halfstep_window_speed_um_s']=float(fine_observed)
   row['original_vs_halfstep_window_difference_um_s']=float(observed-fine_observed)
   checks.append(dict(length_um=L,window=name,maximum_integrated_speed_identity_error_um_s=float(error),
    original_vs_native_halfstep_window_difference_um_s=float(abs(observed-fine_observed))))
   result.append(row)
 return result,checks

MECHANICAL=[('mechanical_speed_endpoint_um_s','Endpoint-driven advection','#D55E00','-'),('mechanical_speed_rho_um_s','Rho-driven advection','#0072B2','-'),
 ('mechanical_speed_membrane_um_s','Membrane-driven advection','#777777','-'),('speed_diffusion_um_s','Bound diffusion','#56B4E9',':'),
 ('speed_exchange_um_s','Binding exchange','#009E73','-'),('speed_uniform_replenishment_um_s','Uniform replacement','#CC79A7','-'),
 ('centroid_speed_um_s','Observed speed','#202020','--')]
ROUTES=[('speed_interior_advection_um_s','Interior advection','#D55E00','-'),('speed_old_end_exit_um_s','Old-end loss','#8C564B','-'),
 ('speed_new_end_exit_um_s','New-end loss','#0072B2','-'),('speed_diffusion_um_s','Bound diffusion','#56B4E9',':'),
 ('speed_exchange_um_s','Binding exchange','#009E73','-'),('speed_uniform_replenishment_um_s','Uniform replacement','#CC79A7','-'),
 ('centroid_speed_um_s','Observed speed','#202020','--')]

def render(summary):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig=plt.figure(figsize=(14.5,9.5))
 grid=fig.add_gridspec(2,3,left=.08,right=.985,bottom=.13,top=.82,hspace=.75,wspace=.16)
 for row,specs in enumerate([MECHANICAL,ROUTES]):
  values=[r[k]*1000 for r in summary for k,*_ in specs]
  lo,hi=min(values),max(values);pad=.07*(hi-lo)
  for col,(window,title) in enumerate([('early','First 60 s of cue'),('late','Final 60 s before withdrawal'),('whole','Whole selected cue')]):
   ax=fig.add_subplot(grid[row,col]);rr=[r for r in summary if r['window']==window]
   for key,label,color,ls in specs:
    ax.plot(ALL_LENGTHS,[r[key]*1000 for r in rr],color=color,ls=ls,lw=2.3 if key=='centroid_speed_um_s' else 1.65,marker='o' if key=='centroid_speed_um_s' else None,ms=3.5)
   ax.axhline(0,color='#AAAAAA',lw=.7,zorder=0);ax.grid(axis='y',alpha=.15)
   ax.set(xlim=(17.5,39.5),ylim=(lo-pad,hi+pad),xticks=[18,21,25,29,33,37,39],xlabel='Cell length (µm)')
   if col==0:ax.set_ylabel('Contribution to centroid speed\n(10⁻³ µm/s; positive = leftward)')
   else:ax.tick_params(labelleft=False)
   if row==0:ax.set_title(title,fontsize=12,pad=12)
   ax.text(-.09,1.04,'ABCDEF'[row*3+col],transform=ax.transAxes,fontsize=13,weight='bold')
 fig.suptitle('MCA redistribution speed across length',fontsize=19,y=.99)
 fig.text(.53,.948,'Relative-20% cue · all twelve eligible lengths · bound-MCA centroid in nominal coordinates',ha='center',fontsize=10)
 for specs,y,title in [(MECHANICAL,.912,'Advection separated by mechanical forcing'),(ROUTES,.465,'The same speed separated into interior transport and boundary loss')]:
  fig.text(.08,y,title,fontsize=12,weight='bold')
  handles=[Line2D([],[],color=c,ls=ls,lw=2 if k=='centroid_speed_um_s' else 1.6,label=label) for k,label,c,ls in specs]
  fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,y-.008),ncol=4,frameon=False,fontsize=9,columnspacing=1.8)
 fig.text(.08,.055,'Contributions in each panel sum to observed speed. Positive shifts MCA toward the new left rear; negative shifts it toward the challenged right front.',fontsize=9)
 fig.text(.08,.032,'Early windows share elapsed time; late and whole-cue windows follow each cell’s own selected duration. Mechanical contributions are instantaneous source allocations.',fontsize=9)
 fig.text(.08,.011,'Both rows describe the same centroid budget. They are alternative decompositions and must not be added together.',fontsize=9)
 for ext in ['.png','.pdf','.svg']:fig.savefig(output(STEM+ext),dpi=220,facecolor='white')
 plt.close(fig)

def publish():
 source(Path(__file__))
 manifest=read(RUN/'manifest.json')
 rows,audits,trajectories=[],[],{}
 for slot,L in enumerate(ALL_LENGTHS):
  path=RUN/'workers'/f'{slot:03d}';meta=json.loads(path.with_suffix('.json').read_text());assert meta['passed']
  assert sha(path.with_suffix('.csv'))==meta['csv_sha256'] and sha(path.with_suffix('.npz'))==meta['npz_sha256']
  for name,h in meta['source_sha256'].items():assert sha(ROOT/name)==h;SOURCES[name]=h
  source(path.with_suffix('.json'));source(path.with_suffix('.csv'));source(path.with_suffix('.npz'))
  with path.with_suffix('.csv').open() as f:
   for r in csv.DictReader(f):
    rows.append({k:(v=='True' if k=='native_field_tail_pass' else float(v)) for k,v in r.items()})
  audits.extend(meta['checks'])
  with np.load(path.with_suffix('.npz')) as z:trajectories.update({k:z[k].copy() for k in z.files})
 summary,checks=summarize(rows,trajectories)
 render(summary)
 csvwrite(output(STEM+'_window_summary.csv'),summary)
 csvwrite(output(STEM+'_native_interval_budgets.csv'),rows)
 generated=now()
 caption=('Bound-MCA centroid-speed decomposition across all twelve eligible corrected deBelly lengths (18,19,21,23,25,27,29,31,33,35,37,39 µm). '
 'Columns show the first 60 s, final 60 s and whole selected relative-20% cue. Positive speed is leftward toward the new rear; units are 10^-3 µm/s. '
 'A–C partition advection by endpoint forcing, direct Rho stress-gradient forcing and membrane-motion drag, alongside bound diffusion, binding exchange, spatially uniform mass correction and observed speed. '
 'D–F partition that same advection into interior transport, old-left-end loss and new-right-end loss; diffusion, exchange and uniform replacement are unchanged. The two rows are alternative decompositions, not additive budgets. '
 'Every contribution is evaluated using independently reconstructed native pre-closure bound/unbound pools. The changing-bound-mass centroid identity is exact at each native step. '
 'Cortical transport uses the preceding cortical solution with the current endpoint commands. Mechanical components use identical donors chosen by the full native velocity, so their flux contributions add exactly. '
 'Binding exchange is the bound-pool on/off term; it is not the numerical split-exchange residual in a total-MCA budget. Uniform replacement is the original native mass correction. '
 'The initial stationary mechanical reconstruction has a negligible exported remainder included in closure checks. '
 'Exact original checkpoints, eleven-field phase continuity, shortest tested scalar-1%-strength selection and independent final-window scalar gates are verified. Field-tail flags are retained separately. '
 'Half-second cue histories are reused where available; L=18,37,39 are replayed at the original 0.5-s native step and checked at every original sample for every recorded field. '
 'To reproduce Figure Z exactly, a partial original saved interval is assigned its mean native-step contribution. The difference from using the recovered half-step interpolation is explicitly exported and audited. '
 'Early/late/whole totals match Figure Z. Late and full-cue windows use length-dependent cue durations and are descriptive, not independent predictions of those durations. '
 'Lines connect sampled lengths without fitting or declaring length-regime boundaries. These are source allocations on the actual coupled trajectories, not independently evolved causal ablations. '
 'Diffusion reflects at both domain endpoints. Native simulation evidence and the main report are unchanged.\n')
 output(STEM+'_caption.txt').write_text(caption)
 audit=dict(generated_at=generated,producer='finished',passed=True,lengths_um=ALL_LENGTHS,intervals=len(rows),jobs=audits,window_checks=checks,
  maximum_step_centroid_budget_error_um_s=max(a['errors']['centroid_budget'] for a in audits),
  maximum_integrated_speed_identity_error_um_s=max(c['maximum_integrated_speed_identity_error_um_s'] for c in checks),
  maximum_original_vs_halfstep_window_difference_um_s=max(c['original_vs_native_halfstep_window_difference_um_s'] for c in checks),
  parameters=meta['parameters'])
 write(output(STEM+'_audit.json'),audit)
 html='<!doctype html><html><head><meta charset="utf-8"><title>MCA speed decomposition</title><style>body{font:16px system-ui;max-width:1600px;margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.5}</style></head><body><h1>MCA redistribution speed across length</h1><p>Generated '+generated+' · producer finished · all twelve lengths verified.</p><p><a href="../vectors/'+STEM+'.pdf">PDF</a> · <a href="../data/'+STEM+'_window_summary.csv">Window contributions</a> · <a href="../metadata/'+STEM+'_audit.json">Audit</a></p><img src="../../'+STEM+'.png" alt="Six panels decompose MCA centroid speed across cell length"><p>'+caption+'</p></body></html>'
 output(STEM+'.html').write_text(html)
 (OUT/'README.md').write_text('# MCA redistribution speed across length\n\n[All development experiments](../README.md)\n\nRelative-20% cue at all twelve eligible lengths. Columns compare early, late and whole-cue averages. The top row separates advection by mechanical forcing; the bottom row separates interior transport from boundary loss. Each budget sums to the actual bound-MCA centroid speed.\n\n- [PNG]('+STEM+'.png) · [PDF](supporting/vectors/'+STEM+'.pdf) · [reader](supporting/viewers/'+STEM+'.html)\n- [Window summary](supporting/data/'+STEM+'_window_summary.csv) · [Native-step contributions](supporting/data/'+STEM+'_native_interval_budgets.csv)\n- [Caption](supporting/metadata/'+STEM+'_caption.txt) · [Audit](supporting/metadata/'+STEM+'_audit.json)\n')
 artifacts=[output(STEM+ext) for ext in ['.png','.pdf','.svg','.html','_window_summary.csv','_native_interval_budgets.csv','_caption.txt','_audit.json']]
 write(output(STEM+'.manifest.json'),dict(generated_at=generated,passed=True,producer='finished',artifact_root='../..',lengths_um=ALL_LENGTHS,
  sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},source_sha256=SOURCES,reproduction='python polarity1d/manuscript/code/mechanisms/length_dev_mca_speed_decomposition.py render'))
 write(RUN/'status.json',dict(generated_at=generated,producer='finished',complete=True,audit_passed=True))
 print(json.dumps({k:v for k,v in audit.items() if k.startswith('maximum')},indent=2));print(output(STEM+'.png'),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','worker','render']);ap.add_argument('--slot',type=int)
 a=ap.parse_args()
 if a.action=='launch':launch()
 elif a.action=='worker':worker(a.slot)
 else:publish()
