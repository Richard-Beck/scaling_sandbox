"""Read-only native snapshots for the corrected linear deBelly width reversal."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'
DEFAULT=ROOT/'data/runs/revision_debelly_width_snapshots_20261019'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,a):p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def csvwrite(p,a):
    keys=list(dict.fromkeys(k for r in a for k in r))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(a)


def initialize(out):
    assert not out.exists();out.mkdir();(out/'frozen_sources').mkdir();(out/'results').mkdir();cases=[]
    m=read(BASE/'manifest.json');prov=read(BASE/'provenance.json')
    sources={str(BASE/'manifest.json'):sha(BASE/'manifest.json'),str(Path(__file__)):sha(Path(__file__))}
    for rel,h in prov['frozen_sources'].items():sources[str(BASE/'frozen_sources'/rel)]=h
    helper=ROOT/'code/mechanisms/revision_compensation_excess_analysis_20261019.py';sources[str(helper)]=sha(helper);shutil.copy2(helper,out/'frozen_sources'/helper.name)
    for c in m['cases']:
        if c['model']!='debelly' or c['length'] not in [18.,39.] or c['width'] not in [1.5,5.] or c['case_kind']!='constant_on':continue
        r=read(BASE/'workers'/f'{c["slot"]:04d}'/'result.json');assert r['case']==c
        w=next(w for w in r['withdrawals'] if w['duration']==2000.)
        cases.append(dict(case=c,on=r['on'],off=w))
    assert len(cases)==4
    write(out/'manifest.json',dict(native_base=str(BASE),cases=cases,source_hashes=sources,
        on_snapshot_times=[0.,50.,100.,150.,200.,300.,400.,500.,800.,2000.],off_snapshot_times=[0.,100.,1000.,10000.],
        native_state_law='Instantaneous conditional reaction/diffusion at saved states. Not a decomposition of staggered finite updates or an evolved pathway intervention.',
        cue_geometry='Input S is fixed on nominal solver x; physical x+displacement and endpoints exported separately.',no_native_integrations=True))
    frozen=out/'frozen_sources'/Path(__file__).name;shutil.copy2(Path(__file__),frozen)
    shell=out/'worker.sh';shell.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(frozen)+' worker --out '+str(out)+'\n');shell.chmod(0o750)
    cmd='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:10:00 --job-name=F2-DB-width-snapshots --output='+str(out/'%j.out')+' --error='+str(out/'%j.err')+' '+str(shell)
    (out/'submission_command.txt').write_text(cmd+'\n');print(sha(out/'manifest.json'));print(cmd)


def worker(out):
    m=read(out/'manifest.json')
    for path,h in m['source_hashes'].items():assert sha(Path(path))==h,path
    base=Path(m['native_base']);sys.path.insert(0,str(base/'frozen_sources/code/mechanisms'))
    import revision_five_model_compensation_20261019 as p
    import revision_compensation_excess_analysis_20261019 as excess
    parameters=p.n.model_for('debelly').parameters;assert parameters.center_boundary_motion is False
    profiles=[];regional=[];sources={};max_cue_error=0.;max_telescoping_error=0.
    for ref in m['cases']:
        c=ref['case'];rp=base/c['regional_input'];assert sha(rp)==c['regional_input_sha256']
        with np.load(rp) as zz:z={k:zz[k].copy() for k in zz.files}
        edge,orgw,full,initial=excess.initial_core(z,c);dx=c['length']/c['cells'];old=orgw;cuew=full*z['initial_cue_footprint'];distal=full*(z['material_fraction']>=.5)
        for phase,item,times in [('cue_on',ref['on'],m['on_snapshot_times']),('after_withdrawal',ref['off'],m['off_snapshot_times'])]:
            archive=base/item['archive'];assert sha(archive)==item['archive_sha256'];sources[str(archive)]=item['archive_sha256']
            with np.load(archive) as a:
                t=a['time'];ids=[]
                for tt in times:
                    hit=np.flatnonzero(t==tt);assert len(hit)==1;ids.append(int(hit[0]))
                f={name[6:]:a[name][ids].copy() for name in a.files if name.startswith('field_')}
                x=a['x'].copy()
            for j,tt in enumerate(times):
                u=f['rac'][j];rho=f['rho'][j];S=z['initial_signal'] if phase=='cue_on' else np.zeros_like(u)
                alpha=f['rac_activation'][j];delta=f['delta_alpha_mca'][j]
                err=float(np.max(abs(alpha-parameters.rac_basal_activation-delta-S)));max_cue_error=max(max_cue_error,err);assert err<2e-13
                cue_source=parameters.rac_deactivation*S/(1+rho*rho)
                mca_source=parameters.rac_deactivation*delta/(1+rho*rho)
                basal_source=parameters.rac_deactivation*parameters.rac_basal_activation/(1+rho*rho)
                sink=parameters.rac_deactivation*u
                # Exact cell-centred Neumann FV faces; conditional instantaneous law.
                faces=np.r_[0.,parameters.rac_deactivation*parameters.rac_diffusion*np.diff(u)/dx,0.]
                diffusion=np.diff(faces)/dx
                left=float(f['left_boundary_displacement'][j,0]);right=float(f['right_boundary_displacement'][j,0])
                corecontrast=float(max(u[:][orgw>0].max()-u.min(),0)/initial['initial_profile_contrast'])
                for name,w in [('organization_initial_core',old),('actual_nominal_cue_footprint',cuew),('fixed_distal_native_cells',distal)]:
                    # Fractional-weight telescoping sum, without dropping edge contributions.
                    diff_integral=float(diffusion@w);telescoped=float(faces[1:]@(w/dx)-faces[:-1]@(w/dx))
                    max_telescoping_error=max(max_telescoping_error,abs(diff_integral-telescoped))
                    regional.append(dict(length=c['length'],width=c['width'],phase=phase,native_time=tt,region=name,
                        raw_rac_amount=float(u@w),current_minimum_excess_amount=float(np.maximum(u-u.min(),0)@w),
                        native_peak=float(u[w>0].max()),mean_rho=float(rho@w/w.sum()),mean_mca_bound=float(f['mca_bound'][j]@w/w.sum()),
                        mean_tension=float(f['tension'][j]@w/w.sum()),mean_cortical_velocity=float(f['velocity'][j]@w/w.sum()),
                        conditional_cue_source=float(cue_source@w),conditional_MCA_source=float(mca_source@w),conditional_basal_source=float(basal_source@w),
                        conditional_rac_deactivation=float(sink@w),conditional_rac_diffusion=diff_integral,
                        conditional_net_rac_rate=float((cue_source+mca_source+basal_source-sink+diffusion)@w),
                        current_profile_floor=float(u.min()),core_peak_contrast_over_basal=corecontrast,
                        left_endpoint_displacement=left,right_endpoint_displacement=right,actual_physical_domain_length=c['length']+right-left,
                        actual_cue_sample_physical_left=float((x+f['displacement'][j])[z['initial_cue_footprint']].min()),
                        actual_cue_sample_physical_right=float((x+f['displacement'][j])[z['initial_cue_footprint']].max()),
                        source_archive=str(archive),source_archive_sha256=item['archive_sha256']))
                for k,xx in enumerate(x):
                    profiles.append(dict(length=c['length'],width=c['width'],phase=phase,native_time=tt,nominal_x=float(xx),
                        actual_physical_x=float(xx+f['displacement'][j,k]),nominal_cue_signal=float(S[k]),
                        **{name:float(f[name][j,k]) for name in ['rac','rho','mca_bound','mca_unbound','tension','velocity','displacement','rac_activation','rho_activation','delta_alpha_mca','delta_beta_tension']}))
    csvwrite(out/'results/native_snapshot_profiles.csv',profiles);csvwrite(out/'results/native_snapshot_regional_law.csv',regional)
    write(out/'results/validation.json',dict(manifest_sha256=sha(out/'manifest.json'),native_source_hashes=sources,
        artifact_hashes={str(q.relative_to(out)):sha(q) for q in (out/'results').glob('*.csv')},profile_rows=len(profiles),regional_rows=len(regional),
        maximum_cue_decomposition_error=max_cue_error,maximum_conditional_diffusion_telescoping_error=max_telescoping_error,
        numerical_scope='Saved-state conditional law, not an exact sequential update budget or causal pathway isolation.',no_native_integrations=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','worker']);parser.add_argument('--out',type=Path,default=DEFAULT)
    a=parser.parse_args()
    if a.action=='init':initialize(a.out.resolve())
    else:worker(a.out.resolve())
