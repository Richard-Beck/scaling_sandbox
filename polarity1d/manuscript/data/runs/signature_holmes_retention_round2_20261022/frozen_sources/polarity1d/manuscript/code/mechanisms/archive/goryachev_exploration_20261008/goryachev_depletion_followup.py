"""Tight feedback availability and length-matched inactive transport controls."""
from pathlib import Path
from dataclasses import dataclass
import argparse, copy, json, shlex, shutil, subprocess, sys, traceback
import numpy as np
import goryachev_depletion as g
from models.goryachev import Goryachev

OUT=g.n.DATA/'runs/goryachev_depletion_followup_20261008'
BASE=g.n.DATA/'runs/goryachev_depletion_20261008_v2'
KINDS=['transport_matched','both_initial_transport','feedback_tight']


@dataclass(frozen=True)
class ExtendedControl(g.AvailabilityControl):
    @property
    def diffusivities(self):
        return self.du, self.dv*((self.length/9.)**2 if self.kind in KINDS[:2] else 1.)

    def reaction(self,time,x,state,stimulus):
        if self.kind=='transport_matched':return Goryachev.reaction(self,time,x,state,stimulus)
        u,v=state;cv,fv=v,v
        if self.kind=='both_initial_transport':
            factor=self.reference_v/self.initial_v
            cv=factor*v;fv=v*(1.+self.challenger_weight(x)*(factor-1.))
        elif self.kind=='feedback_tight':
            fv=v+self.challenger_weight(x)*(self.effective(v)-v)
        else:raise ValueError(self.kind)
        transfer=self.a*self.ec*u*u*fv+self.b*self.ec*u*fv-self.c*u
        transfer+=stimulus*cv
        return np.vstack((transfer,-transfer))


def submit(out,action,array=None,dependency=None):
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),action,'--out',str(out)])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00',
          '--job-name=gory-depletion-followup-'+action,'--chdir='+str(g.n.ROOT),
          '--output='+str(out/'logs'/f'{action}_%A_%a.out'),'--error='+str(out/'logs'/f'{action}_%A_%a.err')]
    if array:opts+=['--array='+array]
    if dependency:opts+=['--dependency=afterany:'+dependency]
    return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]


def initialize(out):
    if out.exists():raise ValueError('Use fresh destination')
    assert g.r.read(BASE/'results/status.json')['complete']
    base=g.verify(BASE);m=copy.deepcopy(base);jobs=[]
    for old in base['jobs']:
        if old['kind']!='native':continue
        for kind in KINDS:
            job=copy.deepcopy(old);job.update(slot=len(jobs),kind=kind,epsilon_fraction=1e-5)
            jobs.append(job)
    m.update(jobs=jobs,kinds=KINDS,created_at=g.n.now(),base=str(BASE.relative_to(g.n.ROOT)),
             intervention='During cue only: transport_matched sets inactive diffusion to 10*(L/9)^2, matching L^2/Dv to L=9. both_initial_transport adds initial-substrate compensation of cue and challenger feedback. feedback_tight uses reference-substrate feedback with epsilon/reference=1e-5 and native cue/diffusion. Native recovery is restored after every pulse.',
             epsilon_fraction=1e-5)
    for folder in ['logs','workers','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    g.r.write(out/'manifest.json',m)
    hashes=g.r.read(BASE/'provenance.json')['sources']
    hashes[str(Path(__file__).relative_to(g.n.ROOT))]=g.n.sha(Path(__file__))
    for name in hashes:
        dest=out/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(g.n.ROOT/name,dest)
    g.r.write(out/'provenance.json',dict(sources=hashes,manifest_sha256=g.n.sha(out/'manifest.json')))
    g.r.write(out/'results/status.json',dict(generated_at=g.n.now(),producer='running',complete=False,jobs=len(jobs)))
    array=submit(out,'worker',f'0-{len(jobs)-1}');final=submit(out,'finish',dependency=array)
    g.r.write(out/'launch.json',dict(generated_at=g.n.now(),array=array,finalizer=final))
    print(json.dumps(dict(jobs=len(jobs),array=array,finalizer=final)))


def render(out,m,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True)
    import csv
    baseline=list(csv.DictReader((BASE/'results/endpoints.csv').open()))
    allrows=rows+[dict(x,length=float(x['length']),endpoint=float(x['endpoint']) if x['endpoint'] else None,field_tail_pass=x['field_tail_pass']=='True') for x in baseline if x['kind'] in ['native','both_initial']]
    for i,gate in enumerate(g.r.DEFINITIONS):
        for j,geometry in enumerate(g.t.TESTS):
            ax=axes[i,j]
            for kind in ['native','both_initial']+KINDS:
                rs=sorted([x for x in allrows if x['gate']==gate and x['geometry']==geometry and x['kind']==kind],key=lambda x:x['length'])
                line=ax.plot([x['length'] for x in rs],[x['endpoint'] if x['endpoint'] else np.nan for x in rs],'-',label=kind)[0]
                for x in rs:
                    ax.plot(x['length'],x['endpoint'] or 5120.,'o' if x['endpoint'] else '^',color=line.get_color(),mfc=line.get_color() if x['field_tail_pass'] and x['endpoint'] else 'white')
            ax.set(title=geometry+' / '+gate,yscale='log',xlabel='Length',ylabel='Shortest passing tested cue duration');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'results/duration_comparison.png',dpi=170);fig.savefig(out/'results/duration_comparison.pdf');plt.close(fig)
    (out/'results/README.md').write_text('# Goryachev follow-up controls\n\nGenerated '+g.n.now()+'. Producer finished.\n\n'+m['intervention']+'\n\nAll start exact native checkpoints, conserve actual mass, restore native recovery, and retain separate orientation/strength gates and field flags. Matched no-cue controls are audited. Searches refine every sampled fail-to-pass boundary to 1%; narrower unsampled intervals remain possible. Open symbols retain failed field flags; triangles indicate censoring at 5120. Diffusion matching changes physical transport and the cue-only no-cue trajectory; it does not equalize pole material or supply concentration. Tight feedback compensation still needs verification against actual substrate approaching zero.\n')
    (out/'results/index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Goryachev follow-up</title></head><body><h1>Goryachev follow-up</h1><p>Generated '+g.n.now()+'; producer finished.</p><img style="width:100%;max-width:1400px" src="duration_comparison.png"><p><a href="README.md">Design</a> · <a href="endpoints.csv">Endpoints</a> · <a href="independent_audit.json">Audit</a></p></body></html>')


def self_check():
    x=(np.arange(201)+.5)*39/201;s=np.vstack([np.linspace(.001,100,201),np.linspace(.02,.4,201)]);cue=np.linspace(0,6.4,201)
    a=ExtendedControl(kind='transport_matched',length=39.,reference_v=.3,initial_v=.09)
    np.testing.assert_array_equal(a.reaction(0,x,s,cue),Goryachev().reaction(0,x,s,cue))
    assert a.diffusivities==(0.0025,10*(39/9)**2)
    for kind,basekind in [('both_initial_transport','both_initial'),('feedback_tight','feedback_clamp')]:
        a=ExtendedControl(kind=kind,length=39.,reference_v=.3,initial_v=.09)
        b=g.AvailabilityControl(kind=basekind,length=39.,reference_v=.3,initial_v=.09)
        np.testing.assert_array_equal(a.reaction(0,x,s,cue),b.reaction(0,x,s,cue))
    print('Transport coefficient and reaction equivalence checks passed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker','finish','self-check']);p.add_argument('--out',type=Path,default=OUT);p.add_argument('--slot',type=int)
    args=p.parse_args();out=args.out.resolve()
    try:
        if args.action=='self-check':self_check()
        elif args.action=='init':initialize(out)
        else:
            g.AvailabilityControl=ExtendedControl;g.render=render
            if args.action=='worker':g.worker(out,args.slot)
            else:g.finish(out)
    except Exception:
        if out.exists():g.r.write(out/'results'/('worker_'+str(args.slot)+'.error.json' if args.action=='worker' else 'error.json'),dict(generated_at=g.n.now(),error=traceback.format_exc()))
        raise
