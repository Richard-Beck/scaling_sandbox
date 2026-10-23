"""Documented numerical retry for one stalled constant-mean preparation.
Native equations, inputs, pulse amplitudes and measured conditions are identical.
Only the long unobserved preparation uses original campaign ODE tolerances and a
bounded max_step; all response measurements retain the tighter declared settings.
"""
from pathlib import Path
import shutil,sys,os
import figure3_round5_edge_20261023 as edge
BASE=edge.RUN
RETRY=BASE.with_name('figure3_round5_edge_numeric_retry_20261023')
def prepare():
    RETRY.mkdir()
    for name in ('inputs','frozen_sources'):shutil.copytree(BASE/name,RETRY/name)
    for name in ('workers','logs','results'):(RETRY/name).mkdir()
    m=edge.read(BASE/'manifest.json');m['numerical_retry']=dict(slot=71,reason='Single long mean-input integration stalled under stricter than native tolerances; original output retained. No model or physical input change.',preparation_solver=dict(rtol=1e-8,atol=1e-10,max_step=5.),response_solver=dict(rtol=2e-10,atol=2e-12),original_run=str(BASE))
    edge.write(RETRY/'manifest.json',m);edge.write(RETRY/'freeze.json',dict(manifest_sha256=edge.sha(RETRY/'manifest.json')))
    shutil.copy2(Path(__file__),RETRY/'frozen_sources'/Path(__file__).name)
def worker():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('SLURM required')
    original=edge.solve_ivp
    def solve(fun,t_span,y0,*args,**kwargs):
        if t_span[1]==1600.:
            kwargs.update(rtol=1e-8,atol=1e-10,max_step=5.)
        return original(fun,t_span,y0,*args,**kwargs)
    edge.solve_ivp=solve;edge.RUN=RETRY
    edge.worker(71)
if __name__=='__main__':
    prepare() if sys.argv[1]=='prepare' else worker()
