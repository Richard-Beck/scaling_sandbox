"""Immutable runtime for late precision additions; see runtime_manifest.json."""
from common import *
import argparse,time,traceback
import numpy as np
import native_campaign as native
def worker(slot,topup=False):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native integration requires SLURM')
    m=read(RUN/'manifest.json');ci,within=divmod(slot,72 if topup else 24)
    rep=within+24 if topup else within;c=m['conditions'][ci]
    if digest(SOURCE/'manifest.json')!=m['source_sha256']:raise RuntimeError('Source manifest changed')
    folder=RUN/'workers'/f'c{ci:03d}'/f'seed{rep:03d}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'status.json').exists() and read(folder/'status.json').get('complete'):return
    started=time.monotonic();status=dict(condition=c,replicate=rep,complete=False,job_id=os.environ['SLURM_JOB_ID'],runtime_manifest_sha256=digest(Path(__file__).with_name('runtime_manifest.json')))
    try:
        seed=202610230600000+c['condition_index']*1000+rep
        t,p,C,basal,checks,audit=native.solve(SOURCE,c,seed)
        n=p.shape[1];theta=2*np.pi*(np.arange(n)+.5)/n
        # Positive Fourier convention matches existing M1 and heading exactly.
        modes=(p@np.exp(1j*np.outer(theta,np.arange(min(33,n//2+1))))/n).astype('complex64')
        P=np.abs(modes[:,1]/modes[:,0]).astype('float32')
        temp=folder/'modes.tmp.npz'
        np.savez_compressed(temp,time=t.astype('float64'),circumference=C.astype('float32'),modes=modes,P_abs=P)
        temp.replace(folder/'modes.npz')
        audit['mode_only_output']=True;write(folder/'noise_audit.json',audit)
        status.update(complete=True,seed=seed,validation=checks,sha256=digest(folder/'modes.npz'),
                      bytes=(folder/'modes.npz').stat().st_size)
    except Exception as e:status.update(error=str(e),traceback=traceback.format_exc());raise
    finally:status['wall_seconds']=time.monotonic()-started;write(folder/'status.json',status)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--slot',type=int,required=True);p.add_argument('--topup',action='store_true');a=p.parse_args();worker(a.slot,a.topup)
