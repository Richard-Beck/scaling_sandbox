"""Ordered multiprocessing preserves all tracked scientific output bytes."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import figure4_round3_analysis_20261022 as A
import figure4_round3_estimators_20261022 as E


class ParallelTracked(unittest.TestCase):
    def test_serial_and_parallel_tracked_exports_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);groups={};contexts={};trackgroups={};t=np.arange(201)*5.
            for condition,size in enumerate([20.,40.]):
                key=('wave_pinning',size,.03);moments=[];records=[]
                for seed in range(2):
                    theta=np.cumsum(np.random.default_rng(81+condition*2+seed).normal(0,.03,len(t)));theta[100:]+=.7
                    P=(.8*np.exp(1j*theta)).astype('complex64');ids=(np.arange(len(t))>=100).astype('int32');path=root/f'c{condition}_s{seed}.npz'
                    np.savez_compressed(path,time=t.astype('float32'),circumference=np.full(len(t),size,dtype='float32'),P_real=P.real,P_imag=P.imag,
                                        region_offset=np.arange(len(t)+1),region_id=ids,region_heading=np.mod(theta,2*np.pi).astype('float32'),dominant_track_id=ids)
                    target,lags=E.lag_grid(key[0],t[-1]-t[0],5.);moments.append(E.seed_moments(P.astype(complex),abs(P)>=.1,5.,lags));records.append(dict(compact=str(path)))
                groups[key]=E.summarize(moments,[size]*2,('compact',key),target,lags);contexts[key]=dict(source='parallel_equivalence_fixture');trackgroups[key]=records
            with patch.dict(os.environ,{'SLURM_JOB_ID':'','SLURM_JOBID':'','SLURM_CPUS_PER_TASK':'4'}):self.assertEqual(A._compact_tracked_worker_limit(4),1)
            with patch.dict(os.environ,{'SLURM_JOB_ID':'fixture','SLURM_CPUS_PER_TASK':'4'}):
                serial=A.compact_tracked(root/'serial',groups,contexts,trackgroups,workers=1)
                parallel=A.compact_tracked(root/'parallel',groups,contexts,trackgroups,workers=4)
            self.assertEqual(serial['workers'],1);self.assertEqual(parallel['workers'],2)
            self.assertEqual(parallel['BLAS_threads'],1)
            for name in ['drift_decomposition.csv','within_tenure_msd.csv','detector_free_mle.csv','seed0_switch_events.csv']:
                self.assertEqual((root/'serial/supporting/data'/name).read_bytes(),(root/'parallel/supporting/data'/name).read_bytes(),name)


if __name__=='__main__':unittest.main()
