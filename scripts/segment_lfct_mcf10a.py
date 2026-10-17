#!/usr/bin/env python3
"""GPU-only CPSAM inference, independent phase-cell and RFP-nucleus masks."""
import argparse,hashlib,json,os,socket,time,importlib.metadata
from pathlib import Path
import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter
import torch
from cellpose import models
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/lfct_mcf10a_fld102'
MODEL=Path('/home/4473331/.cellpose/models/cpsam')
PARAMS=dict(diameter=None,normalize=True,flow_threshold=0.4,cellprob_threshold=0.0,min_size=15,max_size_fraction=0.4,do_3D=False,augment=False,batch_size=8,bsize=256,tile_overlap=0.1,resample=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('--frames',type=str);p.add_argument('--rfp-tile',type=int,default=0);p.add_argument('--rfp-background',action='store_true');p.add_argument('--output-subdir',default='cpsam');p.add_argument('--channels',default='phase,rfp');p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=1);a=p.parse_args()
    assert torch.cuda.is_available(),'CUDA unavailable; refusing slow CPU fallback'
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
    torch.manual_seed(0);np.random.seed(0)
    manifest=json.loads((DATA/'source_manifest.json').read_text())
    indices=[int(x) for x in a.frames.split(',')] if a.frames else list(range(a.shard,len(manifest['frames']),a.shards))
    output=DATA/a.output_subdir;output.mkdir(exist_ok=True)
    model=models.CellposeModel(gpu=True,pretrained_model=str(MODEL))
    model_hash=hashlib.sha256(MODEL.read_bytes()).hexdigest()
    info=dict(cellpose=importlib.metadata.version('cellpose'),torch=torch.__version__,model=str(MODEL),model_sha256=model_hash,parameters=PARAMS,seed=0,host=socket.gethostname(),gpu=torch.cuda.get_device_name(0),slurm_job_id=os.environ.get('SLURM_JOB_ID'))
    print(json.dumps(info),flush=True)
    for frame in indices:
        row=manifest['frames'][frame]
        for channel in a.channels.split(','):
            params=dict(PARAMS)
            if channel=='rfp' and a.rfp_tile:
                params['normalize']={'normalize':True,'tile_norm_blocksize':a.rfp_tile,'percentile':[1,99]}
            folder=output/channel;folder.mkdir(exist_ok=True)
            dest=folder/f'{frame:04d}.tif';metadata=folder/f'{frame:04d}.json'
            if dest.exists() and metadata.exists():
                previous=json.loads(metadata.read_text())
                assert previous['parameters']==params and previous['model_sha256']==model_hash and previous.get('rfp_background',False)==(channel=='rfp' and a.rfp_background)
                continue
            before=time.monotonic();image=tifffile.imread(ROOT/row[channel]['path'])
            if channel=='rfp' and a.rfp_background:
                raw=image.astype(np.float32)
                image=np.maximum(gaussian_filter(raw,1.0)-gaussian_filter(raw,50.0),0)
            masks,flows,styles=model.eval(image,**params)
            assert masks.shape==image.shape and masks.max()<65536
            tifffile.imwrite(dest,masks.astype(np.uint16),compression='zlib',metadata={'axes':'YX'})
            stats={**info,'parameters':params,'rfp_background':channel=='rfp' and a.rfp_background,'frame':frame,'channel':channel,'input':row[channel]['path'],'labels':int(masks.max()),'foreground_pixels':int(np.count_nonzero(masks)),'seconds':time.monotonic()-before}
            metadata.write_text(json.dumps(stats,indent=2))
            print(json.dumps({k:stats[k] for k in ['frame','channel','labels','foreground_pixels','seconds']}),flush=True)
            del masks,flows,styles
if __name__=='__main__':main()
