"""CPSAM inference to two lossless multipage TIFF stacks with read-back validation."""
import importlib.metadata,os,socket,time
from pathlib import Path
import numpy as np
import tifffile
from common import config,sha,write_json,read_json,digest

def run(prepared,out,fingerprint):
    cfg=config();frames=prepared['frames'];out=Path(out);out.mkdir(parents=True,exist_ok=True)
    marker=out/'segmentation.json'
    if marker.exists():
        previous=read_json(marker);assert previous['fingerprint']==fingerprint,'Changed inputs/configuration: choose a new run profile/output'
        for c in ['cells','nuclei']:assert sha(out/(c+'.tif'))==previous['files'][c]['sha256']
        return previous
    import torch
    from cellpose import models
    assert torch.cuda.is_available(),'GPU unavailable; refusing CPU fallback'
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK',4)));torch.manual_seed(0);np.random.seed(0)
    start=time.monotonic();model=models.CellposeModel(gpu=True,pretrained_model=cfg['model']);load_seconds=time.monotonic()-start
    info=dict(fingerprint=fingerprint,frames=frames,cellpose=importlib.metadata.version('cellpose'),torch=torch.__version__,model_sha256=sha(cfg['model']),gpu=torch.cuda.get_device_name(0),host=socket.gethostname(),job_id=os.environ.get('SLURM_JOB_ID'),model_load_seconds=load_seconds,parameters={},timings=[],files={})
    before=time.monotonic()
    for channel,name in [('phase','cells'),('rfp','nuclei')]:
        params=dict(cfg['segmentation'])
        if channel=='rfp':params['normalize']=cfg['rfp_normalization']
        info['parameters'][channel]=params
        dest=out/(name+'.tif');tmp=out/(name+'.tif.part');pixel_hashes=[]
        with tifffile.TiffWriter(tmp,bigtiff=False) as writer:
            for frame in frames:
                image=tifffile.imread(prepared['images'][channel][str(frame)])
                then=time.monotonic();masks,flows,styles=model.eval(image,**params);elapsed=time.monotonic()-then
                assert masks.shape==image.shape and masks.max()<65536
                masks=masks.astype(np.uint16)
                import hashlib
                pixel_hashes.append(hashlib.sha256(masks.tobytes()).hexdigest())
                writer.write(masks,compression='deflate',metadata=None,description=str(frame),photometric='minisblack')
                info['timings'].append(dict(frame=frame,channel=channel,seconds=elapsed,labels=int(masks.max())))
                print(f'{prepared["key"]} {channel} frame={frame} masks={masks.max()} seconds={elapsed:.3f}',flush=True)
                del masks,flows,styles
        with tifffile.TiffFile(tmp) as check:
            assert len(check.pages)==len(frames)
            for frame,page,expected in zip(frames,check.pages,pixel_hashes):
                assert int(page.description)==frame and page.dtype==np.uint16
                assert hashlib.sha256(page.asarray().tobytes()).hexdigest()==expected
        tmp.replace(dest);info['files'][name]=dict(sha256=sha(dest),bytes=dest.stat().st_size,pixel_sha256=pixel_hashes)
    info['inference_and_writing_seconds']=time.monotonic()-before
    info['max_gpu_memory_bytes']=torch.cuda.max_memory_allocated()
    write_json(marker,info);return info
