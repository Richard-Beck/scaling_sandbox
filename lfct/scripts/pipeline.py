#!/usr/bin/env python3
"""One field per invocation: resumable compact outputs, owned scratch cleanup after validation."""
import argparse,os,resource,shutil,tempfile,time,traceback
from pathlib import Path
from common import BASE,config,field_row,profile_fields,frame_selection,output_dir,digest,code_hash,sha,read_json,write_json,csv_read
from prepare import prepare


def validate(out,prepared):
    import numpy as np,tifffile
    from analyze import load_spots
    frames=prepared['frames'];source=load_spots(prepared['annotation'],set(frames));expected={(f,s['id']):s for f,spots in source.items() for s in spots}
    rows=csv_read(out/'detections.csv.gz');assert len(rows)==len(expected)
    lookup={(int(r['frame']),int(r['spot_id'])):r for r in rows};assert set(lookup)==set(expected)
    for key,s in expected.items():
        r=lookup[key];assert float(r['x'])==s['x'] and float(r['y'])==s['y'];assert r['track_id']==str(s['track']) if s['track'] is not None else not r['track_id'];assert r['track_name']==(s['name'] or '')
    for kind,col in [('cells','cell_label'),('nuclei','nucleus_label')]:
        with tifffile.TiffFile(out/(kind+'.tif')) as stack:
            assert len(stack.pages)==len(frames)
            for f,page in zip(frames,stack.pages):
                assert int(page.description)==f;mask=page.asarray();assert mask.dtype==np.uint16 and mask.shape==(prepared['height'],prepared['width'])
                assigned=[int(r[col]) for r in rows if int(r['frame'])==f and r[col]];assert len(assigned)==len(set(assigned));assert set(assigned)<=set(np.unique(mask))
    return dict(frames=len(frames),source_detections_preserved=len(rows),duplicate_assignments=0,status='passed')


def main():
    p=argparse.ArgumentParser();p.add_argument('--profile',choices=['smoke','pilot','full'],default='smoke');p.add_argument('--field');p.add_argument('--run-name');p.add_argument('--index',type=int);p.add_argument('--scratch-root',default=os.environ.get('SLURM_TMPDIR','/tmp'));a=p.parse_args()
    fields=profile_fields(a.profile);key=a.field if a.field else fields[a.index or 0];assert key in fields
    run_name=a.run_name or a.profile
    assert run_name and all(c.isalnum() or c in '_-' for c in run_name),'Use a simple run name'
    row=field_row(key);frames=frame_selection(row,a.profile);out=output_dir(run_name,key);out.mkdir(parents=True,exist_ok=True)
    cfg=config();fingerprint=digest(dict(source=row,frames=frames,segmentation=cfg['segmentation'],rfp=cfg['rfp_normalization'],model_sha256=sha(cfg['model']),**({'source_mapping':cfg['source_mapping_overrides'][key]} if key in cfg.get('source_mapping_overrides',{}) else {})))
    full_fingerprint=digest(dict(segmentation=fingerprint,config=cfg,code=code_hash()))
    completed=out/'completed.json'
    if completed.exists():
        prior=read_json(completed)
        if prior['full_fingerprint']==full_fingerprint:
            for name,h in prior['output_hashes'].items():assert sha(out/name)==h,'Completed output checksum mismatch'
            for name,h in prior['qc_asset_hashes'].items():assert sha(BASE/'reports'/run_name/name)==h,'QC PNG missing or changed'
            print(f'{key}: validated resume, no download or inference',flush=True);return
    lock=out/'.running';lock.mkdir() # fail rather than race two writers for one field
    scratch=Path(tempfile.mkdtemp(prefix='lfct-',dir=a.scratch_root));write_json(scratch/'owner.json',dict(key=key,pid=os.getpid(),output=str(out)))
    start=time.monotonic()
    try:
        prepared=prepare(row,frames,scratch)
        from segment import run as segment
        info=segment(prepared,out,fingerprint)
        from analyze import run as analyze
        summary=analyze(prepared,out,a.profile,run_name);summary['total_movie_frames']=int(row['nframes']);write_json(out/'analysis.json',summary)
        verified=validate(out,prepared);write_json(out/'validation.json',verified)
        qc=read_json(out/'qc.json');assets={image:sha(BASE/'reports'/run_name/image) for item in qc for f in item['frames'] for image in f['images'].values()}
        outputs=['cells.tif','nuclei.tif','segmentation.json','detections.csv.gz','masks.csv.gz','runs.csv.gz','frame_qc.csv','analysis.json','qc.json','validation.json']
        provenance=dict(key=key,profile=a.profile,run_name=run_name,complete_field=len(frames)==int(row['nframes']),frames=frames,full_fingerprint=full_fingerprint,config=cfg,source=row,downloads=prepared['sources'],job_id=os.environ.get('SLURM_JOB_ID'),prepare_seconds=prepared['seconds'],total_seconds=time.monotonic()-start,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,output_hashes={name:sha(out/name) for name in outputs},qc_asset_hashes=assets,scratch_cleanup='validated owned scratch removed; pre-existing inputs retained')
        write_json(completed,provenance)
        # Only delete this invocation's mkdtemp directory; never external reused images.
        assert read_json(scratch/'owner.json')['output']==str(out) and scratch.name.startswith('lfct-')
        shutil.rmtree(scratch)
        print(f'COMPLETED {key}: {summary}; seconds={provenance["total_seconds"]:.1f}',flush=True)
    except Exception:
        write_json(out/'failure.json',dict(key=key,scratch=str(scratch),error=traceback.format_exc()))
        raise
    finally:lock.rmdir()
if __name__=='__main__':main()
