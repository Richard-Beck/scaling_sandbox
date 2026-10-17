"""Assignments, exact pixel-square clearances, uninterrupted runs, and bounded QC assets."""
import gzip,math,time,xml.etree.ElementTree as ET
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np,tifffile
from PIL import Image
from common import BASE,config,read_json,write_json,csv_write
from geometry import mask_objects,assign_masks,add_nuclear_overlap_matches,clearances,point_square_gap,pixel_square_gap,bounds


def load_spots(annotation,selected):
    with gzip.open(annotation,'rb') as f:doc=ET.parse(f)
    meta=doc.find('Settings/ImageData').attrib
    assert float(meta['pixelwidth'])==float(meta['pixelheight'])==1,'Physical-unit XML requires explicit conversion'
    names={};membership={};inc=defaultdict(list);out=defaultdict(list)
    for t in doc.findall('Model/AllTracks/Track'):
        tid=int(t.get('TRACK_ID'));names[tid]=t.get('name')
        for e in t.findall('Edge'):
            a,b=int(e.get('SPOT_SOURCE_ID')),int(e.get('SPOT_TARGET_ID'))
            for sid in (a,b):
                assert sid not in membership or membership[sid]==tid
                membership[sid]=tid
            out[a].append(b);inc[b].append(a)
    allspots={};byframe=defaultdict(list)
    for node in doc.findall('Model/AllSpots/SpotsInFrame/Spot'):
        sid=int(node.get('ID'));f=int(node.get('FRAME'));x=float(node.get('POSITION_X'));y=float(node.get('POSITION_Y'))
        original=[]
        if f in selected:
            vals=np.fromstring(node.text or '',sep=' ')
            if len(vals)>=6 and len(vals)%2==0:original=np.round(vals.reshape(-1,2)+[x,y],2).tolist()
        s=dict(id=sid,frame=f,x=x,y=y,track=membership.get(sid),name=names.get(membership.get(sid)),original=original)
        allspots[sid]=s;byframe[f].append(s)
    for f in sorted(byframe):
        for s in byframe[f]:
            parents=inc[s['id']]
            if len(parents)==1 and len(out[parents[0]])==1 and allspots[parents[0]]['frame']==f-1:s['branch']=allspots[parents[0]]['branch']
            else:s['branch']=s['id']
    return {f:byframe[f] for f in sorted(selected)}


def eligible_runs(rows,min_frames):
    by=defaultdict(list)
    for r in rows:by[r['branch_start_spot_id']].append(r)
    results=[];longest=0
    def finish(chain):
        nonlocal longest
        longest=max(longest,len(chain))
        if len(chain)<min_frames:return
        finite=[r['clearance_px'] for r in chain if r['clearance_px'] is not None]
        results.append(dict(run_id=f'{chain[0]["branch_start_spot_id"]}:{chain[0]["frame"]}',branch_start_spot_id=chain[0]['branch_start_spot_id'],track_id=chain[0]['track_id'],track_name=chain[0]['track_name'],first_frame=chain[0]['frame'],last_frame=chain[-1]['frame'],n_frames=len(chain),first_spot_id=chain[0]['spot_id'],last_spot_id=chain[-1]['spot_id'],min_clearance_px=min(finite) if finite else None,median_clearance_px=float(np.median(finite)) if finite else None))
    for branch,items in by.items():
        chain=[]
        for r in sorted(items,key=lambda r:r['frame']):
            if not r['eligible'] or (chain and r['frame']!=chain[-1]['frame']+1):finish(chain);chain=[]
            if r['eligible']:chain.append(r)
        finish(chain)
    return results,longest


def frame_analysis(cells,nuclei,spots):
    cfg=config();iso=cfg['isolation'];match=cfg['matching']
    objs={c:mask_objects(a) for c,a in [('cells',cells),('nuclei',nuclei)]}
    cm=assign_masks(cells,objs['cells'],spots,match['max_gap_px'],match['margin_px'])
    nm=assign_masks(nuclei,objs['nuclei'],spots,match['max_gap_px'],match['margin_px'])
    add_nuclear_overlap_matches(cells,nuclei,objs['nuclei'],cm,nm)
    gap=clearances(objs['cells']);bycell={o['label']:o for o in objs['cells']}
    orphan=[]
    for o in objs['nuclei']:
        labs,counts=np.unique(cells[nuclei==o['label']],return_counts=True)
        if max([int(n) for l,n in zip(labs,counts) if l] or [0])/o['area']<.5:orphan.append(o)
    missing=[(i,s) for i,s in enumerate(spots) if cm[i]['label'] is None]
    rows=[]
    for i,s in enumerate(spots):
        c,n=cm[i],nm[i];reasons=[];distance=None;nearest=None;border=None;unknown=[]
        if c['status']!='contained':reasons.append('cell_'+c['status'])
        if c['label'] is not None:
            o=bycell[c['label']];distance,nearest=gap[c['label']];box=bounds(o)
            border=float(min(box[0]+.5,box[1]+.5,cells.shape[1]-.5-box[2],cells.shape[0]-.5-box[3]))
            if distance is not None and distance<iso['clearance_px']:reasons.append('neighbor_too_close')
            if border<iso['clearance_px']:reasons.append('field_edge')
            for j,other in missing:
                if j!=i and point_square_gap(o,[other['x'],other['y']])<=iso['unknown_neighbor_guard_px']:unknown.append('spot:'+str(other['id']))
            for nucleus in orphan:
                if pixel_square_gap(o,nucleus)<=iso['unknown_neighbor_guard_px']:unknown.append('nucleus:'+str(nucleus['label']))
            if unknown:reasons.append('uncertain_neighbor')
        if n['cellOverlap'] is not None and n['cellOverlap']<.5:reasons.append('low_nucleus_cell_overlap')
        rows.append(dict(frame=s['frame'],spot_id=s['id'],track_id=s['track'],track_name=s['name'],branch_start_spot_id=s['branch'],x=s['x'],y=s['y'],original_outline=bool(s['original']),cell_label=c['label'],cell_status=c['status'],cell_gap_px=c['gap'],nucleus_label=n['label'],nucleus_status=n['status'],nucleus_cell_overlap=n['cellOverlap'],clearance_px=distance,nearest_cell_label=nearest,edge_clearance_px=border,eligible=not reasons,exclusions=';'.join(reasons),uncertain_neighbors=';'.join(unknown)))
    public={}
    for key,matches in [('cells',cm),('nuclei',nm)]:
        assigned={m['label']:i for i,m in enumerate(matches) if m['label'] is not None}
        public[key]=[{**{k:v for k,v in o.items() if k!='tree'},'spot_index':assigned.get(o['label'])} for o in objs[key]]
    return rows,public


def png(source,dest,crop):
    a=tifffile.imread(source);lo,hi=np.percentile(a,[1,99.5]);x,y,w,h=crop
    a=(np.clip((a[y:y+h,x:x+w].astype(np.float32)-lo)/max(1,hi-lo),0,1)*255).astype('uint8')
    dest.parent.mkdir(parents=True,exist_ok=True);Image.fromarray(a).save(dest,optimize=True)
    return dict(low=float(lo),high=float(hi),bytes=dest.stat().st_size)


def choose_clips(rows,runs,maxclips,length):
    clips=[]
    for run in sorted(runs,key=lambda r:(-r['n_frames'],r['run_id']))[:maxclips]:
        clips.append(dict(kind='isolated',branch=run['branch_start_spot_id'],frames=list(range(run['first_frame'],min(run['last_frame']+1,run['first_frame']+length))),run_id=run['run_id']))
    if not clips and maxclips:
        # A real contiguous negative control exercises targeted QC even if no positive run exists.
        by=defaultdict(list)
        for r in rows:by[r['branch_start_spot_id']].append(r)
        for branch,items in sorted(by.items(),key=lambda kv:-len(kv[1])):
            fs=sorted(r['frame'] for r in items)
            for j in range(len(fs)-length+1):
                part=fs[j:j+length]
                if part==list(range(part[0],part[0]+length)):
                    return [dict(kind='control_not_isolated',branch=branch,frames=part,run_id=None)]
    return clips


def run(prepared,out,profile,run_name=None):
    start=time.monotonic();cfg=config();out=Path(out);frames=prepared['frames'];spots=load_spots(prepared['annotation'],set(frames))
    allrows=[];inventory=[];fq=[];geometry={}
    with tifffile.TiffFile(out/'cells.tif') as cstack,tifffile.TiffFile(out/'nuclei.tif') as nstack:
        assert len(cstack.pages)==len(nstack.pages)==len(frames)
        for page,f in enumerate(frames):
            cells,nuclei=cstack.pages[page].asarray(),nstack.pages[page].asarray()
            rows,public=frame_analysis(cells,nuclei,spots[f]);allrows.extend(rows);geometry[f]=public
            for kind in ['cells','nuclei']:
                for o in public[kind]:
                    idx=o['spot_index'];inventory.append(dict(frame=f,kind=kind,label=o['label'],area_px=o['area'],x=o['x'],y=o['y'],edge=o['edge'],spot_id=spots[f][idx]['id'] if idx is not None else None))
            fq.append(dict(frame=f,detections=len(rows),cells=len(public['cells']),nuclei=len(public['nuclei']),assigned_cells=sum(r['cell_label'] is not None for r in rows),assigned_nuclei=sum(r['nucleus_label'] is not None for r in rows),eligible=sum(r['eligible'] for r in rows)))
    runs,longest=eligible_runs(allrows,cfg['isolation']['min_frames'])
    csv_write(out/'detections.csv.gz',allrows);csv_write(out/'masks.csv.gz',inventory)
    csv_write(out/'runs.csv.gz',runs,['run_id','branch_start_spot_id','track_id','track_name','first_frame','last_frame','n_frames','first_spot_id','last_spot_id','min_clearance_px','median_clearance_px'])
    csv_write(out/'frame_qc.csv',fq)
    # Deterministic QC selection independent of segmentation success.
    n=cfg['profiles'][profile].get('overview_per_field',cfg['qc']['overview_per_field'])
    if profile=='full':
        import hashlib
        seed=cfg['qc']['seed']+int(hashlib.sha256(prepared['key'].encode()).hexdigest()[:8],16)
        rng=np.random.default_rng(seed);thirds=np.array_split(np.array(frames),3);chosen=rng.choice(3,size=min(n,3),replace=False)
        overview=sorted(int(rng.choice(thirds[i])) for i in chosen)
    else:overview=sorted({frames[int(i)] for i in np.linspace(0,len(frames)-1,min(n,len(frames)))})
    items=[dict(kind='overview',frames=overview,crop=[0,0,prepared['width'],prepared['height']],branch=None)]
    clips=choose_clips(allrows,runs,cfg['qc']['target_clips_per_field'],cfg['qc']['clip_frames'])
    for clip in clips:
        points=[]
        for f in clip['frames']:
            rr=[r for r in allrows if r['frame']==f and r['branch_start_spot_id']==clip['branch']]
            if not rr:continue
            r=rr[0];points.append([r['x'],r['y']])
            for o in geometry[f]['cells']:
                if o['label'] in (r['cell_label'],r['nearest_cell_label']):
                    for loop in o['paths']:points.extend(loop)
        pts=np.array(points);pad=cfg['qc']['crop_margin_px'];lo=np.maximum(0,np.floor(pts.min(0)-pad)).astype(int);hi=np.minimum([prepared['width'],prepared['height']],np.ceil(pts.max(0)+pad)).astype(int)
        clip['crop']=[int(lo[0]),int(lo[1]),int(hi[0]-lo[0]),int(hi[1]-lo[1])];items.append(clip)
    qc=[];assetdir=BASE/'reports'/(run_name or profile)/'assets'/prepared['key']
    for number,item in enumerate(items):
        entry={**item,'key':prepared['key'],'frames':[]}
        for f in item['frames']:
            frame_data=dict(frame=f,crop=item['crop'],spots=spots[f],rows=[r for r in allrows if r['frame']==f],**geometry[f],images={},display={})
            for channel in ['phase','rfp']:
                dst=assetdir/f'{number}_{channel}_{f:04d}.png'
                frame_data['display'][channel]=png(prepared['images'][channel][str(f)],dst,item['crop'])
                frame_data['images'][channel]=str(dst.relative_to(BASE/'reports'/(run_name or profile)))
            entry['frames'].append(frame_data)
        qc.append(entry)
    write_json(out/'qc.json',qc)
    summary=dict(key=prepared['key'],profile=profile,processed_frames=len(frames),total_movie_frames=None,detections=len(allrows),assigned_cells=sum(r['cell_label'] is not None for r in allrows),assigned_nuclei=sum(r['nucleus_label'] is not None for r in allrows),eligible_observations=sum(r['eligible'] for r in allrows),isolated_runs=len(runs),longest_eligible_run=longest,cell_masks=len([m for m in inventory if m['kind']=='cells']),nuclear_masks=len([m for m in inventory if m['kind']=='nuclei']),exclusions=dict(Counter(reason for r in allrows for reason in r['exclusions'].split(';') if reason)),analysis_qc_seconds=time.monotonic()-start)
    write_json(out/'analysis.json',summary)
    return summary
