#!/usr/bin/env python3
"""Match CPSAM masks to immutable TrackMate spots and build an audited full-stack viewer."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,base64,csv,io,json,math,xml.etree.ElementTree as ET
import cv2
import numpy as np
import tifffile
from scipy.spatial import cKDTree
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/lfct_mcf10a_fld102'
OUT=ROOT/'reports/lfct_mcf10a_fld102_cpsam_assets'
PALETTE=['#56B4E9','#E69F00','#009E73','#CC79A7','#F0E442','#D55E00','#8172B3','#80CDC1']
MAX_GAP=15.0
MARGIN=3.0


def read_tracks():
    doc=ET.parse(ROOT/'data/lfct_trackmate_20x/mcf10a/FLD_102.xml')
    names={};membership={};incoming=defaultdict(list);outgoing=defaultdict(list)
    for t in doc.findall('Model/AllTracks/Track'):
        tid=int(t.get('TRACK_ID'));names[tid]=t.get('name')
        for e in t.findall('Edge'):
            a,b=int(e.get('SPOT_SOURCE_ID')),int(e.get('SPOT_TARGET_ID'))
            for sid in (a,b):
                assert sid not in membership or membership[sid]==tid
                membership[sid]=tid
            outgoing[a].append(b);incoming[b].append(a)
    frames=[[] for _ in range(385)];all_spots={}
    for el in doc.findall('Model/AllSpots/SpotsInFrame/Spot'):
        sid=int(el.get('ID'));f=int(el.get('FRAME'));x=float(el.get('POSITION_X'));y=float(el.get('POSITION_Y'))
        vals=np.fromstring(el.text or '',sep=' ')
        old=[]
        if len(vals)>=6 and len(vals)%2==0:
            old=np.round(vals.reshape(-1,2)+[x,y],2).tolist()
        spot=dict(id=sid,frame=f,x=x,y=y,track=membership.get(sid),name=names.get(membership.get(sid)),original=old)
        frames[f].append(spot);all_spots[sid]=spot
    for f in frames:
        for spot in sorted(f,key=lambda x:x['id']):
            sid=spot['id'];parents=incoming[sid]
            if len(parents)==1 and len(outgoing[parents[0]])==1:
                parent=all_spots[parents[0]]
                assert parent['frame']<spot['frame']
                spot['branch']=parent['branch']
            else:spot['branch']=sid
    return frames,names


def mask_objects(labels):
    objects=[]
    # Exact integer pixel masks are retained on disk; viewer paths simplify <=0.65 px.
    for lab in np.unique(labels):
        if lab==0:continue
        yy,xx=np.where(labels==lab)
        xmin,xmax,ymin,ymax=int(xx.min()),int(xx.max()),int(yy.min()),int(yy.max())
        binary=(labels[ymin:ymax+1,xmin:xmax+1]==lab).astype('uint8')
        contours,_=cv2.findContours(binary,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
        paths=[];boundaries=[]
        for contour in contours:
            c=contour[:,0,:]+[xmin,ymin]
            if len(c)>0:boundaries.extend(c.tolist())
            ap=cv2.approxPolyDP(contour,.65,True)[:,0,:]+[xmin,ymin]
            if len(ap)>=3:paths.append(ap.tolist())
        # Dense boundary pixels for exact Euclidean distances, not simplified vertices.
        eroded=cv2.erode(binary,np.ones((3,3),np.uint8),borderType=cv2.BORDER_CONSTANT,borderValue=0)
        by,bx=np.where(binary-eroded)
        tree=cKDTree(np.column_stack([bx+xmin,by+ymin]))
        objects.append(dict(label=int(lab),area=len(xx),x=float(xx.mean()),y=float(yy.mean()),edge=xmin==0 or ymin==0 or xmax==labels.shape[1]-1 or ymax==labels.shape[0]-1,paths=paths,tree=tree))
    return objects


def assign_masks(labels,objects,spots):
    bylab={o['label']:o for o in objects};owners=defaultdict(list);answers=[]
    h,w=labels.shape
    for i,s in enumerate(spots):
        x,y=s['x'],s['y'];inside=0
        valid=0<=x<w and 0<=y<h
        if valid:inside=int(labels[min(h-1,int(math.floor(y+.5))),min(w-1,int(math.floor(x+.5)))])
        answer=dict(label=None,status='no_candidate',gap=None,containing=inside)
        if not valid:answer['status']='out_of_bounds'
        elif inside:owners[inside].append(i)
        answers.append(answer)
    occupied=set()
    for label,indices in owners.items():
        occupied.add(label)
        for i in indices:
            if len(indices)==1:answers[i].update(label=label,status='contained',gap=0.0)
            else:answers[i]['status']='shared_mask'
    proposals=defaultdict(list)
    for i,s in enumerate(spots):
        if answers[i]['status']!='no_candidate':continue
        candidates=sorted((float(o['tree'].query([s['x'],s['y']])[0]),o['label']) for o in objects)
        if not candidates or candidates[0][0]>MAX_GAP:continue
        gap,label=candidates[0];answers[i]['gap']=round(gap,3)
        if label in occupied:answers[i]['status']='nearest_mask_occupied';continue
        if len(candidates)>1 and candidates[1][0]-gap<MARGIN:
            answers[i]['status']='ambiguous_nearby';continue
        proposals[label].append(i)
    for label,indices in proposals.items():
        for i in indices:
            if len(indices)==1:answers[i].update(label=label,status='nearby')
            else:answers[i]['status']='competing_centroids'
    assigned=[a['label'] for a in answers if a['label'] is not None]
    assert len(assigned)==len(set(assigned))
    assert all(l in bylab for l in assigned)
    return answers


def add_nuclear_overlap_matches(cell_labels,nuc_labels,nuc_objects,cell_matches,nuc_matches):
    # A nucleus may be displaced from the tracked whole-cell centroid. Permit a
    # fallback only if exactly one nucleus has >=50% of its pixels in that cell,
    # is unassigned, and is not directly claimed by multiple tracked centroids.
    overlap={}
    for obj in nuc_objects:
        labs,counts=np.unique(cell_labels[nuc_labels==obj['label']],return_counts=True)
        overlap[obj['label']]={int(l):int(n)/obj['area'] for l,n in zip(labs,counts) if l}
    used={m['label'] for m in nuc_matches if m['label'] is not None}
    blocked={m['containing'] for m in nuc_matches if m['status']=='shared_mask'}
    for i,(cm,nm) in enumerate(zip(cell_matches,nuc_matches)):
        if nm['label'] is not None or cm['label'] is None or nm['status']=='shared_mask':continue
        candidates=[label for label,ratios in overlap.items() if ratios.get(cm['label'],0)>=.5]
        if len(candidates)==1 and candidates[0] not in used and candidates[0] not in blocked:
            label=candidates[0];nm.update(label=label,status='cell_overlap',gap=None);used.add(label)
    for cm,nm in zip(cell_matches,nuc_matches):
        nm['cellOverlap']=round(overlap.get(nm['label'],{}).get(cm['label'],0),4) if cm['label'] and nm['label'] else None
    assigned=[m['label'] for m in nuc_matches if m['label'] is not None]
    assert len(assigned)==len(set(assigned))


def color_graph(frames,key):
    graph=defaultdict(Counter)
    for spots in frames:
        for i,s in enumerate(spots):
            ident=s[key]
            if ident is None:continue
            graph[ident]
            neighbors=sorted((j for j in range(len(spots)) if j!=i),key=lambda j:(s['x']-spots[j]['x'])**2+(s['y']-spots[j]['y'])**2)[:4]
            for j in neighbors:
                other=spots[j][key]
                if other is not None and other!=ident:graph[ident][other]+=1;graph[other][ident]+=1
    colors={};usage=Counter()
    while len(colors)<len(graph):
        node=max((n for n in graph if n not in colors),key=lambda n:(len({colors[v] for v in graph[n] if v in colors}),sum(graph[n].values()),-n))
        color=min(range(len(PALETTE)),key=lambda c:(sum(w for v,w in graph[node].items() if colors.get(v)==c),usage[c],c))
        colors[node]=color;usage[color]+=1
    return {str(k):PALETTE[v] for k,v in colors.items()}


def dump_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def background(path,dest):
    image=tifffile.imread(path)
    lo,hi=np.percentile(image,[1,99.5])
    display=(np.clip((image.astype(np.float32)-lo)/max(hi-lo,1),0,1)*255).astype('uint8')
    # Lossy JPEG is solely a viewing background. Inference uses original uint16 TIFFs.
    Image.fromarray(display).save(dest,quality=55,optimize=True)
    return dict(low=float(lo),high=float(hi),saturationFraction=float(np.mean(image>=65504)),std=float(image.std()))


def main():
    p=argparse.ArgumentParser();p.add_argument('--frames',default=None);p.add_argument('--rfp-dir',default='cpsam_rfp_tile128');p.add_argument('--pilot',action='store_true');a=p.parse_args()
    OUT.mkdir(exist_ok=True);(OUT/'backgrounds').mkdir(exist_ok=True)
    source=json.loads((DATA/'source_manifest.json').read_text());spots_by_frame,names=read_tracks()
    frames=[int(x) for x in a.frames.split(',')] if a.frames else list(range(385))
    results=[];spot_rows=[];mask_rows=[];frame_rows=[]
    for f in frames:
        spots=spots_by_frame[f]
        labels={c:tifffile.imread(DATA/(a.rfp_dir if c=='rfp' else 'cpsam')/c/f'{f:04d}.tif') for c in ['phase','rfp']}
        assert all(l.shape==(904,1224) for l in labels.values())
        objects={c:mask_objects(labels[c]) for c in labels}
        matches={c:assign_masks(labels[c],objects[c],spots) for c in labels}
        add_nuclear_overlap_matches(labels['phase'],labels['rfp'],objects['rfp'],matches['phase'],matches['rfp'])
        frame_data=dict(frame=f,spots=[],cells=[],nuclei=[],flags=[])
        if f==0:frame_data['flags'].append('Severely defocused first image pair; no source tracked detections.')
        for i,s in enumerate(spots):
            cell,nuc=matches['phase'][i],matches['rfp'][i]
            flags=[]
            if cell['label'] is None:flags.append('cell:'+cell['status'])
            if nuc['label'] is None:flags.append('nucleus:'+nuc['status'])
            if nuc.get('cellOverlap') is not None and nuc['cellOverlap']<.5:flags.append('low_nucleus_cell_overlap')
            s={**s,'cell':cell,'nucleus':nuc,'flags':flags}
            frame_data['spots'].append(s)
            spot_rows.append(dict(frame=f,image_number=f+1,spot_id=s['id'],track_id=s['track'],track_name=s['name'],derived_branch_start_spot_id=s['branch'],x=s['x'],y=s['y'],original_outline=bool(s['original']),cell_mask=cell['label'],cell_status=cell['status'],cell_gap_px=cell['gap'],nucleus_mask=nuc['label'],nucleus_status=nuc['status'],nucleus_gap_px=nuc['gap'],nucleus_cell_overlap=nuc['cellOverlap'],issues=';'.join(flags)))
        for channel,key in [('phase','cells'),('rfp','nuclei')]:
            assigned={ans['label']:i for i,ans in enumerate(matches[channel]) if ans['label'] is not None}
            for obj in objects[channel]:
                idx=assigned.get(obj['label']);public={k:v for k,v in obj.items() if k!='tree'};public['spotIndex']=idx
                frame_data[key].append(public)
                mask_rows.append(dict(frame=f,channel=channel,mask_label=obj['label'],area_px=obj['area'],mask_centroid_x=obj['x'],mask_centroid_y=obj['y'],touches_image_edge=obj['edge'],assigned_spot_id=spots[idx]['id'] if idx is not None else None,assignment_status=matches[channel][idx]['status'] if idx is not None else 'unassigned'))
        display_stats={}
        for c in ['phase','rfp']:
            dest=OUT/'backgrounds'/f'{c}_{f:04d}.jpg'
            display_stats[c]=background(ROOT/source['frames'][f][c]['path'],dest)
        frame_data['displayStats']=display_stats
        nc=sum(m['label'] is not None for m in matches['phase']);nn=sum(m['label'] is not None for m in matches['rfp'])
        frame_data['summary']=dict(tracked=len(spots),cells=nc,nuclei=nn,unassignedCells=len(objects['phase'])-nc,unassignedNuclei=len(objects['rfp'])-nn)
        frame_rows.append(dict(frame=f,image_number=f+1,tracked=len(spots),cell_masks=len(objects['phase']),matched_cells=nc,nuclear_masks=len(objects['rfp']),matched_nuclei=nn,original_missing_outlines=sum(not s['original'] for s in spots),cell_shared_centroids=sum(m['status']=='shared_mask' for m in matches['phase']),nucleus_shared_centroids=sum(m['status']=='shared_mask' for m in matches['rfp']),low_nucleus_cell_overlap=sum('low_nucleus_cell_overlap' in s['flags'] for s in frame_data['spots']),rfp_saturated_fraction=display_stats['rfp']['saturationFraction'],frame_flags=';'.join(frame_data['flags'])))
        results.append(frame_data)
        if f%25==0 or a.pilot:print(json.dumps(frame_rows[-1]),flush=True)
    total=len(spot_rows);cell_status=Counter(r['cell_status'] for r in spot_rows);nuc_status=Counter(r['nucleus_status'] for r in spot_rows)
    summary=dict(frames=len(results),tracked_detections=total,matched_cells=sum(r['cell_mask'] is not None for r in spot_rows),matched_nuclei=sum(r['nucleus_mask'] is not None for r in spot_rows),cell_status=dict(cell_status),nucleus_status=dict(nuc_status),missing_original_outlines=sum(not r['original_outline'] for r in spot_rows),unassigned_cell_masks=sum(r['channel']=='phase' and r['assigned_spot_id'] is None for r in mask_rows),unassigned_nuclear_masks=sum(r['channel']=='rfp' and r['assigned_spot_id'] is None for r in mask_rows),low_nucleus_cell_overlap=sum('low_nucleus_cell_overlap' in r['issues'] for r in spot_rows))
    model_meta={c:json.loads((DATA/(a.rfp_dir if c=='rfp' else 'cpsam')/c/f'{frames[-1]:04d}.json').read_text()) for c in ['phase','rfp']}
    payload=dict(width=1224,height=904,frames=results,summary=summary,palette=PALETTE,trackColors=color_graph(spots_by_frame,'track'),branchColors=color_graph(spots_by_frame,'branch'),source=source,model=model_meta)
    prefix='pilot_' if a.pilot else ''
    (OUT/f'{prefix}viewer_data.json').write_text(json.dumps(payload,separators=(',',':')))
    (OUT/f'{prefix}summary.json').write_text(json.dumps(summary,indent=2))
    dump_csv(OUT/f'{prefix}assignments.csv',spot_rows);dump_csv(OUT/f'{prefix}mask_inventory.csv',mask_rows);dump_csv(OUT/f'{prefix}frame_qc.csv',frame_rows)
    print(json.dumps(summary,indent=2),flush=True)
    if not a.pilot:
        from render_lfct_mcf10a_viewer import render
        render(payload)
if __name__=='__main__':main()
