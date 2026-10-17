"""Mask geometry and conservative one-to-one assignment, adapted from the validated MCF10A prototype."""
import math
from collections import defaultdict
import cv2
import numpy as np
from scipy.spatial import cKDTree

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


def assign_masks(labels,objects,spots,max_gap=15.0,margin=3.0):
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
        if not candidates or candidates[0][0]>max_gap:continue
        gap,label=candidates[0];answers[i]['gap']=round(gap,3)
        if label in occupied:answers[i]['status']='nearest_mask_occupied';continue
        if len(candidates)>1 and candidates[1][0]-gap<margin:
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



def pixel_square_gap(a,b):
    """Exact distance between unions of closed unit pixel squares (touching = 0).

    Boundary pixel centers give an initial candidate. Any better square pair
    must have center distance within sqrt(2) of that candidate's center distance.
    The second query evaluates every such possible improving pair.
    """
    pa=a['tree'].data;pb=b['tree'].data
    distances,indices=b['tree'].query(pa)
    d=float(distances.min());bound=d+math.sqrt(2)+1e-9
    possible=pa[distances<=bound]
    best=float('inf')
    for point,neighbors in zip(possible,b['tree'].query_ball_point(possible,bound)):
        if not neighbors:continue
        delta=np.maximum(np.abs(pb[neighbors]-point)-1,0)
        best=min(best,float(np.sqrt((delta*delta).sum(axis=1)).min()))
        if best==0:break
    return best

def point_square_gap(obj,point):
    pts=obj['tree'].data;d,_=obj['tree'].query(point)
    indices=obj['tree'].query_ball_point(point,float(d)+math.sqrt(.5)+1e-9)
    delta=np.maximum(np.abs(pts[indices]-point)-.5,0)
    return float(np.sqrt((delta*delta).sum(axis=1)).min())

def bounds(obj):
    pts=obj['tree'].data
    return np.r_[pts.min(axis=0)-.5,pts.max(axis=0)+.5]

def box_gap(a,b):
    return float(np.linalg.norm(np.maximum(np.maximum(a[:2]-b[2:],b[:2]-a[2:]),0)))

def clearances(objects):
    by={o['label']:o for o in objects};boxes={k:bounds(o) for k,o in by.items()};cache={};result={}
    for lab,obj in by.items():
        best=float('inf');neighbor=None
        candidates=sorted((box_gap(boxes[lab],boxes[k]),k) for k in by if k!=lab)
        for lower,k in candidates:
            if lower>best:break
            pair=tuple(sorted((lab,k)))
            if pair not in cache:cache[pair]=pixel_square_gap(obj,by[k])
            gap=cache[pair]
            if gap<best:best,neighbor=gap,k
            if best==0:break
        result[lab]=(None if neighbor is None else best,neighbor)
    return result
