#!/usr/bin/env python3
"""Regenerate conservative raster digitization; paths are relative to this script."""
import csv, hashlib, json, platform
from pathlib import Path
import numpy as np
import scipy
from scipy import ndimage
from PIL import Image, ImageDraw, __version__ as pillow_version
BASE=Path(__file__).resolve().parent
C=json.loads((BASE/'readings.json').read_text())
IM=Image.open(BASE/C['source_figure']).convert('RGB')
A=np.asarray(IM).astype(int)
YEAR=C['year']; SCALE=C['ymax']/(C['y0']-C['ymax_pixel'])
def value(y): return (C['y0']-y)*SCALE
def calibration_u(y):
    # Worst-case first-order bound from independently +/-0.5 px calibration endpoints.
    t=(C['y0']-y)/(C['y0']-C['ymax_pixel'])
    return SCALE*.5*(abs(1-t)+abs(t))
def write_csv(name,rows,fields):
    with (BASE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
fields=['study','publication_version','figure','panel','cell_line','condition','clone','series','point_id','x_pixel','y_pixel','value','units','evidence_type','pixel_reading_uncertainty_px','reading_uncertainty_units','calibration_uncertainty_units','reported_n','reported_n_definition','error_bar_definition','qc_status','source_file']
rows=[];components=[]
overlay=IM.crop(C['crop']).resize(((C['crop'][2]-C['crop'][0])*4,(C['crop'][3]-C['crop'][1])*4),Image.Resampling.NEAREST)
draw=ImageDraw.Draw(overlay)
def xy(x,y):return ((x-C['crop'][0])*4,(y-C['crop'][1])*4)
def add(clone,j,series,x,y,ev,pid):
    rows.append(dict(study=f'jemaa_{YEAR}',publication_version=C['version'],figure=C['figure'],panel=C['panel'],cell_line=C['cell_line'],condition='diploid' if j<2 else 'tetraploid',clone=clone,series=series,point_id=pid,x_pixel=round(x,3),y_pixel=round(y,3),value=round(value(y),5),units='um/hour',evidence_type=ev,pixel_reading_uncertainty_px=1,reading_uncertainty_units=round(SCALE,5),calibration_uncertainty_units=round(calibration_u(y),5),reported_n=C['n'],reported_n_definition='caption n; does not establish number of displayed cells',error_bar_definition='SEM' if series!='isolated_dot' else 'not applicable',qc_status='reviewed',source_file=C['source_figure']))
for j,(clone,x,mean) in enumerate(zip(C['clones'],C['xs'],C['means'])):
    # Strict color mask; full components are kept before shape filtering.
    if j<2:mask=(A[:,:,1]>A[:,:,0]+22)&(A[:,:,1]>A[:,:,2]+15)&(A[:,:,1]>80)
    else:mask=(A[:,:,0]>A[:,:,1]+45)&(A[:,:,0]>A[:,:,2]+45)&(A[:,:,0]>120)
    region=np.zeros(mask.shape,bool);region[max(0,int(C['ymax_pixel'])+2):int(C['y0'])+4,x-C['radius']:x+C['radius']+1]=True
    labels,n=ndimage.label(mask&region,structure=np.ones((3,3)))
    accepted=0
    for k in range(1,n+1):
        yy,xx=np.where(labels==k);area=len(xx);w=int(xx.max()-xx.min()+1);h=int(yy.max()-yy.min()+1);cx=float(xx.mean());cy=float(yy.mean())
        ok=C['dot_area'][0]<=area<=C['dot_area'][1] and C['dot_width'][0]<=w<=C['dot_width'][1] and C['dot_height'][0]<=h<=C['dot_height'][1] and abs(cy-mean)>6 and area/(w*h)>.5 and xx.min()>x-C['radius'] and xx.max()<x+C['radius']
        status='isolated_candidate' if ok else 'unresolved_or_fragment'
        pid=f'{clone}_component_{k:03d}'
        components.append(dict(clone=clone,component_id=pid,area_pixels=area,x_min=int(xx.min()),x_max=int(xx.max()),y_min=int(yy.min()),y_max=int(yy.max()),x_centroid=round(cx,3),y_centroid=round(cy,3),status=status))
        if ok:
            add(clone,j,'isolated_dot',cx,cy,'digitized individual point',pid);accepted+=1
            q=xy(cx,cy);draw.ellipse([q[0]-12,q[1]-12,q[0]+12,q[1]+12],outline='blue',width=2)
            draw.text((q[0]+13,q[1]-5),str(k),fill='blue')
        elif area>=C['dot_area'][0]:
            q=xy(xx.min()-1,yy.min()-1);r=xy(xx.max()+1,yy.max()+1);draw.rectangle([*q,*r],outline='#888888',width=1)
    add(clone,j,'mean',x,mean,'digitized summary',clone+'_mean')
    q=xy(x-10,mean);r=xy(x+10,mean);draw.line([q,r],fill='blue',width=2)
    for series,y in [('sem_upper_endpoint',C['sem_top'][j]),('sem_lower_endpoint',C['sem_bottom'][j])]:
        if y is not None:
            add(clone,j,series,x,y,'digitized summary',clone+'_'+series);draw.line([xy(x-5,y),xy(x+5,y)],fill='#aa00ff',width=2)
    print(clone,'isolated',accepted,'mean',round(value(mean),3))
# Calibration marks at original plotted axis endpoints, with explicit value labels.
axis_x=214 if YEAR==2017 else 533
for y,v in [(C['y0'],0),(C['ymax_pixel'],C['ymax'])]:
    q=xy(axis_x,y);draw.line([(q[0]-12,q[1]),(q[0]+12,q[1])],fill='blue',width=2)
overlay.save(BASE/'digitization_overlay.png')
IM.crop(C['crop']).save(BASE/'panel_crop.png')
write_csv('digitized_points.csv',rows,fields)
write_csv('component_qc.csv',components,list(components[0]))
measure=[]
for j,clone in enumerate(C['clones']):
    m=value(C['means'][j]);t=C['sem_top'][j];b=C['sem_bottom'][j]
    measure.append(dict(study=f'jemaa_{YEAR}',figure=C['figure'],panel=C['panel'],cell_line=C['cell_line'],condition='diploid' if j<2 else 'tetraploid',clone=clone,metric='mean individual cell speed',value=round(m,5),units='um/hour',evidence_type='digitized summary',summary_definition=('mean horizontal line; caption mean +/- SEM' if YEAR==2017 else 'mean horizontal line; caption states SEM, endpoints unresolved'),sem_upper=round(value(t)-m,5) if t is not None else '',sem_lower=round(m-value(b),5) if b is not None else '',sem_status='resolvable; +/-1 pixel endpoint reading' if t is not None else 'not resolvable; no inferred SEM',reported_n=C['n'],observation_window=C['window'],reading_uncertainty_units=round(SCALE,5),calibration_uncertainty_units=round(calibration_u(C['means'][j]),5),derivation='digitized_points.csv:'+clone+'_mean',isolated_points=sum(r['clone']==clone and r['series']=='isolated_dot' for r in rows)))
write_csv('measurements.csv',measure,list(measure[0]))
inputs={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in [C['source_figure'],C['source_article'],'readings.json']}
protocol=dict(input_sha256=inputs,publication_version=C['version'],figure=C['figure'],panel=C['panel'],context_panel=C['context_panel'],crop_bounds_xyxy=C['crop'],coordinates='original full figure pixels; origin top left; y downward',y_axis=dict(scale='linear',endpoints=[dict(y_pixel=C['y0'],value=0),dict(y_pixel=C['ymax_pixel'],value=C['ymax'])],endpoint_uncertainty_px=.5,units='um/hour'),x_axis='categorical clone identities; jitter has no physical meaning',timing_evidence=C['window'],reported_summary='2017 means +/-SEM n=3; 2023 caption states SEM n=5, central horizontal lines treated as means',uncertainty=dict(reading_px=1,calibration='worst-case propagated +/-0.5px endpoints; separate from biological SEM',component_centroid='color mask centroid, limited by rasterization'),software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pillow=pillow_version),extraction_settings=C,numerical_source_status='No clone numerical migration dataset in downloaded article or supplements; figures used',interpretation='isolated points are incomplete biased subset; do not calculate population means from these or use displayed cell count as independent replicate n')
(BASE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
