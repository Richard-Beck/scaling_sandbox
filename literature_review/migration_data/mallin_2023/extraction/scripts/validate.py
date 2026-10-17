#!/usr/bin/env python
"""Source hashes, coordinate round-trip, summary-bar readings and threshold sensitivity."""
import csv,json,hashlib
from pathlib import Path
import cv2,numpy as np
ROOT=Path(__file__).resolve().parents[1];STUDY=ROOT.parent
p=json.loads((ROOT/'protocol.json').read_text());cfg=json.loads((ROOT/'calibration.json').read_text())
for f,h in p['source_hashes'].items():assert hashlib.sha256((STUDY/f).read_bytes()).hexdigest()==h
rows=list(csv.DictReader((ROOT/'digitized_points.csv').open()));im=cv2.imread(str(STUDY/'figures/mallin2023_fig2.png'),0)
sens=[]
for c in cfg:
    rr=[r for r in rows if r['figure_panel']==c['panel'] and r['condition']==c['condition']]
    x0,y0,x1,y1=c['crop']; scatter='xmax'in c
    for r in rr:
        yzero=c['axis'][3];ymax=c['axis'][4];ytop=c['axis'][5]
        assert abs(yzero-float(r['value'])/ymax*(yzero-ytop)-float(r['y_px']))<.001
    for t in [80,180]:
        mask=(im[y0:y1,x0:x1]<t).astype('uint8')
        if scatter:mask=cv2.morphologyEx(mask,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
        n,lab,stats,cent=cv2.connectedComponentsWithStats(mask,8);pts=[]
        for j in range(1,n):
            x,y,w,h,a=stats[j];ok=(35<=a<=105 and 7<=w<=12 and 7<=h<=12) if scatter else (3<=a<=13 and 2<=w<=5 and 2<=h<=5)
            if ok:pts.append(cent[j]+[x0,y0])
        matched=[]
        for r in rr:
            dist=min([float(np.linalg.norm(q-[float(r['x_px']),float(r['y_px'])])) for q in pts],default=999)
            if dist<=1.5:matched.append(dist)
        sens.append(dict(panel=c['panel'],condition=c['condition'],threshold=t,nominal_resolved_points=len(rr),alternate_resolved_points=len(pts),nominal_matched_within_1_5_px=len(matched),max_matched_shift_px=round(max(matched,default=0),4)))
    if 'summary'in c:
        x,y=c['summary'];lo=int(y)-4;hi=int(y)+5;dark=(im[lo:hi,x-25:x+26]<128).sum(axis=1)
        assert abs(float(np.mean(np.flatnonzero(dark==max(dark))+lo))-y)<=1
with (ROOT/'threshold_sensitivity.csv').open('w',newline='')as f:
    w=csv.DictWriter(f,fieldnames=list(sens[0]));w.writeheader();w.writerows(sens)
(ROOT/'validation.json').write_text(json.dumps(dict(source_hashes_match=True,pixel_value_roundtrip_max_px=.001,summary_bar_repeat_reading_agrees_within_px=1,nominal_points=len(rows),native_image_dimensions_px=[im.shape[1],im.shape[0]],scope='Internal consistency and threshold sensitivity; visual review of overlay and contact sheets performed separately'),indent=2)+'\n')
print('Validated',len(rows),'points and six central bars; wrote threshold_sensitivity.csv')
