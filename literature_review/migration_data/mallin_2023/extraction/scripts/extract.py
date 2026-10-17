#!/usr/bin/env python
"""Conservative connected-glyph figure digitization and streaming-video audit."""
import csv, hashlib, json, math
from pathlib import Path
import cv2
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; STUDY=ROOT.parent
FIG=STUDY/'figures/mallin2023_fig2.png'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def csvout(name, rows):
    if rows:
        with (ROOT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
# Global native PNG pixel coordinates, origin upper left; axis endpoints read manually.
CFG=[
 dict(panel='2A',condition='non_PACC',crop=[182,190,218,377],axis=[0,0,142,382,1250,77],metric='accumulated_distance',unit='um',summary=[200,304]),
 dict(panel='2A',condition='PACC',crop=[269,99,306,378],axis=[0,0,142,382,1250,77],metric='accumulated_distance',unit='um',summary=[287,341.5]),
 dict(panel='2B',condition='non_PACC',crop=[535,230,570,386],axis=[0,0,495,389,300,83],metric='euclidean_distance',unit='um',summary=[553,349.5]),
 dict(panel='2B',condition='PACC',crop=[622,104,658,386],axis=[0,0,495,389,300,83],metric='euclidean_distance',unit='um',summary=[640,320]),
 dict(panel='2D',condition='non_PACC',crop=[1789,251,1826,413],axis=[0,0,1749,416,1,110],metric='directness',unit='1',summary=[1807,375.5]),
 dict(panel='2D',condition='PACC',crop=[1876,152,1913,413],axis=[0,0,1749,416,1,110],metric='directness',unit='1',summary=[1895,289.5]),
 dict(panel='2E',condition='non_PACC',crop=[96,675,412,916],axis=[76,463,0,920,150,646],xmax=8000,metric='euclidean_distance',unit='um'),
 dict(panel='2E',condition='PACC',crop=[586,663,947,919],axis=[582,968,0,925,500,649],xmax=150000,metric='euclidean_distance',unit='um'),
 dict(panel='2F',condition='non_PACC',crop=[1099,681,1429,917],axis=[1085,1479,0,923,.6,659],xmax=8000,metric='directness',unit='1'),
 dict(panel='2F',condition='PACC',crop=[1561,682,1943,917],axis=[1559,1969,0,920,1,646],xmax=150000,metric='directness',unit='1'),
]
(ROOT/'calibration.json').write_text(json.dumps(CFG,indent=2)+'\n')
a=np.asarray(cv2.imread(str(FIG),cv2.IMREAD_GRAYSCALE)); overlay=cv2.imread(str(FIG))
rows=[]; unresolved=[]; summaries=[]; counts=[]
for c in CFG:
    x0,y0,x1,y1=c['crop'];scatter='xmax' in c
    mask=(a[y0:y1,x0:x1]<128).astype('uint8')
    if scatter:
        # Remove thin regression line without deconvolving any merged glyphs.
        mask=cv2.morphologyEx(mask,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
    n,lab,stats,cent=cv2.connectedComponentsWithStats(mask,8)
    accepted=[]; reject=0
    for j in range(1,n):
        bx,by,bw,bh,area=stats[j];x,y=cent[j]+[x0,y0]
        ok=(35<=area<=105 and 7<=bw<=12 and 7<=bh<=12) if scatter else (3<=area<=13 and 2<=bw<=5 and 2<=bh<=5)
        manual_ambiguous=(c['panel']=='2A' and c['condition']=='PACC' and any(abs(x-px)<1 and abs(y-py)<1 for px,py in [(290.5,296.5),(290.4545,370.3636)]))
        if manual_ambiguous:ok=False
        if not ok:
            if area>=3:
                reject+=1
                unresolved.append(dict(panel=c['panel'],condition=c['condition'],x_min_px=int(bx+x0),y_min_px=int(by+y0),width_px=int(bw),height_px=int(bh),dark_area_px=int(area),status='manual_review_possible_merged_glyph' if manual_ambiguous else 'unresolved_merged_or_line_connected_glyph',inferred_cell_count=''))
                cv2.rectangle(overlay,(bx+x0,by+y0),(bx+x0+bw,by+y0+bh),(0,128,255),1)
            continue
        accepted.append((float(x),float(y),int(area)))
    for k,(x,y,area) in enumerate(sorted(accepted,key=lambda q:q[1]),1):
        ax0,ax1,unused,yzero,ymax,ytop=c['axis'];sy=ymax/(yzero-ytop);value=(yzero-y)*sy
        denom=yzero-ytop
        calerr=abs(sy)*math.sqrt(((y-ytop)/denom)**2+((yzero-y)/denom)**2)
        area_value=(x-ax0)/(ax1-ax0)*c['xmax'] if scatter else ''
        sx=c['xmax']/(ax1-ax0) if scatter else 0
        xcal=abs(sx)*math.sqrt(((x-ax1)/(ax1-ax0))**2+((x-ax0)/(ax1-ax0))**2) if scatter else ''
        pid=f"{c['panel']}_{c['condition']}_{k:03d}"
        rows.append(dict(study='mallin_2023',publication_version='2023_version_of_record',figure_panel=c['panel'],condition=c['condition'],point_id=pid,evidence_type='digitized_individual_point',metric=c['metric'],x_px=round(x,4),y_px=round(y,4),cell_area_um2=round(area_value,3) if scatter else '',value=round(value,6),units=c['unit'],reading_uncertainty_px=1,reading_uncertainty_value=round(sy,6),axis_endpoint_uncertainty_px=1,calibration_uncertainty_value=round(calerr,6),area_reading_uncertainty_um2=round(sx,3) if scatter else '',area_calibration_uncertainty_um2=round(xcal,3) if scatter else '',glyph_area_px=area,qc='resolved_glyph_partial_distribution',source='figures/mallin2023_fig2.png',source_sha256=sha(FIG)))
        cv2.circle(overlay,(round(x),round(y)),7 if scatter else 4,(0,180,0),1)
        if scatter: cv2.putText(overlay,str(k),(round(x)+6,round(y)-4),cv2.FONT_HERSHEY_SIMPLEX,.25,(150,30,200),1,cv2.LINE_AA)
    counts.append(dict(panel=c['panel'],condition=c['condition'],resolved_points=len(accepted),unresolved_regions=reject,reported_n='50' if c['condition']=='non_PACC' else '49',n_scope='caption explicitly attaches n only to panel 2C; extrapolation to this panel unverified'))
    if 'summary' in c:
        x,y=c['summary'];ax0,ax1,u,yzero,ymax,ytop=c['axis'];sy=ymax/(yzero-ytop)
        summaries.append(dict(study='mallin_2023',publication_version='2023_version_of_record',figure_panel=c['panel'],condition=c['condition'],metric=c['metric'],units=c['unit'],value=round((yzero-y)*sy,6),x_px=x,y_px=y,reading_uncertainty_px=1,reading_uncertainty_value=round(sy,6),axis_endpoint_uncertainty_px=1,calibration_uncertainty_value=round(sy*math.sqrt(((y-ytop)/(yzero-ytop))**2+((yzero-y)/(yzero-ytop))**2),6),evidence_type='digitized_summary',summary_definition='horizontal_central_bar_statistic_unspecified_in_caption_and_methods',observation_window_h=24,source='figures/mallin2023_fig2.png',source_sha256=sha(FIG)))
        cv2.drawMarker(overlay,(round(x),round(y)),(255,0,0),cv2.MARKER_CROSS,15,2)
    cv2.rectangle(overlay,(x0,y0),(x1,y1),(255,160,0),1)
    # Full un-resampled crops accompany full overlay.
    cv2.imwrite(str(ROOT/'overlays'/f"{c['panel']}_{c['condition']}.png"),overlay[y0:y1,x0:x1])
cv2.imwrite(str(ROOT/'overlays'/'figure2_annotated.png'),overlay)
csvout('digitized_points.csv',rows);csvout('measurements.csv',summaries);csvout('unresolved_regions.csv',unresolved);csvout('coverage.csv',counts)
# Video streaming-frame changes. Adjacent MAD is diagnostic, not proof of acquisition identity.
movies=[]; framediffs=[]
for num in [1,2]:
    p=STUDY/'source'/f'mallin_movie{num}.mp4';cap=cv2.VideoCapture(str(p));fps=cap.get(cv2.CAP_PROP_FPS)
    total=int(cap.get(cv2.CAP_PROP_FRAME_COUNT));prev=None; seen=set();exact_dup=0;ds=[];samples=[]
    wanted=set(np.linspace(0,total-1,6,dtype=int));idx=0
    while True:
        ok,f=cap.read()
        if not ok:break
        digest=hashlib.sha256(f.tobytes()).hexdigest()
        if digest in seen:exact_dup+=1
        seen.add(digest)
        g=cv2.cvtColor(f,cv2.COLOR_BGR2GRAY)
        if prev is not None:
            diff=np.abs(g.astype('float32')-prev.astype('float32'));mad=float(diff.mean());ds.append(mad)
            framediffs.append(dict(movie=p.name,frame_index=idx,playback_time_s=round(idx/fps,6),adjacent_grayscale_MAD=round(mad,6),fraction_pixels_change_gt10=round(float((diff>10).mean()),6),experimental_time_h='',qc='streaming_frame_not_acquisition_frame'))
        if idx in wanted:
            tile=cv2.resize(f,(512,384));cv2.putText(tile,f'frame {idx}, playback {idx/fps:.2f} s',(8,24),cv2.FONT_HERSHEY_SIMPLEX,.6,(0,255,255),1);samples.append(tile)
        prev=g;idx+=1
    cap.release()
    cv2.imwrite(str(ROOT/'overlays'/f'movie{num}_contact_sheet.jpg'),np.vstack([np.hstack(samples[:3]),np.hstack(samples[3:])]))
    movies.append(dict(movie=p.name,sha256=sha(p),decoded_frames=idx,container_frames=total,width_px=int(g.shape[1]),height_px=int(g.shape[0]),playback_fps=fps,playback_duration_s=idx/fps,unique_exact_decoded_images=len(seen),exact_duplicate_decoded_images=exact_dup,adjacent_MAD_min=float(min(ds)),adjacent_MAD_median=float(np.median(ds)),adjacent_MAD_max=float(max(ds)),acquisition_interval_reported_min=30,experiment_duration_reported_h=24,acquisition_frame_count='unresolved_streaming_repetition_and_encoding',spatial_calibration_um_per_px=None))
csvout('movie_frame_diagnostics.csv',framediffs)
protocol=dict(study='mallin_2023',publication_version='2023_version_of_record',doi='10.1007/s10585-023-10216-8',software=dict(python='existing cpose Conda environment',opencv=cv2.__version__,numpy=np.__version__),source_hashes={str(p.relative_to(STUDY)):sha(p) for p in [FIG,STUDY/'source/mallin2023.html',STUDY/'derived/mallin2023.txt']},axis_scale='linear',calibration='calibration.json; native image 2003 x 999; axis endpoints reading uncertainty 1 px',point_algorithm='dark threshold <128; scatter 5x5 elliptical opening removes thin regression line; connected components; strict single glyph size exclusion; no merged-point deconvolution',point_uncertainty='reading +/-1 pixel separately from propagated independent 1-pixel axis endpoint sensitivity; these are sensitivity bounds, not biological variability or confidence intervals',movies=movies,summary_definition='unresolved: original caption and methods do not specify whether central bars are means or medians',sample_size='50 non-PACC and 49 PACC stated for 2C only',cross_panel_pairing='not recoverable; paired cell area and motility only within an individual scatter panel')
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
# Stable large image changes distinguish holds from candidate acquisition transitions.
groups=[]
for m in movies:
    d=[r for r in framediffs if r['movie']==m['movie']]
    starts=[0]+[r['frame_index'] for r in d if r['adjacent_grayscale_MAD']>.2]
    m['candidate_image_states_MAD_gt_0_2']=len(starts)
    m['candidate_image_states_MAD_gt_0_5']=1+sum(r['adjacent_grayscale_MAD']>.5 for r in d)
    m['candidate_image_states_MAD_gt_1_0']=1+sum(r['adjacent_grayscale_MAD']>1 for r in d)
    for k,st in enumerate(starts):
        en=starts[k+1]-1 if k+1<len(starts) else m['decoded_frames']-1
        groups.append(dict(movie=m['movie'],candidate_state_index=k,start_streaming_frame=st,end_streaming_frame=en,experimental_time_h='',qc='candidate_acquisition_state_threshold_0_2_requires_original_movie_verification'))
csvout('movie_candidate_states.csv',groups)
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps(dict(coverage=counts,movies=movies),indent=2))
