"""Reproduce conservative nuclear tracks and figure digitization; run in cpose."""
from pathlib import Path
import csv, hashlib, json
import cv2
import numpy as np
from scipy.optimize import linear_sum_assignment
from PIL import Image, ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
STUDY=OUT.parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,rows):
    if not rows: return
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
cv2.setNumThreads(1)
cap=cv2.VideoCapture(str(STUDY/'source/xuan2018_movie.avi'))
frames=[]
while True:
    ok,f=cap.read()
    if not ok: break
    frames.append(cv2.cvtColor(f,cv2.COLOR_BGR2GRAY))
cap.release()
scale=25/78 # white bar occupies columns 910..987 inclusive; endpoints 909.5,987.5
REVIEWED={10,15,16,17,18,21,24,25}
tracks={}; active={}; nextid=1; detections=[]
for t,im in enumerate(frames):
    smooth=cv2.GaussianBlur(im,(9,9),0)
    mask=(smooth>45).astype('uint8'); mask[1280:1350,880:1020]=0
    n,labels,stats,cents=cv2.connectedComponentsWithStats(mask)
    ds=[]
    for k in range(1,n):
        x,y,w,h,area=stats[k]
        if area<220: continue
        weights=np.maximum(smooth.astype(float)-30,0)*(labels==k)
        yy,xx=np.nonzero(labels==k); ww=weights[yy,xx]
        ds.append(dict(x=float(np.average(xx,weights=ww)),y=float(np.average(yy,weights=ww)),area=int(area),edge=bool(x==0 or y==0 or x+w==1040 or y+h==1392)))
    detections.append(ds)
    old=list(active); matched={}
    if old and ds:
        cost=np.full((len(old),len(ds)),1e6)
        for a,i in enumerate(old):
            p=active[i]
            for b,d in enumerate(ds):
                dist=np.hypot(p['x']-d['x'],p['y']-d['y']); ratio=d['area']/p['area']
                if dist<60 and .5<ratio<2: cost[a,b]=dist+10*abs(np.log(ratio))
        aa,bb=linear_sum_assignment(cost)
        for a,b in zip(aa,bb):
            if cost[a,b]<1e5:
                # Reject competing identities near a merge/crossing. No gap closure.
                alternatives=cost[:,b][cost[:,b]<1e5]
                if len(alternatives)>1 and sorted(alternatives)[1]-sorted(alternatives)[0]<10: continue
                matched[b]=old[a]
    active={}
    for b,d in enumerate(ds):
        i=matched.get(b)
        if i is None: i=nextid; nextid+=1; tracks[i]=[]
        d=dict(d,frame=t,cell_id=f'N{i:03d}')
        tracks[i].append(d); active[i]=d

rows=[]; qc=[]
for i,ts in tracks.items():
    large=i==10
    cls='apparent_PGCC' if large else 'smaller_nucleus_unconfirmed_ploidy'
    status='edge_censored' if any(p['edge'] for p in ts) else ('short_segment' if len(ts)<10 else 'candidate_usable')
    if i in REVIEWED: status='reviewed_full_window'
    for p in ts:
        rows.append(dict(study='xuan_2018',movie='xuan2018_movie.avi',cell_id=p['cell_id'],frame_index=p['frame'],experimental_time_h='',x_px=p['x'],y_px=p['y'],x_um=p['x']*scale,y_um=p['y']*scale,nuclear_mask_area_px=p['area'],classification=cls,qc_status=status,localization_uncertainty_px=3,evidence_type='retracked_representative_movie'))
    qc.append(dict(cell_id=ts[0]['cell_id'],classification=cls,classification_evidence='large bright irregular nucleus, first-frame mask 6565 px; ploidy not measured' if large else 'smaller fluorescent nuclear object; ploidy not measured',first_frame=ts[0]['frame'],last_frame=ts[-1]['frame'],n_frames=len(ts),gaps=0,qc_status=status,exclusions='boundary localization' if status=='edge_censored' else ('fewer than 10 frames' if status=='short_segment' else ''),correction='none',identity_note='terminated at unmatched/ambiguous detection; newly detected segments may represent the same cell'))
write('tracks.csv',rows);write('track_qc.csv',qc)
review=OUT/'review';review.mkdir(exist_ok=True)
writer=cv2.VideoWriter(str(OUT/'tracking_overlay.mp4'),cv2.VideoWriter_fourcc(*'mp4v'),8,(1040,1392))
for t,im in enumerate(frames):
    rgb=cv2.cvtColor(im,cv2.COLOR_GRAY2BGR)
    for i,ts in tracks.items():
        past=[p for p in ts if p['frame']<=t]
        now=[p for p in ts if p['frame']==t]
        if not now: continue
        col=(0,0,255) if i==10 else (0,255,255)
        pts=np.array([[p['x'],p['y']] for p in past],dtype=np.int32)
        if len(pts)>1:cv2.polylines(rgb,[pts],False,col,2)
        x,y=pts[-1];cv2.circle(rgb,(x,y),5,col,2);cv2.putText(rgb,f'N{i:03d}',(x+6,y),cv2.FONT_HERSHEY_SIMPLEX,.5,col,1)
    cv2.putText(rgb,f'frame {t}; experimental time unresolved',(20,25),cv2.FONT_HERSHEY_SIMPLEX,.6,(0,255,255),1)
    writer.write(rgb)
    sheetidx=t//6; slot=t%6
    if slot==0:sheet=Image.new('RGB',(1560,1392))
    thumb=Image.fromarray(cv2.cvtColor(rgb,cv2.COLOR_BGR2RGB)).resize((520,696))
    sheet.paste(thumb,((slot%3)*520,(slot//3)*696))
    if slot==5 or t==len(frames)-1:sheet.save(review/f'frames_{sheetidx*6:02d}_{t:02d}.jpg',quality=90)
writer.release()
cal=Image.fromarray(frames[0]).convert('RGB');draw=ImageDraw.Draw(cal)
draw.line((909.5,1291,987.5,1291),fill='red',width=2);draw.text((850,1250),'25 um / 78 px',fill='red');cal.save(OUT/'scale_calibration.png')
metrics=[];msd=[]
for i,ts in tracks.items():
    if i not in REVIEWED:continue
    xy=np.array([[p['x'],p['y']] for p in ts])*scale
    steps=np.linalg.norm(np.diff(xy,axis=0),axis=1);path=steps.sum();net=np.linalg.norm(xy[-1]-xy[0])
    for metric,value,unit in [('path_length',path,'um'),('net_displacement',net,'um'),('directness',net/path if path else 0,'dimensionless'),('sampled_step',steps.mean(),'um/frame')]:
        metrics.append(dict(cell_id=f'N{i:03d}',metric=metric,value=float(value),units=unit,first_frame=ts[0]['frame'],last_frame=ts[-1]['frame'],timing_scenario='frame_index',evidence_type='retracked_representative_movie',definition='path=sum consecutive distances; net=end-start; directness=net/path; sampled_step=path/(n-1)'))
    for label,dt in [('methods_10min',1/6),('12h_endpoints',12/73),('supplement_16h_endpoints',16/73)]:
        metrics.append(dict(cell_id=f'N{i:03d}',metric='sampled_speed',value=float(steps.mean()/dt),units='um/h',first_frame=ts[0]['frame'],last_frame=ts[-1]['frame'],timing_scenario=label,evidence_type='retracked_representative_movie',definition='conditional scenario only; constant spacing assumed'))
    for lag in range(1,len(ts)):
        vals=np.sum((xy[lag:]-xy[:-lag])**2,axis=1)
        msd.append(dict(cell_id=f'N{i:03d}',lag_frames=lag,n_pairs=len(vals),msd_um2=float(vals.mean()),first_frame=ts[0]['frame'],last_frame=ts[-1]['frame'],evidence_type='retracked_representative_movie'))
write('measurements.csv',metrics);write('msd.csv',msd)
fig,ax=plt.subplots(figsize=(9,5))
for i,ts in tracks.items():
    if i not in REVIEWED:continue
    xy=np.array([[p['x'],p['y']] for p in ts])*scale
    ax.plot([p['frame'] for p in ts],np.linalg.norm(xy-xy[0],axis=1),label=f'N{i:03d}',linewidth=3 if i==10 else 1)
ax.set(xlabel='Original frame index',ylabel='Displacement from segment start (um)');ax.legend(ncol=4,fontsize=6);fig.tight_layout();fig.savefig(OUT/'displacement_vs_frame.png');plt.close(fig)

# Figure 4: original full-image coordinates, manual axis anchors and bar readings.
fpath=STUDY/'figures/xuan2018_fig4.png';im=np.array(Image.open(fpath).convert('RGB')); overlay=Image.fromarray(im);draw=ImageDraw.Draw(overlay)
points=[]
for panel,xs,y0,ymax,vmax,means,lows,highs,units,metric in [
 ('B',[755,964],469,52,30,[147,264],[216,296],[79,231],'um/h','average_speed'),
 ('C',[1331,1540],469,52,.8,[235,105],[261,130],[211,80],'fraction','slow_mode_fraction'),
 ('E',[758,966],1002,586,1.6,[725,619],[739,650],[711,590],'dimensionless','MSD_logarithmic_slope')]:
    for j,condition in enumerate(['Non-PGCC','PGCC']):
        for kind,y in [('mean',means[j]),('SEM_lower',lows[j]),('SEM_upper',highs[j])]:
            value=(y0-y)/(y0-ymax)*vmax
            points.append(dict(study='xuan_2018',version='2018_version_of_record',figure='4',panel=panel,condition=condition,series=metric,point_kind=kind,x_px=xs[j],y_px=y,x_value='',y_value=value,x_units='',y_units=units,pixel_reading_uncertainty=1.5,axis_calibration_uncertainty_px=1,value_uncertainty=value*0+vmax/(y0-ymax)*1.803,evidence_type='digitized_summary',source_sha256=sha(fpath)))
            draw.ellipse((xs[j]-4,y-4,xs[j]+4,y+4),outline='lime',width=2)
    draw.line((xs[0]-20,y0,xs[-1]+20,y0),fill='lime');draw.line((xs[0]-20,ymax,xs[-1]+20,ymax),fill='lime')
for condition,col in [('Non-PGCC','blue'),('PGCC','red')]:
    for x in range(174,433):
        rgb=im[610:930,x,:].astype(float)
        mask=(rgb[:,2]>150)&(rgb[:,0]<110)&(rgb[:,1]<110) if col=='blue' else (rgb[:,0]>150)&(rgb[:,1]<110)&(rgb[:,2]<110)
        yy=np.where(mask)[0]+610
        # exclude legend and keep trajectory band; trace starts below legend
        yy=yy[yy>650] if x<300 else yy
        if len(yy)==0:continue
        y=float(np.median(yy)); xv=10**(3+(x-206)/(347-206));yv=10**((1002-y)/(1002-586)*4)
        points.append(dict(study='xuan_2018',version='2018_version_of_record',figure='4',panel='D',condition=condition,series='MSD',point_kind='plotted_curve_sample',x_px=x,y_px=y,x_value=xv,y_value=yv,x_units='s',y_units='um2',pixel_reading_uncertainty=1.5,axis_calibration_uncertainty_px=1,value_uncertainty=yv*(10**(1.803*4/416)-1),evidence_type='digitized_curve',source_sha256=sha(fpath)))
        draw.point((x,int(y)),fill='lime')
write('digitized_points.csv',points);overlay.save(OUT/'figure4_extraction_overlay.png')
for p in points:
    if p['panel']=='D':continue
    metrics.append(dict(cell_id=p['condition'],metric=p['series']+'_'+p['point_kind'],value=p['y_value'],units=p['y_units'],first_frame='',last_frame='',timing_scenario='published_population_summary',evidence_type='digitized_summary',definition=f"Figure 4{p['panel']} {p['point_kind']}; source pixel/calibration uncertainty in digitized_points.csv; mean ± SEM"))
write('measurements.csv',metrics)
protocol=dict(study='xuan_2018',software=dict(opencv=cv2.__version__,numpy=np.__version__),inputs={str(p.relative_to(STUDY)):sha(p) for p in [STUDY/'source/xuan2018_movie.avi',fpath,STUDY/'derived/xuan2018.txt',STUDY/'derived/xuan2018_supp.txt']},frames=len(frames),dimensions=[1040,1392],duplicate_frames=[i for i in range(1,len(frames)) if np.array_equal(frames[i],frames[i-1])],scale=dict(bar_endpoints=[[909.5,1291],[987.5,1291]],length_um=25,um_per_px=scale,endpoint_uncertainty_px=1),timing=dict(status='unresolved',main_methods='10 minute intervals for 12 hours',supplement='16 hour time lapse',frame_count_note='74 unique frames implies 73 intervals, not an exact match for 12h at 10min',physical_tracks_time='blank',scenarios=['10min per frame','12h / 73 intervals','16h / 73 intervals']),tracking=dict(smoothing='9x9 Gaussian',threshold=45,min_area_px=220,max_step_px=60,area_ratio=[.5,2],ambiguity_cost_margin=10,centroid='intensity weighted above background 30',gap_closure=False,scale_bar_exclusion=[880,1280,1020,1350],stage_drift_correction='none; no defensible fixed image features'),figure_calibration=dict(D=dict(x_scale='log10',x_anchors=[[206,1000],[347,10000]],y_scale='log10',y_anchors=[[1002,1],[586,10000]],crop=[133,586,485,1003]),summary='mean ± SEM; Fig4D caption n=5, 50+ cells; general methods at least 3 experiments',review_status='candidate tracks; see METHODS and review sheets'))
protocol['review_status']='all 74 frames inspected in contact sheets; 8 full-window tracks accepted; other segments require further identity review'
frame_seen={}; repeated=[]
for index,frame in enumerate(frames):
    frame_hash=hashlib.sha256(frame.tobytes()).hexdigest()
    if frame_hash in frame_seen: repeated.append({'frame':index,'identical_to':frame_seen[frame_hash]})
    else:frame_seen[frame_hash]=index
protocol['all_exact_repeated_frames']=repeated
protocol['reviewed_full_window_ids']=[f'N{i:03d}' for i in sorted(REVIEWED)]
protocol['localization_uncertainty_note']='3 px is a conservative working sensitivity allowance, not a measured confidence interval'
protocol['figure_panel_crops_xyxy']={'B':[550,0,1075,510],'C':[1110,0,1650,510],'D':[0,510,510,1070],'E':[545,510,1080,1070]}
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
comparison=[]
for i in sorted(REVIEWED):
    m={r['metric']:r['value'] for r in metrics if r['cell_id']==f'N{i:03d}' and r['timing_scenario']=='frame_index'}
    comparison.append(dict(cell_id=f'N{i:03d}',classification='apparent_PGCC' if i==10 else 'smaller_nucleus_unconfirmed_ploidy',first_frame=0,last_frame=73,n_frames=74,**m))
write('equal_window_comparison.csv',comparison)
print('tracks',len(tracks),'rows',len(rows),'MSD curve points',sum(p['panel']=='D' for p in points))
