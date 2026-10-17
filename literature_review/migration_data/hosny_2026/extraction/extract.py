#!/usr/bin/env python3
"""Regenerate reviewed Hosny preprint plot readings. Python: scaling_sandbox env."""
from pathlib import Path
import json,csv,hashlib,re,html,platform
import numpy as np
from PIL import Image,ImageDraw
from scipy.ndimage import label,find_objects
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parent
VERSION='bioRxiv 2025.09.03.673887v1 (posted 2025-09-05)'
def write(name,rows):
 if rows:
  fields=list(dict.fromkeys(k for r in rows for k in r))
  with (ROOT/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# Full figure pixel coordinates. Entries: fig,panel,metric,axis_y_low,axis_y_high,value_low,value_high,[(x,median_y,Q1_y,Q3_y)]
P=[
(2,'I','accumulated_distance',1038,860,0,500,[(117,976,989,959),(153,1003,1010,987),(189,963,977,936)]),
(2,'J','euclidean_distance',1038,860,0,300,[(375,990,1006,961),(412,1007,1017,989),(448,1023,1029,1011)]),
(2,'K','directness',1038,860,0,1,[(620,946,972,906),(657,983,998,957),(694,1002,1016,992)]),
(3,'C','parallel_FMI',921,693,-1,1,[(679,779,795,751),(726,796,823,773),(775,803,831,791)]),
(3,'D','perpendicular_FMI',921,693,-1,1,[(1033,806,844,782),(1082,807,833,768),(1129,808,835,789)]),
(5,'J','accumulated_distance',1151,930,0,400,[(169,1041,1077,1009),(221,1058,1088,1035),(271,1083,1108,1054),(323,1105,1127,1077)]),
(5,'K','euclidean_distance',1151,930,0,250,[(614,1074,1097,1044),(666,1103,1115,1070),(717,1098,1126,1081),(768,1126,1138,1114)]),
(5,'L','directness',1151,930,0,1,[(1030,1042,1077,1011),(1081,1081,1101,1045),(1133,1070,1107,1030),(1184,1094,1124,1072)]),
(6,'B','parallel_FMI',324,60,-1,1,[(657,163,191,127),(718,183,210,155),(780,175,199,154),(841,187,211,172)]),
(6,'C','perpendicular_FMI',321,56,-1,1,[(1051,184,216,150),(1113,193,219,162),(1174,185,210,161),(1236,184,203,153)]),
(7,'E','accumulated_distance',638,380,0,500,[(211,528,555,500),(285,539,570,512)]),
(7,'F','euclidean_distance',634,377,0,150,[(673,534,578,517),(748,553,588,506)]),
(7,'G','directness',633,374,0,1,[(1136,586,603,543),(1211,577,611,553)]),
(7,'I','parallel_FMI',1056,780,-1,1,[(725,918,943,901),(804,928,954,891)]),
(7,'J','perpendicular_FMI',1056,779,-1,1,[(1126,914,943,886),(1201,917,949,896)]),
(8,'J','accumulated_distance',1151,930,0,400,[(170,1050,1075,1018),(222,1059,1083,1040),(273,1077,1089,1028),(325,1099,1113,1072)]),
(8,'K','euclidean_distance',1151,930,0,250,[(614,1074,1097,1044),(666,1103,1115,1070),(717,1098,1126,1081),(768,1126,1138,1114)]),
(8,'L','directness',1151,930,0,1,[(1030,1039,1084,1009),(1082,1092,1120,1056),(1133,1064,1101,1035),(1184,1096,1120,1071)]),
(9,'B','parallel_FMI',319,55,-1,1,[(655,160,188,140),(716,182,217,160),(778,171,194,151),(839,186,221,173)]),
(9,'C','perpendicular_FMI',317,53,-1,1,[(1051,183,206,157),(1113,181,215,152),(1174,182,215,159),(1236,184,207,158)]),
(10,'E','accumulated_distance',650,424,0,500,[(203,557,576,534),(270,573,594,548)]),
(10,'F','euclidean_distance',648,424,0,200,[(652,605,626,573),(717,614,628,595)]),
(10,'G','directness',648,423,0,1,[(1095,604,633,572),(1162,615,635,594)]),
(10,'I','parallel_FMI',1081,820,-1,1,[(712,955,988,919),(787,948,976,918)]),
(10,'J','perpendicular_FMI',1079,817,-1,1,[(1109,942,968,914),(1194,948,989,920)])]
rows=[];cal=[]
(ROOT/'overlays').mkdir(exist_ok=True)
for fig in sorted({p[0] for p in P}):
 source=STUDY/f'figures/hosny_preprint_fig{fig}.jpg';im=Image.open(source).convert('RGB');d=ImageDraw.Draw(im)
 for f,p,m,yl,yh,vl,vh,points in [p for p in P if p[0]==fig]:
  if fig in [2,3]:groups=[('PACC','baseline_hypoxia'),('Non-PACC','baseline_hypoxia'),('PC3','normoxia')]
  elif fig in [5,6,8,9]:groups=[('PACC','siNTC_hypoxia'),('PACC','siHIF_hypoxia' if fig in [5,6] else 'siRhoA_hypoxia'),('Non-PACC','siNTC_hypoxia'),('Non-PACC','siHIF_hypoxia' if fig in [5,6] else 'siRhoA_hypoxia')]
  else:groups=[('PC3','siNTC_normoxia'),('PC3','siHIF_normoxia' if fig==7 else 'siRhoA_normoxia')]
  bounds=[min(q[0] for q in points)-40,yh-8,max(q[0] for q in points)+40,yl+8]
  cal.append(dict(figure=f,panel=p,source=str(source.relative_to(STUDY)),crop_bounds=bounds,y_axis={'pixel_low':yl,'pixel_high':yh,'value_low':vl,'value_high':vh,'scale':'linear'},reading_uncertainty_px=3,axis_anchor_uncertainty_px=1))
  d.rectangle(bounds,outline='green',width=2)
  for (x,med,q1,q3),(group,cond) in zip(points,groups):
   for stat,y in [('median',med),('Q1',q1),('Q3',q3)]:
    value=vl+(yl-y)/(yl-yh)*(vh-vl);u=3*abs((vh-vl)/(yl-yh));cu=abs((vh-vl)/(yl-yh))*(abs(y-yl)+abs(yh-y))/abs(yl-yh)
    rows.append(dict(study='hosny_2026',version=VERSION,figure=f,panel=p,condition=cond,group=group,metric=m,statistic=stat,x_pixel=x,y_pixel=y,value=round(value,6),units='um' if 'distance' in m else 'dimensionless',n_cells=60,n_experiments=3 if fig in [2,3] else '',observation_hours=16,observation_definition='distances rescaled to 16 h by 16/T for shorter tracks; T>=1 h',acquisition_interval_hours=0.5,evidence_type='digitized_summary',reading_uncertainty=round(u,6),calibration_uncertainty=round(cu,6),uncertainty_definition='absolute sensitivity bounds; not biological uncertainty',source=str(source.relative_to(STUDY)),source_sha256=sha(source)))
    d.line((x-14,y,x+14,y),fill='lime' if stat=='median' else 'orange',width=2)
   d.text((x-15,yl+10),group,fill='green')
 im.save(ROOT/f'overlays/F{fig}_summaries.png')
write('digitized_points.csv',rows)
write('measurements.csv',[{k:v for k,v in r.items() if k not in ['x_pixel','y_pixel']} for r in rows])
# Printed correlation values take precedence over trying to reconstruct full clouds.
assoc=[('2','L','PACC','directness',-0.02492),('2','M','Non-PACC','directness',-0.02701),('3','E','PACC','parallel_FMI',0.03887),('3','F','Non-PACC','parallel_FMI',0.07138),('S3','A','PACC','accumulated_distance',-0.07611),('S3','B','Non-PACC','accumulated_distance',-0.04518),('S3','C','PACC','euclidean_distance',-0.05301),('S3','D','Non-PACC','euclidean_distance',0.07708),('S4','A','PACC','perpendicular_FMI',0.01156),('S4','B','Non-PACC','perpendicular_FMI',-0.08941)]
write('size_associations.csv',[dict(study='hosny_2026',version=VERSION,figure=f,panel=p,group=g,x_metric='cell_area_um2',y_metric=m,spearman_rho=r,significance='ns (printed)',evidence_type='author_reported_numerical_value',source=(f'figures/hosny_preprint_fig{f}.jpg' if not f.startswith('S') else f'derived/hosny_supp_image{f[1:]}.tiff')) for f,p,g,m,r in assoc])
# Visually audited isolated dot components only, never deblend overlapping dots.
scatter=[('S3','A',107,387,66,253,1000,5000,0,500,'PACC','accumulated_distance',[]),('S3','B',573,854,66,253,0,2000,0,400,'Non-PACC','accumulated_distance',[]),('S3','C',107,387,423,610,1000,5000,0,300,'PACC','euclidean_distance',[]),('S3','D',573,854,423,610,0,2000,0,250,'Non-PACC','euclidean_distance',[]),('S4','A',89,404,88,298,1000,5000,-1,1,'PACC','perpendicular_FMI',[(152,195,298,255)]),('S4','B',594,908,86,295,0,2000,-1,1,'Non-PACC','perpendicular_FMI',[(640,193,825,253)])]
paired=[];qc=[]
for f in ['S3','S4']:
 source=STUDY/f'derived/hosny_supp_image{f[1:]}.tiff';im=Image.open(source).convert('RGB');gray=np.asarray(im.convert('L'));d=ImageDraw.Draw(im)
 for _,p,x0,x1,y0,y1,xlo,xhi,ylo,yhi,g,m,masks in [s for s in scatter if s[0]==f]:
  binary=gray<100
  keep=np.zeros(binary.shape,bool);keep[y0+5:y1-5,x0+6:x1-5]=True
  if f=='S4':keep[188:197,:]=False # axis/zero line; do not fabricate points intersecting it
  for a,b,c,e in masks:keep[b:e,a:c]=False;d.rectangle((a,b,c,e),outline='red')
  lab,n=label(binary&keep);cnt=0;reject=0
  for sl in find_objects(lab):
   if sl is None:continue
   yy,xx=sl;region=lab[sl];identifier=region.max();coords=np.argwhere(region==identifier);area=len(coords);w=xx.stop-xx.start;h=yy.stop-yy.start
   if not(8<=area<=42 and 3<=w<=8 and 3<=h<=8 and .55<=w/h<=1.8):reject+=1;continue
   x=xx.start+coords[:,1].mean();y=yy.start+coords[:,0].mean();cnt+=1
   paired.append(dict(study='hosny_2026',version=VERSION,figure=f,panel=p,group=g,condition='baseline_hypoxia',point_id=f'{f}{p}_{cnt:02d}',x_pixel=round(x,3),y_pixel=round(y,3),cell_area_um2=round(xlo+(x-x0)/(x1-x0)*(xhi-xlo),3),metric=m,value=round(yhi+(y-y0)/(y1-y0)*(ylo-yhi),6),units='um' if 'distance' in m else 'dimensionless',evidence_type='digitized_individual_point',qc_status='isolated_visible_marker; subset_only',reading_uncertainty_px=1.5,cell_area_reading_uncertainty_um2=round(1.5*(xhi-xlo)/(x1-x0),3),value_reading_uncertainty=round(1.5*(yhi-ylo)/(y1-y0),6),calibration_anchor_uncertainty_px=1,cell_area_calibration_uncertainty_um2=round((xhi-xlo)/(x1-x0),3),value_calibration_uncertainty=round((yhi-ylo)/(y1-y0),6),source=str(source.relative_to(STUDY)),source_sha256=sha(source)))
   d.ellipse((x-5,y-5,x+5,y+5),outline='lime',width=1)
  qc.append(dict(figure=f,panel=p,isolated_markers=cnt,rejected_components=reject,expected_n=60,status='partial; overlapping markers and text/axis intersections unresolved; do not refit correlation'))
  cal.append(dict(figure=f,panel=p,source=str(source.relative_to(STUDY)),crop_bounds=[x0,y0,x1,y1],x_axis=[x0,x1,xlo,xhi],y_axis=[y0,y1,yhi,ylo],scale='linear',masked_rectangles=masks,reading_uncertainty_px=1.5,calibration_anchor_uncertainty_px=1))
 im.save(ROOT/f'overlays/{f}_isolated_points.png')
write('size_metric_points.csv',paired);write('scatter_qc.csv',qc)
protocol={'study':'hosny_2026','version':VERSION,'generated_by':'extract.py','python':platform.python_version(),'numpy':np.__version__,'calibration':cal,'input_hashes':{str(p.relative_to(STUDY)):sha(p) for p in list((STUDY/'figures').glob('*.jpg'))+list((STUDY/'derived').glob('hosny_supp_image[34].tiff'))+[STUDY/'source/hosny_preprint.html',STUDY/'source/hosny_supp.docx',STUDY/'source/hosny_preprint_v1.pdf']},'timing':'Published acquisition every 0.5 h over 16 h; distances rescaled by 16/T for shorter tracks, excluding fewer than 3 frames (T<1 h). Source trajectories unavailable.','drift':'No trajectories retrieved; no drift correction inferred.','final_publication':'Metadata/abstract accessible; full text and supplement inaccessible in this session. All readings belong to preprint v1.','uncertainty':'Pixel sensitivity bounds separate from plotted biological IQR; no bootstrap or biological inference.'}
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print('summary rows',len(rows),'isolated scatter markers',len(paired));print(qc)

rayleigh=[(3,'B',[('PC3','normoxia',.673),('PACC','baseline_hypoxia',.0254),('Non-PACC','baseline_hypoxia',.0390)]),(6,'A',[('PACC','siNTC_hypoxia',.00352),('PACC','siHIF_hypoxia',.114),('Non-PACC','siNTC_hypoxia',.0152),('Non-PACC','siHIF_hypoxia',.416)]),(7,'H',[('PC3','siNTC_normoxia',.246),('PC3','siHIF_normoxia',.664)]),(9,'A',[('PACC','siNTC_hypoxia',.0264),('PACC','siRhoA_hypoxia',.286),('Non-PACC','siNTC_hypoxia',.0351),('Non-PACC','siRhoA_hypoxia',.457)]),(10,'H',[('PC3','siNTC_normoxia',.195),('PC3','siRhoA_normoxia',.467)])]
write('rayleigh_tests.csv',[dict(study='hosny_2026',version=VERSION,figure=f,panel=p,group=g,condition=c,p_value=v,n_cells=60,evidence_type='author_reported_numerical_value',source=f'figures/hosny_preprint_fig{f}.jpg') for f,p,groups in rayleigh for g,c,v in groups])
write('reported_values.csv',[dict(figure=f,panel=p,group=g,metric=m,value=v,SEM=sem,units='um2',evidence_type='author_reported_numerical_value',source=source) for f,p,g,m,v,sem,source in [('2','D','Non-PACC','mean_cell_area',487,3,'source/hosny_preprint.html'),('2','D','PACC','mean_cell_area',1490,48,'source/hosny_preprint.html'),('S2','A','normoxic_parent','mean_nuclear_area',165,.64,'source/hosny_supp.docx'),('S2','A','normoxic_parent','mean_cell_area',403,3.2,'source/hosny_supp.docx')]])
