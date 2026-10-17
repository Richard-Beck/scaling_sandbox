"""Freeze reviewed scalar readings and full published MSD curves; never fit movie tracks."""
import csv, hashlib, json, math
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'literature_review/migration_data'
def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
def num(row,key,default=0.):return float(row[key]) if row.get(key) else default
old=read(ROOT/'wgd_langevin/data/observations.csv')
xuan=read(SOURCE/'xuan_2018/extraction/digitized_points.csv')
curves=[]
for i,r in enumerate(xuan):
    if r['point_kind']!='plotted_curve_sample':continue
    value=num(r,'y_value');unc=num(r,'value_uncertainty');lag=num(r,'x_value')/3600
    curves.append(dict(id=f'xuan_curve_{i:03}',study='Xuan2018',group='reference' if r['condition']=='Non-PGCC' else 'polyploid',metric='msd',lag_h=lag,value=value,low=value-unc,high=value+unc,reading_uncertainty=value*(10**(1.5*4/416)-1),calibration_uncertainty=value*(10**(1*4/416)-1),combined_raster_uncertainty=unc,pixel_reading_uncertainty=1.5,axis_calibration_uncertainty_px=1,lag_uncertainty_h=lag*(10**(math.hypot(1.5,1)/141)-1),log_scale=math.log(3000/2000)/2 if r['condition']=='Non-PGCC' else math.log(5000/4000)/2,statistic='published_population_MSD_curve',x_px=r['x_px'],y_px=r['y_px'],source='literature_review/migration_data/xuan_2018/extraction/digitized_points.csv',source_sha256=r['source_sha256'],role='fit',uncertainty_definition='raster reading/calibration sensitivity, not biological error',scale_definition='legacy late-MSD half-log tolerance; balanced log-lag curve contribution'))
for group in ('reference','polyploid'):
    part=[r for r in curves if r['group']==group];a=min(math.log(r['lag_h']) for r in part);b=max(math.log(r['lag_h']) for r in part)
    for r in part:r['lag_bin']=min(2,int(3*(math.log(r['lag_h'])-a)/(b-a)))
    counts={j:sum(r['lag_bin']==j for r in part) for j in range(3)}
    for r in part:r['objective_weight']=1/(3*counts[r['lag_bin']])
write('xuan_msd_curve.csv',curves)
scalars=[]
for original in old:
    r=original.copy();study=r['study'];metric=r['metric'];group=r['group'];value=num(r,'value');unc=cal=0.;semlo=semhi='';stat='legacy_approximation';src=r['source'];sha='';notes=r['notes']
    if study=='Xuan2018' and metric in ('speed','msd_slope','slow_fraction'):
        series={'speed':'average_speed','msd_slope':'MSD_logarithmic_slope','slow_fraction':'slow_mode_fraction'}[metric]
        points=[p for p in xuan if p['series']==series and p['condition']==('Non-PGCC' if group=='reference' else 'PGCC')]
        mean=next(p for p in points if p['point_kind']=='mean');value=num(mean,'y_value');unc=num(mean,'value_uncertainty')*1.5/1.803;cal=num(mean,'value_uncertainty')/1.803;stat='mean';sha=mean['source_sha256'];src='literature_review/migration_data/xuan_2018/extraction/digitized_points.csv'
        semlo=value-num(next(p for p in points if p['point_kind']=='SEM_lower'),'y_value');semhi=num(next(p for p in points if p['point_kind']=='SEM_upper'),'y_value')-value
    elif study=='Xuan2018' and metric=='msd':
        part=sorted((p for p in curves if p['group']==group),key=lambda p:p['lag_h']);t=num(r,'lag_h');point=min(part,key=lambda p:abs(math.log(p['lag_h']/t)));value=point['value'];unc=point['combined_raster_uncertainty']*1.5/1.803;cal=point['combined_raster_uncertainty']/1.803;stat='published_population_MSD_curve';src=point['source'];sha=point['source_sha256']
        if r['role']=='fit':r['role']='replaced';notes+='; late scalar replaced by balanced full population curve, not additionally scored'
        else:
            # Interpolate in plotted log-log coordinates to the same original 600-second check.
            for p,q in zip(part[:-1],part[1:]):
                if p['lag_h']<=t<=q['lag_h']:
                    f=math.log(t/p['lag_h'])/math.log(q['lag_h']/p['lag_h']);value=math.exp(math.log(p['value'])*(1-f)+math.log(q['value'])*f);combined=value*(p['combined_raster_uncertainty']/p['value']*(1-f)+q['combined_raster_uncertainty']/q['value']*f);unc=combined*1.5/1.803;cal=combined/1.803;break
            notes+='; log-log interpolated at original lag; overlaps fitted curve, correlated check, not independent holdout'
    elif study=='Mallin2023':
        table=read(SOURCE/'mallin_2023/extraction/measurements.csv');mapped={'path_length':'accumulated_distance','endpoint':'euclidean_distance','directness':'directness'}[metric]
        p=next(p for p in table if p['condition']==('non_PACC' if group=='reference' else 'PACC') and p['metric']==mapped);value=num(p,'value');unc=num(p,'reading_uncertainty_value');cal=num(p,'calibration_uncertainty_value');stat='central_bar_statistic_unspecified';src='literature_review/migration_data/mallin_2023/extraction/measurements.csv';sha=p['source_sha256'];notes+='; statistic remains unspecified, legacy model mean readout is an approximation'
    elif study in ('Jemaa2017','Jemaa2023'):
        folder='jemaa_2017' if study=='Jemaa2017' else 'jemaa_2023';clone=r['id'].split('_')[-1];p=next(p for p in read(SOURCE/folder/'extraction/measurements.csv') if p['clone']==clone);value=num(p,'value');unc=num(p,'reading_uncertainty_units');cal=num(p,'calibration_uncertainty_units');semlo=p['sem_lower'];semhi=p['sem_upper'];stat='mean';src=f'literature_review/migration_data/{folder}/extraction/measurements.csv';notes+='; '+p['sem_status']
    # Recenter old descriptive tolerance preserving its two relative endpoints, without claiming raster bounds are biological uncertainty.
    raster=math.hypot(unc,cal) if study=='Xuan2018' else unc+cal
    legacy=num(original,'value');r.update(value=value,low=value*num(original,'low')/legacy,high=value*num(original,'high')/legacy,statistic=stat,reading_uncertainty=unc,calibration_uncertainty=cal,combined_raster_uncertainty=raster,raster_low=value-raster,raster_high=value+raster,sem_lower=semlo,sem_upper=semhi,scale_definition='legacy relative heuristic tolerance recentered on reviewed reading; not SEM/CI',source=src,source_sha256=sha,notes=notes)
    scalars.append(r)
write('observations.csv',scalars)
hosny=[]
metricmap={'accumulated_distance':'path_length','euclidean_distance':'endpoint','directness':'directness','parallel_FMI':'parallel','perpendicular_FMI':'perpendicular'}
for p in read(SOURCE/'hosny_2026/extraction/measurements.csv'):
    if p['condition']!='baseline_hypoxia':continue
    if p['metric'] not in metricmap:raise ValueError(p['metric'])
    metric=metricmap[p['metric']];p=dict(p,id='P' if p['group']=='PACC' else 'N',study='Hosny_baseline',group='polyploid' if p['group']=='PACC' else 'reference',metric=metric,source='literature_review/migration_data/hosny_2026/extraction/measurements.csv',role='fit' if metric in ('path_length','directness','parallel') and p['statistic']=='median' else 'check')
    hosny.append(p)
write('hosny_observations.csv',hosny)
paths=[ROOT/'wgd_langevin/data/observations.csv',ROOT/'wgd_langevin/hosny/data.py',SOURCE/'OU_FITTING_DISPOSITION.md']+[SOURCE/f'{s}/extraction/{file}' for s,file in [('xuan_2018','digitized_points.csv'),('mallin_2023','measurements.csv'),('jemaa_2017','measurements.csv'),('jemaa_2023','measurements.csv'),('hosny_2026','measurements.csv')]]
(HERE/'source_provenance.json').write_text(json.dumps({'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'output_rows':{'observations.csv':len(scalars),'xuan_msd_curve.csv':len(curves),'hosny_observations.csv':len(hosny)},'campaign_studies':['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023','Hosny_baseline']},indent=2)+'\n')
assert len(curves)==510
assert len(hosny)==30
assert all(abs(sum(p['objective_weight'] for p in curves if p['group']==g)-1)<1e-12 for g in ('reference','polyploid'))
print('Frozen inputs:',len(scalars),'scalar/context rows;',len(curves),'curve points;',len(hosny),'Hosny baseline quantiles')
