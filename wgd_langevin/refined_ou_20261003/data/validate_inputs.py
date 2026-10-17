"""Check frozen input mapping, scopes, statistic distinctions and curve normalization."""
import hashlib, importlib.util, json, math
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('reviewed_adapter',HERE/'adapter.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
provenance=json.loads((HERE/'source_provenance.json').read_text())
assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest for path,digest in provenance['source_sha256'].items())
assert len(m.OBS)==25 and len(m.CURVE)==510 and len(m.HOSNY)==30
assert [r['id'] for r in m.OBS if r['study']=='Jemaa2017' and r['role']=='fit']==['sarcoma_Cl1','sarcoma_Cl8']
assert [r['id'] for r in m.OBS if r['study']=='Jemaa2023' and r['role']=='fit']==['rko_D1','rko_T1']
assert len(m.observations('Xuan2018'))==2 and all(r['metric']=='speed' for r in m.observations('Xuan2018'))
assert len(m.observations('Xuan2018',True))==2
assert all(r['statistic']=='central_bar_statistic_unspecified' for r in m.OBS if r['study']=='Mallin2023')
assert all(not r['sem_lower'] and not r['sem_upper'] for r in m.OBS if r['study']=='Jemaa2023')
assert all(r['statistic']=='mean' and float(r['sem_lower'])>0 and float(r['sem_upper'])>0 for r in m.OBS if r['id'].startswith('xuan_speed'))
assert all(float(r['low'])<float(r['value'])<float(r['high']) for r in m.OBS)
curves=m.curve_observations()
for g in (0,1):
    part=[r for r in curves if r['group']==g]
    assert all(r['obs']>0 and r['lag_h']>0 and r['low']<r['obs']<r['high'] for r in part)
    assert len(set(r['lag_h'] for r in part))==len(part)
    assert abs(sum(r['weight'] for r in part)-1)<1e-12
    for b in range(3):assert abs(sum(r['weight'] for r in part if r['lag_bin']==b)-1/3)<1e-12
hosny=m.observations('Hosny_baseline');assert len(hosny)==18 and sum(r['role']=='fit' for r in hosny)==6
assert len(m.observations('Hosny_baseline',True))==12
for g in (0,1):
    for metric in ('path_length','directness','parallel','endpoint','perpendicular'):
        rows=[r for r in m.HOSNY if r['group']==('polyploid' if g else 'reference') and r['metric']==metric];q={r['statistic']:float(r['value']) for r in rows}
        assert q['Q1']<=q['median']<=q['Q3']
record={'passed':True,'source_hashes_verified':len(provenance['source_sha256']),'scalar_rows':len(m.OBS),'curve_points':len(curves),'hosny_quantile_rows':len(m.HOSNY),'fitted_scalar_counts':{s:sum(s!='Hosny_baseline' or r['quantile']==1 for r in m.observations(s)) for s in m.STUDIES},'curve_weight_per_group':1,'curve_bin_weight':1/3,'checks':['original clone selections','five-study scope','SEM separate from tolerance','Mallin statistic unspecified','Hosny median mask and IQR','full curve counts and balanced bins','source-table hashes']}
(HERE/'input_validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
