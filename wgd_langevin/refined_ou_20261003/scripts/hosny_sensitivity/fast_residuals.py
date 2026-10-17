"""Equivalent compiled residual vectors for sensitivity objectives only.

Metadata rows retain observation identity and masks; current predictions are
returned only through the residual vector. Never use this adapter for fit-row
exports. Frozen fit/model scripts remain unchanged.
"""
import numpy as np

def compile_residuals(fit):
    original=fit.data_residuals
    records=[]
    for study in fit.STUDIES:
        obs=fit.FIT_OBS[study]
        if study=='Xuan2018':obs=obs+fit.CURVE_OBS
        records.extend([dict(study=study,**o) for o in obs])
    observed=np.array([o['obs'] for o in records])
    is_hosny=np.array([o['study']=='Hosny_baseline' for o in records])
    scale=np.array([o['width'] if o['study']=='Hosny_baseline' else o.get('log_scale',np.log(o['high']/o['low'])/2) for o in records])
    weight=np.sqrt([o.get('objective_weight',1.) for o in records])
    curves={g:np.array([i for i,o in enumerate(records) if o['metric']=='msd_curve' and o['group']==g],int) for g in (0,1)}
    curve_lags={g:np.array([records[i]['lag_h'] for i in curves[g]]) for g in (0,1)}
    scalar=[i for i,o in enumerate(records) if o['metric']!='msd_curve']
    def residuals(pred,checks=False):
        if checks:return original(pred,checks=True)
        value=np.empty(len(records))
        for i in scalar:
            o=records[i];p=pred[o['study']][o['group']]
            value[i]=p['quantiles'][o['metric']][o['quantile']] if 'quantile' in o else p[o['metric']]
        for g,ids in curves.items():
            p=pred['Xuan2018'][g];value[ids]=np.interp(curve_lags[g],p['curve_lags_h'],p['msd_curve'])
        r=np.empty_like(value)
        r[is_hosny]=(value[is_hosny]-observed[is_hosny])/scale[is_hosny]
        r[~is_hosny]=np.log(np.maximum(value[~is_hosny],1e-300)/observed[~is_hosny])/scale[~is_hosny]
        return r*weight,records
    return residuals
