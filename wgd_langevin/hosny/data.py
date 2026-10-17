"""Transcription of literature_review/hosny_2026.md only.
Quantiles are [25th, median, 75th]; missing quartiles remain missing.
"""
import numpy as np
Q=(.25,.5,.75)
METRICS=('path','directness','parallel','endpoint','perpendicular')

def row(id,label,group,block,treatment,path,endpoint,directness,parallel,perp,p):
    return dict(id=id,label=label,group=group,block=block,treatment=treatment,
                path=path,endpoint=endpoint,directness=directness,parallel=parallel,perpendicular=perp,rayleigh=p)
DATA=[
 row('O','Normoxia','normoxia','baseline','none',[170,210,285],[20,42,67],[.12,.20,.26],[-.29,0,.17],[-.18,0,.20],.673),
 row('N','Hypoxic non-PACC','non-PACC','baseline','none',[85,105,150],[34,55,95],[.26,.31,.49],[-.17,.10,.39],[-.21,.02,.31],.0390),
 row('P','Hypoxic PACC','PACC','baseline','none',[140,175,223],[60,82,124],[.37,.51,.74],[.05,.25,.50],[-.33,.03,.21],.0254),
 row('HN0','HIF non-PACC control','non-PACC','HIF','control',[86,121,175],[39,60,78],[.25,.36,.50],[-.06,.13,.30],[np.nan,.04,np.nan],.0152),
 row('HN1','HIF non-PACC KD','non-PACC','HIF','KD',[52,81,141],[17,30,50],[.13,.26,.36],[-.14,.04,.17],[np.nan,.02,np.nan],.416),
 row('HP0','HIF PACC control','PACC','HIF','control',[154,198,258],[62,87,120],[.33,.49,.62],[0,.21,.49],[np.nan,.03,np.nan],.00352),
 row('HP1','HIF PACC KD','PACC','HIF','KD',[121,166,207],[43,54,92],[.21,.32,.49],[-.14,.06,.25],[np.nan,-.03,np.nan],.114),
 row('RN0','RhoA non-PACC control','non-PACC','RhoA','control',[111,132,215],[39,63,85],[.22,.38,.52],[-.08,.13,.27],[np.nan,0,np.nan],.0351),
 row('RN1','RhoA non-PACC KD','non-PACC','RhoA','KD',[75,94,133],[16,38,57],[.15,.25,.42],[-.24,.03,.15],[np.nan,0,np.nan],.457),
 row('RP0','RhoA PACC control','PACC','RhoA','control',[148,182,240],[52,87,115],[.27,.51,.65],[0,.22,.36],[np.nan,0,np.nan],.0264),
 row('RP1','RhoA PACC KD','PACC','RhoA','KD',[117,164,210],[32,59,93],[.18,.28,.45],[-.23,.04,.19],[np.nan,0,np.nan],.286),
]
BYID={r['id']:r for r in DATA}
PAIRS=[('N','P','Baseline PACC/non-PACC'),('HN0','HN1','HIF · non-PACC'),('HP0','HP1','HIF · PACC'),('RN0','RN1','RhoA · non-PACC'),('RP0','RP1','RhoA · PACC')]
RADIUS=np.array([2000,4000,4800,5600,6400,7500,8650,11400])
RADIAL_FMI=np.array([.145,.075,.100,.108,.350,.312,.085,.105])
RADIAL_SE=np.array([.065,.070,.062,.067,.052,.066,.073,.074])
