"""Untuned transfer to the existing De Belly Figure 2 establishment/countercues."""
import sys
import numpy as np
from numba import njit
from design import *
from reduced import drift
from scaling import parameters
sys.path.insert(0,str(ROOT/'wgd_langevin/logistic_feedback/debelly_matching'))
from campaign import load_data

@njit(cache=True)
def forcing(t,pulses):
    u=0.
    for start,stop,taper,amplitude in pulses:
        if start<=t<stop:
            w=1.
            if taper>0 and t>stop-taper:w=.5*(1+np.cos(np.pi*(t-stop+taper)/taper))
            u+=amplitude*w
    return u

@njit(cache=True)
def simulate(theta,times,pulses,family,refine=1.):
    out=np.zeros((len(times),2));p=q=0.;t=times[0]
    cap=min(.1,.25*theta[0],.25*theta[1],.25/(theta[6]+1e-6))/refine
    for i in range(1,len(times)):
        while t<times[i]-1e-10:
            h=min(cap,times[i]-t)
            for start,stop,taper,amplitude in pulses:
                for event in (start,stop-taper,stop):
                    if event>t+1e-10:h=min(h,event-t)
            a,b=drift(p,q,theta,family,forcing(t+1e-11,pulses))
            c,d=drift(p+h*a/2,q+h*b/2,theta,family,forcing(t+h/2,pulses))
            e,f=drift(p+h*c/2,q+h*d/2,theta,family,forcing(t+h/2,pulses))
            g,j=drift(p+h*e,q+h*f,theta,family,forcing(t+h-1e-11,pulses))
            p+=h*(a+2*c+2*e+g)/6;q+=h*(b+2*d+2*f+j)/6;t+=h
        out[i,0],out[i,1]=p,q
    return out

def run():
    cases,_,paths=load_data();rows=[]
    for family in FAMILIES:
        source=OUT/'scaling'/f'debelly_{family}.json'
        record=json.loads(source.read_text());law=next(x for x in record['laws'] if x['law']==record['selected_law'])
        for c in cases:
            if c['length'] not in lengths('debelly'):continue
            if (c['panel']=='A' and c['variant']!='full') or (c['panel']=='C' and c['variant']!='beta_1x'):continue
            theta=parameters(law,c['length']);pulses=c['pulses'].copy();pulses[:,:3]/=1000
            pred=simulate(theta,c['times']/1000,pulses,FAMILIES.index(family))[:,0]
            y=c['values'][:,0]
            tail=c['times']>=c['times'][-1]-1000
            final=float(np.median(pred[tail]));actual=float(np.median(y[tail]))
            rows.append(dict(family=family,length=c['length'],panel=c['panel'],task_id=c['task_id'],
                             rmse=float(np.sqrt(np.mean((pred-y)**2))),final_actual=actual,final_predicted=final,
                             final_class_agrees=bool((abs(final)>=.02)==(abs(actual)>=.02) and
                                                    (abs(actual)<.02 or np.sign(final)==np.sign(actual))),
                             physical_range=bool(np.max(np.abs(pred))<=.5)))
    save_json(OUT/'retained_transfer.json',dict(rows=rows,
              note='No fitting or cue-gain adjustment. Rac only; latent q is not equated to MCA.',
              source_hashes=digest([Path(__file__)]+paths)))

if __name__=='__main__':run()
