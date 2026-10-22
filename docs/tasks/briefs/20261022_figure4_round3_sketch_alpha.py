import csv,math,collections,json,sys
D=sys.argv[1]
rows=[r for r in csv.DictReader(open(D+'/condition_lag_metrics.csv')) if r['arm']=='common' and r['view']=='main']
by=collections.defaultdict(dict)
for r in rows:
    try: C=float(r['C']); cov=float(r['coverage'])
    except: continue
    by[(r['model'],float(r['f']))].setdefault(float(r['C_label']),{})[float(r['lag'])]=(C,cov,r['forced_combined'],r['regime'])
out=[]
for (m,f),sz in sorted(by.items()):
    sizes=sorted(sz)
    if len(sizes)<2: continue
    for i,s in enumerate(sizes):
        nb=sizes[max(0,i-1):i+2]
        lags=set.intersection(*[set(sz[x]) for x in nb])
        best=None
        for L in lags:
            vals=[sz[x][L] for x in nb]
            if any(not(0.03<v[0]<0.995) or v[1]<0.5 for v in vals): continue
            score=max(abs(math.log(-math.log(v[0]))-math.log(-math.log(0.5))) for v in vals)
            if best is None or score<best[0]: best=(score,L)
        if best is None: continue
        L=best[1]
        xs=[math.log(x) for x in nb]; ys=[math.log(-math.log(sz[x][L][0])/L) for x in nb]
        mx=sum(xs)/len(xs); my=sum(ys)/len(ys)
        sl=sum((a-mx)*(b-my) for a,b in zip(xs,ys))/sum((a-mx)**2 for a in xs)
        forced=any(sz[x][L][2]=='True' for x in nb)
        out.append(dict(model=m,f=f,size=s,alpha=-sl,lag=L,n=len(nb),forced=forced,regime=sz[s][L][3]))
json.dump(out,open('alpha_data.json','w'))
for o in out:
    if o['model']!='spring_ring' or True: pass
import itertools
for m in ['wave_pinning','holmes_model3','spring_ring','goryachev']:
    print(m)
    for o in out:
        if o['model']==m: print('  f=%g C=%g alpha=%.2f lag=%g n=%d forced=%s %s'%(o['f'],o['size'],o['alpha'],o['lag'],o['n'],o['forced'],o['regime']))
