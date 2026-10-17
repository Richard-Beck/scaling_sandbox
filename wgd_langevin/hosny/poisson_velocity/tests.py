import numpy as np
from model import simulate,metrics,moments,msd,angle_table

def run_tests():
    checks=[]
    def check(name,actual,expected,tolerance):
        checks.append(dict(test=name,actual=float(actual),expected=float(expected),tolerance=float(tolerance),passed=bool(abs(actual-expected)<=tolerance)))
    tau,lam,A,kappa=.8,1.3,4.,.7
    p,v=simulate(tau,lam,A,kappa,n=32768,duration=4,return_velocity=True,seed=47)
    u,cv=moments(tau,lam,A,kappa)
    d=p[:,-1];sq=(d*d).sum(axis=1)
    check('Stationary mean displacement',d[:,1].mean(),u*4,5*d[:,1].std()/np.sqrt(len(d)))
    check('Analytic MSD',sq.mean(),msd(4,tau,lam,A,kappa),5*sq.std()/np.sqrt(len(sq)))
    vv=((v[:,0]-[0,u])**2).sum(axis=1)
    check('Stationary velocity variance trace',vv.mean(),cv,5*vv.std()/np.sqrt(len(vv)))
    cov=((v[:,0]-[0,u])*(v[:,1]-[0,u])).sum(axis=1)
    check('Velocity covariance at 0.5 h',cov.mean(),cv*np.exp(-.5/tau),5*cov.std()/np.sqrt(len(cov)))
    fine=simulate(tau,lam,A,kappa,n=128,cadence=.1,duration=4,seed=47)
    check('Same events under finer sampling',np.max(abs(fine[:,::5]-p[:128])),0,1e-10)
    scaled=simulate(tau,lam,2*A,kappa,n=128,duration=4,seed=47)
    check('Exact kick amplitude scaling',np.max(abs(scaled-2*p[:128])),0,1e-10)
    m=metrics(p);check('Directional fractions identity',np.max(abs(m['directness']**2-m['parallel']**2-m['perpendicular']**2)),0,1e-12)
    from scipy.special import i1e,i0e
    table=angle_table(kappa);check('Von Mises angle inverse CDF',np.cos(table).mean(),i1e(kappa)/i0e(kappa),1e-4)
    unbiased=simulate(.1,.6,4,0,n=32768,duration=2,seed=93)
    y=unbiased[:,-1,1];check('Uniform kicks have zero mean',y.mean(),0,5*y.std()/np.sqrt(len(y)))
    if not all(x['passed'] for x in checks):raise AssertionError(checks)
    return checks
if __name__=='__main__':print(run_tests())
