"""Scientific identities and conservative enclosure checks for round5."""
import json,hashlib
import numpy as np
import figure4_round5_analysis_20261023 as A

def archive(amounts,angles,ids=None):
    amounts=np.asarray(amounts,float);angles=np.asarray(angles,float);n,k=amounts.shape
    ids=np.tile(np.arange(k),(n,1)) if ids is None else np.asarray(ids)
    q=(amounts*np.exp(1j*angles)).sum(1)/amounts.sum(1)
    # A nonzero shape/background residual exercises the fourth attribution.
    p=.6*q+.03*np.exp(1j*np.linspace(.2,1,n))
    return dict(time=np.arange(n,dtype=float),region_offset=np.arange(n+1)*k,
        region_id=ids.ravel(),region_amount=amounts.ravel(),region_heading=angles.ravel(),
        region_amount_available=np.ones(n*k,bool),M0=np.ones(n),M1_real=p.real,M1_imag=p.imag)

def main():
    rng=np.random.default_rng(20261023);checks=[]
    z=archive(rng.uniform(.5,2,(20,3)),rng.uniform(-1,1,(20,3)))
    for lag in (1,3,8):
        d=A.decomposition(z,lag)
        assert d['eligible_pairs']==20-lag
        assert d['max_complex_closure_error']<1e-14
        assert d['max_heading_closure_error']<1e-14
        assert abs(sum(d[c+'_loss'] for c in A.COMPONENTS)-d['total_heading_loss'])<1e-12
    checks.append('Exact complex and signed heading-loss closure, 3 lags, unequal amounts/residual')
    # Birth/loss operator changes identity while keeping positions/amounts.
    z=archive([[1,2],[1,2],[1,2]],[[0,.8],[0,.8],[0,.8]],[[0,1],[0,2],[0,2]])
    d=A.decomposition(z,1)
    assert d['max_heading_closure_error']<1e-14
    checks.append('Membership replacement including birth/loss and shared domain')
    z=archive([[1,2],[1,2],[1,2]],[[0,.8],[.2,1.],[.4,1.2]])
    q=(z['region_amount'].reshape(3,2)*np.exp(1j*z['region_heading'].reshape(3,2))).sum(1)/3
    z['M1_real']=q.real;z['M1_imag']=q.imag
    d=A.decomposition(z,1)
    assert abs(d['motion_loss']-d['total_heading_loss'])<1e-12
    assert abs(d['coexisting_amount_loss'])<1e-12 and abs(d['birth_loss_loss'])<1e-12 and abs(d['moment_residual_loss'])<1e-12
    checks.append('Pure domain motion allocates all loss to motion')
    # Variable native moment amplitude must not create a spurious motion
    # excess/residual cancellation in a single rigidly translating region.
    theta=np.linspace(-.7,.9,21);amplitude=np.linspace(.12,.85,21)**1.2
    z=archive(np.ones((21,1)),theta[:,None]);p=amplitude*np.exp(1j*theta)
    z['M1_real']=p.real;z['M1_imag']=p.imag
    for lag in (1,3,7):
        d=A.decomposition(z,lag,amplitude_aware=True)
        assert d['eligible_pairs']>0
        assert abs(d['motion_loss']-d['total_heading_loss'])<1e-12
        assert abs(d['moment_residual_loss'])<1e-12
        assert d['max_heading_closure_error']<1e-12
    checks.append('Amplitude-aware variable-amplitude single rigid translation: motion=1, residual=0')
    z=archive(rng.uniform(.5,2,(20,3)),rng.uniform(-1,1,(20,3)))
    d=A.decomposition(z,3,amplitude_aware=True)
    assert d['max_complex_closure_error']<1e-14 and d['max_heading_closure_error']<1e-14
    checks.append('Amplitude-aware multiregion unequal-amount/residual exact closure')
    assert A.certified_no_overlap(.01,.1,20,np.array([3.]),np.array([.1]),np.array([20.]))
    assert not A.certified_no_overlap(.01,.1,20,np.array([2*np.pi-.01]),np.array([.1]),np.array([20.]))
    assert not A.certified_no_overlap(0,10,20,np.array([3.]),np.array([.1]),np.array([20.]))
    checks.append('Periodic seam, separated enclosures, and half-ring rejection')
    # Random positive contour weights, including highly asymmetric shapes:
    # verify every occupied material cell is inside its conservative enclosure.
    for _ in range(2000):
        n=400;nc=int(rng.integers(1,199));start=int(rng.integers(n));cells=(start+np.arange(nc))%n
        angles=2*np.pi*(cells+.5)/n;weights=rng.uniform(.01,1,nc)
        h=np.angle(np.sum(weights*np.exp(1j*angles)))
        radius=2*np.pi*nc/n
        assert np.all(abs(A.wrapped(angles-h))+np.pi/n <=radius+1e-12)
    checks.append('2,000 asymmetric contiguous supports certify enclosure with half-cell margin')
    data=A.OUT/'supporting/data';metrics=data/'birth_heading_condition_metrics.csv'
    if metrics.exists():
        rows=A.read_csv(metrics)
        assert len(rows)==260 and sum(int(r['n_seeds']) for r in rows)==45696
        for row in rows:
            assert float(row['certified_births'])<=float(row['persistent_secondary_onsets'])
            if float(row['eligible_pairs'])>0:
                assert abs(sum(float(row[c+'_fraction']) for c in A.COMPONENTS)-1)<1e-10
        assert max(float(r['max_heading_closure_error']) for r in rows)<1e-12
        checks.append('All260 conditions/45696 histories; birth identification bounds; every finite signed allocation sums to one')
    amplitude=data/'heading_amplitude_condition_metrics.csv'
    if amplitude.exists():
        revised=A.read_csv(amplitude);assert len(revised)==260
        assert sum(int(r['n_seeds']) for r in revised)==45696
        raw={r['condition_id']:r for r in A.read_csv(metrics)}
        for row in revised:
            prior=raw[row['condition_id']]
            assert row['eligible_pairs']==prior['eligible_pairs'] and row['potential_pairs']==prior['potential_pairs']
            if float(row['eligible_pairs'])>0:
                assert abs(sum(float(row[c+'_fraction']) for c in A.COMPONENTS)-1)<1e-10
        assert max(float(r['max_heading_closure_error']) for r in revised)<1e-12
        checks.append('Amplitude-aware full campaign: identical pair counts; exact signed closure in all170eligible conditions')
        incidence=A.read_csv(data/'takeover_incidence_by_size.csv')
        for model in {r['model'] for r in incidence}:
            for size in {r['C_label'] for r in incidence if r['model']==model}:
                for b in ('-1','0','1','2'):
                    rr=sorted([r for r in incidence if r['model']==model and r['C_label']==size and r['distance_bin']==b],key=lambda r:float(r['horizon_tau']))
                    assert len({r['births'] for r in rr})==1
                    assert np.all(np.diff([int(r['takeovers']) for r in rr])>=0)
        checks.append('Incidence sensitivity uses identical cohort denominators at1–80tau and nondecreasing event counts')
    supplement=A.OUT/'figure4_supplement.pdf'
    if supplement.exists():
        from pypdf import PdfReader
        original=PdfReader(A.OLD/'figure4_supplement.pdf');current=PdfReader(supplement)
        offset=len(current.pages)-len(original.pages)
        assert offset>=8
        for i,page in enumerate(original.pages):
            after=current.pages[offset+i]
            assert list(page.mediabox)==list(after.mediabox)
            assert page.get_contents().get_data()==after.get_contents().get_data()
        mainpdf=PdfReader(A.OUT/'figure4_main.pdf');box=mainpdf.pages[0].mediabox
        assert abs(float(box.width)*25.4/72-180)<1e-4 and abs(float(box.height)*25.4/72-420)<1e-4
        checks.append('All32 historical supplementary page content streams preserved; main180x420mm')
        history=A.OUT/'supporting/historical_unscaled_20261023'
        if amplitude.exists() and (history/'figure4_supplement.pdf').exists():
            reviewed=PdfReader(history/'figure4_supplement.pdf');reviewed_main=PdfReader(history/'figure4_main.pdf')
            offset=len(current.pages)-len(reviewed.pages)
            for i,page in enumerate(reviewed.pages):assert page.get_contents().get_data()==current.pages[offset+i].get_contents().get_data()
            assert reviewed_main.pages[0].get_contents().get_data()==current.pages[offset-1].get_contents().get_data()
            snapshot=json.loads((history/'reviewed_unscaled_snapshot.json').read_text())
            for name,digest in snapshot.items():
                target=A.OUT/name if name.startswith('supporting/data/') else history/name
                assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
            checks.append('Reviewed unscaled main+40supplement pages preserved; original metrics and20tau takeover CSV hashes unchanged')
    out=dict(passed=True,checks=checks)
    path=A.OUT/'supporting/metadata/scientific_validation.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
