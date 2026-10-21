"""Review tables from the completed E0 analysis; never starts new experiments."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as s:
        for b in iter(lambda:s.read(4<<20),b''):h.update(b)
    return h.hexdigest()


def load(path):return json.loads(Path(path).read_text())
def table(path):return list(csv.DictReader(Path(path).open()))


def write_csv(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as s:
        w=csv.DictWriter(s,fieldnames=fields);w.writeheader();w.writerows(rows)


def main(root):
    result=root/'results/analysis';gate=load(result/'status.json');audit=load(result/'independent_estimand_audit.json')
    if not(gate['passed'] and gate['complete'] and audit['passed'] and audit['condition_groups']==110):
        raise ValueError('Complete E0 and independent estimand audit are required')
    for name,digest in gate['outputs'].items():
        if sha(result/name)!=digest:raise ValueError('E0 CSV changed')
    manifest=load(root/'analysis_manifest.json')
    repository=next(p for p in root.parents if (p/'docs/tasks/briefs/20261020_figure4_e1_reviewed_selection.json').is_file())
    selection=repository/'docs/tasks/briefs/20261020_figure4_e1_reviewed_selection.json'
    decisions=[r for r in table(result/'decisions.csv') if r['view']=='main']
    e4=[]
    for g in manifest['groups']:
        worker=root/f"analysis_workers/{g['task']:04d}";r=load(worker/'result.json');m=r['metrics']['main']
        if sha(worker/'arrays.npz')!=r['arrays_sha256']:raise ValueError('Worker arrays changed')
        with np.load(worker/'arrays.npz',allow_pickle=False) as z:
            cov=z['main__coverage'];C=z['main__C'];lag=z['main__lag']
            usable=np.flatnonzero((cov>=.5)&np.isfinite(C))
            if len(usable):
                j=int(usable[-1]);status='trigger' if C[j]>.8 else 'not_triggered'
                details=dict(last_usable_lag=float(lag[j]),last_usable_C=float(C[j]),last_usable_coverage=float(cov[j]),
                    last_usable_C_ci_lower=float(z['main__C_ci_lower'][j]),last_usable_C_ci_upper=float(z['main__C_ci_upper'][j]))
            else:
                status='no_lag_meets_coverage';details=dict(last_usable_lag='',last_usable_C='',last_usable_coverage='',last_usable_C_ci_lower='',last_usable_C_ci_upper='')
            max_C_coverage=float(np.max(cov))
        e4.append(dict(model=g['model'],C_label=g['C_label'],f=g['f'],status=status,
            trigger=status=='trigger',forced=m['forced'],existing_noisy_seeds=len(g['noise']),
            proposed_E4_seeds=12 if status=='trigger' else 0,original_horizon=g['noise'][0]['native_horizon'],
            proposed_horizon=4*g['noise'][0]['native_horizon'] if status=='trigger' else '',
            maximum_lag_coverage=max_C_coverage,tau_status=m['tau_status'],**details))
    recommendations=[
        dict(model='wave_pinning',transition_f=.03,strong_f=.1,
             evidence='At .03, C20 tau=1849.86 while C60 is censored; C60 multiple-front fraction=.23. At .1 it is .95, with non-forced means.',
             alternative_transition_f=.01,alternative_strong_f=.1),
        dict(model='holmes_model3',transition_f=.03,strong_f=.1,
             evidence='At .03, C60 multiple-front fraction=.72 and persistence ratio=2.77; at .1 fraction=.92 and ratio=2.38. Matched means are non-forced.',
             alternative_transition_f=.01,alternative_strong_f=.1),
        dict(model='spring_ring',transition_f=.03,strong_f=.1,
             evidence='At .03 the C3 multiple-front fraction is close to .5; at .1 it is .79. Both means are non-forced; .3 and 1 are forced.',
             alternative_transition_f=.01,alternative_strong_f=.1),
        dict(model='goryachev',transition_f=.02,strong_f=.03,
             evidence='Root-reviewed scientific bracket: C60 stays polarized through .02 and first becomes low-contrast dominated at .03 (fraction=.86). E2 .02 sigma-one small/mid companions supply missing references; existing C60 refinement seeds202610195800100..111 must pair exactly. The .03 mean is forced, so E1 tests fluctuation versus mean.',
             alternative_transition_f='',alternative_strong_f='')
    ]
    for row in recommendations:
        row['status']='root_reviewed_selection_pending_user_interim_delivery'
    paths={'provisional_main_decisions.csv':decisions,'e4_trigger_inventory.csv':e4,'e1_recommendations.csv':recommendations}
    out=root/'results/interim_handoff'
    out.mkdir(exist_ok=False)
    for name,rows in paths.items():write_csv(out/name,rows)
    trigger=[r for r in e4 if r['trigger']]
    memo=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,scope='E0 review handoff, no E1/E2/E4 execution',
        E0_status_sha256=sha(result/'status.json'),independent_estimand_audit_sha256=sha(result/'independent_estimand_audit.json'),
        source_sha256=sha(__file__),condition_groups=110,E4_trigger_groups=len(trigger),E4_no_usable_lag_groups=sum(r['status']=='no_lag_meets_coverage' for r in e4),
        reviewed_E1_selection_path=str(selection),reviewed_E1_selection_sha256=sha(selection),
        E4_requested_native_integrations=12*len(trigger),
        E4_triggers_by_model={model:sum(r['model']==model for r in trigger) for model in sorted({r['model'] for r in e4})},
        E4_caveat='Every trigger is shown, including forced conditions. The brief count of 100-200 runs is not a cap. A fourfold horizon does not guarantee uncensoring at the weakest inputs.',
        provisional_H1='Inconclusive in all four models; weakest-input persistence remains censored, fhalf intervals overlap, and empirical persistence ratios at weak input are not identified.',
        provisional_H2='Holmes satisfies the model peak-drop criterion in the main view, but it becomes inconclusive in the second-half-of-analysis and late-quarter windows. No non-forced evaluable condition has a ratio CI wholly below one. Overall H2 remains pending A6 and evidence in a second model.',
        physical_scope='Simulated polarity heading and qualified site mobility, not observed migration or an assumed cell-velocity law.',
        alternatives=[
            dict(framing='Model-dependent attenuation of the size advantage',support='Holmes ratio falls from 2.77 at f=.03 to1.24 at .3 with large-cell multiplicity above .5; Spring ratios also decrease across .01-.1.',weakness='No overall H2 decision; weak-input ratios and fhalf remain censored.',next='E2 and A6; retain full regime maps.'),
            dict(framing='Enlargement increases persistence over a measured intermediate-input range',support='Non-forced ratio CIs exceed one in WP, Holmes and Spring at several strengths.',weakness='H1 is not established at every evaluable f and fhalf cannot be resolved from weak baselines.',next='E4 plus adequate weak-input replication.'),
            dict(framing='Multiplicity versus contrast loss under stronger input',support='Large WP/Holmes/Spring become multiple-front dominated; large Gory becomes low contrast.',weakness='Gory strong-input means are forced, so H3 requires E1 and A5 fingerprints.',next='Matched-mean E1 sigma sweep and optical A5/A6.')],
        outputs={name:sha(out/name) for name in paths})
    (out/'memo.json').write_text(json.dumps(memo,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in memo.items() if k not in ('alternatives','outputs')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();main(args.out.resolve())
