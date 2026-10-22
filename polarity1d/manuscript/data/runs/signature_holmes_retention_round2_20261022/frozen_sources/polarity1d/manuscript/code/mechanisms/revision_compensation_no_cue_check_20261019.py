"""Export existing no-cue regional loss checks; no native reads/integration."""
from pathlib import Path
import csv
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/revision_compensation_excess_baseline_20261019'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/five_model_compensation_20261019/selection'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    (OUT/'supporting/data').mkdir(parents=True,exist_ok=True);(OUT/'supporting/metadata').mkdir(parents=True,exist_ok=True)
    m=json.loads((RUN/'manifest.json').read_text());cs={x['case']['slot']:x['case'] for x in m['cases']};rr=[]
    sources={str((RUN/'manifest.json').relative_to(ROOT)):sha(RUN/'manifest.json'),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}
    for model in m['models']:
        horizon=max(max(c['pulse_durations']) for c in cs.values() if c['model']==model and c['case_kind']=='constant_on')
        p=RUN/'workers'/model/'sampled_loss_events.csv';sources[str(p.relative_to(ROOT))]=sha(p)
        for e in csv.DictReader(p.open()):
            c=cs[int(e['case_slot'])]
            if c['case_kind']!='zero_control' or e['observable']!='excess_fraction' or float(e['loss_fraction']) not in [.25,.5,.75]:continue
            entry=float(e['upper']) if e['upper'] else None
            rr.append(dict(model=model,length=c['length'],case_slot=c['slot'],roi=e['roi'],background_rule=e['background_rule'],loss_fraction=float(e['loss_fraction']),
                matched_native_ON_horizon=horizon,full_no_cue_observation_horizon=float(e['observation_end']),
                status=e['status'],entry_lower=float(e['lower']),entry_upper=entry,
                sampled_loss_entry_within_matched_ON_horizon=bool(entry is not None and entry<=horizon),
                sampled_loss_entry_through_full_10000=bool(entry is not None),
                initial_checkpoint_sha256=c['checkpoint_sha256']))
    assert len(rr)==1746
    path=OUT/'supporting/data/no_cue_localized_excess_loss_check.csv'
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
    proof=dict(controls=97,convention_records=len(rr),matched_horizon_entries=sum(r['sampled_loss_entry_within_matched_ON_horizon'] for r in rr),
        full_10000_entries=sum(r['sampled_loss_entry_through_full_10000'] for r in rr),source_hashes=sources,
        artifact_hashes={str(path.relative_to(OUT)):sha(path)},native_integrations=0)
    (OUT/'supporting/metadata/no_cue_check.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
