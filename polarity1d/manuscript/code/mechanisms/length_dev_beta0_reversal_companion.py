"""Export original-cutoff diagnostics and verify beta=0 development artifacts."""
from pathlib import Path
import argparse,csv,hashlib,json,sys
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(Path(__file__).parent)]
import length_dev_output_layout as layout
DEFAULT=ROOT/'data/runs/debelly_beta0_reversal_20261010'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/beta0_reversal'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def export(run):
    assert read(run/'results/status.json')['complete']
    assert read(run/'results/validation.json')['passed'] and read(run/'results/independent_audit.json')['passed']
    import debelly_beta0_reversal as campaign
    manifest,groups,pending,errors=campaign.r.collect(run)
    assert not pending and not errors
    rows=[]
    for p in sorted((run/'work/preparations').glob('[0-9][0-9][0-9].json')):
        d=read(p);rows.append({k:d.get(k) for k in ['length','original_eligible','verified','eligible','original_cutoff_eligible','asymmetric_eligible','initial_orientation','initial_contrast','elapsed','checkpoint','checkpoint_sha256']})
    diagnostics=[]
    for case in manifest['cases']:
        ts=groups[case['case_index']]
        for gate in ['orientation','full_polarity']:
            native_passes=[t for t in ts if t['duration']>0 and t['original_orientation_success'] and (gate=='orientation' or t['maximum_relative_polarity_error_last_500']<=.01)]
            # The original assay also excludes starts below the original cutoff.
            eligible=case['original_cutoff_eligible']
            selected=min(native_passes,key=lambda t:t['duration']) if native_passes and eligible else None
            diagnostics.append(dict(length=case['length'],condition=case['condition'],definition=gate,original_cutoff_eligible=eligible,original_gate_status='ineligible' if not eligible else 'success' if selected else 'censored',shortest_sampled_original_gate_pulse=selected['duration'] if selected else None,original_gate_bracket_refined=False,note='Original-gate diagnostic on the variant duration grid; only the development gates drove local refinement.'))
    artifacts={}
    for name,data in [('beta0_reversal_preparations.csv',rows),('beta0_reversal_original_cutoff_diagnostics.csv',diagnostics)]:
        path=layout.dev_output_path(OUT,name)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
        artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    checked=0
    for p in (OUT/'supporting/metadata').glob('*.manifest.json'):
        if p.name=='beta0_reversal_companion.manifest.json':continue
        m=read(p)
        for key in ['sha256','figure_sha256','data_sha256']:
            for name,digest in m.get(key,{}).items():
                path=OUT/name;assert path.exists() and sha(path)==digest,(p,name);checked+=1
    links=0
    import re
    for p in [OUT/'README.md',OUT/'supporting/viewers/beta0_reversal_status.html']:
        text=p.read_text()
        targets=re.findall(r'\]\(([^)]+)\)',text) if p.suffix=='.md' else re.findall(r'(?:href|src)="([^"]+)"',text)
        for target in targets:
            if target.startswith(('https:','http:','#')):continue
            assert (p.parent/target).exists(),(p,target);links+=1
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,artifact_hashes_checked=checked,local_links_checked=links,verified_starts=sum(bool(r['verified']) for r in rows),original_cutoff_eligible_starts=sum(bool(r['original_cutoff_eligible']) for r in rows))
    path=layout.dev_output_path(OUT,'beta0_reversal_artifact_audit.json');path.write_text(json.dumps(audit,indent=2)+'\n');artifacts[layout.dev_artifact_key(OUT,path)]=sha(path)
    layout.dev_output_path(OUT,'beta0_reversal_companion.manifest.json').write_text(json.dumps(dict(generated_at=audit['generated_at'],producer='finished',artifact_root='../..',sha256=artifacts,script_sha256=sha(Path(__file__)),campaign_manifest_sha256=sha(run/'manifest.json'),independent_audit_sha256=sha(run/'results/independent_audit.json'),reproduction=f'{sys.executable} {__file__} --run {run}'),indent=2)+'\n')
    print(json.dumps(audit,indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,default=DEFAULT);a=ap.parse_args();export(a.run.resolve())
