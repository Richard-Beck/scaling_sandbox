"""Freeze remaining strength series with exact accepted high-dose native reuses.

The native worker/engine are unchanged. This initializer only proves source and
protocol identities and removes already-executed groups from the submit array.
"""
from pathlib import Path
import argparse
import copy
import shutil
import revision_five_model_compensation_campaign_20261019 as campaign
import revision_five_model_compensation_20261019 as p

ROOT=Path(__file__).resolve().parents[2]
SCREEN=ROOT/'data/runs/revision_five_model_high_strength_screen_20261019'
DEFAULT=ROOT/'data/runs/revision_five_model_compensation_strength_v2_20261019'


def screen_reuses():
    manifest=p.verify(SCREEN); provenance=p.read(SCREEN/'provenance.json')
    for rel,h in p.n.sources().items():assert provenance['frozen_sources'][rel]==h,rel
    engine='code/mechanisms/revision_five_model_compensation_20261019.py'
    assert provenance['frozen_sources'][engine]==p.n.sha(ROOT/engine)
    reused={}
    for c in manifest['cases']:
        folder=SCREEN/'workers'/f'{c["slot"]:03d}'
        r=p.read(folder/'result.json');screen=p.read(folder/'numerical_screen.json')
        assert r['case']==c and r['manifest_sha256']==p.n.sha(SCREEN/'manifest.json')
        assert screen['status'] in ['passes','native_high_dose_restart_checks_pass_no_step_convergence_claim']
        assert screen['primary_native_result_sha256']==p.n.sha(folder/'result.json')
        assert p.n.sha(folder/'cue.npz')==r['cue_archive_sha256']
        proof=dict(source_manifest=str((SCREEN/'manifest.json').relative_to(ROOT)),source_manifest_sha256=p.n.sha(SCREEN/'manifest.json'),
            source_result=str((folder/'result.json').relative_to(ROOT)),source_result_sha256=p.n.sha(folder/'result.json'),
            numerical_screen=str((folder/'numerical_screen.json').relative_to(ROOT)),numerical_screen_sha256=p.n.sha(folder/'numerical_screen.json'),
            numerical_screen_status=screen['status'],cue_archive=str((folder/'cue.npz').relative_to(ROOT)),cue_archive_sha256=r['cue_archive_sha256'],withdrawals=[])
        for w in r['withdrawals']:
            if w['duration']==0:continue
            item=copy.deepcopy(w)
            for name,hname in [('archive','archive_sha256'),('restart','restart_sha256')]:
                path=SCREEN/w[name];assert p.n.sha(path)==w[hname];item[name]=str(path.relative_to(ROOT))
            proof['withdrawals'].append(item)
        assert [w['duration'] for w in proof['withdrawals']]==p.PULSES[c['model']]
        reused[campaign.scientific_key(c)]=(c,proof)
    assert len(reused)==10
    return reused


def initialize(out):
    reuses=screen_reuses()
    campaign.initialize(out,'strength')
    m=p.read(out/'manifest.json');provenance=p.read(out/'provenance.json');matched=[]
    for c in m['cases']:
        pair=reuses.get(campaign.scientific_key(c))
        if not pair:continue
        source,proof=pair
        for field in ['model','length','cells','checkpoint_sha256','response_amplitude','width','fixed_scales','background_stimulus']:
            assert c.get(field)==source.get(field),(c['slot'],field)
        assert source['pulse_durations']==c['pulse_durations']
        if c['model']=='debelly':assert source['center_boundary_motion']==c['center_boundary_motion']==False
        with p.np.load(out/c['regional_input']) as current, p.np.load(SCREEN/source['regional_input']) as original:
            for name in ['initial_active','initial_core_mask','initial_weights','initial_cue_footprint','initial_halfcontrast_threshold']:
                p.np.testing.assert_array_equal(current[name],original[name])
        c['exact_reuse']=proof;matched.append(c['slot'])
    assert len(matched)==10
    jobs=[dict(task=i,case_slot=c['slot']) for i,c in enumerate([c for c in m['cases'] if not c['exact_reuse']])]
    m.update(jobs=jobs,exact_reuse_case_groups=10,new_native_workers=len(jobs),
        reused_high_strength_case_slots=matched,reused_logical_protocol_rows=50,
        high_strength_screen_manifest_sha256=p.n.sha(SCREEN/'manifest.json'),
        strength_levels=[.03,.1,.3,3.],baseline_one_times_in_separate_baseline=True)
    assert len(jobs)==378 and m['logical_protocol_rows']==1940
    rel=str(Path(__file__).relative_to(ROOT));provenance['frozen_sources'][rel]=p.n.sha(Path(__file__))
    target=out/'frozen_sources'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(Path(__file__),target)
    p.write(out/'manifest.json',m);provenance['manifest_sha256']=p.n.sha(out/'manifest.json');p.write(out/'provenance.json',provenance)
    p.write(out/'results/latest_status.json',dict(generated_at=p.n.now(),producer='frozen_not_submitted',complete=False,
        total_case_groups=388,total_native_workers=378,logical_protocol_rows=1940,exact_high_strength_reuse_groups=10))
    command='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 --array=0-377 --job-name=F2-compensation-strength --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(out/'worker.sh')
    (out/'submission_command.txt').write_text(command+'\n')
    print(dict(out=str(out),manifest_sha256=p.n.sha(out/'manifest.json'),new_workers=378,reused_groups=10,protocol_rows=1940))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=DEFAULT)
    initialize(parser.parse_args().out.resolve())
