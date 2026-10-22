from copy import deepcopy
from pathlib import Path
import sys
import json

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure4_regions_20261020 import track_profiles
from signature_figure4_regions_audit_20261020 import audit_case, frame_checks, takeover_checks
from signature_figure4_regions_campaign_20261020 import analyze_case, sha


def example():
    q=(np.arange(80)+.5)/80
    p=.2+np.exp(-.5*(((q-.2+.5) % 1-.5)/.07)**2)
    p=np.tile(p,(3,1))
    result=track_profiles(np.arange(3.),p,20.,basal_contrast=1.)
    return p,result["frames"]


def test_independent_frame_audit_accepts_complete_activity_and_physical_widths():
    p,frames=example()
    assert frame_checks(p,20.,1.,frames)["frames"]==3


@pytest.mark.parametrize("mutation",["count","width","amount","contrast"])
def test_independent_frame_audit_rejects_measurement_and_roi_errors(mutation):
    p,frames=example();frames=deepcopy(frames)
    if mutation=="count":frames[1]["region_count"]+=1
    if mutation=="width":frames[1]["regions"][0]["width"]*=.5
    if mutation=="amount":frames[1]["regions"][0]["excess"]*=.5
    if mutation=="contrast":frames[1]["contrast"]*=.5
    with pytest.raises(ValueError):frame_checks(p,20.,1.,frames)


def test_takeover_audit_rejects_pairwise_win_without_global_observed_dominance():
    event=dict(time=2.,first_qualifying_frame_time=1.,newborn_track_id=1,
               previous_dominant_excess=1.,newborn_excess=4.,globally_dominant_required=True)
    result=dict(frames=[dict(time=1.,dominant_track_id=1),dict(time=2.,dominant_track_id=1)],takeovers=[event])
    takeover_checks(result)
    result["frames"][1]["dominant_track_id"]=2
    with pytest.raises(ValueError):takeover_checks(result)


def test_native_adapter_to_all_view_audit_checks_identity_and_every_imaging_variant(tmp_path):
    native=tmp_path/"native";native.mkdir()
    p,_=example();basal=p[0];time=np.arange(0.,401.,40.)
    np.savez(native/"profiles.npz",time=time,primary=np.tile(basal,(len(time),1)),
             circumference=20.,material_positions=(np.arange(len(basal))+.5)/len(basal))
    np.savez(native/"basal.npz",primary=basal,circumference=20.)
    np.savez(native/"checkpoint.npz",state=np.array([basal,np.ones_like(basal)]))
    for name in ("scalar.json","status.json"):(native/name).write_text("{}\n")
    case=dict(id="example",model="wave_pinning",C_label=20.,native_cells=len(basal),
              profile_cadence=40.,native_horizon=400.,actualC0=20.,basal_contrast=float(np.ptp(basal)))
    for key,name in (("profile","profiles.npz"),("scalar","scalar.json"),("basal_profile","basal.npz"),
                     ("basal_checkpoint","checkpoint.npz"),("source_status","status.json")):
        case[key+"_path"]=str(native/name);case[key+"_sha256"]=sha(native/name)
    workers=tmp_path/"workers"
    analyze_case(case,workers/case["id"],inventory_sha256="test_inventory")
    result=audit_case(case,workers,"test_inventory")
    assert result["views_checked"]==24 and result["all_view_native_frames"]>result["frames"]
    assert result["measurement_sha256"]==sha(workers/case["id"]/"result.json")
    receipt=workers/case["id"]/"result.json"
    missing=json.loads(receipt.read_text());missing["sensitivities"].pop()
    receipt.write_text(json.dumps(missing))
    with pytest.raises(ValueError):audit_case(case,workers,"test_inventory")
