"""Explicit Spring counterpart to the historical chemical-only descriptor test."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
RUN=ROOT/"data/runs/revision_spring_distance_duration_20261019"
DATA=TOPIC/"spring_distance/supporting/data"
META=TOPIC/"spring_distance/supporting/metadata"


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    source=DATA/"spring_distance_duration_summary.csv";frame=pd.read_csv(source)
    long=frame[np.isclose(frame.cue_duration,1000.)].copy();assert len(long)==14
    rows=[]
    for _,r in long.iterrows():
        path=Path(r.actual_recovery_path);assert sha(path)==r.actual_recovery_sha256
        with np.load(path) as z:y=z["states"][-1];q=z["material_positions"];n=len(q)
        peak=int(np.argmax(y[:n]));gap=float(abs((q[peak]-r.cue_center_fraction+.5)%1-.5)*y[-1])
        associated=gap<=.1+y[-1]/(2*n)
        assert associated and r.resolved_peaks_final==1 and r.half_contrast_components_final==1 and r.cue_aligned_final_15deg
        rows.append(dict(resting_circumference=r.resting_circumference,basal_circumference=r.basal_circumference,
            fraction=r.fraction,physical_distance_at_onset=r.physical_distance_at_onset,prescribed_normal_angle_deg=r.normal_angle_deg,
            actual_peak_to_cue_distance_final=gap,final_physical_circumference=r.physical_C_final,
            final_activity_angle_deg=r.angle_final,final_R=r.R_final,final_q1=r.q1_final,final_q2=r.q2_final,
            final_original_core_excess_fraction=r.old_core_excess_final,resolved_peaks=1,half_contrast_regions=1,
            finite_fate_category="cue-associated single peak",full_field_tail_stationary=bool(r.tail_stationary),
            actual_recovery_path=str(path),actual_recovery_sha256=r.actual_recovery_sha256))
    out=pd.DataFrame(rows);table=DATA/"spring_standard_long_cue_descriptor_coverage.csv";out.to_csv(table,index=False)
    ranges=[]
    for C0,g in out.groupby("resting_circumference"):
        ranges.append(dict(resting_circumference=float(C0),eligible_geometries=len(g),finite_fate_classes=g.finite_fate_category.nunique(),
            final_angle_min=float(g.final_activity_angle_deg.min()),final_angle_max=float(g.final_activity_angle_deg.max()),
            original_core_excess_min=float(g.final_original_core_excess_fraction.min()),original_core_excess_max=float(g.final_original_core_excess_fraction.max())))
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,model="Spring mean-force global-extension ring analogue",
        standard_cue_duration=1000.,cue_free_readout=13000.,requested_geometries=15,eligible_geometries=14,
        initial_overlap_exclusions=[dict(resting_circumference=1.,fraction=.3)],
        categorical_fate_invariant=True,observed_category="cue-associated single peak",
        descriptor_ranking_status="uninformative_constant_endpoint_class",
        classifier_accuracy_reported=False,
        ranking_limit="This constant-class T1000 sample cannot discriminate physical, angular or other descriptors for class prediction. It does not establish class invariance for shorter cues, other doses, untested placements or sizes.",
        continuous_readouts_retained=True,within_size_distance_angle_collinear=True,within_size_ranges=ranges,
        interval_correspondence=dict(uniform_field_native_operator_carryover="verified by retained Spring operator and numerical validation",
            exact_nonuniform_interval_reflection="not_applicable_to_declared_mean_force_closure",
            reason="Native interval extension samples the two endpoint Rac force laws; the declared ring averages the nonlinear force law around its perimeter. These generally differ for nonuniform fields. Uniform material dilution/diffusion carryover is preserved, but does not establish an exact nonuniform reflected-interval trajectory.",
            native_2d_mechanical_equivalence=False),
        current_counterparts=dict(basal_organization_gain_recovery=str((TOPIC.parent/"spring_function/README.md").relative_to(REPO)),
            challenged_regional_function=str((TOPIC.parent/"spring_challenged_function/README.md").relative_to(REPO)),
            three_size_distance_duration=str((TOPIC/"spring_distance/spring_distance_duration.png").relative_to(REPO))),
        older_comparison_scope="The held-out-size descriptor fit and exact chemical interval/ring reflection maps cover Wave-pinning, Goryachev and Holmes only; they are historical chemical-only support.",
        native_units="Uncalibrated author-model length and time; no cross-model perimeter or clock conversion",
        source_sha256={str(source):sha(source),str(Path(__file__)):sha(__file__)},output_table_sha256=sha(table))
    (META/"descriptor_coverage.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(dict(passed=True,eligible_geometries=14,endpoint_classes=1,ranking="uninformative",table=str(table)),indent=2))


if __name__=="__main__":main()
