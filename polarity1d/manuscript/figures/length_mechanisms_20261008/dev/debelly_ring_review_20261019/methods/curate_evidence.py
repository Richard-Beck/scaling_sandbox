#!/usr/bin/env python3
"""Retain a compact, source-pinned assessment of the stopped ring attempt.

This copies existing small evidence products and selects existing scalar records.
It does not evolve a model, recompute an analysis, or choose a new closure.
The retained report reads its internal evidence; source paths are provenance.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
MANUSCRIPT = REPO / "polarity1d/manuscript"
RUNS = MANUSCRIPT / "data/runs"
OUT = MANUSCRIPT / "figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name: str):
    return json.loads((RUNS / name).read_text())


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main():
    now = datetime.now(timezone.utc).isoformat()
    OUT.mkdir(parents=True, exist_ok=True)
    records = []

    def retain(source: str, dest: str, purpose: str):
        src = REPO / source
        dst = OUT / dest
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        records.append({"source_path": source, "source_sha256": sha(src),
                        "retained_path": dest, "retained_sha256": sha(dst),
                        "bytes": dst.stat().st_size, "transformation": "byte-for-byte copy",
                        "purpose": purpose})

    def run_copy(source: str, dest: str, purpose: str):
        retain("polarity1d/manuscript/data/runs/" + source, "evidence/" + dest, purpose)

    # Four small evidence groups; the separate bounded event viewer is curated
    # by the implementation agent. No full trajectories or sparse matrices here.
    copies = [
        ("debelly_protrusive_closure_design_20261019/mechanics_evidence.json", "01_closure/native_source_findings.json", "Primary-source locators distinguishing truncated native geometry and prescribed membrane/cortex actuation."),
        ("debelly_ring_literature_20261019/literature_review.json", "01_closure/literature_review.json", "Complementary models, primary links, inspected-code status and biological differences; no replacement selected."),
        ("debelly_gap_on_face_challenge_C50_20261019/analysis/comparison.svg", "02_alignment/comparison.svg", "The same qualitative traveling/retargeting behavior survives fixed grid-face gap placement and aligned cue placement."),
        ("debelly_gap_on_face_challenge_C50_20261019/analysis/comparison.json", "02_alignment/comparison.json", "Compact numerical values and limits behind the alignment comparison."),
        ("debelly_gap_on_face_challenge_C50_20261019/analysis/validation_analysis.json", "02_alignment/analysis_validation.json", "Independent agreement audit of the retained alignment analysis."),
        ("debelly_paired_mca_20261019/analysis/four_branch_comparison.svg", "03_numerics/comparison.svg", "Native versus conservative paired MCA responses at two timesteps in fixed physical coordinates."),
        ("debelly_paired_mca_20261019/analysis/four_branch_comparison.json", "03_numerics/comparison.json", "Four response summaries, complete observed intervals, pool budgets and interpretation limits."),
        ("debelly_paired_mca_20261019/analysis/four_branch_comparison.csv", "03_numerics/comparison.csv", "Tabular version of the four-branch numerical comparison."),
        ("debelly_paired_mca_20261019/analysis/validation_analysis.json", "03_numerics/analysis_validation.json", "Independent analysis agreement audit."),
        ("debelly_paired_mca_20261019/campaign/validation_result.json", "03_numerics/terminal_validation.json", "Independent completed raw/protocol/checkpoint/pool audit, distinct from physical adequacy."),
        ("debelly_gap_factorial_C50_20261019/analysis/impl_four_corner_table.json", "03_numerics/native_factorial.json", "Native two-grid/two-timestep sensitivity, with one complete target-trough return per corner."),
        ("debelly_gap_factorial_C50_20261019/validation_result.json", "03_numerics/native_factorial_validation.json", "Independent native four-corner terminal and selected-cycle budget audit."),
        ("debelly_gap_oscillation_C50_20261019/impl_intervention_comparison.svg", "04_mechanism/retarget_control.svg", "Suppressing the complete retarget operation removes the selected return while retaining polarization."),
        ("debelly_gap_oscillation_C50_20261019/impl_intervention_analysis.json", "04_mechanism/retarget_control.json", "Exact fixed-topology fork, return window and finite-time fates."),
        ("debelly_gap_oscillation_C50_20261019/interventions/fixed_topology/validation_result.json", "04_mechanism/retarget_validation.json", "Independent fixed-topology control audit."),
        ("debelly_gap_linear_mode_C50_20261019/mechanics_mode_growth.svg", "04_mechanism/onset_growth.svg", "Finite-time growth of the observed physical reflection-antisymmetric mode in all seven native fields."),
        ("debelly_gap_linear_mode_C50_20261019/mechanics_observed_mode.json", "04_mechanism/onset_summary.json", "Observed growth, field parity, selected one-step allocations and nonsmooth-map limitations."),
        ("debelly_gap_linear_mode_C50_20261019/validation_followup.json", "04_mechanism/onset_selected_validation.json", "Selected observed-mode/controller/production checks; scope is narrower than a full nonlinear explanation."),
        ("debelly_gap_travel_clock_C50_20261019/mechanics_travel_profiles.svg", "04_mechanism/travel_shape.svg", "Ten saved states show evolving width, edge chemistry, tension and declining old-gap actuation."),
        ("debelly_gap_travel_clock_C50_20261019/mechanics_travel_profiles.json", "04_mechanism/travel_summary.json", "Selected actual states and exact one-step Rac-production/flow allocations; no evolved ablation."),
        ("debelly_gap_travel_clock_C50_20261019/validation_result.json", "04_mechanism/travel_validation.json", "Independent selected-state geometry, source-accounting and passage-time audit."),
    ]
    for args in copies:
        run_copy(*args)

    retain("polarity1d/manuscript/DEBELLY_BOUNDARIES.md",
           "evidence/01_closure/corrected_linear_boundary_scope.md",
           "Protected corrected linear model scope; the stopped ring attempt does not replace it.")

    alignment = read("debelly_gap_on_face_challenge_C50_20261019/analysis/comparison.json")
    numerics = read("debelly_paired_mca_20261019/analysis/four_branch_comparison.json")
    factorial = read("debelly_gap_factorial_C50_20261019/analysis/impl_four_corner_table.json")
    onset = read("debelly_gap_linear_mode_C50_20261019/mechanics_observed_mode.json")
    travel = read("debelly_gap_travel_clock_C50_20261019/mechanics_travel_profiles.json")
    control = read("debelly_gap_oscillation_C50_20261019/impl_intervention_analysis.json")
    angles_name = "debelly_gap_initial_angles_C50_20261019/analysis/angle_comparison.json"
    angles = read(angles_name)
    angle_selection = {
        "scope": "Selection of the existing five-angle table, without recomputation.",
        "source_path": "polarity1d/manuscript/data/runs/" + angles_name,
        "source_sha256": sha(RUNS / angles_name),
        "table": angles["table"], "definitions": angles["definitions"],
        "independent_raw_validation": angles["independent_raw_validation"],
        "limit": "Tiny cardinal offsets increase during finite apparently stationary preparations. These tests do not establish asymptotic stability or uniquely identify the seed of drift.",
    }
    selected_path = OUT / "evidence/02_alignment/initial_angle_table.json"
    dump(selected_path, angle_selection)
    records.append({"source_path": angle_selection["source_path"],
                    "source_sha256": angle_selection["source_sha256"],
                    "retained_path": str(selected_path.relative_to(OUT)),
                    "retained_sha256": sha(selected_path), "bytes": selected_path.stat().st_size,
                    "transformation": "Selected existing table, definitions and audit pointer; no numerical recomputation",
                    "purpose": "Retain fresh-angle results without the original lengthy per-frame summaries."})

    literature = read("debelly_ring_literature_20261019/literature_review.json")
    lit_ids = ["ivsic_2026", "mietke_active_surface_family", "werner_2022", "chen_2023"]
    lit = [{k: x.get(k) for k in ["id", "title", "citation", "doi", "geometry", "biology", "feasibility", "use", "code", "sources", "companions"]}
           for item in lit_ids for x in literature["candidates"] if x["id"] == item]

    evidence_groups = [
        {"id": "closure", "title": "The native model does not provide the missing closed protrusive cap",
         "path": "evidence/01_closure/native_source_findings.json", "figure": None,
         "why_retain": "Explains why changing a boundary condition alone cannot preserve the complete native protrusive mechanism.",
         "limit": "Existing source inspection and design requirements, not an implemented finite-cap solution."},
        {"id": "alignment", "title": "Alignment changes do not remove traveling/retargeting dynamics",
         "path": "evidence/02_alignment/comparison.json", "figure": "evidence/02_alignment/comparison.svg",
         "why_retain": "Rules out the tested grid-face/cue-alignment explanation as a sufficient remedy.",
         "limit": "Three prepared C50 branches and one complete observed return each; not every possible discretization artifact is excluded."},
        {"id": "numerics", "title": "Timing changes with timestep and MCA update",
         "path": "evidence/03_numerics/comparison.json", "figure": "evidence/03_numerics/comparison.svg",
         "csv": "evidence/03_numerics/comparison.csv",
         "why_retain": "The most direct evidence against using the current finite return time as a robust assay measurement.",
         "limit": "Own-method basal preparations and only two timestep/grid levels; no convergence or unique numerical-source causality claim."},
        {"id": "mechanism", "title": "Mechanistic progress is useful but incomplete",
         "path": "evidence/04_mechanism/retarget_control.json", "figure": "evidence/04_mechanism/retarget_control.svg",
         "additional_figures": ["evidence/04_mechanism/onset_growth.svg", "evidence/04_mechanism/travel_shape.svg"],
         "why_retain": "Preserves the bounded retargeting control, reproducible onset mode and evolving travel-state evidence.",
         "limit": "The complete retarget operation is tested jointly; instantaneous source allocation does not identify a unique oscillation clock or establish biological necessity."},
    ]
    assessment = {
        "generated_utc": now, "producer_state": "finished_stopped_development_assessment",
        "objective_state": "Stopped by user; no further model development or campaigns authorized by this report.",
        "title": "Periodic deBelly: significant unresolved hurdles",
        "verdict": "This attempt did not establish a faithful, numerically robust periodic deBelly model preserving its full Rac-driven protrusion mechanism.",
        "executive_summary": "The point-gap construction retains both native Rac-driven flank arms and produces nontrivial polarized mechanics, but it introduces an omitted-cap/interface interpretation whose material, geometry and actuator-reference assumptions are not a demonstrated closed-cell model. Its polarized domain travels and retargets the gap; tested alignment changes do not remove that behavior, and finite return times change substantially with numerical treatment. Useful onset and travel diagnostics were obtained, but a complete physical explanation and converged assay counterpart were not achieved.",
        "is_promising": "Mechanistically interesting as a future research project; not ready as a periodic counterpart in the current cross-model assay comparison. The shortest publication path is to omit it from the current ring assays while retaining corrected linear deBelly.",
        "scope": {"geometry": "Fixed circular assay coordinate with mechanically cut flanks; chemistry remains periodic.",
                  "main_protocol": numerics["protocol"],
                  "MCA": "Membrane–cortex attachment, represented by bound and unbound pools.",
                  "do_not_infer": ["A requirement that every challenged state become stationary.", "Loss of polarization merely because a domain moves.", "Cell migration speed from the Rac first-moment phase.", "A continuous material membrane, resolved normal protrusion or force-free cell torque law from the point-gap construction."]},
        "five_findings": [
            {"id": "native_closure", "headline": "A faithful closed protrusive closure is missing.",
             "finding": "The native one-dimensional model prescribes cortex velocity and separately relaxed membrane displacement at a front. The published axisymmetric formulation also truncates the protrusive cap rather than providing a regular uncut cap law. A fully periodic bulk control removes the direct Rac-to-membrane flank arm; a cortex-only seam does too. Neither meets the requested preservation of direct Rac-driven protrusion.",
             "status": "Source-grounded structural limitation.", "evidence": "closure",
             "remaining": "A declared driven interface can be studied as an analogue, but its port reactions/work, absolute actuator reference and material reservoir must remain explicit. A continuous cap instead requires new geometric, constitutive and material-return choices; those were not adopted."},
            {"id": "gap_travel", "headline": "The protrusive gap construction retains polarity but develops traveling and retargeting dynamics.",
             "finding": "In the retained C50 long-cue response, the Rac domain remains polarized while moving relative to the fixed-at-birth gap. Healing and new opening reset the gap location; Rac motion continues in its previous direction initially, with the later turn outside the dense 30-second event window. Suppressing the complete retarget operation removes the selected return without removing polarization.",
             "evidence": "mechanism", "measured_return_window_s": control["preregistered_return_window"],
             "return_window": control["return_window"],
             "limit": "This fork implicates the combined lifecycle/mechanical/material operation, not a uniquely isolated controller, transport term or biological oscillator."},
            {"id": "alignment", "headline": "The tested grid-alignment repair does not cure the motion.",
             "finding": "Fixed existing-face gap births, with either the original cue center or a grid-face-aligned cue, retain the same qualitative traveling/retargeting response. Fresh cardinal preparations show tiny growing offsets despite passing finite basal gates, so the original zero-angle preparation is not evidence of asymptotic stability.",
             "rows": [{"key": x["key"], "point_trough_return_intervals_s": x["point_observed_trough_intervals"], "cue_location_um": x.get("own_cue_location")} for x in alignment["rows"]],
             "evidence": "alignment", "limit": "A negative result for the tested placement/alignment fixes; not proof that all discrete effects are absent."},
            {"id": "numerical_dependence", "headline": "Finite return times depend on the numerical pipeline.",
             "native_factorial": [{"N": x["N"], "dt_s": x["dt"], "point_trough_return_intervals_s": x["point_observed_trough_intervals"]} for x in factorial["rows"]],
             "native_paired": [{"key": x["key"], "dt_s": x["dt"], "point_trough_return_intervals_s": x["point_observed_trough_intervals"], "own_basal_qualified": x["basal_qualified"]} for x in numerics["rows"]],
             "finding": "Halving the native timestep shortens the observed return by about 13%; doubling the mesh changes it by about 2–3% at these levels. Conservative pairing of bound/unbound exchange removes the known numerical creation term and retains oscillatory responses, but changes their return intervals to about 11,000 seconds.",
             "evidence": "numerics", "limit": "Only one complete native interval per corner and two paired intervals per timestep are observed. Requalification is part of each pipeline. Pairing alone is not an accuracy or convergence proof, and period changes do not isolate the numerical exchange residual as the cause."},
            {"id": "mechanism_progress", "headline": "Onset and travel were characterized, but the mechanism is not fully resolved.",
             "finding": "The saved near-symmetric preparation develops a common physical reflection-antisymmetric mode across all seven fields and four controllers, with an e-folding time near 2,840 seconds. At finite amplitude, the active width, MCA, Rho and tension evolve as the domain travels; the old gap's direct command progressively switches off before its connectivity threshold is crossed. Rac has no direct advection: movement reflects spatially unequal reaction and diffusion, coupled indirectly to mechanics through MCA and Rho.",
             "observed_mode_growth": onset["growth_fits"]["Rac_phase"],
             "travel_interpretation": travel["interpretation"],
             "evidence": "mechanism", "remaining": "Instantaneous field/source accounts and a complete-retarget fork do not establish the unique feedback setting growth, reversal or cycle timing. No continuum instability, limit cycle or physical cap model was established."},
        ],
        "recommendations": [
            {"rank": 1, "action": "OMIT from current ring assays; keep corrected linear deBelly.",
             "rationale": "Protect the supported linear results and the interpretable ring comparison. Label this retained report as an unsuccessful extension attempt with useful diagnostics, rather than include an unresolved counterpart as if it were the same model.",
             "scope": "No change to corrected linear figures, other supported ring models or independent Spring work."},
            {"rank": 2, "action": "MODIFY only as a bounded future redesign if membrane–cortex attachment biology is essential.",
             "rationale": "Choose explicitly between an externally referenced driven omitted-cap/interface analogue and a finite active-cap model. The former is closer to the prescribed-boundary native assay; the latter can support stronger continuous-material claims but requires new geometry, material measures, recycling and mechanical closure. Preserve both direct Rac arms and separate physical fluxes from numerical correction.",
             "minimum_future_gate": "Declare the biological interpretation first, then verify work/material/reference accounting and numerical robustness before comparative assays. No new construction is selected here."},
            {"rank": 3, "action": "REPLACE or complement it if genuinely closed mechanics is the scientific requirement.",
             "rationale": "Inspect a published closed-cell/active-surface model with appropriate code and biology. Treat it as a different model, not a validated deBelly port. Existing literature offers closure ideas and complementary mechanisms, not a drop-in replacement.",
             "candidates": lit},
        ],
        "literature_judgment": {"assessment": "No inspected primary candidate preserves the complete deBelly Rac/Rho–MCA–tension feedback while supplying demonstrated uncut deforming mechanics.",
                               "implementation_caution": "Public code was inspected, not run or independently validated here. A published geometry or source repository does not guarantee that it implements the intended assays or has suitable stability, conservation and convergence properties.",
                               "recommended_inspection_order": ["Ivšić 2026 for a complementary closed-cell polarity implementation", "Mietke active-surface family for normal/tangential mechanics and surface transport", "Werner 2022 if explicit membrane–cortex linkers are indispensable", "Chen 2023 for an economical closed-curve mechanical template, with different adhesion biology"],
                               "evidence": "evidence/01_closure/literature_review.json"},
        "evidence_groups": evidence_groups,
        "recommended_display": ["evidence/02_alignment/comparison.svg", "evidence/03_numerics/comparison.svg", "evidence/04_mechanism/retarget_control.svg"],
        "optional_diagnostic_figures": ["evidence/04_mechanism/onset_growth.svg", "evidence/04_mechanism/travel_shape.svg"],
        "bounded_viewer": {"path": "viewer/index.html", "role": "Inspect all seven actual fields before/after each retained gap event; dense windows have native 0.5-second resolution, with separate same-clock topology snapshots. It is evidence for the attempted construction, not a functioning-model claim."},
        "retention": {"large_histories_copied": False, "new_analysis_or_trajectory": False, "source_artifacts_modified": False,
                      "provenance": "evidence_manifest.json", "asset_count": len(records),
                      "data_note": "Copied summaries retain original source paths and hashes as historical provenance. Those paths are not runtime dependencies of this compact assessment."},
    }
    dump(OUT / "executive_assessment.json", assessment)
    methods = OUT / "methods"
    methods.mkdir(exist_ok=True)
    shutil.copyfile(Path(__file__), methods / "curate_evidence.py")
    manifest = {"generated_utc": now, "producer_state": "finished_small_evidence_curation",
                "scope": "Existing-output curation only; no model selected and no new analysis or dynamics.",
                "curation_source": str(Path(__file__).relative_to(REPO)), "curation_source_sha256": sha(Path(__file__)),
                "assets": records, "total_asset_bytes": sum(x["bytes"] for x in records),
                "executive_assessment_sha256": sha(OUT / "executive_assessment.json"),
                "runtime_dependencies": "Only the copied internal evidence and separately retained viewer are needed to read the final report."}
    dump(OUT / "evidence_manifest.json", manifest)
    print(json.dumps({"out": str(OUT), "assets": len(records), "bytes": manifest["total_asset_bytes"],
                      "executive_sha256": manifest["executive_assessment_sha256"]}))


if __name__ == "__main__":
    main()
