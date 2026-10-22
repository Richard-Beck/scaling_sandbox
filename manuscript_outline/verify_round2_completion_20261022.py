"""Requirement-by-requirement acceptance of the completed round-2 package.

Scientific campaign audits are prerequisites, not substitutes for checking the
actual plotted values, retained pages, document contents and delivery assets.
"""
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "manuscript_outline/revision_round2_20261022"
FIG = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022"
RUN = ROOT / "polarity1d/manuscript/data/runs"
csv.field_size_limit(50_000_000)


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(4 << 20),b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def check(result, requirement, condition, evidence):
    assert condition, requirement
    result.append(dict(requirement=requirement,status="proved",evidence=evidence))


def main():
    checks=[]
    document=read(ROOT/"manuscript_outline/outline_v5.json")
    baseline=read(ROOT/"manuscript_outline/outline_v4.json")
    section={s["id"]:s for s in document["sections"]}
    integration=read(DEST/"completion_integration.json")
    v2=read(FIG/"figure2/supporting/metadata/verification.json")
    a2=read(FIG/"figure2/supporting/metadata/brief_acceptance.json")
    v3=read(FIG/"figure3/verification.json")
    v4=read(FIG/"figure4/supporting/metadata/completion_scientific_audit.json")
    check(checks,"Complete scientific execution, not a successful agent exit",v2["complete"] and a2["complete"] and v3["passed"] and v4["passed"],
        [relative(FIG/"figure2/supporting/metadata/verification.json"),relative(FIG/"figure3/verification.json"),relative(FIG/"figure4/supporting/metadata/completion_scientific_audit.json")])
    # Figure 2: inspect actual frozen design and every final measurement cell.
    campaign=RUN/"signature_holmes_retention_round2_20261022"
    manifest=read(campaign/"manifest.json")
    audit=read(campaign/"results/independent_native_readout_audit.json")
    check(checks,"3.1 All 32 Holmes trajectories independently audited",audit["passed"] and audit["workers"]==32 and len(audit["trajectories"])==32,relative(campaign/"results/independent_native_readout_audit.json"))
    cases=manifest["cases"]
    deterministic=[c for c in cases if c["kind"]=="train"]
    check(checks,"3.1 L31/L35/L39 × six prescribed gaps",{(c["length"],c["gap"]) for c in deterministic}=={(length,gap) for length in (31,35,39) for gap in (0,.1,.5,1,4,16)},relative(campaign/"manifest.json"))
    reference=manifest["reference_calibration"]
    check(checks,"3.1 Four pulses each one quarter of 1.5×L27 calibration",reference["duration"]==18125 and reference["pulse_duration"]==1.5*reference["duration"]/4 and all(c["pulse_count"]==4 and c["calibration_duration"]==reference["duration"] for c in cases),relative(campaign/"manifest.json"))
    check(checks,"3.1 L39 S3-matched bounded-input check and controls",Counter(c["kind"] for c in cases)==dict(train=18,noise=12,mean=1,zero=1) and all(c["length"]==39 and c["gap"]==16 for c in cases if c["kind"]!="train"),relative(campaign/"manifest.json"))
    for row in audit["trajectories"]:
        assert sha(row["result_path"])==row["result_sha256"]
        assert all(row[k] for k in ("all_artifact_hashes_verified","all_readouts_independently_recomputed","all_native_phases_agree","field_tail_independently_recomputed"))
    endpoint=rows(FIG/"figure2/holmes_large_L_endpoints.csv")
    check(checks,"3.1 D at every pulse and inter-pulse gap endpoint",len(endpoint)==32*7 and Counter(r["kind"] for r in endpoint)==dict(on=128,gap=96),relative(FIG/"figure2/holmes_large_L_endpoints.csv"))
    for row in endpoint:
        if row["D"] and row["share"]:
            assert math.isclose(float(row["D"]),2*float(row["share"])-1,rel_tol=1e-10,abs_tol=1e-10)
    noise=rows(FIG/"figure2/holmes_L39_bounded_noise_retention.csv")
    assert len(noise)==9
    for summary in noise:
        selected=[r for r in endpoint if r["forcing_kind"]==summary["kind"] and r["kind"]=="gap"]
        tasks=sorted({r["task"] for r in selected},key=int)
        sample=[sorted((r for r in selected if r["task"]==task),key=lambda r:float(r["time"]))[int(summary["gap_index"])-1] for task in tasks]
        eligible=[r for r in sample if r["S3_denominator_eligible"]=="True"]
        assert len(sample)==int(summary["seeds"]) and len(eligible)==int(summary["eligible"])
        assert sum(r["S3_success"]=="True" for r in eligible)==int(summary["successes"])
        retention=[float(r["share"])/float(r["pulse_share"]) for r in sample]
        assert math.isclose(min(retention),float(summary["minimum_raw_share_retention"]),rel_tol=1e-10)
        assert math.isclose(max(retention),float(summary["maximum_raw_share_retention"]),rel_tol=1e-10)
    noisy=[r for r in noise if r["kind"]=="noise"]
    check(checks,"3.1 Noise successes and failures independently counted by gap without pseudoreplication",all(int(r["seeds"])==12 for r in noisy) and any(int(r["successes"])<int(r["eligible"]) for r in noisy) and min(float(r["minimum_raw_share_retention"]) for r in noisy)<.1,relative(FIG/"figure2/holmes_L39_bounded_noise_retention.csv"))
    classes=rows(FIG/"figure2/class_membership_check.csv")
    check(checks,"3.3 Every class-table cell has numerical evidence",len(classes)==20 and len({(r["model"],r["signature"]) for r in classes})==20 and all(r["decision"] and r["qualification"] for r in classes),relative(FIG/"figure2/class_membership_check.csv"))
    check(checks,"Author class-membership stop remains untriggered",not v2["class_membership_changed"] and all(r["membership_changed"]=="False" for r in classes),relative(FIG/"figure2/class_membership_check.csv"))
    fates=rows(FIG/"figure2/train_fates_and_retention.csv")
    check(checks,"3.2D Complete maps distinguish 18 large-L retention endpoints",len(fates)==180 and sum(r["status"]=="observed_retention_without_target_calibration" for r in fates)==18,relative(FIG/"figure2/train_fates_and_retention.csv"))
    check(checks,"7:03 correction: failed deBelly field tails remain visible",sum(r["model"]=="debelly" and r["tail_pass"]=="False" for r in fates)==31 and v2["debelly_failed_field_tails_visible"]==31 and "31" in " ".join(section["figure2"]["legend"]),relative(FIG/"figure2/train_fates_and_retention.csv"))
    initiation=rows(FIG/"figure2/initiation_normalized.csv")
    check(checks,"3.2B deBelly first-update ratios remain unavailable",all(not r["ratio_estimate"] for r in initiation if r["model"]=="debelly" and float(r["length"])!=19),relative(FIG/"figure2/initiation_normalized.csv"))
    traces=rows(FIG/"figure2/gap_retention_traces.csv")
    check(checks,"3.2E All five L39 histories include new Holmes measurements",{r["model"] for r in traces}=={"goryachev","wave_pinning","holmes_model3","debelly","spring"},relative(FIG/"figure2/gap_retention_traces.csv"))
    check(checks,"3.1 No continuous commitment by 100,000 explicitly reported","No continuous commitment by 100,000" in " ".join(section["figure2"]["legend"]),relative(FIG/"figure2/legend.json"))
    for ident,value in a2["requirements"].items():
        check(checks,"Figure 2 acceptance: "+ident,value["passed"],value.get("evidence"))
    # Figure 3: new evidence actually replaces the nine original horizon bounds.
    evidence=rows(FIG/"figure3/B_complete_with_extensions.csv")
    extensions=[r for r in evidence if r["extension_status"]]
    check(checks,"2.4 Nine approved exact-endpoint continuations completed",len(extensions)==9 and v3["approved_extensions_completed"]==9 and len(v3["native_field_validation"])==9,relative(FIG/"figure3/verification.json"))
    check(checks,"2.4 Six resolved and three censored at 80,000",sum(r["resolution_status"]=="observed_saved_sample_crossing" for r in extensions)==6 and sum(r["resolution_status"]=="right_censored_at_completed_horizon" and float(r["resolution_lower"])==80000 for r in extensions)==3,relative(FIG/"figure3/B_complete_with_extensions.csv"))
    check(checks,"2.1 All 48 decomposition terms and interpolation queries checked",v3["decomposition_terms"]==48 and v3["interpolation_queries"]==48 and v3["finite_product_identities"]==12,relative(FIG/"figure3/verification.json"))
    terms=rows(FIG/"figure3/decomposition_terms.csv")
    for row in terms:
        if not math.isfinite(float(row["upper"])):
            assert not row["point"] or not math.isfinite(float(row["point"]))
    check(checks,"2.1 Censoring never substituted as a ratio point",v3["bounds_have_no_point_values"],relative(FIG/"figure3/decomposition_terms.csv"))
    message="Separation slows the resolution of coexisting regions at short range; at matched separation, cell size dominates both driven takeover and resolution."
    check(checks,"7:03 correction: exact Figure 3 message",message in section["figure3"]["title"] and section["figure3"]["legend"][0].startswith(message),relative(FIG/"figure3/content.json"))
    check(checks,"7:03 correction: measured short-range gap ranges reported",(FIG/"figure3/short_range_evidence.csv").is_file() and "short" in " ".join(section["figure3"]["legend"]).lower(),relative(FIG/"figure3/short_range_evidence.csv"))
    check(checks,"2.2 Native kymograph values and source hashes checked",sum(v3["native_field_csv_rows"].values())==3183488 and v3["native_source_checksums_verified"]==6,relative(FIG/"figure3/verification.json"))
    check(checks,"2.2 Holmes bounds stated after reanalysis",v3["holmes_censored"]==12 and v3["holmes_eligible"]==23 and section["figure3"]["coverage"]["B_resolution_lower_bounds"]==46,relative(FIG/"figure3/verification.json"))
    check(checks,"2.3 Full former main/supplement pages retained",v3["preserved_original_page_content_streams"]==9,relative(FIG/"figure3/verification.json"))
    # Figure 4 audit covers all original raw scalar data and lag rules.
    check(checks,"1.1 Existing-output correlation, rate, coverage and whole-seed uncertainty",v4["native_scalar_hashes"]==2364 and v4["condition_lag_rows"]==15540 and v4["exported_seed_lags"]==156660,relative(FIG/"figure4/supporting/metadata/completion_scientific_audit.json"))
    data=FIG/"figure4/supporting/data"
    lags=rows(data/"all_lag_size_ratios.csv")
    check(checks,"1.4 Every usable lag and lag-consistency flag displayed",v4["all_lag_pages"]==56 and all(r["non_exponential"] in ("True","False") for r in lags),relative(FIG/"figure4/figure4_all_lag_ratios.pdf"))
    tolerance=rows(data/"noise_tolerance.csv")
    designed=set(float(f) for r in tolerance for f in json.loads(r["designed_inputs"]))
    check(checks,"1.1 Fixed-reference-lag tolerance includes all added designed inputs",{.005,.015,.02,.05}<=designed and {float(r["reference_lag"]) for r in tolerance}=={200,2000},relative(data/"noise_tolerance.csv"))
    h1=rows(data/"S0_posthoc_H1.csv")
    main=[r for r in h1 if r["model"]=="all" and r["view"]=="main" and r["forced_inclusive"]=="False"]
    check(checks,"Author-authorized Case A: amended pass with frozen calls alongside",len(main)==2 and all(r["status"]=="pass" and r["frozen_rate_status"]=="inconclusive" and r["post_hoc"]=="True" for r in main),relative(data/"S0_posthoc_H1.csv"))
    check(checks,"1.2 Both observation arms, all views, both forcing choices",len(h1)==140 and len({(r["arm"],r["view"],r["forced_inclusive"]) for r in h1})==28,relative(data/"S0_posthoc_H1.csv"))
    check(checks,"Author correction: H2/H3/A6 gates retired",v4["retired_gate_rows"]==196,relative(data/"S0_retired_gates.csv"))
    check(checks,"1.4 Old tau and seven-view robustness retained",(data/"old_tau_vs_new_rate.csv").is_file() and (FIG/"figure4/supporting/vectors/legacy_seven_view_robustness.pdf").is_file(),relative(FIG/"figure4/supporting/metadata/main_render_receipt.json"))
    oldfour=FIG.parent/"signature_figure4_20261020/supporting/vectors"
    previous=[*PdfReader(oldfour/"figure4_revision_main.pdf").pages,*PdfReader(oldfour/"figure4_revision_complete_supplement.pdf").pages]
    combined=PdfReader(FIG/"figure4/figure4_supplement.pdf").pages
    check(checks,"1.4 Complete superseded main and original supplement preserved",all(a.get_contents().get_data()==b.get_contents().get_data() for a,b in zip(combined[-len(previous):],previous)),dict(original_pages=len(previous),supplement=relative(FIG/"figure4/figure4_supplement.pdf")))
    current=rows(data/"main_size_ratios.csv")
    check(checks,"Author correction: Goryachev excluded from Figure 4 main",all(r["model"]!="goryachev" for r in current),relative(data/"main_size_ratios.csv"))
    legend=" ".join(section["figure4"]["legend"])
    check(checks,"Author correction: geometry and multiplicity descriptive, weak inputs model-specific","post hoc" in legend.lower() and "model-specific" in legend and "do not establish multiple fronts as the cause" in legend,relative(FIG/"figure4/supporting/metadata/main_legend.md"))
    # Artifact checks apply to every requested output, not only the integrated PDF.
    palette={"wave_pinning":"#2878b5","goryachev":"#c83b3b","debelly":"#82599c","holmes_model3":"#36944c","spring":"#c59a25"}
    four_render=read(FIG/"figure4/supporting/metadata/main_render_receipt.json")
    palette4={k.replace("spring_ring","spring"):v for k,v in four_render["palette"].items()}
    check(checks,"0 Fixed model hues across all main figures",v2["palette"]==palette and v3["model_palette_verified"] and all(palette[k]==v for k,v in palette4.items()),relative(FIG/"figure4/supporting/metadata/main_render_receipt.json"))
    mainstems={2:"figure2_round2",3:"figure3_main",4:"figure4_main"}
    artifact=[]
    for number,stem in mainstems.items():
        folder=FIG/f"figure{number}"
        for extension in ("png","pdf","svg"):
            path=folder/(stem+"."+extension)
            check(checks,f"4 Figure {number} main {extension.upper()}",path.is_file() and path.stat().st_size>1000,relative(path))
            artifact.append(dict(path=relative(path),sha256=sha(path),bytes=path.stat().st_size))
        supplement=ROOT/section[f"figure{number}"]["supplement_pdf"]
        check(checks,f"4 Figure {number} complete supplement",supplement.is_file(),relative(supplement))
        artifact.append(dict(path=relative(supplement),sha256=sha(supplement),bytes=supplement.stat().st_size))
        check(checks,f"4 Figure {number} README and plotted CSVs",(folder/"README.md").is_file() and bool(list(folder.rglob("*.csv"))),relative(folder/"README.md"))
        check(checks,f"0 Figure {number} message-led title and legend",section[f"figure{number}"]["legend"][0] and not any(word in section[f"figure{number}"]["title"] for word in ("AUDITED FINAL","IN PROGRESS","pending")),relative(folder))
        assert " ".join(section[f"figure{number}"]["legend"]).lower().count("native units")<=1
    check(checks,"0 Main panel counts at most six",a2["requirements"]["3.2_six_main_panels_A_to_F"]["passed"] and v3["main_panels"]==5 and four_render["main_panels"]==5,dict(figure2=6,figure3=5,figure4=5))
    original_inputs={r["path"]:r["sha256"] for r in read(ROOT/"manuscript_outline/outline_v4.build_manifest.json")["inputs"]}
    for ident in ("figure1","figure5","figure6"):
        original=next(s for s in baseline["sections"] if s["id"]==ident)
        unchanged=section[ident]==original and all(sha(ROOT/a["path"])==original_inputs[a["path"]] for a in original["assets"])
        check(checks,f"4 {ident} content and artwork unchanged",unchanged,relative(ROOT/"manuscript_outline/outline_v4.build_manifest.json"))
    check(checks,"7:03 correction: cross-figure narrative and Spring exception",integration["cross_figure_narrative"]==document["overview"][1] and "mainly because of size itself" in document["overview"][1] and "Spring" in document["overview"][1],relative(DEST/"completion_integration.json"))
    outline_check=read(ROOT/"manuscript_outline/outline_v5.validation.json")
    current_build=read(ROOT/"manuscript_outline/outline_v5.build_manifest.json")
    assert all(sha(ROOT/r["path"])==r["sha256"] for r in current_build["inputs"])
    actual_pdf=sha(ROOT/"manuscript_outline/outline_v5.pdf")
    check(checks,"4 Integrated outline_v5 renders and exact current source hashes pass",integration["final_outline_v5_integrated"] and outline_check["passed"] and outline_check["stage"]==document["stage"] and outline_check["pdf_sha256"]==current_build["pdf"]["sha256"]==actual_pdf,relative(ROOT/"manuscript_outline/outline_v5.validation.json"))
    for path in [DEST/"completion.html",FIG/"figure3/index.html",FIG/"figure4/index.html"]:
        content=path.read_text()
        assert not re.search(r'http-equiv\s*=\s*[\"\x27]?refresh|setInterval\s*\([^;]*location\.reload',content,re.I)
        assert "Generated" in content or "generated" in content
    check(checks,"HTML reports show generation time without timed self-refresh",True,relative(DEST/"completion.html"))
    # Extract the actual final manuscript chapters for convenient standalone review.
    build=read(ROOT/"manuscript_outline/outline_v5.build_manifest.json")
    reader=PdfReader(ROOT/"manuscript_outline/outline_v5.pdf");writer=PdfWriter()
    for page in build["page_inventory"]:
        if page["section_id"] in ("figure2","figure3","figure4"):
            writer.add_page(reader.pages[page["number"]-1])
    with (DEST/"revised_figures_2_4.pdf").open("wb") as stream:
        writer.write(stream)
    report=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),requirements_checked=len(checks),
        checks=checks,artifacts=artifact,all_required_work_complete=True,remaining_scientific_qualifications=section["figure2"]["limitations"]+section["figure3"]["limitations"]+section["figure4"]["limitations"],
        delivery_commit_and_remote_verification="Recorded separately after commit/push; this report covers scientific work and artifacts.")
    (DEST/"completion_acceptance.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(dict(passed=True,requirements_checked=len(checks),chapters_pages=len(writer.pages))))


if __name__=="__main__":
    main()
