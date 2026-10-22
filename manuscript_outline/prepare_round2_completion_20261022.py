"""Integrate the completed round-2 figures and the author's final corrections.

The Figure 4 hypothesis amendment is explicitly post hoc. Figures 1, 5 and 6
are copied from outline_v4 without any change to their content or artwork.
"""
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
OUTLINE = ROOT / "manuscript_outline"
DEST = OUTLINE / "revision_round2_20261022"
FIG = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022"


def read(path):
    return json.loads(path.read_text())


def relative(path):
    return str(path.relative_to(ROOT))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main():
    baseline = read(OUTLINE / "outline_v4.json")
    document = deepcopy(baseline)
    two = read(FIG / "figure2/legend.json")
    three = read(FIG / "figure3/content.json")
    four = read(FIG / "figure4/supporting/metadata/main_summary.json")
    assert "pending" not in two["status"].lower(), two["status"]
    assert four["case"] == "A" and four["post_hoc"]
    narrative = (
        "Enlarged cells resist imposed replacement through distinct initiation, spatial and temporal response signatures (Figure 2). "
        "Distance between competing regions matters only at short range; enlarged cells resist replacement and resolve coexistence slowly mainly because of size itself (Figure 3). "
        "Spring's middle size is the supported exception in which extra accessible separation adds a substantial resolution-time contribution. "
        + four["message"] + " (Figure 4). Geometry and multiple-front fractions are descriptive: erosion below the single-front geometric expectation can also occur while a single front persists. "
        "Observed migration summaries, conditional OU velocity-memory estimates and simulated polarity-heading readouts remain separate evidence."
    )
    document.update(
        date="22 October 2026",
        stage="Completed round-2 Figures 2–4, including approved Holmes follow-up",
        subtitle="Working manuscript outline with revised figures, complete retained supplements and explicit scientific qualifications",
    )
    document["overview"][1] = narrative
    for section in document["sections"]:
        ident = section["id"]
        if ident == "figure2":
            section.update(two)
            section["title"] = "Figure 2. " + two["title"].removeprefix("Figure 2. ")
            section["main_pdf"] = relative(FIG / "figure2/figure2_round2.pdf")
            section["supplement_pdf"] = relative(FIG / "figure2/figure2_round2_supplement.pdf")
            section["reader"] = relative(FIG / "figure2/README.md")
            section["role"] = "Distinguish replacement-resistance mechanisms using initiation, equal-area cue widths and measured pulse/gap retention."
            section["narrative"] = [two["legend"][0],
                "The targeted Holmes measurements report activity at every pulse and gap endpoint, independently of eventual replacement.",
                "No continuous commitment by 100,000 at Holmes L31, L35 and L39 is a finite-horizon finding of extreme replacement resistance."]
        elif ident == "figure3":
            section.update(three)
            section["narrative"] = [three["legend"][0],
                "At matched separation, size dominates both driven takeover and cue-free resolution. Extra separation contributes substantially in the supported Spring middle-size exception.",
                "The nine approved Holmes continuations extend the original exact endpoints to 80,000; unresolved conditions remain lower bounds and ratios of two censored times remain indeterminate.",
                "The Holmes large-size matched-gap lower bound is below one and does not establish a positive large-size effect; the full qualification is retained in the legend and decomposition CSV."]
            section["measurements"] = [
                "Measure complete-region continuous-cue handover against initial edge gap at fixed circumference.",
                "Withdraw at challenger/incumbent excess ratio 1 and measure first single-region resolution against the actual withdrawal gap.",
                "Decompose size at matched gap and extra separation using interpolation only within measured support; propagate censoring and saved-sample brackets."]
            with (FIG / "figure3/B_complete_with_extensions.csv").open(newline="") as stream:
                complete=list(csv.DictReader(stream))
            eligible=[r for r in complete if r["primary_claim_eligible"]=="True"]
            section["coverage"]["B_resolution_lower_bounds"]=sum(r["resolution_status"]=="right_censored_at_completed_horizon" for r in eligible)
            section["coverage"]["approved_Holmes_extensions"]=9
            section["coverage"]["Holmes_ratio1_censored"]=three["provenance"]["holmes_B_censored"]
            section["actual_horizons"]["Holmes_approved_B_extensions"]=dict(cases=9,cue_free_horizon=80000,original_horizon=20000)
            section["sources"].extend([
                dict(label="Current B evidence including all nine extensions",path=relative(FIG / "figure3/B_complete_with_extensions.csv")),
                dict(label="Matched-gap decomposition and censoring",path=relative(FIG / "figure3/decomposition_terms.csv")),
                dict(label="Measured short-range evidence",path=relative(FIG / "figure3/short_range_evidence.csv")),
                dict(label="Independent native-field and decomposition validation",path=relative(FIG / "figure3/verification.json")),
            ])
        elif ident == "figure4":
            section.update(
                title="Figure 4. " + four["message"],
                status="round2_complete_Case_A_posthoc",
                role="Measure the heading decay-rate advantage and describe geometric normalization and spatial organization under fluctuating input.",
                assets=[dict(path=relative(FIG / "figure4/figure4_main.png"), caption=four["message"],
                    page_title="Figure 4. Enlarged cells hold heading longer under fluctuating input",
                    page_size=[16.5, 21], minimum_printed_width_inches=14)],
                main_pdf=relative(FIG / "figure4/figure4_main.pdf"),
                supplement_pdf=relative(FIG / "figure4/figure4_supplement.pdf"),
                reader=relative(FIG / "figure4/index.html"),
                legend=[four["legend"]],
                narrative=[four["message"],
                    "Post hoc H1 passes in both main-view observation arms; frozen rate H1 remains inconclusive. All seven views, both arms and both forcing choices remain alongside in S0.",
                    "Geometry and multiple-front fractions are descriptive. The data do not establish multiple fronts as the cause of the declining advantage."],
                measurements=[
                    "Use endpoint coverage R ≥ 0.1 and a whole-seed bootstrap. Select the shortest usable lag at which small-cell autocorrelation is at most 0.98; report every usable lag and flag drift greater than 1.5-fold.",
                    "Measure fixed-reference-lag noise tolerance at autocorrelation 0.5 across every designed input, using log-input interpolation and explicit bounds.",
                    "Normalize the rate ratio by actual mean-circumference geometry, relative to the single-front expectation; this is not site diffusivity."],
                limitations=[
                    "The H1 domain restriction to inputs below the largest size's estimated noise tolerance is post hoc.",
                    "Weak-input geometric comparisons describe point estimates; Holmes and Spring confidence intervals also admit sub-geometric gains.",
                    "Non-exponential lags, forcing, off-grid lags and nonpositive autocorrelations remain explicit. Goryachev is complete in the supplement.",
                    "H2, H3 and A6 tau-based gates are retired; their frozen historical calls are preserved for provenance.",
                    "Simulated activity-heading persistence is distinct from observed migration and conditional OU persistence estimates."],
                equations=[
                    dict(math=r"C(\ell)=\langle\cos[\theta(t+\ell)-\theta(t)]\rangle,\qquad k(\ell)=-\frac{\ln C(\ell)}{\ell}",number="1",
                        description="Both endpoints must satisfy the declared magnitude rule. The logarithmic rate is defined only for positive autocorrelation."),
                    dict(math=r"\rho_k=\frac{k_{\mathrm{small}}}{k_{\mathrm{large}}},\qquad \rho_{\mathrm{geo}}=\frac{\rho_k}{(\mathcal{C}_{\mathrm{large}}/\mathcal{C}_{\mathrm{small}})^2}",number="2",
                        description="Calligraphic C denotes actual mean circumference. Geometric normalization is a relative single-front expectation."),
                    dict(math=r"C(\ell_{\mathrm{ref}};f_{1/2})=\frac{1}{2}",number="3",
                        description="The fixed reference lag is 2,000 for chemical models and 200 for Spring; crossing interpolation is logarithmic in input fraction.")],
            )
            section["sources"].extend([
                dict(label="Frozen and post hoc H1 decisions",path=relative(FIG / "figure4/supporting/data/S0_posthoc_H1.csv")),
                dict(label="Every lag and lag-consistency flags",path=relative(FIG / "figure4/supporting/data/all_lag_size_ratios.csv")),
                dict(label="Complete lag-by-lag ratio supplement",path=relative(FIG / "figure4/figure4_all_lag_ratios.pdf")),
                dict(label="Source and rendering receipt",path=relative(FIG / "figure4/supporting/metadata/main_render_receipt.json")),
            ])
            section["main_source_data"]=[dict(panel=panel,measurement=measurement,
                path=relative(FIG / "figure4/supporting/data" / name),
                sha256=hashlib.sha256((FIG / "figure4/supporting/data" / name).read_bytes()).hexdigest())
                for panel,measurement,name in [
                    ("A","Native f=0.003 paired activity histories and f=0.1 large Wave-pinning","main_native_kymographs.csv.gz"),
                    ("B","Heading autocorrelation and whole-seed confidence bands","main_autocorrelation.csv"),
                    ("C","Primary decay-rate size ratio against input","main_size_ratios.csv"),
                    ("D","Descriptive multiple-front fractions","main_multiple_fronts.csv"),
                    ("E","Fixed-reference-lag noise tolerance","main_noise_tolerance.csv")]]
            section["supplement_scope"]="New rate, geometry, fixed-lag tolerance, all-lag consistency and post hoc/frozen S0 analyses, followed by the complete retained 20 October supplement and superseded main artwork."
    for ident in ("figure1", "figure5", "figure6"):
        assert next(s for s in document["sections"] if s["id"] == ident) == next(s for s in baseline["sections"] if s["id"] == ident)
    for section in document["sections"]:
        if section["id"] not in ("figure2", "figure3", "figure4"):
            continue
        for asset in section["assets"]:
            width, height = Image.open(ROOT / asset["path"]).size
            page_width = max(asset["page_size"][0],asset["minimum_printed_width_inches"]+1.5)
            asset["page_size"] = [page_width,max(asset["page_size"][1],3.5+(page_width-1.4)*height/width)]
    write(OUTLINE / "outline_v5.json", document)
    writer = PdfWriter()
    for name in ("figure2/figure2_round2.pdf", "figure3/figure3_main.pdf", "figure4/figure4_main.pdf"):
        writer.append(FIG / name)
    writer.add_metadata({"/Title":"Completed round-2 Figures 2–4"})
    with (DEST / "figures_2_4.pdf").open("wb") as stream:
        writer.write(stream)
    generated = datetime.now(timezone.utc).isoformat()
    write(DEST / "completion_integration.json",dict(generated_utc=generated,
        final_outline_v5_integrated=True,protected_sections=["figure1","figure5","figure6"],
        cross_figure_narrative=narrative,figure4_case="A",figure4_H1_posthoc=True,
        figure2_status=two["status"],figure3_status=three["status"],source_brief=relative(ROOT / "docs/tasks/briefs/20261022_figures_round2.txt"),
        steering=relative(DEST / "urgent_user_steering_20261022T110301Z.txt")))
    body=["<!doctype html><html lang='en'><meta charset='utf-8'><title>Completed Figures 2–4 round 2</title>",
        "<style>body{max-width:1450px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}img{width:100%}a{color:#166aa2}</style>",
        "<h1>Completed manuscript Figures 2–4, round 2</h1>",
        "<p>Generated "+generated+". Simulation producers and figure generation complete.</p>",
        '<p><a href="../outline_v5.pdf">Integrated outline_v5</a> · <a href="figures_2_4.pdf">Three main plates</a> · <a href="revised_figures_2_4.pdf">Figures with legends</a></p>',
        "<p>"+html.escape(narrative)+"</p>"]
    for number in (2,3,4):
        section=next(s for s in document["sections"] if s["id"]==f"figure{number}")
        link=lambda path:html.escape(os.path.relpath(ROOT/path,DEST),quote=True)
        body.extend(["<h2>"+html.escape(section["title"])+"</h2>",
            '<p><a href="'+link(section["main_pdf"])+'">Main PDF</a> · <a href="'+link(section["supplement_pdf"])+'">Complete supplement</a> · <a href="'+link(relative(FIG/f"figure{number}/README.md"))+'">Changes, plotted values and qualifications</a></p>',
            '<img alt="Figure '+str(number)+'" src="'+link(section["assets"][0]["path"])+'">',
            "<p>"+html.escape(section["legend"][0])+"</p>"])
    body.append("</html>")
    (DEST / "completion.html").write_text("\n".join(body))
    print(OUTLINE / "outline_v5.json")


if __name__ == "__main__":
    main()
