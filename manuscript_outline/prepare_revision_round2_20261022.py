"""Integrate the round-2 Figure 4 review interim without composing a new main.

Run after the scalar-history analysis, S0 replay and interim rendering. The
20 October Figure 2/3 artwork remains pending its later revision. Protected
Figures 1, 5 and 6 are copied exactly. No model integrations are launched.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
OUTLINE = ROOT / "manuscript_outline"
FIGURE = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4"
DEST = OUTLINE / "revision_round2_20261022"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def main():
    baseline_path = OUTLINE / "outline_v4.json"
    baseline = json.loads(baseline_path.read_text())
    summary_path = FIGURE / "supporting/metadata/interim_summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["case"] in ("A", "B"), summary
    analysis_path = FIGURE / "supporting/metadata/analysis_receipt.json"
    audit_path = FIGURE / "supporting/metadata/S0_decision_replay_audit.json"
    analysis = json.loads(analysis_path.read_text())
    audit = json.loads(audit_path.read_text())
    assert analysis["existing_outputs_only"] and analysis["no_simulations"]
    # S0 renderer validates all original calls before replacing either estimand.
    assert audit.get("passed") is True, audit
    supplement_parts = [
        FIGURE / "figure4_S0_round2_old_new.pdf",
        FIGURE / "supporting/vectors/figure4_reanalysis_supplement.pdf",
        FIGURE / "supporting/vectors/figure4_old_new.pdf",
        FIGURE / "supporting/vectors/legacy_seven_view_robustness.pdf",
    ]
    DEST.mkdir(parents=True, exist_ok=True)
    supplement = DEST / "figure4_interim_complete_supplement.pdf"
    writer = PdfWriter()
    for part in supplement_parts:
        writer.append(part)
    writer.add_metadata({"/Title": "Figure 4 round-2 review interim: S0, rate estimands and retained robustness"})
    with supplement.open("wb") as stream:
        writer.write(stream)
    document = deepcopy(baseline)
    document.update(date="22 October 2026",
                    stage="Round-2 Figure 4 interim — awaiting author review; later revisions pending",
                    subtitle="Working manuscript outline with the round-2 persistence reanalysis and explicit review stop")
    narrative = ("Enlarged cells resist imposed replacement through model-specific response signatures (Figure 2). "
                 "Larger rings can support more widely separated active regions, whose cue-free resolution can be slow in mass-conserved models (Figure 3). "
                 + summary["summary"] + " These are polarity-heading readouts; observed migration summaries and conditional OU velocity-memory estimates remain separate evidence.")
    document["overview"][1] = narrative
    document["overview"].append(
        "Round-2 delivery stops at the Figure 4 reanalysis, as requested in the brief. "
        "The new estimand, geometric normalization and old/new S0 tables are integrated as review-interim supplementary sheets. "
        "The Figure 4 main is not recomposed; Figures 2 and 3 retain their 20 October artwork pending review and the subsequent work package. "
        "No Holmes retention runs or optional Figure 3 horizon extensions have been launched.")
    section = next(s for s in document["sections"] if s["id"] == "figure4")
    section["status"] = "round2_interim_pending_author_review"
    section["title"] = "Figure 4. Heading persistence under fluctuating input — round-2 reanalysis pending review"
    section["role"] = ("Review the rate-based persistence estimand before selecting the revised main composition. "
                       "The following main artwork is retained from 20 October; the new interim sheets follow in the supplementary section.")
    for asset in section["assets"]:
        asset["caption"] = "20 October main artwork retained pending round-2 review; see the new interim sheets"
        asset["page_title"] = "Figure 4. Previous main artwork — round-2 composition awaits review"
    section["legend"].insert(0, summary["summary"] + " The retained artwork below uses the previous estimand; new rate-based values and decisions appear in the round-2 interim supplement.")
    section["narrative"] = [summary["summary"],
        "Geometric normalization is relative to the single-front geometric expectation; it is not a site diffusivity.",
        "S0 substitutes the new size ratio and fixed-lag noise tolerance. Its explicitly direct-tau A6 predictor gate and H3 fingerprints remain historical, as recorded in the decision audit."]
    section["limitations"] = list(summary.get("qualifications", [])) + section["limitations"]
    section["limitations"].append("Author review is required before Figure 4 main composition and the later Figure 3/2 revisions. The complete round-2 work package is not yet delivered.")
    section["sources"].extend([
        {"label": "Combined round-2 interim supplement", "path": relative(supplement)},
        {"label": "Round-2 interim summary and qualifications", "path": relative(summary_path)},
        {"label": "Rate estimand source and scalar-history receipt", "path": relative(analysis_path)},
        {"label": "Exact old-S0 replay and new-estimand audit", "path": relative(audit_path)},
    ])
    previews = sorted(FIGURE.glob("*.png"))
    assert previews, "No rendered review-interim sheets"
    interim = dict(id="figure4_round2_interim", title="Figure 4 round-2 review interim",
                   role="New persistence estimand and old/new decision table; these sheets are review material, not a recomposed main figure.",
                   assets=[dict(path=relative(p), caption="Round-2 review interim: " + p.stem.replace("_", " "),
                                page_size=[14.4, 12.2] if p.stem.startswith("S0_") else [16.4, 8.0],
                                minimum_printed_width_inches=12.5 if p.stem.startswith("S0_") else 14.5,
                                page_title="Figure 4 round-2 review interim") for p in previews],
                   legend=[summary["summary"], "Model colours are fixed across the new sheets. Bounds use open triangles and arrows; unavailable values carry explicit labels. Native units are retained.",
                           "The S0 tables pair old and revised calls across both observation arms, all seven views and both forcing selections. H3 and the A6 predictor gate retain their stated historical tau definitions."],
                   narrative=[summary["summary"]], measurements=[],
                   limitations=summary.get("qualifications", []) + summary.get("pending", []),
                   sources=[{"label": "Complete interim reader", "path": relative(FIGURE / "index.html")},
                            {"label": "Scientific methods and qualifications", "path": relative(FIGURE / "supporting/metadata/scientific_readme.md")},
                            {"label": "Lag-by-lag condition measurements", "path": relative(FIGURE / "supporting/data/condition_lag_metrics.csv")},
                            {"label": "Primary size ratios and geometric normalization", "path": relative(FIGURE / "supporting/data/size_ratios.csv")},
                            {"label": "Fixed-reference-lag noise tolerance", "path": relative(FIGURE / "supporting/data/noise_tolerance.csv")}])
    document.setdefault("supplementary", []).append(interim)
    protected = ("figure1", "figure5", "figure6")
    for ident in protected:
        assert next(s for s in document["sections"] if s["id"] == ident) == next(s for s in baseline["sections"] if s["id"] == ident)
    write(OUTLINE / "outline_v5.json", document)
    generated = datetime.now(timezone.utc).isoformat()
    # Measure this new output group without relabelling older storage inventory.
    inodes, allocated = set(), 0
    for path in FIGURE.rglob("*"):
        if path.is_file() and not path.is_symlink():
            stat = path.stat()
            inode = (stat.st_dev, stat.st_ino)
            if inode not in inodes:
                inodes.add(inode)
                allocated += stat.st_blocks * 512
    record = dict(generated_utc=generated, status="interim_pending_author_review", no_simulations_launched=True,
                  generator_sha256=sha(Path(__file__)),
                  new_output_storage=dict(path=relative(FIGURE), measured_utc=generated,
                                          allocated_bytes=allocated, allocated_GiB=allocated / 2**30,
                                          native_inputs_modified=False),
                  source_outline={"path": relative(baseline_path), "sha256": sha(baseline_path)},
                  protected_sections_verified=list(protected), previous_artwork_retained=["figure2", "figure3", "figure4"],
                  new_main_composition=False, case=summary["case"], cross_figure_narrative=narrative,
                  inputs=[dict(path=relative(p), sha256=sha(p)) for p in [summary_path, analysis_path, audit_path, *previews]],
                  combined_supplement=dict(path=relative(supplement), sha256=sha(supplement),
                                           pages=len(PdfReader(supplement).pages),
                                           parts=[dict(path=relative(p), sha256=sha(p)) for p in supplement_parts]),
                  outline={"path": "manuscript_outline/outline_v5.json", "sha256": sha(OUTLINE / "outline_v5.json")})
    write(DEST / "integration.json", record)
    reader = DEST / "index.html"
    href = lambda p: html.escape(os.path.relpath(p, reader.parent), quote=True)
    lines = ['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
             '<title>Round-2 Figure 4 interim</title><style>body{max-width:1280px;margin:2rem auto;padding:0 1rem;font:17px/1.55 system-ui;color:#203348}img{width:100%;height:auto}a{color:#126478}</style>',
             '<h1>Round-2 Figure 4 interim — author review required</h1>',
             '<p>Generated ' + html.escape(generated) + '. Static interim; analysis producer has finished. No timed refresh.</p>',
             '<p>' + html.escape(summary["summary"]) + '</p>',
             '<p><a href="../outline_v5.pdf">Integrated outline_v5</a> · <a href="figure4_interim_complete_supplement.pdf">Combined interim supplement</a> · <a href="' + href(FIGURE / 'index.html') + '">Full interim reader and CSVs</a></p>',
             '<p>Figures 2–3 and the Figure 4 main retain their previous artwork pending review. Figures 1, 5 and 6 are unchanged.</p>']
    for p in previews:
        lines += ['<figure><img src="' + href(p) + '" alt="' + html.escape(p.stem, quote=True) + '"></figure>']
    lines += ['<h2>Cross-figure narrative</h2><p>' + html.escape(narrative) + '</p>', '</html>']
    reader.write_text("\n".join(lines) + "\n")
    print(json.dumps({"outline": "manuscript_outline/outline_v5.json", "case": summary["case"], "review_interim_sheets": len(previews), "protected_sections_verified": list(protected)}))


if __name__ == "__main__":
    main()
