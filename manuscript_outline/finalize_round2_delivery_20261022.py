"""Refresh manuscript navigation and measured storage after full acceptance."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/"manuscript_outline/revision_round2_20261022"
MARKER="## Round-2 completion retained outputs (22 October 2026)"


def main():
    acceptance=json.loads((DEST/"completion_acceptance.json").read_text())
    assert acceptance["passed"] and acceptance["all_required_work_complete"]
    generated=datetime.now(timezone.utc).isoformat()
    narrative=json.loads((DEST/"completion_integration.json").read_text())["cross_figure_narrative"]
    notes=(
        "Figures 2–4 are complete, including the 32 Holmes train/noise histories, all nine approved Figure 3 extensions, "
        "the corrected short-range interpretation and the Figure 4 Case A main. Six Figure 3 extensions resolve; three remain censored at 80,000. "
        "Holmes class membership is unchanged. Numerical and finite-horizon qualifications remain explicit. "
        "Figure 4 H1 is explicitly **post hoc**, with frozen calls alongside; geometry and multiple-front fractions are descriptive. "
        "Figures 1, 5 and 6 retain their exact preceding content and artwork."
    )
    for path,prefix in [(ROOT/"manuscript_outline/README.md",""),(ROOT/"current manuscript figs/README.md","../manuscript_outline/")]:
        old=path.read_text()
        if old.startswith("# Completed Figures 2–4 round 2"):
            old=old.split("<!-- preserved earlier editions -->\n",1)[1]
        old=old.replace("# Figures 2–4 round-2 follow-up (22 October 2026)","# Historical 7 a.m. round-2 follow-up (22 October 2026)",1)
        old=old.replace("# Current manuscript figure copies","# Historical 20 October manuscript figure copies",1)
        # Historical interim links must continue to open the archived interim.
        if prefix:
            old=old.replace("](../manuscript_outline/outline_v5.pdf)","](../manuscript_outline/revision_round2_20261022/historical_interim/outline_v5.pdf)")
        else:
            old=old.replace("](outline_v5.pdf)","](revision_round2_20261022/historical_interim/outline_v5.pdf)")
        intro=("# Completed Figures 2–4 round 2 (22 October 2026)\n\n"
            f"[**Integrated outline_v5**]({prefix}outline_v5.pdf) · "
            f"[**Figures with legends**]({prefix}revision_round2_20261022/revised_figures_2_4.pdf) · "
            f"[**Three main plates**]({prefix}revision_round2_20261022/figures_2_4.pdf) · "
            f"[**Reader and complete supplements**]({prefix}revision_round2_20261022/completion.html).\n\n"
            +notes+"\n\n"+narrative+"\n\n"
            f"[Requirement-by-requirement acceptance]({prefix}revision_round2_20261022/completion_acceptance.json). "
            f"Generated {generated}. Earlier editions and their dated status notes are preserved below.\n\n---\n\n<!-- preserved earlier editions -->\n")
        path.write_text(intro+old)
    readme=("# Completed round-2 manuscript Figures 2–4\n\n"
        "[Integrated outline_v5](../outline_v5.pdf) · [Figures with legends](revised_figures_2_4.pdf) · "
        "[Three main plates](figures_2_4.pdf) · [Reader and supplements](completion.html) · "
        "[Acceptance evidence](completion_acceptance.json).\n\n"+notes+"\n\n"+narrative+"\n\n"
        "The scientific producers have stopped and every required trajectory is accounted for. "
        "Original evidence and complete supplements are retained; raw native arrays and the whole-seed bootstrap archive remain on HPC. "
        "The previous interim and in-progress review remain historical artifacts.\n\n"
        "Rebuild with the existing scaling_sandbox environment:\n\n"
        "```sh\npython manuscript_outline/prepare_round2_completion_20261022.py\n"
        "python scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v5.json --output manuscript_outline/outline_v5.pdf\n"
        "python scripts/verify_manuscript_outline.py manuscript_outline/outline_v5.json --report manuscript_outline/outline_v5.validation.json\n"
        "python manuscript_outline/verify_round2_completion_20261022.py\n```\n\n"
        f"Generated UTC: {generated}. Commit and remote verification are recorded in completion_delivery.json after publication.\n")
    (DEST/"COMPLETION_README.md").write_text(readme)
    # Keep the directory's default entry point current while retaining the
    # earlier interim text as explicitly dated history.
    directory_readme=DEST/"README.md"
    historical=directory_readme.read_text()
    marker="<!-- preserved interim directory notes -->\n"
    if marker in historical:
        historical=historical.split(marker,1)[1]
    historical=historical.replace("](../outline_v5.pdf)","](historical_interim/outline_v5.pdf)")
    historical=historical.replace("](../outline_v5.validation.json)","](historical_interim/outline_v5.validation.json)")
    directory_readme.write_text(readme+"\n---\n\n## Historical interim notes\n\n"+marker+historical)
    # Close the live status page only after the complete scientific audit.
    progress=dict(generated_utc=generated,producer_state="stopped",completed_workers=32,total_workers=32,
                  final_integration="complete",reader="completion.html",acceptance="completion_acceptance.json")
    (DEST/"completion_progress.json").write_text(json.dumps(progress,indent=2)+"\n")
    (DEST/"completion_progress.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Round-2 work complete</title>'
        '<style>body{max-width:1000px;margin:3rem auto;font:18px/1.5 system-ui}</style>'
        '<h1>Figures 2–4 complete</h1><p>Data generated '+generated+'. All scientific producers have stopped; all 32 Holmes histories and all nine Figure 3 extensions are accounted for.</p>'
        '<p><a href="completion.html">Open the final figures, integrated outline and complete supplements.</a></p>')
    groups=[
        ("polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/",
         "Completed and independently audited native evidence: 18 Holmes L31–L39 four-pulse trains, 12 L39 bounded-input seeds and two controls. Exact checkpoints, frozen sources, all native phase fields and joined trajectories are scientifically authoritative; retain these raw outputs."),
        ("polarity1d/manuscript/data/runs/signature_figure3_holmes_extensions_round2_20261022/",
         "Nine approved exact-endpoint Holmes cue-free continuations through 80,000; six resolutions and three lower bounds, with native fields, conservation and resolution audits. Required evidence for the revised decomposition."),
        ("polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/",
         "Completed derived Figure 2–4 artwork, complete retained supplements, plotted-value CSVs and scientific/display acceptance. Includes the HPC-only 164 MiB whole-seed bootstrap archive; do not mistake derived presentation for the underlying native trajectories."),
        ("manuscript_outline/revision_round2_20261022/",
         "Final manuscript reader, main plates and legends, exact-source acceptance, archived review interim and historical recovery bundles. Final integration is complete; previous in-progress artifacts are retained history."),
    ]
    allocated=[]
    for path,purpose in groups:
        value=int(subprocess.check_output(["du","-sB1",str(ROOT/path)],text=True).split()[0])
        allocated.append(dict(path=path,allocated_bytes=value,allocated_GiB=value/(1<<30),purpose=purpose,status="complete; producers stopped"))
    storage=dict(measured_utc=datetime.now(timezone.utc).isoformat(),method="du -sB1 allocated blocks; group totals may include small receipts written immediately after measurement",groups=allocated,scientific_producers_stopped=True)
    (DEST/"completion_storage.json").write_text(json.dumps(storage,indent=2)+"\n")
    manifest=ROOT/"DATA_MANIFEST.md";text=manifest.read_text()
    if MARKER in text:text=text.split(MARKER,1)[0].rstrip()+"\n"
    text+="\n"+MARKER+"\n\nMeasured "+storage["measured_utc"]+" using allocated blocks (`du -sB1`), after the scientific producers stopped. "
    text+="This updates the round-2 groups below; the older whole-workspace total remains an explicitly dated historical measurement. Original inputs were preserved and no scientific data were deleted.\n\n"
    text+="| Root-relative group | Allocated storage | Purpose and scientific status |\n| --- | ---: | --- |\n"
    for row in allocated:
        text+=f"| `{row['path']}` | {row['allocated_GiB']:.4f} GiB ({row['allocated_bytes']:,} allocated bytes) | {row['purpose']} |\n"
    manifest.write_text(text)
    print(json.dumps(dict(generated_utc=generated,storage_groups=allocated)))


if __name__=="__main__":
    main()
