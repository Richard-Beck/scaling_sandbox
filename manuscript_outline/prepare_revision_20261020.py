#!/usr/bin/env python3
"""Assemble the completed October 20 revisions and preserve Figures 1, 5 and 6.

This script reads only reviewed figure exports. It does not launch simulations,
choose a Figure 4 framing, or infer scientific completion from scheduler status.
The integration plan must contain the exact reviewed identities for all three
figures before any report is written.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUTLINE = ROOT / "manuscript_outline"
PLAN = ROOT / "docs/tasks/briefs/20261020_outline_integration_plan.json"


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def checked(record):
    path = (ROOT / record["path"]).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise ValueError("Missing repository export: " + str(path))
    if digest(path) != record["sha256"]:
        raise ValueError("Reviewed export changed: " + relative(path))
    return path


def section_digest(section):
    return hashlib.sha256(json.dumps(section, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":")).encode()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def load_reviewed_sections(plan):
    records = plan.get("reviewed_figures", {})
    if set(records) != {"2", "3", "4"}:
        raise ValueError("All three completed and independently reviewed figures are required")
    sections, receipts = {}, {}
    for number in (2, 3, 4):
        record = records[str(number)]
        content_path = checked(record["content"])
        review_path = checked(record["completion_review"])
        review = load(review_path)
        if review.get("independent_review_passed") is not True:
            raise ValueError(f"Figure {number} independent completion review has not passed")
        section = load(content_path)
        section = copy.deepcopy(section.get("section", section))
        if section.get("id") != f"figure{number}" or section.get("status") != "final":
            raise ValueError(f"Figure {number} is not a final reviewed section")
        for field in ("title", "role", "assets", "legend", "narrative", "measurements", "limitations", "sources"):
            if not section.get(field):
                raise ValueError(f"Figure {number} lacks {field}")
        for field in ("main_pdf", "supplement_pdf", "reader"):
            target = checked(record[field])
            if field in section and relative(target) != section[field]:
                raise ValueError(f"Figure {number} {field} disagrees with its reviewed identity")
            section.setdefault(field, relative(target))
        volumes = section.get("supplement_volumes", [])
        reviewed_volumes = record.get("supplement_volumes", [])
        if volumes != reviewed_volumes:
            raise ValueError(f"Figure {number} supplement volumes differ from the reviewed complete set")
        if number == 4 and not volumes:
            raise ValueError("Figure 4 requires its complete native-atlas appendix volumes")
        volume_paths = [volume["path"] for volume in volumes]
        if len(set(volume_paths)) != len(volume_paths) or section["supplement_pdf"] in volume_paths:
            raise ValueError(f"Figure {number} supplement volumes are duplicated")
        for volume in volumes:
            checked(volume)
            if not volume.get("label") or not isinstance(volume.get("pages"), int) or volume["pages"] <= 0:
                raise ValueError(f"Figure {number} supplement volume lacks a label or page count")
        reviewed_assets = {row["path"]: row for row in record["assets"]}
        if set(reviewed_assets) != {asset["path"] for asset in section["assets"]}:
            raise ValueError(f"Figure {number} asset coverage differs from reviewed content")
        for asset in section["assets"]:
            path = checked(reviewed_assets[asset["path"]])
            with Image.open(path) as image:
                width, height = image.size
            asset.setdefault("minimum_printed_width_inches", 10.)
            asset.setdefault("page_title", section["title"])
            # The outline adds its own header, source caption and margins.
            # Native figure sheet dimensions therefore need extra space here;
            # preserve every pixel and the requested minimum printed width.
            original_size = asset.get("page_size", [11.8, 8.])
            minimum = asset["minimum_printed_width_inches"]
            page_width = max(original_size[0], minimum + 1.5)
            page_height = max(original_size[1], 3.5 + (page_width - 1.4) * height / width)
            if page_width > 24 or page_height > 24:
                raise ValueError("Figure needs continued readable outline sheets: " + asset["path"])
            asset["page_size"] = [page_width, page_height]
        for source in section["sources"]:
            if source.get("path") and not (ROOT / source["path"]).is_file():
                raise ValueError("Missing figure source: " + source["path"])
        sections[number] = section
        receipts[str(number)] = copy.deepcopy(record)
    return sections, receipts


def assemble(plan_path, check_only=False):
    plan = load(plan_path)
    baseline_path = checked({"path": plan["baseline"]["path"], "sha256": plan["baseline"]["sha256"]})
    baseline = load(baseline_path)
    sections, receipts = load_reviewed_sections(plan)
    for number in (1, 5, 6):
        original = baseline["sections"][number - 1]
        expected = plan["baseline"]["protected_section_sha256"][f"figure{number}"]
        if section_digest(original) != expected:
            raise ValueError(f"Protected Figure {number} baseline changed")
    if check_only:
        return {"passed": True, "written": False, "reviewed_figures": [2, 3, 4]}

    document = copy.deepcopy(baseline)
    document.update(date="20 October 2026", stage="Completed Figures 2–4 revisions with stated scientific qualifications")
    document["overview"] = [
        baseline["overview"][0],
        "Figures 2–4 test how challenge shape and timing separate response mechanisms, how physical separation affects competition at fixed size, and how replicated fluctuating input changes polarity-direction persistence and spatial organization. Figure 2 reports equal-area widths and pulse trains. Figure 3 distinguishes cue-driven handover, cue-free resolution and directly measured activity direction. Figure 4 follows its completed hypothesis, coverage and sensitivity decisions. Detailed native histories and every finite-horizon or numerical qualification remain in the linked supplements.",
        baseline["overview"][-1],
    ]
    additions = []
    for number, section in sections.items():
        document["sections"][number - 1] = section
        volumes = section.get("supplement_volumes", [])
        supplement_text = "The comprehensive supplement retains the full designed condition grid, unavailable measurements, finite-horizon bounds, native histories, controls and the stated numerical and observational sensitivities. Earlier figure versions remain development history."
        if volumes:
            supplement_text += " " + section["supplement_scope"] + ". Every numbered appendix below is part of the complete supplement."
        additions.append({
            "title": f"Figure {number}: complete revised evidence",
            "text": [supplement_text],
            "path": receipts[str(number)]["supplement_pdf"]["path"],
            "sources": [
                {"label": "Figure reader and source data", "path": receipts[str(number)]["reader"]["path"]},
                {"label": "Independent completion review", "path": receipts[str(number)]["completion_review"]["path"]},
                *[{"label": volume["label"], "path": volume["path"]} for volume in volumes],
            ],
        })
    for number in (1, 5, 6):
        if document["sections"][number - 1] != baseline["sections"][number - 1]:
            raise ValueError(f"Protected Figure {number} changed during assembly")
    # Existing supplementary sections are preserved as explicitly labelled history.
    historical = copy.deepcopy(baseline.get("supplementary", []))
    for item in historical:
        if isinstance(item, dict) and item.get("title"):
            item["title"] = "Earlier supporting material: " + item["title"]
    document["supplementary"] = additions + historical
    document.setdefault("references", []).extend(
        {"label": f"Figure {n} October 20 revision brief", "path": plan[f"figure{n}"]["brief"]}
        for n in (2, 3, 4))
    write_json(OUTLINE / "outline_v4.json", document)
    review = copy.deepcopy(document)
    review.update(title="Revised Figures 2–4: challenge signatures, separation and noisy polarity", subtitle="Completed revisions, detailed legends and linked comprehensive evidence")
    review["sections"] = document["sections"][1:4]
    review["supplementary"] = additions
    write_json(OUTLINE / "revision_20261020/revised_figures_2_4.json", review)
    generated = datetime.now(timezone.utc).isoformat()
    record = dict(generated_utc=generated, integration_launches_simulations=False,
                  source_outline={"path": relative(baseline_path), "sha256": digest(baseline_path)},
                  outline={"path": "manuscript_outline/outline_v4.json", "sha256": digest(OUTLINE / "outline_v4.json")},
                  integration_plan={"path": relative(plan_path), "sha256": digest(plan_path)},
                  generator_sha256=digest(Path(__file__)), reviewed_figures=receipts,
                  protected_sections_verified=["figure1", "figure5", "figure6"])
    write_json(OUTLINE / "revision_20261020/integration.json", record)
    reader = OUTLINE / "revision_20261020/index.html"
    href = lambda path: html.escape(os.path.relpath(ROOT / path, reader.parent), quote=True)
    lines = ["<!doctype html><html lang=\"en\"><meta charset=\"utf-8\">",
             '<meta name="viewport" content="width=device-width, initial-scale=1">',
             "<title>Revised Figures 2–4</title>",
             "<style>body{max-width:1344px;margin:2rem auto;padding:0 1rem;font:17px/1.55 system-ui;color:#203348}img{max-width:100%;height:auto}.figure-scroll{overflow-x:auto}section{margin:3rem 0}a{color:#126478}</style>",
             "<h1>Revised Figures 2–4</h1>",
             f"<p>Generated {html.escape(generated)}. Static report assembled from completed, independently reviewed figure exports. Scientific and numerical qualifications accompany each figure.</p>",
             f'<p><a href="{href("manuscript_outline/outline_v4.pdf")}">Integrated manuscript outline</a> · <a href="revised_figures_2_4.pdf">Figures and legends</a></p>']
    for number, section in sections.items():
        rec = receipts[str(number)]
        lines.extend([f'<section id="figure{number}"><h2>{html.escape(section["title"])}</h2>',
                      f'<p>{html.escape(section["role"])}</p>'])
        for paragraph in section["narrative"]:
            lines.append("<p>" + html.escape(paragraph) + "</p>")
        for asset in section["assets"]:
            minimum_px = int(asset.get("minimum_display_width_px", 0))
            image_style = f' style="min-width:{minimum_px}px"' if minimum_px else ""
            lines.append(f'<div class="figure-scroll"><a href="{href(rec["main_pdf"]["path"])}"><img{image_style} src="{href(asset["path"])}" alt="{html.escape(asset.get("caption", section["title"]), quote=True)}"></a></div>')
        supplement_label = "Supplement statistical core and atlas index" if section.get("supplement_volumes") else "Complete supplement"
        lines.append(f'<p><a href="{href(rec["supplement_pdf"]["path"])}">{supplement_label}</a> · <a href="{href(rec["reader"]["path"])}">Figure reader and data</a> · <a href="{href(rec["completion_review"]["path"])}">Independent review</a></p>')
        if section.get("supplement_volumes"):
            lines.append('<p>' + html.escape(section["supplement_scope"]) + '.</p><ul>')
            for volume in section["supplement_volumes"]:
                lines.append(f'<li><a href="{href(volume["path"])}">{html.escape(volume["label"])}</a> ({volume["pages"]} pages)</li>')
            lines.append('</ul>')
        lines.append("<details><summary>Detailed legend and limitations</summary>")
        for paragraph in section["legend"] + section["limitations"]:
            lines.append("<p>" + html.escape(paragraph) + "</p>")
        lines.append("</details></section>")
    lines.append("</html>")
    reader.write_text("\n".join(lines) + "\n")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=PLAN)
    parser.add_argument("--check", action="store_true", help="Verify the complete reviewed intake without writing reports")
    args = parser.parse_args()
    print(json.dumps(assemble(args.plan.resolve(), args.check), indent=2))
