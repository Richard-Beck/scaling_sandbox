#!/usr/bin/env python3
"""Render the manuscript mechanism crosswalk from its single retained JSON input.

No scientific data, source report, simulator or native history is opened. The
source JSON remains the authoritative text; wrapping changes whitespace only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
WIDTH = 10.0
FONT = 9.5
LEADING = 1.23
INK = "#193847"
TEAL = "#087c87"
MUTED = "#465c68"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wrap(text, width_inches, renderer, figure, size=FONT, weight="normal"):
    """Wrap only at word boundaries using actual output-font widths."""
    prop = FontProperties(family="DejaVu Sans", size=size, weight=weight)
    limit = width_inches * figure.dpi
    lines = []
    line = ""
    for word in text.split():
        trial = f"{line} {word}".strip()
        measured = renderer.get_text_width_height_descent(trial, prop, False)[0]
        if measured > limit and line:
            lines.append(line)
            line = word
        else:
            line = trial
    if line:
        lines.append(line)
    assert " ".join(lines) == " ".join(text.split())
    return "\n".join(lines), len(lines) * size * LEADING / 72


def draw_sheet(data, rows, sheet, stem, outdir, generated):
    # Measurement canvas uses the same DPI, font and physical width as output.
    probe = plt.figure(figsize=(WIDTH, 8), dpi=180)
    probe.canvas.draw()
    renderer = probe.canvas.get_renderer()
    margin = 0.38
    table_width = WIDTH - 2 * margin
    widths = [2.96, 3.14, table_width - 6.10]
    padding = 0.15
    title, title_h = wrap(data["title"], table_width, renderer, probe, 15, "bold")
    scope, scope_h = wrap(data["scope"], table_width, renderer, probe, 9.5)
    header_h = 0.50
    laid_rows = []
    for row in rows:
        mechanism, mh = wrap(row["candidate_mechanism"], widths[0] - 2 * padding,
                             renderer, probe)
        evidence, eh = wrap(row["evidence"], widths[1] - 2 * padding, renderer, probe)
        measurement, qh = wrap(row["measurement"], widths[2] - 2 * padding,
                               renderer, probe)
        height = max(0.24 + mh, eh, qh) + 2 * padding + 0.08
        laid_rows.append((row, mechanism, evidence, measurement, height))
    scopes = []
    for row in rows:
        text = row["model"] + ": " + row["inference_scope"]
        paragraph, ph = wrap(text, table_width - 2 * padding, renderer, probe, 9.5)
        scopes.append((paragraph, ph))
    scope_box_h = 0.37 + sum(ph + 0.09 for _, ph in scopes) + 0.12
    table_h = header_h + sum(row[-1] for row in laid_rows)
    top_h = 0.30 + 0.24 + title_h + 0.13 + scope_h + 0.25
    height = top_h + table_h + 0.23 + scope_box_h + 0.55
    plt.close(probe)

    fig = plt.figure(figsize=(WIDTH, height), dpi=180)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(height, 0)
    ax.axis("off")
    artists = []
    bound_checks = []

    def text(x, y, value, size=FONT, weight="normal", color=INK, box=None):
        art = ax.text(x, y, value, va="top", ha="left", fontsize=size,
                      fontfamily="DejaVu Sans", fontweight=weight, color=color,
                      linespacing=LEADING, parse_math=False)
        artists.append(art)
        if box is not None:
            bound_checks.append((art, box))
        return art

    y = 0.28
    text(margin, y, f"FIGURE 2 CONTINUED · EXISTING MECHANISTIC EVIDENCE · SHEET {sheet}/2",
         9.4, "bold", TEAL)
    y += 0.24
    text(margin, y, title, 15, "bold")
    y += title_h + 0.13
    text(margin, y, scope, 9.5, color=MUTED)
    y += scope_h + 0.25
    table_top = y
    ax.add_patch(Rectangle((margin, y), table_width, header_h, color=INK, zorder=0))
    x = margin
    for label, width in zip(["Model / candidate mechanism", "Observed evidence", "Corresponding measurements"], widths):
        text(x + padding, y + 0.13, label, 10.7, "bold", "white")
        x += width
    y += header_h
    for i, (row, mechanism, evidence, measurement, rh) in enumerate(laid_rows):
        ax.add_patch(Rectangle((margin, y), table_width, rh,
                              facecolor="#f1f6f7" if i % 2 == 0 else "#ffffff",
                              edgecolor="#cbdadd", linewidth=0.6, zorder=0))
        x = margin
        text(x + padding, y + padding, row["model"], 11.1, "bold", TEAL,
             box=(x + 0.05, y + 0.04, x + widths[0] - 0.05, y + rh - 0.04))
        text(x + padding, y + padding + 0.24, mechanism,
             box=(x + 0.05, y + 0.04, x + widths[0] - 0.05, y + rh - 0.04))
        x += widths[0]
        text(x + padding, y + padding, evidence,
             box=(x + 0.05, y + 0.04, x + widths[1] - 0.05, y + rh - 0.04))
        x += widths[1]
        text(x + padding, y + padding, measurement,
             box=(x + 0.05, y + 0.04, x + widths[2] - 0.05, y + rh - 0.04))
        y += rh
    for offset in [widths[0], widths[0] + widths[1]]:
        ax.plot([margin + offset] * 2, [table_top + header_h, y], color="#cbdadd", lw=0.6)
    y += 0.23
    ax.add_patch(Rectangle((margin, y), table_width, scope_box_h,
                          facecolor="#edf5f3", edgecolor="#aac7c2", linewidth=0.7))
    text(margin + padding, y + 0.12, "Inference scope", 11.1, "bold", TEAL)
    sy = y + 0.37
    for paragraph, ph in scopes:
        text(margin + padding, sy, paragraph, 9.5,
             box=(margin + 0.04, y + 0.04, WIDTH - margin - 0.04, y + scope_box_h - 0.02))
        sy += ph + 0.09
    footer_y = height - 0.38
    footer = "Source: mechanism_crosswalk_v2.json · Exact assays, gates and source citations accompany the full legend."
    text(margin, footer_y, footer, 9.0, color=MUTED)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    violations = []
    for art, box in bound_checks:
        bbox = art.get_window_extent(renderer).transformed(ax.transData.inverted())
        x0, x1 = min(bbox.x0, bbox.x1), max(bbox.x0, bbox.x1)
        y0, y1 = min(bbox.y0, bbox.y1), max(bbox.y0, bbox.y1)
        left, top, right, bottom = box
        if x0 < left - 0.005 or x1 > right + 0.005 or y0 < top - 0.005 or y1 > bottom + 0.005:
            violations.append({"text": art.get_text(), "bounds_inches": [x0, y0, x1, y1],
                               "allowed_inches": list(box)})
    assert not violations, violations
    # Text boxes in a column may share a containing cell; test actual bounding
    # boxes independently so model labels cannot collide with their paragraphs.
    actual_bounds = [a.get_window_extent(renderer) for a in artists]
    overlaps = []
    for i, a in enumerate(actual_bounds):
        for j in range(i):
            b = actual_bounds[j]
            if min(a.x1, b.x1) - max(a.x0, b.x0) > 1 and min(a.y1, b.y1) - max(a.y0, b.y0) > 1:
                overlaps.append([artists[j].get_text(), artists[i].get_text()])
    assert not overlaps, overlaps
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=220, facecolor="white")
    fig.savefig(pdf, facecolor="white", metadata={
        "Title": data["title"] + f" — sheet {sheet} of 2",
        "Subject": "Existing mechanism evidence and corresponding measurements",
        "Creator": "render_mechanism_crosswalk_v2.py",
    })
    plt.close(fig)
    return {"sheet": sheet, "models": [r["model"] for r in rows],
            "size_inches": [WIDTH, height], "minimum_font_points_at_10inch_width": 9.0,
            "text_bounds_violations": violations, "text_overlaps": overlaps,
            "png": {"path": str(png.relative_to(ROOT)), "sha256": digest(png)},
            "pdf": {"path": str(pdf.relative_to(ROOT)), "sha256": digest(pdf)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "manuscript_outline/mechanism_crosswalk_v2.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "manuscript_outline/assets")
    args = parser.parse_args()
    source = args.input.resolve()
    outdir = args.output_dir.resolve()
    data = json.loads(source.read_text())
    assert len(data["rows"]) == 5
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    outdir.mkdir(parents=True, exist_ok=True)
    generated = datetime.now(timezone.utc).isoformat()
    outputs = [draw_sheet(data, data["rows"][:2], 1, "mechanism_crosswalk_v2", outdir, generated),
               draw_sheet(data, data["rows"][2:], 2, "mechanism_crosswalk_v2_continued", outdir, generated)]
    assets = [{"path": r["png"]["path"],
               "caption": f"Figure 2 continued, mechanism-to-measurement companion, sheet {r['sheet']} of 2. "
                          + ", ".join(r["models"]) + ". Existing evidence and inference scope remain distinct."}
              for r in outputs]
    (outdir / "mechanism_crosswalk_v2.assets.json").write_text(json.dumps(assets, indent=2) + "\n")
    manifest = {"generated_utc": generated, "passed": True,
                "source_json": {"path": str(source.relative_to(ROOT)), "sha256": digest(source)},
                "renderer": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": digest(Path(__file__))},
                "scientific_source_reads": [str(source.relative_to(ROOT))],
                "new_integrations": 0, "all_five_models": True,
                "source_wording_preserved": True, "text_bounds_and_overlaps_checked": True,
                "sheets": outputs}
    (outdir / "mechanism_crosswalk_v2.provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (outdir / "mechanism_crosswalk_v2.validation.json").write_text(json.dumps({
        "generated_utc": generated, "passed": True, "sheets": len(outputs),
        "all_five_models": True, "no_overlapping_or_out_of_cell_text": True,
        "minimum_font_points_at_10inch_width": 9.0,
        "scope": "Exact source text and layout checks; scientific inference is reviewed separately."
    }, indent=2) + "\n")
    print(json.dumps({"passed": True, "assets": assets}, indent=2))


if __name__ == "__main__":
    main()
