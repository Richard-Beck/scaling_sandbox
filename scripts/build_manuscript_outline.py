#!/usr/bin/env python3
"""Build an offline manuscript outline PDF from a JSON document.

Usage::

    python scripts/build_manuscript_outline.py --input outline.json \
        --output "current manuscript figs/manuscript_outline.pdf"

Dependencies: matplotlib, Pillow and pypdf (the existing scaling_sandbox Conda
environment supplies these). No browser, network, TeX installation or simulation
is used. Run --example to print the input schema. All local input paths are
relative to --repo-root, which defaults to this script's parent repository.

Figures are displayed without cropping on dedicated pages. Prose is wrapped
using measured font widths and automatically continued on subsequent pages.
The PDF contains searchable text and clickable, URL-encoded relative links to
the source files. An adjacent build manifest records the exact input/source
hashes, image dimensions, page count and layout validation results.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import quote, urlparse

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans", "pdf.fonttype": 42,
    "mathtext.fontset": "dejavusans", "axes.unicode_minus": False,
})
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf.annotations import Link


INK = "#172d40"
ACCENT = "#147d83"
MUTED = "#536575"
LIGHT = "#d9e3e8"
PAGE = (8.27, 11.69)
LANDSCAPE = PAGE[::-1]
EXAMPLE = {
    "title": "Manuscript outline",
    "subtitle": "Scientific argument and proposed main figures",
    "date": "2026-10-17",
    "stage": "Checkpoint integration",
    "overview": ["A paragraph explaining the scientific argument."],
    "sections": [{
        "id": "figure_1", "title": "Figure 1. A precise scientific claim",
        "role": "Main manuscript figure",
        "assets": [{"path": "figures/example.png", "caption": "Original figure, panels A–D."}],
        "legend": ["(A) A detailed explanation of the first panel.", "(B–D) Remaining panels."],
        "narrative": ["The inference supported by the figure."],
        "measurements": ["The proposed experiment and its decision criterion."],
        "limitations": ["The boundary of the inference."],
        "equations": [{"math": "\\tau_p = k^{-1}", "number": "1",
                       "description": "Display equations are centered and numbered."}],
        "sources": [{"label": "Figure source", "path": "figures/example.png"}],
    }],
    "supplementary": ["A description of proposed supplementary material."],
    "references": [{"label": "A bibliographic reference", "url": "https://example.org"}],
}


def paragraphs(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [p.strip() for p in value.split("\n\n") if p.strip()]
    if isinstance(value, list):
        return [str(p).strip() for p in value if str(p).strip()]
    raise TypeError(f"Expected a string or list of strings, received {type(value).__name__}")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


class OutlineBuilder:
    def __init__(self, doc, repo, output, input_path):
        self.doc, self.repo, self.output = doc, repo, output
        self.generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        self.pages, self.links, self.files = [], [], {}
        self.fig = self.canvas = None
        self.pdf = None
        self.page_links, self.bounds = [], []
        self.w, self.h = PAGE
        self.y, self.bottom = 0, 0
        self.section_id = None
        self.page_title = None
        self.track(input_path, "document")
        self.track(Path(__file__), "renderer")

    def resolve(self, name):
        path = Path(name)
        return (self.repo / path).resolve() if not path.is_absolute() else path.resolve()

    def track(self, path, role):
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Missing {role}: {path}")
        rel = os.path.relpath(path, self.repo)
        if rel not in self.files:
            self.files[rel] = {"path": rel, "bytes": path.stat().st_size,
                               "sha256": digest(path), "roles": []}
        if role not in self.files[rel]["roles"]:
            self.files[rel]["roles"].append(role)
        return rel

    def link_for(self, source):
        if urlparse(str(source)).scheme in {"https", "http", "doi"}:
            return str(source)
        return quote(os.path.relpath(self.resolve(source), self.output.parent), safe="/.-_~")

    def width(self, text, size, weight="normal"):
        prop = FontProperties(family="DejaVu Sans", size=size, weight=weight)
        # Matplotlib's PDF backend uses the same font metrics. The slight wrap
        # allowance covers backend rounding without shrinking the printed text.
        return self.canvas.get_renderer().get_text_width_height_descent(
            text, prop, ismath=bool(re.search(r"(?<!\\)\$", text)))[0] * 72 / self.fig.dpi

    def wrap(self, text, size=10.5, width=None, weight="normal"):
        available = (width or self.w - 1.4) * 72 - 1.5
        result = []
        # Keep complete math expressions together so a line never has unmatched $.
        for raw_line in str(text).splitlines() or [""]:
            tokens = re.findall(r"\$(?:\\.|[^$])*\$|\S+", raw_line)
            line = ""
            for token in tokens:
                trial = f"{line} {token}".strip()
                if self.width(trial, size, weight) <= available:
                    line = trial
                    continue
                if line:
                    result.append(line)
                    line = ""
                if self.width(token, size, weight) <= available:
                    line = token
                    continue
                if token.startswith("$") and token.endswith("$"):
                    raise ValueError(f"Equation exceeds printable width; shorten it: {token}")
                # Long source paths and URLs may break at any character.
                fragment = ""
                for char in token:
                    if fragment and self.width(fragment + char, size, weight) > available:
                        result.append(fragment)
                        fragment = ""
                    fragment += char
                line = fragment
            result.append(line)
        return result

    def text(self, text, x, y, size=10.5, weight="normal", color=INK, url=None):
        artist = self.fig.text(x / self.w, 1 - y / self.h, text, fontsize=size,
                               weight=weight, color=color, ha="left", va="top")
        box = artist.get_window_extent(self.canvas.get_renderer())
        point_box = [box.x0 * 72 / self.fig.dpi, box.y0 * 72 / self.fig.dpi,
                     box.x1 * 72 / self.fig.dpi, box.y1 * 72 / self.fig.dpi]
        if (point_box[0] < -0.8 or point_box[1] < -0.8 or
                point_box[2] > self.w * 72 + 0.8 or point_box[3] > self.h * 72 + 0.8):
            raise ValueError(f"Text leaves page bounds: {text[:100]!r}: {point_box}")
        for previous in self.bounds:
            px0, py0, px1, py1 = previous["rect_pt"]
            intersection_w = min(px1, point_box[2]) - max(px0, point_box[0])
            intersection_h = min(py1, point_box[3]) - max(py0, point_box[1])
            if intersection_w > .8 and intersection_h > .8:
                raise ValueError(f"Overlapping text: {previous['text']!r} and {text!r}")
        self.bounds.append({"text": text, "rect_pt": [round(v, 2) for v in point_box]})
        if url:
            self.page_links.append({"rect": point_box, "url": url})

    def finish_page(self):
        if self.fig is None:
            return
        self.text(str(len(self.pages) + 1).zfill(2), self.w - .95, self.h - .48,
                  size=9, color=MUTED)
        self.text(f"Generated {self.generated.replace('T', ' ').replace('+00:00', ' UTC')}",
                  .7, self.h - .48, size=7.1, color=MUTED)
        self.pdf.savefig(self.fig, dpi=300)
        self.pages.append({"number": len(self.pages) + 1, "kind": self.page_kind,
                           "section_id": self.page_section_id, "title": self.page_title,
                           "size_in": [self.w, self.h], "text_blocks": len(self.bounds),
                           "assets": self.page_assets, "bounds_valid": True,
                           "text_overlap_check": "passed"})
        self.links.append(self.page_links)
        self.fig.clear()
        self.fig = self.canvas = None

    def new_page(self, title, kind="prose", landscape=False, continued=False):
        self.finish_page()
        self.w, self.h = LANDSCAPE if landscape else PAGE
        self.fig = Figure(figsize=(self.w, self.h), dpi=100, facecolor="white")
        self.canvas = FigureCanvasAgg(self.fig)
        self.canvas.draw()
        self.bounds, self.page_links = [], []
        self.page_kind, self.page_title = kind, title
        self.page_section_id, self.page_assets = self.section_id, []
        self.fig.add_artist(Rectangle((.7 / self.w, 1 - .47 / self.h),
                                     .55 / self.w, .045 / self.h,
                                     transform=self.fig.transFigure, color=ACCENT))
        stage = str(self.doc.get("stage", "Manuscript outline")).upper()
        if self.width(stage, 8) > (self.w - 2.1) * 72:
            while stage and self.width(stage + "…", 8) > (self.w - 2.1) * 72:
                stage = stage[:-1]
            stage += "…"
        self.text(stage, 1.4, .36, size=8, color=MUTED)
        self.y = .94
        heading = title + (" · continued" if continued else "")
        for line in self.wrap(heading, 20, weight="bold"):
            self.text(line, .7, self.y, size=20, weight="bold")
            self.y += .355
        self.y += .23
        self.bottom = self.h - .86

    def ensure(self, height):
        if self.y + height > self.bottom:
            self.new_page(self.page_title, kind=self.page_kind, continued=True)

    def paragraph(self, value, size=10.5, color=INK, bullet=False, url=None, gap=.12):
        indent = .18 if bullet else 0
        lines = self.wrap(value, size=size, width=self.w - 1.4 - indent)
        step = size * 1.40 / 72
        # Keep the first two lines together where possible; long paragraphs flow.
        self.ensure(min(2, len(lines)) * step + .07)
        first = True
        for line in lines:
            self.ensure(step)
            if bullet and first:
                self.text("•", .7, self.y, size=size, color=ACCENT)
            self.text(line, .7 + indent, self.y, size=size, color=color, url=url)
            self.y += step
            first = False
        self.y += gap

    def heading(self, title):
        self.ensure(.62)
        self.y += .10
        for line in self.wrap(title, size=12.3, weight="bold"):
            self.text(line, .7, self.y, size=12.3, weight="bold", color=ACCENT)
            self.y += .25
        self.y += .05

    def source(self, source):
        if isinstance(source, str):
            self.paragraph(source, size=9.2, color=MUTED, gap=.035)
            return
        target = source.get("path") or source.get("url")
        label = source.get("label") or source.get("title") or target or "Source"
        if target:
            if not urlparse(str(target)).scheme:
                self.track(self.resolve(target), "source")
            self.paragraph(label, size=9.2, color=ACCENT, url=self.link_for(target), gap=.035)
        else:
            self.paragraph(label, size=9.2, color=MUTED, gap=.035)

    def equation(self, equation):
        if isinstance(equation, str):
            equation = {"math": equation}
        formula = equation["math"].strip().strip("$")
        formula = f"${formula}$"
        size = 14
        while size > 10 and self.width(formula, size) > (self.w - 2.6) * 72:
            size -= .5
        width = self.width(formula, size) / 72
        if width > self.w - 2.6:
            raise ValueError(f"Display equation is too wide: {formula}")
        prop = FontProperties(family="DejaVu Sans", size=size)
        _, height, _ = self.canvas.get_renderer().get_text_width_height_descent(formula, prop, ismath=True)
        height = height / self.fig.dpi + .22
        self.ensure(height + .14)
        self.y += .08
        self.text(formula, (self.w - width) / 2, self.y, size=size)
        if equation.get("number") is not None:
            label = f"({equation['number']})"
            self.text(label, self.w - .7 - self.width(label, 10.5) / 72, self.y,
                      size=10.5)
        self.y += height
        if equation.get("description"):
            self.paragraph(equation["description"], size=10)

    def figure(self, section, asset, index, count):
        path = self.resolve(asset["path"])
        rel = self.track(path, "figure")
        with Image.open(path) as im:
            im.load()
            # Preserve the original pixel grid; transparent images composite white.
            if im.mode in {"RGBA", "LA"} or "transparency" in im.info:
                rgba = im.convert("RGBA")
                image = Image.new("RGB", rgba.size, "white")
                image.paste(rgba, mask=rgba.getchannel("A"))
            else:
                image = im.convert("RGB")
        pw, ph = image.size
        self.files[rel]["pixels"] = [pw, ph]
        title = section["title"]
        if count > 1:
            title += f" · image {index + 1} of {count}"
        self.new_page(title, kind="figure", landscape=pw / ph > 1.24)
        caption = asset.get("caption", "")
        caption_lines = self.wrap(caption, size=9.2) if caption else []
        caption_height = len(caption_lines) * .18 + (.12 if caption else 0)
        if caption_height > 1.2:
            raise ValueError("Asset captions must be short (put detailed text in legend).")
        top = self.y
        available_h = self.bottom - top - caption_height - .08
        available_w = self.w - 1.4
        if available_h <= 0:
            raise ValueError(f"No space remains for figure: {title}")
        scale = min(available_w / pw, available_h / ph)
        iw, ih = pw * scale, ph * scale
        left = .7 + (available_w - iw) / 2
        image_top = top + (available_h - ih) / 2
        axis = self.fig.add_axes([left / self.w, 1 - (image_top + ih) / self.h,
                                 iw / self.w, ih / self.h])
        axis.imshow(image, interpolation="none", aspect="equal")
        axis.set_axis_off()
        axis.set_position([left / self.w, 1 - (image_top + ih) / self.h,
                           iw / self.w, ih / self.h])
        self.y = top + available_h + .10
        for line in caption_lines:
            self.text(line, .7, self.y, size=9.2, color=MUTED,
                      url=self.link_for(asset["path"]))
            self.y += .18
        self.current_asset = rel
        self.page_assets.append(rel)

    def section(self, section):
        self.section_id = section.get("id", section["title"])
        assets = section.get("assets", [])
        for i, asset in enumerate(assets):
            self.figure(section, asset, i, len(assets))
        self.new_page(section["title"], kind="figure_notes")
        if section.get("role"):
            self.paragraph(section["role"], size=10, color=MUTED)
        blocks = [("Detailed legend", "legend", False),
                  ("Narrative points and conclusions", "narrative", True),
                  ("Proposed measurements and experimental follow-up", "measurements", False),
                  ("Scope, limitations and pending validation", "limitations", False)]
        for heading, key, bullet in blocks:
            content = paragraphs(section.get(key))
            if content:
                self.heading(heading)
                for para in content:
                    self.paragraph(para, bullet=bullet)
        if section.get("equations"):
            self.heading("Model definitions and readouts")
            for equation in section["equations"]:
                self.equation(equation)
        sources = section.get("sources", [])
        if sources:
            self.heading("Sources and reproducibility")
            for source in sources:
                self.source(source)

    def build(self, destination):
        with PdfPages(destination, metadata={
                "Title": self.doc["title"],
                "Author": self.doc.get("author", "Scaling sandbox manuscript working collection"),
                "Subject": self.doc.get("subtitle", "Vision-aligned manuscript outline"),
                "Keywords": "manuscript outline; original figures; detailed legends; reproducibility",
                "Creator": "scripts/build_manuscript_outline.py",
                "CreationDate": datetime.now(timezone.utc),
        }) as self.pdf:
            self.new_page("Manuscript outline", kind="cover")
            self.y += .25
            for line in self.wrap(self.doc["title"], size=27, weight="bold"):
                self.text(line, .7, self.y, size=27, weight="bold")
                self.y += .49
            self.y += .15
            if self.doc.get("subtitle"):
                self.paragraph(self.doc["subtitle"], size=13, color=MUTED)
            self.paragraph(f"{self.doc.get('date', self.generated[:10])}  ·  "
                           f"{self.doc.get('stage', 'Working outline')}", size=10, color=ACCENT)
            self.heading("Scientific argument and manuscript plan")
            for para in paragraphs(self.doc.get("overview", self.doc.get("overview_paragraphs"))):
                self.paragraph(para)
            self.heading("Proposed main figures")
            for section in self.doc["sections"]:
                self.paragraph(section["title"], size=10, bullet=True)
            for section in self.doc["sections"]:
                self.section(section)
            supplements = self.doc.get("supplementary", [])
            if supplements:
                self.section_id = None
                for item in supplements:
                    if isinstance(item, dict) and "assets" in item:
                        self.section(item)
                    else:
                        if self.page_kind != "supplementary":
                            self.section_id = None
                            self.new_page("Supplementary material and supporting analyses", kind="supplementary")
                        if isinstance(item, dict):
                            if item.get("title"):
                                self.heading(item["title"])
                            for para in paragraphs(item.get("text", item.get("description", ""))):
                                self.paragraph(para)
                            if item.get("path"):
                                self.source({"label": item.get("source_label", item["path"]),
                                             "path": item["path"]})
                            for source in item.get("sources", []):
                                self.source(source)
                        else:
                            self.paragraph(item, bullet=True)
            references = self.doc.get("references", [])
            if references:
                self.section_id = None
                self.new_page("References and evidence provenance", kind="references")
                for ref in references:
                    self.source(ref)
            self.finish_page()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, help="Manuscript JSON document")
    parser.add_argument("--output", type=Path, help="Output PDF")
    parser.add_argument("--manifest", type=Path, help="Defaults to OUTPUT.build_manifest.json")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--example", action="store_true", help="Print an example JSON document")
    args = parser.parse_args()
    if args.example:
        print(json.dumps(EXAMPLE, indent=2, ensure_ascii=False))
        return
    if not args.input or not args.output:
        parser.error("--input and --output are required unless --example is used")
    repo = args.repo_root.resolve()
    input_path = args.input.resolve()
    output = args.output.resolve()
    manifest_path = (args.manifest or output.with_suffix(".build_manifest.json")).resolve()
    doc = json.loads(input_path.read_text())
    if not doc.get("title") or not doc.get("sections"):
        raise ValueError("The document requires a title and at least one section.")
    for section in doc["sections"]:
        if not section.get("title"):
            raise ValueError("Every section needs a title.")
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    builder = OutlineBuilder(doc, repo, output, input_path)
    # Use temporary output so missing inputs or malformed text never replace an
    # earlier usable manuscript outline.
    with tempfile.TemporaryDirectory(prefix="manuscript_outline_") as tmp:
        intermediate = Path(tmp) / "rendered.pdf"
        builder.build(intermediate)
        reader = PdfReader(intermediate)
        if len(reader.pages) != len(builder.pages):
            raise RuntimeError("Renderer and PDF page counts disagree.")
        writer = PdfWriter()
        writer.append(reader)
        bookmarked_sections = set()
        for i, page in enumerate(builder.pages):
            section_id = page["section_id"]
            if section_id is not None and section_id not in bookmarked_sections:
                writer.add_outline_item(page["title"], i)
                bookmarked_sections.add(section_id)
            elif section_id is None and page["kind"] in {"cover", "references", "supplementary"}:
                # Repeated continuation pages share a title and need one bookmark.
                key = (page["kind"], page["title"])
                if key not in bookmarked_sections:
                    writer.add_outline_item(page["title"], i)
                    bookmarked_sections.add(key)
        annotations = 0
        for i, links in enumerate(builder.links):
            for link in links:
                writer.add_annotation(i, Link(rect=tuple(link["rect"]), url=link["url"], border=[0, 0, 0]))
                annotations += 1
        writer.add_metadata({key: str(value) for key, value in reader.metadata.items()})
        final = Path(tmp) / "final.pdf"
        with final.open("wb") as handle:
            writer.write(handle)
        verify = PdfReader(final)
        text = [page.extract_text() or "" for page in verify.pages]
        if any(not page.strip() for page in text):
            raise RuntimeError("A PDF page lacks extractable text.")
        if not all(math.isfinite(float(v)) for page in verify.pages for v in page.mediabox):
            raise RuntimeError("Invalid PDF page bounds.")
        output.write_bytes(final.read_bytes())
    manifest = {
        "schema_version": 1, "generated_utc": builder.generated,
        "document_stage": doc.get("stage"), "title": doc["title"],
        "pdf": {"path": os.path.relpath(output, repo), "bytes": output.stat().st_size,
                "sha256": digest(output), "pages": len(builder.pages)},
        "inputs": sorted(builder.files.values(), key=lambda item: item["path"]),
        "page_inventory": builder.pages,
        "validation": {"all_local_sources_present": True, "all_text_within_page_bounds": True,
                       "all_text_nonoverlapping": True,
                       "all_pages_have_extractable_text": True, "uncropped_figures": True,
                       "clickable_source_annotations": annotations,
                       "extracted_characters": sum(map(len, text)),
                       "renderer": f"matplotlib {matplotlib.__version__}",
                       "figure_save_dpi": 300},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"pdf": str(output), "manifest": str(manifest_path),
                      "pages": len(builder.pages), "bytes": output.stat().st_size,
                      "annotations": annotations}, indent=2))


if __name__ == "__main__":
    main()
