#!/usr/bin/env python3
"""Build literature/index.csv from card frontmatter and check the folder.

Checks:
- every card has the required fields, with allowed values;
- the key matches the filename;
- every card key has an entry in references.bib, and every bib entry has a card;
- every [@key] cited in syntheses/ and positioning.md has a card;
- every theme named in a card has a synthesis file.

Usage: python literature/tools/build_index.py
Exits non-zero if any check fails. Writes index.csv only when all checks pass.
Standard library only.
"""

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CARDS = ROOT / "cards"
SYNTH = ROOT / "syntheses"
BIB = ROOT / "references.bib"
INDEX = ROOT / "index.csv"

REQUIRED = [
    "key", "title", "authors", "year", "venue", "type", "readout",
    "size_varied", "kinetics_fixed", "noise", "role", "themes", "sections",
    "access", "status", "checked",
]
OPTIONAL = ["doi", "url"]
LISTS = {"readout", "role", "themes", "sections"}
ALLOWED = {
    "type": {"theory", "model", "experiment", "review", "method", "data"},
    "readout": {"observed_migration", "ou_inference", "simulated_polarity",
                "observed_polarity", "none"},
    "size_varied": {"none", "persistence", "reorientation", "threshold",
                    "number_of_sites", "copy_number", "shape", "other"},
    "kinetics_fixed": {"yes", "no", "na"},
    "noise": {"intrinsic", "extrinsic", "both", "none", "na"},
    "role": {"background", "precedent", "competitor", "counterpoint",
             "evidence", "method"},
    "sections": {"fig1", "fig2", "fig3", "fig4", "fig5", "fig6",
                 "methods", "discussion"},
    "access": {"full_text", "summary", "abstract", "secondary"},
    "status": {"stub", "read", "verified"},
}
KEY_RE = re.compile(r"^[a-z]+[0-9]{4}[a-z]+$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
CITE_RE = re.compile(r"\[@([a-z0-9]+)\]")


def parse_frontmatter(text):
    """Parse a flat YAML subset: 'key: value' and 'key: [a, b]'."""
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError("no frontmatter")
    meta = {}
    for line in lines[1:]:
        if line.strip() == "---":
            return meta
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"bad frontmatter line: {line!r}")
        name, value = line.split(":", 1)
        name, value = name.strip(), value.strip()
        if value.startswith("[") and value.endswith("]"):
            meta[name] = [v.strip() for v in value[1:-1].split(",") if v.strip()]
        else:
            if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
                value = value[1:-1]
            meta[name] = value
    raise ValueError("frontmatter not closed")


def main():
    errors = []
    rows = []
    themes_available = {p.stem for p in SYNTH.glob("*.md")}

    for path in sorted(CARDS.glob("*.md")):
        try:
            meta = parse_frontmatter(path.read_text(encoding="utf-8"))
        except ValueError as exc:
            errors.append(f"{path.name}: {exc}")
            continue
        for field in REQUIRED:
            if field not in meta or meta[field] in ("", []):
                errors.append(f"{path.name}: missing {field}")
        for field in meta:
            if field not in REQUIRED and field not in OPTIONAL:
                errors.append(f"{path.name}: unknown field {field}")
        key = meta.get("key", "")
        if key != path.stem:
            errors.append(f"{path.name}: key {key!r} does not match filename")
        if key and not KEY_RE.match(key):
            errors.append(f"{path.name}: key {key!r} is not <author><year><word>")
        for field, allowed in ALLOWED.items():
            values = meta.get(field, [])
            values = values if isinstance(values, list) else [values]
            for v in values:
                if v not in allowed:
                    errors.append(f"{path.name}: {field}={v!r} not allowed")
        for field in LISTS:
            if field in meta and not isinstance(meta[field], list):
                errors.append(f"{path.name}: {field} must be a [list]")
        for theme in meta.get("themes", []):
            if theme not in themes_available:
                errors.append(f"{path.name}: no synthesis for theme {theme!r}")
        if not DATE_RE.match(meta.get("checked", "")):
            errors.append(f"{path.name}: checked must be YYYY-MM-DD")
        rows.append(meta)

    card_keys = {r.get("key") for r in rows}

    bib_keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", BIB.read_text(encoding="utf-8")))
    for k in sorted(card_keys - bib_keys):
        errors.append(f"references.bib: no entry for card {k}")
    for k in sorted(bib_keys - card_keys):
        errors.append(f"references.bib: entry {k} has no card")

    cited_files = list(SYNTH.glob("*.md")) + [ROOT / "positioning.md", ROOT / "queue.md"]
    cited_files += list(CARDS.glob("*.md"))
    for path in cited_files:
        if not path.exists():
            continue
        for k in sorted(set(CITE_RE.findall(path.read_text(encoding="utf-8")))):
            if k not in card_keys:
                errors.append(f"{path.relative_to(ROOT)}: cites [@{k}] with no card")

    if errors:
        print(f"{len(errors)} problem(s):", file=sys.stderr)
        for e in errors:
            print("  " + e, file=sys.stderr)
        return 1

    columns = ["key", "year", "title", "authors", "venue", "doi", "type",
               "readout", "size_varied", "kinetics_fixed", "noise", "role",
               "themes", "sections", "access", "status", "checked", "card"]
    with INDEX.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for meta in sorted(rows, key=lambda m: m["key"]):
            out = {c: meta.get(c, "") for c in columns}
            for c in LISTS:
                out[c] = ";".join(meta.get(c, []))
            out["card"] = f"cards/{meta['key']}.md"
            writer.writerow(out)
    print(f"index.csv: {len(rows)} cards; all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
