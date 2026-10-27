"""Storage and release evidence for the R1/R2 ring figures (no simulations).

The two figure rebuild entry points own extraction and rendering. This command
records independently measured allocation, deliverable hashes and scope checks.
It does not promote a running campaign to complete or delete any records.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "polarity1d/manuscript").is_dir())
OUT = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/ring_r1r2_20261027"
BASE = "802a0ce"


def now():
    return datetime.now(timezone.utc).isoformat()


def command(*args):
    return subprocess.run(args, cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temp.replace(path)


def measure(paths):
    rows = []
    for relative in paths:
        path = ROOT / relative
        if not path.exists():
            rows.append(dict(path=relative, exists=False))
            continue
        allocated = int(command("du", "-sx", "-B1", str(path)).split()[0])
        rows.append(dict(path=relative, exists=True, allocated_bytes=allocated,
                         allocated_gib=allocated / 2**30))
    receipt = dict(measured_utc=now(), measurement="du -sx -B1 (allocated bytes)",
                   atomic=False, groups=rows,
                   filesystem=command("df", "-B1", str(ROOT)),
                   personal_quota="quota -vs produced no quota information; filesystem free space is not personal quota headroom")
    write(OUT / "supporting/metadata/storage_measurement.json", receipt)
    print(json.dumps(receipt, indent=2))


def audit():
    # The base already contains the established figures. No revision of them is
    # part of this release, including the integrated outline or migration text.
    protected = ["manuscript_outline", "current manuscript figs",
                 "wgd_langevin/refined_ou_20261003",
                 "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026",
                 "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023",
                 "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022"]
    changes = command("git", "diff", "--name-only", BASE, "--", *protected).splitlines()
    records, missing, formats = [], [], []
    for label in ("R1", "R2"):
        folder = OUT / f"figure{label}"
        required = [folder / f"figure{label}_main.{ext}" for ext in ("png", "pdf", "svg")]
        required += [folder / "supporting/metadata" / name for name in ("main_legend.md", "scientific_readme.md")]
        required += sorted((folder / "supporting/data").glob("*.csv"))
        required += sorted((folder / "supporting/data").glob("*.csv.gz"))
        if not list((folder / "supporting/data").glob("*.csv")):
            missing.append(str((folder / "supporting/data/*.csv").relative_to(ROOT)))
        for path in required:
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(str(path.relative_to(ROOT)))
            else:
                records.append(dict(path=str(path.relative_to(ROOT)), bytes=path.stat().st_size, sha256=sha(path)))
        svg = folder / f"figure{label}_main.svg"
        pdf = folder / f"figure{label}_main.pdf"
        png = folder / f"figure{label}_main.png"
        width = None
        if svg.is_file():
            value = ET.parse(svg).getroot().attrib["width"]
            match = re.fullmatch(r"([\d.]+)(pt|mm|in|px)", value)
            if match:
                width = float(match[1]) * {"pt": 25.4 / 72, "mm": 1, "in": 25.4, "px": 25.4 / 96}[match[2]]
        pdf_width = None
        if pdf.is_file():
            info = command("pdfinfo", str(pdf))
            match = re.search(r"Page size:\s*([\d.]+) x ([\d.]+) pts", info)
            if match:
                pdf_width = float(match[1]) * 25.4 / 72
        pixels = None
        if png.is_file():
            header = png.open("rb").read(24)
            if header.startswith(b"\x89PNG\r\n\x1a\n"):
                pixels = list(struct.unpack(">II", header[16:24]))
        formats.append(dict(figure=label, svg_width_mm=width, pdf_width_mm=pdf_width,
                            png_pixels=pixels, width_180mm_pass=all(w is not None and abs(w - 180) < .2 for w in (width, pdf_width))))
    wrong_width = [row["figure"] for row in formats if not row["width_180mm_pass"]]
    receipt = dict(generated_utc=now(), base_commit=command("git", "rev-parse", BASE),
                   current_commit=command("git", "rev-parse", "HEAD"),
                   branch=command("git", "branch", "--show-current"),
                   protected_scope_changes=changes, missing=missing, artifacts=records,
                   format_checks=formats, wrong_width=wrong_width,
                   structural_checks_pass=not (changes or missing or wrong_width),
                   qualification="Structural checks and hashes do not prove scientific correctness, campaign completion or visual quality; those require the panel and campaign audits and fresh visual reviews.")
    write(OUT / "supporting/metadata/release_audit.json", receipt)
    print(json.dumps({key: value for key, value in receipt.items() if key != "artifacts"}, indent=2))
    return not (changes or missing or wrong_width)


def portable(python):
    """Run both renderers in a temporary source/data tree with no HPC records."""
    source = ROOT / "polarity1d/manuscript/code/mechanisms"
    rows = []
    with tempfile.TemporaryDirectory(prefix="ring-r1r2-portable-") as temporary:
        tree = Path(temporary)
        code = tree / source.relative_to(ROOT)
        code.mkdir(parents=True)
        for path in source.glob("*.py"):
            shutil.copy2(path, code / path.name)
        destination = tree / OUT.relative_to(ROOT)
        # Inputs to the portable renderer include no native histories, worker
        # outputs or campaign roots. Generated outputs are removed before test.
        shutil.copytree(OUT, destination)
        for label in ("R1", "R2"):
            folder = destination / f"figure{label}"
            for suffix in ("png", "pdf", "svg"):
                (folder / f"figure{label}_main.{suffix}").unlink(missing_ok=True)
            entry = code / f"rebuild_figure{label}_20261027.py"
            args = [python, "-B", str(entry)] + (["render"] if label == "R2" else [])
            result = subprocess.run(args, cwd=tree, capture_output=True, text=True, timeout=600)
            expected = [folder / f"figure{label}_main.{suffix}" for suffix in ("png", "pdf", "svg")]
            rows.append(dict(figure=label, returncode=result.returncode,
                             stdout=result.stdout[-4000:], stderr=result.stderr[-4000:],
                             all_exports_present=all(p.is_file() and p.stat().st_size > 0 for p in expected),
                             png_identical_to_workspace=expected[0].is_file() and sha(expected[0]) == sha(OUT / f"figure{label}" / expected[0].name)))
        no_campaign_records = not (tree / "polarity1d/manuscript/data/runs").exists()
    passed = no_campaign_records and all(row["returncode"] == 0 and row["all_exports_present"] for row in rows)
    receipt = dict(generated_utc=now(), python=python, passed=passed,
                   no_campaign_records_in_test_tree=no_campaign_records, tests=rows,
                   method="Copy executable Python sources and compact deliverables only into a fresh temporary root; remove all figure exports, rebuild R1 and R2, then remove temporary test tree.")
    write(OUT / "supporting/metadata/portable_rebuild_audit.json", receipt)
    print(json.dumps(receipt, indent=2))
    return passed


def scientific_checks():
    """Independent algebra checks on plotted records, not new estimators."""
    def rows(path):
        with path.open(newline="") as stream:
            return list(csv.DictReader(stream))

    first = OUT / "figureR1/supporting"
    second = OUT / "figureR2/supporting"
    checks = []

    def record(name, errors, tolerance=1e-10):
        maximum = max(errors, default=math.inf)
        checks.append(dict(check=name, comparisons=len(errors), maximum_absolute_error=maximum,
                           tolerance=tolerance, passed=bool(errors) and maximum <= tolerance))

    record("R1B closed-form drift null", [abs(float(r["C"]) - math.exp(-(2 * math.pi / float(r["C_label"]))**2 * float(r["Ds"]) * float(r["lag"]))) for r in rows(first / "data/B_analytic_drift.csv")])
    record("R2B closed-form jump null", [abs(float(r["C"]) - math.exp(-float(r["lambda0"]) * (float(r["size"]) / 40)**float(r["beta"]) * float(r["lag"]))) for r in rows(second / "data/B_null_heading_acf.csv")])
    front = json.loads((first / "metadata/isolated_front.json").read_text())
    errors = []
    for row in rows(first / "data/G_conduction_points.csv"):
        circumference, diffusivity, theta = [float(row[key]) for key in ("C", "D", "theta")]
        crossover = 2 * diffusivity / (front["kappa"] * front["delta_u"] * theta * (1 - theta))
        predicted = front["lambda_push"] / 2 / (1 + circumference / crossover)
        errors.extend((abs(crossover - float(row["Cstar"])), abs(predicted - float(row["predicted_gain"]))))
    record("R1G isolated-front prediction with exact-ring theta", errors)
    record("R1H closed-form native exponent prediction", [abs(float(r["alpha"]) - 2 - 2 * float(r["C"]) / (float(r["C"]) + float(r["Cstar"]))) for r in rows(first / "data/H_native_prediction.csv")])
    rates = rows(second / "data/mode1_rates.csv")
    keyed = {(r["arm"], r["model"], r["rate_level"], r["amplitude_level"], float(r["size"])): r for r in rates}
    alpha = rows(second / "data/EF_local_alpha.csv")
    discrepancies, pending, verified = [], 0, 0
    for row in alpha:
        if row["status"] == "execution_pending":
            pending += 1
            continue
        sizes = (1., 2., 3.) if row["model"] == "spring_ring" else (20., 40., 60.)
        selected = [keyed[("per_perimeter" if size == sizes[1] else row["arm"], row["model"], row["rate_level"], row["amplitude_level"], size)] for size in sizes]
        x = [math.log(float(r["mean_physical_C"])) for r in selected]
        mean = sum(x) / 3
        denominator = sum((value - mean)**2 for value in x)
        weights = [-(value - mean) / denominator for value in x]
        lower, upper = 0., 0.
        for weight, rate in zip(weights, selected):
            if abs(weight) < 1e-14:
                continue
            if rate["status"] == "resolved":
                low = high = math.log(float(rate["gamma"]))
            elif rate["status"] == "beyond_analysis_lag":
                low, high = -math.inf, math.log(float(rate["bound"]))
            elif rate["status"] == "below_frame_resolution":
                low, high = math.log(float(rate["bound"])), math.inf
            else:
                raise ValueError("Unrecognized rate status: " + rate["status"])
            lower += weight * (low if weight > 0 else high)
            upper += weight * (high if weight > 0 else low)
        for expected, key in ((lower, "alpha_lower"), (upper, "alpha_upper")):
            observed = float(row[key])
            discrepancies.append(0. if expected == observed else abs(expected - observed))
        verified += 1
    record("R2E/F signed OLS rate-bound propagation using actual circumference", discrepancies)
    passed = all(check["passed"] for check in checks)
    receipt = dict(generated_utc=now(), passed=passed, checks=checks,
                   verified_alpha_conditions=verified, pending_alpha_conditions=pending,
                   full_campaign_completion_proven=pending == 0,
                   scope="Independent algebra checks of plotted analytic curves and OLS interval propagation only. Original rate/heading estimators remain unchanged. Pending execution is not censoring and is not validated by this audit.")
    write(OUT / "supporting/metadata/scientific_algebra_audit.json", receipt)
    print(json.dumps(receipt, indent=2))
    return passed


def main():
    parser = argparse.ArgumentParser(__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    storage = sub.add_parser("measure")
    storage.add_argument("paths", nargs="+")
    sub.add_parser("audit")
    sub.add_parser("science")
    rebuild = sub.add_parser("portable")
    rebuild.add_argument("--python", default="/home/4473331/.conda/envs/scaling_sandbox/bin/python")
    args = parser.parse_args()
    if args.command == "measure":
        measure(args.paths)
    elif args.command == "portable":
        if not portable(args.python):
            raise SystemExit(1)
    elif args.command == "science":
        if not scientific_checks():
            raise SystemExit(1)
    elif not audit():
        raise SystemExit(1)


if __name__ == "__main__":
    main()
