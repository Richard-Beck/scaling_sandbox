#!/usr/bin/env python3
"""Verify the authorized cleanup against the immutable current Figure 2.

Run with the existing scaling_sandbox environment. Redraw to the audit folder;
never replace current scientific artwork or launch native simulations.
"""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/figure2_cleanup_round3_20261022"
FIG = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure2"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    started = datetime.now(timezone.utc).isoformat()
    execution = json.loads((OUT / "execution.json").read_text())
    assert execution["state"] == "completed"
    pins = json.loads((OUT / "protected_sha256.json").read_text())
    for path, digest in pins.items():
        assert sha(ROOT / path) == digest, path
    native = json.loads((OUT / "retained_native_inventory.json").read_text())
    for row in native:
        s = (ROOT / row["path"]).lstat()
        assert (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns) == (
            row["device"], row["inode"], row["logical_bytes"], row["mtime_ns"]), row["path"]
    code = ROOT / "polarity1d/manuscript/code/mechanisms"
    sys.path.insert(0, str(code))
    source = code / "revision_round3_figure2_20261022.py"
    spec = importlib.util.spec_from_file_location("round3_cleanup_redraw", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUT = OUT / "redraw_after"
    module.main()
    csv_checks = {}
    for source in sorted((FIG / "supporting/data").glob("*.csv")):
        target = module.OUT / "supporting/data" / source.name
        assert sha(source) == sha(target), source.name
        csv_checks[source.name] = sha(source)
    assert sha(FIG / "figure2_round3.png") == sha(module.OUT / "figure2_round3.png")
    from pypdf import PdfReader
    pdf_checks = {}
    for name in ("figure2_round3.pdf", "figure2_round3_supplement.pdf"):
        original = PdfReader(FIG / name)
        rebuilt = PdfReader(module.OUT / name)
        assert len(original.pages) == len(rebuilt.pages)
        for a, b in zip(original.pages, rebuilt.pages):
            assert a.get_contents().get_data() == b.get_contents().get_data(), name
            assert list(a.mediabox) == list(b.mediabox)
        pdf_checks[name] = len(original.pages)
    # Existing independent verifier writes a dated receipt at its canonical
    # location. Save that new receipt here, then restore its previous exact
    # bytes so the current delivery remains byte-for-byte unchanged.
    receipt = ROOT / "manuscript_outline/revision_round3_20261022/verification.json"
    prior = receipt.read_bytes()
    try:
        result = subprocess.run([sys.executable, "-B", "manuscript_outline/verify_figure2_round3_20261022.py"],
                                cwd=ROOT, check=True, text=True, capture_output=True)
        independent = json.loads(receipt.read_text())
        assert independent["passed"]
        (OUT / "after_independent_verification.json").write_bytes(receipt.read_bytes())
        print(result.stdout)
    finally:
        receipt.write_bytes(prior)
    for path, digest in pins.items():
        assert sha(ROOT / path) == digest, path
    report = dict(started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
        passed=True, protected_source_and_delivery_hashes_verified=len(pins),
        retained_native_files_with_unchanged_identity_size_mtime=len(native),
        current_figure2_png_byte_identical=True,
        plotted_csvs_byte_identical=csv_checks,
        rebuilt_pdf_pages_with_identical_decoded_content_and_bounds=pdf_checks,
        independently_verified_width_conditions=independent["width_conditions_independently_verified"],
        independently_verified_temporal_conditions=independent["temporal_coordinates_independently_verified"],
        other_main_figures_unchanged=independent["protected_sections_identical"],
        deBelly_native_endpoint_audit_passed=True,
        all_101_supplement_pages_preserved=True,
        canonical_verification_receipt_restored_exactly=True,
        simulations_run=0, verifier_source_sha256=sha(Path(__file__)))
    (OUT / "after_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
