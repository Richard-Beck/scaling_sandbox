#!/usr/bin/env python3
"""Re-render the supplied mechanism checkpoint from a frozen Git revision.

Only the declared compact source closure is reconstructed. Large CSVs are
expanded from verified lossless gzip copies in temporary storage. The original
renderer bytes and all scientific functions remain unchanged; a narrow runtime
adapter resolves historical repository paths inside the temporary checkout.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time


BASE = "polarity1d/manuscript"
TOPIC = BASE + "/figures/length_mechanisms_20261008/dev"
CHECKPOINT = TOPIC + "/vision_aligned_checkpoint_20261017"
MANIFEST = CHECKPOINT + "/supporting/metadata/vision_aligned_checkpoint.manifest.json"
COMPRESSION = TOPIC + "/vision_noise_20261017/supporting/metadata/publication_compression_20261017.json"
RENDERER = BASE + "/code/mechanisms/vision_aligned_checkpoint_20261017.py"
LAYOUT = BASE + "/code/mechanisms/length_dev_output_layout.py"
ANALYSIS_AUDIT = TOPIC + "/vision_noise_20261017/supporting/metadata/vision_noise_analysis_audit.json"
DEFAULT_REVISION = "6199cdc"


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_bytes(repo, revision, name):
    return subprocess.check_output(["git", "-C", str(repo), "show", f"{revision}:{name}"])


def safe_name(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or str(path) != name:
        raise ValueError(f"Expected a repository-relative source path: {name}")
    return name


def extract_files(repo, revision, names, destination):
    """Stream only named regular Git blobs; never consult worktree inputs."""
    command = ["git", "-C", str(repo), "archive", "--format=tar", revision, "--", *sorted(names)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    assert process.stdout is not None
    extracted = set()
    try:
        with tarfile.open(fileobj=process.stdout, mode="r|") as archive:
            for member in archive:
                if member.isdir():
                    continue
                name = safe_name(member.name)
                if name not in names or not member.isfile():
                    raise ValueError(f"Unexpected or nonregular Git archive input: {name}")
                target = destination / name
                target.parent.mkdir(parents=True, exist_ok=True)
                incoming = archive.extractfile(member)
                assert incoming is not None
                with incoming, target.open("wb") as outgoing:
                    shutil.copyfileobj(incoming, outgoing, length=2**20)
                extracted.add(name)
    finally:
        process.stdout.close()
        code = process.wait()
    if code or extracted != names:
        raise RuntimeError(f"Git input extraction failed ({code}); missing {sorted(names - extracted)}")


def reconstruct(repo, revision, destination):
    manifest_bytes = git_bytes(repo, revision, MANIFEST)
    manifest = json.loads(manifest_bytes)
    if not manifest["checks"].get("checkpoint_complete"):
        raise ValueError("Selected revision has no completed mechanism checkpoint")
    compression = json.loads(git_bytes(repo, revision, COMPRESSION))
    compressed = {row["raw_path"]: row for row in compression["files"]}
    tree = subprocess.check_output(
        ["git", "-C", str(repo), "ls-tree", "-r", "--name-only", revision], text=True
    ).splitlines()
    tracked = set(tree)
    png_names = {name for name in manifest["sha256"] if name.endswith(".png")}
    names = {MANIFEST, COMPRESSION, LAYOUT} | {CHECKPOINT + "/" + name for name in png_names}
    expanded = []
    for name, expected in manifest["source_sha256"].items():
        safe_name(name)
        if name in tracked:
            names.add(name)
        elif name in compressed:
            row = compressed[name]
            if row["raw_sha256"] != expected or row["decompressed_sha256"] != expected:
                raise ValueError(f"Compression map disagrees with original raw-source identity: {name}")
            names.add(safe_name(row["compressed_path"]))
            expanded.append(row)
        else:
            raise FileNotFoundError(f"Required source is absent from Git and lossless packaging: {name}")
    extract_files(repo, revision, names, destination)
    expansion_records = []
    for row in expanded:
        source = destination / row["compressed_path"]
        if sha(source) != row["compressed_sha256"]:
            raise ValueError(f"Compressed source hash mismatch: {row['compressed_path']}")
        target = destination / row["raw_path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256()
        size = 0
        with gzip.open(source, "rb") as incoming, target.open("wb") as outgoing:
            for block in iter(lambda: incoming.read(2**20), b""):
                digest.update(block)
                size += len(block)
                outgoing.write(block)
        if digest.hexdigest() != row["raw_sha256"] or size != row["raw_bytes"]:
            raise ValueError(f"Decompressed raw-source identity mismatch: {row['raw_path']}")
        expansion_records.append({key: row[key] for key in (
            "raw_path", "raw_sha256", "raw_bytes", "compressed_path",
            "compressed_sha256", "compressed_bytes", "decompressed_sha256")})
    for name, expected in manifest["source_sha256"].items():
        if sha(destination / name) != expected:
            raise ValueError(f"Frozen source hash mismatch: {name}")
    for name in png_names:
        if sha(destination / CHECKPOINT / name) != manifest["sha256"][name]:
            raise ValueError(f"Published PNG identity mismatch: {name}")
    audit = json.loads((destination / ANALYSIS_AUDIT).read_text())
    marker = "/polarity1d/manuscript/"
    generator = audit["generator_path"]
    if marker not in generator or not generator.startswith("/"):
        raise ValueError("Cannot identify the historical repository prefix from the frozen analysis audit")
    historical_root = generator.split(marker, 1)[0]
    return dict(
        manifest=manifest, manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
        extracted_files=len(names), extracted_bytes=sum((destination / name).stat().st_size for name in names),
        expanded_csvs=expansion_records, historical_repository_root=historical_root,
        layout_sha256=sha(destination / LAYOUT),
        original_source_hashes_verified=len(manifest["source_sha256"]),
    )


def worker(config_path):
    """Run immutable render functions with declared-input-only path lookup."""
    config = json.loads(Path(config_path).read_text())
    scratch = Path(config["scratch"]).resolve()
    output = Path(config["output"]).resolve()
    original = Path(config["historical_repository_root"])
    expected = config["manifest"]["source_sha256"]
    blocked = []
    relocated = set()
    read_count = 0
    start = time.perf_counter()
    title_changes = []

    def under(path, parent):
        return path == parent or parent in path.parents

    def guard(event, arguments):
        # Reading source files from the original worktree is forbidden. Python
        # libraries outside it remain available, as does this fresh output.
        if event != "open" or not isinstance(arguments[0], (str, bytes, os.PathLike)):
            return
        path = Path(os.fsdecode(arguments[0])).absolute()
        if under(path, original) and not under(path, output) and not under(path, scratch):
            blocked.append(str(path))
            raise PermissionError(f"Portable render attempted original-worktree access: {path}")

    sys.addaudithook(guard)
    sys.dont_write_bytecode = True
    code = scratch / BASE / "code/mechanisms"
    sys.path.insert(0, str(code))
    spec = importlib.util.spec_from_file_location("frozen_manuscript_mechanisms", scratch / RENDERER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_source = module.Evidence.source

    def portable_source(self, path, expected_digest=None):
        nonlocal read_count
        path = Path(path)
        if path.is_absolute() and under(path, original):
            relative = path.relative_to(original)
            relocated.add(relative.as_posix())
            path = scratch / relative
        elif not path.is_absolute():
            path = scratch / path
        path = path.resolve()
        if not under(path, scratch):
            raise FileNotFoundError(f"Undeclared external renderer input: {path}")
        name = path.relative_to(scratch).as_posix()
        if name not in expected:
            raise FileNotFoundError(f"Renderer requested a source outside the frozen input closure: {name}")
        if expected_digest is not None and expected_digest != expected[name]:
            raise ValueError(f"Nested source identity disagrees with frozen checkpoint: {name}")
        read_count += 1
        return original_source(self, path, expected[name])

    # The sole runtime replacement is file lookup/verification. No figure,
    # selection, transformation, model, assay or statistics function is changed.
    module.Evidence.source = portable_source
    if config["manuscript_numbering"]:
        from matplotlib.figure import Figure
        original_text = Figure.text

        def numbered_text(figure, x, y, text, *args, **kwargs):
            if isinstance(text, str):
                for working, manuscript in ((1, 2), (2, 3), (3, 4)):
                    prefix = f"Figure {working} ·"
                    if text.startswith(prefix):
                        replacement = f"Figure {manuscript} ·" + text[len(prefix):]
                        title_changes.append(dict(original=text, displayed=replacement, x=x, y=y))
                        text = replacement
                        break
            return original_text(figure, x, y, text, *args, **kwargs)

        Figure.text = numbered_text
    module.OUT = output
    module.render_full_checkpoint()
    generated = json.loads((output / "supporting/metadata/vision_aligned_checkpoint.manifest.json").read_text())
    if generated["source_sha256"] != expected:
        raise ValueError("Re-render did not use exactly the original frozen source identities")
    comparison = []
    for name, original_hash in config["manifest"]["sha256"].items():
        path = output / safe_name(name)
        actual = sha(path)
        comparison.append(dict(path=name, original_sha256=original_hash, regenerated_sha256=actual,
                               byte_identical=actual == original_hash, bytes=path.stat().st_size))
    pngs = [row for row in comparison if row["path"].endswith(".png")]
    data = [row for row in comparison if Path(row["path"]).suffix in (".csv", ".npz")]
    pixel_comparison = []
    for row in pngs:
        reference = module.plt.imread(scratch / CHECKPOINT / row["path"])
        actual = module.plt.imread(output / row["path"])
        if reference.shape != actual.shape:
            raise ValueError(f"Regenerated PNG dimensions changed: {row['path']}")
        mask = module.np.any(reference != actual, axis=-1)
        ys, xs = module.np.nonzero(mask)
        title_bottom = int(module.np.ceil(.06 * len(reference)))
        box = [int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1] if len(xs) else None
        body_equal = bool(module.np.array_equal(reference[title_bottom:], actual[title_bottom:]))
        pixel_comparison.append(dict(path=row["path"], shape=list(actual.shape),
                                     changed_pixels=int(mask.sum()), changed_pixel_box=box,
                                     title_band_rows=[0, title_bottom],
                                     pixels_outside_title_band_identical=body_equal,
                                     original_body_array_sha256=hashlib.sha256(reference[title_bottom:].tobytes()).hexdigest(),
                                     regenerated_body_array_sha256=hashlib.sha256(actual[title_bottom:].tobytes()).hexdigest()))
    label_only = all(row["pixels_outside_title_band_identical"] for row in pixel_comparison)
    if config["manuscript_numbering"] and (len(title_changes) != 4 or not label_only):
        raise ValueError("Manuscript numbering changed more than the four title labels")
    for row in data:
        row["equivalent_after_path_normalization"] = row["byte_identical"]
        if row["byte_identical"]:
            continue
        # The unchanged renderer records the actual input checkout in these
        # five metadata cells. Numeric values and serialized arrays must match
        # the published bytes after only this documented prefix normalization.
        if row["path"] != "supporting/data/vision_aligned_linear_longest_examples.csv":
            raise ValueError(f"A scientific data export differs: {row['path']}")
        raw = (output / row["path"]).read_bytes()
        cells = []
        for index, record in enumerate(csv.DictReader(io.StringIO(raw.decode("utf-8")))):
            for key, value in record.items():
                if str(scratch) in value:
                    if key != "source_archive_path":
                        raise ValueError(f"Unexpected portable metadata change: {row['path']}:{key}")
                    cells.append(dict(row=index, column=key, value=value))
        normalized = raw.replace((str(scratch) + "/").encode(), (str(original) + "/").encode())
        digest = hashlib.sha256(normalized).hexdigest()
        if len(cells) != 5 or digest != row["original_sha256"]:
            raise ValueError("Longest-example export changed beyond its five source-path metadata cells")
        row.update(equivalent_after_path_normalization=True,
                   original_prefix_normalized_sha256=digest, relocated_metadata_cells=cells)
    result = dict(
        schema=1, generated_at=now(), revision=config["revision"],
        frozen_checkpoint_manifest_sha256=config["manifest_sha256"],
        frozen_renderer_sha256=sha(scratch / RENDERER),
        portable_wrapper_sha256=sha(__file__),
        input_source_hashes_verified=len(expected), input_source_read_calls=read_count,
        original_source_bytes_unchanged=True, scientific_functions_unchanged=True,
        adapter="Runtime Evidence.source lookup maps the recorded historical repository prefix to the isolated checkout; all raw SHA-256 checks remain exact.",
        historical_paths_relocated=sorted(relocated), original_worktree_access_blocked=True,
        original_worktree_open_attempts=blocked, native_integrations=0,
        temporary_input_files=config["extracted_files"], temporary_input_bytes=config["extracted_bytes"],
        expanded_csvs=config["expanded_csvs"], layout_sha256=config["layout_sha256"],
        runtime_replaced_methods=["Evidence.source"] + (["matplotlib.figure.Figure.text"] if config["manuscript_numbering"] else []),
        original_render_entrypoint="render_full_checkpoint (the original --full branch)",
        source_array_archive_sha256={name:digest for name,digest in expected.items() if name.endswith(".npz")},
        environment=dict(python=platform.python_version(), numpy=module.np.__version__,
                         matplotlib=module.matplotlib.__version__, blas_threads=1),
        output_comparisons=comparison, pngs_byte_identical=all(row["byte_identical"] for row in pngs),
        pngs_compared=len(pngs),
        manuscript_numbering=config["manuscript_numbering"], title_changes=title_changes,
        pixel_comparisons=pixel_comparison, png_pixels_outside_title_band_identical=label_only,
        scientific_data_exports_byte_identical=all(row["byte_identical"] for row in data),
        scientific_data_exports_equal_after_path_normalization=all(row["equivalent_after_path_normalization"] for row in data),
        scientific_data_export_sha256={row["path"]:row["regenerated_sha256"] for row in data},
        render_wall_seconds=time.perf_counter() - start,
        limits=["No native simulations or 243-worker statistical analysis are rerun.",
                "PDF dates and SVG metadata/identifiers may differ after rendering.",
                "Only the source closure declared by the supplied checkpoint is reconstructed."],
    )
    (output / "portable_render_provenance.json").write_text(json.dumps(result, indent=2) + "\n")
    if len(pngs) != 4 or (not config["manuscript_numbering"] and not result["pngs_byte_identical"]):
        raise ValueError("Scientific PNG comparison failed; inspect portable_render_provenance.json")
    print(json.dumps(dict(passed=True, input_sources=len(expected),
                          pngs_byte_identical=sum(row["byte_identical"] for row in pngs),
                          title_only_numbering=config["manuscript_numbering"],
                          original_worktree_open_attempts=len(blocked), output=str(output))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--revision", default=DEFAULT_REVISION, help="Frozen source Git revision; defaults to the first published outline")
    parser.add_argument("--output", type=Path, help="Fresh destination; existing paths are refused")
    parser.add_argument("--scratch-parent", type=Path, help="Optional parent for the automatically removed temporary input checkout")
    parser.add_argument("--manuscript-numbering", action="store_true", help="Change only the four displayed figure titles: working 1/2/3 become manuscript 2/3/4")
    parser.add_argument("--_worker", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args._worker:
        worker(args._worker)
        return
    if args.output is None:
        parser.error("--output is required")
    repo = args.repo.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite an existing output: {output}")
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", f"{args.revision}^{{commit}}"], text=True).strip()
    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="portable_manuscript_mechanisms_", dir=args.scratch_parent) as folder:
        scratch = Path(folder).resolve()
        config = reconstruct(repo, revision, scratch)
        config.update(revision=revision, scratch=str(scratch), output=str(output),
                      manuscript_numbering=args.manuscript_numbering)
        driver = scratch / "scripts/render_manuscript_mechanisms.py"
        driver.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(__file__, driver)
        config_path = scratch / "portable_render_config.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
                           OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONPATH="",
                           MPLCONFIGDIR=str(scratch / ".matplotlib"))
        subprocess.run([sys.executable, str(driver), "--_worker", str(config_path)],
                       cwd=scratch, env=environment, check=True)


if __name__ == "__main__":
    main()
