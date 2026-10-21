"""Read-only per-seed A5 measurement adapter; no native integrations."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path

import numpy as np
import signature_figure4_regions_20261020 as regions


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verified(path, expected):
    path = Path(path).resolve()
    if not expected or sha(path) != expected:
        raise ValueError("Native Figure 4 source identity failed: "+str(path))
    return path


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+".tmp")
    if path.suffix == ".gz":
        # The same native observation may be selected in both common and
        # longest arms. Its byte identity must not depend on wall-clock time
        # or the destination filename.
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8") as stream:
                    json.dump(value, stream, separators=(",", ":"), allow_nan=False)
                    stream.write("\n")
    else:
        temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def analyze_case(case, output, *, inventory_sha256):
    output = Path(output).resolve()
    if (output/"result.json").exists():
        raise FileExistsError("Refusing to overwrite a completed region measurement")
    if any(output == Path(case[key]).parent or output.is_relative_to(Path(case[key]).parent)
           for key in ("profile_path", "scalar_path")):
        raise ValueError("Region measurements cannot modify original native worker directories")
    files = {key:verified(case[key+"_path"], case[key+"_sha256"])
             for key in ("profile", "scalar", "basal_profile", "basal_checkpoint", "source_status")}
    with np.load(files["profile"]) as archive:
        time, primary = archive["time"], archive["primary"]
        C = np.broadcast_to(archive["circumference"], time.shape).copy()
        q = (archive["material_positions"] if "material_positions" in archive else
             archive["x"]/float(np.asarray(archive["circumference"])))
    if (primary.shape != (len(time), case["native_cells"])
            or not np.allclose(q, (np.arange(primary.shape[1])+.5)/primary.shape[1],
                               rtol=0, atol=1e-12)):
        raise ValueError("Native periodic profile coordinates do not match declared cell quadrature")
    if not np.allclose(np.diff(time), case["profile_cadence"], rtol=0, atol=1e-9):
        raise ValueError("Native region profile cadence differs from the case inventory")
    if abs(time[-1]-case["native_horizon"]) > 1e-8:
        raise ValueError("Native region profile does not reach its declared horizon")
    with np.load(files["basal_profile"]) as archive:
        basal = archive["primary"]
        actualC0 = float(archive["circumference"])
    if (basal.ndim != 1 or len(basal) != primary.shape[1]
            or not np.isclose(np.ptp(basal), case["basal_contrast"], rtol=1e-10, atol=1e-12)
            or not np.isclose(actualC0, case["actualC0"], rtol=1e-12, atol=1e-12)):
        raise ValueError("Exact basal profile differs from the canonical region inventory")
    keep = time >= .5*case["native_horizon"]-1e-10
    time, primary, C = time[keep], primary[keep], C[keep]
    spring = case["model"] in ("spring", "spring_ring")
    kernel, main_cadence = (.15, 8.) if spring else (1.5, 80.)
    main_factor = int(round(main_cadence/case["profile_cadence"]))
    if (main_factor < 1 or not np.isclose(main_factor*case["profile_cadence"], main_cadence,
                                        rtol=0, atol=1e-10)):
        raise ValueError("Native profiles cannot provide the predeclared common imaging cadence")
    variants, main = [], None
    for variant in regions.imaging_variants(time, primary, C,
            basal_contrast=case["basal_contrast"], kernel_length=kernel,
            basal_profile=basal, basal_circumference=actualC0):
        name = f"contour_{variant['threshold']:g}_every_{variant['downsample_factor']}_{'optical' if variant['smoothed'] else 'raw'}"
        path = output/(name+".json.gz")
        save(path, variant)
        row = {k:variant[k] for k in ("threshold", "downsample_factor", "smoothed", "observed_basal_contrast")}
        row.update(path=str(path), sha256=sha(path), summary=variant["result"]["summary"])
        variants.append(row)
        if variant["threshold"] == .5 and variant["downsample_factor"] == main_factor and variant["smoothed"]:
            main = row
    if main is None:
        ix = np.unique(np.r_[np.arange(0, len(time), main_factor), len(time)-1])
        observed_basal = float(np.ptp(regions.gaussian_filter1d(basal,
            kernel*len(basal)/actualC0, mode="wrap")))
        result = regions.track_profiles(time[ix], primary[ix], C[ix],
            basal_contrast=observed_basal, threshold=.5, smoothing=kernel)
        path = output/f"primary_common_cadence_{main_cadence:g}.json.gz"
        variant = dict(threshold=.5, downsample_factor=main_factor, smoothed=True,
                       observed_basal_contrast=observed_basal, result=result)
        save(path, variant)
        main = {k:variant[k] for k in ("threshold", "downsample_factor", "smoothed", "observed_basal_contrast")}
        main.update(path=str(path), sha256=sha(path), summary=result["summary"])
    record = dict(generated_utc=datetime.now(timezone.utc).isoformat(), case=case,
        inventory_sha256=inventory_sha256, source_sha256=sha(Path(__file__)),
        tracking_source_sha256=sha(Path(regions.__file__)),
        native_input_sha256={key:sha(path) for key, path in files.items()},
        raw_profile_cadence=case["profile_cadence"], common_imaging_cadence=main_cadence,
        actual_analysis_window=[float(time[0]), float(time[-1])],
        profile_cadence_limits_lifetimes=True, primary=main, sensitivities=variants,
        native_integrations=0, passed=True,
        scientific_scope="A5 per-realization observed-region statistics only; condition uncertainty, forced-state joins and A6 prediction remain separate.")
    save(output/"result.json", record)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--task", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Complete native-profile region analysis requires SLURM")
    inventory_path = verified(args.inventory, args.inventory_sha256)
    cases = json.loads(inventory_path.read_text())["cases"]
    case = cases[args.task]
    result = analyze_case(case, args.output/case["id"], inventory_sha256=args.inventory_sha256)
    print(json.dumps(dict(id=case["id"], passed=result["passed"], primary=result["primary"]["summary"])))


if __name__ == "__main__":
    main()
