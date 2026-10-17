#!/usr/bin/env python3
"""Compare the body-fitted solver CSV with the authors' saved VTU trajectory.

Run this in the repository's existing ``vtk-render`` Conda environment.
The reference directory is deliberately supplied by the caller: the tracked
implementation does not depend on ``dev/debelly`` or on generated output.
"""

import argparse
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


FIELDS = {
    "rac": ("rac", "rac"),
    "rho": ("rho", "rho"),
    "mcab": ("Ezrin bound", "mcab"),
    "tension": ("tension", "tension"),
}


def probe_vtu(path, arclength, radius):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    theta = arclength / radius
    points = vtk.vtkPoints()
    for angle in theta:
        points.InsertNextPoint(-radius * np.cos(angle), radius * np.sin(angle), 0.0)
    target = vtk.vtkPolyData()
    target.SetPoints(points)
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(target)
    probe.SetSourceData(reader.GetOutput())
    probe.Update()
    data = probe.GetOutput().GetPointData()
    valid = vtk_to_numpy(data.GetArray("vtkValidPointMask")).astype(bool)
    values = {name: vtk_to_numpy(data.GetArray(name)) for name in
              ["rac", "rho", "Ezrin bound", "tension", "uₕ", "x"]}
    return valid, values


def metrics(reference, candidate):
    error = candidate - reference
    span = np.ptp(reference)
    scale = span if span > 1e-12 else max(np.max(np.abs(reference)), 1e-12)
    return np.sqrt(np.mean(error**2)) / scale, np.linalg.norm(error) / np.linalg.norm(reference)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("candidate_csv", type=Path)
    parser.add_argument("--times", nargs="+", type=int, default=[75, 100, 150, 200, 300, 350])
    parser.add_argument("--samples", type=int, default=240)
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--max-nrmse", type=float, default=0.05,
                        help="fail if any profile exceeds this span-normalized RMSE")
    args = parser.parse_args()

    table = np.genfromtxt(args.candidate_csv, delimiter=",", names=True)
    smax = args.radius * np.arccos(-0.9)
    # Avoid exact cut-domain corners, where VTK interpolation is ambiguous.
    sample_s = np.linspace(0.01, smax - 0.01, args.samples)
    print("time field       NRMSE(span)  relative_L2  reference_range       candidate_range")
    worst = (0.0, None, None)
    for time in args.times:
        matches = list(args.reference_dir.rglob(f"_sur_{time}.vtu"))
        if len(matches) != 1:
            raise RuntimeError(f"expected one _sur_{time}.vtu below {args.reference_dir}, found {len(matches)}")
        valid, ref = probe_vtu(matches[0], sample_s, args.radius)
        rows = table[np.isclose(table["time"], time)]
        if not len(rows):
            raise RuntimeError(f"candidate CSV has no time {time}")
        order = np.argsort(rows["arclength"])
        source_s = rows["arclength"][order]
        for label, (vtk_name, csv_name) in FIELDS.items():
            candidate = np.interp(sample_s, source_s, rows[csv_name][order])
            reference = ref[vtk_name]
            use = valid & np.isfinite(reference) & np.isfinite(candidate)
            nrmse, rel = metrics(reference[use], candidate[use])
            worst = max(worst, (nrmse, time, label), key=lambda item: item[0])
            print(f"{time:4d} {label:10s} {nrmse:11.5f} {rel:12.5f}  "
                  f"[{reference[use].min():.6g},{reference[use].max():.6g}]  "
                  f"[{candidate[use].min():.6g},{candidate[use].max():.6g}]")
        for label, vtk_name, cx, cy in [
            ("velocity", "uₕ", "velocity_x", "velocity_y"),
            ("displace", "x", "displacement_x", "displacement_y"),
        ]:
            candidate = np.column_stack([
                np.interp(sample_s, source_s, rows[cx][order]),
                np.interp(sample_s, source_s, rows[cy][order]),
            ])
            reference = ref[vtk_name][:, :2]
            use = valid & np.all(np.isfinite(reference), axis=1) & np.all(np.isfinite(candidate), axis=1)
            refmag = np.linalg.norm(reference[use], axis=1)
            canmag = np.linalg.norm(candidate[use], axis=1)
            nrmse, rel = metrics(refmag, canmag)
            worst = max(worst, (nrmse, time, label), key=lambda item: item[0])
            print(f"{time:4d} {label:10s} {nrmse:11.5f} {rel:12.5f}  "
                  f"[{refmag.min():.6g},{refmag.max():.6g}]  [{canmag.min():.6g},{canmag.max():.6g}]")
    print(f"worst NRMSE(span)={worst[0]:.5f} at time={worst[1]} field={worst[2]}")
    if worst[0] > args.max_nrmse:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
