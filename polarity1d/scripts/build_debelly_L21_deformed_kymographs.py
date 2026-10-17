#!/usr/bin/env python3
"""Display the retained L=21 trajectory in displacement-mapped coordinates."""

from __future__ import annotations

import json
import argparse
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from build_debelly_figure2_kymograph_supplement import (
    ROOT, CAMPAIGN, FIELDS, _rows, _trajectory, _edges, _half_mass_position,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--uncentered', action='store_true', help='Plot the L=21 independent-endpoint run')
    args = parser.parse_args()
    rows = [r for r in _rows(CAMPAIGN / "task_summary.csv")
            if r["assay"] == "single_stimulus_length_test"
            and r["coupling_variant"] == "full" and float(r["length_um"]) == 21]
    if len(rows) != 1:
        raise ValueError("Expected exactly one full-model L=21 trajectory")
    source_campaign = CAMPAIGN
    if args.uncentered:
        source_campaign = ROOT / 'results/campaigns/debelly_L21_uncentered'
        with np.load(source_campaign/'trajectory.npz', allow_pickle=False) as raw:
            payload = {key:raw[key] for key in raw.files}
    else:
        payload = _trajectory(rows[0])
    time = payload["time"] - 200.0
    x = payload["x"]
    reference_edges = _edges(x)
    dx = np.diff(reference_edges)
    displacement = payload["field__displacement"]
    edge_displacement = np.column_stack((
        payload["field__left_boundary_displacement"][:, 0],
        0.5 * (displacement[:, :-1] + displacement[:, 1:]),
        payload["field__right_boundary_displacement"][:, 0],
    ))
    physical_edges = reference_edges[None, :] + edge_displacement
    stretch = np.diff(physical_edges, axis=1) / dx
    if not np.all(np.isfinite(stretch)) or np.min(stretch) <= 0:
        raise ValueError("The displacement map folds or contains nonfinite values")
    local_change = 100 * (stretch - 1)
    length = physical_edges[:, -1] - physical_edges[:, 0]
    np.testing.assert_allclose(
        length, 21 + payload["field__boundary_separation"][:, 0], atol=1e-10,
    )
    np.testing.assert_allclose(np.sum(stretch * dx, axis=1), length)

    saved_fields = list(FIELDS) + [("velocity", "Cortical velocity (µm/s)", "RdBu_r")]
    fields = saved_fields + [("local_length_change", "Local length change (%)", "coolwarm")]
    values = {f: payload[f"field__{f}"] for f, _, _ in FIELDS}
    values["velocity"] = payload["field__velocity"]
    values["local_length_change"] = local_change
    norms = {}
    for field, _, _ in fields:
        a = values[field]
        if field in {"tension", "local_length_change", "velocity"}:
            limit = max(abs(a.min()), abs(a.max()), 1e-12)
            norms[field] = matplotlib.colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        else:
            norms[field] = matplotlib.colors.Normalize(a.min(), a.max())

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "pdf.fonttype": 42})
    fig = plt.figure(figsize=(21, 10), layout="constrained")
    grid = fig.add_gridspec(3, len(fields), height_ratios=[1, 1, 0.48])
    axes = [[fig.add_subplot(grid[r, c]) for c in range(len(fields))] for r in range(2)]
    ylimits = (physical_edges.min() - 0.25, physical_edges.max() + 0.25)
    mappables = {}
    for row, (end, symlog) in enumerate(((time[-1], True), (800.0, False))):
        indices = np.flatnonzero(time <= end)
        # Retain all early samples and thin only the long, slowly varying tail.
        if symlog:
            indices = np.unique(np.r_[indices[time[indices] <= 800], indices[::4], indices[-1]])
        t = time[indices]
        t_edges = _edges(t)
        t_edges[0], t_edges[-1] = t[0], t[-1]
        mapped = np.column_stack([
            np.interp(t_edges, time, physical_edges[:, j])
            for j in range(physical_edges.shape[1])
        ])
        for col, (field, label, cmap) in enumerate(fields):
            ax = axes[row][col]
            mesh = ax.pcolormesh(
                np.broadcast_to(t_edges, mapped.T.shape), mapped.T,
                values[field][indices].T, shading="flat", cmap=cmap,
                norm=norms[field], rasterized=True,
            )
            mappables[field] = mesh
            for boundary in (0, -1):
                ax.plot(time, physical_edges[:, boundary], color="black", lw=0.8)
            for reference in (0, 21):
                ax.axhline(reference, color="gray", ls=":", lw=0.8)
            for event in (0, 160):
                ax.axvline(event, color="white", ls="--", lw=0.8)
            if field in {"rac", "mca_bound"}:
                reference_x50 = 21 * _half_mass_position(payload, field)
                mapped_x50 = np.array([
                    np.interp(pos, reference_edges, edges)
                    for pos, edges in zip(reference_x50, physical_edges)
                ])
                ax.plot(time, mapped_x50, color="white", lw=2)
                ax.plot(time, mapped_x50, color="black", lw=0.7)
            if symlog:
                ax.set_xscale("symlog", linthresh=200)
            ax.set_xlim(time[0], end)
            ax.set_ylim(*ylimits)
            ax.set_xlabel("Time from cue onset (s)" + ("; symlog" if symlog else ""))
            if col == 0:
                ax.set_ylabel(("Full trajectory" if row == 0 else "Early response") + "\nPosition (µm)")
            if row == 0:
                ax.set_title(label)
    for col, (field, label, _) in enumerate(fields):
        cb = fig.colorbar(mappables[field], ax=[axes[0][col], axes[1][col]],
                          orientation="horizontal", fraction=0.045, pad=0.025)
        cb.set_label(label)

    length_ax = fig.add_subplot(grid[2, :3] if args.uncentered else grid[2, :])
    length_ax.plot(time, length, color="#245e89", lw=1.6)
    length_ax.axhline(21, color="gray", ls=":", lw=1)
    length_ax.axvspan(0, 160, color="#f2b134", alpha=0.25, label="Rac cue (160 s)")
    length_ax.set_xscale("symlog", linthresh=200)
    length_ax.set_xlim(time[0], time[-1])
    length_ax.set_xlabel("Time from cue onset (s; symlog)")
    length_ax.set_ylabel("Total length (µm)")
    length_ax.legend(loc="upper right", frameon=False)
    if args.uncentered:
        endpoint_ax = fig.add_subplot(grid[2, 3:])
        for field, label in (("left_boundary_displacement", "Left (stimulated)"),
                             ("right_boundary_displacement", "Right")):
            endpoint_ax.plot(time, payload[f'field__{field}'][:,0], label=label)
        endpoint_ax.axvspan(0,160,color='#f2b134',alpha=.25,lw=0)
        endpoint_ax.set_xscale('symlog',linthresh=200)
        endpoint_ax.set_xlim(time[0],time[-1])
        endpoint_ax.set_xlabel('Time from cue onset (s; symlog)')
        endpoint_ax.set_ylabel('Endpoint displacement (µm)')
        endpoint_ax.legend(frameon=False)
    title = ('L₀ = 21 µm · symmetric endpoint rules · NO CENTERING · left Rac cue · no substrate clutch\n'
             if args.uncentered else 'L₀ = 21 µm · full model · deformed physical coordinates\n')
    fig.suptitle(title +
                 "Color shows the saved fields or local extension/compression; black boundaries show cell extent",
                 fontsize=14)

    stem = "Figure2_supplement_L21_deformed_kymographs"
    if args.uncentered:
        stem += '_uncentered'
    output = ROOT / "results/figures"
    output.mkdir(exist_ok=True, parents=True)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"{stem}.{extension}", dpi=250)
    plt.close(fig)
    data = ROOT / "results/data"
    data.mkdir(exist_ok=True, parents=True)
    np.savez_compressed(data / f"{stem}.npz", time_from_cue_s=time,
                        reference_edges_um=reference_edges, physical_edges_um=physical_edges,
                        local_length_change_percent=local_change, total_length_um=length,
                        **{f: values[f] for f, _, _ in saved_fields})
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_task": rows[0]["task_id"], "source_campaign": str(source_campaign),
        "center_boundary_motion": not args.uncentered,
        "length_range_um": [float(length.min()), float(length.max())],
        "local_length_change_range_percent": [float(local_change.min()), float(local_change.max())],
        "coordinate_map": "Physical position = reference position + saved membrane displacement. Interior edge displacement is the mean of adjacent cell centers; boundary displacements are exact saved endpoints.",
        "local_length_change": "100 * (deformed bin width / reference bin width - 1). Positive means extension; negative means compression.",
        "cortical_velocity": "Saved field__velocity in the model reference frame (micrometres/second). Blue is negative (leftward); red is positive (rightward). Symmetric color limits centered on zero are shared across both time views. No substrate-friction or translation post-processing is applied.",
        "interpretation": "Visualization of deformation on the fixed reference-domain model, not a new moving-domain simulation. Field amplitudes are unchanged; no concentration Jacobian correction is applied. Rac/MCA lines map the original reference-coordinate half-mass positions into physical space.",
        "display": "White dashed lines delimit the cue (0–160 s). Gray dotted horizontal lines mark reference endpoints. Full duration above, linear early-time view below; identical colors and spatial limits in both rows.",
    }
    (data / f"{stem}.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    print(output / f"{stem}.png")


if __name__ == "__main__":
    main()
