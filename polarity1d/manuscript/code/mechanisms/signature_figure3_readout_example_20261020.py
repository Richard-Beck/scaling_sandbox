"""Render an actual retained native-field example for the Figure 3 ROI decision.

Read-only diagnostic. This never changes frozen evidence or calibrates A2/B.
"""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import signature_figure3_campaign_20261020 as campaign
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_native_20261020 as native
import signature_figure3_regions_20261020 as regions


def build(run, destination):
    plan = campaign.verify(run)
    slot = 107
    status_path = run/f"workers/primary/{slot:04d}/status.json"
    record = native.read(status_path)
    placement = record["placement"]
    assert placement["model"] == "holmes_model3"
    assert placement["edge_gap_10"] == 1.
    assert record["handover_time"] is None
    assert record["final_winner"]["winner"] == "challenger"
    adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
    C = placement["basal_circumference"]
    old_fraction = placement["initial_10_percent_contour"]["origin"]/C
    cue_fraction = placement["cue_center_fraction"]
    with np.load(run/placement["input"]) as z:
        basal = adapter.primary(z["state"])
    samples = []
    archives = []
    for phase in record["phases"]:
        archive = run/phase["path"]
        assert native.sha(archive) == phase["sha256"]
        archives.append(dict(path=str(archive), sha256=native.sha(archive)))
        with np.load(archive) as z:
            for t, state in zip(z["time"], z["states"]):
                samples.append((phase["absolute_start"]+float(t), state.copy()))
    selected = []
    for target in (0., 480., 485., samples[-1][0]):
        time, state = min(samples, key=lambda item: abs(item[0]-target))
        active = adapter.primary(state)
        old = campaign.scalar(adapter, state, basal, placement)
        proposed = regions.readout(active, C, old_fraction, cue_fraction)
        selected.append(dict(time=time, state=state, active=active, old=old, proposed=proposed))

    # Coordinates start opposite the old peak so that both sites are visible.
    xx = geometry.arc(adapter.q, old_fraction, 1.)*C
    order = np.argsort(xx)
    xx = xx[order]
    old_edge = placement["initial_50_percent_contour"]
    old_mask = (xx >= old_edge["left"]) & (xx <= old_edge["right"])
    cue_center = placement["center_distance"]
    cue_mask = abs(xx-cue_center) <= placement["support"]/2
    colors = dict(incumbent="#20639b", challenger="#e88322")
    before, after = selected[1:3]
    assert before["time"] == 480. and after["time"] == 485.
    assert before["proposed"]["regional_dominance"] < 0 < after["proposed"]["regional_dominance"]
    figure, axes = plt.subplots(4, 2, figsize=(13, 11.5), sharex=True, sharey=True)
    max_y = max(float(v["active"].max()) for v in selected)
    for row, item in enumerate(selected):
        y = item["active"][order]
        old, proposed = item["old"], item["proposed"]
        for col, ax in enumerate(axes[row]):
            ax.plot(xx, y, color="#17212b", lw=1.7)
            ax.set_ylim(0, max_y*1.13)
            ax.set_xlim(-C/2, C/2)
            ax.axvline(0, color=colors["incumbent"], ls=":", lw=.8)
            ax.axvline(cue_center, color=colors["challenger"], ls=":", lw=.8)
            ax.spines[["top", "right"]].set_visible(False)
            ax.text(.015, .93, f"t = {item['time']:g}", transform=ax.transAxes, va="top", fontsize=11)
            if col == 0:
                floor = float(basal.min())
                ax.axvspan(old_edge["left"], old_edge["right"], color=colors["incumbent"], alpha=.07)
                ax.axvspan(cue_center-placement["support"]/2, cue_center+placement["support"]/2,
                           color=colors["challenger"], alpha=.12)
                for mask, color in ((old_mask, colors["incumbent"]), (cue_mask, colors["challenger"])):
                    ax.fill_between(xx, floor, y, where=mask, color=color, alpha=.3)
                ax.axhline(floor, color="#555555", ls="--", lw=.8)
                caption = (f"Incumbent excess = {old['incumbent_excess']:.3f}\n"
                           f"Cue-footprint excess = {old['challenger_excess']:.3f}\n"
                           f"Frozen basal floor = {floor:.3f}")
            else:
                floor = proposed["common_current_floor"]
                ax.axhline(floor, color="#555555", ls="--", lw=.8)
                if proposed["two_region_available"]:
                    for role in ("incumbent", "challenger"):
                        basin = proposed[role+"_basin"]
                        mask = ((adapter.q*C-basin["basin_left_physical"]) % C < basin["basin_width_physical"])[order]
                        ax.fill_between(xx, floor, y, where=mask, color=colors[role], alpha=.3)
                    detail = f"Challenger/incumbent = {proposed['challenger_to_incumbent_ratio']:.4f}"
                else:
                    role = "challenger" if "cue_associated" in proposed["state_class"] else "incumbent"
                    ax.fill_between(xx, floor, y, color=colors[role], alpha=.3)
                    detail = "One resolved region; no two-region ratio"
                caption = (f"Incumbent excess = {proposed.get('incumbent_excess', 0.):.3f}\n"
                           f"Challenger excess = {proposed.get('challenger_excess', 0.):.3f}\n"
                           f"Common current floor = {floor:.3f}\n{detail}")
            ax.text(.985, .94, caption, transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
        axes[row, 0].set_ylabel("Active Rac concentration\n(native model units)")
    axes[0, 0].set_title("Original: fixed incumbent window and cue footprint", fontsize=12)
    axes[0, 1].set_title("Proposal: complete active regions between profile minima", fontsize=12)
    for ax in axes[-1]:
        ax.set_xlabel("Arc distance from initial incumbent peak (native length units)")
    figure.suptitle("Figure 3 measurement decision — identical native profiles, different counting windows\n"
                   "Holmes 3, circumference 37.8, external edge gap 1; native cue held continuously", fontsize=14)
    figure.text(.5, .025, "Blue: incumbent-associated activity. Orange: challenger-associated activity. Dashed line: subtracted background.\n"
                "Original handover is unobserved by 20,000; proposed two-region excess handover lies in saved bracket [480, 485].\n"
                "Amounts integrate active Rac above background (concentration × length); basin boundaries use a continuous subcell saddle fit.",
                ha="center", va="bottom", fontsize=10)
    figure.subplots_adjust(top=.89, bottom=.12, left=.075, right=.98, hspace=.25, wspace=.12)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination.with_suffix(".png"), dpi=180)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)

    # Show whether the ROI result persists with the same floor on both methods.
    profiles = []
    for item in selected:
        proposed = item["proposed"]
        floor_matched = {}
        if proposed["two_region_available"]:
            for role in ("incumbent", "challenger"):
                b = proposed[role+"_basin"]
                floor_matched[role+"_excess_at_frozen_basal_floor"] = regions.excess_integral(
                    item["active"], C, b["basin_left_physical"], b["basin_width_physical"], float(basal.min()))
        profiles.append(dict(time=item["time"], original_fixed_footprint_readout=item["old"],
            proposed_complete_region_readout=proposed, floor_matched_complete_region_diagnostic=floor_matched,
            primary_field=item["active"].tolist()))
    evidence = dict(generated_at=native.now(), purpose="Concrete readout-definition proposal; no A2/B calibration or native rerun",
        source_run=str(run), source_manifest_sha256=native.sha(run/"manifest.json"),
        source_status=str(status_path), source_status_sha256=native.sha(status_path), native_archives=archives,
        basal_input_sha256=native.sha(run/placement["input"]),
        render_source=str(Path(__file__).resolve()), render_source_sha256=native.sha(Path(__file__)),
        region_source=str(Path(regions.__file__).resolve()), region_source_sha256=native.sha(Path(regions.__file__)),
        placement=placement, active_field="Native Holmes 3 active Rac concentration, first native state field",
        amount_units="Native concentration times native arc length; no absolute molecule-count calibration",
        original_handover=record["handover_time"], original_handover_status=record["handover_status"],
        proposed_two_region_handover_saved_bracket=[480., 485.], native_grid_physical_positions=(adapter.q*C).tolist(),
        handover_saved_bracket_evidence=dict(before_time=before["time"], after_time=after["time"],
            before_dominance=before["proposed"]["regional_dominance"],
            after_dominance=after["proposed"]["regional_dominance"],
            resolution="Saved native fields bracket the event; no interpolated event time is asserted"),
        proposed_rule="Two resolved peaks: continuous local-parabolic subcell saddle boundaries define each basin; integrate the actual native piecewise-linear profile above the common instantaneous minimum. Sum partitions all active excess. Handover is challenger excess first exceeding incumbent excess while two peaks remain resolved.",
        ambiguity_rule="One moving/merged region, extra peaks, ambiguous site association or low contrast get explicit separate outcomes; do not invent a continuous two-region handover or eligible B seed.",
        proposed_definition_approved=True,
        approval_authority="User: yes obviously you should not restrict challenger activity counting to the area under the cue.", profiles=profiles,
        outputs={str(destination.with_suffix(s)):native.sha(destination.with_suffix(s)) for s in (".png", ".pdf")})
    native.write(destination.with_suffix(".json"), evidence)
    return dict(png=str(destination.with_suffix(".png")), evidence=str(destination.with_suffix(".json")),
                native_integrations=0, handover_bracket=[480., 485.])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=campaign.DEFAULT)
    parser.add_argument("--destination", type=Path, default=native.ROOT.parents[1]/"docs/tasks/briefs/20261020_figure3_readout_example")
    args = parser.parse_args()
    print(json.dumps(build(args.run.resolve(), args.destination.resolve()), indent=2))


if __name__ == "__main__":
    main()
