"""Existing-data geometry feasibility for shared initial physical cue arcs."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/arc_review/supporting"
NAMES = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes", "spring":"Spring"}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def intersection(d, width, basal, level):
    c = float(basal["circumference"])
    lo, hi = float(basal[f"edge_left_{level}"]), float(basal[f"edge_right_{level}"])
    return sum(max(0., min(hi, d+width/2+j*c)-max(lo, d-width/2+j*c)) for j in (-1,0,1))


def main():
    basal_path = OUT / "data/initial_basal12_edges.csv"
    placement_path = OUT / "data/initial_placement48.csv"
    basals = list(csv.DictReader(basal_path.open()))
    placements = list(csv.DictReader(placement_path.open()))
    assert len(basals) == 12 and len(placements) == 48
    rows, replay = [], []
    for model, name in NAMES.items():
        all_sizes = sorted([b for b in basals if b["model"] == model], key=lambda b:float(b["size_parameter"]))
        widths = {float(p["full_physical_cue_support"]) for p in placements if p["model"] == model}
        assert len(widths) == 1
        width = next(iter(widths))
        for scope, selected in (("all_three", all_sizes), ("upper_pair", all_sizes[1:])):
            for level in ("50", "10"):
                ceiling = min(float(b["circumference"])/2 for b in selected)
                limiting = max(selected,key=lambda b:float(b[f"edge_right_{level}"]))
                lower = float(limiting[f"edge_right_{level}"]) + width/2
                # The forward support must also remain short of the periodic copy
                # of the negative incumbent edge. This may limit a very wide cue.
                wrap_limit = min(float(b["circumference"])+float(b[f"edge_left_{level}"])-width/2 for b in selected)
                upper = min(ceiling,wrap_limit)
                feasible = lower < upper
                interval = f"({lower:.9g}, {upper:.9g}{']' if ceiling < wrap_limit else ')'}" if feasible else "empty"
                row = dict(model=model,model_label=name,size_subset=scope,
                           size_labels="; ".join(b["size_label"] for b in selected),
                           initial_contrast_edge_percent=int(level),full_physical_cue_support=width,
                           common_shortest_arc_ceiling=ceiling,
                           minimum_centre_arc_for_positive_gap=lower,
                           lower_bound_exclusive=True,lower_bound_limiting_basal=limiting["size_label"],
                           periodic_copy_clearance_upper_bound=wrap_limit,
                           available_upper_bound=upper,upper_bound_inclusive=ceiling < wrap_limit,
                           feasible=feasible,available_interval=interval,
                           available_interval_width=max(0.,upper-lower))
                rows.append(row)
                if feasible:
                    # Direct periodic support intersections at an interior point
                    # and the allowed endpoint reproduce the analytic interval.
                    for d in ((lower+upper)/2,upper if ceiling < wrap_limit else upper-1e-10):
                        overlap=max(intersection(d,width,b,level) for b in selected)
                        assert overlap < 1e-10
                        replay.append(dict(model=model,subset=scope,level=level,d=d,max_overlap=overlap))
                else:
                    assert upper <= lower
    path=OUT/"data/shared_arc_feasibility.csv"
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table="\n".join(f'| {r["model_label"]} | {r["size_labels"]} | {r["initial_contrast_edge_percent"]}% | {r["common_shortest_arc_ceiling"]:.7g} | {r["minimum_centre_arc_for_positive_gap"]:.7g} | {r["available_interval"]} |' for r in rows)
    note=f'''# Shared initial-arc feasibility from the existing basals

Generated {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}. This is an initial-geometry calculation only. No new native trajectory, parameter change, size selection or experimental grid is chosen.

The lower bound is the largest positive initial operational edge among the selected basals plus half the native full cue support (3 for chemical models, 0.2 for Spring). A strictly positive gap requires the centre arc to exceed that bound. The common shortest-arc ceiling is half the smallest actual basal circumference. Full support is treated periodically: near an antipode it can wrap across the plotting seam. The periodic copy of the negative edge is also checked. It does not further reduce any interval in this table.

| Model | Selected sizes | Initial edge | Shared d ceiling | d must exceed | Available interval |
|---|---|---:|---:|---:|---|
{table}

Spring's three-size interval is empty at both definitions. With resting C₀ = 2/3 only, the half-contrast interval is (0.821007819, 1.02209036] and the 10%-contrast interval is narrower, (0.962025358, 1.02209036]. Thus narrowing to the upper pair makes an initially separated full footprint geometrically possible, while changing the original size comparison. This does not establish a distinct challenger response or predict its kinetics.

For the chemical three-size range, Wave-pinning has a nonempty interval at both edge definitions; Holmes has a narrow half-contrast interval and an empty 10%-contrast interval. The upper pair enlarges the available range for both. Goryachev permits nonempty intervals in both size subsets and definitions.

Edges are linear-interpolated crossings of min(u) + level × [max(u) − min(u)] in the initial connected region containing the declared preparation axis. The 10% level is a sensitivity definition, not a confidence interval or a biological boundary. Native basal tails need not vanish beyond it. Endpoint contact at the lower bound has zero intersection length but is excluded here because it has no positive physical gap.

Changing a shared centre arc changes the cue's normal angle according to its actual initial circumference. An edge-relative cue instead changes centre distance across sizes. Spring uses exact stretched basal circumference, rather than resting C₀; its material-anchored cue centre has a time-varying physical arc during subsequent mechanics. These calculations concern prescribed initial positions, not evolved site separations, ancestry or cell motion.

Sources: [native basal edges](../data/initial_basal12_edges.csv), [native cue widths and placements](../data/initial_placement48.csv), [full feasibility table](../data/shared_arc_feasibility.csv). Existing native runs, report HTML and frozen review PDFs remain unchanged.
'''
    (OUT/"metadata/SHARED_ARC_FEASIBILITY.md").write_text(note)
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),status="passed",rows=len(rows),basals=12,native_work="none",
               feasible_intervals=sum(r["feasible"] for r in rows),empty_intervals=sum(not r["feasible"] for r in rows),
               periodic_intersection_replays=len(replay),maximum_replayed_overlap=max(r["max_overlap"] for r in replay),
               source_sha256={str(p):sha(p) for p in (basal_path,placement_path,Path(__file__))},
               artifact_sha256={str(p):sha(p) for p in (path,OUT/"metadata/SHARED_ARC_FEASIBILITY.md")})
    (OUT/"metadata/shared_arc_feasibility_validation.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(note);print(json.dumps(audit,indent=2))


if __name__=="__main__":main()
