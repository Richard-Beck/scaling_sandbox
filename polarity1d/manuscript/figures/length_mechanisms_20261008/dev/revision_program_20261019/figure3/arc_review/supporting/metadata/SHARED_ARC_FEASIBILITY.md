# Shared initial-arc feasibility from the existing basals

Generated 2026-10-19 19:39:50 UTC. This is an initial-geometry calculation only. No new native trajectory, parameter change, size selection or experimental grid is chosen.

The lower bound is the largest positive initial operational edge among the selected basals plus half the native full cue support (3 for chemical models, 0.2 for Spring). A strictly positive gap requires the centre arc to exceed that bound. The common shortest-arc ceiling is half the smallest actual basal circumference. Full support is treated periodically: near an antipode it can wrap across the plotting seam. The periodic copy of the negative edge is also checked. It does not further reduce any interval in this table.

| Model | Selected sizes | Initial edge | Shared d ceiling | d must exceed | Available interval |
|---|---|---:|---:|---:|---|
| Wave-pinning | C=20; C=30; C=40 | 50% | 10 | 7.640032 | (7.64003217, 10] |
| Wave-pinning | C=20; C=30; C=40 | 10% | 10 | 8.80077 | (8.80076969, 10] |
| Wave-pinning | C=30; C=40 | 50% | 15 | 7.640032 | (7.64003217, 15] |
| Wave-pinning | C=30; C=40 | 10% | 15 | 8.80077 | (8.80076969, 15] |
| Goryachev | C=20; C=30; C=40 | 50% | 10 | 2.216895 | (2.2168945, 10] |
| Goryachev | C=20; C=30; C=40 | 10% | 10 | 2.957948 | (2.95794819, 10] |
| Goryachev | C=30; C=40 | 50% | 15 | 2.216895 | (2.2168945, 15] |
| Goryachev | C=30; C=40 | 10% | 15 | 2.957179 | (2.95717901, 15] |
| Holmes | C=20; C=30; C=40 | 50% | 10 | 9.768941 | (9.76894112, 10] |
| Holmes | C=20; C=30; C=40 | 10% | 10 | 10.78902 | empty |
| Holmes | C=30; C=40 | 50% | 15 | 9.768941 | (9.76894112, 15] |
| Holmes | C=30; C=40 | 10% | 15 | 10.78902 | (10.7890228, 15] |
| Spring | C0=1; C0=2; C0=3 | 50% | 0.5207016 | 0.8210078 | empty |
| Spring | C0=1; C0=2; C0=3 | 10% | 0.5207016 | 0.9620254 | empty |
| Spring | C0=2; C0=3 | 50% | 1.02209 | 0.8210078 | (0.821007819, 1.02209036] |
| Spring | C0=2; C0=3 | 10% | 1.02209 | 0.9620254 | (0.962025358, 1.02209036] |

Spring's three-size interval is empty at both definitions. With resting C₀ = 2/3 only, the half-contrast interval is (0.821007819, 1.02209036] and the 10%-contrast interval is narrower, (0.962025358, 1.02209036]. Thus narrowing to the upper pair makes an initially separated full footprint geometrically possible, while changing the original size comparison. This does not establish a distinct challenger response or predict its kinetics.

For the chemical three-size range, Wave-pinning has a nonempty interval at both edge definitions; Holmes has a narrow half-contrast interval and an empty 10%-contrast interval. The upper pair enlarges the available range for both. Goryachev permits nonempty intervals in both size subsets and definitions.

Edges are linear-interpolated crossings of min(u) + level × [max(u) − min(u)] in the initial connected region containing the declared preparation axis. The 10% level is a sensitivity definition, not a confidence interval or a biological boundary. Native basal tails need not vanish beyond it. Endpoint contact at the lower bound has zero intersection length but is excluded here because it has no positive physical gap.

Changing a shared centre arc changes the cue's normal angle according to its actual initial circumference. An edge-relative cue instead changes centre distance across sizes. Spring uses exact stretched basal circumference, rather than resting C₀; its material-anchored cue centre has a time-varying physical arc during subsequent mechanics. These calculations concern prescribed initial positions, not evolved site separations, ancestry or cell motion.

Sources: [native basal edges](../data/initial_basal12_edges.csv), [native cue widths and placements](../data/initial_placement48.csv), [full feasibility table](../data/shared_arc_feasibility.csv). Existing native runs, report HTML and frozen review PDFs remain unchanged.
