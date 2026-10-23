# deBelly: corrected independent endpoint extension

**The old centered extension is flawed and deprecated. Do not use simulations with `center_boundary_motion=True` for scientific inference or current manuscript figures.**

The old two-end extension computed local Rac/tension protrusion commands, then subtracted their average from both commands and endpoint displacements. This forced equal-and-opposite motion even when only one end was active. Because the interior history and relaxation were not transformed consistently, it changed strain, tension, drag and MCA transport. It was not simply removing whole-cell translation.

The corrected extension retains the same reflected local protrusion law at both ends, but imposes each endpoint's own command and relaxed displacement directly (`center_boundary_motion=False`). Reflection symmetry means swapping ends swaps the result; it does not impose simultaneous equal protrusion. An inactive endpoint remains constrained. This is a prescribed-boundary assay, not a stress-free or free-migration cell. Quasistatic mechanics and uniform MCA replenishment retain their declared original-model assumptions.

Use `code/debelly_corrected_launch.py launch --out data/runs/FRESH_NAME` (with an absolute destination or from the manuscript directory) for fresh campaigns; it includes both preparation-search stages before reversal. `code/debelly_corrected_reversal.py` supplies the frozen scientific runner. Its `model_for('debelly')` factory explicitly chooses the corrected boundary rule. Historical core defaults and frozen campaigns are preserved for provenance; direct use of their centered default does not implement the correction.

The fresh campaign is `data/runs/debelly_corrected_three_test_20261008/`. At all 22 original deBelly lengths, independently run homogeneous rest → finite nominal left relative-20% patch → verified basal release, select preparation by the established strength/duration rules, and share that exact checkpoint across fixed width 3, relative width 20%, and finite constant slope 0.24. Both scalar gates, basal recovery, field flags and local nonmonotonic duration refinement retain their definitions. No legacy deBelly checkpoint or trial is reused. Unchanged four-model native archives are reused only with exact input identities and undergo complete independent audit.

Corrected mechanism figures are grouped in `figures/deBelly_mechanism/`. The flow decomposition separates instantaneous prescribed-endpoint, direct Rho-stress and membrane-drag forcing at a fixed state; it is not a complete causal pathway partition. Converted MCA first-moment speeds use nominal coordinates and must not be called laboratory centroid velocities.

Old centered figures and data remain immutable historical records, explicitly deprecated. Current canonical figures may be replaced only after complete corrected execution and independent audit.
