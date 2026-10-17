# Four-representative presentation

Review PNG: [hosny_four_representatives.png](results/hosny_four_representatives.png).

This is a presentation revision, with no refitting or new mechanism assumptions.
The four choices use exactly the current exhaustive population-specific parameter
fits, paired under identical HIF/RhoA subsets across PACC and non-PACC:

- K=0: no KD-specific parameters.
- K=2: HIF concentration and RhoA concentration (the representative used previously).
- K=2: HIF concentration and RhoA memory.
- K=8: all four parameters may change after either KD.

Both K=2 choices allow the selected parameters to increase or decrease. In
particular, the HIF-concentration/RhoA-memory model is the unrestricted-sign
exhaustive fit, not the earlier decrease-constrained reference. Parameter values
can differ from that reference because the baseline and both KDs were refitted
jointly under the relaxed restrictions. K is per population; across both
populations the total scalar parameter count is eight baseline values plus 2K.

The old mechanism table has been removed, and panels are renumbered:

- A: the same shared Pareto frontier, highlighting the four chosen models.
- B: a 4x4 grid: representative models are rows, parameters are columns.
  Each baseline/HIF/RhoA category has adjacent blue non-PACC and orange PACC bars.
  Every facet has its own linear scale starting at zero, with numerical bar labels.
  Common normoxia values are listed above the grid so they do not dominate the
  separate facet scales. No parameter uncertainty bars are implied. The canvas
  is taller to accommodate this grid while preserving the other panels' physical sizes.
- C: all five reported trajectory summary metrics, displaying all four models'
  quartiles alongside each source observation. Controls remain separate
  observations on shared baselines. Endpoint/perpendicular metrics remain held out.
- D: baseline spaghetti plots for the four choices, each at 1 h and 16 h, with
  common axis limits across choices for each duration. Blue/orange encode
  non-PACC/PACC here; titles identify the model using the global model colours.
  The short trajectories are exact prefixes of the long trajectories, sampled at
  30-minute camera intervals. These plots show baseline populations explicitly.
- E: a 2x2 model facet grid of analytical MSD predictions for all six hypoxic
  baseline/KD states. Population is colour; treatment is line style. Curves overlap
  when predictions coincide. Shading marks times beyond the observed 16 h.
- F: condensed 2x2 baseline ratio facets across the shared frontier, with all four
  selected models highlighted, including the off-frontier HIF-kappa/RhoA-tau choice.
- G: condensed 2x2 within-model ratio profiles, reused from the HIF-kappa/RhoA-kappa
  model only. Other selected models are not silently assigned these profiles.
  All facets use the same error scale; these are not confidence intervals.

All data comparisons are regenerated with 65,536 cells per state using the existing
exact Poisson-event simulator. The frontier scores retain the exhaustive sweep's
32,768-cell independent evaluations. Input hashes, sample size and scope are in
results/metadata.json. Complete bar values, model quartiles and MSD values are in
CSV; PNG is the only image export. No external data were used.

From the repository root, reproduce using the existing environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/hosny/poisson_velocity/shared_mechanisms/presentation/prepare.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/hosny/poisson_velocity/shared_mechanisms/presentation/plot.py

These serial simulation/plotting steps fit within the repository's local resource
limits. The previous compact parameter-comparison figure is retained separately.
