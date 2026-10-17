# Holmes Model 3 regime search

This note records the parameter audit and bounded regime search for the
two-GTPase, no-phosphoinositide version of Holmes et al. Model 3. It is not a
claim that the selected activation rate is a published default.

## Primary-source correction

The original paper's reaction-diffusion equations (1) activate the active
pool by `I_G * cbar(L) * G_inactive / G_total` and transfer the equal opposite
amount out of the inactive pool. We hold `cbar(L)=1`, as the shared adapter
does, because the paper does not report the membrane-cycling constants needed
to calculate its length dependence. Original Model 3 Eqs. (5) and (6), with
PI feedback removed, are:

- Rac activation rate: `Ihat_R / (2 * (1 + (Rho/a1)^n))`.
- Rho activation rate: `Ihat_Rho / (2 * (1 + (Rac/a2)^n))`.

The first-pass regime diagnostic omitted the factor 2 in both denominators. It
also described the Model 4 Table 1 Rac activation value 0.2 as a Model 3
default. Both choices were wrong. `holmes_homogeneous_audit.json` and all
outputs under `spatial/` are retained as a record of that discarded diagnostic
and must not be used to choose a regime or included in the comparison.
Corrected outputs use
`equations_revision=original_PLOS_Eq5_Eq6_factor2_no_PI` and live in
`holmes_eq56_homogeneous_audit.json`, `spatial_eq56/`, or `pulse_eq56/`.

The primary source is Holmes, Lin, Levchenko, and Edelstein-Keshet (2012),
“Modelling Cell Polarization Driven by Synthetic Spatially Graded Rac
Activation,” PLOS Computational Biology 8(6):e1002366,
[doi:10.1371/journal.pcbi.1002366](https://doi.org/10.1371/journal.pcbi.1002366).
Model 3 is described on printed page 6; Eq. (1) and Table 1 appear on printed
pages 3–4. Table 1 does not identify a unique Model 3 basal Rac activation
value. The model class's inherited 0.2 is therefore treated only as a low-end
audit point. The search varies that one activation rate and holds the
reciprocal Hill structure, Rho activation, totals, thresholds, decay rates,
diffusion values, and no-PI structure fixed. No species or feedback loops are
added.

The corrected reaction is in `holmes_regime.py` as `HolmesModel3Eq56`. Its
source-check reports an independently calculated expected rate at Rac=1 and
Rho=1.25: Eq. (5) yields `Ihat_R/4`, and Eq. (6) yields `6.6/4=1.65` for the
current Rho activation setting. The spot check also verifies both reaction
terms sum to zero within each conserved pool. The homogeneous audit stores
nullcline residuals and linear stability for each root.

## Search performed

The cheap homogeneous-root audit tested Rac basal activation values
0.2, 0.5, 1, 2, 3, 4, 6, and 8, with all other parameters fixed. The corrected
audit found one stable homogeneous root at each value. Active Rac at those
roots was approximately 0.032, 0.080, 0.159, 0.320, 0.512, 1.391, 2.125,
and 2.605 μM, respectively. This nullcline scan does not by itself exclude a
finite-amplitude spatial polarization.

The reference-length spatial screen used L=20 μm, 201 cells, no continuing
cue, and a finite conservative local Rac perturbation in a 20%-length patch
near x/L=0.15. The perturbation transferred up to 25% or 75% of the total Rac
concentration from inactive to active in the patch; it did not alter either
spatially integrated Rac or Rho pool. Each trajectory ran for 4,000 seconds.
The control criterion required persistent active-Rac contrast, a stable
nonzero orientation in the final quarter, and a negative spatial correlation
between active Rac and Rho. This tests the spatial front pattern as well as
its persistence after withdrawal.

| Rac activation | 25% perturbation | 75% perturbation |
|---:|:---:|:---:|
| 0.2 | no | no |
| 2.0 | no | no |
| 3.0 | no | yes |
| 3.25 | yes | yes |
| 3.5 | yes | yes |
| 3.75 | yes | yes |
| 4.0 | no | no |

The 3.25–3.75 interval is a successful interior band rather than a boundary
point selected to maximize assay endpoints. For the successful controls,
Rac/Rho spatial correlations were about -0.9999 and maximum total-pool mass
drift was below 1.3e-11 μM·μm (Rac) and 2.0e-12 μM·μm (Rho). The exact
measurements are in `spatial_eq56/`.

A second control starts from the lowest stable homogeneous state and applies a
finite 40-second local Gaussian Rac-activation cue, then checks whether the
pattern persists after cue removal. At Rac activation 3.5, amplitudes 0.2,
0.5, 1, and 2 each produced the same stable front, with final orientation
-0.325, minimum tail Rac contrast 1.119 μM, Rac/Rho correlation -0.99987,
and total-pool drift below 5e-12 μM·μm. The smallest tested local cue (0.2)
therefore suffices to demonstrate stimulus-triggered autonomous polarity.
These outputs are in `pulse_eq56/`.

The common comparison uses the corrected no-PI equations and one activation
value, 3.5, at every length. This is selected from the interior of the
reference-length success band 3.25–3.75; it is not a published default and was
not selected by maximizing four-assay endpoint counts. Assay-specific finite
cues will be calibrated at L=20 μm before the all-length comparison and then
held fixed across lengths. No length-specific tuning is allowed.

## Reference-length assay calibration

The L=20 μm pilots start the direct-reversal assay with a continuous leftward
linear gradient, and start finite-countercue and periodic-following from the
same finite-pulse-trained, no-cue-confirmed left state. The trialed gradient
amplitudes were 0.2, 0.5, 1, and 2. Competition used the fixed two-pulse
protocol (amplitude 0.2, width 0.2L, 40 s, 2% imbalance).

| Test | 0.2 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| Direct reversal latency | censored | 35 s | 22.5 s | 17.5 s |
| Finite countercue shortest successful duration | none through 5120 s | 20 s | 10 s | 10 s |
| Periodic switches, P=100 s | 0/12 | 0/12 | 12/12, 15 s median | 12/12, 10 s median |
| Periodic switches, P=400 s | 0/12 | 12/12, 67.5 s median | 12/12, 30 s median | 12/12, 27.5 s median |

The common gradient amplitude is 1: it is the lowest candidate that gives
stable reversals, finite reversals, and 12/12 periodic switches at both
reference periods. Competition retains its separately specified 0.2 paired
pulse. Finite reversal at 1 is stable for the whole 500-s confirmation tail;
the periodic test succeeds in every scored switch at P=100 and P=400. The full
pilot measurements are saved under `reference_eq56/`.

A representative source-positive control was rerun at 101 and 201 cells. The
final orientations (-0.325437 and -0.325442), minimum tail contrasts (1.119281
and 1.119284 μM), and Rac/Rho correlations (-0.999867 at both resolutions)
agree closely. This confirms the chosen source state does not depend on the
201-cell mesh.

The final fixed-regime, all-length Holmes campaign completed as SLURM job
20631903. It ran 35 tasks over lengths 10, 15, 20, 25, 30, 40, and 50 μm,
with competition, direct reversal, finite countercue, and periodic periods 100
and 400 s. Each task records source hashes before and after the run; all 35
tasks report identical before/after hashes and the same protocol version.
Results and the source-hash manifest are in `holmes_data_final/` and
`holmes_campaign_manifest.json`.

An initial 35-task run (job 20631814) was retained in `data/` for history. A
scoring review then tightened finite confirmation to require every sample in
the final 500-s window to remain beyond the polarity threshold. The source was
changed after that initial run, so its outputs are superseded. The final
campaign reran all 35 tasks with the stricter criterion and per-task hashes.
The initial code/results remain available and are not mixed into the final
report.

| Length (μm) | Competition | Direct reversal (s) | Finite countercue (s) | P=100 latency (s) | P=400 latency (s) |
|---:|---|---:|---:|---:|---:|
| 10 | 3127.5 s | 20 | 5 | 12.5 | 30 |
| 15 | censored | 20 | 10 | 12.5 | 30 |
| 20 | censored | 22.5 | 10 | 15 | 30 |
| 25 | censored | 22.5 | 10 | 15 | 32.5 |
| 30 | censored | 22.5 | 10 | 15 | 32.5 |
| 40 | censored | 25 | 10 | 17.5 | 37.5 |
| 50 | censored | 25 | 10 | 17.5 | 40 |

Competition establishes both peaks at withdrawal for all lengths; only L=10
resolves within the 20,000-s post-cue horizon. All direct preparations meet
the stable left-polarity criterion. Every finite result meets the full
confirmation-window rule. The smallest conservative positive-polarity margin
among finite successes is 0.0374 above the +0.05 threshold at L=50; all
minimum contrasts also exceed 0.05. Direct-preparation orientation remains at
least 0.219 below the -0.05 cutoff. Every periodic task has 12/12 successful
switches, and each reported median latency matches the recorded per-switch
latencies.

## Artifacts

- `holmes_regime.py`: corrected Model 3 kinetics, source spot check,
  homogeneous roots, conservative patch control, and cue-driven control.
- `holmes_eq56_homogeneous_audit.json`: corrected homogeneous-root results.
- `holmes_spatial_search.slurm`: reference-length spatial patch campaign.
- `holmes_pulse_search.slurm`: finite local-cue persistence campaign.
- `holmes_adapter.py`: fixed-regime `run_case(assay,length,period=None,cells=201)` interface.
- `holmes_reference_calibration.py` and `holmes_reference_calibration.slurm`: L=20 μm stimulus calibration.
- `holmes_campaign.slurm` and `holmes_run.py`: final 35-task all-length campaign.
- `reference_eq56/`: calibration outputs; `mesh_eq56/`: representative 101/201-cell refinement.
- `holmes_data_final/` and `holmes_campaign_manifest.json`: final outputs and per-task source hashes.
- `data/`: preserved superseded first run, not used in the figure.
- `spatial/`: superseded pilot with missing Eq. (5)/(6) factors; preserve only
  for audit history, never as valid comparison output.
