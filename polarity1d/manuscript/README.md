> **Do not use the flawed centered deBelly extension (2026-10-08).** The old reflection extension with `center_boundary_motion=True` is deprecated: subtracting endpoint averages inside the dynamics creates artificial opposite-end motion and changes tension. Earlier centered-model preparations, reversal thresholds and mechanism figures are historical evidence only, not results for the corrected model. The active correction keeps reflection symmetry with independent local prescribed endpoint commands (`center_boundary_motion=False`). See [the correction and executable entry point](DEBELLY_BOUNDARIES.md). Fresh preparations and three-test reruns cover all 22 deBelly lengths in `data/runs/debelly_corrected_three_test_20261008/`; mechanism figures are regenerated in `figures/deBelly_mechanism/`. The corrected campaign is complete and independently audited; canonical figures now use corrected deBelly results. Previous canonical figures are archived and explicitly deprecated.

# Polarity manuscript workflow

The biological motivation is the evidence for greater migration persistence in
WGD-associated enlarged cells, summarized in the
[working manuscript migration figure](../../current%20manuscript%20figs/main_03_migration_OU.png)
and its [refined OU source analysis](../../wgd_langevin/refined_ou_20261003/README.md).
There is no completed manuscript. Increasing the 1D domain length represents the
size change associated with WGD; it is not an experiment on transient shortening
before a cell turns. The present question is whether size-dependent polarity
mechanisms can help explain the migration phenotype with local kinetic
parameters unchanged. Finite opposing-cue thresholds, incumbent repair and
transport budgets supply mechanistic evidence, while spontaneous directional
persistence and whole-cell migration are distinct outcomes still requiring a
connection. See the [figure collection overview](../../current%20manuscript%20figs/README.md)
for the measured summaries, inferred persistence parameters and source limits.

The current [five-model, three-test comparison](THREE_TESTS.md) uses one exact
nominal relative-20% basal polarized checkpoint at each model/length for fixed
physical width 3 (Spring native 0.2), relative width 20%, and a finite gradient
with constant physical slope. Goryachev's two patches use a single audited
K=64 multiplier (peak 6.4) across every length; their amplitude never scales
with length. Other models retain nominal patch peaks.

The completed replacement has 74,183 audited trials (2,886 new deBelly trials), all 22 original deBelly lengths, and unchanged four-model results. All sampled reversal brackets are resolved. Panel C now shows Spearman correlations using all audited scalar successes, with at least four distinct lengths and nonconstant durations. All 30 cells are available. Native-field stability remains a separate diagnostic, retained in the supplement and supporting counts; it does not filter these correlations. See [the inclusion rule](THREE_TESTS.md#current-main-figure-statistic-2026-10-10).

The canonical [main heatmap](figures/main_polarity.png) and
[length-versus-cue-duration supplement](figures/supplement_width_response.png)
use exactly those three tests. The current [mechanism report](figures/length_mechanisms_20261008/index.html)
and [vector PDF](figures/length_mechanisms_20261008/length_mechanisms.pdf)
connect the length responses to chemical reorganization distance,
Goryachev/Wave-pinning reaction kinetics, incumbent repair, sink separation,
retained domain expansion and corrected deBelly mechanics.

The 14 October author-draft revision has twelve vector figure pages. Figure 1 combines
three cue colors within each model column under the canonical 1% gate, then
sums reflected-distribution Wasserstein distances across total chemical pools.
Active and inactive forms are added before measuring each pool; this spatial
distance increases at every adjacent eligible length in all five models.
Figures 2–7 use grouped facets and direct panel-by-panel prose.

The deBelly sequence uses relative-20% illustrative cues throughout. Figures 8–9
closely follow the author's draft, introducing the model equations alongside
Rac/Rho activation, cortical mechanics and bound-MCA redistribution. The eight
equations use centred publication-style display mathematics, right-aligned
numbers and separate parameter descriptions. Figure 9
puts mechanical illustrations before the local MCA budgets and shows common-time
outward bound-MCA loss across length. Figure 10 presents the full nonmonotone
duration curves, basal endpoint MCA, common-time endpoint changes and incumbent
MCA/Rac trajectories. The intermediate-length interpretation combines a more
depleted starting incumbent with broadly similar first-180-s accumulation over
23–35 µm. At 37–39 µm, accumulation also slows; no common MCA shutdown
concentration or proportional outflux law is inferred. Figure 11 connects basal
profiles and mechanics to that starting depletion, including rightward interior
flow at 19–21 µm. Figure 12 explains persistent incumbent outflow and later
challenged-front MCA rebound at 39 µm.

Figures 1–7 retain their established content and assets. New deBelly figures,
technical notes, primary-paper checks, source hashes and reviews are under
`figures/length_mechanisms_20261008/revision_20261014_draft_alignment/`.
Native trajectories and development originals are preserved. The preceding
reader, PDF and renderer are hash-preserved in
`archive_20261014_before_draft_alignment/`; earlier full-field, recovery,
centroid, equal-command and feedback-ablation companions remain linked from
the archives.

Use the existing scaling_sandbox environment (the vector reader also uses pypdf):

```bash
python code/mechanisms/length_revision_20261014_draft_alignment_figures.py
python code/mechanisms/length_revision_20261014_draft_report.py
```

Commands above run from this manuscript directory. Native replay execution uses
SLURM. Figure/reader rendering reuses retained histories; budget reproduction
performs only the existing half-second diagnostic steps. See the figure
manifests for exact inputs, archived execution sources and audits.
The original [historical reader](figures/three_test_mechanisms/historical_index_20261006.html)
is preserved; its centered deBelly evidence is deprecated.
Current comparison generation time, producer state and audits are in
[data/runs/debelly_corrected_three_test_20261008/](data/runs/debelly_corrected_three_test_20261008/README.md).

The [initial-polarity length report](figures/initial_polarity_length/initial_polarity_length.pdf)
shows all stored pre-cue fields at short, middle and long eligible lengths for
all five models, with pole concentration, half-contrast width, integrated material
and opposite-end inactive concentration across all 97 qualifying checkpoints.
Spring extension and its native global elastic term, and deBelly bound-MCA
polarity, have separate summary panels. deBelly has no explicit inactive Rac
field. Definitions and exclusions are visible in the report. Regenerate with
`code/mechanisms/initial_polarity_report.py`; compact CSVs and checkpoint-hash
validation are in `data/runs/initial_polarity_length_corrected_20261008/`.

The [single manuscript PNG](figures/initial_polarity_length/initial_polarity_length.png)
condenses pages 2–3 into five model columns and four metric rows, emphasizing
increasing pole material through width or concentration. Its
[caption](figures/initial_polarity_length/initial_polarity_length_caption.txt),
single-page PDF and SVG are beside the PNG. Regenerate with
`code/mechanisms/initial_polarity_manuscript_figure.py`; the original report is
preserved, and all 97 input checkpoint hashes are rechecked before rendering. The 12 deBelly starting states now use corrected independent endpoints; the 85 other-model metric rows are exactly unchanged. Historical initial-profile outputs are archived under `figures/archive/initial_polarity_length_centered_deprecated_20261008/`.

The [cue-and-recovery trajectory report](figures/polarity_recovery/polarity_recovery.pdf)
shows active/Rac polarity (plus bound-MCA polarity for deBelly) on logarithmic
elapsed-time axes for all five models, all three cue classes, and three lengths
per model. Each page pairs polarity with old/new pole width, peak, mean and
integrated active material. Its 45 native trajectories preserve censored outcomes
and field-tail flags. Regenerate with `code/mechanisms/polarity_recovery_report.py`;
measurements and provenance are in `data/runs/polarity_recovery_report_20261008/`.
Both new reports live only within the manuscript folder.

The [deBelly boundary-resolution audit](data/runs/debelly_boundary_resolution_20261008/README.md)
identifies in-loop centering as a mechanical boundary-condition change and tests
an independent-endpoint correction at five lengths. All corrected MCA-zero
withdrawals retain reversal, but timing and rate trends differ from the legacy
adapter. The [comparison figure](figures/debelly_boundary_resolution/boundary_resolution.pdf)
and exact source budgets document the correction and its prescribed-endpoint
limits. Earlier deBelly mechanism figures below describe the legacy adapter.

The revised [deBelly advection figure](figures/debelly_advection_v2/debelly_advection.pdf)
retains MCA polarity/rate, physical centroid speed and boundary/Rho/membrane flow
decomposition across nine lengths (23–39 µm). Physical Rho profiles and an
identical-profile domain-extension calculation explain how the same local Rho
forcing generates faster flow in a longer cortex. Reproduction and checks are
in [the revision notes](data/runs/debelly_advection_revision_20261008/README.md).
The [original figure](figures/debelly_advection/debelly_advection.pdf) retains the
complete MCA term budget. Regenerate the revision with
`code/mechanisms/debelly_advection_figure_v2.py`.

The alternative [Goryachev inactive-pool figure](figures/Goryachev_length_mechanism/inactive_pool_perturbation/Goryachev_length_mechanism.png) and [Wave-pinning comparison](figures/Goryachev_wave_pinning_trial/inactive_pool_perturbation/goryachev_wave_pinning_trial.png) remove 0.1%, 0.5% or 1% of the entire inactive profile briefly, then return the exact removed material. Native kinetics remain unchanged. The 0.1-unit primary pulse and 1-unit sensitivity test completed all 288 independently audited histories. Both [reports and methods](data/runs/inactive_pool_figure_exports_20261008/README.md) retain morphology, fractional damage, restoring reaction and the new recovery curves; earlier activation-suppression versions remain separate. Regenerate with `code/mechanisms/inactive_pool_mechanism_figures.py render`.

The [single Goryachev length-mechanism report](figures/Goryachev_length_mechanism/index.html)
centers two supported effects: lower inactive concentration weakens initial
activation by an equal countercue, and the incumbent has a stronger immediate
restoring response per unit damage in longer cells. The
[main figure](figures/Goryachev_length_mechanism/Goryachev_length_mechanism.png)
and [four-page PDF](figures/Goryachev_length_mechanism/Goryachev_length_mechanism_report.pdf)
show full steady active/inactive profiles, the spatial restoring reaction and
a short damage/recovery trace, alongside the suppression length trends.
Supporting pages retain native relative-20% snapshots, late recovery and turnover context. Slower late recovery
remains visible; the two effects are coupled, not exhaustive causal fractions.
Regenerate with `code/mechanisms/goryachev_length_mechanism.py render`.
Compact hash-checked inputs and the
[cleanup inventory](data/runs/goryachev_length_mechanism_20261008/cleanup_inventory.json)
are under `data/runs/goryachev_length_mechanism_20261008/`. Earlier exploratory
Goryachev presentations are archived; the native comparison and audited raw
simulation evidence remain preserved. Fresh restoring-response experiments use
`code/mechanisms/goryachev_restoring_response.py` with a fresh output directory.

The [inactive-pool sink report](figures/inactive_sink_loss/index.html) compares
Goryachev and Wave-pinning incumbent peak-height and active-amount loss rates
across 1,040 positive-sink conditions. Width-1 sinks sit at the opposite end or
3 native length units beyond the initial half-height pole edge, with fixed
positions and exact retained basal starts. No no-sink simulations or reversal
cues were used. The [main figure](figures/inactive_sink_loss/sink_peak_loss.png)
and [PDF](figures/inactive_sink_loss/sink_peak_loss_report.pdf) show strength
saturation, higher loss plateaus with closer sinks, and residual length effects.
Long Wave-pinning cells have a local peak-rate dip despite faster active-amount
loss. Sharp crossings received 121 dense local refinements; fourteen accuracy
replays and full checkpoint/material-budget audits passed. Protocol, measurements
and reproducible code are linked in the
[run notes](data/runs/inactive_sink_loss_20261008/README.md).

The corrected [Wave-pinning diffusion-contrast experiment](data/runs/wave_pinning_diffusion_contrast_20261008/README.md)
sets active diffusion to 0.01 (tenfold decrease) and inactive diffusion to 100
(tenfold increase). It continues the retained relative-20% polarized branch to
verified basal states under the variant, then shares each exact variant checkpoint
across fixed support 3, relative 20%, and finite constant-slope challenges.
Its [status and three length-response plots](data/runs/wave_pinning_diffusion_contrast_20261008/results/index.html)
replace the superseded run that decreased both diffusivities. Final figures
appear after execution and independent audit. The plotted clock is shortest
passing finite cue duration, excluding preparation and basal recovery.

The following reports and the original matched-width comparison are retained
evidence. The current figures supersede their earlier three-width presentation.

The [continuous-gradient PDF report](figures/continuous_fixed_slope/continuous_gradient_report.pdf)
uses the same fixed slopes and exact relative-20% basal checkpoints across all
five models. Its readout is the first recorded time that rightward orientation
meets or exceeds the original leftward magnitude while the gradient stays on.
All 97 eligible integrations and the independent field audit completed:
deBelly and Spring reach the target at every eligible length; Wave-pinning,
Goryachev and Holmes remain below it through 20,000 native units. Histories,
crossing brackets and provenance are in `data/runs/continuous_fixed_slope_20261006/`.

An additional five-model [finite fixed-slope countergradient comparison](GRADIENT_REVERSAL.md)
was added on 2026-10-06. It starts from the same corrected basal-release
checkpoints, retains the existing outcomes and duration search, and uses one
physical slope per model across all lengths and all three preparation conditions.
Its bounded slope pilot and full campaign live under `data/runs/fixed_slope_*_20261006/`;
the added figures are separate from the retained patch comparison.

The [single-PDF mechanism report](figures/fixed_slope_mechanisms/length_duration_mechanisms.pdf)
examines the nonmonotonic duration curves in Wave-pinning, deBelly and Holmes 3.
It uses one preparation, spatial snapshots at three lengths per model, neighboring
passing/failing pulses, extended native recovery and labeled causal diagnostics.
The three preparations have identical duration thresholds at all 55 eligible
model-length pairs examined. Reproduction scripts are in `code/mechanisms/`;
native histories, audits and provenance are in
`data/runs/fixed_slope_mechanisms_20261006/`.

For the minimal Git codeset, start with [CODESET.md](CODESET.md): environment,
file map, verification, figure regeneration and fresh SLURM runs. Git retains
the 291 exact reversal input checkpoints and compact figure inputs; bulk trial
histories and the historical personal explorer remain local.

This is the single working folder for the polarity manuscript. The comparison
uses five models: Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring. Otsuji is
excluded because its initial polarization was not clear; pure LEGI is excluded
because proper basal withdrawal removes its polarity.

The retained 2026-10-03 reversal comparison starts from qualifying **corrected finite-patch
establishment endpoints**, at nominal cue amplitude. Its right challenge mirrors
the left preparation patch, using fixed widths 1.5 and 3 or 20% of instantaneous
length (Spring fixed widths 0.1 and 0.2 in native units). Historical
linear-conditioned reversal results do not supply this figure's length trends.

| Location | Contents |
| --- | --- |
| [code/](code/) | Native solvers, corrected experiment, figure builders, tests and independent audits |
| [figures/index.html](figures/index.html) | Main figure, detailed reversal supplement, corrected establishment heatmaps and personal field explorer |
| [data/establishment_basal/](data/establishment_basal/) | Corrected establishment records and selected native checkpoints/histories |
| [data/runs/](data/runs/) | Fresh five-model matched reversal campaign and compact trial evidence |
| [data/reversal_fields/](data/reversal_fields/) | Retained personal diagnostic only, including its original experiment identity |
| [data/provenance/](data/provenance/) | Source/copy hashes, validation and explicit path mappings |

Read [PROTOCOL.md](PROTOCOL.md) and [ASSUMPTIONS.md](ASSUMPTIONS.md) before
interpreting the results. [REPRODUCE.md](REPRODUCE.md) gives the commands.
[RETENTION.md](RETENTION.md) marks the older seven-model suite and exploratory
campaigns deprecated. No unrelated workstream is included in this cleanup.

The original main explanation and length summary, and the original width supplement,
were generated from the completed fresh campaign after all 73,975 trials passed
independent audit. Across 660 readouts, 375 succeeded, 207 were censored within
the tested range, and 78 were preparation-ineligible. All 291 qualifying no-cue
controls passed. The landing page links final data-generation time and state.
The corrected establishment heatmaps are retained independently. The interactive
old field explorer is for personal inspection; it is outside the manuscript
comparison.
