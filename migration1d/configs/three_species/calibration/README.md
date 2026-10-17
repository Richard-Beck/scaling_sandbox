# Three-species calibration

The first standardized screen follows Tier 1 of the Holmes parameter-mapping
recommendation. It uses the six-field, GTPase-only reduction and varies

1. physical domain size through `Lambda = L*sqrt(d/D_m)`;
2. basal Rac activation through `kappa_R = I_R1/(d_R*R_total)`;
3. inhibitory Hill cooperativity `n`.

All other Holmes kinetics, the three conserved concentrations, and the
diffusivities are fixed at their mapped values. In particular, `D_m=0.1`,
`D_i=50`, and changing `Lambda` changes physical `L`; it does not silently
change transport or turnover.

Each condition is initialized at its relaxed homogeneous rest state. A
Gaussian Rac-activation pulse centered at 10% of cell length acts from 10--30 s,
and the model is observed until 300 s. The primary outcome is now the active-Rac
state defined exactly as in the single-species sweep: normalized amplitude,
connected high-domain count, and late stationarity distinguish homogeneous-low,
homogeneous-high, single-polarized, dual-polarized, multi-peak, coarsening, and
dynamic states.

The analysis also retains a three-species composite annotation. Coherent
polarity requires Cdc42 and Rac range contrasts of at least 0.2, positive
Cdc42--Rac correlation, negative Cdc42--Rho correlation, and oppositely directed
first moments. Linearization of the conserved homogeneous circuit and nonzero
Neumann modes then distinguishes finite-cue polarization (stable rest, persistent
post-cue pattern) from spontaneous linear instability. These are secondary
mechanistic annotations, not substitutes for the common active-Rac outcome.

The production grid contains 13 Lambda values by 13 kappa values at `n=2,3,4`
(507 PDE runs). A matched exact Holmes-baseline versus isometric twofold-volume
comparison is additionally written to `matched_wgd_null.csv`; this applies the
linear factor `2^(1/3)` while holding kinetics and concentrations fixed.

The production WGD screen is configured separately in
`../testing/wgd_null_screen.json`. It evaluates every Tier-1 point at 0.5x, 1x,
and 2x volume as a 1,521-task unthrottled SLURM array with one simulation per
task. Every task writes one disjoint JSON result. A dependent merge job validates
all 1,521 assignments, schemas, keys, outcomes, and conservation errors before
it writes the WGD-specific report, exact transition maps, and transition-count
matrices.

The full experimental extension is configured by
`../testing/wgd_null_extended.json`. It exactly matches the latest
single-species suite's task shape: 6,084 two-seeded-domain runs (three volumes,
three cooperativity slices, four asymmetries, and a 13x13 grid) plus 2,601
polarity-reversal runs (three volumes, three cooperativity slices, and a 17x17
grid). Its 8,685 unthrottled one-simulation tasks write independently, followed
by an `afterok` validation/merge job. Active Rac uses the same domain-count,
eligibility, endpoint-contrast, and sustained-reversal definitions as the
single-species suite; Cdc42/Rac/Rho coherence remains secondary.

The six-field screen enables the solver's local reaction--diffusion Jacobian
sparsity pattern. On the canonical 81-node condition this reproduces the dense
BDF result to machine precision while avoiding dense finite-difference
factorizations that dominate the multi-field runtime.
