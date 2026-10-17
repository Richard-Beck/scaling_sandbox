# Descriptions and methods for the four main figures

This document describes the four retained figures below at the level of their lettered panels and internal facets. Descriptions were checked against the PNGs, plotting scripts, campaign manifests, and retained result tables. They describe the current versions of these figures; older figures and superseded protocols elsewhere in the repository are not substituted for them. Equations and parameter relationships are written in plain text.

1. [Polarity Figure 1: cross-model comparison](polarity1d/results/figures/Figure1_WP1_cross_model.png)
2. [Polarity Figure 2: De Belly mechanochemical model](polarity1d/results/figures/Figure2_standalone_v2_julia_activation_defaults_linear50.png)
3. [Gaussian velocity model for literature migration data](wgd_langevin/gaussian_literature/results/literature_gaussian_overview.png)
4. [Hosny migration: four representative knockdown models](wgd_langevin/hosny/poisson_velocity/shared_mechanisms/presentation/results/hosny_four_representatives.png)

## 1. Polarity Figure 1: competition, reversal, and dynamic tracking across four models

**Purpose.** This figure compares three distinct demands on intracellular polarity: selecting one front after stimulation at both ends, replacing an established front when a gradient reverses, and following a continuously changing directional signal. The comparison separates basic responsiveness to stimulation from the ability to coordinate polarity across a cell of increasing length.

### Models and common methods

Four one-dimensional models are compared: mass-conserved wave-pinning, the Otsuji redistribution model, the reduced Goryachev Cdc42 feedback model, and local excitation/global inhibition (LEGI). The plotted observables are active GTPase for wave-pinning and Otsuji, active Cdc42 for Goryachev, and the response species for LEGI. LEGI is a signal-sensing comparator; its behavior is not interpreted as autonomous competition between self-maintaining polarity fronts.

The reaction–diffusion systems are solved with conservative, cell-centred finite volumes and an adaptive sparse BDF integrator. The competition and reversal length scans cover 1–39 µm in increments of 2 µm. The retained competition runs end at 20,000 s; reversal runs have 20,000 s of follow-up after reversal begins. These common observation horizons are important because an absent event can reflect a delay beyond observation rather than an impossible transition.

A spatial profile is classed as polarized when its contrast, defined as the maximum minus the minimum divided by the sum of their absolute values, reaches 0.1. Front detection also imposes peak-prominence and separation criteria. Competition resolution requires an established multifront state followed by a sustained single-front state; the competition persistence requirement is 50 s. The curves summarize the median resolution latency among runs with observed resolution, rather than assigning a numerical event time to unresolved runs. Competition uses five small, mass-conserving initial perturbations per length for the three polarity models; LEGI is deterministic.

For directional measurements, orientation is the normalized difference between the mean response in the rightmost and leftmost 20% of the domain. Reversal in panel F uses the stricter front-based reversal output: an initially polarized, oriented state must be replaced by a qualifying opposite-side polarized state, with no more than one detected maximum and at least 20 s of persistence. The orientation cutoff is 0.05. Runs that were not appropriately polarized before the switch are distinguished in the data from conditioned runs that did not reverse.

The profile panels use relative position along the cell. Within each model facet, the displayed snapshots share a common minimum-to-maximum rescaling. Their heights therefore compare the evolution of that model's response, not absolute concentrations across different models.

### Panel A — Input for two-front competition

The grayscale image shows the imposed stimulus as a function of time and relative position, with localized stimulation at both cell ends. Each cue occupies the outer 10% of the domain. The nominal peak activation rate is 0.2 per second, and a small prescribed left–right asymmetry breaks exact symmetry. The cues remain at full strength through 20 s and taper to zero by 40 s. The red dashed line marks 20 s, the first profile time in panel B.

Because cue width is a fraction of cell length, increasing length also increases the physical width of the stimulated region. The length scan therefore describes this specific stimulus geometry.

### Panel B — Representative spatial responses during competition

Four facets show wave-pinning, Otsuji, Goryachev, and LEGI at a cell length of 15 µm. Curves are snapshots at 20, 200, and 20,000 s. The polarity models initially respond at both ends, but one side can subsequently win: wave-pinning retains a relatively broad high-activity region, whereas Otsuji and Goryachev develop more sharply localized profiles. Their rates of eliminating the competing response differ markedly. LEGI provides the corresponding sensing response and is not included in the intrinsic competition-time comparison.

### Panel C — Competition-resolution time versus cell length

The three colored curves show resolution latency for wave-pinning, Otsuji, and Goryachev on a logarithmic time axis. The vertical dotted line marks 15 µm, linking the scan to panel B. Horizontal red guides mark 200 and 20,000 s. The translucent regions above each curve indicate the length–time region after the observed resolution event, up to the observation horizon; they are not uncertainty bands.

Wave-pinning and Otsuji show strong increases in resolution time with length and eventually have no plotted resolution within the observation window. Goryachev resolves over a wider length range, with a shallower increase at larger sizes. Missing points are not plotted as explicit censoring symbols and should not be read as zero latency or proof that competition can never resolve.

### Panel D — Input for gradient reversal

The grayscale image shows a persistent linear gradient initially favoring the left end and subsequently favoring the right end. Time zero is the start of a 20-s cross-fade between the original and reflected gradients. Red dashed lines mark the profile times immediately before reversal, 20 s after reversal starts, and 100 s after it starts.

The actual simulations condition the cell for at least 2,000 s before reversal, with longer conditioning at larger lengths. The plot shows only the nearby part of that history.

### Panel E — Representative responses to reversal

The four model facets again use a 15-µm cell. Curves show the response immediately before reversal and at 20 and 100 s afterward. Wave-pinning visibly relocates its broad response to the newly favored side. LEGI also follows the new gradient. Otsuji remains concentrated at the original side over these times, while Goryachev retains a strong old-side response and reverses more slowly. These profiles illustrate why responding to an initial cue does not guarantee rapid replacement of an established front.

### Panel F — Reversal time versus cell length

Successful reversal times are plotted against length on a logarithmic time axis. Wave-pinning and LEGI reverse relatively quickly across much of the scan. Goryachev reversal slows sharply with increasing length and is no longer observed at the larger lengths within the common horizon. Otsuji contributes no successful reversal curve under this protocol.

The 15-µm guide and 20- and 100-s horizontal lines connect the scan to panel E. As in panel C, shaded regions indicate times after successful events, not uncertainty. Missing results include both insufficient initial conditioning and failure to reverse after conditioning; these categories remain distinct in the source table.

### Panel G — Input for dynamic tracking

The grayscale image shows half a cycle of a nonnegative linear gradient whose direction changes sinusoidally. The illustrated period is 200 s. Red dashed lines identify 25, 50, and 75 s within the displayed cycle, corresponding to the snapshots in panel H. Thus this assay tests continuing adaptation to a moving directional preference, rather than response to a single abrupt switch.

### Panel H — Spatial profiles during a tracking cycle

The four facets show a 15-µm cell driven with a 200-s period, sampled at the three indicated times in the final simulated cycle. Wave-pinning and LEGI substantially redistribute their response as the cue changes. Otsuji remains strongly localized at the original side. Goryachev changes its spatial profile but displays weaker directional modulation under this condition. The snapshots provide spatial context for the amplitude and phase summaries in panel I.

### Panel I — Tracking amplitude and phase across lengths and periods

Each model has two heatmaps: response amplitude on the left and phase lag, expressed as a fraction of the forcing period, on the right. Rows are lengths of 5, 15, 25, and 39 µm; columns are periods of 50, 100, 200, and 400 s. Red boxes identify the 15-µm, 200-s examples in panel H.

Each run spans six cycles, with the first cycle excluded from scoring. Orientation is analyzed at the forcing frequency to quantify amplitude and phase. Colored cells must pass both the original response criteria and a waveform audit: amplitude at least 0.05, absolute lag no greater than one-quarter period, positive cue–response correlation, first-harmonic R-squared at least 0.8, cycle-amplitude coefficient of variation no greater than 0.25, and circular phase dispersion no greater than 0.1. Gray cells with crosses fail the audit and are excluded from the color scales; they do not necessarily have literally zero response.

Wave-pinning generally tracks more strongly when the forcing is slow relative to the cell's response, especially at smaller lengths. Otsuji fails the displayed tracking criteria throughout. Goryachev has relatively weak tracking amplitudes, and LEGI loses reliable tracking for some combinations of large size and short period.

**Interpretive scope.** These are protocol-specific comparisons of qualitative equation-level reconstructions. The numerical transition lengths should not be presented as measured biological size thresholds or exact reproductions of the original publications.

**Local provenance:** [figure builder](polarity1d/scripts/build_wp1_manuscript_figures.py), [common-horizon campaign](polarity1d/campaigns/wp1_figure1_common_horizon_v1.json), [tracking campaign](polarity1d/campaigns/retained_phase1_tracking_v1.json), and [panel source data](polarity1d/results/data/manuscript/).

## 2. Polarity Figure 2: establishment and reversal in the De Belly mechanochemical model

**Purpose.** This figure examines how cell length and mechanochemical feedback affect the establishment, maintenance, and replacement of intracellular polarity. It follows both Rac and bound membrane–cortex attachment components (MCA), so a transient Rac response can be distinguished from a lasting redistribution of the coupled system.

### Model, geometry, and readouts

The model combines local Rac–Rho antagonism with mechanics. Rac promotes protrusion; mechanical tension feeds back onto Rho; Rho-generated active stress drives cortical flow; and redistribution of MCA feeds back onto Rac. Bound and unbound MCA exchange and are transported by the mechanical velocities, while Rac and Rho undergo reaction and diffusion.

The current figure uses the explicitly defined **reflection-symmetric extension** of the one-dimensional model. The same protrusion law acts at both ends with opposite outward directions, and mean displacement and velocity are removed to suppress rigid translation. Consequently, these experiments measure intracellular redistribution on a deforming reference cell, not the turning or migration of a whole cell. This distinction also separates this figure from older one-sided, anchored-boundary assays in the repository.

The activation settings are alpha0 = beta0 = 1, alpha = 0.2, and baseline beta = 0.2. Alpha controls MCA-dependent Rac feedback; beta controls tension-dependent Rho feedback. The beta sweep uses 0, 0.2, and 0.4, labeled 0x, 1x, and 2x. Simulations start from a homogeneous algebraic fixed point and use the finite-volume mechanochemical solver with a 0.5-s timestep, conservative transport, and positive MCA binding exchange. The length scan covers 1–39 µm in 2-µm increments.

The vertical coordinate in the trajectory panels is x50/L: the interpolated location containing half of the total field mass to its left, divided by cell length. A value of 0.5 is centered; lower and higher values indicate leftward and rightward mass biases. This is a mass-position readout, not a peak location or the edge-orientation score used in Figure 1. Color encodes cell length, from purple for short cells to yellow for long cells. Time axes use a linear region near zero and logarithmic scaling farther away.

A field is considered settled when its x50/L enters and subsequently remains within 0.002 of the median over the final 1,000 s. At least 1,000 s of confirmation must remain. A settled state is called polarized only if its displacement from 0.5 is at least 0.02. These are operational criteria for the plotted polarity coordinate, rather than a proof that every model variable has reached an exact fixed point.

### Panel A — Response to one 160-s stimulus under four model variants

The eight plots form four rows and two columns. Columns show Rac and bound MCA. Rows show the full model, beta set to zero, alpha set to zero, and a local-inhibition control. The local-inhibition control removes both mechanical feedbacks, protrusion, and active stress while retaining the local biochemical interactions.

A cue stimulates the leftmost 20% of the cell for 160 s, including an 80-s taper; the pale shaded interval marks the cue. The actual simulation runs to 10,000 s, and the displayed time is relative to stimulus onset at simulation time 200 s. Every colored trace corresponds to a different length. Markers identify the first qualifying polarized steady state.

All variants can show a transient leftward Rac response during stimulation. In the full model and beta-zero variant, sufficiently large cells retain a displaced Rac distribution and an oppositely biased bound-MCA distribution. Removing alpha eliminates lasting polarity in these readouts: Rac returns toward the center and the MCA redistribution decays. The local-inhibition control similarly relaxes toward the centered state, with bound MCA remaining centered. Thus immediate cue responsiveness and self-maintained polarity are distinct outcomes.

### Panel B — Time to polarized steady state after the single cue

This plot summarizes only the full-model results from panel A that reach a qualifying polarized steady state. Green and purple curves give the Rac and bound-MCA settling times measured from stimulus onset.

Both curves are nonmonotonic: the smallest lengths that retain polarity settle relatively slowly, intermediate lengths settle more quickly, and larger cells take progressively longer. Bound MCA generally settles later than Rac. Lengths without a qualifying polarized steady state are omitted. “Polarization time” here means sustained settling of x50/L, not the first detectable response during the cue.

### Panel C — Competing-cue responses with different beta values

The six plots show Rac and bound MCA in columns, with beta = 0x, 1x, and 2x in rows. Unlike panel A, these runs receive **two competing end cues** over the outer 20% at each end, with prescribed asymmetry, from 0 to 80 s and a 40-s taper. They are followed to 20,000 s. The shaded region marks this shorter competition stimulus.

After the immediate response, the trajectories can spend a long period near the centered state before developing a lasting Rac/MCA asymmetry. Shorter cells remain near the center, while larger cells develop displaced steady states. Circles mark polarized steady states and crosses mark settled but unpolarized states. Changing beta alters both the final displacement and the delay before settling, with particularly visible effects on bound MCA.

### Panel D — Settling time versus length after competing cues

Separate plots summarize Rac and bound MCA from panel C. Orange, purple, and green denote beta = 0x, 1x, and 2x. Only polarized, uncensored settling events are shown.

The curves again have a broad minimum at intermediate lengths. The beta-zero cases tend to settle sooner, and stronger beta feedback generally increases the settling time over the larger-length range. The bound-MCA times are typically longer than the Rac times. These values should not be directly treated as repetitions of panel B, because the single-cue and competing-cue protocols differ.

### Panel E — Trajectories during and after finite opposing cues

The six plots use the full model at beta = 1x. Columns show Rac and bound MCA; rows show opposing-cue durations of 150, 300, and 600 s. A 600-s left-end conditioning pulse first establishes a state, and the mirrored right-end counter-cue is applied only after both Rac and bound MCA satisfy the settling criterion. Time zero is the start of that counter-cue, and orange shading marks its duration. Each simulation continues for 10,000 s after counter-cue removal.

Rac responds strongly during the opposing stimulus, but this immediate deflection does not always persist afterward. Bound MCA often changes more slowly. Short counter-cues commonly lead to recovery toward the previous polarity; longer counter-cues can produce sustained reversal for some lengths. At other lengths, even a strong temporary response relaxes back or fails the final reversal criterion. The plots make the difference between forced displacement and lasting repolarization visible.

### Panel F — Shortest tested counter-cue producing sustained reversal

The two plots show the shortest successful duration based on Rac and bound MCA, respectively, for each length and beta setting. Durations were tested on a grid of 50 equally spaced values between 10 and 1,000 s. These are therefore discrete tested thresholds, with a spacing of about 20.2 s, rather than continuously optimized switching times.

Success requires a settled final state on the opposite side of the center, lying within 0.05 of the reflected pre-cue x50/L. Both pre-cue and final states must be at least 0.02 from the center, preventing small sign changes in nearly unpolarized states from being counted. The final state must also pass the 0.002 settling tolerance and confirmation rule. Rac and MCA are assessed separately in their respective plots.

The required duration rises with length over much of the plotted range, especially at beta = 1x and 2x. Stronger tension-to-Rho feedback therefore tends to increase the resistance to replacement of an established state in larger cells under this protocol. Cases with no successful duration are omitted, and lines break across missing values; absence does not establish impossibility for durations longer than the tested range.

**Local provenance:** [figure-specific readme](polarity1d/results/data/figure2/README.md), [figure builder](polarity1d/scripts/build_debelly_standalone_figure2_v2.py), [campaign generator](polarity1d/scripts/create_debelly_standalone_figure2_julia_activation_defaults.py), [model assumptions](polarity1d/MODEL_ASSUMPTIONS.md), and [panel source data](polarity1d/results/data/figure2/).

## 3. Gaussian velocity model: literature comparisons of reference and polyploid migration

**Purpose.** This figure asks how much of the reported migration differences can be described by changing two quantities: velocity memory and stochastic forcing strength. Xuan and Mallin provide pairs of fitting summaries; the Jemaà observations provide speed-only constraints. Blue denotes the reference population and orange/red the polyploid-associated population. These labels cover different experimental comparisons, including selected tetraploid clones, PGCCs, and treatment-associated PACCs, and should not be interpreted as interchangeable clean manipulations of genome doubling.

### Model and fitting procedure

Cells move in two dimensions by integrating a velocity that relaxes toward zero with memory time tau and is continuously perturbed by Gaussian white noise of strength eta. In plain-text notation: dx = v dt; dv = -(v/tau) dt + eta dW, with independent noise on each axis. The stationary RMS speed is eta times sqrt(tau), so it is derived from the two fitted parameters. It is not a third free parameter.

Within each condition, every simulated cell has the same parameters. Initial velocities are stationary Gaussian, and there is no directional drift, permanent cell-to-cell heterogeneity, or constant-speed constraint. Position and velocity are sampled using the exact integrated Gaussian process, including the camera interval used by each experiment. This matters because frame-to-frame speed can be much lower than instantaneous speed when velocity memory is shorter than the interval between images.

For Xuan, fitting uses mean sampled speed and one late-lag mean-squared displacement (MSD), near 40,000 s. The simulated observation schedule is 10-min sampling over 12 h. For Mallin, fitting uses mean accumulated path length and mean per-track directness over 24 h with 30-min sampling. Directness is endpoint displacement divided by path length, calculated for each track before averaging.

The objective is a sum of squared log residuals scaled by the squared half-log widths of the supplied reading intervals. These intervals represent tolerances in reading reported figures, not standard errors or confidence intervals. At each tau, eta is profiled analytically. Gaussian moments supply mean speed, path, endpoint, and MSD. Mallin directness uses an interpolated bank of 16,384 simulated cells at each of 161 log-spaced memory values; final checks use 65,536 fresh cells. The memory search spans 0.0005–48 h.

### Panel A — Fitted parameter locations and compatible regions

Separate Xuan and Mallin plots show memory time horizontally and stationary RMS speed vertically, both on logarithmic axes. Stars mark optimal fits. Colored shading marks parameter combinations that satisfy both fitted reading intervals simultaneously; it is a compatibility region, not a statistical confidence region.

For Xuan, the polyploid fit has longer memory and lower RMS speed than the reference fit. For Mallin, the polyploid fit also has longer memory and lower speed, but the reference optimum lies at the lower memory bound. No compatible reference region is found on the plotted grid: its path and directness targets cannot both be brought within their reading intervals with this model.

### Panel B — Restricted comparisons of which parameters differ

The heatmap compares four joint fits for each study: no population difference, only tau differs, only eta differs, or both differ. Numbers are reading-scaled discrepancy scores, with lower values indicating closer agreement. Shared parameters are reoptimized jointly under each restriction.

Allowing both parameters to differ gives the lowest score for both studies. Xuan reaches essentially zero discrepancy for its fitted central values, whereas Mallin retains residual error. This is a descriptive comparison conditional on the tau/eta parameterization. Holding eta fixed does not hold RMS speed fixed, and the scores are not likelihood ratios or significance tests.

### Panel C — Memory, forcing strength, and derived RMS speed

Three grouped bar plots show tau, eta, and stationary RMS speed for the two populations in each study. Xuan's fitted memory increases from approximately 0.107 to 0.725 h, while RMS speed decreases from approximately 32.6 to 17.3 µm/h. Mallin's polyploid fit has a memory of approximately 2.12 h and RMS speed of 8.06 µm/h.

The dagger on the Mallin reference bars is essential: its memory is at the lower bound of 0.0005 h, accompanied by very large eta and instantaneous RMS speed. Those values are unresolved, rather than reliable biological estimates. Very short memory and very large instantaneous fluctuations can compensate to preserve the observed camera-sampled step sizes.

### Panel D — Fitted measurements and a held-out endpoint check

Five small plots compare observed summaries with model predictions: Xuan sampled speed, Xuan late MSD, Mallin accumulated path length, Mallin directness, and Mallin endpoint displacement. Open diamonds and their bars show reported values and reading intervals; colored dots show predictions.

The first four quantities were used for fitting. Mallin endpoint displacement was held out. Xuan's two fitting targets are reproduced, as are the Mallin path means and polyploid directness. Mallin reference directness remains slightly high, approximately 0.144 versus a supplied interval of 0.130–0.140. Predicted Mallin endpoints, approximately 44 and 69 µm, fall inside the supplied intervals. Because endpoint and fitted summaries come from the same experimental tracks, this is a held-out-summary check rather than independent experimental validation.

### Panel E — Xuan MSD curves, early-lag behavior, and crossover

The full predicted MSD curves are shown on logarithmic axes, together with the early- and late-lag observations. Only the late-lag points were fitted. The polyploid prediction starts below the reference at short times but overtakes it later, illustrating how lower instantaneous motility can coexist with larger long-time displacement when memory is longer.

The gray vertical band marks the reported crossover interval of roughly 2–3 h. The model crosses at about 1.05 h, earlier than that interval. At the nominal 600-s early lag, the predicted polyploid MSD is about 7.7 µm², below the supplied 9–15 µm² interval; the source's early time coordinate is only approximately known. These mismatches remain visible checks on the model. With the same tau in both populations, changing eta alone only rescales MSD and cannot produce crossing curves.

### Panel F — Sensitivity of the slow fraction to an unknown threshold

The horizontal axis is a common distance threshold for classifying a one-hour window as slow. The vertical axis is the fraction of time classified as slow across overlapping one-hour windows. Solid curves use net displacement; dashed curves use accumulated camera-sampled path. Colored horizontal bands show the reported slow-fraction intervals.

The panel explores two interpretations of a source criterion described in terms of moving less than one average cell length. The numerical cell length and exact measurement definition are unavailable, so no threshold is fitted and no separate population-specific thresholds are asserted. The result demonstrates sensitivity to the operational definition rather than establishing a unique slow fraction from the fitted parameters alone.

### Panel G — Jemaà speed-only parameter families

Curves show the combinations of memory and stationary RMS speed compatible with the selected clone-speed observations. Solid curves correspond to the 2023 RKO comparison and dashed curves to the 2017 sarcoma comparison; blue and orange retain the reference/polyploid distinction. Thin bands reflect the speed-reading intervals.

For each assumed memory, a different forcing strength can reproduce the sampled speed. The camera cadence is included in this calculation. Consequently, speed alone defines a family of solutions and cannot identify a unique memory or noise strength. These curves are constraints, not uniquely fitted parameter estimates.

### Panel H — Simulated migration tracks at short and long durations

Four trajectory plots show Xuan at 1 and 12 h and Mallin at 1 and 24 h. Each displays twelve cells per population, started at a common origin. The short tracks are exact prefixes of the long tracks, and positions are sampled at the study-specific camera intervals.

The tracks illustrate the consequences of the fitted memory and forcing differences, including more persistent excursions for the longer-memory populations. Differences between individual tracks arise entirely from stochastic motion under common population parameters. These are illustrative model realizations, not measured trajectories or additional fitting data.

**Local provenance:** [model and fitting readme](wgd_langevin/gaussian_literature/README.md), [results report](wgd_langevin/gaussian_literature/REPORT.md), [figure builder](wgd_langevin/gaussian_literature/plot.py), and [observations and result tables](wgd_langevin/gaussian_literature/results/).

## 4. Hosny migration: four representative models with shared knockdown mechanisms

**Purpose.** This figure compares alternative explanations of baseline and knockdown migration in PACC and non-PACC populations. The allowed parameter changes caused by each knockdown are shared between populations, but the fitted parameter values and effect sizes remain population-specific. It examines fit quality, parameter tradeoffs, and predicted motion across four representative mechanism choices.

### Stochastic model and observation process

Position integrates a velocity that decays exponentially between randomly timed kicks. Four parameters describe each state: tau, the velocity-decay memory in hours; lambda, the Poisson kick rate per hour; A, the scale of Rayleigh-distributed kick magnitudes in µm/h; and kappa, the concentration of kick angles around the oxygen-rich direction. Kappa = 0 gives uniform directions; increasing kappa concentrates kicks toward the oxygen source. The directed mean motion emerges from these biased kicks, rather than from an additional fitted drift parameter.

All cells in a state share the same parameters. The simulator integrates velocity decay exactly between events and camera frames, with approximately stationary initialization from past kicks. Simulated tracks last 16 h and are sampled every 30 min. A is a velocity-jump scale, not a stationary speed or velocity standard deviation, and lambda is a physical event rate, not a numerical timestep parameter.

Fits use the 25th percentile, median, and 75th percentile of three trajectory quantities: accumulated path length, directness, and displacement along the oxygen gradient divided by path length. Squared residuals are scaled by observed interquartile ranges, with numerical floors, and averaged across quantiles, metrics, and source datasets. Untreated, HIF control-siRNA, and RhoA control-siRNA observations each constrain one common baseline within a population, while remaining separate equally weighted source datasets. Baseline and both knockdown states are fitted jointly.

The underlying search considered all 15 nonempty subsets of the four parameters for each knockdown, giving 225 HIF/RhoA subset pairs per population, plus a no-change model. It did not include models in which only one knockdown has an empty subset. Matching subset pairs were then combined across populations. A simulation bank with 2,048 cells per grid point supported optimization, using repeated differential-evolution searches and local refinement. Fresh 32,768-cell simulations supplied the frontier scores; the displayed distribution comparisons were regenerated with 65,536 cells per state. The presentation itself introduces no new fitting.

The four highlighted choices are:

- **K = 0:** neither knockdown changes any parameter.
- **K = 2, directionality/directionality:** HIF and RhoA knockdown can each change kappa only.
- **K = 2, directionality/memory:** HIF knockdown can change kappa and RhoA knockdown can change tau.
- **K = 8:** each knockdown can change all four parameters.

K counts additional knockdown parameters per population. Across both populations, the hypoxic comparison therefore has eight baseline scalar parameters plus 2K knockdown parameters. Shared mechanism subsets do not mean shared numerical parameters. Selected parameters may increase or decrease; the directionality/memory example is the unrestricted-sign exhaustive fit, not an older decrease-constrained version.

### Panel A — Fit error versus mechanism complexity

Gray points show the combined fitting error of the shared-mechanism models against K. The gray line connects the empirical Pareto frontier: choices for which another equally or less complex model does not give an equal or better error. Four larger colored symbols identify the selected representatives.

The no-change model has a score of about 0.0585. Allowing each knockdown to change directionality reduces it to about 0.0364. The equally complex HIF-directionality/RhoA-memory model scores about 0.0500 and lies above the frontier. Allowing all parameters to change reduces the score to about 0.0176. These are descriptive discrepancy scores. The figure does not apply a likelihood-based criterion, establish statistical significance, or identify one uniquely proven mechanism.

### Panel B — Parameter values for each representative

The 4-by-4 grid has one representative model per row and one parameter per column: memory, kick rate, kick magnitude, and kick directionality. Within each facet, baseline, HIF-knockdown, and RhoA-knockdown categories contain adjacent blue non-PACC and orange PACC bars. Numerical labels give the fitted values. Every facet has its own linear scale beginning at zero, so bar heights should be compared within a facet and labels used for comparisons across facets.

The K = 0 row repeats the baseline parameters under both knockdowns. The directionality/directionality row changes only kappa; the directionality/memory row changes HIF kappa and RhoA tau. The K = 8 row allows all columns to change under both knockdowns. Common normoxia reference parameters are printed separately above the grid.

The panel also shows how different mechanism assumptions produce different baseline estimates even though they fit the same source summaries. For example, larger memory can trade off against smaller kick magnitude. All bars are fitted point estimates; no parameter uncertainty intervals are shown.

### Panel C — Full distribution-summary comparison with the observations

Five metric columns show accumulated distance, directness, along-gradient displacement fraction, endpoint displacement, and across-gradient displacement fraction. Rows show normoxia and the baseline, HIF-knockdown, and RhoA-knockdown states for each population. Colored bands span model Q25–Q75, with ticks marking medians, in a consistent four-model order. Black symbols and bars show source medians and population interquartile ranges. Distinct symbols keep untreated and the two control-siRNA observations visible even though they share a fitted baseline.

The first three columns are fitting targets. Endpoint displacement and the across-gradient fraction are held-out checks. Along-gradient and across-gradient quantities use displacement projected onto the corresponding axis divided by accumulated path length; the source notes leave the exact normalization uncertain. The panel exposes whether each mechanism reproduces the spread as well as the central tendency, and whether improvement in fitted metrics also carries over to held-out summaries. The bands are population variability, not confidence intervals for fitted parameters or means.

### Panel D — Baseline trajectory realizations

Eight plots show the four representative models in columns, with 1-h trajectories above 16-h trajectories. Blue denotes non-PACC and orange PACC, with twelve cells per population. An arrow marks the oxygen-rich direction. The 1-h tracks are exact prefixes of the 16-h tracks and use the same 30-min sampling interval.

All representatives share the same spatial limits within a duration, allowing a visual comparison of their predicted spread and directional bias. These plots show untreated/shared-baseline populations only. They illustrate how different parameter combinations can produce broadly similar observed-scale movement despite quite different inferred short-time dynamics.

### Panel E — Predicted MSD for baseline and knockdown states

A 2-by-2 grid gives one MSD plot per representative. Each facet contains all six hypoxic states: two populations, each with baseline, HIF-knockdown, and RhoA-knockdown conditions. Population is encoded by blue/orange color and treatment by line style. Curves are calculated analytically from the stationary Poisson-kick model, including its directed motion.

Curves coincide where the assumed mechanism gives identical predictions. For example, all three treatments overlap within each population in the no-change model. Shading marks times beyond the observed 16-h window. No MSD observations were supplied as fitting targets, so the entire panel is a model prediction, with the shaded portion additionally extending beyond the measurement duration.

### Panel F — Baseline PACC/non-PACC ratios across mechanism choices

Four facets show the baseline ratios of tau, lambda, A, and kappa. The horizontal axis is K. Lines follow the shared frontier, and symbols highlight all four representatives, including the off-frontier directionality/memory model. A horizontal reference at one denotes equal baseline parameters between populations.

These ratios reveal which inferred contrasts remain similar across model structures and which depend strongly on the allowed knockdown effects. In the plotted choices, PACC memory is longer and kick rate lower, while the inferred magnitude ratio changes substantially between structures. Directionality is also greater in the PACC baseline fits. Agreement across several structures is descriptive robustness, not proof that a parameter is uniquely identified within any one structure.

### Panel G — Within-model sensitivity to fixed baseline ratios

The four facets examine the same baseline ratios, but only within the HIF-directionality/RhoA-directionality representative. At each horizontal-axis value, one PACC/non-PACC baseline ratio is fixed and the remaining baseline and knockdown parameters are jointly refitted. The vertical axis reports extra fitting error, with the same error scale in all facets. The ratio grid spans 0.03–30 and includes equality at one.

Broad, shallow profiles indicate ratios that can change substantially while other parameters compensate. Sharper increases indicate ratios that are harder to alter without degrading agreement with the supplied summaries. These are descriptive optimization profiles, not confidence intervals or posterior distributions. They do not apply automatically to the other three representative models.

**Interpretive scope.** The analysis uses the supplied Hosny summary observations rather than raw trajectories or replicate-level data. It does not model cell area, spatially varying oxygen response, changing oxygen conditions, or permanent cell heterogeneity. The four mechanisms are alternative conditional explanations of those summaries; improvements in discrepancy do not establish the molecular action of a knockdown.

**Local provenance:** [presentation readme](wgd_langevin/hosny/poisson_velocity/shared_mechanisms/presentation/README.md), [Poisson model and observation process](wgd_langevin/hosny/poisson_velocity/README.md), [exhaustive-search methods](wgd_langevin/hosny/poisson_velocity/exhaustive/README.md), [shared-mechanism and ratio-profile methods](wgd_langevin/hosny/poisson_velocity/shared_mechanisms/README.md), and [presentation result tables](wgd_langevin/hosny/poisson_velocity/shared_mechanisms/presentation/results/).
