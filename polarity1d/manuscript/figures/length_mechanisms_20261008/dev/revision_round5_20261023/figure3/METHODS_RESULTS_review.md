# Figure 3: size-dependent boundary response and shared translation

Generated 2026-10-23T19:03:33.862615+00:00. All native producers finished; this derivative uses existing data only.

Cue placement separates sliding of the incumbent from takeover by a new region; enlargement slows both routes. Weak boundary susceptibility changes differently across models: terminal response at the largest size is 70% of the smallest-size response in Wave-pinning, 9% in Goryachev, 91% in Holmes and 99% in Spring.

The deterministic boundary response predicts physical domain drift without fitting when shared translation is retained. Wave-pinning’s C60/C20 drift ratio at f = 0.001 is 0.494 predicted and 0.473 measured. Holmes and Spring show much weaker susceptibility changes. The independent-edge approximation is refuted by boundary co-motion; Goryachev’s drift prediction remains quantitatively inconsistent and does not explain its full heading persistence.

## Assays and size matching

The nominal weak-pulse peaks are 10⁻⁵ Wave-pinning, 10⁻⁴ Holmes, 4 × 10⁻⁴ Spring and 6.4 × 10⁻⁴ Goryachev, each 10⁻⁴ of its model’s fixed native input reference. Gaussian SD and pulse duration are 1.5 and 50 chemical units, or 0.15 and 5 Spring units. Each starts from the stationary polarized state under the same condition’s mean input, without noise.

A–B use finite-strength continuous cues placed across or away from the incumbent boundary, measuring time to half-turn. C uses a separate weak brief boundary pulse. B uses chemical sizes 30/37.8/47.6 (Spring 2.5/3.15/3.97), whereas C–D use 20/30/40/50/60 (Spring 1/1.5/2/2.5/3). A matched-size link between strong-cue reorientation and small-signal susceptibility is therefore not established; the two assays are complementary and do not uniquely sort the Figure 2 mechanism classes.

## Derivation

The supplementary PDF renders the numbered equations. For a stationary native state, linear perturbations satisfy

$$
\delta\dot{\mathbf{y}}=J\delta\mathbf{y}+B\delta\mathbf{u},\qquad J\mathbf{q}=0,\quad J^{\mathsf{T}}\mathbf{p}=0,\quad\mathbf{p}^{\mathsf{T}}\mathbf{q}=1. \tag{1}
$$

Here **q** is the derivative of the translated complete profile. Conserved-pool adjoints have zero overlap with it. The native spatial shape for the velocity response to input cell j is

$$
\dot{s}=\sum_j g_j\delta u_j,\qquad g_j=(\mathbf{p}^{\mathsf{T}}B)_j. \tag{2}
$$

The native gain under the actual edge-pulse shape $h_e$ is $Q_e=\sum_j g_jh_{e,j}$. The measured finite-window terminal displacement per dose is $M_e$. The kernel used for D is therefore $\widehat g_{e,j}=(M_e/Q_e)g_j$: its amplitude comes from the pulse measurement, and only its unfitted spatial shape comes from the native Jacobian. The raw native gain $Q_e$ is an independent numerical check. A localized pulse alone cannot determine every spatial component, so the native linearization is an explicit additional assumption, not hidden identification from one pulse. The covariance of the actual lognormal forcing is

$$
K_u(j,k;t,t+v)=\mu^2\left\{\exp\left[\sigma^2 C_x(j-k)C_t(t,t+v)\right]-1\right\}. \tag{3}
$$

The periodic Gaussian spatial filter is discretely L2 normalized exactly as in the frozen campaign. Cx is its circular autocorrelation, including the zero spatial mode. Ct sums the four OU-knot covariances weighted by the two variance-normalized interpolation coefficients at each time. We integrate v from zero to 20τ at knot/12 spacing and average 24 midpoint knot phases. The omitted tail is bounded by the exponentially small latent covariance; numerical refinement is checked independently. No Gaussian small-sigma replacement is used at sigma=1. With the convention MSD=2Dt, the marginal edge diffusion is

$$
D_e=\sum_{j,k}g_jg_k\int_0^{\infty}\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v. \tag{4}
$$

Both zero-frequency edge kernels coincide because both boundaries share translation. Transient width and shape changes decay and are not assigned a nonzero asymptotic diffusion. The requested approximation and its covariance-complete counterpart are

$$
D_s^{\mathrm{ind}}=\frac{D_L+D_R}{4},\qquad D_s=\frac{D_L+D_R+2D_{LR}}{4}. \tag{5}
$$

The prescribed independence assumption therefore undercounts the shared neutral translation. The main panel shows both the shared-translation prediction and the rejected independent-edge candidate. Retaining the cross term uses the same native model and measured pulse gain; it is not a fitted factor-of-two correction. The integrated-displacement summary is explicitly finite:

$$
\mathcal A_{20\tau}=\frac{1}{a\tau}\int_0^{20\tau}\left[x_e^{\mathrm{pulse}}(t)-x_e^{\mathrm{control}}(t)\right]\,\mathrm{d}t. \tag{6}
$$

It cannot be replaced by an integral to infinity when the pulse leaves a permanent position shift.

The kernel used for each marginal diffusivity is explicitly measurement-scaled:

$$
\widehat g_{e,j}=\frac{M_e}{Q_e}g_j,\qquad M_e=\frac{\Delta x_e(20\tau)}{a\tau},\quad Q_e=\sum_jg_jh_{e,j}. \tag{7}
$$

Equation (4) uses this measured-gain kernel in place of the raw native g. Neither $M_e$ nor the normalized spatial shape uses observed drift.

## Numerical preparation and checks

Original Figure 4 basal states are copied with their hashes. Each is allowed to settle at its exact condition mean for 16000 chemical units or 1600 Spring units. Eligibility requires a unique half-contrast region, contrast ≥0.05, max absolute native RHS ≤1e-6, and tail change ≤1e-4 of max(1, maximum state). Pulses use BDF rtol 2e-10, atol 2e-12, exact restart at withdrawal. Both edges receive independent half/nominal/double pulses; a matched unpulsed continuation removes any preparation drift. Current half-contrast crossings use periodic linear interpolation. Spring positions are material fractions times mean-input circumference; pure extension is excluded.

Native sparse chemical Jacobians use centered numerical differentiation; Spring uses its complete analytic native Jacobian. Translation eigenvectors are selected by overlap with the native profile derivative; adjoint conservation modes are excluded. The spatial input derivative is measured directly from the unchanged native RHS. Every comparison table retains linearity, pulse validation, eigenvalue and stationary residuals. Amplitudes and horizons were fixed before diffusion comparisons.

## Independence diagnostics and precision

The pre-existing float16 exemplars are the only saved profiles; at most 400 bins and 2000 frames each. Current half-contrast edges are identified only when exactly one region exists and current contrast is at least 0.1 of the original native basal contrast. Displacements are unwrapped within uninterrupted unique-region segments; no interpolation joins segments. For every lag we report left/right MSD, cross-product, centered increment correlation, and actual pair count. We retain the empirical cross-product in the variance identity, including systematic increments, rather than replacing it by the centered covariance. Filled supplement points require both MSDs greater than ten times a conservative squared float16 positional error bound, derived from half ULP of profile values and local crossing slope. Low-input short-lag results often fail this precision gate and remain open symbols. Two exemplar seeds describe covariance but do not establish a campaign confidence interval.

## Undefined cases and scope

- goryachev C=20, f=0.01: No unique polarized half-contrast boundary at the declared mean input
- goryachev C=30, f=0.01: No unique polarized half-contrast boundary at the declared mean input

No missing execution is represented as censoring. Goryachev tracked drift can be very small and retains substantial disagreement. None of this equates tracked domain motion with whole-heading decay or with empirical migration persistence. Length represents WGD-associated enlargement; observed migration summaries and conditional OU fits remain distinct from these 1D polarity readouts.

The maximum change in any terminal response over its final 0.5τ window is below 0.026% of the condition’s maximum terminal displacement; the transient is therefore settled within the observation window while a neutral positional offset remains. The covariance quadrature refinement changes the integrated kernel by 0.0026% (chemical) and 0.0007% (Spring). The declared interpolation differs from a continuous OU covariance and is retained exactly. Numerical checks are saved under `supporting/metadata/numerical_checks.json`.

## Fixed-input size test

For each model and input, both observed and predicted diffusivity are divided by their own value at that model’s fixed smallest grid size. The same reference size is retained at every input. The resulting 115 defined prediction pairs are exported with observed bootstrap ratio intervals in `supporting/data/size_normalized_prediction_comparison.csv`. Five Goryachev ratios at f=0.01 are unavailable because the C20 mean-state response is undefined; C40 is not substituted. At f=0.001, Wave-pinning’s C60/C20 ratio is 0.49390049 predicted and 0.47324320 measured. The size-normalized candidate ratios coincide to numerical precision because the rejected independence assumption changes the absolute diffusivity by essentially one common factor; therefore the size subrow shows the covariance-complete prediction, while the absolute subrow retains both candidates.

C insets use the finite integrated response through 20τ divided by the smallest-size response, with the same zero-based 0–1.1 range for all models. Their source is `supporting/data/normalized_integrated_edge_response.csv`. These insets compare integrated displacement; the reported 0.701/0.914/0.989/0.0878 folds compare terminal displacement per pulse dose. They are distinct summaries.

## Relaxation row and historical evidence

For A’s cue-free row, the cue is switched off at 95% heading or at 20,000 if 95% is not reached. Cue-free relaxation ends at the native stationarity tolerance or at four times the cue-on duration. Example i withdraws at 790.5 and relaxes 258.6; example ii withdraws at 227.6 and relaxes 291.3. Both reach the stationarity criterion 2.2683×10⁻⁷ in active-species RHS. Profiles and all A/B numerical coordinates remain unchanged.

The reviewed main and complete supplement at commit 49754ab are preserved unchanged and appended after the new supplementary pages; their earlier single-candidate interpretation is historical. The current main exposes both candidates and the fixed-input size test. All covariance integration, conservation and amplitude checks remain valid because no scientific numerical source changed.
