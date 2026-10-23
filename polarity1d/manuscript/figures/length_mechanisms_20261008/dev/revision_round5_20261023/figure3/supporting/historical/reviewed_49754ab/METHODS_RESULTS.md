# Figure 3: deterministic boundary response and no-retuning drift

Generated 2026-10-23T18:40:00.249716+00:00; deterministic and profile-analysis producers complete.

Figures A/B retain exactly the preceding quantitative source files, with source SHA-256 checks in `supporting/metadata/preserved_AB.json`. C uses Figure 4’s size grid, not B’s three sizes. No native parameters or noise realization changed.

C. Control-subtracted displacement of the pulsed falling half-contrast edge, divided by peak input a, at f = 0.001. Shading denotes Figure 4 circumference: chemical C = 20, 30, 40, 50, 60; Spring resting C = 1, 1.5, 2, 2.5, 3. The stationary polarized state at each declared constant mean is used, without noise or retuning. A Gaussian pulse centered on the boundary has spatial SD ℓ = 1.5 and duration τ = 50 (Spring 0.15 and 5). Peak a is 10⁻⁴ of the native reference: 10⁻⁵ Wave-pinning, 10⁻⁴ Holmes, 4 × 10⁻⁴ Spring and 6.4 × 10⁻⁴ Goryachev. The grey interval is pulse-on. Inset A₂₀τ is displacement integrated through 20τ divided by pulse dose aτ. Both edges, a/2 and 2a are checked in the supplement. A residual positional offset is allowed because translation is neutral; this is not an infinite displacement integral. Largest/smallest terminal displacement per pulse dose is 0.701 Wave-pinning, 0.914 Holmes, 0.989 Spring and 0.0878 Goryachev.

D. The prescribed independent-edge prediction, Dₛ = (Dleft + Dright)/4, against physical tracked domain drift measured in Figure 4 at f ≤ 0.01 (markers). Colors and shading match C. Each measured terminal edge gain per pulse dose scales the native spatial response shape, normalized to give unit response to that exact localized pulse; this extra spatial-shape assumption is explicit. The native linearization supplies its shape only, while measured C amplitudes supply its gain; the imposed discrete spatial covariance and phase-averaged OU-knot interpolation enter the exact lognormal covariance μ²[exp(σ²CₓCₜ) − 1], with σ = 1. There is no diffusion fitting or parameter adjustment. Line: 1:1; horizontal bars: Figure 4 whole-seed-bootstrap 95% CI. This predicts Figure 4’s tracked drift, not total heading decay. The independent prediction/observation medians are Wave-pinning 0.511, Holmes 0.540, Spring 0.534, Goryachev 0.437. Goryachev spans 0.00359–4.01 across conditions, so a single bias does not summarize its failure. Goryachev C20/C30 at f = 0.01 has no polarized boundary at the declared mean, so its prediction is undefined.

The edge-independence supplement measures simultaneous half-contrast displacements from existing Figure 4 exemplar profiles, seeds 0 and 1, using uninterrupted single-region intervals and exposing float16 precision limits. At quarter-window lag, resolved median correlations are Wave-pinning 0.980, Holmes 0.983, Spring 0.997, Goryachev 0.759. Positive common motion violates edge independence. Retaining the native shared translation covariance doubles the prediction; this sensitivity is supplementary and does not adjust any parameter. Goryachev’s disagreement and short valid-region exposure remain explicit.

Cellular counterpart of C: photoactivation straddling a leading-edge boundary, tracking both flanking boundaries and their displacement after withdrawal. Native units are not calibrated between models. deBelly is excluded from rings.

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

The prescribed independence assumption therefore undercounts the shared neutral translation. Keeping covariance is a separate, untuned sensitivity; we retain the failing independence prediction in the main panel. The integrated-displacement summary is explicitly finite:

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
