# All-five linear gain and incumbent-repair assays

This is a bounded design/evidence checkpoint, not a launched campaign. Figure 1 remains linear 1D and includes Wave-pinning, Goryachev, corrected deBelly, Holmes 3 and Spring. The same definitions should subsequently be applied to supported WP/G/H rings. No new integrations or figures were generated for this design.

**Recommendation.** Extend B to all five using cue-weighted incremental native activation gain. Replace the old inactive-pool hold/return assay in main-panel C with a common, small, spatially targeted active-site deactivation followed by paired native recovery. Keep the old WP/G hold/return evidence as support. Perturbing deBelly MCA would ask a different question and should not fill its Rac-repair entry.

## B: incremental response to the model's native activation channel

At an exact qualified basal checkpoint, let a unit-peak spatial cue profile be \(p(x)\), and let \(s\) be its amplitude in that model's native input units. The extra initial production of the active observable is linear in this input in all five models:

\[
\left.\frac{\partial \dot a(x)}{\partial s}\right|_{s=0}=p(x)\,g(x),
\qquad
G=\frac{\int p(x)\,g(x)\,\mathrm dx}{\int p(x)\,\mathrm dx}.
\tag{1}
\]

| Model | Active observable and native input channel | Exact local gain \(g(x)\) | Source |
|---|---|---|---|
| Wave-pinning | Active u; additive activation coefficient | Inactive v | `code/core/models/wave_pinning.py:34–42` |
| Goryachev | Active u; extra inactive-to-active transfer coefficient | Inactive v | `code/core/models/goryachev.py:51–56` |
| Holmes 3 | Active Rac; additive Rac GEF activation | Inactive Rac / 7.5 | `code/core/models/holmes_model3.py:47–58`; Rac total in `code/core/src/polarity1d/holmes.py:32–53` |
| Spring | Active u; additive basal activation b | Inactive v | `code/core/models/spring.py:10–11,103–121` |
| Corrected deBelly | Rac; additive Rac activation parameter alpha | 0.04 / (1 + Rho²) | `code/core/models/debelly.py:344–363,806–815` |

deBelly has no explicit inactive-Rac field or conserved Rac pool. Its gain measures Rho-dependent attenuation of native input, not substrate depletion. Rac/Rho/MCA fields must not be relabelled as interchangeable resources. Spring's input is a conservative u/v transfer on its moving material domain; its current physical domain length is used for cue integration.

Use the same cue-weighted definition for linear and ring results, retaining the native input channels. Ratios to each model's own declared reference size compare size dependence. Absolute gains are not a common biochemical calibration or a prediction of challenger-establishment probability. A gain measured at cue onset is an incremental source, not total initial growth.

### Available evidence: all 97 linear checkpoints

The [algebraic extraction](../../../../code/mechanisms/all_five_linear_design_evidence_20261017.py) checked every retained checkpoint hash and independently verified the native increment identity; maximum discrepancy was 2.67 × 10⁻¹⁵. It uses a fixed end-to-zero raised-cosine support of 3, or 0.2 for Spring. Both last-cell and cue-weighted gains are retained in [the CSV](supporting/data/native_gain_design_evidence.csv), with [verification and source hashes](supporting/metadata/native_gain_design_evidence.json).

| Model | Qualified linear sizes | Highlighted pair | Larger/smaller cue-weighted gain |
|---|---:|---:|---:|
| Wave-pinning | 20 | 19 → 39 | 0.99999999 |
| Goryachev | 21 | 19 → 39 | 0.65361060 |
| Corrected deBelly | 12 | 19 → 39 | 0.84569943 |
| Holmes 3 | 23 | 20 → 40 | 0.99984253 |
| Spring | 21 | L₀ 1 → 2 | 1.00422832 |

The WP/G values agree with the existing point-gain panel to the shown precision because inactive substrate is essentially flat over the patch. deBelly supplies an additional nontrivial size response: approximately 15.4% lower native input gain over 19 → 39, mediated directly by the basal Rho field. This does not prove weaker finite challenger activation or identify its full migration consequences.

### Linear-to-ring comparison must respect the factor of two

The linear cue implementation uses end-to-zero width \(w\) and \(\cos(\pi d/w)\); the ring implementation calls its **total** support \(W\) and uses \(\cos(2\pi d/W)\). Even reflection of a length-L Neumann interval onto a C = 2L ring therefore requires W = 2w, unchanged peak, mirrored complete state and twice the integrated cue. The current ring width-3 campaign corresponds to linear end support 1.5, whereas the linear fixed-support-3 assay reflects to ring total width 6. Source: `code/core/polarity_establishment.py:41–57` and `code/mechanisms/ring_distance_prototype_20261017.py:97–111`.

Do not silently call the present width-3 ring challenges exact equivalents of linear width-3 challenges. Algebraic gain can be re-extracted at either width without integration. Existing outcome histories remain valid at their declared width. An exact reflection test additionally requires matching cell-centered meshes (ring 2N cells), mirrored **all** chemical fields and cue/preparation histories; separately prepared basals are not automatically pointwise identical. This relation applies to the native chemical WP/G/H models, not deBelly's mechanics or Spring's endpoint-force law.

## C: a small incumbent-targeted deactivation and paired recovery

The shared measured question should be: **after the same fractional local loss of incumbent activity, how much of that loss does the native coupled system restore over an early declared window?** This is a perturbation-specific recovery measurement, not generic stability, force, or whole-cell movement.

Freeze the original basal incumbent center, background \(b_0\), half-contrast core \(\Omega_0\), and its radius/endpoint width \(w_0\). On a ring use the shortest distance \(r\) from the frozen incumbent center; on a linear interval use distance from its incumbent endpoint. Define a compact taper and a fixed activity-removal direction from the **original basal** profile:

\[
\chi_0(r)=
\begin{cases}
\tfrac12\left[1+\cos(\pi r/w_0)\right],&0\leq r<w_0,\\
0,&r\geq w_0,
\end{cases}
\qquad
h_0(x)=\chi_0(r)\,[a_0(x)-b_0].
\tag{2}
\]

For a ring with unequal left/right basal half-contrast radii, use the corresponding frozen radius on each side. This preserves the same definition without forcing symmetry. The taper and its first derivative vanish at the core edge. This avoids a hard-mask concentration jump and its artificial high-frequency diffusion transient. The native basal profile and all geometric definitions are recorded, rather than reconstructed from the challenged profile.

At the assay state \(Y_*\), apply one instantaneous kick

\[
a_{\mathrm p}(x,0)=a_*(x)-\varepsilon h_0(x),
\qquad \varepsilon=0.005,
\tag{3}
\]

and add \(\varepsilon h_0\) pointwise to the corresponding inactive pool in WP/G/H/Spring. The complete other native fields and controller state are copied unchanged. Holmes Rho and inactive Rho are untouched, as is Spring's pair of moving endpoints. Corrected deBelly changes Rac only: its Rac equation has no tracked inactive pool; Rho, bound/unbound MCA, displacement, velocity, tension and all four endpoint-controller fields remain the exact same initial values. It is not appropriate to invent a Rac pool or compensate through MCA. The deBelly active-observable loss is a model-specific intervention on a nonconserved field; this difference must be stated in the caption/methods.

This frozen, compact direction is crucial for the later ring heatmaps: a full **current** excess-profile reduction would damage the challenger too. Equation (3) instead targets the original incumbent at cue onset or a specified challenged time. Before any run, require \(a_*\geq\varepsilon h_0\) pointwise and finite nonnegative fields. If it fails, report that state/amplitude as unavailable; do not clip, cap, move the core, or replace the direction. Cue/core overlap exclusions remain explicit.

Evolve the perturbed and exact same-start unperturbed sham under unchanged native coupled dynamics and the **same external cue schedule**. Basal Figure 1 C has no challenger; a challenged-state ring C map can start, for example, 10 native units after cue onset while preserving that cue. No parameter or field is frozen during recovery. No material is subsequently returned or state reset.

Define the integrated deficit within the frozen original core and the finite-window repair coefficient:

\[
D(t)=\int_{\Omega_0}[a_{\mathrm s}(x,t)-a_{\mathrm p}(x,t)]\,\mathrm dx,
\qquad
K_{\tau}=\frac{D(0)-D(\tau)}{\tau D(0)}.
\tag{4}
\]

Positive K means the perturbation-induced deficit shrinks relative to its paired sham; negative K means it grows. Preserve the actual sham and probe amount changes separately. A positive K can coexist with incumbent loss in both trajectories. It is a finite-window coefficient, not an instantaneous Jacobian eigenvalue or a proof of incumbent maintenance against an actual cue.

For WP/G/H/deBelly use native finite-volume fractional edge-cell integration. For Spring freeze the core in material coordinates, integrate each evolving active profile with its own physical length, and include dilution/domain changes rather than fixing the domain artificially. This follows its native concentration-to-amount convention (`code/core/models/spring.py:23–35,45–49,103–121`; incumbent integration in `code/mechanisms/initial_polarity_report.py:27–35`).

### Early window: choose once per model from a declared pilot

Do **not** impose the same numerical time across uncalibrated models. Goryachev can rapidly recapture an active-to-inactive kick through its u²-dependent conversion; a long common window could saturate its recovery and conceal size dependence. deBelly instead has a discrete, staggered native 0.5-s update. Its response must be measured by actual paired native steps, not a fabricated common continuous RHS.

Before the full grid, run a bounded **15-case pilot**: shortest, middle and longest eligible linear basals for each of the five models. Record full short deficit curves for sham and epsilon = 0.001, 0.005, 0.01; use tight paired continuous integrations and the retained native deBelly step. The pilot's purpose is to choose one early window per model that remains fixed across geometry, size, separation and assay stage. Stop after this report and window-selection decision; do not launch the full grid automatically. The subsequent ring port must check the smallest/largest supported rings and a representative early challenged state against the frozen early-window criterion before expanding the maps, since its size range extends beyond the linear campaign.

Prespecify candidate horizons spanning native 10⁻⁶–1 for WP/G/H/Spring, using a logarithmic sequence; for deBelly use exact native times 0.5, 1, 2, 5 and 10 s. Choose the **largest** candidate at which every pilot condition has resolved paired error and \(|1-D(\tau)/D(0)|\leq0.2\). This rule prevents complete repair from being misreported as weak size dependence and also excludes large early amplification from the linear-response scalar. Do not choose separate horizons for different sizes or separations. Require the epsilon-normalized response to agree across the three amplitudes within a declared tolerance, suggested 5%; if the gate fails, retain the failure and do not tune an individual condition. Near-zero responses need numerical bounds rather than a minimum-effect selection rule.

Use the pilot only to fix the clock, not to select a favorable result. Freeze the window and pilot decision record before the full Figure 1/ring campaigns. Display representative complete short recovery curves in support so the scalar's relation to the dynamics is visible. If no candidate meets the gate, report that common early response as unavailable rather than silently moving to a different assay. An unresolved zero or sign change cannot be displayed as a positive fold.

For continuous solvers check paired error against tighter tolerances and K at adjacent early horizons. Preserve the finite-window definition even when its early values converge. For deBelly, halving dt without changing endpoint relaxation is **not** a temporal-convergence check: `debelly.py:621–625` rescales the 0.98 relaxation only for the separate clutch construction. Any explicit corrected-linear refinement would need relaxation \(0.98^{\Delta t/0.5}\), its own declared implementation, and comparison to the retained dt = 0.5 response. No such refinement is required merely to claim the native finite-step result, but that claim must stay discrete.

### Normalization and display

Keep native \(K_\tau\), \(D(0)\), \(D(\tau)\), the chosen window, epsilon, numerical error bounds and raw sham/probe changes in exported data. Comparisons across models concern size trends and behavior under the same intervention question; they do not equate native clocks. Positive within-model reference ratios may be plotted with a shared linear/log axis, depending on the verified range. If the reference is near zero, negative or changes sign, use signed native coefficients with explicit model scales or an additional dimensionless repair fraction; do not force a logarithmic fold panel.

## What can be reused, and what requires fresh integration

| Item | Existing evidence | Required new work |
|---|---|---|
| All-five Figure 1 A and native basal profiles | 97 retained qualified checkpoints and existing organization metrics | No new preparation; preserve exact state/hash and unavailable lengths. |
| All-five B | Algebraic source derivatives verified here for every checkpoint | Render all-five curves; optional exact cue-width/reflection re-extraction uses existing states only. |
| All-five C under the common active-site assay | No existing campaign provides Equation (3)–(4) for all five | Bounded clock/linearity pilot, then native paired short integrations from retained checkpoints. |
| WP/G old inactive hold/return | Fully audited established interval evidence | Retain as support; do not rename its coefficients as the new C. |
| Ring organization and B | Existing native-ring basals and cue-weighted gain maps | Reuse at declared widths; exact reflection requires a specifically matched source state/mesh and width accounting. |
| Ring C size/separation heatmaps | Existing maps use a different inactive-pool hold/return disturbance | Fresh Equation (3)–(4) short paired runs from exact recorded challenged states. Basal C alone has no separation variable. |
| Later challenger fates/kymographs | Existing complete challenge histories | Reuse with exact highlighted state/map links, including a sole cue-associated peak example. These fates are not determined by the new coefficient alone. |

The design does not launch this work. A subsequent implementation should freeze inputs and sources, audit conservation (including the two Holmes pools and Spring physical amount), validate zero-kick identity, preserve all deBelly controllers, record numerical exclusions, and use SLURM for the combined campaign. New large outputs, if they reach approximately 2 GiB, must be measured and added to the root data manifest before handoff.
