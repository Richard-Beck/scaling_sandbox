# Common activity assay: independent design and evidence review

Reviewed 17 October 2026. Scope: the all-five linear assay [design](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017/DESIGN.md), its gain CSV/JSON, the native source equations and retained basal checkpoints. This review performed algebraic and hash checks only. It did not integrate trajectories, generate figures or alter retained evidence.

**Verdict: implement the bounded 15-case pilot now, after freezing the narrow definitions below.** The intervention is scientifically coherent and portable enough to test. A general stability analysis, mesh campaign, parameter search or deBelly ring development is not a prerequisite. Stop after the pilot report and one model-specific clock decision; the pilot does not authorize the full linear or ring grid.

The assay measures the response to a declared local activity-loss direction. It can support an enlargement-dependent change in that response. It cannot, by itself, establish greater incumbent maintenance under a challenger, generic stability, a causal effect of incumbent amount, or a migration-persistence mechanism. The proposed narrative should remain incumbent enlargement, model-specific input/recovery consequences, and location as a modifier.

## Essential definitions to freeze before the pilot

1. **Use the existing basal background and core exactly.** `initial_polarity_report.measure` defines the background as the minimum of the original basal active field, the threshold halfway between its minimum and maximum, and the core edge by the first boundary-connected threshold crossing, interpolated between native sample positions. Reuse its recorded background and width, and independently reproduce them from the checkpoint. Freeze these quantities in the basal coordinates. The tapered direction is nonnegative in these basals; it does not alter tiny inactive-end or active-tail values outside its support. On a challenged ring, freeze the original incumbent center and the two original basal core edges; do not recenter on whichever site subsequently wins.

2. **Describe epsilon accurately.** The perturbation removes epsilon times a tapered basal excess profile. It is approximately a specified peak/local fractional excess attenuation; it is not an identical integrated fraction of core activity across models or sizes. Record initial removed amount, its fraction of basal core activity and its fraction of basal excess activity. Normalization by the actual initial deficit makes the coefficient interpretable without claiming equal whole-core damage.

3. **Specify the numerical and amplitude gates before seeing recovery results.** “Resolved paired error” and “suggested 5%” are not executable rules. A concrete pilot rule is provided below. Near-null responses should be retained with a numerical resolution label, rather than rejected for lacking an effect. An all-five panel must not acquire its membership by selecting the models or sizes with convenient positive coefficients.

4. **Apply clock gates to the entire recorded prefix, not just its endpoint.** Choose the largest candidate for which every sampled time through that candidate passes the amplitude and nonsaturation gates. An endpoint-only rule can accept a late return after an earlier saturation or overshoot. Freeze the exact logarithmic candidate list: the root decision is 25 logarithmically spaced candidates from 10⁻⁶ through 1 for the continuous models. For deBelly retain the stated 0.5, 1, 2, 5 and 10 s candidates and record every native 0.5-s step. No size, separation or challenged state gets its own optimized clock.

5. **Implement Spring as a physical-amount comparison between two evolving domains.** Let the frozen material core end at \(q_0=w_0/L(0)\). Its deficit is

\[
D(t)=L_{\mathrm{s}}(t)\int_0^{q_0}a_{\mathrm{s}}(q,t)\,\mathrm{d}q
     -L_{\mathrm{p}}(t)\int_0^{q_0}a_{\mathrm{p}}(q,t)\,\mathrm{d}q.
\tag{1}
\]

Use the native nodal piecewise-linear/trapezoidal convention, inserting the interpolated value at \(q_0\), as in `initial_polarity_report.amount`. Multiplying a difference array by one trajectory's length is incorrect once the endpoints diverge. Copy the two endpoints unchanged at the kick, preserve pointwise active-plus-inactive concentration there, and check the full physical active-plus-inactive amount during recovery. Do not use finite-volume cell weights for Spring.

6. **Preserve the native restart and update contracts.** deBelly needs all eleven restart entries, including all four endpoint-controller values. Only Rac changes initially; signed velocity, displacement and tension require finiteness, not nonnegativity. Nonnegativity applies to chemical concentrations. Use the corrected independent-endpoint model with its retained 0.5-s stepping. Output times must not split those native steps. A smaller step without rescaling the 0.98 controller relaxation changes the physical model and is not a convergence test. No such refinement is required for this bounded native-step pilot. Continuous paired runs need tight tolerances exposed explicitly; the existing convenience wrappers' defaults are not a substitute for the proposed error check.

## A concrete numerical rule that admits null responses

Define the signed dimensionless local repair fraction and the reported coefficient as

\[
F_\varepsilon(t)=1-\frac{D_\varepsilon(t)}{D_\varepsilon(0)},
\qquad
K_\tau=\frac{F_\varepsilon(\tau)}{\tau}.
\tag{2}
\]

For continuous models, run the same paired problem at two frozen solver settings, for example relative/absolute tolerances 10⁻¹⁰/10⁻¹² and 10⁻¹²/10⁻¹⁴. These are proposed pilot settings, not demonstrated accuracy. Both runs must use identical raw input states and exact output times. For each amplitude and sample, a conservative empirical tolerance-sensitivity estimate can be

\[
E_\varepsilon(t)=\max\!\left(
2\left|F_\varepsilon^{\mathrm{tight}}(t)-F_\varepsilon^{\mathrm{tighter}}(t)\right|,
\frac{2|D_{\mathrm{zero}}(t)|}{D_\varepsilon(0)},
100\epsilon_{\mathrm{mach}}\frac{M_{\mathrm{core},0}}{D_\varepsilon(0)}
\right).
\tag{3}
\]

Here the zero-kick path must reproduce its sham, and the final term supplies a floating-point scale. This is an empirical numerical-error estimate, not a certified mathematical error bound. For the native deBelly update, preserve the exact discrete map and replace the tolerance-sensitivity comparison with a comparison against independent dense solves of the same native matrices at exactly 0.5-s steps. Do not pretend that a continuous-solver tolerance experiment applies to it.

The dense arithmetic reference is acceptable with a narrow validation: independently check band-to-dense indexing on a hand-checkable matrix; copy the matrix and right-hand side before a native solve can overwrite them; check the scaled residual of every native/reference solve, with a declared target such as 10⁻¹²; and record the differences in all eleven per-step fields/controller values as well as the resulting paired repair fractions. Replace only the linear algebra callable and restore it reliably after each reference run. Keep all native thresholds, clipping safeguards, MCA closure and update ordering identical. Hash the wrapper and sources. This checks arithmetic of the same discrete model, not temporal convergence or biological validity.

Require \(E_\varepsilon\leq10^{-5}\) in repair-fraction units. Retain responses compatible with zero; for example, flag \(|F_\varepsilon|\leq3E_\varepsilon\) as sign-unresolved. Do not form a positive fold from such a reference. Across epsilon = 0.001, 0.005 and 0.01, require each pair to satisfy

\[
|F_{\varepsilon_i}-F_{\varepsilon_j}|
\leq0.05\max\bigl(|F_{\varepsilon_i}|,|F_{\varepsilon_j}|\bigr)
     +3(E_{\varepsilon_i}+E_{\varepsilon_j}).
\tag{4}
\]

This tests a 5% response agreement where measurable and permits responses indistinguishable within numerical error. It imposes no minimum biological effect. For all three amplitudes and all 15 pilot cases, require

\[
\max_{0<t\leq\tau}\bigl(|F_\varepsilon(t)|+E_\varepsilon(t)\bigr)\leq0.2
\tag{5}
\]

over the recorded prefix. Choose the largest declared candidate passing these rules, separately for each model. Export every candidate and failure; freeze the decision before the complete grid. The pilot clock does not prove that untested states pass: later linear/ring cases that violate the same gate must remain explicit failures or unavailable entries, without tuning their individual clocks. A common model-specific clock may be reconsidered only as a new declared protocol, not silently changed within a map.

## What the coefficient actually measures

Paired shams are essential. A positive coefficient means that the perturbation-induced deficit in the frozen core shrinks relative to the same evolving, cued sham. The two incumbents can both lose activity while that happens. This is the correct reason to replace the old main-panel inactive-pool hold/return coefficient, while retaining that older intervention as mechanistic support under its own name.

The fixed-core measurement includes native transport and, for Spring, physical length changes. Activity redistribution can reduce a local deficit without restoring the whole-domain active amount. Export full-domain active deficit and outside-core deficit from the same saved states, together with the sham/probe core amount changes. These inexpensive diagnostics prevent a local response from being described as whole-activity restoration. A reaction/transport budget would be useful if a later mechanistic claim needs it; it is not needed to launch this pilot. The basal-width-scaled direction also means that a size trend is a response to a consistently defined family of perturbations, not an isolated causal effect of incumbent amount.

The taper removes the hard-mask jump and has zero value and first derivative at the continuum core edge. It does not abolish all finite-grid or diffusive effects. Use raw initial native arrays, record the discrete removed amount, and avoid relative error tests against nearly zero tails. Preserve native internal positivity safeguards, while introducing no additional clipping of the proposed kick or solver output. At a fractional finite-volume edge, the implemented cell-average perturbation and fractional-core measurement should be recorded explicitly; their full-domain and core initial losses need not be numerically identical. Do not clip or redistribute a failed kick. At a challenged state, test the frozen direction's positivity first and retain the unavailable status if it fails.

## Read-only evidence audit

All eight hashes recorded in `native_gain_design_evidence.json` match the current named files. All 97 CSV checkpoints match both their recorded hashes and the basal metric records. Independent re-evaluation reproduces every cue-weighted gain and the native increment identity, with maximum discrepancy 2.6645352591003757 × 10⁻¹⁵. This verifies the recorded algebraic extraction; the eventual pilot still needs its complete runtime dependency/source freeze.

| Model | Records | Native frozen-state gain | Highlighted larger/smaller ratio |
|---|---:|---|---:|
| Wave-pinning | 20 | Inactive v | 19→39: 0.9999999869 |
| Goryachev | 21 | Inactive v | 19→39: 0.6536106032 |
| Corrected deBelly | 12 | 0.04 / (1 + Rho²) | 19→39: 0.8456994345 |
| Holmes 3 | 23 | Inactive Rac / 7.5 | 20→40: 0.9998425251 |
| Spring | 21 | Inactive v | L₀ 1→2: 1.0042283216 |

The source terms have the correct input identities. In Holmes the external input is added after the endogenous Rho-inhibition factor, so the external gain is inactive Rac / 7.5, with no extra inhibition denominator. deBelly's algebraic value is the frozen-state Rac production susceptibility; its actual next native step also includes staggered mechanics/MCA/Rac/Rho updates. Do not rename the algebraic susceptibility a measured one-step growth rate. Spring's cue is an active/inactive transfer and its gain uses the current physical domain.

Using the proposed basal direction, every retained checkpoint remains chemically nonnegative even at epsilon = 0.01. The smallest Goryachev active tails are approximately 5.59 × 10⁻⁴¹. Its longest L = 39 core contains only four sample centers; preserve that case and mark its representation. Direct rechecking gives 55 centers at L = 3 (width 0.81404), eight at L = 19 (width 0.71419), and four at L = 39 (width 0.79923). The original review mistakenly attached the four-center flag to the shortest size; the attribution is corrected here. This observation does not justify an advance mesh campaign, but a size trend dominated by the coarsely represented long-cell profile would need a targeted follow-up before a continuum-shape interpretation.

The shortest/middle/longest eligible pilot sizes, using the existing lower-middle index convention, are WP 5/19/39, Goryachev 3/19/39, deBelly 18/27/39, Holmes 10/28/50, and Spring L₀ 0.6/1.6/2.6. These are 15 starting states and three probe amplitudes per state, not 15 total trajectories. Freeze the exact list in the pilot manifest.

## Ring portability and figure claims

For the native WP/G/H chemical equations, even reflection of an N-cell Neumann interval gives a 2N-cell periodic ring with circumference twice the interval length. Reflect every field and the cue/preparation history. Whole-incumbent width and amount double; peak, mean, normalized gain and the paired coefficient agree under this exact correspondence. Reflect the deactivation direction and core too, so both the deficit and its initial value double. The one-sided linear cue's end-to-zero support w becomes ring total support 2w with the same peak. Existing ring total-width-3 histories are valid at their declared width, but are not exact counterparts of linear support 3.

One further discrete detail matters: `ring_cue` currently centers its peak on a cell center, whereas an interval end incumbent reflects about a boundary face. For an exact identity test, construct mirrored state/cue arrays by concatenating each interval array with its reverse, or explicitly supply the corresponding face-centered cue. Merely setting circumference = 2L, cells = 2N and width = 2w in the current cell-centered cue helper leaves a half-cell offset. This is an implementation detail for the exact comparison, not a defect in existing ring experiments.

Basal recovery has no future separation variable. A separation map must use the same declared early challenged time and cue schedule, preserve the exact full starting state for each location comparison, and retain the basal coefficient as a reference. The frozen incumbent direction avoids direct damage to a remote challenger, subject to declared overlap exclusions. Check the supported smallest/largest ring states and an early challenged state against the frozen clock before expanding the map. The fate examples and global activity vector remain separate observable evidence; neither follows mathematically from the new scalar.

The sensible all-five main-text output is A: incumbent organization, B: native incremental input gain, C: local conditional recovery, with native coefficients and genuine nulls available in support. Shared within-model normalization is useful only for well-resolved positive references. Very small differences in WP/Holmes gain should remain visibly small, and signed or unresolved recovery must not be forced onto a positive log-fold display. This permits a comparison of mechanisms without promising a universal stability result.

## Anticipated user critique — hypothetical

“Why this metric? I want to know whether a larger incumbent resists or suppresses a challenger. A tiny artificial activity cut is a different experiment. Show me what the paired curves actually do, explain why this cut is comparable across models, and show the actual cue-driven incumbent loss beside the coefficient. If the answer is almost flat, show it flat. Do not tune the time until every model tells the same story, or call local deficit shrinkage stability.”

That critique would be fair. The assay earns its place by supplying a transparent, bounded local-response comparison that the old two-model pool intervention could not provide. Its biological relevance still comes from the displayed native cue histories and activity patterns, rather than from the coefficient alone.

## Bounded implementation pre-launch check

The implementation at `polarity1d/manuscript/code/mechanisms/all_five_incumbent_recovery_20261017.py` was checked without integration. `band_to_dense` implements the correct storage identity, and `AuditedBanded` copies both inputs before the native overwrite-capable call and checks the prescribed scaled residual. The independent nonsymmetric four-by-four fixture passes both solve paths, with residuals 1.44 × 10⁻¹⁷ and 2.89 × 10⁻¹⁷. The native-to-dense substitution is enclosed in `try/finally` restoration. The fifteen-case algebraic preflight passes all source-state, kick, controller and pool checks. Independent Spring core weights agree with direct interpolated-edge trapezoidal integration to 1.33 × 10⁻¹⁵ in the checked fixtures. The gate implementation uses the full sampled prefix and the fixed 25-candidate continuous grid, with all three amplitude pairs and a single clock selected across each model's three sizes.

One launch-blocking packaging omission was identified and sent to the author: the frozen worker imports `native.py`, which immediately reads `code/core/protocol.json`, but the initial source-copy list included only core Python files. The author fixed this by including core JSON files, including `protocol.json` and `environment.json`; the correction was checked in the source. The author will run the frozen CLI's help command after initialization to verify import/path closure without integration. No launch blocker remains in this bounded check: the authorized pilot may be frozen and launched. No integrations were run for this check.
