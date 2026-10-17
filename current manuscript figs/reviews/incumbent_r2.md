# Fresh-context review of main-text candidate 1

Reviewed 17 October 2026. This review reads the six-panel PNG, candidate discussion, technical caption, measurements, native reaction definitions and retained removal/return archives. It does not edit the candidate or perform new integrations. The numerical decomposition is in [incumbent_r2_return_state_decomposition.csv](incumbent_r2_return_state_decomposition.csv).

**Decision:** the measured immediate net recovery is real, but the figure needs a focused conceptual restructure before it is a convincing main-text candidate. The supported mechanism is stronger incumbent **capture of available inactive substrate**, together with reduced remote substrate availability in Goryachev. The present result does not support describing the taller incumbent as having a stronger direct restoring response to active-field loss at unchanged inactive concentration. The retained data already distinguish these interpretations; another simulation campaign is unnecessary for this revision.

## The decisive returned-state result

Write the native basal transfer reaction as

$$
f(u,v)=A(u)v-cu,
\qquad
A_{\mathrm G}(u)=\frac{0.33u^2+0.67u}{101},
\qquad
A_{\mathrm{WP}}(u)=0.067+\frac{u^2}{1+u^2}.
\tag{1}
$$

Here \(c=0.01733\) for Goryachev and \(c=1\) for Wave-pinning. Let \(p\) and \(s\) denote the actual perturbed and matched-sham states immediately after exact inactive-profile return. The following is an exact, ordered finite-difference identity:

$$
\begin{aligned}
f(u^p,v^p)-f(u^s,v^s)
={}&\underbrace{\left[A(u^p)-A(u^s)\right]v^s-c\left(u^p-u^s\right)}_{\text{active-state contribution at sham inactive field}}\\
&+\underbrace{A(u^p)\left(v^p-v^s\right)}_{\text{inactive-state contribution at perturbed active field}}.
\end{aligned}
\tag{2}
$$

Each term is integrated over the unchanged original half-contrast core with the native finite-volume weights and divided by its actual active deficit at return. This ordering matters: the active change is made first at sham inactive concentration, then the inactive change is made at perturbed active concentration. Reversing the order reallocates the nonlinear cross term. These are fixed-state algebraic contributions, not separately evolved interventions or an attribution of reversal time.

| Model | Length | Active-state contribution / damage | Inactive-state contribution / damage | Their sum: chemical rate / damage | Diffusion / damage | Actual net rate / damage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Goryachev | 9 | −0.02438 | 0.99747 | 0.97309 | 0.00456 | 0.97766 |
| Goryachev | 21 | −0.02504 | 5.41906 | 5.39402 | 0.00463 | 5.39865 |
| Goryachev | 39 | −0.02346 | 15.30885 | 15.28539 | 0.00410 | 15.28948 |
| Wave-pinning | 9 | 0.09699 | 0.52035 | 0.61734 | 0.01728 | 0.63462 |
| Wave-pinning | 21 | 0.16566 | 0.58304 | 0.74870 | 0.00699 | 0.75569 |
| Wave-pinning | 39 | 0.18724 | 0.60295 | 0.79019 | 0.00367 | 0.79387 |

The CSV contains all 36 model–length states with full precision. The split reconstructs the archived chemical rate at every state. Across all 18 Goryachev lengths the active-state contribution is negative, between −0.02516 and −0.02346 per native time unit. Including the actual active-diffusion difference still leaves this fixed-sham-inactive hybrid negative, between −0.02047 and −0.01936. Thus Goryachev's positive chemical repair in this experiment comes entirely from the elevated inactive field after return. Its direct active-state contribution opposes recovery slightly. Wave-pinning has positive contributions from both changed fields.

This is more informative than the near coincidence of reaction-only and net rates in the existing final row. That coincidence establishes that diffusion does not erase the effect; it does not identify which changed chemical field produces it.

## Genuine recovery versus trivial external replacement

The active field is continuous at return. It subsequently gains material relative to its sham, and the full native reaction-plus-diffusion budget gives a positive deficit-decay rate. Calling that measured response immediate core recovery is correct.

The return does not add an extra bolus beyond what was withheld. After return the global inactive surplus relative to sham equals the global active loss, within the retained mass balance. At Goryachev lengths 9, 21 and 39, both are respectively 1.159, 1.165 and 1.213 times the core active deficit. The pulse left more material inactive because less was activated during the hold; exact return exposes that residual inactive surplus. The strong response is native conversion of that redistributed material, rather than recovery by an instantaneous active-field reset.

Nevertheless, the sign of rebound after restoring substrate is not by itself an interesting mechanism. The informative result is the length dependence of substrate capture. Goryachev concentrates active material into a much taller pole; its unsaturated \(A_{\mathrm G}(u)\) grows strongly with local active concentration. Wave-pinning mainly broadens a pole with approximately unchanged height, and its saturating \(A_{\mathrm{WP}}(u)\) changes much less. Equation (2) exposes this distinction directly. A useful physical statement is that the enlarged Goryachev incumbent can recapture available inactive material more strongly, while leaving less inactive material available for initial activation at the opposite end.

This formulation connects the two observed effects through a shared pool. It does not establish that recapture defeats a simultaneously growing rival; the global removal experiment never creates that rival.

## Is normalization by tiny unequal damage meaningful?

For this intervention the ratio has an exact operational meaning: the instantaneous logarithmic decay rate of the sham-subtracted core deficit. Equal active damage is not required for that identity. The available amplitude controls also support the small-perturbation interpretation: over 0.1%, 0.5% and 1% removal at the primary hold time, the largest relative range in net gain is 0.0043% across Goryachev lengths and 0.1225% across Wave-pinning lengths. The retained tighter-tolerance endpoint reruns address numerical accuracy at even smaller damage.

The large Goryachev ratio is not created by dividing by a decreasing absolute deficit. From length 9 to 39, actual core damage grows from 0.0001590 to 0.0003717, whereas the excess net rate grows from 0.0001555 to 0.0056833. Fractional damage decreases because the original incumbent contains more active material. Both the raw fractional traces and the normalized rates therefore carry information.

The restriction is response direction. Equal fractional removal does not create the same active and inactive perturbation shape at every length. The ratio compares native responses to the specified global disturbance, not a state-independent restoring coefficient. Tiny fractional losses are suitable for this local diagnostic; they do not establish tolerance of a biologically consequential finite challenge.

## Focused revision that would make this a useful candidate

Retain the two-model profile comparison and the actual early damage/recovery traces. Replace the current reaction-only versus net overplot with the ordered active-state and inactive-state contributions in equation (2), alongside the actual net rate. Keep the diffusion agreement in the technical evidence. Put the native activation laws near the profile comparison, so that taller versus broader poles have a concrete kinetic meaning.

Lead with the coherent message: enlargement reorganizes a conserved pool differently in these two models; in Goryachev, the taller incumbent has stronger substrate capture and lower remote inactive concentration. Treat “stronger maintenance” as the hypothesis that motivated the test, not as a synonym for its entire measured effect. The opposing-cue identity is useful only if it is linked to that same pool competition logic; as a detached annotation it risks looking like inactive concentration plotted twice.

With that restructure, this is a worthwhile single main-text candidate to discuss with the user. The current six-panel form is not ready to carry a broad incumbent-strength claim. No additional native integration is needed to make the narrower chemical mechanism concrete and reviewable. An independently evolved perturbation would become necessary only if the intended claim were direct active-field stability or successful competition against a rival.

## Anticipated user critique

The following is constructed wording in the user's plausible direct style, not an actual quotation:

- You took substrate away and put it back. I expect activation to resume. What have you shown about a larger incumbent, beyond that obvious response?
- The large number is almost entirely the response to extra inactive material left after the pulse. Call that substrate capture and explain why it grows with length. Do not let “repair” quietly become intrinsic stability.
- Why do I need the last row to tell me that diffusion makes almost no difference? Use that space to show what the taller pole actually does differently from the broad pole.
- I can accept this as a possible way enlargement favors the old pole, if you show the shared-pool logic clearly. It still does not tell me whether a challenger forms, persists or wins.

These objections identify one decisive issue—what chemical change causes the rebound—rather than a need for an expansive new campaign. Equation (2) and the retained state data answer that issue now.

## Evidence checked

The candidate and technical caption are under `figures/length_mechanisms_20261008/dev/main_candidate_1_20261017/`; the decomposition reads its `supporting/data/length_measurements.csv` and the exact corresponding `returnedstate` arrays in `data/runs/inactive_pool_restoring_response_20261008/`, together with each `epsilon_0.npz` sham. It uses the native rates in `code/core/models/goryachev.py` and `code/core/models/wave_pinning.py`. Primary protocol is 0.5% whole inactive-profile removal, 0.1 native-unit hold and exact pointwise return. Amplitude checks use `data/runs/inactive_pool_figure_exports_20261008/response_vs_length.csv`; accuracy checks use the original run's `results/accuracy_check.json`. All paths in this paragraph are relative to `polarity1d/manuscript/`.
