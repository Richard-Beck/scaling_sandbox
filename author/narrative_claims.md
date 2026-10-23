# WGD manuscript: narrative and claims supported by the round-3 figure set

Oct 23, 2026 · @Richard

Across five polarity models with unchanged local kinetics, directionally polarized cells resist reorientation upon challenge by opposing cues and, under fluctuating input, maintain their heading for longer. Model analyses predict that resistance to reorientation arises through one of three routes: slower initiation of a competing domain, saturation that forces a challenger to expand, or a slow mechanical memory. They predict that persistence under fluctuating input arises because polarity domains drift no faster, and in some models more slowly, around a larger perimeter. The analyses also suggest experiments to discriminate these possibilities. Regardless of the mechanism, the enhanced directional persistence predicted by all models is consistent with published migration data. Follow-up experiments in a WGD+/− cell pair are planned.

This maps the claims the round-3 figure set (outline\_v9, commit `87cb5f5` on `codex/round3-outline-compact-20261023`) supports as it stands. Each section title is the claim; bullets carry figure callouts. **(not shown)** marks a step the argument needs that no panel shows. Fig 2 and 3 numbers are computed from the committed CSVs; Fig 1 numbers are read off panel C; Fig 4 and 5 numbers come from their legends. Manuscript text and talking points are not yet vetted.

## 1. Larger cells need longer cues to reorient (Fig. 1)

- We polarized model cells of different lengths, with local kinetics unchanged (Fig. 1A).
- Each cell was then challenged with a finite opposing cue of one of three shapes: fixed width, width proportional to length, or a constant-slope ramp (Fig. 1B).
- For each length, we determined the shortest pulse that produced lasting reorientation **(not shown; duration curves in Fig. S1)**.
- Under fixed-width and proportional-width cues, the required pulse duration increased with length in all five models (Fig. 1C). The rank correlation was perfect in four models and weaker in deBelly (+0.29 to +0.55), whose duration curves are non-monotone (Fig. S1).
- Under the ramp, the length dependence was model-specific (Fig. 1C). Because the ramp peak increases with length (Fig. 1B, right), this comparison confounds size with cue dose.
- The magnitude of the increase in required duration **(not shown)**.

## 2. Larger cells resist replacement by three routes, set by how each model stores the extra active material (Fig. 2)

- Larger cells contain more active material, which the models accommodate in one of two ways (Fig. 2A). Wave-pinning, Holmes, deBelly and Spring broaden the active domain at an approximately fixed peak; Goryachev maintains domain width and increases peak height. These two scaling modes impose different requirements on a challenging domain.
- **Missing step: the mechanistic basis of the two scaling modes (not shown; candidate panel: Wave-pinning vs Goryachev net-activation maps, length-mechanisms report Fig. 3).**
  - In Goryachev, the taller incumbent depletes the inactive pool: opposite-end inactive concentration falls three-fold between L9 and L39. An equal cue therefore activates one third as much material, and the local activation threshold for a challenger rises from 18 to 58.
  - In Wave-pinning, active concentration saturates at a fixed plateau and net activation is confined to the domain boundary, at a level independent of length. A challenger therefore cannot prevail by increasing its height and must instead expand.
- Under continuous opposing cues, replacement time increased 5–12-fold from L19 to L39 in every model, at every cue width for which both lengths were measured (Fig. 2C).
- The Goryachev constraint appears as delayed initiation, approximately 2.5-fold slower at L39 (Fig. 2B). Initiation time was independent of length in Wave-pinning, Holmes and Spring, and apparently in deBelly.
- The expansion constraint appears as a cue-width penalty. At equal cue area, narrow cues were 3.8-fold (Wave-pinning) and 45-fold (Holmes) slower than broad cues at L19, and failed within the simulation horizon at larger lengths. In Goryachev, deBelly and Spring, replacement time varied by at most 1.6-fold across widths (Fig. 2C).
- Pulse trains distinguish three response classes (Fig. 2D, E):
  - **Goryachev:** challenger activity collapses during each gap and must be re-initiated by each pulse.
  - **Wave-pinning and Holmes:** challenger activity is retained across gaps.
  - **deBelly and Spring:** visible challenger activity is lost during gaps, yet subsequent pulses still produce gains and the challenger ultimately takes over.
- Outside Goryachev, the step that slows with size is removal of the incumbent, not initiation of the challenger. Initiation was length-independent (Fig. 2B) while replacement slowed 5–12-fold (Fig. 2C), and challengers in Wave-pinning and Holmes retained their gains across gaps (Fig. 2D, E). This asymmetry, in which challengers start readily but incumbents resist removal, recurs under fluctuating input (Section 4).
- **Missing step: the state that deBelly retains across gaps (not shown; candidate panel: deBelly spring/dashpot graphic, length-mechanisms report Fig. 10).** Membrane strain and cortical flow persist after Rac activity decays, indicating a mechanical rather than chemical memory. Spring, also a mechanical model, presumably behaves similarly.
- These signatures predict distinct biosensor readouts in WGD− and WGD+ cells (Fig. 2F; prospective). Because Fig. 2E shows L39 only, the size-dependent differences depicted in F are **(not shown)**.

## 3. Larger rings turn more slowly toward a new cue, because of size rather than distance (Fig. 3)

- On rings of fixed circumference, a continuous cue was placed at varying distances from the edge of the incumbent domain, and we measured the time required to rotate halfway toward the cue (Fig. 3A, B).
- Two turning regimes emerged depending on whether the cue overlapped the incumbent domain (Fig. 3A). An overlapping cue displaced the incumbent toward it, whereas a separated cue nucleated a second domain that subsequently replaced the incumbent.
- Larger rings reoriented more slowly at every cue position in all four models (Fig. 3B). Per volume doubling, half-turn time increased approximately 1.5-fold in Goryachev, 2–3-fold in Wave-pinning and Spring, and 3–7-fold in Holmes **(fold change not plotted)**.
- Once the cue was separated from the incumbent, further separation had little effect in Wave-pinning, Holmes and Spring: half-turn time varied by less than about 10% across gaps. Distance therefore acted mainly by determining the turning regime. Goryachev was the exception, with half-turn time increasing 2–4-fold with gap.
- The two regimes may correspond to the two ways heading is lost under fluctuating input (Section 4). Input at the incumbent's edge would displace the domain, as an overlapping cue does. Input far from it would have to nucleate a domain that then wins, as a separated cue does **(speculative)**.
- Fig. 3C presents the same data as mean angular speed. It could instead show the small-signal counterpart of Fig. 3B: edge displacement after a weak, brief pulse at the incumbent's edge, at each size **(candidate panel; see Section 4)**.
- Resolution of coexisting domains after cue withdrawal **(not shown)**.

## 4. Larger cells maintain heading longer under fluctuating input, an advantage that diminishes as input strengthens (Fig. 4)

- We summarized the size dependence of persistence by an exponent α: the rate at which heading decorrelates, k, scales with circumference C as k ∝ C^−α, so a larger α means a stronger advantage for larger cells. α is estimated locally from three adjacent circumferences, so each size has its own estimate; where a single power law holds, these coincide. Simple null models of polarity-domain wandering predict distinct values of α (Fig. 4A):
  - Drift at a fixed physical diffusivity gives α = 2.
  - Relocation at a size-independent rate gives α = 0.
  - Relocation at a rate proportional to perimeter gives α = −1.
- Combining drift with relocation defines two noise-dependent regimes (Fig. 4B). Drift dominates at weak input, giving α ≈ 2 (the drift regime). Relocations increase more steeply with input and dominate at strong input, driving α toward 0 or −1 (the relocation regime). The transition occurs at lower input in larger rings, because drift costs a larger ring less angle and relocations therefore dominate sooner. The input dependences are illustrative.
- The regimes also set how noise tolerance scales with size. We define noise tolerance as the input at which heading correlation over a fixed reference lag falls to 0.5. In the drift regime it increases in proportion to circumference, three-fold from C20 to C60; relocations reduce this fold (Fig. 4B, right).
- Wave-pinning, Holmes and Spring show the drift regime at weak input (Fig. 4C). α was approximately 2 in Holmes and Spring, and 2.5–3 in Wave-pinning, where larger cells therefore gain more persistence than drift geometry alone predicts. As input strengthened, α declined toward 1 in all three, as expected on entering a second regime, although α did not reach 0 within the tested range. In Goryachev, heading broke down too abruptly for local estimates.
- Noise-tolerance folds agree with this picture (Fig. 4C, right). Holmes (2.2) and Wave-pinning (2.4) fall below the drift-only fold of 3, indicating that the second regime limits tolerance in larger cells. Spring (3.1) matches drift alone. Goryachev (4.1) exceeds every null model.
- At inputs near the transition, the larger ring maintained its heading markedly longer in all four models (Fig. 4D).
- The drift coefficient was measured directly by tracking the centroid of the dominant active domain while it remained dominant and fitting its mean squared displacement (Fig. 4E, row 1; normalized to the smallest ring at the same input). Size reduced drift in Wave-pinning (about two-fold from C20 to C60) and, far more strongly, in Goryachev, but not in Holmes or Spring. This is consistent with Fig. 4C: α above 2 in Wave-pinning, approximately 2 in Holmes and Spring, and a Goryachev tolerance fold above the drift-only expectation. The Holmes and Spring results were surprising given Fig. 3, where larger rings slid toward a continuous cue 1.3–3.5-fold more slowly per volume doubling (computed from Fig. 3 data; not plotted). Cue-driven sliding therefore does not predict noise-driven drift in these models. One reconciliation is that Fig. 3 drives the whole domain with a strong, sustained cue, whereas drift may reflect small, transient responses at the incumbent's two edges **(speculative)**.
- Drift behaved as a small-signal response. For f ≤ 0.01, drift diffusivity scaled as f^1.94–1.99 at every size in Wave-pinning, Holmes and Spring, as expected for a linear response to input whose variance scales with f²; the domain's MSD was close to diffusive (exponent about 1.1) **(computed from Fig. 4 data; not plotted)**. Goryachev did not scale this way (f^1.7–3.2, MSD exponent about 0.75), consistent with a pinned spike. If the incumbent's two edges respond independently to local input, drift is set by how far an edge moves per unit of input at that edge. The edge-pulse response proposed for Fig. 3C would then predict drift and its size dependence without noise simulations **(speculative; candidate panel: predicted versus measured drift)**.
- The strongly suppressed drift in larger Goryachev rings may reflect the same inactive-pool depletion that delays initiation (Fig. 2B) **(hypothesis)**. The slower takeover by separated domains in larger rings (Fig. 3) is consistent with relocation rates that were similar or lower in larger rings (Fig. 4E, row 2).
- In the drift regime, drift plus counted relocations predicted the observed rate of heading decorrelation to within 15% (Fig. 4E, row 3). Counted relocations were absent below f ≈ 0.01, and the regime was nearly pure: multiple domains occupied at most 0.1% of frames in Wave-pinning for f ≤ 0.01, and up to 7% (Holmes) and 13% (Spring) only at f = 0.01 (computed from Fig. 4 data; not plotted). The detector counts only regions above 10% of basal contrast, so sub-threshold satellites, visible near input hot spots, are not counted. Because drift alone accounts for the observed decay here, such satellites appear not to affect heading (candidate supplement: detector-threshold sweep on the exemplar profiles).
- The second regime in the models is not the relocation regime of the nulls. The decline in α coincided with the appearance of multiple active domains, which were more frequent in larger rings (Fig. 4E, row 4). Counted relocation rates were similar or lower in larger rings (Fig. 4E, row 2). Drift diffusivity rose with size at strong input, up to about four-fold in Holmes and two-fold in Spring. Even so, drift plus relocations under-predicted the observed decorrelation where multiple domains predominated: 0 of 33 such conditions were within 15% (Fig. 4E, row 3). Heading is therefore lost through a process the tracked measures do not capture **(not identified)**.
- Larger rings spawned secondary domains far more readily than they lost their incumbents. In Holmes at f = 0.03, the multiple-front fraction rose from 7% at C20 to 69% at C60, while counted takeovers only doubled **(computed from Fig. 4 data; not plotted)**. This is the asymmetry of Fig. 2: in Wave-pinning, Holmes, deBelly and plausibly Spring, a challenger starts readily at any size, but the incumbent becomes harder to remove. Goryachev, in which initiation itself slows with size, showed the reverse trend, with multiple domains becoming rarer in larger rings. If secondary domains persist without winning, heading would be lost through shifts in the balance of activity between coexisting domains rather than through relocation, which the dominant-domain tracking cannot register **(speculative)**.
  - Candidate panel: an exact decomposition of heading change into domain motion, balance between coexisting domains, and domain birth or loss, computable from the stored per-region centroids and amounts.
  - Candidate panel: birth rate per unit perimeter, birth distance from the incumbent's edge, and takeover probability by distance and size. This is the stochastic counterpart of the two turning regimes in Fig. 3.
- The detector-free mixture fit disagrees with the counted-event estimates, over-predicting decorrelation by orders of magnitude at weak input (Fig. 4E, rows 2–3, dashed).
- These quantities suggest observables for the migration experiments **(speculative)**:
  - diffusion of the lateral edges of the leading active or protrusive region;
  - birth rate, position and fate of secondary protrusions;
  - heading autocorrelation.

  If secondary protrusions are rare and edge fluctuations are size-independent, the models predict persistence rising roughly with the square of size, about 1.6-fold per WGD. Frequent secondary protrusions that rarely take over would erode that advantage. Such measurements could place an experiment in a regime, and point to conditions that maximize or remove a WGD advantage. Photoactivation straddling or separated from the leading-edge boundary (Fig. 3) would probe the two regimes directly **(candidate panel: schematic predictions in the style of Fig. 2F)**.

## 5. Published migration data fit greater persistence in enlarged cells (Fig. 5)

- Enlarged or polyploid groups were not consistently faster: speed was lower in Xuan and Mallin and higher in Jemaà. They did, however, show greater directness (Mallin 3.1-fold, Hosny 1.7-fold) or greater long-lag displacement (Xuan 2.1-fold) (Fig. 5A).
- In Xuan, the mean squared displacement of the enlarged group was initially lower but exceeded that of the reference group beyond approximately 2.3 h, the profile expected for slower but more persistent motion (Fig. 5B).
- Ornstein–Uhlenbeck fits inferred longer velocity memory in all five comparisons (pooled fold approximately 16; Fig. 5E), robust to the choice of prior and pooling (Fig. 5F).
- The fit does not reproduce the shape or crossover of the Xuan curves (Fig. 5B), so the inferred fold is conditional on the model.

## 6. New WGD motility experiments (Fig. 6, reserved)

Reserved for the SUM-159 WGD+/− motility assays; no data are yet available.

## What the set doesn't yet support

- **No figure connects imposed replacement (Figs 1–3) to noise-driven persistence (Fig. 4).** The link must be made in the text or shown directly. The revised Fig. 4E strengthens one candidate: where multiple active domains predominate, heading is lost faster than drift and relocation of the dominant domain predict, which suggests that competition among coexisting domains, the process studied in Figs 2–3, drives the decline in α. The excess of secondary-domain births over takeovers in larger rings, which mirrors the readily started but rarely successful challengers of Fig. 2, would be the most direct bridge **(hypothesis; not shown)**.
- **Magnitude (not shown).** At α ≈ 2, a single WGD (a 1.26-fold increase in circumference) would reduce the rate of heading decorrelation only about 1.6-fold, whereas the published OU persistence fold is approximately 16.
- **The connection between simulated heading persistence and migration persistence remains a hypothesis.**
- **The figures do not support** claims that distance matters beyond determining the turning regime, that multiple domains cause the decline in α, or any specific mechanism of heading loss in the multiple-domain regime.

* **Figs 3–4 lack experimental signatures.** Unlike Fig. 2F, neither figure proposes a measurement that would discriminate the mechanisms it describes. The summary's claim that the analyses suggest discriminating experiments currently rests on Fig. 2 alone. Candidate signatures are proposed in Sections 3–4: the edge-pulse response, and the statistics of secondary-domain births and takeovers.

## The thread in one paragraph

Across five models, larger cells resist reorientation (Figs 1–3). Under fluctuating input they maintain their heading for longer, by approximately the amount expected from drift alone at weak input, and this advantage diminishes as input strengthens and multiple domains appear, through a process not captured by drift or counted relocations, plausibly shifting balance among secondary domains that, as under imposed cues, start readily but rarely displace the incumbent (Fig. 4). Published migration data show the increase in persistence this would predict (Fig. 5).
