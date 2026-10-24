# Positioning: manuscript claims against the literature

One row per claim in the narrative doc (round 3, outline_v9, live revision 138). Status:
- **established**: shown in prior work; present as pedagogical, with citations.
- **extends**: known mechanism or method; our result adds size dependence or a new setting.
- **novel**: not found in the literature checked.
- **contradicted**: prior work points the other way; address it.
- **untested**: a hypothesis in the manuscript with no direct test yet.

Based on the literature check of 2026-10-23, which was not exhaustive. Everything in the Note column is our interpretation unless it cites a card.

## Section 1. Larger cells need longer cues to reorient (Fig 1)

| Claim | Status | Literature | Note |
|---|---|---|---|
| Required pulse duration rises with length at fixed kinetics in all five models | novel | [@buttenschoen2022repolarization] (reversal at one size), [@hubatsch2019threshold] (size effects at fixed kinetics, small end) | Check Buttenschön et al. 2020 (queue, top) before claiming |
| Fixed kinetics across sizes is a reasonable assumption | established | [@hubatsch2019threshold] | Measured PAR kinetics adapt little to size. Cite in Methods |

## Section 2. Three routes to resisting replacement (Fig 2)

| Claim | Status | Literature | Note |
|---|---|---|---|
| Two scaling modes: broaden at fixed plateau versus raise the peak | extends | [@chiou2018competition], [@chiou2021number] | The saturation rule framed as a size effect |
| Goryachev: larger incumbent depletes the inactive pool, delaying initiation | extends | [@chiou2021number] | Substrate depletion limiting sites is known; its size dependence is ours |
| Wave-pinning, Holmes: saturation forces a challenger to expand | extends | [@chiou2018competition] | Width penalty at equal area and its growth with size are ours |
| deBelly, Spring: slow mechanical memory carries challenger gains across gaps | extends | [@prentice2016memory] | Cytoskeletal memory is established in cells; the model class and size dependence are ours |
| Outside Goryachev, removing the incumbent is the step that slows with size | novel | | |
| Pulse-train signatures distinguish the routes (Fig 2F) | novel | [@prentice2016memory] (stimulus-switching protocol) | Hypothesis-generating |

## Section 3. Larger rings turn more slowly, because of size rather than distance (Fig 3)

| Claim | Status | Literature | Note |
|---|---|---|---|
| Two turning regimes: an overlapping cue displaces the incumbent; a separated cue nucleates a competitor | extends | [@kryvoruchko2026intrinsic] (overlap dependence, one size), [@buttenschoen2022repolarization] (distinct reversal routes), [@bosgraaf2009pseudopodia] (split versus de novo protrusions) | |
| Larger rings reorient more slowly at every cue position | novel | | |
| Beyond setting the regime, separation matters little | extends | [@brauns2021wavelength] | Fast cytosolic exchange makes competition near-global. Check the coarsening law |

## Section 4. Larger cells hold heading longer under fluctuating input (Fig 4)

| Claim | Status | Literature | Note |
|---|---|---|---|
| A diffusing domain loses heading at a rate ∝ D/C², so α = 2 | established | [@kryvoruchko2026intrinsic], [@walther2012stochastic], [@dyer2013wandering] | Pedagogical |
| Drift is a small-signal response, D ∝ f², near-diffusive | established | [@rocco2001front] | Theory of front diffusion; Goryachev departs |
| Physical drift diffusivity depends on size, by model (Wave-pinning down about 2×, Goryachev strongly, Holmes and Spring flat) | novel | | Consider an α = 3 content-averaged null (`syntheses/noise_driven_domain_motion.md`) |
| Edge-pulse response predicts drift and its size scaling | extends | [@rocco2001front], [@weyer2025interface] | Known method (projection onto the translation mode) applied to size |
| The size advantage diminishes at strong input as multiple domains appear | novel | [@chiou2018competition], [@brauns2021wavelength] (coexistence and coarsening background) | |
| Secondary domains start readily in larger rings but rarely take over; Goryachev shows the reverse | novel | [@chiou2018competition] | The saturation rule predicts the per-model direction |
| Heading loss in the multiple-domain regime is not captured by drift plus relocations | novel | [@mikailu2026noise] (noise correlation controls transient patches) | Mechanism not identified |
| Noise-driven relocation is reachable | contradicted (partly) | [@guan2023relocation] | Intrinsic noise may be too weak in a strongly self-reinforcing circuit. Our f is a scan of extrinsic input, so say so |

## Section 5. Published migration data fit greater persistence in enlarged cells (Fig 5)

| Claim | Status | Literature | Note |
|---|---|---|---|
| Enlarged cells are more direct or travel further at long lags | established | Dossiers in `evidence/` | Source papers attribute it to cytoskeleton, vimentin or hypoxia, not size |
| Enlarged cells can be slower yet more persistent | contradicted (by a general rule) | [@maiuri2015ucsp] | Runs against UCSP. This strengthens a speed-independent explanation; state it |
| OU fits infer longer velocity memory | extends | [@mitterwallner2020gle] | OU is a limit of a GLE; the fold is model-conditional |

## Cross-cutting

| Claim | Status | Literature | Note |
|---|---|---|---|
| Simulated heading persistence maps to migration persistence | untested | [@kryvoruchko2026intrinsic] (τ ≈ 1/D), [@vaidziulyte2022persistent] (persistence from protrusion statistics) | Hypothesis |
| Larger size raises persistence | contradicted (in another model class) | [@thuroff2019bridging] (lower), [@ron2025size] (higher) | Scope the claim to polarity biochemistry on a fixed 1D geometry |
| One WGD gives about 1.6-fold slower decorrelation at α ≈ 2, against an OU fold of about 16 | untested | [@prentice2016memory], [@mitterwallner2020gle] | Candidate multipliers listed in `syntheses/persistence_and_migration_theory.md` |
