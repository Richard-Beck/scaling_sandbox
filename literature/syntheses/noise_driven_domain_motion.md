# Noise-driven motion of a polarity domain

Serves Fig 4 (drift, relocation, multiple domains) and the proposed Fig 3C–D (edge-pulse response, drift prediction). Drafted 2026-10-23.

## Established

- **A single domain under noise performs a random walk, and heading decorrelates by rotational diffusion.** In stochastic wave-pinning the pinning position fluctuates more as molecule number falls, until the domain collapses [@walther2012stochastic]. On a ring, a wave-pinning peak's angular diffusion falls roughly as 1/N with copy number, and 1/D gives persistence times of hours [@kryvoruchko2026intrinsic]. Their heading is the vector sum of membrane-bound positions, the same construction as our P = M1/M0.
- **Low copy number gives multiple transient peaks.** Near 2×10³ molecules no stable peak survives and several transient peaks appear [@kryvoruchko2026intrinsic].
- **Domain wandering is seen in cells and can be useful.** The yeast polarity patch wanders in sub-saturating pheromone; actin-directed traffic contributes, receptors bias it, and blocking it impairs gradient tracking [@dyer2013wandering].
- **Whether noise relocates a site depends on the circuit.** Molecular noise moves clusters with few limiting molecules or short dwell times [@ramirez2021stochastic]. In the core yeast circuit it does not overcome positive feedback; relocation needed a second feedback loop [@guan2023relocation].
- **Front diffusion is a projection of noise onto the translation mode.** Bistable and pushed fronts wander diffusively. The diffusion coefficient follows from splitting motion into a position variable plus fluctuations orthogonal to the translation mode, and at small noise scales with noise variance [@rocco2001front].
- **Noise model matters.** External and internal noise are not interchangeable, and temporal correlation controls whether transient secondary patches form [@mikailu2026noise]. In wave-pinning, polarity loss under extrinsic noise depends on its spatial and temporal correlation [@defranciscis2012bounded].
- **Extrinsic fluctuating fields are a standard way to set persistence** in reaction-diffusion cell models [@camley2014prm].

## Where we extend or depart

- **Size as the variable.** None of the above varies cell size at fixed kinetics and fixed concentration while measuring drift. Kryvoruchko & Camley vary copy number at fixed size [@kryvoruchko2026intrinsic]. Walther et al. may have varied width with copy number; not yet confirmed [@walther2012stochastic].
- **Size-dependent drift diffusivity.** Fig 4E: physical drift diffusivity D_s falls about two-fold from C20 to C60 in Wave-pinning, falls much more in Goryachev, and is flat in Holmes and Spring. Nothing in the papers checked reports this.
- **Drift as a small-signal response across models.** Fig 4 drift scales as f^1.94–1.99 with a near-diffusive MSD in Wave-pinning, Holmes and Spring, as front-diffusion theory predicts [@rocco2001front]. Goryachev departs (f^1.7–3.2, MSD exponent about 0.75), consistent with a pinned spike rather than a front pair.
- **The multiple-domain regime as a separate channel of heading loss.** Covered in `competition_coexistence.md`.

## Our interpretation (I)

- **Two candidate drift nulls, not one.**
  - *Edge-limited drift.* If heading moves only when the two edges move, and each edge responds to local input with a size-independent coefficient, D_s is size-independent and α = 2. This is the current Fig 4 drift null.
  - *Content-averaged drift.* If the domain's position averages independent fluctuations over its content, as in the independent-molecule argument of [@kryvoruchko2026intrinsic], then D_s ∝ 1/content. Wave-pinning, Holmes and Spring broaden their domain with size (Fig 2A; deBelly also does, but is not in the ring figures), so content ∝ C and α = 3.
  - Wave-pinning's α of 2.5–3 sits between the two. Holmes and Spring at α ≈ 2 look edge-limited. A content-averaged null would let Fig 4 say which models are edge-limited, which is a sharper statement than "above or at geometry".
- **Intrinsic noise would add size dependence.** At fixed concentration a WGD cell has about twice the molecules. If intrinsic noise contributes, noise per unit length is unchanged but content grows, pushing broadening models towards the content-averaged null. Our input f is extrinsic, so this is a caveat for the cell comparison, not a model result.
- **The edge-pulse → drift prediction is a known method applied to a new question.** The edge-pulse response is the response function onto the translation mode [@rocco2001front]. Predicting D_s and its size dependence from it, without noise simulations, would be new. Present it that way.
- **Interface theory may supply the edge law.** Effective interfacial tension in mass-conserved patterns [@weyer2025interface] is a candidate description of each edge. Departures from α = 2 would then have to come from coupling between the edges through the shared cytosolic pool, which is exactly where the scaling modes of Fig 2A differ.

## Implications for the manuscript

- **Pedagogical:** drift picture, geometric α = 2, f² scaling. Cite as established; keep short.
- **Novel:** measured D_s(C) by model; α above 2 in Wave-pinning; the edge-pulse prediction of drift and its size scaling.
- **Hypothesis-generating:** edge diffusion of the leading active region and its size dependence, measured in SUM-159 (`measurement_methods.md`).
- Cite [@kryvoruchko2026intrinsic] as the nearest neighbour. Note in Methods that our noise is extrinsic, and state its correlation structure.

## To check

- Walther et al. full text: was width varied with N, and is pinning-position variance reported against width?
- Armero et al. (original front-diffusion result), Meerson et al. 2011 (front diffusion from shot noise), Altschuler et al. 2008, Hecht et al. 2010. See `queue.md`.
- Correlation time and length of our Fig 4 input, and whether satellite frequency depends on them.
