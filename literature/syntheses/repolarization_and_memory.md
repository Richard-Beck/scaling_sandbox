# Reorienting an established polarity, and polarity memory

Serves Fig 1 (cue duration for lasting reorientation), Fig 2 (routes and pulse-train signatures), and Fig 3 (turning regimes). Drafted 2026-10-23.

## Established

- **Reversal mechanism depends on model class.** PDE bifurcation analysis shows that a wave-pinning model and a mutually antagonistic Rac–Rho model reverse an established polarity by distinct routes [@buttenschoen2022repolarization].
- **Input must overlap the domain to move it, in a local-inhibition model.** After a contact that locally raises the off rate, reorientation of a stochastic wave-pinning peak falls as the contact moves from the peak towards the rear, and becomes negligible once it no longer overlaps the domain. An inhibitor that spreads beyond the contact restores a response to rear contacts and drives near-complete reversals [@kryvoruchko2026intrinsic].
- **Cells keep a directional memory in a slow structural layer.** Neutrophil-like cells in 1D channels remember direction through a long-lived cytoskeletal asymmetry, with rear moesin as one long-lived element; perturbing it changes memory [@prentice2016memory].
- **Split versus de novo protrusions differ in their effect on direction.** Split protrusions sustain runs; de novo protrusions interrupt them [@bosgraaf2009pseudopodia].

## Where we extend or depart

- **Size dependence of reorientation.** Fig 1: the shortest pulse giving lasting reorientation lengthens with size in all five models under fixed-width and proportional cues. Fig 2C: replacement time rises 5–12-fold from L19 to L39. Fig 3: half-turn time rises with size at every cue position. The papers checked analyse reorientation at a single size.
- **Three routes, classified across five models.** Slower initiation (Goryachev), a saturation-driven expansion constraint (Wave-pinning, Holmes), and a slow mechanical memory (deBelly, Spring), each with a pulse-train signature (Fig 2D–F). Model-class dependence of reversal is established [@buttenschoen2022repolarization]; this specific classification and its size dependence are ours.
- **Two turning regimes on rings.** Fig 3: a cue overlapping the incumbent displaces it; a separated cue nucleates a competitor that replaces it. The overlap dependence has a single-size precedent [@kryvoruchko2026intrinsic].

## Our interpretation (I)

- The mechanical-memory class (deBelly, Spring) is the model counterpart of the cytoskeletal memory measured by Prentice-Mott et al. Their 1D-channel switching protocol is a template for the pulse-train experiment in Fig 2F.
- Bosgraaf & Van Haastert's split and de novo protrusions map onto our overlapping and separated regimes. That gives the Fig 3 regimes a measurable cellular counterpart.
- Buttenschön & Edelstein-Keshet's two routes may correspond to our displacing and nucleating regimes. Check once read.

## Implications for the manuscript

- **Pedagogical:** that reversal mechanism differs by model class, and that a slow structural layer can store direction.
- **Novel:** size dependence of all three reorientation measures; the three-route classification with its pulse-train signatures.
- **Hypothesis-generating:** pulse trains with biosensor readout at gap ends (Fig 2F); photoactivation straddling or separated from the leading-edge boundary (Fig 3).

## To check

- [@buttenschoen2022repolarization] full text: is domain length varied, and do their routes match our regimes?
- Jilkine & Edelstein-Keshet 2011 (card to be written from `jilkine_ref.md`): gradient-reversal behaviour of the Wave-pinning and Goryachev models.
- Hadjitheodorou et al. 2021 (neutrophil reorientation limited at the rear by myosin II) and Skoge et al. 2014 (cellular memory in chemotaxis). See `queue.md`.
- Optogenetic stimulus tools (Wu 2009 photoactivatable Rac; de Seze 2025 RhoA GEF), cards to be written from the archived notes.
