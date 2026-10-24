# From polarity to migration persistence

Serves Fig 4 (how simulated heading maps to migration), Fig 5 (OU inference from published data), and the Discussion. Drafted 2026-10-23.

The manuscript keeps three readouts distinct: observed migration summaries, conditional OU velocity-memory estimates, and simulated heading persistence. Their connection is the hypothesis. This synthesis collects what the literature says about each link.

## Established

- **Polarity diffusion sets persistence time, in models.** A wave-pinning peak's angular diffusion D_P gives a persistence time of about 1/D_P, of order an hour at plausible copy numbers [@kryvoruchko2026intrinsic]. Fluctuating inhibitory fields are a standard way to set persistence in reaction-diffusion cell models [@camley2014prm].
- **Persistence can be predicted from protrusion statistics.** In RPE1 and HeLa cells, a minimal model reproduces persistence time from the number of simultaneous protrusions and the alignment of the polarity axis with motion [@vaidziulyte2022persistent].
- **Faster cells are more persistent (UCSP).** Across many systems persistence is coupled to speed, through advection of polarity cues by actin flows [@maiuri2015ucsp].
- **Persistence need not reflect a slow molecular memory.** A generalized Langevin fit to cancer-cell tracks gives a memory kernel of about 6 min yet an emergent persistence time near 164 min, from a near-cancellation of instantaneous and delayed friction. The OU model is a limit of this description [@mitterwallner2020gle].
- **Size can raise or lower persistence through shape and mechanics.** Larger contact area eases polarization and raises persistence and speed [@ron2025size]; in a mechanical cell model larger cells are less persistent [@thuroff2019bridging].
- **Slow structural memory adds persistence.** A long-lived cytoskeletal asymmetry gives directional memory beyond signalling timescales [@prentice2016memory].

## Where we extend or depart

- **A polarity-level reason for size-dependent persistence, independent of speed.** Our models give longer heading persistence in larger cells without any change in speed or actin flow.
- **The published enlarged cells do not follow UCSP.** In Xuan 2018 and Mallin 2023 the enlarged groups are slower yet more persistent (Fig 5A). Under UCSP slower cells should be less persistent [@maiuri2015ucsp]. Jemaà shows the enlarged group faster, so there UCSP could explain a persistence gain on its own.

## Our interpretation (I)

- **The departure from UCSP strengthens the Fig 5 argument.** Where speed and persistence move apart with size, a speed-coupled mechanism cannot be the whole story. State this in the Fig 5 text as our inference, and treat Jemaà as the weakest test.
- **The magnitude gap.** At α ≈ 2, one WGD (about 1.26-fold in circumference) slows heading decorrelation about 1.6-fold, while the pooled OU persistence fold is about 16 (narrative, "What the set doesn't yet support"). Candidate multipliers from the literature: slow structural memory [@prentice2016memory]; intrinsic-noise averaging over a larger domain (α up to 3, `noise_driven_domain_motion.md`); and the sensitivity of emergent persistence to friction balance [@mitterwallner2020gle]. None is shown; list them as open in the Discussion.
- **The OU fold is model-conditional.** The OU fit does not reproduce the Xuan crossover (Fig 5B). A generalized Langevin fit [@mitterwallner2020gle] is a natural supplement check.
- **Protrusion statistics bridge the readouts.** The indices of [@vaidziulyte2022persistent] measure in cells the quantities our Fig 4 decomposition tracks in models (number of active regions; heading–domain alignment). Using them would connect simulated heading persistence to observed migration without passing through OU.

## Implications for the manuscript

- **Pedagogical:** τ ≈ 1/D for a diffusing polarity axis; UCSP as the default explanation of persistence differences.
- **Novel:** a speed-independent, size-based contribution to persistence, supported by enlarged cells that are slower but more persistent.
- **Hypothesis-generating:** measure unicity, alignment and edge diffusion in SUM-159 WGD+/− (`measurement_methods.md`).

## To check

- Whether UCSP is claimed within a population or only across conditions.
- Persistence statements in the full text of [@ron2025size] and [@thuroff2019bridging].
- Full author list for [@mitterwallner2020gle].
