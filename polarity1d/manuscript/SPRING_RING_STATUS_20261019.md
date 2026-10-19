# Spring ring development status

**Current comparative scope, 19 October:** the author subsequently requested
Spring throughout the Figures 3–4 analyses. The
[revision results](figures/length_mechanisms_20261008/dev/revision_program_20261019/RESULTS.md)
extend the bounded pilot below to three resting circumferences, matched basal
and challenged functional assays, distance–duration responses and replicated
stochastic dynamics. This page preserves the original one-size validation
record; its former sidecar scope is historical. The mean-force closure remains
explicit in the comparative work.

The [scientific reader](figures/length_mechanisms_20261008/dev/spring_ring_validation_20261019/index.html)
shows the declared circular model, an interactive raw-profile comparison of the
three retained challenger locations, and the latest collected numerical
validation snapshot. The reader states its generation time and scheduler
observation. It has no automatic browser refresh.

The existing Spring ring is a **mean-force global-extension analogue**. Its
periodic active/inactive Rac fields couple to one dynamic circumference. The
native local reaction, diffusion, force switch, extension-feedback and dilution
coefficients remain unchanged; replacing two endpoint forces with the spatial
mean of their local force law is an explicit new closure. It supplies no local
elastic deformation field or force-propagation distance. See the
[closure and equations](SPRING_RING_CLOSURE_20261017.md).

The historical closure record documents agent review of a bounded feasibility
pilot. On 19 October, after the native convention was explained, the user
explicitly agreed: **“Yes—carry over the existing stretching, dilution and
diffusion terms.”** The
[native interval adapter](code/core/models/spring.py) verifies that carryover:
both chemical fields already contain the uniform fractional-length dilution
term and material-coordinate diffusion scales with inverse length squared.
The ring preserves these terms while replacing endpoint force sampling with the
mean-force rule. This records the user's approval of the native geometric
carryover without treating it or historical agent review as acceptance of every
modeling choice. Resting circumference
2 is a native-unit choice, without a micron calibration or an established
conversion from interval length to cell perimeter. The full periodic cue has
peak 4 and physical support 0.2; it is not dose-matched to an endpoint half-patch.

The retained preparation at resting circumference 2 produces one polarized Rac
region after a 10-unit cue and qualifies over a 13,000-unit cue-free
continuation. The retained three-position assay uses a 1,000-unit cue, followed
by 13,000 cue-free units. Near activity bridges to the incumbent; more separated
challengers coexist with a stationary old-site maximum before it decays. Every
tested long cue leaves the same basal profile translated to its target. The
[result note](SPRING_RING_DISTANCE_PILOT_20261017.md) and
[independent scientific review](../../current%20manuscript%20figs/reviews/spring_distance_r1.md)
define the observed readouts and their limits.

New numerical checks are complete and independently audited at the bounded
one-size scope: all nine workers pass finite-field, amount, restart-chain and
final stationarity checks. Fresh preparations qualify at 200 and 400 cells;
solver/tolerance, cyclic-rotation and one small conserved-perturbation checks
also pass. Two finite cues at the same offset retain different sites: 0.1 units
returns near the original incumbent, while 10 units retains activity near the
challenger. Neither a minimum duration nor general stability is established.
Their authoritative
manifest, worker records, summary and final audit belong under
`data/runs/spring_ring_validation_20261019/`. A completed worker or terminal
scheduler state alone does not establish a successful validation verdict. The
reader labels the collected state and links both audits. The
[signed geometry accounting](data/runs/spring_ring_validation_20261019/results/geometry_dilution_20261019/summary.json)
separates the concentration contribution of expansion/contraction from Rac
exchange and global extension feedback. A downloadable SLURM/API restart
example is available in the reader.

These experiments concern Rac organization and site competition. They do not
establish a size-dependent migration-persistence effect, a minimum effective
pulse duration, or equivalence to the authors' deforming 2D mechanics. The
canonical five-model interval comparison remains distinct from the ring
adaptation. The earlier three-chemical-model ring checkpoint is preserved;
the current Figures 3–4 revision adds Spring as described above.

Regenerate the reader after collection:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/spring_ring_reader_20261019.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/organize_length_dev.py --group spring_ring_validation_20261019
```
