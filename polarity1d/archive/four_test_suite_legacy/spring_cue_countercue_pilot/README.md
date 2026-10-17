# Spring cue/countercue reference-length pilot

User requested 2026-10-01: quickly test cue-conditioned spring finite countercue,
without preserving spontaneous polarity. This is a bounded pilot, not a full
campaign or a replacement for revision-6 results. Spring competition stays paused.

Two initial states at resting L0=1 and 201 cells: exactly spatially homogeneous,
and the established default cosine perturbation of amplitude 0.001. Neither
receives spontaneous preconditioning or orientation reflection.

Both receive a left-edge Gaussian activation cue, amplitude 4 and native width
parameter 0.2, until full chemical/mechanical settling, then autonomous settling
without cue. Each phase uses the revision-6 500-unit confirmation window,
0.1% pointwise range/projected drift tests and 10000-unit cap. For this pilot,
spring's RHS check uses the exact nonlinear cue term implemented by its solver;
revision-6's cue*v drift approximation is not used. No model coefficients change.

Prepared states must meet revision-6 left polarity eligibility. Eighteen right
Gaussian durations from 0.02 to 5, including values around the revision-6
reference threshold 0.171, restart independently from the prepared checkpoint.
Recovery remains 1500 native time units and success uses revision-6 scoring.
Raw prepared checkpoints, final pulse/recovery profiles, JSON scores and source
hashes are retained in results/. Job array: 20658283, QoS small, one core and 3 GB.
Native length/time are uncalibrated; no generalization across lengths is implied.

## Result

Both initial states settle under the left Gaussian in 900 units, then without
cue in 700 units. Neither is mirrored. Both yield eligible left polarity and
a shortest tested successful right pulse of 0.171, matching revision 6 at L0=1.
The prepared active profile differs from the revision-6 spontaneous checkpoint
by at most 1.9e-8; physical length differs by at most 2e-10. This reference-length
pilot supports cue conditioning as an alternative; other lengths are untested.
