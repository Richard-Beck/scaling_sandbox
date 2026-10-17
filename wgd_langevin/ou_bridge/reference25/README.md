> Active campaign (2026-09-26): [direct WGD fits with persistence-only differences](../../persistence_only/results/report.html). OU-to-polarity fitting and transfer have been abandoned. The material below is historical.

# Fixed-length polarity representation review

All five spatial models use L = 25 µm. Each of four reduced candidates gets one
parameter set jointly fitted to all 18 protocols. Competition is included.
This campaign ends at the representation review; it does not infer length laws,
fit or predict WGD experiments, or estimate noise from deterministic responses.

Open `results/report.html`. It is rebuilt when fits finish, shows generation time,
and has no browser timer. Existing bridge results are not inputs to these fits.
The completed interpretation is in [results/REVIEW.md](results/REVIEW.md).

## Observable and dynamics

The observable is p = 2*x50/L - 1, where x50 is the interpolated active-field
half-mass position, measured from the left boundary. De Belly uses Rac; the
other models use their registered active field (LEGI: response). Zero total mass
maps to p=0 and is flagged. This is a chemical polarity readout, with hypothetical
translation dX/dt = v0*p; De Belly's internal mechanics remain those of the
reflection-symmetric extension with Julia activation defaults. Its reference
domain length is fixed, while its internal mechanical fields remain active.

The four drift equations are the existing Gaussian, linear memory, nonlinear
coupled and self-activation candidates in `../reduced.py`. Their coordinate is
now on twice the previous scale. Bounds on gain, ps, vs, self_a and self_K are
correspondingly doubled. No clipping or tanh observation map is applied. The
report flags predictions outside [-1,1]. The nonlinear candidates are OU-inspired
extensions, not Gaussian OU processes. The latent q is not measured MCA polarity.

## Inputs

Endpoint stimuli use a cosine profile over the outer 10% of the domain for the
four Jilkine models and 20% for De Belly. Endpoint pulses are rectangular in time.
Existing model-specific amplitudes and clocks are retained: peak 0.2, pulse 40 s,
numerical time unit 100 s, conditioning/relaxation horizon 4000 s for the Jilkine
models (20000 s for Otsuji); De Belly uses peak 3.6, pulse 160 s, unit 1000 s,
and horizon 10000 s. These amplitudes are model inputs, not equivalent physical
ligand concentrations. Exact seconds and segment amplitudes are in manifest.json.

- Pulses: strengths 0.5x, 1x, 2x at the standard duration; durations 0.5x and 2x
  at standard strength; one mirrored standard pulse (six trajectories).
- Weak probes: ±5% directional modulation for one numerical unit after tonic
  gradient conditioning; return to the tonic gradient (two trajectories).
- Countercues: first a standard left pulse, followed by cue-free conditioning;
  a right pulse of 0.5x, 2x or 8x duration, plus 2x strength at 2x duration and
  a mirrored standard case. Follow with cue-free relaxation (five trajectories).
- Tracking: tonic gradient conditioning then four cycles with periods 1, 4 and
  16 numerical units (three trajectories).
- Competition: simultaneous endpoint pulses at 1.01x and 0.99x amplitude, and
  their mirror, followed by relaxation (two trajectories).

Tonic gradients run from 0.9*peak at the left to 0.1*peak at the right. Tracking
rotates this gradient sinusoidally while holding its spatial mean fixed.
The reduced input u is the right-half mean minus the left-half mean of the
spatial stimulus, integrated analytically. Competition therefore carries a small
net u despite strong bilateral input; failure here can expose this input reduction
as well as insufficient latent dynamics. A zero p cannot distinguish symmetric
fronts from a uniform profile, so competition profiles are retained for review.

## Fitting and execution

The least-squares objective equally weights the five stimulus groups, then each
trajectory within a group, then each constant/oscillating stimulus segment. Within
a segment saved samples have equal weight; onset sampling is denser. All complete
conditioning histories are included. No trajectory amplitude normalization or
separately fitted v0 is used. Four seeded starts, maximum 240 evaluations each;
convergence status is shown. Errors are descriptive, not statistical acceptance.

Use the existing scaling_sandbox Conda environment. SLURM uses red/small, one CPU
and 3 GB per task; spatial array 0–89, fit array 0–19 after spatial completion.
No array throttle. All outputs stay in this directory.

`check_results.py` checks analytic input projection, the OU step solution,
reflection of polarity, and final reduced trajectories at half integration steps.
Self-activation uses an additional local step limit near its sharp transition at
p=0; the initial five fits before this numerical correction are archived under
`results/archive_before_integration_fix/` and excluded from the review.
