# Taniguchi PIP2/PIP3 excitable-wave port

This package implements the equations and parameter tables transcribed in
`docs/jobs/02_taniguchi.txt`; no article or supplementary-information lookup
was used. `fixed.py` implements SI Eq. S6 on the stated equilateral triangular
lattice. `phase_field.py` implements SI Eqs. S8--S9 on the stated regular
square lattice and includes every Fig. 4F--I parameter set from Table S1.

## Reconstructed numerical choices

- Stochastic synthesis is an instantaneous, exactly balanced `U -> V`
  impulse. The source specifies the distribution but not whether an impulse
  or one-step source was used.
- The initial phase field is the tanh profile proposed in the job description.
- The default computational box is `[-10,10]^2`, placing its boundary in the
  negligible tail of the radius-5 diffuse cell; the source does not specify
  its box extent.
- A stochastic centre in the deformable model is sampled uniformly from the
  sharp intracellular set `phi >= 0.5`.
- Concentration recovery from conservative variables uses
  `max(phi, 1e-3)`, and conservative chemical fields are zeroed below that
  cutoff. Concentration is undefined as `phi -> 0`; this safeguard prevents
  the exponentially small exterior tail from feeding roundoff into the flux.
- The fixed-circle boundary is the masked triangular-lattice approximation to
  a circle. A missing neighbour is reflected to the centre value, giving the
  requested zero-flux condition.

The deformable solver integrates the conservative variables `phi*U` and
`phi*V` together with `phi` using explicit midpoint RK2. The phase field is
bounded to `[0, 1]` after each complete step to remove floating-point
overshoot.

## Run

Use the project environment:

```bash
conda run -n scaling_sandbox python scripts/validate_taniguchi.py
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=fixed --condition=5.7 --end-time=30
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=phase --condition=4F --end-time=20
```

The full four-condition runs have SLURM entry points in
`workflow/run_taniguchi_fixed.slurm` and
`workflow/run_taniguchi_phase.slurm`. Neither the API nor the command-line
driver contains a cell-size/sandbox experiment.

The job text mentions additional Movie S3 parameterizations but does not list
their numerical values. They are optional in the stated acceptance criteria
and are therefore not implemented; adding them would require the missing
values from the user rather than a literature search.
