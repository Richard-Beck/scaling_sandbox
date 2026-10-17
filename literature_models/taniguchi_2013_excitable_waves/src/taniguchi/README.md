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

The production `fast` backend uses the same equations, stencil, grid, and RK2
update as the original `reference` backend, but fuses scalar work and reuses
15 scratch arrays instead of allocating full-grid temporaries twice per
step. The command-line production timestep is `5e-4`; the former `8e-5`
timestep remains available with `--dt=8e-5`, and the original kernel remains
available with `--backend=reference`. In a controlled ten-time-unit excited
cell comparison (including a more aggressive `1e-3` test), the fast timestep retained a `1.0` outline intersection over
union at every sampled time, with phase-field RMS error below `8.2e-6` versus
the `8e-5` result.

## Run

Use the project environment:

```bash
conda run -n scaling_sandbox python scripts/validate_taniguchi.py
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=fixed --condition=5.7 --end-time=30
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=phase --condition=4F --end-time=20
```

For a strict numerical-reference run:

```bash
conda run -n scaling_sandbox python scripts/run_taniguchi.py \
  --model=phase --condition=4F --end-time=20 \
  --backend=reference --dt=8e-5
```

The full four-condition runs have SLURM entry points in
`workflow/run_taniguchi_fixed.slurm` and
`workflow/run_taniguchi_phase.slurm`. The native reproduction driver does not
silently add a size perturbation. The separate WGD v1 adapter uses the same
optimized phase-field API under the scaling rules in
`experiments/wgd_v1/manifests/taniguchi_phase_field.md`.

## Microscopy-style time lapses

The phase-field trajectories can be rendered without rerunning the model. The
default renderer produces both a PIP3 false-colour movie and a cyan-PIP2 /
magenta-PIP3 composite. Both use the moving `phi=0.5` cell outline, fixed
movie-wide intensity limits, a time label, and a scale bar. Five
activity-aware full-resolution PNG frames (endpoints plus separated PIP3
peaks) accompany each MP4.

```bash
module load FFmpeg/7.1.2-GCCcore-14.3.0
conda run -n scaling_sandbox python scripts/render_taniguchi.py \
  reports/output/taniguchi/phase_{4F,4G,4H,4I}.npz --mode=both
```

For an HPC production render, submit `workflow/render_taniguchi_phase.slurm`
after the four simulation files exist. Render metadata records the intensity
limits and movie settings beside each condition's outputs in
`reports/output/taniguchi/rendered/`.

The simulation workflow accepts `TANIGUCHI_SEED`, `TANIGUCHI_OUTPUT`,
`TANIGUCHI_DT`, `TANIGUCHI_BACKEND`, and `TANIGUCHI_OUTPUT_INTERVAL`; the
render workflow accepts `TANIGUCHI_INPUT` and `TANIGUCHI_RENDER_OUTPUT`. These
environment overrides allow reproducible stochastic realizations to coexist
without overwriting the default seed-1 trajectories.

The job text mentions additional Movie S3 parameterizations but does not list
their numerical values. They are optional in the stated acceptance criteria
and are therefore not implemented; adding them would require the missing
values from the user rather than a literature search.
