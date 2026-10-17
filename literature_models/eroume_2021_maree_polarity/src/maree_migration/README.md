# Eroumé/Marée deformable-cell migration port

This package embeds the exact nine spatial fields and reaction expressions in
`R/vcell_model.cpp` in a Taniguchi-style diffuse domain. Each chemical is
stored and advanced as `q_i = phi*c_i`:

`dq_i/dt = div(phi*D_i*grad(c_i)) + phi*f_i(c, S)`.

The source chemistry and initial values are transcribed from
`scenarios/wgd_spatial_scale.vcml`. In particular, `J5 = Ip1-delta_P1*P1` is
retained even though the represented uniform PI pool does not appear in that
source. No 1x/2x scaling controls are present.

## Scope and deviations

- The nine Eroumé reactions and diffusion constants are unchanged. The Python
  production integrator uses first-order explicit time stepping rather than
  the R solver's Strang/RK4 reaction plus Crank--Nicolson diffusion.
- The smooth circular starting cell replaces the fixed VCML image geometry.
- The phase field uses the standard quartic double well whose stationary flat
  interface is the initialized `tanh(z/(sqrt(2)*epsilon))` profile; this is a
  normalization choice, not a biochemical change.
- The motility term is explicitly **phenomenological**, not part of the
  published Eroumé model: interface-normal drive is proportional to
  `Rac_act + 0.35*Cdc42_act - Rho_act`. Its interface mean is removed and a
  phase-field area constraint supplies the isotropic pressure. This is the
  requested minimal provisional Marée-inspired coupling; no adhesion,
  traction, nucleus, or explicit actin module is claimed.
- Concentration recovery uses `max(phi, 1e-3)`. Saved concentrations outside
  the cell are numerical extension values and should be masked by `phi`.
- A fixed 96 x 96 micrometre laboratory box is used. Stimulus coordinates are
  never cell-recentered. Runs approaching the box boundary should increase
  `box_size`.
- Linear activation is clipped at zero if a requested gradient would otherwise
  make the Cdc42 activation rate negative.

`free` means uniform stimulation with small, seeded, pairwise mass-balanced
spatial perturbations. `--noise-amplitude 0` is the deterministic symmetric
control. `transient` means free migration *after transient polarization* and
must not be described as spontaneous migration. If the source-derived network
does not break symmetry, that negative result is retained rather than repaired
by changing biochemical parameters.

The reported polarity vector is the first spatial moment of active Cdc42
relative to the phase-field centroid, because the cue enters the Cdc42
activation reaction directly. Rac remains the main movie channel.

## Run

```bash
conda run -n scaling_sandbox python scripts/run_maree_migration.py \
  --mode free --duration 600 --seed 1 --render

conda run -n scaling_sandbox python scripts/run_maree_migration.py \
  --mode rotate --gradient-strength 0.05 --gradient-angle 0 \
  --rotate-at 300 --rotate-angle 90 --duration 600 --render

conda run -n scaling_sandbox python scripts/run_maree_migration.py \
  --mode transient --gradient-angle 0 --duration 600 --render
```

Load the cluster FFmpeg module before rendering if `ffmpeg` is not already on
`PATH`. The importable API consists of `SimulationConfig`,
`StimulusProtocol`, and `run_simulation`. Each output directory contains the
float32 sparse trajectory, a CSV of migration observables, a JSON summary,
main and diagnostic MP4s, a seven-frame montage, and the montage source PNGs.
Saved simulations can be rerendered independently with
`scripts/render_maree_migration.py OUTPUT_DIRECTORY`.
