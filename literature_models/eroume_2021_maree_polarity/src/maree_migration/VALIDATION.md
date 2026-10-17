# Validation snapshot

Validated on 2026-08-23 with the repository `scaling_sandbox` Conda
environment. These are plausibility demonstrations, not fitted results.

## Numerical checks

- The Python reaction rates match independently evaluated values from
  `R/vcell_model.cpp` at the VCML initial state to `1e-13` relative tolerance.
- Frozen-boundary uniform states remain spatially uniform, and the legacy
  `2.6 + 0.05*x` stimulus produces the source model's expected Cdc42 direction.
- Active+inactive Cdc42, Rac, and Rho material changes by at most `9e-16` in
  the moving examples. This directly checks the conservative `q=phi*c`
  update rather than reconstructing material from cutoff movie fields.
- Leakage is defined as material in the deep exterior (`phi < 0.05`), not
  material in the diffuse interface itself.
- Area variation was 0.96% in the constant-gradient example and 1.00% in the
  rotation example. The free/noise start has a one-time 2.85% phase-interface
  relaxation before settling.
- The default 96-micrometre, 1-micrometre grid ran 15,000 steps (60 s) in
  16.5 s after Numba compilation on the validation host. A 600 s run therefore
  extrapolates to roughly three minutes on the same CPU.

The focused Eroumé/migration and existing Taniguchi suites pass together:
`18 passed`.

## Biological outcomes (source-derived chemistry)

| protocol | observation | outcome |
|---|---|---|
| uniform + 1% seeded noise, 30 s | 0.040 micrometre displacement; final active-Cdc42 polarity moment 0.045 micrometre | no convincing spontaneous symmetry breaking |
| uniform deterministic control, 30 s | displacement and polarity moment both below `3e-14` micrometre | symmetric control remained symmetric |
| localized pulse for 20 s, then uniform, 90 s | 3.35 micrometre displacement; polarity moment fell from 0.581 at 45 s to 0.021 at 90 s | polarity and translation were not maintained after cue removal |
| 0.05 linear gradient at 0 degrees, 60 s | 9.25 micrometre displacement in +x; final polarity/cue error 0.01 degrees | directional response succeeded |
| gradient 0 to 90 degrees at 45 s, 90 s total | centroid moved from (56.98,48.00) at switching to (57.18,49.59); final polarity/cue error -11.35 degrees | polarity and migration responded to rotation |

Thus the directional and cue-rotation acceptance cases succeed, but persistent
free migration does not. The failure is retained rather than changing the
published biochemical parameters or silently tuning the phenomenological
boundary coefficients. `transient` output is correctly labelled “free after
transient polarization,” not spontaneous migration.

Example trajectories, CSVs, seven-frame sequences/montages, and main plus
diagnostic MP4s are under `reports/output/maree_migration/examples/`.
