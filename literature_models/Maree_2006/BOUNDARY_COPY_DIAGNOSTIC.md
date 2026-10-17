# Boundary-copy diagnostic

## Conclusion

The copy sampler, connectivity check, directional pushing lookup, and Metropolis
draw all operate. The standard reconstruction is immobile because the literal
area penalty charges `0.1 nm^-1` for either a one-pixel extension or retraction
when the cell is exactly at its target area. That dominates the typical contact,
actin, and Rho terms and makes the expected copy count effectively zero.

This exposes a parameter-units/discretization problem in the reconstruction,
not a failure to visit boundary sites. With `dx=100 nm`, a pixel has area
`10,000 nm2`, and the implemented expression gives

```text
lambda * (delta_area)^2 = 1e-9 * 10000^2 = 0.1 nm^-1.
```

The specification's stated units are themselves dimensionally inconsistent:
`lambda [nm^-3] * area^2 [nm^4]` does not have the `nm^-1` units of pushing,
temperature, and yield energy. The historical implementation likely used a
different lattice/physical-area conversion or normalization.

## Standard 20-second diagnostic

The native 500×500, 30,000-site model generated 84,352 cell/medium proposals.
None were rejected by the connectivity check. The integrated Metropolis
expectation was only 0.00553 accepted copies, and zero were observed.

During 10–20 s, an average extension proposal had these contributions:

| Contribution | Mean (`nm^-1`) |
|---|---:|
| Area | +0.10000 |
| Contact | +0.00371 |
| Pushing | −0.01914 |
| Rho | +0.01754 |
| Net effective energy | approximately +0.1021 |

Even the most favorable extension observed by 20 s had `dH=+0.0549`. Including
the `Hb=0.046` yield offset at `T=0.008` makes its acceptance probability tiny.

## Positive controls

All controls use the same sampler and seed-independent model logic.

| Control | Expected copies | Observed | Result |
|---|---:|---:|---|
| Target raised by exactly one pixel, 3 MCS | 1.000 | 1 extension | Correct target response |
| Initial outward `P_total=0.2 nm^-1`, 3 MCS | 3.017 | 3 | Correct pushing response |
| `lambda=1e-11`, 10 MCS | 8.51 | 12 | Consistent stochastic response |

The one-pixel target control accepted the energetically favored extension and
then stopped at the new target. The preloaded-P control produced both extensions
and retractions. These tests rule out a dead proposal path, an always-failing
connectivity condition, and a sign error that universally suppresses pushing.

## Area-coefficient sweep

Ten-MCS native-size controls without the polarization pulse gave:

| `lambda` | Expected copies | Observed copies |
|---:|---:|---:|
| `1e-9` | 0.000034 | 0 |
| `5e-10` | 0.0179 | 0 |
| `2e-10` | 3.26 | 2 |
| `1e-10` | 5.56 | 7 |
| `5e-11` | 8.24 | 8 |
| `1e-11` | 8.44 | 11 |

Thus the onset of ordinary area fluctuations occurs around a five- to ten-fold
reduction, well before a hundred-fold reduction. This does not by itself identify
the historically correct conversion or guarantee the reported migration speed.

## Reproduction

Every MCS now records proposal counts, predicted acceptances, accepted copies,
minimum effective energies, and mean area/contact/pushing/Rho contributions in
`metrics.csv`. Diagnostic-only CLI controls are `--target-area-um2`,
`--area-lambda`, `--contact-j`, `--temperature`, `--yield-energy`, and
`--initial-p-total`.

Raw outputs are under `runs/boundary_copy_diagnostics/` in this working tree.
The next step should be an explicitly labeled area-conversion calibration, with
the standard literal interpretation retained as a reference.
