# Native maintenance at the observed states

This is a direct evaluation of the frozen Goryachev equation on the already saved full fields, with the cue set to zero after withdrawal. It adds no integration or new perturbation. All rates are native activity amount per native time; activity concentration and amount remain distinct.

The native transfer is

$$J=(a e_c u^2+b e_c u+S)v-cu.\tag{1}$$

Recruitment below means endogenous recruitment only (the first two terms); stimulation is absent. Diffusion uses the exact frozen conservative cell-centred no-flux Laplacian. Regional amount is the sum of active concentration times the fixed native cell width. The distal ROI contains cells whose centres satisfy $x_i \ge L/2$, including the entire centre cell for these odd $N=201$ meshes; it therefore begins at the left face of that centre cell, slightly before the ideal geometric half. The CSV records every endpoint index, cell count, $\Delta x$ and face coordinate. The frozen original-core mask is evaluated separately.

| L | Cue width | Off time | Distal active amount | Endogenous recruitment | Deactivation | Net reaction | Native diffusion |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 1.8 | 0 | 12.428 | 0.025158 | 0.21538 | -0.19022 | 9.23e-06 |
| 9 | 1.8 | 100 | 11.267 | 0.21913 | 0.19526 | 0.023876 | 9.28e-06 |
| 9 | 3 | 0 | 13.466 | 0.012981 | 0.23337 | -0.22038 | 6.54e-06 |
| 9 | 3 | 100 | 8.8569 | 0.13126 | 0.15349 | -0.022231 | -3.1e-05 |
| 25 | 3 | 0 | 35.798 | 0.1216 | 0.62037 | -0.49877 | 1.5e-14 |
| 25 | 3 | 100 | 31.383 | 0.52516 | 0.54388 | -0.018712 | 3.96e-15 |
| 25 | 5 | 0 | 43.484 | 0.089902 | 0.75357 | -0.66367 | -2.06e-12 |
| 25 | 5 | 100 | 37.524 | 0.63383 | 0.6503 | -0.016468 | -3.16e-12 |

At L9 withdrawal, the broader response has more distal active amount but less endogenous recruitment and more deactivation. Both responses initially lose distal activity. At 100 cue-free units the narrow response has positive net endogenous reaction, while the broad response remains negative. This makes the concentration-dependent native maintenance account concrete at those states. The small regional active-diffusion terms are instantaneous contributions across this particular ROI boundary; they do not establish that transport or inactive redistribution is unimportant for the preceding or subsequent competition.

The regional balance independently equals endogenous recruitment minus deactivation plus the exact sum of incoming face fluxes. Maximum absolute residual is 2.22e-16 native amount/time. The CSV also records actual neighbouring-frame secants and their intervals. In particular the forward 0.25-unit secant at withdrawal is a finite-time observation, not the instantaneous derivative; no exact agreement with it is assumed.

The L25 pair is contextual, with its distinct pulse 815 and width contrast 3 / 5; it is not a controlled change in length alone. Both L25 distal regions still have slightly negative net reaction at 100 cue-free units, although the broader challenge ultimately wins. These two snapshots therefore do not explain that width-fate difference or define a universal maintenance threshold. The table identifies native production/loss at each saved state without allocating the causal explanation of eventual replacement. The inherited $N=201$ spatial resolution has not undergone an independent convergence study; the native input and clocks remain uncalibrated to an optical actuator or reporter.

[Exact regional values and provenance](supporting/data/native_regional_state_law.csv) · [Source/artifact checks](supporting/metadata/validation.json)
