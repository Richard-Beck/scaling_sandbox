# Spring stimulation regime search

The spring mechanics and biochemical coefficients remain at the author-code defaults: basal activation 4, tension feedback 80, and the original diffusion and reaction parameters. We changed the stimulation schedules to use the spring model's native drive scale. All lengths use the same selected drive amplitude; only the tested finite-countercue duration varies.

## Search and selection

First, the established standalone positive control was reproduced at L0=1 with 201 cells, a right-edge Gaussian cue of amplitude 4 and width 0.2, and the original 1500-unit conditioning and recovery periods. Durations 0.15 and 0.17 did not reverse the front; 0.18 and 0.20 did. The new adapter reproduces that threshold at 0.18. Its initial orientation was -0.1047859, successful final orientation was +0.1047859, and the largest mass error among these duration tests was 2.70e-11. Results are in [spring_positive_control_pilot.json](spring_positive_control_pilot.json).

At L0=1, amplitude 4 also passed one pilot for each other assay. The two-end Gaussian competition cue (width 0.2, 40 units, 2% imbalance) left both endpoint peaks established at withdrawal and resolved competition in 20 units. The edge-to-edge linear-gradient reversal was eligible after 2000 units of leftward conditioning and reversed 12.5 units after the direction crossfade began. Alternating linear gradients at periods 100 and 400 followed in 12/12 switches each, with median response latencies of 3 and 9 units. The per-assay records are `spring_competition_pilot.json`, `spring_gradient_pilot.json`, and `spring_periodic_pilot.json`.

Amplitude 4 is therefore the selected fixed stimulus regime for all spring lengths and assays. No amplitude-8 or amplitude-16 search was needed. The complete, frozen 40-task length sweep uses L0 values 0.6, 0.8, 1, 1.2, 1.4, 1.8, 2.2, and 2.6; periodic following includes both periods. Finite countercue durations are the same at every length: 0.02, 0.05, 0.1, 0.15, 0.17, 0.18, 0.2, 0.3, 0.5, 1, 2, and 5 author-time units. This brackets the reference threshold closely while retaining shorter and longer tests.

## Scoring details

Finite-countercue eligibility requires a stable leftward polarity over the final 500 conditioning units and active-field contrast of at least 0.05 throughout that window. Success requires a stable rightward polarity and the same contrast floor throughout the final 500 recovery units. Competition reports `resolved_during_cue` separately when the loser/winner peak ratio is already below 0.1 at withdrawal; such a case has no post-withdrawal latency. Otherwise both peak excesses must exceed 0.1 at withdrawal, and the resolution ratio must remain below 0.1 for 50 units. The first 100 units after competition withdrawal and the gradient response window use 0.5-unit samples; finite Gaussian challenges are sampled at 0.002 to 0.01 units. Coarser samples are used for settled tails. Reported latencies are therefore quantized to the corresponding sampling interval.

## Fixed length sweep results

The fixed-amplitude campaign completed 40/40 tasks with status `success`. The table gives the observed endpoint for each L0; all values are in native author-model time. Latencies are sampled estimates, with a 0.5-unit interval for the early competition and gradient response windows and 1-unit periodic output for these periods.

| L0 | Competition | Gradient reversal | Finite countercue threshold | Period 100 | Period 400 |
|---:|---:|---:|---:|---:|---:|
| 0.6 | 42.5 | 12.0 | 0.15 | 3 | 9 |
| 0.8 | 27.0 | 12.5 | 0.15 | 3 | 8 |
| 1.0 | 19.5 | 12.5 | 0.18 | 3 | 9 |
| 1.2 | 16.5 | 12.5 | 0.3 | 3 | 9 |
| 1.4 | 15.0 | 12.5 | 0.3 | 4 | 9 |
| 1.8 | 12.0 | 12.5 | 0.5 | 4 | 9 |
| 2.2 | 16.0 | 13.0 | 0.5 | 4 | 9 |
| 2.6 | 25.0 | 13.0 | 0.5 | 4 | 10 |

Periodic following succeeded on 12/12 switches at every L0 and period. The largest mass drift across any spring result was 4.68e-9. The final array was SLURM job 20631850; it supersedes the first, coarser competition-latency run (job 20631761). Each final task JSON records its task and a consistent hash set for the regime adapter, shared mechanics adapter, and spring source files, under `data/task_00.json` through `data/task_39.json`.
