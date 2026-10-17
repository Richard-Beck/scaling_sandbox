# Revision 2: validated stimulation regimes

The revised comparison contains 275 task rows and no missing or duplicate task indices. It reuses the original 160 rows for wave-pinning, LEGI, Otsuji, and Goryachev verbatim, and replaces all 115 spring, Holmes Model 3, and deBelly rows with campaigns using fixed, validated stimulation regimes. Final task data are in [`combined_data/`](combined_data/); [`manifest.json`](manifest.json), [`protocol.json`](protocol.json), and the report manifest record task keys, scoring rules, and model-specific runtime provenance.

The all-model outcome counts are 182 successes, 77 cases with no qualifying endpoint within the tested window, and 16 ineligible starting states. There were no failed executions. Per-model status counts are:

| Model | Success | No endpoint in tested window | Ineligible |
|---|---:|---:|---:|
| wave-pinning | 36 | 4 | 0 |
| LEGI | 19 | 18 | 3 |
| Otsuji | 6 | 30 | 4 |
| Goryachev | 20 | 14 | 6 |
| Holmes Model 3 | 29 | 6 | 0 |
| spring | 40 | 0 | 0 |
| deBelly | 32 | 5 | 3 |

“No endpoint in tested window” describes right-censored observation: it does not claim the system can never reach an endpoint. An ineligible case failed the assay's stated starting-state requirement and is kept separate from a stimulus-response failure. Periodic tests report the median only for qualifying switches; detailed switch outcomes are in [`periodic_switches.csv`](results/periodic_switches.csv).

## Selected regimes and remaining limits

- **Spring:** original coefficients and native drive conventions, with fixed stimulus amplitude 4 for all lengths. The native-scale Gaussian reference reproduction found the shortest tested successful finite cue at 0.18 author-time units at L0=1; the eight-size campaign reports 0.15–0.5 depending on L0. Competition has a 2% initial two-peak imbalance; periodic following succeeds on all tested switches at periods 100 and 400. All 40 tasks succeeded. Regime bounds, response sampling precision, and full measurements are in [SPRING_SEARCH.md](SPRING_SEARCH.md).
- **Holmes Model 3:** corrected no-PI Eq. 5–6 kinetics, fixed Rac activation 3.5, competition cue amplitude 0.2, and gradient/finite/periodic cue amplitude 1. The chosen activation lies inside the reference-length spatial-polarization band; it is not a claimed published default. Across seven lengths, competition resolves only at 10 μm (the other six cases have no qualifying resolution within 20,000 s). All direct reversals, finite countercues, and periodic switches qualify. See [HOLMES_SEARCH.md](HOLMES_SEARCH.md), including the correction to the first-pass kinetic formula and the archived superseded run.
- **deBelly:** restored the reflection-symmetric boundary extension and established Figure 2 activation defaults; fixed cue amplitude 3.6 and the established pulse shape. Finite-countercue durations were screened at L=21 μm with no coefficient or threshold tuning; 250 s was the shortest tested successful duration there. Across eight lengths, finite reversal is ineligible at 9, 12, and 15 μm because the conditioned polarity decays during the fixed settling interval; eligible lengths succeed at tested durations 400, 250, 400, 600, and 800 s. Competition resolves only at 27, 33, and 39 μm, with no qualifying resolution at smaller sizes. See [DEBELLY_SEARCH.md](DEBELLY_SEARCH.md) for the bounded trials and endpoint details.

The three revised campaigns fixed each model's stimulation regime across length; no per-length coefficient or amplitude tuning was used. Finite-countercue results are the shortest successful duration on the tested grid, not interpolated estimates. Their bounds and individual attempts are in [`finite_duration_tests.csv`](results/finite_duration_tests.csv). Latencies are limited by each adapter's saved sampling interval; see the spring search note for its especially fine, short-response sampling details.

## Figures and tables

- [Seven-model overview (PNG)](results/four_test_all_models.png) · [PDF](results/four_test_all_models.pdf)
- [Peak competition](results/competition.png) · [Gradient reversal](results/gradient_reversal.png)
- [Finite countercue](results/finite_countercue.png) · [Periodic following](results/periodic_following.png)
- [Endpoint summary](results/endpoints.csv) · [finite-duration attempts](results/finite_duration_tests.csv) · [periodic switches](results/periodic_switches.csv)

Each plot retains native length and time units. The spring L0 axis remains in the author's dimensionless length unit and has no physical calibration to the micrometre axes used by other models. Numeric y ticks on the compact log-scale panels are sparse display ticks; the underlying endpoint values are retained in the CSV tables.
