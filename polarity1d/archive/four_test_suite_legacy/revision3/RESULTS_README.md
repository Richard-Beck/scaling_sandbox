# Revision 3: endpoint assays and periodic harmonic response

Revision 3 retains the 165 nonperiodic task rows byte-for-byte from revision 2 and replaces its 110 periodic-following rows with the archived first-harmonic amplitude and phase reliability analysis. Revision 2 remains unchanged as the record of its periodic switch-latency scoring; use revision 3 as the latest comparison.

The 165 endpoint cases comprise 106 successes, 43 with no qualifying endpoint in the tested window, and 16 ineligible starts; no execution failures occurred. Of the 110 periodic fits, 92 pass the archived reliability gate and 18 do not. Counts by model are:

| Model | Periodic fits passing reliability / total |
|---|---:|
| Wave-pinning | 16/16 |
| Otsuji | 0/16 |
| Goryachev | 16/16 |
| LEGI | 14/16 |
| Holmes Model 3 | 14/14 |
| Spring | 16/16 |
| deBelly | 16/16 |

All 16 Otsuji fits and two LEGI fits are retained with numeric estimates and rejection reasons in the CSV. The phase/amplitude plots omit them from their numeric scales and mark them in the bottom status strip.

The seven-model overview has three endpoint columns—competition, gradient reversal, and finite countercue—and two periodic readouts: harmonic amplitude and phase lag. A periodic response's phase is reported as a signed fraction of its stimulus period; positive values mean the response lags. The harmonic fit window is P through 6P inclusive, so the initial stimulus cycle is excluded. Periods 100 and 400 are distinguished by color and refer to each model's native time units.

The harmonic measurement uses the archived outer-edge orientation: compare mean active concentration in the outer 20% of each cell end. The shared scorer projects that signal onto the first harmonic, then applies the archived reliability requirements for amplitude, phase lag, cue correlation, harmonic R-squared, cycle-to-cycle amplitude variation, and phase consistency. Unreliable values remain in [`periodic_harmonics.csv`](results/periodic_harmonics.csv), together with their failure reasons and all scorer metrics, but do not appear on either plotted metric scale. A cross in the labeled plot strip marks an unreliable fit. These cases are not missing trajectories. Per-task traces are retained in `diagnostics/`.

The periodic starting states are unchanged from revision 2 and do not form a controlled common preparation. Wave-pinning, Otsuji, LEGI, Goryachev, and deBelly begin from their model's default homogeneous or algebraic fixed-point state; spring and Holmes Model 3 begin from a prepared left-polarized state. That remaining starting-state difference can affect amplitude and phase reliability. The selected stimulation regimes and preparation details are documented in [SPRING_SEARCH.md](../revision2/SPRING_SEARCH.md), [HOLMES_SEARCH.md](../revision2/HOLMES_SEARCH.md), and [DEBELLY_SEARCH.md](../revision2/DEBELLY_SEARCH.md). Preparation duration, cue details, native parameters, sampled traces, and source hashes are also recorded in each periodic task JSON under [`data/`](data/) and in [`report_manifest.json`](results/report_manifest.json).

The archived reliability gate requires amplitude at least 0.05, absolute phase lag no more than 0.25P, cue correlation greater than zero, harmonic R-squared at least 0.8, cycle-amplitude coefficient of variation no more than 0.25, and cycle-phase dispersion no more than 0.1. These are the thresholds in the original [WP1 tracking-reliability audit](../../../polarity1d/scripts/analyze_wp1_tracking_reliability_audit.py), implemented by [`harmonic.py`](harmonic.py). All models except spring save 80 samples per period. Spring retains its 1-author-time-unit output interval and maximum solver step, providing 100 samples per P=100 cycle and 400 samples per P=400 cycle; the scorer uses those native samples directly.

The endpoint table [`endpoints.csv`](results/endpoints.csv) contains only the 165 nonperiodic cases. Periodic cases have no latency endpoint in this report. [`periodic_harmonics.csv`](results/periodic_harmonics.csv) contains all 110 periodic fits, including unreliable raw fit values. The 275-task identity and the byte hashes of all reused rows are in [`collection_manifest.json`](collection_manifest.json).

## Figures and tables

- [Seven-model comparison (PNG)](results/four_test_all_models.png) · [PDF](results/four_test_all_models.pdf)
- [Periodic amplitude and phase (PNG)](results/periodic_harmonics.png) · [PDF](results/periodic_harmonics.pdf)
- [Competition](results/competition.png) · [Gradient reversal](results/gradient_reversal.png) · [Finite countercue](results/finite_countercue.png)
- [Nonperiodic endpoints](results/endpoints.csv) · [all periodic harmonic metrics](results/periodic_harmonics.csv)

The three individual nonperiodic assay figures are copied unchanged from revision 2. The all-model figure plots native length and time units. Spring L0 remains an author-model length unit without a physical calibration to micrometres. Period labels state native time units rather than implying a shared seconds calibration.
