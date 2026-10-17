# deBelly regime search and reproduction notes

## Positive control

The first check uses the established Figure 2 baseline: default model
coefficients, Julia activation defaults (`rac_basal_activation=1`,
`rho_basal_activation=1`), `reflection_symmetric_extension`, algebraic fixed
point initialization, and cell-centred input. At L=21 um and 201 cells, it
applies the left half-cosine pulse (amplitude 3.6, width 0.2 L) from 0 to 600 s,
then the right pulse from 1310 to 1610 s, followed by 10,000 s of recovery.
The previous 105-cell Figure 2 trajectory used this same 1310 s onset; its
published x50 audit reported a mirrored target of 0.552895 and successful
reversal at 300 s.

The 201-cell rerun reproduced the intended polarity change. The pre-challenge
x50 was 0.445341 L (mirrored target 0.554659 L); final x50 was 0.552938 L.
Using the previous x50 criterion of within 0.05 L of the mirrored target for
50 s, the target was reached at 1655 s, 345 s after cue onset. That is after
the 300 s pulse ended, so the response completed during recovery. The final
first-moment orientation was +0.06328, final active-field range was 0.22182,
and the orientation range during the final 500 s was 1.11e-8. MCA conservation
error was 4.44e-16. The common first-moment threshold +0.05 was crossed at
1325 s, 15 s after cue onset; that timing differs from the established x50
endpoint and must be reported as a separate diagnostic, not substituted into
the previous threshold estimate.

The historical onset is kept as a positive-control schedule. The pre-challenge
first-moment range from 810 to 1310 s was 0.01391, above the shared suite's
0.005 settling criterion. A separate fixed-wait audit at 10600 s passed: over
9600–10600 s, Rac x50 range was 3.48e-11 L, MCA x50 range was 7.78e-11 L,
orientation range was 3.34e-11, and final orientation was -0.06328. Thus the
shared finite-cue tests use the same fixed 600 s conditioning pulse followed by
10000 s without cue at every length. There is no length-specific start time.

Raw 201-cell data are in `pilot_debelly_l21.npz`; the reproduction runner is
`debelly_regime.py`. This pilot changes only the boundary formulation to the
established symmetric baseline and uses the activation defaults already used
in the prior Figure 2 campaign. No parameter or success-threshold tuning was
done.

## Bounded assay regime search and selection

All bounded trials used L=21 um and 201 cells. They varied stimulation only;
all model coefficients and the predeclared common thresholds stayed fixed.
The selected amplitude is the same at every length in the full screen.

| Assay | Bounded L21 regimes | Observation | Selected full-screen regime |
|---|---|---|---|
| Finite countercue | Localized amplitude 3.6, width 0.2 L; durations 100, 150, 200, 250, 300, 400, 600, 800, 1000 s | 100 and 150 s briefly visit the mirrored x50 target but return to the original state. 200 s ends near uniform (orientation +0.0024, contrast 0.0083). 250 s is the first tested autonomous reversal (orientation +0.0633, contrast 0.2218); the prior ~232 s threshold lies between 200 and 250 s. | Same pulse shape and amplitude; fixed 600 s conditioning + 10000 s free settling; listed duration grid; 10000 s recovery. |
| Gradient reversal | Full-domain linear gradients with amplitude 0.1, 1, 3.6; 2000 s conditioning and 20 s crossfade | All were eligible. Latencies from crossfade start were 2675, 90, and 45 s, respectively. | Amplitude 3.6; retain 2000 s conditioning, 20 s crossfade, and 20000 s response horizon. |
| Periodic following | Linear sinusoidal gradients at amplitudes 0.1, 1, 3.6, periods 100 and 400 s | At amplitude 0.1 neither period switched. At amplitude 1, P100 had 0/12 and P400 had 12/12 responses. At amplitude 3.6 both periods had 12/12. | Amplitude 3.6 at P=100 and 400 s; six cycles plus the quarter-period confirmation tail. |
| Competition | Dual localized 3.6 pulses, width 0.2 L, 2% left-right imbalance; durations 80 and 600 s | Both produced two eligible peaks at withdrawal. Neither had a qualifying loser/winner ratio below 0.1 for 50 s while the winner remained above 0.1. At 20,000 s after withdrawal, the 80 s case was nearly uniform (peak excesses 0.0185/0.0034); the 600 s case retained only a weak bias (0.0738/0.0060, below the winner threshold). | Use the 80 s protocol and report the outcome as censored if no relative-resolution event occurs. No amplitude or threshold increase was used to force an endpoint. |

The final bounded traces, including full pilot fields for inspection, are in
`pilot_debelly_l21.json`, `pilot_debelly_l21.csv`, and the corresponding
`pilot_debelly_*.npz` files. The common periodic score counts a response only
when the starting orientation has the opposite sign to the upcoming cue and
the new orientation and active contrast exceed 0.05 for 25 s within that
half-cycle. Competition resolution uses the ratio below 0.1 and winner peak
above 0.1 for 50 s; early resolution is kept separate from post-withdrawal
latency. The selected deBelly suite runner is `debelly_regime.py`.

The frozen 40-task length campaign is in `debelly_manifest.json` and writes to
`debelly_data/`. It uses lengths 9, 12, 15, 18, 21, 27, 33, and 39 um with 201
cells. All 40 tasks completed without execution failures: 32 successes, five
censored competition outcomes, and three ineligible finite-countercue
outcomes.

| Length (um) | Competition (s after withdrawal) | Gradient reversal (s) | Finite countercue (s) | Period 100 (successes/12; median s) | Period 400 (successes/12; median s) |
|---:|---:|---:|---:|---:|---:|
| 9 | censored | 40 | ineligible | 12/12; 25 | 12/12; 45 |
| 12 | censored | 40 | ineligible | 12/12; 25 | 12/12; 45 |
| 15 | censored | 45 | ineligible | 12/12; 25 | 12/12; 45 |
| 18 | censored | 45 | 400 | 12/12; 25 | 12/12; 45 |
| 21 | censored | 45 | 250 | 12/12; 25 | 12/12; 45 |
| 27 | 7755 | 45 | 400 | 12/12; 25 | 12/12; 45 |
| 33 | 3675 | 45 | 600 | 12/12; 25 | 12/12; 45 |
| 39 | 6070 | 45 | 800 | 12/12; 25 | 12/12; 45 |

At L=9, 12, and 15 um, the left-conditioning state loses autonomous polarity
during the fixed 10000 s free interval: the pre-challenge orientations are
-1.1e-13, -3.5e-9, and -8.6e-5, respectively. Their contrast is negligible,
so finite reversal is correctly ineligible. No length-specific settling time
or parameter change was used. Gradient reversal and periodic following
succeed at every length. Competition is length-dependent: L=27, 33, and 39 um
resolve after withdrawal with a persistent winner; at smaller lengths both
endpoint peaks decay below the winner threshold before qualifying resolution.
At the reference L=21 pilot, 80 s and 600 s cues established two peaks but did
not produce a scored resolution; the selected 80 s regime is censored there.
