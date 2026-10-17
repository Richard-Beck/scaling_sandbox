# Mechanical-adapter validation

This is a validation record for the spring and deBelly adapters, not a report of
the final seven-model screen. The center-length SLURM checks used the suite
runner at 201 cells: spring L0=1 and deBelly L=15 µm. The later full campaign
will regenerate final rows with the corrected runner.

## Sources, defaults, and units

- Spring uses the preserved Zmurchok, Collette, Rajagopal & Holmes (2020)
  Figure 4 author code, loaded through
  `reports/trial_report_assets/zmurchok2020/run.py`. Its reference resting
  length is L0=1 in author model-length units; time is author model time and is
  uncalibrated. The biochemical coefficients, basal activation, tension
  feedback, and diffusion values are unchanged as L0 is swept.
- deBelly uses `DeBellyParameters()` from `polarity1d/src/polarity1d/debelly.py`.
  Length is micrometers and time is seconds. The run keeps the released default
  `published_one_sided` boundary mode, algebraic fixed-point initialization,
  Rac optogenetic target, and all source-default coefficients. No symmetric
  extension or manuscript-tuned parameter set is substituted.

## Cue comparison

I checked the adapters against `rd_adapter.py` and `protocols.py`:

- Gradient reversal and finite countercue use the same nonnegative linear ramps
  as the reaction-diffusion models: left cue 0.1 times (1 - x/L), right cue
  0.1 times x/L. Their edge-to-edge difference is 0.1 and their spatial mean
  is 0.05.
- Periodic following uses a cosine modulation with mean 0.05 and peak-to-trough
  difference 0.1. The zero crossings occur at P/4 plus k times P/2. Runs now
  continue through 6.25P so the final switch has its full response window.
- Competition uses the same half-cosine edge profiles, 0.1L width, 0.2 peak
  drive, 20-unit plateau, 20-unit half-cosine withdrawal, and 2% relative
  left-to-right cue asymmetry as the reaction-diffusion adapter.
- Finite countercue is a rightward linear cue switched on after 2000 units of
  left cue and 2000 units without cue, then switched off after each tested
  duration, followed by 5000 units of recovery.

The numeric cue values are held fixed across lengths. In deBelly they add to
the Rac activation input. In spring they add to the author model's basal
activation before its original nonlinear exchange law. Each output row retains
native time and length labels; those units are not calibrated across models.

## Numerical checks and first observations

At 101 cells with periodic input for 1000 time units, relative spring total
active-plus-inactive mass drift was 6.61e-10. The deBelly MCA conservation
error was 4.44e-16. A separate n=201 run from the deBelly simulator reports its
MCA conservation error in trajectory metadata; the adapter now passes that
diagnostic into suite results. Spring result rows include relative total-mass
drift.

At 201 cells, complete periodic windows took 10.5 and 12.3 seconds for spring
at P=100 and P=400, and 0.6 and 2.4 seconds for deBelly. After correcting the
periodic scorer to reject states already aligned at switch time, spring L0=1
had no measurable switches at either period; six of twelve switches had been
incorrectly counted as zero-latency responses before that correction. deBelly
L=15 µm also had no measurable switch responses.

The first deBelly finite-countercue screen at L=15 µm was ineligible: the
published one-sided default did not retain a leftward orientation through the
two-thousand-unit stimulus-free interval. Its orientation remained near zero
after both conditioning and the tested countercues. This is an ineligible
reversal test under the shared protocol, not evidence that an eligible reversal
failed. The one-sided boundary setup can bias left/right comparisons and should
be kept visible in interpreting its rows.

The initial SLURM smoke attempt also exposed a suite-adapter bug: nonperiodic
tasks passed `period=None` into the schedule constructor. That has been fixed by
giving unused periods a neutral default; the failed competition/reversal smoke
rows are not valid results. The full suite should be used for final plotted
data.

In the completed 275-task suite, the 40 spring rows contain five competition
resolutions and three ineligible competitions; all eight gradient-reversal and
finite-countercue rows are censored, and all sixteen periodic rows are censored.
The maximum spring mass drift is 1.08e-8. The 40 deBelly rows contain eight
ineligible competitions; gradient reversal has five ineligible and three
censored rows; finite countercue has four ineligible, three censored, and one
success at L=21 µm with a 5120-second tested cue duration; all sixteen periodic
rows are censored. The maximum deBelly MCA conservation error is 4.44e-16.
These are endpoint-screen outcomes at fixed default coefficients, not fits or
claims that the two models share calibrated length or time scales.
