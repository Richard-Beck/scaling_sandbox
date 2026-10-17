# Validation-only amendment: flat cached activity-site labels

The verified native producer, all 243 archives, frozen scientific plan, noise
law, statistical estimators, time windows, bootstrap draws/seed and medoid
selection remain unchanged. Failed statistics job 21126773 and its complete
workflow/source/status are preserved under ../statistics_workflow/.

The original strict 3e-10 centroid identity assertion compared cached labels in
flat states as if each were a current activity-site centroid. In Goryachev C20
high-mean deterministic control slot0223, 697 saved profiles are numerically
flat (contrast <=1e-12, zero current regions). The native last-valid label is
updated at scalar cadence5, while the independent saved-profile reconstruction
uses cadence40. Labels differ by3.0027896791295916e-5, with exactly zero position
difference on all nonflat profiles. These flat labels describe no identifiable
current site and were already ineligible for the site-displacement estimator.

The separate analysis.py retains the original strict tolerance on every
nonflat profile, independently asserts exact cache carry-forward at each
cadence, confirms flat/current-region status, and asserts no flat frame is
eligible for physical site comparisons. Per-case counts/errors are exported
inside each source record of the full analysis audit. Cached positions remain
raw metadata; no heading or current physical site is inferred from them.

REGRESSION.json covers exact slot0223 and resolved Wave-pinning/Goryachev/Holmes
cases. Resolved outputs are array-identical to the failed implementation.
Every analysis function except load_case (validation) and analyze (audit
metadata) is AST-identical. Original correlation/moment/peak/gating fixtures
pass. Sources and failed-status/producer/verifier identities are frozen in
preparation.json before the single statistics/report-only SLURM retry.
