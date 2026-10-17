# Single-species testing

Testing contains two protocols beyond the single-cue calibration experiment.
Both use directly observable polarized-domain counts. A qualifying domain is a
connected interval above the midpoint of the active-GTPase profile, and the
whole profile must first satisfy `P=(max(u)-min(u))/mean(u) >= 0.5`. This prevents
small numerical ripples in an effectively homogeneous profile from being
reported as poles.

## Two seeded domains

`dual_pole_map.json` initializes equal Gaussian active enrichments of amplitude
2 at relative
coordinates 0.2 and 0.8, with matched removal from the cytoplasmic pool so that
`u+h*v` is conserved pointwise. It repeats the `(mu, Lambda, h)` map for
initial asymmetries `epsilon=0,0.001,0.01,0.05`.

The principal output is the number of qualifying polarized domains at 0, 100,
and 300 seconds. For example, `2 -> 2 -> 1` means two seeded domains persisted
through 100 seconds but only one remained at the endpoint. `2 -> 2 -> 2`
demonstrates persistence over the 300-second observation window, not formal
infinite-time stability.

## Polarity reversal

`polarity_reversal.json` runs a matched `(mu, Lambda, h)` sweep. A left cue acts
from 10--30 seconds and a matched right cue from 150--170 seconds. A run is
eligible to test reversal only when it has one polarized domain and left-end
contrast at most -0.5 immediately before the right cue. Eligible endpoints are
classified as reversed, retained original pole, lost polarity, multiple domains,
or ambiguous. Reversal time requires a single right-biased domain to persist for
ten seconds.

Formal continuation and eigenvalue analysis remain deferred until these regime
maps identify informative parameter regions.
