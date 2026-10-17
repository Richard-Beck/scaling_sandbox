# Unused diagnostic variant — not the operative analysis

This folder preserves an independently developed observability-only adapter
and a complete 243-case/194643-profile regression. It was **not submitted** and
must not be confused with the operative retry.

The actual retry is SLURM **21126788**, using
`../statistics_centroid_amendment_20261017/analysis.py`, SHA256
`61622d8d21102397d38e9dc4050cc83852a60169c0e08fdedff494e1ab415dc7`.
Its separately frozen workflow preserves the failed 21126773 sources/logs.

The unused variant provides corroborating diagnostic evidence: only Goryachev
C20 high-mean control slot0223 has unobserved centers (697 saved profiles and
5572 scalar frames), no such frame meets the unchanged physical/direction
eligibility gates, and every identifiable saved center agrees exactly. No
native trajectory, parameter, stimulus law or statistic was changed. Its
source hash and complete regression are preserved in `REGRESSION.json`.
