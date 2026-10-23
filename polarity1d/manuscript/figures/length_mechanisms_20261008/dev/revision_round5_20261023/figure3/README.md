# Figure 3 round 5

[Main PDF](figure3_round5.pdf) · [complete supplement](figure3_round5_supplement.pdf) · [legend](legend.md) · [methods/results](METHODS_RESULTS.md) · [reader](index.html).

Cue placement relative to the incumbent boundary separates sliding from takeover, and enlargement slows half-turn reorientation. In the weak-input deterministic response the two boundaries share a long-lived translation: the prescribed independent-edge approximation gives median prediction/observation ratios 0.511, 0.540 and 0.534 in Wave-pinning, Holmes and Spring. Existing noise exemplars confirm correlated boundary motion. Goryachev does not give a general drift-based explanation of heading persistence.

Deterministic campaign: `polarity1d/manuscript/data/runs/figure3_round5_edge_pulses_v2_20261023`. SLURM array 21332598, 120 independent one-core tasks, 4 GB each, small QoS, 45-minute limit. Existing-profile edge analysis: SLURM 21332996 (cohort-corrected; first pass 21332981), one core / 4 GB. No new noise campaign. First failed API attempt is preserved under `figure3_round5_edge_pulses_20261023`, plus the incomplete preparation under `_prepare_failed`; scientific v2 inputs are frozen and valid. An unused solver-only retry in `figure3_round5_edge_numeric_retry_20261023` was canceled after the original strict run completed; no retry output supplies a figure datum.

Rebuild using the existing scaling_sandbox Python:

```sh
python -B polarity1d/manuscript/code/mechanisms/figure3_round5_render_20261023.py
python -B polarity1d/manuscript/code/mechanisms/figure3_round5_checks_20261023.py
```

Every preceding supplement and main page is appended unchanged after the new supplement pages. Source hashes establish A/B preservation. No timed refresh. All data tables retain undefined statuses.
