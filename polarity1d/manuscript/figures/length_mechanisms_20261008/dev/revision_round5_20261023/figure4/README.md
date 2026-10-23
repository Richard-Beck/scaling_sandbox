# Figure 4: enlargement and heading memory

[Main figure](figure4_main.pdf) · [Legend/source contract](legend.json) · [Complete supplement](figure4_supplement.pdf) · [Methods](METHODS.md) · [Results](RESULTS.md)

Enlargement extends heading memory and raises noise tolerance, while the size advantage weakens as input strengthens. The weak-input drift in Wave-pinning, Holmes and Spring is close to quadratic in input, connecting Figure 4 to the deterministic edge-response test in Figure 3.

Direct event accounting distinguishes persistent secondary onsets from certifiable no-parent births. Certified takeover probability declines with size in each model at a common finite follow-up, but different certification fractions and pooled inputs limit population inference. Per-perimeter birth direction varies by model; it is not a universal consequence of enlargement.

Signed heading-loss accounting closes numerically while showing substantial cancellation between point-region changes and the moment residual. This exposes the limitations of a point-domain reconstruction; it does not establish that one mechanism causes a specified fraction of native heading loss.

The complete existing campaign was analyzed: 260 conditions, 45,696 independent histories. No new noisy model runs were generated. Birth rates are certified lower bounds; exact no-parent birth and split classification was not retained in the compact schema. The requested strict all-birth/takeover measures therefore remain partly unidentifiable. `requested_measure_availability.csv` lists all cases. Finite-horizon takeovers use a common 20-noise-block follow-up. Every historical supplementary page is appended unchanged after the eight new pages.

Rebuild with the existing `scaling_sandbox` Conda Python:

```bash
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_validation_20261023.py
sbatch polarity1d/manuscript/code/mechanisms/figure4_round5_analysis_20261023.slurm
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_analysis_20261023.py collect
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_render_20261023.py
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_documents_20261023.py
```

Analysis used SLURM job 21332720 (260 single-core condition tasks, small QoS, no array throttle). Data are generated once and opening/manually refreshing the reader loads those saved results; there is no browser timer. The analysis does not mutate raw campaign data, basals, historical figures, reviews or author text.

Key CSVs: `birth_heading_condition_metrics.csv`, `birth_heading_seed_metrics.csv`, `takeover_by_condition_distance.csv`, `takeover_by_size_distance.csv`, `takeover_seed_counts.csv`, `weak_input_drift_exponents.csv`, `kymograph_centering.csv`, `plotted_D.csv`, and `requested_measure_availability.csv`. Per-condition retained event rows live under `supporting/condition_analysis/` and cover the entire campaign.
