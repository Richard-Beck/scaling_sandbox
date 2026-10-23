# Figure 4: enlargement and heading memory

[Main figure](figure4_main.pdf) · [Legend/source contract](legend.json) · [Complete supplement](figure4_supplement.pdf) · [Methods](METHODS.md) · [Results](RESULTS.md)

Enlargement extends heading memory and raises noise tolerance, while the size advantage weakens as input strengthens. The weak-input drift in Wave-pinning, Holmes and Spring is close to quadratic in input, connecting Figure 4 to the deterministic edge-response test in Figure 3.

Direct event accounting distinguishes persistent secondary onsets from certifiable no-parent births. Certified takeover probability declines with size in each model at a common finite follow-up, but different certification fractions and pooled inputs limit population inference. Per-perimeter birth direction varies by model; it is not a universal consequence of enlargement.

At native moment amplitude, eligible weak-input heading loss is allocated almost entirely to domain motion. As secondary domains appear, membership changes and the amplitude/shape/background remainder can grow; the pattern varies with model and size. This exact accounting remains descriptive, and the original unit-normalized point-region artefact is retained as a supplementary diagnostic.

The complete existing campaign was analyzed: 260 conditions, 45,696 independent histories. No new noisy model runs were generated. Birth rates are certified lower bounds; exact no-parent birth and split classification was not retained in the compact schema. The requested strict all-birth/takeover measures therefore remain partly unidentifiable. `requested_measure_availability.csv` lists all cases. Finite-horizon takeovers use a common 20-noise-correlation-time follow-up. Six new supplementary pages are followed by the reviewed unscaled main and all 40 reviewed supplementary pages, unchanged (47 pages total).

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

The main heading allocation now uses native endpoint amplitude (SLURM 21333382, heading-only). Original unscaled metrics and reviewed main/supplement are preserved under `supporting/historical_unscaled_20261023/`. Main metric tables: `heading_amplitude_condition_metrics.csv`, `heading_amplitude_seed_metrics.csv`; fixed-cohort incidence sensitivity: `takeover_incidence_by_size.csv`.
