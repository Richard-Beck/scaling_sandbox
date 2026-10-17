# Refined-data manuscript figures

- `main_03_migration_OU.png/.pdf/.svg`: updated integrated manuscript OU figure. A is the compact nine-summary comparison; B is the full population MSD and fitted OU curves beside A, with logarithmic axes and an annotated population crossover; C contains the enlarged simulated trackplots alone. D–F retain the relative loss surfaces, paired parameters, and five prior/pooling sensitivity scenarios. The summary MSD fold uses the last overlapping population lag, interpolated from both native traces; it is a display summary rather than a replacement fitting target. Crossover interpolation is descriptive raster-derived information, not a biological uncertainty estimate.
- `supp_04_Hosny_sensitivity.png/.pdf/.svg`: updated existing six-panel Hosny sensitivity supplement. Its supporting broad/focused surfaces, multistart drift profiles and leave-Hosny-out comparison are retained in `../results/hosny_sensitivity/`.
- `xuan_population_MSD_diagnostics.png/.pdf/.svg`: added full-curve shape/residual diagnostic, including equal log-lag fitting regimes and explicitly labeled raster extraction bounds.

Each figure has a separate caption and provenance JSON. `main_03_observation_folds.csv` retains central-value statistic labels and group-specific quartile, raster and SEM fields without inventing ratio confidence intervals. `figure3_persistence.csv` and `figure3_primary_noise.csv` retain the established parameter table names; values are the new fits. `xuan_population_MSD_residuals.csv` contains every fitted population-curve residual and weight.

Regenerate from the repository root in the existing `scaling_sandbox` environment:

```sh
python wgd_langevin/refined_ou_20261003/scripts/main_surface.py
python wgd_langevin/refined_ou_20261003/scripts/build_manuscript_figures.py
python wgd_langevin/refined_ou_20261003/scripts/curve_diagnostics.py
sbatch wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/run.slurm
```

Read the campaign README, results summary and movie-diagnostic documentation together with these captions. Track illustrations in the main figure are simulated population-protocol examples. The representative movie diagnostics are separate and their physical timing remains conditional. Loss surfaces and prior/pooling/profile contrasts are heuristic sensitivity analyses, not likelihoods or confidence intervals. Historical manuscript artifacts outside this folder are unchanged.
