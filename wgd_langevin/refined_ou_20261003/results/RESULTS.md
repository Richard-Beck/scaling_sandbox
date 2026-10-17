# Refined OU fit results

Generated 2026-10-03T17:13:46.416721+00:00. All 40 bounded starts and frozen trajectory validation finished.

|Scenario|Start|Data score|Pooling|Prior|Total|Converged starts|
|---|---:|---:|---:|---:|---:|---:|
|primary|1|17.289|7.4612|15.139|39.889|8/8|
|prior5|6|16.873|6.8643|18.718|42.455|8/8|
|prior20|1|17.936|8.3275|13.949|40.212|8/8|
|strong|2|20.89|7.5849|21.462|49.937|8/8|
|weak|7|13.764|5.9771|10.498|30.239|8/8|

Primary estimates (hours; noise in micrometers/hour^(3/2)):

|Study|tau reference|tau WGD|tau fold|noise reference|noise WGD|noise fold|
|---|---:|---:|---:|---:|---:|---:|
|Xuan2018|0.077231|0.71415|9.247|129.77|19.269|0.14848|
|Mallin2023|0.049694|2.0375|41.001|159.23|5.7048|0.035827|
|Jemaa2017|0.040063|0.5543|13.836|238.63|49.075|0.20565|
|Jemaa2023|0.03427|0.5193|15.153|61.568|13.037|0.21174|
|Hosny_baseline|0.043024|0.53269|12.381|98.359|16.132|0.16401|

Hosny signed drifts: 1.5532, 3.3173 micrometers/hour.

Xuan curve deviations (log ratio RMS; each regime contributes equally):

|Group|Regime|Points|RMS log(predicted/observed)|Mean log ratio|Inside raster envelope|
|---|---|---:|---:|---:|---:|
|0|early|87|0.27131|0.267|0/87|
|0|intermediate|78|0.13484|-0.066579|21/78|
|0|late|87|0.25871|-0.25652|0/87|
|1|early|86|0.12036|-0.008791|18/86|
|1|intermediate|86|0.13959|0.13641|0/86|
|1|late|86|0.12534|-0.053889|28/86|

Primary regularized Jacobian rank: 22/22. Active numerical bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]. These are numerical sensitivity diagnostics, not uncertainty intervals.

Validation files retain two frozen 16,384-track Sobol ensembles per scenario; all fitting curves are analytical. Sampled Xuan MSD validation uses exact Gaussian marginal displacements at every curve lag with common standardized draws, so agreement tests the marginal formula, not independent full-path correlations. Movie trajectory diagnostics are separate. Hosny nonmedian path predictions are undefined under the analytical approximation and only sampled trajectories are used for those checks. Mallin central statistic remains unspecified; the historical mean approximation is retained. Early 600-second Xuan checks overlap the fitted curve and are correlated checks.

Objective tolerances are unchanged historical descriptive scales recentered on reviewed values; raster bounds and SEM are retained separately. Xuan native digitized points are normalized within three equally spaced log-lag regimes; each regime contributes one third of a group curve block, and each group curve block has total weight one, matching one speed target. Raster envelopes are not population confidence intervals.

Data-only Jacobian ranks (relative singular-value threshold 1e-6; central finite differences at three relative steps):

|Scenario|step 1e-4|step 1e-5|step 1e-6|
|---|---:|---:|---:|
|primary|18/22|18/22|18/22|
|prior5|18/22|18/22|18/22|
|prior20|18/22|18/22|18/22|
|strong|18/22|18/22|18/22|
|weak|18/22|18/22|18/22|

The remaining rank is supplied by regularization. Each Jemaà assay observes only speed, leaving a persistence/noise tradeoff per group. Fitted times are therefore conditional on prior and pooling assumptions.

The Xuan reference curve shows systematic early overprediction and late underprediction (mean log ratios +0.267 and -0.257): the fitted OU curve shape does not reproduce the full digitized curve within its raster envelope. This discrepancy is larger than digitization precision and survives the completed optimization; completion is not evidence of model adequacy.
