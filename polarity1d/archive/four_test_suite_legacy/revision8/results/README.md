# Revision 8 outcome modes and length regressions

Generated 2026-10-01T16:16:56.494964+00:00. Scientific simulations are complete; no new simulations were run.

[Combined preparation/response PNG](length_response_and_settling.png) · [PDF](length_response_and_settling.pdf)
[Gradient heatmap PNG](length_response_gradients.png) · [PDF](length_response_gradients.pdf)
[Response PDF](four_test_all_models.pdf) · [Settling PDF](time_to_settle.pdf)

Goryachev competition: 19 successes reclassified as censored with a distinct
resolved-at/near-withdrawal mode. Eighteen satisfy the resolution criterion at
withdrawal; one resolves within the first 2.5-unit response sample. This is an
analysis classification requested by the user, not conventional right-censoring
or evidence of failure to resolve. The original successes and endpoints remain
in raw data, before_mode_revision/, and original_status/original_endpoint CSV
columns. Colors distinguish this mode, unobserved endpoint, ineligible state,
unsettled preparation, unreliable periodic fit and paused assay. Periodic bands
are lower half P100 / upper half P400. Settled preparation times remain valid
even when the prepared state is assay-ineligible.

OLS lines use >3 valid points, include an intercept and stay within their valid
length extent. Periods and amplitude/phase outputs are fitted separately. Linear
fits describe valid subsets; they neither estimate censored times nor demonstrate
a mechanistic scaling law. Fits can span invalid lengths between observations.
The slope heatmap includes all response and settling metrics, including spring's
additional spontaneous gradient preparation. Raw slopes are labeled and use
separate symmetric scales by metric, because time, amplitude and phase differ.
Native model/length/time units differ across models. Gray means no eligible fit
(≤3 valid observations or no applicable phase), not a zero slope.

[Regression table](length_regressions.csv) · [Goryachev flags](goryachev_withdrawal_flags.csv)
[Endpoints and modes](endpoints.csv) · [Settling records](settling_times.csv)
[Provenance](visualization_provenance.json). Original exports are in before_mode_revision/.
Frozen revision-8 simulation inputs remain intact. Use ../..//revision8_visualization/render.py
for these exports; the frozen campaign renderer recreates the original analysis.
