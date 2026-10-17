# Revision 8 completed — 2026-10-01

Both requested PDFs and PNG previews are in results/. All 441 preparations and
4096 response jobs finished with exit code 0. Final collection verified complete
duration coverage; no missing outputs or execution errors. There are 770 case
records: 517 success, 142 censored, 48 unreliable, 35 ineligible, 7 unsettled,
and 21 paused spring-competition records. All frozen source checks passed.

Cue strengths/shapes were checked against revision 6 by exact AST comparison of
the Cue class. Seven-model preparation and checkpoint smoke checks passed;
periodic common sampling was checked for spring and Goryachev. Both final PNGs
were visually inspected for layout and labels. PDFs/SVGs share the same renderer.
Artifact SHA-256 hashes are in results/artifact_hashes.json.

Compared with revision 6, the only categorical change is Otsuji competition at
25 micrometres: success becomes unsettled under the longer confirmation window
and unchanged 10000-unit cap. All 21 spring countercue thresholds are identical
to revision 6 despite cue conditioning and longer recovery. Some measured event
times and four wave-pinning/two Holmes finite thresholds differ; the per-assay
counts and spring comparison are in results/comparison_to_revision6.json.
Successful execution does not certify scientific interpretation.

Remaining protocol differences are explicitly retained in README.md: native
strengths/units/coupling, other cue shapes and preparation histories, periodic
initialization, competition withdrawal, and state-specific numerical mechanics
checks. Spring competition remains paused. Revision 6 is preserved as baseline;
revision 7 remains rejected and was not used as an input.

## Revision-8 visualization update — 2026-10-01

User requested mode-colored outcome ranges, length regressions, a merged plot
and model-by-metric slope heatmap. Nineteen Goryachev competition successes
(18 at zero, one within the first 2.5-unit sample) are classified as censored
with mode resolved_at_withdrawal; original statuses/endpoints remain in raw
records, original CSV columns and results/before_mode_revision/. This is not
conventional right-censoring and does not indicate failure to resolve. Updated
exports have 498 success, 161 censored, 48 unreliable, 35 ineligible, 7 unsettled
and 21 paused cases. All 66 eligible series have OLS lines (>3 valid points);
periods and periodic outputs are fitted separately. Combined plot:
revision8/results/length_response_and_settling.png; slope heatmap:
revision8/results/length_response_gradients.png. Sources and provenance are in
revision8_visualization/render.py and results/visualization_provenance.json.
Use this renderer for current exports; frozen build_report.py restores the
original classification. No simulations were rerun or modified.
