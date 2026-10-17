# Generating trails for the two trial reports

> Current-output note, 2026-09-28: `reports/trial_report.pdf` and its entry-point
> builder now produce **Overall cell size and resistance to repolarization**.
> The historical two-report provenance below describes the files before that
> replacement. The new standalone report's build instructions and source hashes
> are in `reports/trial_report_assets/standalone/REPRODUCE.md` and
> `standalone/provenance.json`. Its content, figures, and renderer are the `v3`
> source files; the `v2` sources and numerical outputs remain available. To rebuild
> the historical revised report, stage its asset files and call
> `build_figures_v2.figure1/figure2/figure3` followed by `render_v2.main`, rather
> than using the current entry-point wrapper. `docs/trial_report.pdf` is untouched.

Investigated and verified 2026-09-28. Both PDFs were rebuilt in isolated
temporary directories using the existing `scaling_sandbox` Conda environment.
Each rebuilt PDF is byte-for-byte identical to its corresponding saved PDF
after replacing only the PDF `/CreationDate` value. No existing report or
analysis output was overwritten during verification.

## Original version: now `docs/trial_report.pdf`

- Title: **Cell length, polarity memory, and persistent cancer-cell migration**.
- 10 pages, 667,734 bytes; embedded creation time 2026-09-26 17:00:58 EDT.
- SHA-256: `bcb658008d7303cdf9d04d00518bd6e13ab8c2ab0814cb2fdf7524f8a398c08c`.
- Originally built at `reports/trial_report.pdf`, as confirmed by the recorded
  successful command and contemporaneous `pdfinfo` output.
- Original builder recovered as
  `reports/trial_report_assets/original_source/build_trial_report.py`.
  It contains the report text, Figure 3 construction, PDF rendering, and
  provenance export. Its original location was `reports/build_trial_report.py`.
- Recovery record: `reports/trial_report_assets/original_source/recovery.json`.

Inputs and assembly:

1. Figure 1 embedded unchanged from
   `polarity1d/results/figures/Figure1_WP1_cross_model.png`.
2. Figure 2 embedded unchanged from
   `polarity1d/results/figures/Figure2_standalone_v2_julia_activation_defaults_linear50.png`.
3. Figure 3 built from the best stored Gaussian-hierarchy fits in
   `wgd_langevin/gaussian_hierarchy/results/`: `fit_05.json` (primary),
   `fit_09.json` (prior5), `fit_18.json` (prior20), `fit_26.json` (strong),
   and `fit_33.json` (weak). The builder selects by minimum `total_score`
   within each scenario.
4. Matplotlib `PdfPages` assembled the manuscript. The successful original
   command was `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/build_trial_report.py`.

Surviving outputs in `reports/trial_report_assets/` include `figure3_wgd.png`,
`figure3_wgd.pdf`, `figure3_observations.csv`, `figure3_persistence.csv`,
`pdf_text.txt`, `page-01.png` through `page-10.png`, and `contact_sheet.png`.
`provenance.json` was generated at 2026-09-26 21:01:04 UTC. All nine recorded
source hashes still match. Fresh plain `pdftotext` output from the PDF in
`docs/` matches `pdf_text.txt` exactly.

The original source and successive edits are recorded in
`/home/4473331/.codex/sessions/2026/09/26/rollout-2026-09-26T16-45-06-01a0df76-e574-71b2-9281-11bb54cb55a9.jsonl`:
line 116 contains its initial source; lines 129, 138, 149, and 163 contain
subsequent edits; line 169 records the successful build; line 179 records
the original PDF metadata. Replaying those edits recovered the final builder.

The recovered script retains its original relative-directory assumptions and
output paths. **Do not execute it directly at its archive location.** For a
separate rebuild, create a staging root, copy it to
`STAGING_ROOT/reports/build_trial_report.py`, and link or copy the repository's
`polarity1d/` and `wgd_langevin/` into that root. Its outputs will then be under
`STAGING_ROOT/reports/`. This is the arrangement used for verification.

## Revised version: `reports/trial_report.pdf`

- Title: **Does cell size set polarity memory?**
- 13 pages, 1,028,728 bytes; embedded creation time 2026-09-26 19:28:59 EDT.
- SHA-256: `01c96bfe833dd7269c5442527f82741e4f72f5ef9554bde62b0bd80aa0217857`.
- Entry point: `reports/build_trial_report.py`, now a wrapper for
  `trial_report_assets/build_figures_v2.py` and `trial_report_assets/render_v2.py`.
- Manuscript source: `reports/trial_report_assets/report_content_v2.py`;
  exported text: `manuscript.md` and `pdf_text_v2.txt` in the same directory.

Inputs and assembly:

1. Figure 1 combines retained canonical-model tables under
   `polarity1d/results/data/manuscript/` with new transport/pool interventions
   under `reports/trial_report_assets/size_mechanism/`.
2. `size_mechanism.py` generated and ran a 432-task campaign. Recorded SLURM
   job **20597822** used `size_mechanism.slurm` (54 array elements, eight tasks
   per element, one CPU, 2 GB, QoS `small`).
3. `refine_size_mechanism.py` added 48 matched-initial-profile refinement
   runs, recorded SLURM job **20597876** (six array elements, one CPU,
   2 GB, QoS `small`). Task JSONs, `summary.csv`, the manifest, and logs survive.
4. Figure 2 replots retained mechanochemical tables and trajectories in
   `polarity1d/results/data/figure2/`, with cue durations recovered from
   `polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response.json`.
5. Figure 3 uses the same five selected Gaussian-hierarchy fits as the original,
   with revised presentation distinguishing observations from fitted quantities.
6. `render_v2.py` typesets the pages with Matplotlib and atomically replaces
   `reports/trial_report.pdf`.

Generated figures are `figure1_size_mechanism`, `figure2_mechanical_test`, and
`figure3_evidence_levels` (PNG/PDF) under `reports/trial_report_assets/`.
Numerical checks are in `mechanism_statistics.json`.
`provenance_v2.json` was generated at 2026-09-26 23:29:04 UTC and records the
two SLURM jobs, selected fits, page count, and 25 source hashes. All 25 hashes
still match. Fresh plain `pdftotext` output matches `pdf_text_v2.txt` exactly.

The same session log records the primary submission at line 281, refinement
submission at line 375, final build at line 474, and final PDF verification
at line 485.

To rebuild the revised report from completed analyses in the repository:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/build_trial_report.py
```

This command replaces the current PDF and regenerates its figure/text/provenance
outputs. For verification we instead copied the wrapper and top-level assets
into a staging root and linked the retained data directories there.

## Upstream trails and remaining limits

- `polarity1d/REPRODUCE_RETAINED_RESULTS.md` describes the retained spatial-model
  workflow. `polarity1d/workflow/reproduce_retained_results.sh` orchestrates
  campaigns, collection, and plotting. Figure builders include
  `scripts/build_wp1_manuscript_figures.py` and
  `scripts/build_debelly_standalone_figure2_v2.py` within `polarity1d/`.
- `wgd_langevin/gaussian_hierarchy/README.md` documents the fit and validation
  campaign, including its SLURM scripts. The builders here consume completed
  fits; no model fitting was needed for PDF verification.
- Neither PDF nor builder has a matching tracked Git history in this checkout.
  The recovered original depends on the local session log, rather than Git.
- The surviving original provenance and exact reconstruction establish that
  `docs/trial_report.pdf` is the earlier generated report. The specific later
  operation that placed that file in `docs/` was not established.
- PDF assembly and figure generation were verified; the upstream scientific
  simulation and fitting campaigns were not rerun during this investigation.
