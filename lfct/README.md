# LFCT compact resegmentation

Entry points: [execution plan and measured smoke results](EXECUTION_PLAN.md), [original design](PLAN.md), [smoke viewer](reports/smoke/index.html), [portable viewer ZIP](reports/smoke.zip).

Run commands from the repository root. The existing `cpose` environment and cached CPSAM model are used; no environment installation is needed.

```bash
mkdir -p lfct/results/logs
OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/cpose/bin/python -m unittest discover -s lfct/tests -v
sbatch lfct/workflow/run_smoke.slurm
# After the job succeeds:
/home/4473331/.conda/envs/cpose/bin/python lfct/scripts/report.py --profile smoke --zip
```

The smoke processes 21 full-resolution phase/RFP pairs across four fields, including 12 consecutive MCF10A frames. The report keeps 12 full-field overview pairs and a 12-frame cropped negative-control clip. Open `index.html` locally, or unzip the portable bundle first. Keep the assets directory beside the HTML. Wheel/slider/arrow keys change frame; Ctrl-wheel or buttons zoom; drag pans. Layers include original XML contours, new cells/nuclei, tracked centroids, and unassigned masks. Original Track IDs, names and Spot IDs remain inspectable.

## Full run

The user accepted HEK293 segmentation limitations and visually checked U87 alignment. Full processing was submitted as field array **20236390** (25 tasks) followed by CPU audit/report job **20236391**, with `afterany:20236390`. The finalizer checks every completed field and all output/PNG hashes; failed or missing fields are listed in `results/full/final_audit.json` and prevent final-report publication.

One-shot entry point: `python3 lfct/workflow/submit_full.py` (`--dry-run` previews commands). Submission IDs are saved in `results/full/submission.json`; the launcher refuses accidental duplicate submissions. Each field requests one A30, four CPUs, 16G and one hour. No array throttle; SLURM schedules available GPUs. The finalizer requests two CPUs, 4G and 15 minutes without a GPU.

Monitor: `squeue -j 20236390,20236391`. Final outputs: `reports/full/index.html`, `reports/full.zip`, and `results/full/final_audit.json`. If fields fail, resolve their recorded issue and resubmit those array indices, then submit a new finalizer dependent on that retry array. Successful fields can validate and exit without inference.

## Durable files

- `config.json`, `scripts/`, `workflow/`, `tests/`: configuration and implementation.
- `data/sources.csv`: 25 original XML hashes and public image URLs. Regenerate deliberately with `scripts/prepare.py --discover`; this accesses the network.
- `data/annotations/`: lossless compressed original XML.
- `results/<run>/fields/<cell>/<field>/`: complete cell/nuclear label stacks, compressed detection/mask/run tables, frame QC, provenance, timings and checksums. TIFF page descriptions contain original XML frame numbers; labels are frame-local, not Track IDs.
- `reports/<run>/`: small HTML, PNG assets and QC selection manifest; `<run>.zip` is a rebuildable delivery copy.

Generated data/results/reports are ignored by Git. Back up these outputs separately. Raw ZIPs and selected source TIFFs live in an invocation-owned scratch directory, removed after successful validation. Existing MCF10A source TIFFs are reused without modification. On failure, `failure.json` records the error; node scratch may disappear after the allocation ends. A stale `.running` lock after an interrupted job requires checking that no writer remains before removing the lock.

A repeat with identical code/settings checks output and PNG hashes and exits without downloads or inference. Changes to analysis can reuse verified segmentation; changes to segmentation/source mapping require a new run name (existing masks are protected). Use `--run-name NAME` for a separate experiment and `LFCT_CONFIG=/absolute/config.json` for alternate settings. Do not edit an active run's configuration. For browser verification, use an existing Playwright installation, e.g. `NODE_PATH=/tmp/lfct_browser/node_modules node lfct/tests/browser_smoke.cjs` (the temporary module location may need updating in future sessions). The corresponding report command accepts `--run-name NAME` too.

## Scientific interpretation

CPSAM runs independently on phase (cells) and RFP (nuclei). Missing, ambiguous and shared-mask assignments remain in the tables. No original IDs are renumbered. Isolation requires a containing unique cell assignment, at least 10 consecutive nonbranching frames, at least 10 px clearance from every other cell mask and the field edge, and no nearby unmatched centroid/orphan nucleus within the provisional 50 px uncertainty guard. Missing nuclei alone do not fail a cell. Splits, merges, gaps and invalid observations break runs. Thresholds use pixels/frames because physical calibration has not been established.

Masks are the exact durable scientific representation; viewer contours are derived display geometry. PNGs are contrast-normalized 8-bit QC images, with bounds recorded in `qc.json`, and cannot replace quantitative original fluorescence data. Candidate isolation is only established at sampled image times and remains subject to segmentation/QC errors.

U87 FLD_1 has an explicit **provisional** source mapping exception: omit numeric source suffix 1 from both 289-image archives to match 288 XML frames. XML 0 maps to suffix 0; XML f≥1 maps to f+1. RFP suffix 1 is anomalous, and visual checks at XML 0,1,144,287 support this mapping. The user subsequently checked and accepted the alignment; see `data/qc_decisions.json`. Other unexpected index mismatches fail rather than guessing. See `results/smoke/alignment/` for evidence.
