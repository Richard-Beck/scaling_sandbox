# Reproducing the corrected polarity core

Use the existing `scaling_sandbox` Conda environment. The complete manuscript
folder can be moved intact; runtime imports stay inside `code/core/`.
The minimal Git checkout is described in [CODESET.md](CODESET.md); optional
historical diagnostic data and exhaustive trial evidence are not shipped.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
conda activate scaling_sandbox
PY=python
$PY -m pytest -q polarity1d/manuscript/code/tests
$PY polarity1d/manuscript/code/core/smoke.py
# Initialize a fresh campaign first; the Git snapshot has no worker chunks.
$PY polarity1d/manuscript/code/core/corrected_reversal.py init \
  --out polarity1d/manuscript/data/runs/smoke_example
$PY polarity1d/manuscript/code/core/corrected_reversal.py smoke \
  --out polarity1d/manuscript/data/runs/smoke_example
```

The corrected matched reversal uses five models: wave-pinning, Goryachev,
deBelly, Holmes 3 and Spring. LEGI and Otsuji are excluded from this rerun. Each
condition starts from its exact qualifying corrected basal-establishment
checkpoint at nominal amplitude factor 1. The same added patch is reflected
to the opposite end, with the same peak and fixed/relative width. It is then
withdrawn to the same basal conditions.

There are 330 model-length-width conditions, 291 qualifying original preparations
and 39 original ineligible preparations. Every qualifying preparation first
receives an independent 10000-unit no-cue check; failures are reported as invalid
preparations and skip the pulse grid. The original checkpoint is never replaced.
Separate orientation and 1% scalar-strength gates use every final-500-unit sample.
Native-field tail variation is a separate flag. The initial pulse grid combines
the original dense duration grid, zero controls and 17 logarithmic durations
within the corrected bounds; near-equal durations are deduplicated. Every sampled
fail-to-pass transition is refined without assuming global monotonicity. Unresolved
numerics and execution errors stop finalization rather than becoming censoring.

The active runner is `code/core/corrected_reversal.py`. It uses SLURM QoS `small`,
one CPU/3 GB per worker, ten pulse durations per chunk, and no array concurrency
throttle. Source, experiment and input-checkpoint hashes are frozen. The campaign
lives at `data/runs/corrected_reversal_20261003/`. The SLURM monitor updates data and
its generation time during simulation; controllers also publish every stage.
For final audit, stop the periodic writer after the last worker array completes
and collect one snapshot before releasing the final controller, so the auditor
reads a stable evidence snapshot. HTML pages have no timed browser refresh.

```bash
# The retained campaign is complete. Use a fresh destination for a rerun.
$PY polarity1d/manuscript/code/core/corrected_reversal.py init \
  --out polarity1d/manuscript/data/runs/corrected_reversal_example
$PY polarity1d/manuscript/code/core/corrected_reversal.py launch \
  --out polarity1d/manuscript/data/runs/corrected_reversal_example
# Collect progress for that new campaign (updates its summary files).
$PY polarity1d/manuscript/code/core/corrected_reversal.py status \
  --out polarity1d/manuscript/data/runs/corrected_reversal_example
```

Compact trial NPZ files retain all scalar pulse/recovery samples, final native
states, and field-tail extrema/last-two-sample statistics. No-cue trials additionally
retain full-recovery field extrema and exact starting profiles for independent
preparation audits. Full spatial histories are reconstructed only for the selected
wave-pinning L15 fixed-20% passing/failing explanation example after finalization.
Figures are finalized only after complete execution and independent audits.

Regenerate the retained figures without launching simulations:

```bash
$PY polarity1d/manuscript/code/build.py render
$PY polarity1d/manuscript/code/build.py verify
# Optional, only in the full local archive with historical diagnostic data:
# $PY polarity1d/manuscript/code/build.py diagnostic
```

To regenerate the corrected establishment heatmaps:

```bash
$PY polarity1d/manuscript/code/core/polarity_establishment_plot.py \
  --out polarity1d/manuscript/data/establishment_basal \
  --figures polarity1d/manuscript/figures/supplement_establishment
```

A fresh establishment rerun requires an explicit destination:

```bash
$PY polarity1d/manuscript/code/core/polarity_establishment.py launch \
  --out polarity1d/manuscript/data/runs/establishment_example
```

`native.py` supplies native solvers and full-state checkpoints; `settling.py`
supplies native-field convergence diagnostics. Neither contains legacy assay
workflows. `protocol.json` holds 154 native model-length/grid settings, amplitudes
and verified-release criteria. Historical experimental descriptions are provenance.

`data/provenance/core_sources.json` records original and adapted source identities;
original corrected-establishment records and frozen sources retain their identities.
Selected corrected basal checkpoints and histories are local. Exhaustive nonselected
establishment histories remain outside the compact core; repeating their previous
full audits requires deprecated source data or a fresh campaign. The retained
interactive reversal-field report keeps its original experiment identity and
serves as a personal diagnostic.

Numerical completion and scalar polarity recovery alone establish neither a
length scaling law nor a causal mechanism.
