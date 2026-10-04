# Minimal polarity figure codeset

Run commands from the repository root. This package reproduces the polarity
artifacts delivered in `current manuscript figs/`: `main_polarity.png`,
`supplement_width_response.png`, and the seven `model_*.png` establishment
heatmaps. It retains the completed 2026-10-03 five-model reversal comparison
and corrected seven-model basal-withdrawal establishment results.

The migration, MSD and Hosny figures are separate workstreams. The
`reversal_length_effect_mechanisms.pdf` explains an older linear-conditioned
campaign; it is a historical diagnostic, outside this current polarity codeset.

## Files

| Path | Purpose |
| --- | --- |
| `code/core/models/`, `code/core/src/polarity1d/` | Native model definitions, parameters, solvers, observables and integration |
| `code/core/native.py`, `settling.py` | Exact full-state restart, coordinate conventions and field convergence |
| `code/core/protocol.json`, `data/inputs/` | Model-length grids, fixed field scales and reversal duration grids |
| `code/core/polarity_establishment.py` | Rest → left cosine patch → basal withdrawal; duration selection and measurement |
| `code/core/corrected_reversal.py` | Exact selected checkpoint → mirrored right patch → basal recovery; both gates and local duration refinement |
| `code/audit_reversal.py` | Independent checkpoint, sample, gate, field flag and shortest-tested-duration audit |
| `code/figure_builder.py`, `code/core/polarity_establishment_plot.py` | Main/width figures, descriptive fits and establishment heatmaps |
| `code/build.py`, `code/check_codeset.py`, `code/tests/` | Rendering entry point, compact evidence verification and scientific tests |
| `data/establishment_basal/` | Full establishment manifest/table, 291 qualifying nominal-amplitude restart checkpoints and one preparation trace |
| `data/runs/corrected_reversal_20261003/results/` | Final readouts, recorded passing/failing example, status and completed audit evidence |
| `CODESET_SHA256.json` | Explicit shipped-file inventory and SHA-256 hashes |

`PROTOCOL.md` defines cues, clocks, thresholds and exclusions; `ASSUMPTIONS.md`
records native conventions and limitations. Five-model reversal excludes LEGI
and Otsuji. Field stability remains separate from scalar reversal success.

## Environment and quick reproduction

Check `conda env list` first; on the original HPC use `scaling_sandbox`. The
verified packaging environment is Python 3.12.14 with versions in `requirements.txt`.
On another machine, install those dependencies in an appropriate environment.
No compiled model extensions, R, Julia, external datasets or repository-root
Python installation are required. Runtime imports are internal to this folder.

```bash
conda activate scaling_sandbox
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
P=polarity1d/manuscript
python "$P/code/build.py" verify
python -m pytest -q "$P/code/tests"
python "$P/code/core/smoke.py"
python "$P/code/build.py" render
python "$P/code/core/polarity_establishment_plot.py" \
  --out "$P/data/establishment_basal" \
  --figures "$P/figures/supplement_establishment"
```

Open `figures/index.html` for generated PNG/PDF/SVG exports. Pages show the
source data-generation time and require manual refresh. Rendering may change
timestamps and PDF metadata; reproduction concerns the scientific content.
`verify` checks shipped hashes, exact checkpoints, all 660 readouts against the
recorded independent audit, and independently recomputes the 18 available fits.
It does not rerun the 73,975 native trial-array audits from compact summaries.

## Fresh simulations on SLURM

Use a fresh destination; the retained summaries have no worker chunks and
cannot be resumed. Initialization and smoke checks are small local operations.

```bash
python "$P/code/core/corrected_reversal.py" init \
  --out "$P/data/runs/reversal_example"
python "$P/code/core/corrected_reversal.py" smoke \
  --out "$P/data/runs/reversal_example"
python "$P/code/core/corrected_reversal.py" launch \
  --out "$P/data/runs/reversal_example"
python "$P/code/core/corrected_reversal.py" status \
  --out "$P/data/runs/reversal_example"
```

The runners generate `sbatch` commands: QoS `small`, one CPU and 3 GB per
worker, no array concurrency throttle. Reversal automatically schedules no-cue
controls, duration chunks, local refinement and the independent final audit.
Sources and inputs are hash-frozen; do not edit them while jobs run. Inspect
`launch.json` and `work/logs/` if a stage fails. Keep the final audit snapshot
stable: stop the periodic monitor after the last workers finish before the
final controller audits, as described in `REPRODUCE.md`.
To render a completed new campaign, pass `--campaign` to `code/figure_builder.py`.

The establishment supplement can also be recomputed from native rest:

```bash
python "$P/code/core/polarity_establishment.py" launch \
  --out "$P/data/runs/establishment_example"
```

That scans all 3,234 combinations and performs checkpoint, persistent-arrival
and endpoint-selection audits before final heatmaps. To feed its endpoints into
a new reversal campaign, work in a **separate fresh checkout**, wait for
establishment `results/status.json` to report complete, preserve its shipped
input directory, and copy the new campaign into the runner's declared input path:

```bash
mv "$P/data/establishment_basal" "$P/data/establishment_snapshot"
cp -a "$P/data/runs/establishment_example" "$P/data/establishment_basal"
# Then init/launch a new reversal destination as above.
```

Bulk chunks, all establishment trial histories, frozen duplicate source trees,
legacy campaigns, historical field explorers and rendered figures are omitted
from Git. Nothing is deleted locally. The supplied native input checkpoints
preserve exact restart identity; summaries alone cannot independently repeat
every historical trial audit. Fresh simulations regenerate that evidence.
