# Portable rendering of the supplied mechanism figures

`scripts/render_manuscript_mechanisms.py` rebuilds the three supplied mechanism
figures, including the two continued ring sheets, from the first published
checkpoint in Git commit `6199cdc`. It uses the declared compact inputs, the
five selected linear histories, four selected noise profiles and lossless
summary tables. It does not rerun simulations or the 243-worker statistical
analysis.

From a Git checkout containing that commit, use the existing
`scaling_sandbox` environment:

```bash
conda run -n scaling_sandbox python scripts/render_manuscript_mechanisms.py \
  --revision 6199cdc --output manuscript_outline/mechanisms_rebuilt
```

The output directory must be new. The script refuses to overwrite an existing
path, so the first outline, its JSON and supplied figures are preserved. The
temporary input checkout is removed after rendering; `--scratch-parent` can
place it on another suitable filesystem.

The wrapper obtains inputs from Git rather than the current worktree. It
verifies every original raw-source SHA-256, verifies compressed identities
separately, and expands only the three oversized CSVs needed by the original
renderer. It imports the original renderer bytes unchanged and invokes the
same function as its `--full` branch. A runtime file-lookup adapter resolves
the recorded historical repository prefix inside the temporary checkout.
The child process blocks attempts to open files in the original worktree.
No bulk native archive is available to that child.

The destination includes `portable_render_provenance.json`, with source and
gzip identities, environment versions, output comparisons and any relocated
metadata. The four original PNGs must match the supplied PNG hashes exactly.
PDF dates and SVG metadata or identifiers can change. The unchanged renderer
also records five temporary `source_archive_path` values in the longest-example
CSV; the wrapper verifies that replacing only that repository prefix recovers
the original CSV hash. Numeric values and all serialized scientific arrays
remain unchanged.

For the manuscript's final figure numbering, a separate fresh destination can
change only the four displayed figure-title labels: working Figures 1, 2 and 3
become manuscript Figures 2, 3 and 4.

```bash
conda run -n scaling_sandbox python scripts/render_manuscript_mechanisms.py \
  --revision 6199cdc --manuscript-numbering \
  --output manuscript_outline/mechanisms_manuscript_numbering
```

This mode preserves the renderer and scientific functions. Its provenance
records the four text replacements, original source-archive hashes and output
data hashes. A pixel comparison requires every pixel outside the top title
band to match the supplied PNGs exactly. This changes labels without changing
the scientific panels.

The isolated original and numbered validation runs are recorded in
`portable_mechanisms_validation_v2.json`. Full campaign reruns remain a
separate SLURM/archive workflow; this command reproduces the supplied rendered
evidence from its committed input closure.
