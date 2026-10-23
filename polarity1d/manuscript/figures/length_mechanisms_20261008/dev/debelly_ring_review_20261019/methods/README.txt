Retained deBelly viewer
======================

Open ../viewer/index.html in a browser. The viewer works directly from local
files: its data scripts use relative paths and require no server or network.
The executive review is ../index.html.

Three retained branches show the original 1,000-s challenger assay, the
existing-face gap variant with the original cue at 20 µm (A, the default), and
the existing-face variant with the cue at 20.061728395061703 µm (B).

Each branch retains all seven fields, native widths, gap identities, positions,
all four flank controllers and next_id at every native 0.5-s state in the event
windows. Same-clock pre-transition and post-transition frames are separate.
The next frame contains the following complete native half-step. There are
18 event clocks across the three branches, comprising 30 opening/closing
records. Each branch also retains 125 exact full-precision native states at
250-s cadence over the full 0–31,000-s observation. The full original histories
are not needed to use or regenerate this viewer.

Regenerate from the retained local data only:

    python methods/render_viewer.py

Run that command from the package root with Python and NumPy installed. The
script also works when invoked by its absolute path from another directory.
It reads only viewer/raw and methods/viewer_ui.py within this package and
writes viewer/index.html, viewer/data/*.js and viewer/report_manifest.json.
No simulation or replay is run. Browser payloads use 10 significant digits;
the downloadable NPZ archives retain the original full precision.

viewer/raw/manifest.json checksums all runtime raw inputs. The generated
viewer/report_manifest.json records local source/input/output checksums.
retained_provenance.json maps retained files to their original identifiers
and hashes. Files under provenance are historical validation receipts, not
runtime dependencies. Repository paths inside those receipts document where
the evidence originated; they are not links or instructions to read files
outside this package. Original replay validation passed before consolidation.
The separate packaging audit verifies local retention and relocation.

Optional browser-control check (Node with built-in WebSocket and Chromium):

    CHROME_BIN=/path/to/chromium node methods/check_viewer.mjs viewer/index.html

This tests all 21 segments, all seven fields, the same-clock stage buttons,
time stepping, reference/scale controls and desktop/mobile layouts. Its local
screenshots and recorded result are retained in methods/browser_checks.
relocation_check.json records a second render from a temporary copied package
while access to the original repository was explicitly forbidden. All 21 data
scripts were byte-identical after relocation and to the pre-cleanup viewer.

No scientific law, parameter, event rule or stored state was changed during
packaging. Development and experiments have stopped.
