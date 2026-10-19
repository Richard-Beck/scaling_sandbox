# Compact figure reconstruction

The [exact source closure](supporting/metadata/compact_reconstruction.json) lists the current main plate and comprehensive supplement inputs with SHA-256 hashes. Both renderers use exported numerical readouts; no new integration or ignored native archive is needed to redraw these figures.

From `polarity1d/manuscript`, use the existing `scaling_sandbox` Python:

```bash
python code/mechanisms/revision_figure3_main_plate_20261019.py
python code/mechanisms/revision_arc_narrow_supplement_20261019.py
```

The main renderer reads seven native-cell, sector, moment, rate and field CSVs. The comprehensive renderer reads all 645 exact condition rows and final active profiles, selected native rates, all 129 longest-cue histories, and the twelve compact initial basal profiles. It produces four five-page core supplements and four separate fifteen-page native-profile appendices. The source tables and manifests retain their native lineage and finite-field-tail flags.

Both renderers resolve the same repository-relative exported sources after a checkout moves. The supplement matches the complete source suffix and requires its identical SHA-256, rejecting ambiguous source identities. Preserved absolute native archive paths remain lineage rather than redraw dependencies. Redrawing may change PDF metadata and generation timestamps; numerical sources and plotted values are the reproducibility target.

The larger native archives remain preserved for independent scientific replay. They are not substitutes for, or dependencies of, the compact plotted-data closure. The [scientific report](../arc_review/index.html), [main legend](main_plate/CAPTION.md) and [comprehensive supplementary legends](supplement/LEGENDS.md) provide interpretation and methods.
