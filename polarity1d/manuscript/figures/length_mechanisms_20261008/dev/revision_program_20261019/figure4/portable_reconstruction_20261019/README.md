# Figure4 exact redraw closure

This committed compact bundle redraws the final one-page main, matching native-profile sheet and the complete32-page experiment supplement without opening native simulation archives. Full scientific replay remains a separate operation using preserved local native archives.

The export contains112 exact field matrices:108 declared condition medoids plus the four main matched-mean controls. Each matrix contains the exact basal row and every native second-half profile, native coordinates/circumference and scalar context at those same profile times. The112 CSV.gz chunks retain30,103,368 binary64 values with zero round-trip error (`%.17g` serialization). The largest compressed chunk is below5MB. No field normalization or spatial/time decimation is introduced.

All case/condition/full-half/quarter/CI tables and representative-selection inputs are compact CSV/JSON. Original absolute native paths and hashes in `FIELD_INDEX.json` and selection metadata are lineage, not redraw dependencies. The redraw process actively blocks opens beneath native data/runs paths. It rebuilds temporary NPZ plot intermediates from committed CSV.gz fields, and needs no ignored native NPZ.

From a fresh repository checkout with NumPy, pandas, Matplotlib and pypdf available:

```bash
python polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/portable_reconstruction_20261019/helpers/redraw.py --out /tmp/figure4_redraw --recompute-ci
```

`--recompute-ci` also reproduces the2000 joint whole-seed intervals from exported per-case sufficient statistics and checks the three CI CSVs byte-for-byte. Omit that flag for redraw from the stored exact intervals. No solver/native adapter is imported and no new integration runs.

The verification reproduced34 PNGs byte-for-byte: the compact main, matching four-model profile-support sheet, all12 native atlas pages,9 temporal pages,9 organization pages and2 late-quarter pages. The comprehensive32-page PDF is reconstructed from those source plots. The earlier expanded signature matrix/two-sheet candidate remains preserved in the source assembly but is not part of the compact redraw target. PDF creation timestamps may differ; figure values and PNG rasters are exact.

[Export and round-trip receipt](metadata/EXPORT_RECEIPT.json) · [Verified redraw receipt](metadata/PORTABLE_REDRAW_CHECK.json) · [Fields and native lineage](metadata/FIELD_INDEX.json) · [Complete reconstruction closure](metadata/RECONSTRUCTION_CLOSURE.json)

Replaying all936 native profile moments, region geometry or numerical/forcing checks requires the local raw archives. This bundle retains their exact identities and all plotted outcomes; it does not replace the raw scientific audit. Spring’s original failed sampled-duration receipt and the separately accepted ten-metric scope remain unchanged.
