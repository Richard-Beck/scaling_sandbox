# Final round-3 supporting package

The current manuscript is [outline_v9.pdf](../outline_v9.pdf), with its editable JSON and exact-input build manifest beside it. It includes the final 22 October Figure 2 revisions, the final Figure 3 single-C revision, and the completed Figure 4 campaign.

Final main plates, legends, supplementary PDFs, and the outline's direct build inputs remain ordinary repository files. supporting_data.zip holds the remaining compact plotted-data tables, bounded illustration inputs, and scientific/provenance receipts for those three final figure versions. It preserves repository-relative paths and original bytes. manifest.json lists every member's size and SHA-256. Full campaign histories, superseded outline PDFs/layouts, and verbose data-retirement ledgers remain on HPC.

To restore the additional supplementary files, from the repository root run:

    python -m zipfile -e manuscript_outline/supplement_round3_20261023/supporting_data.zip .

The extraction restores the original supporting paths used by the figure readers and scientific verification scripts; existing outline build inputs are not overwritten. The archive contains only new supplementary paths. Scientific and integration verification passed on the complete HPC working tree before packaging. On a fresh checkout, extract this archive before running the full scientific verifier.

Recheck the current outline with:

    python scripts/verify_manuscript_outline.py manuscript_outline/outline_v9.json --require-tracked
    python manuscript_outline/verify_figure4_round3_20261022.py

The second verifier also checks the retained Figure 2/3 sections against outline_v8.json. This small editable baseline is preserved; superseded outline PDFs are omitted.
