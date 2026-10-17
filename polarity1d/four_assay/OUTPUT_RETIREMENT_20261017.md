# Deprecated output retirement — 17 October 2026

The user explicitly requested removal of large outputs associated with this
deprecated four-assay tree. Execution status, measured allocation, exact output
roots, archive checks and protected-file checks are recorded in
[`reports/storage_cleanup_four_assay_20261017.json`](../../reports/storage_cleanup_four_assay_20261017.json).
The [initial storage audit](../../reports/storage_audit_20261017.md) preceded this
broader authorization; its six-root NPZ-only plan is not the complete execution
scope.

The cleanup retires pilot and revision-11 worker payloads, older exploratory
continuation arrays and duplicate historical native-viewer arrays. It preserves
compact scalar results, endpoint CSVs, figures, protocols, source code, frozen
inputs and provenance. Nonempty worker metadata and logs are losslessly
consolidated into independently verified `.tar.gz` files in
`retained_metadata_20261017/`; paths inside these archives are relative to this
four-assay directory. Empty worker files carry no content and are omitted.

The accepted v8 `baseline/`, final revision-11
`reference/chained_v11_20261001/` profile/provenance bundle, and saved
`results/assay_explanation_fields/` remain. Five native inputs retain executable
support for the legacy wave-pinning assay explanation:

- `work/checkpoints/task_0024.npz`
- `work/checkpoints/task_0025.npz`
- `work/results/trace_0024.npz`
- `work/histories/conditioning_0024.npz`
- `work/histories/release_0025.npz`

The approved historical `assay_explanation_v2` keeper and its sources remain
outside this tree, under `../archive/four_test_suite_legacy/`. No manuscript data,
current report native arrays, canonical model sources or new circular experiments
are removed.

Before cleanup, all 4,952 canonical selected establishment NPZ files and 54 copied
diagnostic input NPZ files passed their recorded SHA-256 checks. Data-only
regeneration of the current corrected main and width-response supplement from
660 endpoint rows reproduced both canonical PNGs byte for byte, without any
integration or mutation of manuscript inputs.

Historical raw coverage is deliberately reduced. Old worker checkpoint/trace
paths in scalar records are provenance identities, not promises of locally
retained arrays. The six-model symlink reuse in the historical corrected-basal
pilot also loses its original raw targets. Canonical selected establishment
states and selected personal-diagnostic inputs survive independently in
`../manuscript/data/`. Full re-audit or arbitrary retrospective analysis of the
retired native trials requires rerunning their frozen scientific inputs. Neither
archived scalar records nor a regeneration recipe substitutes for those arrays.
