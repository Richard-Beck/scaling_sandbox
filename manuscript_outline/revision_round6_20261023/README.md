# Round6 manuscript revision

[Integrated outline_v12](../outline_v12.pdf) · [Figures 1–4 and legends](figures1_4_with_legends.pdf) · [Reader](index.html).

The starting point is the current outline_v11 PDF. Its exact hash, the editable
outline, the verbatim brief, the supplied narrative and the task prompt are
archived in [request/](request/). Round5 remains intact. Figures 1–2 change layout
and lettering only; Figures 3–4 report the requested new analyses. Figure 5 and
reserved Figure 6 are preserved. The prospective pilot directions are retained;
panel references and model-conditioning notes match the revised figures.

Each figure folder exposes one rebuild entry point and a `scientific_readme.md`
beginning with Unexpected. Main figures, compact plotted values and complete
supplements live beside round5 in `revision_round6_20261023/figure1/` through
`figure4/`. Scientific qualifications and campaign completion are recorded there.

After those figure rebuilds, integrate using the existing environment:

```sh
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/prepare_round6_20261023.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v12.json --output manuscript_outline/outline_v12.pdf
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v12.json --report manuscript_outline/outline_v12.validation.json
```

At the end of the native campaign, refresh allocated storage and require the
terminal allocation audit:

```sh
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/measure_round6_storage_20261023.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/verify_round6_20261023.py --require-complete-mode-campaign
```

[Claude review 1, preserved verbatim](../../reviews/claude/20261023_manuscript_round6.md) · [Integrity record](review_preservation.json) · [Point-by-point response](claude_response.md)
