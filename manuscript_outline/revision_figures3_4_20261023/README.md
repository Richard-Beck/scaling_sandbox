# Figures 3–4: edge response, drift and secondary-domain dynamics

[Integrated outline_v11](../outline_v11.pdf) · [Figures, legends and pilot table](figures3_4_with_legends_and_pilots.pdf) · [Reader](index.html).

Starting point: the committed outline_v10 PDF and corrected standalone Figure 4, recorded by exact hashes in [the starting checkpoint](request/starting_checkpoint.json). The [revision brief](request/revision_brief_20261023.md) and [task prompt](request/user_prompt_20261023.md) are preserved verbatim.

[Delivery and identification limits](revision_acceptance.md) · [Point-by-point Claude responses](claude_response.md) · [Verbatim review hashes](review_preservation.json). Claude reviewed the starting checkpoint and revision checkpoint `49754ab`; the final deliverables incorporate responses to both reviews. The post-review changes have not received another Claude review.

Figures 1, 2, 5 and the reserved Figure 6 preserve their outline_v10 scientific content and artwork. The Figure 6 and bibliography vision citations are repaired from the deleted root MANUSCRIPT_VISION.md to the existing author/vision.md; no author text changes. No new noisy histories are generated. New deterministic edge measurements and derived compact-record analyses are documented with each figure; missing measures and failed assumptions remain explicit. Prospective observables are in [the pilot table](pilot_observables.json).

Rebuild using the existing scaling_sandbox environment after both figure producers finish:

```sh
python -B manuscript_outline/prepare_figures3_4_20261023.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v11.json --output manuscript_outline/outline_v11.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v11.json --report manuscript_outline/outline_v11.validation.json
python -B manuscript_outline/verify_figures3_4_20261023.py
```
