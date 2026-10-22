# Completed round-2 manuscript Figures 2–4

[Integrated outline_v5](../outline_v5.pdf) · [Figures with legends](revised_figures_2_4.pdf) · [Three main plates](figures_2_4.pdf) · [Reader and supplements](completion.html) · [Acceptance evidence](completion_acceptance.json).

Figures 2–4 are complete, including the 32 Holmes train/noise histories, all nine approved Figure 3 extensions, the corrected short-range interpretation and the Figure 4 Case A main. Six Figure 3 extensions resolve; three remain censored at 80,000. Holmes class membership is unchanged. Numerical and finite-horizon qualifications remain explicit. Figure 4 H1 is explicitly **post hoc**, with frozen calls alongside; geometry and multiple-front fractions are descriptive. Figures 1, 5 and 6 retain their exact preceding content and artwork.

Enlarged cells resist imposed replacement through distinct initiation, spatial and temporal response signatures (Figure 2). Distance between competing regions matters only at short range; enlarged cells resist replacement and resolve coexistence slowly mainly because of size itself (Figure 3). Spring's middle size is the supported exception in which extra accessible separation adds a substantial resolution-time contribution. Enlarged cells hold heading longer under fluctuating input — at or above the geometric expectation at weak input — and the advantage declines as input strengthens. (Figure 4). Geometry and multiple-front fractions are descriptive: erosion below the single-front geometric expectation can also occur while a single front persists. Observed migration summaries, conditional OU velocity-memory estimates and simulated polarity-heading readouts remain separate evidence.

The scientific producers have stopped and every required trajectory is accounted for. Original evidence and complete supplements are retained; raw native arrays and the whole-seed bootstrap archive remain on HPC. The previous interim and in-progress review remain historical artifacts.

Rebuild with the existing scaling_sandbox environment:

```sh
python manuscript_outline/prepare_round2_completion_20261022.py
python scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v5.json --output manuscript_outline/outline_v5.pdf
python scripts/verify_manuscript_outline.py manuscript_outline/outline_v5.json --report manuscript_outline/outline_v5.validation.json
python manuscript_outline/verify_round2_completion_20261022.py
```

Generated UTC: 2026-10-22T13:51:12.001470+00:00. Commit and remote verification are recorded in completion_delivery.json after publication.
