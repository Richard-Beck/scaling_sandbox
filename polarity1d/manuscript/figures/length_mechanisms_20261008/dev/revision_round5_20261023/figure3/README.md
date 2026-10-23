# Figure3: size-dependent boundary susceptibility and drift

[Main PDF](figure3_round5_review.pdf) · [supplement](figure3_round5_review_supplement.pdf) · [legend](legend.md) · [methods/results](METHODS_RESULTS_review.md) · [reader](index.html).

Cue placement separates sliding of the incumbent from takeover by a new region; enlargement slows both routes. Weak boundary susceptibility changes differently across models: terminal response at the largest size is 70% of the smallest-size response in Wave-pinning, 9% in Goryachev, 91% in Holmes and 99% in Spring.

The deterministic boundary response predicts physical domain drift without fitting when shared translation is retained. Wave-pinning’s C60/C20 drift ratio at f = 0.001 is 0.494 predicted and 0.473 measured. Holmes and Spring show much weaker susceptibility changes. The independent-edge approximation is refuted by boundary co-motion; Goryachev’s drift prediction remains quantitatively inconsistent and does not explain its full heading persistence.

No new native integrations, stochastic histories or parameter changes. The main displays the untuned shared-translation prediction and keeps the rejected independent-edge candidate. C insets use a shared normalized zero-based scale; D includes a fixed-input size comparison. The reviewed edition remains intact and is appended in the supplement.

Rebuild:

```sh
python -B polarity1d/manuscript/code/mechanisms/figure3_round5_review_20261023.py
python -B polarity1d/manuscript/code/mechanisms/figure3_round5_review_checks_20261023.py
```
