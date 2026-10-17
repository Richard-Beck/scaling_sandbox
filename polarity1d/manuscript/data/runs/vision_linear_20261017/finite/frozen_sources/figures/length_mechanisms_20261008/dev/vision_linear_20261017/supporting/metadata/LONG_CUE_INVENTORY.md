# Established linear long-cue evidence and design inventory

Measured 17 October 2026. This bounded inventory used existing files only; no integration or SLURM job was launched. Paths below are relative to `polarity1d/manuscript/`. The Figure 1 campaign owner subsequently relayed the root-approved fixed-support design below. The separate recovery pilot is not long-cue evidence.

## Approved finite exposure

Restart every eligible size from its exact corrected basal checkpoint. Apply the native right-end raised-cosine patch with **fixed support 3**, or **0.2 native length units for Spring**, and fixed within-model peaks Wave-pinning 0.1, Goryachev 6.4, deBelly 3.6, Holmes 1, Spring 4. The extra-input background is zero. Withdraw after **1000 native units**, or **50 for Spring**. Sensitivity branches at the existing three recovery-pilot sizes per model use 500/2000, or 25/100 for Spring. Preserve native coupled dynamics and deBelly's corrected `center_boundary_motion=False` independent endpoints.

The primary cue fixes local support and peak across sizes and matches the local-input gain geometry. It is a prolonged finite exposure, **not a common physical time calibration, an asymptotic stimulus, or a duration chosen to make every model reverse**. In particular, Holmes at reference length 20 needs a sampled 2360-unit pulse for the original scalar strength gate. The finite result may retain an incumbent, produce multiple active regions, lose polarity, or replace it. Show the full size series and both original gates rather than relabeling every response a pole replacement.

The newly approved withdrawal horizon is 13,000 cue-free native units, with spatial/scalar readouts at 3000, 8000, 10,000 and 13,000. Original archives end at 10,000. The final-500 window and its separate native-field stability diagnostic must be labeled by the horizon to which they apply.

## Exact available coverage and clocks

The corrected canonical campaign is `data/runs/debelly_corrected_three_test_20261008/`. Its other four models exactly preserve the older audited three-test evidence. All 97 checkpoint paths and SHA256 values in `data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv` were matched to canonical relative-preparation cases and rehashed successfully. All three response geometries share these checkpoints. The 13 ineligible preparations remain unavailable, rather than reversal failures.

| Model | Eligible sizes | Fixed-support reference size: orientation / strength clock | Successful fixed-support clock range: orientation / strength | Censored orientation / strength |
|---|---|---|---|---|
| Wave-pinning | 20, L=5–39 on retained grid | L19: 73.5 / 73.5 | 5.90636–3700 / 5.90636–4540 | L39 / L33,35,37,39 |
| Goryachev | 21, L=3–39 | L19: 377.5 / 377.5 | 33.7263–3540 / 33.7263–3540 | none / none |
| deBelly | 12, L=18,19,21,…,39 | L19: 306.5306 / 335.7908 | 280.2296–565.6122 / 287.8061–565.6122 | none / none |
| Holmes 3 | 23, L=10–50 on retained grid | L20: 440 / 2360 | 15.75–3800 / 15.75–2360 | 13 sizes L26–50 / 16 sizes L22–50 |
| Spring | 21, L=0.6–2.6 in steps 0.1 | L1: 0.310142 / 0.310142 | 0.152–5.81632 / 0.152–5.81632 | none / none |

These are shortest **sampled** passing finite pulses, not cue-on reorganization times. The search caps are 5120 for Wave-pinning/Goryachev/Holmes, 2000 for corrected deBelly and 50 for Spring. Nonmonotonic responses preclude inferring any trial from a shorter passing pulse.

All **97 primary scalar histories** already exist at exactly 1000/50 with 10,000-unit withdrawal. They retain pulse/recovery orientation and contrast, final native states, and pointwise final-window extrema/last samples; they do **not** retain complete cue-on or recovery spatial fields.

| Primary archived trial | Orientation passes | Strength passes | Separate field-tail passes |
|---|---:|---:|---:|
| Wave-pinning, 1000 | 16/20 | 14/20 | 19/20 |
| Goryachev, 1000 | 14/21 | 14/21 | 21/21 |
| deBelly, 1000 | 12/12 | 11/12 | 12/12 |
| Holmes, 1000 | 8/23 | 6/23 | 23/23 |
| Spring, 50 | 21/21 | 21/21 | 21/21 |

The 500-unit sensitivity is already scalar-sampled for all 64 Wave-pinning/Goryachev/Holmes sizes, but not for corrected deBelly. The 2000-unit trial exists for all 76 non-Spring sizes. Spring 25/100 are new durations. These all-size counts exceed the planned sensitivity anchors; exact anchor provenance should still be frozen in the new manifest.

All 97 old zero-cue full-field controls pass. Maximum scaled whole-history excursions by model are Wave-pinning 7.97×10⁻¹⁵, Goryachev 6.67×10⁻¹⁵, deBelly 7.3854×10⁻⁵, Holmes 1.28×10⁻¹⁴, Spring 1.77×10⁻¹². deBelly's maximum is L18. Retained pointwise extrema establish these finite-window baseline bounds; they do not supply missing time-resolved sham histories after 10,000.

## Native spatial histories and reuse limits

* `length_quality_debelly_traces_20261009/workers/000.npz`–`011.npz` contain corrected fixed-support uninterrupted cue histories to 2000 at every eligible deBelly size, output interval 1, native step 0.5 and all 11 restart fields. Slots 012–023 are relative-20% cues. The matching manifests retain checkpoint/cue identities and frozen source hashes. No 1000-unit-withdrawal full history is supplied by these uninterrupted records.
* `length_revision_20261009/workers/00.npz`–`05.npz` contain six exact relative-20% uninterrupted histories to 5120: Goryachev then Wave-pinning at L9,21,39, output interval 0.5. They do not implement the approved fixed-support exposure.
* `three_test_mechanisms_20261006/trajectories/` contains selected passing/failing pulses and controls at selected sizes. The other four models remain eligible after exact checkpoint/cue checks; **its centered deBelly histories are deprecated**. Holmes fixed-support/native cap examples cannot fill the all-size primary branches.
* The old `data/reversal_fields/*_full.npz` files belong to the preceding width campaign, with different checkpoint conventions and Goryachev peak 0.1. They are unsuitable substitutes for this campaign.

Fresh exact branches are the straightforward way to obtain all-size spatial histories and the longer withdrawal. Reusing a cue prefix requires an identical complete restart state and sampling/native integration convention; an interpolated profile is not an exact withdrawal checkpoint. No prior spatial archive establishes the approved all-97 cue-plus-13,000 protocol.

## Readout and resource guidance

Freeze the basal central background and incumbent half-contrast connected core before stimulation, and mirror that material-coordinate core for the challenger. Report raw primary-activity integrals and signed excess above the frozen background separately; state whether negative excess occurs instead of silently clipping it. Normalize both sites by the initial incumbent amount and retain absolute values. Report local peaks against the fixed initial half-contrast threshold, coherent widths/site count, incumbent first 50% loss and challenger activation clocks. Preserve unavailable crossings and full nonmonotone histories. Mark cue/core overlap, especially the fixed-width Goryachev L3 case. Spring core membership follows its material coordinate and its integrals use its own instantaneous physical length.

Cue-on and postwithdrawal fate are distinct observations. Use onset, declared cue-on times, exact withdrawal, early cue-free samples, and 3000/8000/10,000/13,000 endpoints. Retain signed global first moment and contrast as activity descriptors, alongside local profiles; neither defines velocity. Keep the original orientation, scalar strength and native-field tail flag separate.

Historical full-history SLURM array 20728794 had 235 completed workers: elapsed median 16 s, maximum 43 s, summed elapsed 3471 s. Corrected deBelly array 20769453 had 36 workers: median 9 s, maximum 21 s, sum 453 s. These support an initial planning allowance of approximately 0.5–2 summed worker-hours for 97 fresh primary branches; 30 anchor sensitivities add work. Reserve one CPU/3 GB under QoS `small`, with no concurrency throttle unless a concrete constraint arises. Actual usage must be measured; historical MaxRSS accounting is not a reliable Python peak estimate.

Dense early/withdrawal samples with modest later native output should produce roughly 2–4 GiB of uncompressed arrays, or approximately 0.8–2 GiB compressed for the primary 97 histories, depending on deBelly field duplication and sampling. Sensitivities and exports can take the group beyond 2 GiB. Measure allocated storage after execution and update the root manifest through its owner. This inventory itself creates only this small note.

The evaluated relative-20% alternative used twice the reference-size scalar strength clock, rounded upward to an existing sampled duration: WP104, G760, deB575.7142857, Holmes1040, Spring0.5954856. Those exact scalar trials exist at all 97 sizes and retain size resistance. That option changes local support/dose with size and was not selected as the primary protocol.

## Source identity

SHA256 measured 17 October 2026. Current `native.py`, `polarity_establishment.py`, and all five model implementations exactly match the corrected canonical campaign's frozen-source digests. Every selected checkpoint also matches the canonical digest. Freeze the new campaign's full executable closure, not only these entry points.

| Path | SHA256 |
|---|---|
| `data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv` | `653800607d4260b7474de98a04947afe9f8240eee41857eee50769387e81ef80` |
| `data/runs/debelly_corrected_three_test_20261008/manifest.json` | `0b7271bbb7b65e3107c86fee4b52da8c7c33bff9a1c6b93fe7d99fd6cbcdfdde` |
| `…/results/records.json.gz` | `8edd8b5d2d35ddbdc2f685cf6ae96b31af4f60bc1daee75c23b09db962b378f4` |
| `…/results/independent_audit.json` | `4c85e842d95f3e00fedd9b6d3e8e6b0d1fa466b8b718f0c71f2f46546e6d73a7` |
| `…/results/replacement_comparison_audit.json` | `c9e5fa8de1a31c89654f82c6fd49d6741f4f066c125cfec86e49b42b3aaa3dcd` |
| `code/core/native.py` | `f94a71a1bff6eba0be19d00a0f498bf6eb7994645cb85323c42c60bb7084c8a4` |
| `code/core/polarity_establishment.py` | `1e5834d0127c31e185af4de14d2bea1789081ac43ff9530e0ee357d2c310d429` |
| `code/core/models/wave_pinning.py` | `3167366db522833c7de3ac7b7c0285d23d27a9bc53409ccb0f0798618a60b54b` |
| `code/core/models/goryachev.py` | `3729ec29ed0861b47470b8acf5bc2e5c8d59c8cfc7940aa9377166296c547887` |
| `code/core/models/debelly.py` | `2c55b50e03a74abaaa55eeccbb4b18bb1b8594d14dc716a901c51fe678793f96` |
| `code/core/models/holmes_model3.py` | `c130081420e002c510a121811339f92576484b7287d339f971193c79ba95587e` |
| `code/core/models/spring.py` | `874871bfd60cb9cc94a8390625f0a57260a6578b9f310c69935dcac5b235b93f` |

