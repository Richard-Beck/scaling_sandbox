# Independent review of the bounded Spring early-window amendment

Design reviewed 17 October 2026 before initialization or launch. Actual outputs independently reviewed later that day, after completion of SLURM 21126520.

**Design verdict: APPROVE launch of the two declared Spring cases only. No essential launch blocker remains.** This approval covers the computational amendment, not acceptance of its unmeasured result. Independent review of actual outputs remains required before the parent campaign substitutes its target prefix or supersedes the original live task.

**Actual-result verdict before the original task finished: ACCEPT the L₀ = 0.8 prefix for substitution at the unchanged reviewed clock. Both cases qualify and all 636 neighbor comparisons pass; no essential result blocker remains.** The independently checked evidence below satisfied the design's result-audit requirement and authorized the prospective supersession/cancellation route while preserving original records. Later target response remained unmeasured at that moment.

**Final handoff status: the original L₀ = 0.8 task completed before any cancellation or substitution. The main 97-case collection uses its complete original archive; all original later candidates are measured there.** The amendment remains separate accepted end-time sensitivity evidence, and only its own later candidates are unmeasured. The final update below supersedes the prospective substitution/cancellation route discussed in the earlier review.

## Scope and scientific identity

I read the root and manuscript agent instructions, the active manuscript protocol and assumptions, the current status, the amendment, its algebraic preflight, the new wrapper, the shared assay, and the original frozen campaign definitions. I ran only read-only source, hash and metadata checks. I did not integrate trajectories, initialize the amendment, submit jobs, modify scientific sources or alter the original task.

The target is the exact original Spring L₀ = 0.8 case (parent slot 78); the prospectively declared control is adjacent completed L₀ = 0.9 (parent slot 79). Both retain the fixed reviewed clock τ = 0.01778279410038923. The wrapper retains all 212 original sample times from zero through τ, including the exact endpoint. It retains all 18 original candidate clocks at or before τ and marks the remaining seven later candidates unmeasured. It does not choose a new clock from the target response.

The equations, retained full-state checkpoints, frozen initial material core and weights, perturbation direction, active-to-inactive transfer and raw field treatment are unchanged. The common worker still executes ten independent trajectories per case: paired sham and ε = 0, 0.001, 0.005, 0.01 at each of two tolerance settings. The independent zero trajectory is checked against the sham in every archived field. Tight LSODA settings remain rtol = 10⁻¹⁰, atol = 10⁻¹², maximum step 0.01; tighter settings remain rtol = 10⁻¹², atol = 10⁻¹⁴, maximum step 0.005. The sole integration change is the final native time, from 1 to τ.

The smaller final time can change adaptive step selection. The declared neighbor comparison therefore tests a real numerical consequence of the amendment. Its allowance uses the original empirical sensitivity definition and the same factor of three already declared for amplitude agreement:

$$
\left|F_{\mathrm{new}}(\varepsilon,t)-F_{\mathrm{old}}(\varepsilon,t)\right|
\leq 3\left[E_{\mathrm{new}}(\varepsilon,t)+E_{\mathrm{old}}(\varepsilon,t)\right].
\tag{1}
$$

The collector evaluates this at all 212 retained times, including zero, for each of the three nonzero amplitudes: 636 comparisons. It asserts that every comparison passes. This is empirical numerical agreement at the declared sampling, not a certified error bound or a spatial convergence test.

Independently, both target and neighbor must pass the original selected-clock gates over every nonzero sampled prefix time: E ≤ 10⁻⁵; pairwise amplitude agreement within 5% plus the original 3(E₁ + E₂) allowance; and |F| + E ≤ 0.2. The existing evaluator enforces these unchanged gates. Near-null and negative values remain eligible, with resolution above 3E separately recorded. A neighbor comparison pass alone cannot qualify the target.

## Execution and provenance checks

The original manifest and all its frozen-source hashes passed verification; the present shared assay and native code have no hash differences from the original campaign's frozen code. Both target input hashes match the original cases. The completed neighbor archive matches its original status digest.

The initializer requires a fresh destination, copies the two exact input checkpoints and both harnesses, copies native Python/configuration dependencies, and records their hashes alongside the parent manifest/provenance and amendment. The frozen wrapper imports its sibling frozen assay, whose ROOT resolves inside the frozen source tree; native imports therefore resolve to the copied core. The worker verifies the manifest, frozen dependencies and copied checkpoints before executing. SLURM QoS `small`, one CPU, one BLAS thread, 3 GB, 20 minutes and an unthrottled two-element array match this bounded computation. Missing or timed-out execution must remain an execution failure.

The following prelaunch anchors allow the later independent result review to verify that collection compared against the same original neighbor, even though the collector reads its status and archive from the preserved parent folder:

| Record | SHA-256 |
| --- | --- |
| Amendment wrapper | `4ad48317635727b55960f1cae86bae847500573430eb08cf05cccf9bcee71988` |
| Shared assay | `78f472edbb2fffa65e8d0d24de2df0549b8f237d625941f90706b8d8feb3ad14` |
| Original parent manifest | `91950b1f882e491483509f6440e1f755dddc81553383bf0eae22cd38cf969f97` |
| Original parent provenance | `616c8487dd936d6f8cdda2bb2b9d91a11b683a88abed75ab5b16387becad19f0` |
| Original L₀ = 0.9 status (`workers/79/status.json`) | `07e2e81933668ba90ce4b56d3c0619957e2bd3812dd9a86f9297b875d04e5b00` |
| Original L₀ = 0.9 archive (`workers/79/curves.npz`) | `b900eb10a1f1df4948fca9638f42dbdccf275946c8299424eca7db8f20be7aef` |

## Required result handoff

After execution, independently verify the frozen closure and exact input hashes, all twenty complete trajectories and 212-point time axes, full-field zero identity, pool and raw-field checks, both selected-clock gates, and every neighbor comparison from archived amounts. Confirm the original-neighbor anchors above. Record maximum differences and margins, rather than relying only on the collector's Boolean result.

Only after that audit passes may the parent accept the L₀ = 0.8 prefix and explicitly supersede/cancel its unused t = 1 task if still live. Preserve original input, source, status, logs, launch records and any original output. The root table must retain an explicit computational-span/provenance field and represent later L₀ = 0.8 candidate responses as **unmeasured**, never extrapolated, failed or censored. The original task's long runtime is not evidence of numerical failure. This amendment supplies no new scientific interpretation beyond the unchanged local-core finite-window recovery readout.

## Independent actual-output audit

The audited run is `data/runs/vision_linear_20261017/recovery/early_window_amendment`. Scheduler accounting reports task 21126520_0 completed with exit 0:0 in 9 s and task 21126520_1 completed with exit 0:0 in 7 s, each with one CPU, 3 GB and QoS `small`. At the audit, the original 21126305_78 task was still `RUNNING`, with its source/input/status/log records preserved. I submitted no job, performed no native integration and did not cancel that task.

I independently verified the amendment manifest/provenance, all 32 frozen-source digests, both input checkpoint digests, the reviewed wrapper/shared-assay anchors and the original parent/neighbor anchors above. The new native and shared code hashes match the original frozen campaign. Both worker status records are complete and their archive hashes match:

| Case | Archive SHA-256 |
| --- | --- |
| L₀ = 0.8 | `73e249a3e67fed96b65461983f7b33215c1af5a468aee1fa2e5c807586b4db9c` |
| L₀ = 0.9 | `0842994cf31b6e3e7fa263c7258a32c7eccb084cdb4f9f0ca64face8bd8da86d` |

All twenty declared trajectories contain exactly the original 212-point time axis, from zero to τ inclusive. Their five native fields are finite, physical lengths stay positive, and length equals the stored endpoint separation. Independent geometry reconstruction reproduces both frozen material cores, perturbation directions and integration weights. Archived initial chemical samples agree with the exact prescribed solver initial states to at most 2.22 × 10⁻¹⁶, consistent with floating-point evaluation of the LSODA output. At each tolerance, independent zero and sham agree **exactly in all five fields** at all times. No clipping is reported or needed: the smallest raw chemical concentration is 0.7084654381128734.

I reconstructed core and whole active amounts from every archived field, matching the archived amounts within 1.11 × 10⁻¹⁵. The independently reconstructed total active-plus-inactive pool has maximum relative drift 2.78 × 10⁻¹⁵, far below the unchanged 10⁻⁸ check. Its tiny difference from the originally reported reduction, at most 1.33 × 10⁻¹⁵, is floating-point reduction order. These raw-field checks precede the response calculation.

Using the archived amounts, I independently recalculated F, the full original E definition, all 36 candidate-prefix gate rows, all 1,272 short-curve rows, and all 636 neighbor-comparison rows. They match their exported values. Both selected-clock results pass the numerical, amplitude and saturation gates over every nonzero prefix sample and every declared nonzero amplitude:

| Quantity at the unchanged τ | Target L₀ = 0.8 | Neighbor L₀ = 0.9 |
| --- | ---: | ---: |
| Primary restored fraction F (ε = 0.005) | 0.06466203463557607 | 0.07561770096413589 |
| Conditional F/τ, native time⁻¹ | 3.63621342464741 | 4.252295816801971 |
| Primary empirical sensitivity E | 9.261906797064512 × 10⁻¹⁰ | 1.6484744680411723 × 10⁻⁹ |
| Maximum prefix E, all amplitudes | 7.543404567655898 × 10⁻⁹ | 9.830903202612262 × 10⁻⁹ |
| Maximum prefix absolute F plus E, all amplitudes | 0.06487520428737914 | 0.07584108810744083 |
| Maximum amplitude difference minus declared allowance | −3.6967312026528 × 10⁻⁸ | −3.913499882000052 × 10⁻⁸ |

Every neighbor comparison satisfies equation (1). The largest absolute F difference is 4.249639085074364 × 10⁻⁹. The largest difference divided by its allowed 3(E_new + E_old) is 0.16757148508445543; it occurs at ε = 0.001 and native time 0.00530884444230989, with difference 1.982451314219702 × 10⁻⁹ and allowance 1.1830481261299042 × 10⁻⁸. The smallest absolute allowance-minus-difference margin across all 636 comparisons is positive, 1.6792856654207753 × 10⁻¹⁰. Thus qualification does not rest only on collector status flags.

The amendment outputs contain no time or candidate gate beyond τ. Its manifest and verification explicitly list the seven later candidates as unmeasured in that bounded run; neither their F nor their gate outcome is inferred from its prefix. The then-planned substitution route would have required the parent to preserve that status and the changed computational-span provenance. That route was not used, as the final update records. This acceptance concerns the fixed local-core finite-window measurement only.

The amendment occupies 9,412,608 allocated bytes at audit. This review created only this small Markdown record; storage accounting for the combined parent output group remains the producer/root manifest owner's responsibility at handoff.

## Final handoff update: original target completed without substitution

Before any cancellation, original task 21126305_78 completed naturally with scheduler state `COMPLETED`, exit 0:0 and elapsed time 11 min 43 s. No cancellation or supersession occurred. Its final status and archive hashes verify:

| Original L₀ = 0.8 record | SHA-256 |
| --- | --- |
| `workers/78/status.json` | `9a885b3ff0c935317ac172cf7b37ef6eb0250bca95def16e1ac19ffaebf765ec` |
| `workers/78/curves.npz` | `b69d9babbb867ae1a783b7a872986ba967982ebee29cc4e639930645cd293e01` |

The original archive contains all 282 original samples through native time 1, including every later candidate clock. The main `recovery_summary.csv` points to this original archive and hash, with `required_window_amendment=False` and `later_candidates_unmeasured=False`. Its selected-clock result is F = 0.06466203407446192, conditional F/τ = 3.6362133930936418 native time⁻¹, and E = 6.735357382936513 × 10⁻⁹. These original values supply the main result; the amendment values in the preceding table remain separate sensitivity results.

As an additional retrospective diagnostic after the original target completed, I independently compared all 212 shared times at all three nonzero amplitudes between its original archive and the amendment, using exactly equation (1) without changing the allowance. **All 636 comparisons pass.** The largest absolute F difference is 2.306086765990756 × 10⁻⁹, the maximum fraction of the empirical allowance is 0.07094342904944452, and the minimum absolute allowance-minus-difference margin is 1.726599800511571 × 10⁻¹⁰. At the selected clock and primary amplitude, the absolute F difference is 5.611141551398191 × 10⁻¹⁰. This target comparison corroborates the predeclared neighbor check; it was not used to select that neighbor or retune any gate.

The final accepted evidence therefore comprises the complete original main-grid measurement plus the separately retained two-case end-time sensitivity experiment. Later responses are absent only from the bounded amendment, not from the original target or the completed main 97-case collection. No new jobs or integrations were run for this update.
