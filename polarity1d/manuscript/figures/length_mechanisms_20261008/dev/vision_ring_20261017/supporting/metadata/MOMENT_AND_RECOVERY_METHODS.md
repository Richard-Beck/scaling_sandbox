# Activity-moment and active-loss definitions

The global extraction imports the reviewed `moment()` and `analytic_checks()` unchanged from `global_polarity_vector_20261017.py`, with a preserved copy in the new extraction run. For primary active field values at native periodic cell centres,

$$
M_n(t)=h\sum_i u_i(t)\exp(\mathrm{i}n\theta_i),\qquad
\theta_i=\frac{2\pi(x_i-x_0)}{C},\qquad h=\frac{C}{N}.
\tag{1}
$$

The activity vector, its magnitude and the mass-normalized diagnostics are

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad
q_1(t)=\frac{|M_1(t)|}{M_0(t)},\qquad q_2(t)=\frac{|M_2(t)|}{M_0(t)}.
\tag{2}
$$

Complex division puts the initial incumbent at $P(0)=1$. Primary angle is $\arg P$, withheld when $R<0.1$; sensitivities use 0.05 and 0.2 without changing the vector. Uniform activity cancels from the first/second nonzero angular moment numerators. It changes $M_0$ and therefore mass-normalized $q_1, q_2$, which are retained alongside the initial-moment normalization. A small first moment can arise from cancellation between strong opposed sites. Neither a nonzero vector nor its angle uniquely identifies a protrusion, tracked site, velocity or direction of migration.

The new common active-loss perturbation uses the frozen basal half-contrast core, background and cosine taper. The removed primary profile is $\varepsilon$ times taper times basal excess above background, with primary$\varepsilon=0.005$. In these three native chemical models the same material is transferred pointwise to the explicit inactive pool. All native fields evolve under the same actual cue as their sham. The perturbation direction stays the original basal direction even at the challenged t10 snapshot; there is no renewed inactive-withdrawal/return step.

$$
D(t)=A_{\mathrm{sham}}(t)-A_{\mathrm{probe}}(t),\qquad
F=1-\frac{D(\tau)}{D(0)},\qquad K=\frac{F}{\tau}.
\tag{3}
$$

All amounts in this assay integrate the fixed original basal core. It does not follow a translating/broadening site. Fixed clocks are WP/Holmes 0.1778279410038923 and Goryachev 0.005623413251903491 native units, identical to the matched linear assay and fixed across sizes/placements within each model. They are not a calibration across models. The main map displays $F_{10}(C, d)/F_0(C)$ with each exact same-size basal control. Raw$F, K, D(0), D(\tau)$, errors, whole-active and outside-core deficits, actual sham/probe changes and rates remain linked.

All 702 ring states pass the numerical, amplitude and prefix-saturation gates; 7020 native trials were independently audited. Positive paired-deficit damping is conditional recovery along this disturbance, not actual net maintenance or generic stability. A black dot flags the eight t10 states lacking a resolved distinct old-site association; those remain regional responses. The independent linear/ring 13-pair audit concerns independently prepared states with different native meshes, so close agreement is not exact mirrored-state identity or agreement within numerical$E$.

The complete original matrix has 867 planned cells:652 complete histories, 198 geometry exclusions and 17 unqualified Holmes C80 preparations. No excluded state receives a calculated vector. Phase-local recovery time is shifted by 1000 before joining; the common withdrawal state is retained once. The four readouts are onset 1000, 4000, 9000, 14000, equivalent to 0, 3000, 8000, 13000 after withdrawal, all exactly present in native samples. Minima and maximum unmasked angles are sampled summaries, not interpolated event times. Later source sampling differs (5 or 50 native units); no greater temporal precision is inferred.

The lean NPZ stores concatenated native-time traces and offsets by complete condition: $M_0$, real/imaginary$M_1, M_2$, $P_x, P_y, R, q_1, q_2$. The 867-row matrix preserves unavailable statuses; 2608 readout rows retain the primary and sensitivity angle gates. Native sources are referenced by immutable path/hash. All extraction is read-only; no solver or deterministic persistence fit is involved.

