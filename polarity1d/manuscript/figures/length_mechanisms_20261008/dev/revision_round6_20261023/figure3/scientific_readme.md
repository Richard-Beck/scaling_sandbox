# Unexpected

The requested map tests support the native linear response away from the edge: the maximum error is 1.69% of the condition's map amplitude. Symmetric centre/antipode pulses give essentially zero shared translation, while off-centre pulses shift the complete domain. This does not repair Goryachev's absolute drift disagreement.

The covariance contribution maps localize 99.59–100.00% of predicted drift within two input correlation lengths of the two edges at the endpoint sizes. Small signed plateau/basal contributions and all negative bins are retained.

The expanded drift-regime alpha bridge contains 53 qualified matched-input pairs, including 15 Wave-pinning pairs at f=0.003, 0.005 and 0.01. Wave-pinning predicts α=2.48–2.88, above 2, and can now be compared with the nontrivial observed size advantage. The narrow predicted range primarily tests the excess α above 2 rather than the within-model trend. Goryachev still has no qualified pair. Every excluded row remains with its explicit reason.

The matched-size comparison is now explicit. The operational separated-cue plateau is the median of the furthest three finite observed two-region half-turns, with every selected gap and range exported. It is not a fitted asymptote.

# Scientific results

Terminal largest/smallest edge-gain folds are Wave-pinning 0.701, Holmes 0.914, Spring 0.989, Goryachev 0.088. Wave-pinning's f=0.001 size-normalized drift ratio remains 0.489645 predicted and 0.473243 measured. The full spatial prediction is approximately the squared terminal edge-gain fold; it is not mathematically identical when the susceptibility shape and covariance change with size.

The linear range is bracketed by the last tested amplitude within 10% of nominal gain and the first stopping amplitude. All 20 native-size sweeps stop on gain departure, before the continuous-cue peak and without losing the unique half-contrast domain:

- Wave-pinning: first >10% departures 0.02048 at C=20 and 0.01024 at C=60, size fold 0.5; last tested linear amplitudes 0.01024 and 0.00512
- Holmes: first >10% departures 0.0512 at C=20 and 0.0512 at C=60, size fold 1; last tested linear amplitudes 0.0256 and 0.0256
- Spring: first >10% departures 0.1024 at C=1 and 0.1024 at C=3, size fold 1; last tested linear amplitudes 0.0512 and 0.0512
- Goryachev: first >10% departures 0.32768 at C=20 and 1.31072 at C=60, size fold 4; last tested linear amplitudes 0.16384 and 0.65536

Thresholds are geometric-grid brackets, not precisely located bifurcations. Fixed-dose durations test 12.5, 25, 50, 100, 200 chemical units and Spring/10 at the unchanged nominal dose.

# Methods and provenance

Starting reference: `manuscript_outline/outline_v11.pdf` and `.json`. A–B source coordinate files are byte-identical to round 5 and all previous editions remain intact. C–H use deterministic assays and existing observed drift/heading tables. No stochastic realization, kinetic parameter fit or retuning was introduced.

The 28 conditions include all 20 original f=0.001 sizes and 8 added B-matched sizes. Original states, adjoints and covariance kernels are reused from the frozen round 5 campaign. Added states periodically resample the nearest original basal onto unchanged native cell density and settle at the exact mean input for 16000 chemical or 1600 Spring units. Eligibility uses the unchanged unique-region, contrast≥0.05, native RHS≤10⁻⁶ and tail-change≤10⁻⁴ gates. All 28 pass. BDF uses rtol 2×10⁻¹⁰, atol 2×10⁻¹², exact pulse withdrawal restart, and the unchanged native sparse Jacobian. Every integration ran under SLURM 21334289 on one CPU, small QoS, 4 GB, 30-minute ceiling; actual workers took under 36 seconds.

For the full native state y, translation mode q, and normalized adjoint p:

$$
Jq=0,\qquad J^{\mathsf T}p=0,\qquad p^{\mathsf T}q=1,\qquad g_j=(p^{\mathsf T}B)_j. \tag{1}
$$

The plotted standard-pulse map is the convolution of raw native g with the exact periodic Gaussian pulse shape h, per peak-amplitude×duration dose:

$$
G(x_0)=\sum_j g_j h(x_j-x_0),\qquad h(d)=\exp[-d^2/(2\ell^2)]. \tag{2}
$$

Raw cell-input g appears as thin lines in the supplement; it is not substituted for G on the main map. Points measure the terminal common shift of the two half-contrast edges divided by pulse dose, subtracting matched unpulsed controls. Every test retains left/right edges, activity centroid, width and region count through 20 nominal pulse durations. Spring uses material fractions multiplied by its mean circumference, excluding uniform extension.

Main C/E/G/H use raw native g without any empirical amplitude normalization. The historical round 5 prediction used the mean of the two independently measured edge/native gain factors s, close to one; that variant is exported separately and retained in the historical supplement. Observed stochastic drift is never used. The raw-native contribution partition uses

$$
c_j=g_j\sum_k g_k\int_0^\infty\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v,\qquad D_s=\sum_j c_j. \tag{3}
$$

The input covariance is the exact native lognormal image of the discrete spatial Gaussian filter with variance-normalized OU-knot interpolation, zero spatial mode retained. Round5's deterministic covariance quadrature is reused. Edge zones use two input correlation lengths; where zones overlap the nearest edge assigns the cell. Plateau and basal remainders exclude edge zones, so the four signed fractions sum to one. No negative contribution is clipped. Every map bin and contribution is in `supporting/data/susceptibility_and_contributions.csv`.

At the right edge, amplitudes double from nominal/2 until gain departs by >10%, topology changes, or the native continuous-cue peak is reached. The first stopping point is retained. Duration tests hold nominal aτ constant and retain the same final 20τ nominal observation horizon. No response threshold or amplitude was tuned to improve agreement.

$$
\alpha_{\mathrm{pred}}=2-\frac{\mathrm{d}\ln D_s}{\mathrm{d}\ln \bar C}. \tag{4}
$$

The local derivative uses the fixed adjacent-size neighborhoods of Figure 4E (two members at grid endpoints, three at interior sizes) and mean-state physical circumference. Each prediction is evaluated at the observation's own exact input using the completed deterministic native mean-state, adjoint and covariance calculations at that input. The main map remains at f=0.001. No assumption of input-independent susceptibility or extrapolation from that map is made. The drift subset requires f≤0.01 and multiple-domain fraction≤0.05 at every size in the fixed neighbourhood; linearity error must be≤0.1 and all original stationarity gates must pass. The observed exponent retains its existing coverage and shared-lag eligibility rules. Ineligible or unavailable estimates remain rows with their reason.

Length denotes WGD-associated enlargement. Observed migration summaries, conditional OU persistence estimates and simulated 1D polarity-heading readouts remain distinct. Domain drift is not identified with total heading decay or empirical migration persistence.

# Outputs and rebuild

Main: `figure3_round6.pdf`, `.png`, `.svg`. Complete supplement appends the unchanged round 5 reviewed supplement after 11 new pages. Compact CSVs include every plotted source point, raw/smoothed maps, signed contribution bins, pulse tests/time courses, amplitude brackets, matched-size comparison and alpha bridge. `supporting/metadata/provenance.json` records hashes and exact sources.

Run `./rebuild.sh` to verify/collect/render from completed workers. Run `./rebuild.sh --native` only to prepare/submit any missing deterministic workers and then render; it requires SLURM and does not repeat completed conditions. The output tree contains only small derivative/native scalar archives; no large full-field campaign was created.
