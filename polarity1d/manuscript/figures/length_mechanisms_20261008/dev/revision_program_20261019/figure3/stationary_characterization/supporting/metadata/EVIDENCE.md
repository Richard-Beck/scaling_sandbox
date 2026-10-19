# Cue-on response and cue-free retention in the same 80 histories

Generated: 2026-10-19T16:28:14.131151+00:00. This is a complementary analysis of the existing 80 standardized-angle histories. No new integrations, sizes, cue protocols or moving cues are introduced. The exact existing initial overlaps, physical circumference histories, continuous target phases, native archive provenance and independent field-tail flags remain attached to every row.

The primary sheet directly pairs withdrawal and retained angles at each tested duration. The paired magnitude/amount sheet shows why activity-vector magnitude and total active amount must remain distinct. The supplied tables contain sampled extrema, fixed-region depletion/rebound and threshold brackets so these quantities need not be inferred from a single endpoint.

## Observed examples

- Wave-pinning, target 90°, duration 10: withdrawal 8.516° → retained 0.071°; R 0.978 → 1.000.
- Goryachev, target 90°, duration 300: withdrawal 33.479° → retained -0.028°; R 0.721 → 1.000.
- Goryachev, target 45°, duration 300: withdrawal 29.592° → retained 44.100°; R 0.924 → 1.000.
- Spring, target 90°, duration 1: withdrawal 37.706° → retained 29.745°; R 1.076 → 1.000.

These examples show that cue-on angular deflection is not a retained-angle surrogate. Post-withdrawal evolution may return toward the original direction or continue toward the input normal. A first moment alone does not identify site lineage: the existing native final-profile checks show all 60 non-antipodal endpoints are one translated full-strength basal region, while Holmes antipodal cancellation requires its separate local-field evidence.

## Measured quantities

The raw physical first moment M₁ has real and imaginary components along fixed ring reference axes. Normalizing by the complex initial M₁ both rescales the vector and rotates its reference axes once by the initial M₁ phase, placing P(0) at 1 + 0i. P therefore measures orientation relative to the initial resultant, with no subsequent rotating reference frame. Its magnitude depends on both total activity and spatial alignment. M₀ is physical total active amount; the native API's M₁/M₀ differs from P = M₁/M₁(0). We export physical M₀ and M₁ directly from every native active field, plus both normalizations. Spring uses instantaneous C(t), and its prescribed material cue offset remains fixed while its physical offset changes as fraction × C(t).

For every history the compact summary exports withdrawal/final angles and R, physical M₀ ratios at withdrawal/final and sampled min/max, sampled minimum R, cue-on and release fixed-original-region minima, final regional value, net release change and rebound from the observed minimum. “Minimum” here always means minimum over saved samples, not a proven continuous-time extremum. Raw angles are preserved separately; plotted angles are masked below R = 0.1.

The fixed-original-region observable is the physical active amount in the frozen basal half-contrast material interval, minus its transported initial background amount, divided by initial excess. Pure affine dilution leaves it unchanged. It can change through profile translation, redistribution or local chemistry, so depletion is not equivalent to destruction or death of an independently tracked incumbent. Overlapping inputs cannot be relabelled as a separate challenger. Rebound is quantified as final excess minus its sampled minimum and as final minus withdrawal; threshold recrossing is explicit.

The crossing table has 1440 records: 80 conditions × 3 levels (0.9/0.5/0.1 of initial regional excess) × 3 phases (cue-on/release/full) × 2 events (first saved entry or entry retained through all later saved phase frames). Times use cue-onset origin. Each first/retained entry has a saved-sample bracket or explicit finite-horizon censoring, and the table reports whether the region rises above that threshold after its first sampled entry. Saved sampling can miss shorter excursions; these brackets do not establish a minimum effective cue duration.

The independent existing whole-grid trace replay agrees to 6.66e-16. Actual native phase hashes and exact primary restart/continuity are checked again. Prior full-field tail flags remain 75 passes / 1 failure / 4 unavailable; activity-angle behaviour is not used to overwrite them. Native physical input full widths stay 3 for chemical models and 0.2 for Spring. Shared centre angles do not match physical or angular footprint, dose or native clocks across models.

## Outputs

- `stationary_endpoint_pairs.csv`: compact paired observables across all 80 conditions.
- `stationary_whole_grid_summary.csv`: all paired and transient observables, exact native provenance, anchor/overlap/geometry and field-tail flags.
- `stationary_whole_grid_traces.npz`: keys `c{condition_id}_{observable}`, including direct physical M0, M1_real, M1_imag, P_real, P_imag, R, q1, q2, C, fixed-region amount/excess and masked/raw angles.
- `stationary_fixed_region_crossings.csv`: all censoring-aware saved-sample brackets and recrossing flags.

The selected local-field/vector sketch produced by the other figure worker is a separate view of these same histories. Existing main/prototype art remains preserved; this checkpoint does not choose the final Figure 3 composition.


$$
M_n(t)=C(t)\int_0^1u(q,t)e^{2\pi i n q}\,\mathrm{d}q,
\quad P(t)=\frac{M_1(t)}{M_1(0)},
\quad R(t)=|P(t)|,
\quad q_1(t)=\frac{|M_1(t)|}{M_0(t)}.
\tag{1}
$$

$$
E_{\Omega_0}^{\mathrm{tr}}(t)=
\frac{C(t)\int_{\Omega_0}u(q,t)\,\mathrm{d}q-C_{\mathrm b}u_{\mathrm{bg},0}|\Omega_0|}
     {C_{\mathrm b}\int_{\Omega_0}[u(q,0)-u_{\mathrm{bg},0}]\,\mathrm{d}q}.
\tag{2}
$$
