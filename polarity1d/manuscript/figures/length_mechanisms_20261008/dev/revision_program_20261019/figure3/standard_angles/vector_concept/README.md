# Spatial activity and the physical complex vector

This separate reviewable sketch uses the same four suggested cases from the completed 80-condition stationary-cue checkpoint. It adds no native integration and leaves all previous art, content and outline files intact.

- [Complex-vector concept sheet](stationary_vector_concept.png) · [PDF](supporting/vectors/stationary_vector_concept.pdf) · [SVG](supporting/vectors/stationary_vector_concept.svg)
- [Actual signed spatial changes](stationary_vector_native_changes.png)
- [Selected measurements](supporting/data/selected_case_measurements.csv) · [Native profile snapshots](supporting/data/selected_native_profiles.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Source and validation record](supporting/metadata/validation.json)

For outward normal angle θⱼ and native active concentration uⱼ, the physical ring quadrature is

$$
M_1(t)=\sum_j u_j(t)\,\Delta s_j(t)e^{i\theta_j},\qquad
\Delta s_j(t)=\frac{C(t)}{N},\qquad
P(t)=\frac{M_1(t)}{M_1(0)}.
\tag{1}
$$

$$
M_0(t)=\sum_j u_j(t)\,\Delta s_j(t),\qquad R(t)=|P(t)|.
\tag{2}
$$

Each complex contribution combines activity, physical segment length and outward-normal direction. Opposed contributions can cancel despite substantial total activity. P is normalized by the initial complex physical moment: division rotates the initial direction to the positive real axis and scales by the initial magnitude, giving P(0)=1. It does not normalize by the total active amount at each frame, so R may exceed 1. The physical M₁/M₀ shape-normalized helper is not substituted. Spring uses its actual circumference at every sample. No profile is divided by an independent history maximum.

## What the actual four cases show

- Wave-pinning, 180°, T=300: the measured real-axis path changes sign near the origin and retains a reversed final full-strength vector. The minimum sampled magnitude is about .00262 at t=205; this sampling does not assert an exactly sampled zero.
- Goryachev, 90°, T=300: the cue produces a substantial excursion to 33.48° with R=.721 at withdrawal, followed by return to −.0278° and R≈1. Its raw active peak changes from about 64.95 to 37.31 at withdrawal. The near-original-axis endpoint does not imply an absent response.
- Holmes 3, 180°, T=1000: the final two-peak profile has R=.03993 while physical total active amount is .99746 of initial. Direction is masked, and loss of the resultant does not mean loss of nearly all activity.
- Spring, 90°, T=10: the cue produces 76.09° with R=1.281 at withdrawal and retains 72.44° with R≈1 after release. The final native profile is one shifted full-strength activity region. Circumference is included in both physical moments. Its sampled minimum R differs from 1 only at numerical roundoff, rather than marking a weak-resultant episode. About 35.6% of the initial cue support overlaps the original half-contrast core, so the case does not start with two separate activity sites.

The antipodal chemical paths have imaginary components at numerical roundoff; equal complex-plane scales prevent numerical asymmetry from being magnified. The R<.1 disk concerns heading interpretation, not validity of the measured weak complex vector. The displayed profile time at minimum R is the exact stored native frame, with no nearest-profile substitution. The signed Δu history companion exposes spatial changes in native concentration units; its colour limits are explicitly numerical rather than a peak normalization.

The physical total active amount M₀ is distinct from total conserved active+inactive pool. Near-constant M₀ with weak M₁ supports cancellation or redistribution of activity rather than global activity disappearance; it does not assign a particular site lineage or functional competition. Shared cue-centre angles do not match native cue dose, footprint or timescale. A finite 13,000-unit withdrawal does not prove asymptotic stability. Activity orientation is not a measured movement vector.

## Source and validation

The selected native challenge/recovery archives are identified and hashed in the measurement table and source JSON. Physical P is recomputed directly from every selected native profile and matches the retained standard-angle trace. All four suggested cases were retained after inspecting their actual response. The complete 80-condition experiment remains in the parent standard_angles directory.
