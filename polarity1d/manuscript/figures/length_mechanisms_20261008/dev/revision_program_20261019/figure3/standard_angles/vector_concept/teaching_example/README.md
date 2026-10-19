# A two-state worked example of the activity vector

[PNG](holmes_two_state_vector.png) · [PDF](supporting/vectors/holmes_two_state_vector.pdf) · [SVG](supporting/vectors/holmes_two_state_vector.svg)

This small teaching sketch uses one already selected native Holmes history. It adds no integration and leaves the preceding concept sheets and Figure 3 artwork unchanged. The measured profile is shown as coloured activity on a schematic circle. Activity itself is real; complex numbers package its two Cartesian directional sums:

$$
M_1=\sum_j u_j\Delta s_j\cos\theta_j+i\sum_j u_j\Delta s_j\sin\theta_j,
\qquad M_0=\sum_j u_j\Delta s_j.
\tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad\Delta s_j=\frac{40}{800}=0.05.
\tag{2}
$$

Initial complex normalization rotates and scales the vector so P(0)=1. Opposing activity can cancel in M₁ while adding positively in M₀. The final R≈0.040 and M₀/M₀(0)≈0.9975 demonstrate that distinction in the actual saved profile. Raw u is neither complex nor independently normalized between states. This is an activity descriptor, without a movement law.

[Caption and limits](supporting/metadata/CAPTION.md) · [Two-state values](supporting/data/two_state_measurements.csv) · [Every physical cell contribution](supporting/data/native_cell_contributions.csv) · [Exact sector sums](supporting/data/exact_bin_sums.csv) · [Sources and validation](supporting/metadata/validation.json)
