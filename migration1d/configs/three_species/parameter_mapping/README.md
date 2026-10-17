# Rac–Rho–Cdc42 polarity parameterization for a coupled membrane–cytoplasm strip

## Scope

A “three-species” Rho-family polarity model is most naturally interpreted as a model containing the three small GTPases

\[
\mathrm{Cdc42},\qquad \mathrm{Rac},\qquad \mathrm{Rho},
\]

but not literally only three dynamical fields. Each GTPase cycles between a slowly diffusing, membrane-associated active form and a more rapidly diffusing inactive pool. The minimal mass-conserving representation therefore has **six fields**:

\[
(C,C_i),\qquad (R,R_i),\qquad (\rho,\rho_i),
\]

where \(C,R,\rho\) denote active Cdc42, Rac and Rho, and \(C_i,R_i,\rho_i\) denote their corresponding inactive pools.

A more mechanistically explicit representation separates inactive membrane-bound and inactive cytosolic GTPase. It therefore has **nine GTPase fields**:

\[
(C_a,C_{mi},C_{ci}),\qquad
(R_a,R_{mi},R_{ci}),\qquad
(\rho_a,\rho_{mi},\rho_{ci}).
\]

This distinction is important for cell-size scaling. The six-field model is the natural analogue of the two-field wave-pinning model and is the most convenient starting point for analysis. The nine-field bulk–surface model is preferable when membrane/cytosol exchange and the changing surface-to-volume ratio are themselves part of the question.

The most useful canonical lineage is:

1. **Marée et al. 2006**: an early multiscale Rac–Rho–Cdc42 / phosphoinositide / actin model for keratocyte polarization and motility.
2. **Jilkine, Marée & Edelstein-Keshet 2007**: systematic analysis of Rac–Rho–Cdc42 spatial segregation and inhibitory crosstalk.
3. **Holmes et al. 2012**: a reduced, nearly one-dimensional mammalian-cell model developed against direct Rac-activation experiments in HeLa cells, with explicit analysis of geometry and parameter regimes.
4. **Marée, Grieneisen & Edelstein-Keshet 2012**: a two-dimensional motile-cell implementation coupling the GTPases to phosphoinositides and the cytoskeleton.
5. **Spill et al. 2016**: an explicit three-dimensional membrane–cytosol implementation of the Rac–Rho–Cdc42 signaling network, including inactive membrane/cytosolic exchange.
6. **Comer et al. 2021**: a useful modern reimplementation of the Marée model with a clear parameter table and cell-shape sensitivity analysis.

There is **no single universal Rac–Rho–Cdc42 wiring diagram**. The crosstalk depends on cell type and stimulus. For the present strip model, the Holmes 2012 network is a particularly useful default because it was designed for a nearly 1D mammalian-cell experiment, keeps the GTPase layer relatively compact, and was explicitly constructed to reproduce polarization under direct Rac stimulation.

---

## Geometry

We consider a one-dimensional representation of a cell of length \(L\), consisting of two adjacent compartments:

* a membrane strip running from \(x=0\) to \(x=L\);
* an adjacent cytoplasmic strip over the same axial coordinate.

For a physically explicit reduction of a bulk–surface cell geometry, define

* \(s_m\): membrane area represented per unit axial length \([\mu\mathrm{m}]\);
* \(A_c\): cytoplasmic volume represented per unit axial length \([\mu\mathrm{m}^2]\).

The effective cytoplasmic depth associated with the modeled membrane is

\[
h=\frac{A_c}{s_m}.
\]

For geometrically similar cells,

\[
h\sim \frac{V}{A_{\rm mem}}
\]

scales with the characteristic cell size.

This creates the same geometric issue as in the one-GTPase wave-pinning model, but now separately for Cdc42, Rac and Rho.

---

# 1. The canonical six-field formulation

## Generic equations

Let

\[
G\in\{C,R,\rho\}
\]

denote active Cdc42, Rac or Rho, and let \(G_i\) denote the corresponding inactive pool.

The reduced Holmes/Marée-type equations have the form

\[
\frac{\partial G}{\partial t}
=
D_m\frac{\partial^2G}{\partial x^2}
+
Q_G(C,R,\rho,\ldots)\frac{G_i}{G_{\rm tot}}
-d_GG,
\]

\[
\frac{\partial G_i}{\partial t}
=
D_i\frac{\partial^2G_i}{\partial x^2}
-
Q_G(C,R,\rho,\ldots)\frac{G_i}{G_{\rm tot}}
+d_GG.
\]

Here

* \(D_m\) is the diffusion coefficient of active, membrane-associated GTPase;
* \(D_i\) is the effective diffusion coefficient of the inactive pool;
* \(Q_G\) is a GEF-mediated activation rate;
* \(d_G\) is the GAP-mediated inactivation rate;
* \(G_{\rm tot}\) is a reference total concentration for that GTPase.

Each GTPase has its own conserved amount:

\[
M_C=\int_0^L(C+C_i)\,dx,
\]

\[
M_R=\int_0^L(R+R_i)\,dx,
\]

\[
M_\rho=\int_0^L(\rho+\rho_i)\,dx.
\]

Thus the important difference from a single-species wave-pinning model is that there are **three separately conserved pools**, coupled through the activation functions.

The active forms diffuse slowly, while the inactive pools redistribute rapidly. This fast inactive redistribution is not merely a numerical convenience: Jilkine et al. showed that it is central to stabilizing stationary spatial segregation rather than allowing activation fronts simply to propagate through the domain.

---

# 2. Canonical GTPase crosstalk

## Jilkine–Marée lineage

Jilkine, Marée & Edelstein-Keshet systematically examined several proposed Rac–Rho–Cdc42 crosstalk schemes. Their key conclusion was that a cooperative, nonlinear inhibitory network can generate the experimentally observed spatial segregation of front-associated Cdc42/Rac from rear-associated Rho.

They considered linear, Michaelis–Menten and Hill-type interaction functions. The Hill formulation is the most useful precedent for a compact polarity model because it permits a stable resting state together with a threshold for rapid polarization.

A representative form of the older Marée/Jilkine network is

\[
Q_C^0
=
\frac{I_C}
{1+\left(\rho/a_1\right)^n},
\]

\[
Q_R^0
=
I_R+\alpha C,
\]

\[
Q_\rho^0
=
\frac{I_\rho+\beta R}
{1+\left(C/a_2\right)^n}.
\]

The biological interpretation is

\[
\rho \dashv C,
\qquad
C\rightarrow R,
\qquad
R\rightarrow\rho,
\qquad
C\dashv\rho.
\]

Thus Cdc42 and Rho form an antagonistic front–rear pair, while Rac is recruited downstream of Cdc42 and feeds into Rho activation.

This topology was used in the Marée motility-model lineage and remains a useful canonical alternative. It should not, however, be treated as a universal mammalian signaling diagram.

---

## Holmes 2012 network

Holmes et al. revisited the network using HeLa cells in narrow microfluidic channels in which a Rac GEF was recruited directly to the membrane by a graded synthetic stimulus.

A key result was that the pre-existing Marée-type GTPase/PI network could polarize in the **wrong direction** under direct Rac stimulation. Holmes et al. therefore rebuilt the signaling network stepwise and adopted a minimal three-GTPase circuit compatible with those experiments.

Ignoring phosphoinositides, the core activation functions can be written as

\[
Q_C
=
\frac{I_C}
{1+\left(\rho/a_1\right)^n},
\]

\[
Q_R
=
I_{R1}+\alpha C+S_R(x,t),
\]

\[
Q_\rho
=
\frac{I_\rho}
{1+\left(R/a_2\right)^n}.
\]

The core crosstalk is therefore

\[
\rho\dashv C,
\qquad
C\rightarrow R,
\qquad
R\dashv\rho.
\]

This circuit creates an effective front–rear antagonism:

* high Cdc42 promotes Rac;
* high Rac suppresses Rho;
* high Rho suppresses Cdc42.

The direct experimental stimulus enters naturally as an additive Rac-activation term

\[
S_R(x,t).
\]

For example, a graded Rac stimulus can be represented as

\[
S_R(x,t)
=
A(t)
\left[
1+g\left(\frac{x}{L}-\frac12\right)
\right],
\]

or a localized pulse as

\[
S_R(x,t)
=
A
\exp\left[
-\frac{(x-x_0)^2}{2\sigma_x^2}
\right]
\mathbf{1}_{t_0<t<t_1}.
\]

This makes the Holmes network especially attractive for controlled polarity-reversal and rapidly changing-cue experiments.

---

# 3. Phosphoinositide feedback

The full Holmes/Marée models also include PIP, PIP\(_2\) and PIP\(_3\). PIP\(_3\) feeds back onto Cdc42 and Rac activation.

A Holmes/Spill-type Rac activation function is

\[
Q_R
=
I_{R1}
+
I_{R2}f_{\rm PI}\frac{P_3}{P_{3b}}
+
\alpha C
+
S_R(x,t).
\]

The Cdc42 activation rate can likewise be multiplied or augmented by a PIP\(_3\)-dependent term.

The PI layer is biologically meaningful and strongly affects competition between multiple fronts. Marée et al. found that PI feedback can accelerate communication between competing Rac/Cdc42-rich regions by driving stronger depletion gradients in the fast inactive GTPase pools.

However, the PI layer adds three more PDEs and several additional kinetic parameters.

For a first strip model aimed at understanding **size, polarity stability and reversal**, the most defensible hierarchy is therefore:

1. start with the six GTPase fields alone;
2. establish the polarized, non-polarized and multipolar regimes;
3. add PIs only if the GTPase-only model cannot reproduce the desired competition/reorientation dynamics.

---

# 4. Published parameterizations

## Holmes et al. 2012

Holmes et al. provide the cleanest baseline for a 1D Rac–Rho–Cdc42 strip model.

Their baseline cell length is

\[
L_0=20~\mu\mathrm{m}.
\]

The total GTPase concentrations are

\[
C_{\rm tot}=2.4~\mu\mathrm{M},
\]

\[
R_{\rm tot}=7.5~\mu\mathrm{M},
\]

\[
\rho_{\rm tot}=3.1~\mu\mathrm{M}.
\]

Their reported activation parameters are approximately

| Parameter | Value | Meaning |
|---|---:|---|
| \(I_C\) | \(2.95~\mu\mathrm{M\,s^{-1}}\) | basal Cdc42 activation |
| \(I_{R1}\) | \(0.2~\mu\mathrm{M\,s^{-1}}\) | basal Rac activation |
| \(I_{R2}\) | \(0.2~\mu\mathrm{M\,s^{-1}}\) | PI-dependent Rac activation coefficient |
| \(I_\rho\) | \(6.6~\mu\mathrm{M\,s^{-1}}\) | basal Rho activation |
| \(a_1\) | \(1.25~\mu\mathrm{M}\) | Rho level for half-maximal Cdc42 inhibition |
| \(a_2\) | \(1~\mu\mathrm{M}\) | Rac level for half-maximal Rho inhibition |
| \(n\) | 3 | inhibitory Hill coefficient |
| \(\alpha\) | \(0.65~\mathrm{s^{-1}}\) | Cdc42-dependent Rac activation |
| \(d_C,d_R,d_\rho\) | \(1~\mathrm{s^{-1}}\) | active-GTPase inactivation |
| \(D_m\) | \(0.1~\mu\mathrm{m}^2\,\mathrm{s^{-1}}\) | active/membrane diffusion |
| \(D_i\) | \(50~\mu\mathrm{m}^2\,\mathrm{s^{-1}}\) | effective inactive diffusion |

The resulting transport ratio is

\[
\boxed{
\frac{D_i}{D_m}=500.
}
\]

The inactivation timescale is

\[
\boxed{
\tau_G=d_G^{-1}\approx1~\mathrm{s}.
}
\]

These parameters were assembled from previous literature and adjusted to the HeLa polarization system rather than fitted as a uniquely identifiable biochemical parameter set. They should therefore be treated as a well-established **model regime**, not as precise molecular constants.

---

## Marée et al. 2012 / Comer et al. 2021 parameterization

The related Marée model uses the older crosstalk topology

\[
\rho\dashv C,
\qquad
C\rightarrow R,
\qquad
R\rightarrow\rho,
\qquad
C\dashv\rho.
\]

A later VCell reimplementation reports the Marée parameter set clearly:

| Parameter | Value |
|---|---:|
| \(C_{\rm tot},R_{\rm tot},\rho_{\rm tot}\) | \(2.4,\,7.5,\,3.1~\mu\mathrm{M}\) |
| \(I_C,I_R,I_\rho\) | \(2.95,\,0.5,\,3.3~\mu\mathrm{M\,s^{-1}}\) |
| \(a_1\) | \(1.25~\mu\mathrm{M}\) |
| \(a_2\) | \(1~\mu\mathrm{M}\) |
| \(n\) | 3 |
| \(\alpha\) | \(4.5~\mathrm{s^{-1}}\) |
| \(\beta\) | \(0.3~\mathrm{s^{-1}}\) |
| \(d_C,d_R,d_\rho\) | \(1~\mathrm{s^{-1}}\) |
| \(D_m\) | \(0.1~\mu\mathrm{m}^2\,\mathrm{s^{-1}}\) |
| \(D_i\) | \(50~\mu\mathrm{m}^2\,\mathrm{s^{-1}}\) |

Again,

\[
D_i/D_m=500.
\]

The agreement in the totals, diffusion coefficients, inhibitory thresholds and inactivation rate across this model lineage is more informative than the differences in individual activation coefficients, which depend strongly on the assumed network topology.

---

## Jilkine et al. 2007

Jilkine et al. are especially important for interpreting the model rather than for supplying a single canonical parameter vector.

They compared several functional forms for GTPase crosstalk:

* linear;
* Michaelis–Menten-like;
* cooperative Hill functions.

Their principal dynamical conclusions were:

* cooperative nonlinearity is required for the desired stable-rest / rapid-polarization behavior;
* bistability in the local GTPase kinetics provides the threshold for segregation;
* rapid redistribution of inactive GTPase stabilizes the spatial fronts;
* crosstalk topology matters, and different plausible biological networks are not dynamically equivalent.

Thus the most portable result from Jilkine et al. is not a single number but the expectation

\[
\boxed{n>1}
\]

with

\[
\boxed{n\approx3}
\]

a natural default in the later Holmes/Marée models.

---

# 5. The nine-field membrane–cytosol formulation

## Why the six-field model is a reduction

For each GTPase, a more explicit biochemical picture contains three states:

1. active membrane-bound GTPase \(G_a\);
2. inactive membrane-bound GTPase \(G_{mi}\);
3. inactive cytosolic GTPase \(G_{ci}\).

GEFs activate the membrane-bound inactive form,

\[
G_{mi}\rightarrow G_a,
\]

GAPs inactivate the active form,

\[
G_a\rightarrow G_{mi},
\]

and inactive GTPase exchanges with the cytosol,

\[
G_{mi}
\rightleftharpoons
G_{ci}.
\]

The six-field models assume that exchange between \(G_{mi}\) and \(G_{ci}\) is sufficiently rapid that the two can be replaced by one composite inactive pool with an effective diffusion coefficient.

Holmes et al. explicitly derive this quasi-steady-state reduction.

---

## Spill et al. 2016

Spill et al. implement the three-state picture in a full three-dimensional cell.

For each GTPase,

\[
\frac{\partial G_a}{\partial t}
=
D_M\nabla_S^2G_a
-
\delta_GG_a
+
I_G\frac{G_{mi}}{G_{\rm tot}},
\]

\[
\frac{\partial G_{mi}}{\partial t}
=
D_M\nabla_S^2G_{mi}
+
\delta_GG_a
-
I_G\frac{G_{mi}}{G_{\rm tot}}
+
k_{\rm on}L_I G_{ci}
-
k_{\rm off}G_{mi},
\]

while the cytosolic form satisfies bulk diffusion,

\[
\frac{\partial G_{ci}}{\partial t}
=
D_C\nabla_V^2G_{ci},
\]

with the membrane boundary condition

\[
-D_C\mathbf n\cdot\nabla G_{ci}
=
k_{\rm on}L_I G_{ci}
-
k_{\rm off}G_{mi}.
\]

This model is conceptually important for size scaling because it makes the membrane and cytoplasmic compartments physically distinct.

In particular, changing cell size changes the relationship between

\[
\text{membrane binding capacity}
\]

and

\[
\text{cytoplasmic reservoir volume}
\]

even if all microscopic rate constants and concentrations are unchanged.

---

# 6. Mapping onto the strip geometry

There are three useful levels of abstraction.

## Option 1: six-field line-density model

This is the closest direct analogue of the standard 1D Holmes equations.

For each GTPase define

\[
U_G(x,t)
=
\text{active amount per unit axial length},
\]

\[
V_G(x,t)
=
\text{inactive amount per unit axial length}.
\]

Then

\[
\frac{\partial U_G}{\partial t}
=
D_mU_{G,xx}
+
F_G,
\]

\[
\frac{\partial V_G}{\partial t}
=
D_iV_{G,xx}
-
F_G.
\]

The conserved amount is

\[
M_G
=
\int_0^L
(U_G+V_G)\,dx.
\]

This is the safest representation for reproducing the native Holmes/Jilkine dynamics.

Its weakness is that membrane area and cytoplasmic volume are hidden in the one-dimensional projection.

---

## Option 2: six-field membrane–cytoplasm strip

For cell-size scaling, define

\[
G(x,t)
=
\text{active membrane-associated concentration},
\]

\[
G_i(x,t)
=
\text{inactive cytoplasmic concentration}.
\]

If the active membrane variable is represented as an equivalent concentration in a membrane reaction layer of thickness \(\ell_m\), then the amount represented by a slice is

\[
dM_G
=
s_m
\left(
\ell_mG+hG_i
\right)dx.
\]

The total amount is therefore

\[
M_G
=
s_m\int_0^L
\left(
\ell_mG+hG_i
\right)dx.
\]

If

\[
F_G
=
Q_G\frac{G_i}{G_{\rm tot}}
-
d_GG
\]

is written in concentration-per-time units for the membrane reaction layer, then

\[
\frac{\partial G}{\partial t}
=
D_mG_{xx}
+
F_G,
\]

\[
\frac{\partial G_i}{\partial t}
=
D_iG_{i,xx}
-
\frac{\ell_m}{h}F_G.
\]

The geometric exchange factor is therefore

\[
\boxed{
\eta_h=\frac{\ell_m}{h}.
}
\]

This is the direct six-field analogue of the membrane–volume correction in the one-GTPase strip model.

The key point is that the published Holmes concentrations cannot be transferred literally into a true surface-density / volume-density model without specifying a conversion length such as \(\ell_m\). The original 1D equations have already absorbed this geometric conversion into their projected variables and effective rates.

---

## Option 3: nine-field strip

The most physically explicit strip model retains

\[
G_a,\qquad G_{mi},\qquad G_{ci}
\]

for each of the three GTPases.

The membrane equations are

\[
\frac{\partial G_a}{\partial t}
=
D_MG_{a,xx}
-
d_GG_a
+
Q_G\frac{G_{mi}}{G_{\rm tot}},
\]

\[
\frac{\partial G_{mi}}{\partial t}
=
D_MG_{mi,xx}
+
d_GG_a
-
Q_G\frac{G_{mi}}{G_{\rm tot}}
+
k_{\rm on}\ell_I G_{ci}
-
k_{\rm off}G_{mi},
\]

and the adjacent cytoplasmic strip satisfies

\[
\frac{\partial G_{ci}}{\partial t}
=
D_CG_{ci,xx}
-
\frac{1}{h}
\left(
k_{\rm on}\ell_I G_{ci}
-
k_{\rm off}G_{mi}
\right).
\]

For each GTPase,

\[
M_G
=
s_m
\int_0^L
\left(
G_a+G_{mi}+hG_{ci}
\right)dx
\]

up to the precise units chosen for membrane densities.

This is the natural 1D reduction of the Spill bulk–surface formulation.

For a size-scaling study, this version is the most mechanistically interpretable because \(k_{\rm on}\) and \(k_{\rm off}\) remain molecular exchange parameters while the geometry enters separately through \(h\).

---

# 7. Useful nondimensionalization

The Holmes parameterization uses approximately equal inactivation rates for the three GTPases,

\[
d_C=d_R=d_\rho=d\approx1~\mathrm{s^{-1}},
\]

and a common active diffusion coefficient

\[
D_m\approx0.1~\mu\mathrm{m}^2\,\mathrm{s^{-1}}.
\]

Choose

\[
\tau=dt
\]

and

\[
\ell_G
=
\sqrt{\frac{D_m}{d}}.
\]

For the Holmes baseline,

\[
\ell_G
=
\sqrt{0.1}
\approx0.316~\mu\mathrm{m}.
\]

Define

\[
X=\frac{x}{\ell_G}.
\]

The dimensionless cell length is

\[
\Lambda
=
L\sqrt{\frac{d}{D_m}}.
\]

For

\[
L=20~\mu\mathrm{m},
\]

\[
\boxed{
\Lambda\approx63.2.
}
\]

---

## Diffusion ratio

Define

\[
\mathcal D
=
\frac{D_i}{D_m}.
\]

For the Holmes/Marée baseline,

\[
\boxed{
\mathcal D=500.
}
\]

This strong transport separation is one of the best-supported structural features of the model class.

---

## Dimensionless activation strengths

Normalize each GTPase by its own total concentration,

\[
c=\frac{C}{C_{\rm tot}},
\qquad
r=\frac{R}{R_{\rm tot}},
\qquad
q=\frac{\rho}{\rho_{\rm tot}}.
\]

For the Holmes network define

\[
\kappa_C
=
\frac{I_C}{d_CC_{\rm tot}},
\]

\[
\kappa_R
=
\frac{I_{R1}}{d_RR_{\rm tot}},
\]

\[
\kappa_\rho
=
\frac{I_\rho}{d_\rho\rho_{\rm tot}}.
\]

Using the published baseline,

\[
\boxed{
\kappa_C\approx1.23,
}
\]

\[
\boxed{
\kappa_R\approx0.0267,
}
\]

\[
\boxed{
\kappa_\rho\approx2.13.
}
\]

These numbers make clear that the apparently large differences among the dimensional activation rates partly reflect differences in the total concentrations assigned to the three GTPases.

---

## Cdc42-to-Rac coupling

Define

\[
\chi_{CR}
=
\frac{\alpha C_{\rm tot}}
{d_RR_{\rm tot}}.
\]

For the Holmes baseline,

\[
\boxed{
\chi_{CR}
\approx0.208.
}
\]

---

## Inhibitory thresholds

The dimensionless Rho level for half-maximal inhibition of Cdc42 is

\[
\theta_{\rho C}
=
\frac{a_1}{\rho_{\rm tot}}.
\]

For Holmes,

\[
\boxed{
\theta_{\rho C}
\approx0.403.
}
\]

The dimensionless Rac level for half-maximal inhibition of Rho is

\[
\theta_{R\rho}
=
\frac{a_2}{R_{\rm tot}},
\]

giving

\[
\boxed{
\theta_{R\rho}
\approx0.133.
}
\]

Thus Rho activation is strongly sensitive to Rac occupying even a modest fraction of the available Rac pool.

---

## Cooperativity

The core inhibitory Hill coefficient is

\[
\boxed{n=3}
\]

in the canonical Holmes/Marée parameterizations.

Jilkine et al. show why \(n>1\) matters dynamically: cooperative inhibition creates the nonlinear threshold structure needed to combine a stable resting state with rapid polarization.

A practical exploration range is therefore

\[
\boxed{n=2-4}
\]

with \(n=3\) as the baseline.

---

## Stimulus strength

For Rac stimulation define

\[
\Sigma_R(x,t)
=
\frac{S_R(x,t)}
{d_RR_{\rm tot}}.
\]

This is preferable to comparing dimensional stimulus amplitudes across parameterizations.

The mean and spatial contrast can then be varied independently:

\[
\Sigma_R(x,t)
=
\Sigma_0(t)
+
\Sigma_g(t)
\left(
2x/L-1
\right).
\]

Holmes et al. explicitly found that stimulus mean and gradient can play distinct roles: the mean shifts the cell's overall sensitivity while the gradient provides the spatial asymmetry that selects the polarization direction.

---

# 8. Geometry and abundance parameters

A single-species wave-pinning model has one conserved abundance parameter. A Rac–Rho–Cdc42 model has three.

For a line-density model define

\[
\mu_C
=
\frac{M_C}{LC_{\rm ref}},
\]

\[
\mu_R
=
\frac{M_R}{LR_{\rm ref}},
\]

\[
\mu_\rho
=
\frac{M_\rho}{L\rho_{\rm ref}}.
\]

For a membrane–cytoplasm strip, the corresponding amount per unit membrane area is

\[
m_G
=
\ell_m\bar G+h\bar G_i.
\]

Thus for each GTPase the geometry enters through

\[
\boxed{h=V/A_{\rm mem}}
\]

or, more precisely, through

\[
\boxed{\ell_m/h}.
\]

For geometrically similar enlargement,

\[
h\propto L_{\rm char}.
\]

At fixed cytoplasmic concentration, a larger cell therefore contains a larger inactive reservoir per unit membrane area.

This is likely to be one of the most important size-dependent control variables in the three-GTPase model.

---

# 9. What should scale under a WGD / cell-size null?

For an isometric concentration-preserving enlargement, the clean null is:

### Fixed molecular parameters

Keep unchanged

\[
D_m,\quad
D_i,\quad
d_G,\quad
I_G,
\quad
a_1,\quad
a_2,\quad
n,\quad
\alpha,
\]

and, in the explicit exchange model,

\[
k_{\rm on},\quad k_{\rm off}.
\]

These are treated as local molecular kinetic or transport properties.

### Fixed concentrations

Keep the baseline concentrations of Cdc42, Rac and Rho approximately unchanged:

\[
C_{\rm tot},\quad
R_{\rm tot},\quad
\rho_{\rm tot}.
\]

Total cellular molecule number then increases with cell volume.

### Change geometry

Increase

\[
L
\]

and, in the physically explicit strip,

\[
h=V/A_{\rm mem}.
\]

For geometrically similar 3D cells scaled by a linear factor \(s\),

\[
L\rightarrow sL,
\]

\[
A\rightarrow s^2A,
\]

\[
V\rightarrow s^3V,
\]

so

\[
h\rightarrow sh.
\]

For a twofold increase in cell volume,

\[
s=2^{1/3}\approx1.26.
\]

Thus both

\[
\Lambda=L\sqrt{d/D_m}
\]

and

\[
h
\]

increase by approximately \(26\%\).

This is the cleanest geometric null.

---

# 10. Supported parameter space for the strip model

The literature supports a canonical **region** more strongly than a unique parameter set.

## Transport

A very well-established baseline in the Holmes/Marée lineage is

\[
\boxed{
D_m\approx0.1~\mu\mathrm{m}^2\,\mathrm{s^{-1}}
}
\]

and

\[
\boxed{
D_i\approx50~\mu\mathrm{m}^2\,\mathrm{s^{-1}}.
}
\]

The corresponding ratio is

\[
\boxed{
D_i/D_m\approx500.
}
\]

For sensitivity analysis, the ratio itself is more important than separately tuning both diffusivities.

A reasonable exploratory range is approximately

\[
\boxed{
D_i/D_m=10^2-10^3,
}
\]

while retaining \(500\) as the literature-centered baseline.

---

## Turnover

The canonical active-GTPase inactivation rate is

\[
\boxed{
d_C=d_R=d_\rho\approx1~\mathrm{s^{-1}}.
}
\]

This should be the default reaction timescale.

For model sensitivity rather than empirical uncertainty, a logarithmic perturbation around this baseline is more interpretable than treating \(0.01-10~\mathrm{s^{-1}}\) as an equally supported biological interval.

---

## Cooperativity

Use

\[
\boxed{n=3}
\]

as baseline and

\[
\boxed{n=2-4}
\]

for the main exploration.

---

## Conserved concentrations

For the Holmes baseline,

\[
\boxed{
(C_{\rm tot},R_{\rm tot},\rho_{\rm tot})
=
(2.4,7.5,3.1)~\mu\mathrm{M}.
}
\]

These should initially be treated as three independent abundance controls rather than collapsed to one “total GTPase” parameter.

---

## Domain size

The Holmes baseline is

\[
\boxed{
L=20~\mu\mathrm{m}.
}
\]

The more useful dimensionless quantity is

\[
\boxed{
\Lambda
=
L\sqrt{d/D_m}.
}
\]

For the baseline,

\[
\Lambda\approx63.
\]

Cell-size scaling should therefore be implemented by varying \(\Lambda\) through the physical geometry, not by altering reaction kinetics to compensate for length.

---

# 11. Recommended baseline model

For the present membrane–cytoplasm strip, the most useful first implementation is the **Holmes six-field GTPase-only model**:

\[
C_t
=
D_mC_{xx}
+
\frac{I_C}{1+(\rho/a_1)^n}
\frac{C_i}{C_{\rm tot}}
-
d_CC,
\]

\[
(C_i)_t
=
D_i(C_i)_{xx}
-
\frac{I_C}{1+(\rho/a_1)^n}
\frac{C_i}{C_{\rm tot}}
+
d_CC,
\]

\[
R_t
=
D_mR_{xx}
+
\left[
I_{R1}
+
\alpha C
+
S_R(x,t)
\right]
\frac{R_i}{R_{\rm tot}}
-
d_RR,
\]

\[
(R_i)_t
=
D_i(R_i)_{xx}
-
\left[
I_{R1}
+
\alpha C
+
S_R(x,t)
\right]
\frac{R_i}{R_{\rm tot}}
+
d_RR,
\]

\[
\rho_t
=
D_m\rho_{xx}
+
\frac{I_\rho}{1+(R/a_2)^n}
\frac{\rho_i}{\rho_{\rm tot}}
-
d_\rho\rho,
\]

\[
(\rho_i)_t
=
D_i(\rho_i)_{xx}
-
\frac{I_\rho}{1+(R/a_2)^n}
\frac{\rho_i}{\rho_{\rm tot}}
+
d_\rho\rho.
\]

Use no-flux boundary conditions,

\[
\partial_xC
=
\partial_xC_i
=
\partial_xR
=
\partial_xR_i
=
\partial_x\rho
=
\partial_x\rho_i
=
0
\]

at \(x=0,L\).

Baseline parameters:

| Parameter | Baseline |
|---|---:|
| \(L\) | \(20~\mu\mathrm m\) |
| \(D_m\) | \(0.1~\mu\mathrm{m}^2/s\) |
| \(D_i\) | \(50~\mu\mathrm{m}^2/s\) |
| \(C_{\rm tot}\) | \(2.4~\mu\mathrm M\) |
| \(R_{\rm tot}\) | \(7.5~\mu\mathrm M\) |
| \(\rho_{\rm tot}\) | \(3.1~\mu\mathrm M\) |
| \(I_C\) | \(2.95~\mu\mathrm{M}/s\) |
| \(I_{R1}\) | \(0.2~\mu\mathrm{M}/s\) |
| \(I_\rho\) | \(6.6~\mu\mathrm{M}/s\) |
| \(a_1\) | \(1.25~\mu\mathrm M\) |
| \(a_2\) | \(1~\mu\mathrm M\) |
| \(n\) | \(3\) |
| \(\alpha\) | \(0.65~s^{-1}\) |
| \(d_C,d_R,d_\rho\) | \(1~s^{-1}\) |

This is deliberately the GTPase-only core. PI feedback can be added later through \(I_{R2}\) and the PIP equations.

---

# 12. Recommended parameter sweep

The model nominally has many dimensional parameters, but the main dynamical questions can be organized around a much smaller set of combinations:

\[
\boxed{
\mathcal D
=
\frac{D_i}{D_m},
}
\]

\[
\boxed{
\Lambda
=
L\sqrt{\frac{d}{D_m}},
}
\]

\[
\boxed{
\kappa_C,\kappa_R,\kappa_\rho,
}
\]

\[
\boxed{
\chi_{CR},
}
\]

\[
\boxed{
\theta_{\rho C},
\theta_{R\rho},
}
\]

\[
\boxed{n,}
\]

together with the three conserved abundances and, for the physical strip,

\[
\boxed{h\ \text{or}\ \ell_m/h.}
\]

A useful hierarchy is:

### Tier 1: establish the basic dynamical regime

Sweep

\[
\Lambda,
\qquad
\kappa_R,
\qquad
n,
\]

while keeping the remaining Holmes parameters fixed.

This determines whether the model is

* insensitive/unpolarized;
* polarizable only above a finite perturbation;
* spontaneously unstable;
* persistently polarized.

### Tier 2: front–rear balance

Sweep

\[
\theta_{\rho C},
\qquad
\theta_{R\rho},
\qquad
\chi_{CR}.
\]

These tune the strength of the Cdc42–Rac front module relative to the Rho rear module.

### Tier 3: conserved-pool effects

Sweep the three GTPase abundances independently:

\[
M_C,\qquad M_R,\qquad M_\rho.
\]

This is the natural analogue of the total-mass sweep in a one-species wave-pinning model.

### Tier 4: size scaling

At a parameter point that supports robust polarization, vary

\[
L
\]

and

\[
h
\]

according to the geometric scaling rule while holding molecular kinetics and concentrations fixed.

This directly tests whether enlargement changes:

* polarization existence;
* front width;
* front/rear contrast;
* polarization onset time;
* susceptibility to reversal;
* recovery after removal of a stimulus;
* coexistence or competition of two fronts.

---

# 13. Polarization criteria

A three-GTPase model needs a clearer polarization definition than a single activator model.

For active Cdc42 and Rac, define front–back contrasts

\[
P_C
=
\frac{C_{\max}-C_{\min}}
{C_{\max}+C_{\min}},
\]

\[
P_R
=
\frac{R_{\max}-R_{\min}}
{R_{\max}+R_{\min}}.
\]

For Rho, whose high-activity region should be opposite the Cdc42/Rac front, define

\[
P_\rho
=
\frac{\rho_{\max}-\rho_{\min}}
{\rho_{\max}+\rho_{\min}}.
\]

A useful composite polarity condition is:

1. \(P_C\) and \(P_R\) exceed a threshold;
2. Cdc42 and Rac maxima are colocated;
3. the Rho maximum is spatially separated from the Cdc42/Rac maximum;
4. the pattern persists after a transient initiating stimulus is removed.

For example,

\[
\mathrm{corr}(C,R)>0
\]

and

\[
\mathrm{corr}(C,\rho)<0.
\]

A directional polarity vector can be defined from the first spatial moment:

\[
p_C
=
\frac{
\int_0^L
(x-L/2)C(x)\,dx
}{
(L/2)\int_0^L C(x)\,dx
},
\]

with analogous definitions for Rac and Rho.

A coherent polarized state should satisfy approximately

\[
\operatorname{sign}(p_C)
=
\operatorname{sign}(p_R)
=
-\operatorname{sign}(p_\rho).
\]

---

# 14. Dual polarization and multiple fronts

The three-GTPase system creates additional possible states beyond ordinary one-front wave pinning.

A dual-front state would have two separated Cdc42/Rac-rich regions with Rho enriched between or opposite them.

Useful diagnostics are:

* number of local maxima of \(C\) and \(R\);
* peak separation;
* lifetime of the weaker peak;
* integrated active mass associated with each peak;
* time to winner-take-all competition;
* whether the state survives after the initiating stimulus is removed.

The canonical literature gives two competing expectations:

1. mass conservation and fast inactive diffusion encourage long-range competition between active zones;
2. additional feedback layers, especially phosphoinositides, can substantially accelerate competition and eliminate multiple fronts.

Thus a GTPase-only model may support long-lived dual fronts in parameter regimes where the full GTPase–PI model resolves them more rapidly.

For the size-scaling question, this is particularly interesting because enlargement increases both the physical separation between fronts and the size of the inactive cytoplasmic reservoir.

---

# 15. Recommended analysis strategy

The most efficient workflow is not a brute-force sweep over every dimensional parameter.

## Step 1: reproduce a canonical polarized state

Implement the Holmes six-field model at

\[
L=20~\mu\mathrm m
\]

with its native parameters.

Use a transient localized or graded Rac stimulus and verify

* Cdc42/Rac enrichment at the stimulated front;
* Rho enrichment at the opposite side;
* persistence after stimulus removal.

## Step 2: local perturbation analysis / bifurcation analysis

Holmes et al. introduced Local Perturbation Analysis specifically for this class of mass-conserved fast–slow reaction–diffusion models.

The useful bifurcation parameter is initially

\[
I_{R1}
\]

or its dimensionless equivalent

\[
\kappa_R.
\]

The expected regimes are

* stable insensitive homogeneous state;
* stable but finitely polarizable state;
* Turing/hypersensitive regime;
* overstimulated homogeneous high-activity state.

For migration-like polarity, the **finitely polarizable** regime is the most relevant because it has a stable resting state but can break symmetry after a sufficiently strong local cue.

## Step 3: continuation in cell size

Once a polarizable baseline is identified, continue the polarized branch in

\[
\Lambda.
\]

This directly answers whether increasing size

* destroys the polarized branch;
* stabilizes it;
* introduces additional patterned branches;
* changes hysteresis or reversal thresholds.

## Step 4: two-parameter diagrams

The most informative first diagrams are likely

\[
(\Lambda,\kappa_R),
\]

\[
(\Lambda,\theta_{R\rho}),
\]

\[
(\Lambda,\theta_{\rho C}),
\]

and

\[
(\Lambda,\mu_G).
\]

These ask whether cell size interacts with overall activation, inhibitory crosstalk or conserved abundance.

## Step 5: explicit geometry

Only after understanding the six-field model should the physical strip parameter

\[
h
\]

or the explicit membrane-binding kinetics be introduced.

Then compare:

1. changing \(L\) alone;
2. changing \(h\) alone;
3. geometric scaling in which both change together.

This separates axial-distance effects from membrane-to-volume effects.

---

# 16. Recommended size-scaling experiments

For a \(1\times\) versus \(2\times\)-volume comparison, let

\[
s=2^{1/3}\approx1.26.
\]

Then use

\[
L_{2\times}=sL_{1\times}.
\]

If the strip is intended to preserve a geometrically similar membrane/bulk relationship, also use

\[
h_{2\times}=sh_{1\times}.
\]

Run the following matched simulations.

### A. Spontaneous/noise perturbation

Start near the homogeneous state with small spatial noise.

Measure whether either size spontaneously polarizes.

This identifies entry into a hypersensitive/Turing-like regime.

### B. Standardized local Rac pulse

Apply the same **local concentration-level** Rac activation stimulus over a fixed fractional region of the cell.

Measure

* polarization probability;
* onset time;
* final contrast;
* front width;
* persistence after stimulus withdrawal.

### C. Competing opposite stimuli

Apply two opposite Rac pulses of equal amplitude.

Measure

* whether one front wins;
* time to symmetry breaking;
* frequency of long-lived dual fronts.

### D. Slightly asymmetric competition

Apply two opposite pulses with a small amplitude difference.

Measure the minimum asymmetry required to select the stronger side.

This is an especially useful test of whether larger cells are better or worse at global front competition.

### E. Reversal

Establish a polarized state, then apply a Rac stimulus at the opposite end.

Measure the amplitude–duration reversal boundary.

This can be represented as a phase diagram

\[
A_{\rm reverse}
\quad\text{versus}\quad
T_{\rm pulse}.
\]

Compare the boundary between \(1\times\) and enlarged cells.

### F. Rapidly alternating cues

Alternate the favored Rac-activation side at a defined period.

Measure the transfer function between stimulus orientation and polarity orientation.

This directly connects the biochemical model to experiments involving rapidly reorienting external polarity cues.

---

# 17. Interpretation

The most canonical Rac–Rho–Cdc42 polarity literature does not support a literal three-field model. It supports a hierarchy:

\[
\boxed{
3\ \text{GTPases}
\rightarrow
6\ \text{active/inactive fields}
\rightarrow
9\ \text{fields if membrane/cytosol exchange is explicit}.
}
\]

For a first model, the Holmes 2012 six-field formulation is the strongest default because it is

* mass conserving;
* explicitly Rac–Rho–Cdc42;
* derived for a nearly 1D mammalian-cell geometry;
* experimentally motivated by direct Rac stimulation;
* equipped with a published parameter set;
* already analyzed by local perturbation and bifurcation methods.

Its canonical dimensional regime is approximately

\[
\boxed{
D_m=0.1~\mu\mathrm{m}^2/s,
\qquad
D_i=50~\mu\mathrm{m}^2/s,
\qquad
d_G=1~s^{-1},
\qquad
n=3,
\qquad
L=20~\mu m.
}
\]

The most useful corresponding dimensionless quantities are

\[
\boxed{
\mathcal D=D_i/D_m\approx500,
}
\]

\[
\boxed{
\Lambda=L\sqrt{d/D_m}\approx63,
}
\]

together with

\[
\boxed{
(\kappa_C,\kappa_R,\kappa_\rho)
}
\]

for activation strength,

\[
\boxed{
(\theta_{\rho C},\theta_{R\rho})
}
\]

for inhibitory thresholds,

\[
\boxed{
\chi_{CR}
}
\]

for Cdc42-to-Rac coupling, and the three conserved GTPase abundances.

For cell-size scaling, the crucial addition is the geometric reservoir parameter

\[
\boxed{
h=\frac{V}{A_{\rm mem}},
}
\]

or, in an explicit exchange model, the combination of \(h\) with \(k_{\rm on}\) and \(k_{\rm off}\).

The recommended progression is therefore

\[
\boxed{
\text{Holmes 6-field native model}
\rightarrow
\text{map polarization regimes}
\rightarrow
\text{vary }\Lambda
\rightarrow
\text{introduce }h
\rightarrow
\text{upgrade to explicit 9-field exchange if needed}.
}
\]

This preserves direct comparability with the canonical literature while isolating the geometry-dependent effects that are central to the cell-size problem.

---

# References

* Marée AFM, Jilkine A, Dawes A, Grieneisen VA, Edelstein-Keshet L. *Polarization and movement of keratocytes: a multiscale modelling approach.* Bulletin of Mathematical Biology. 2006;68:1169–1211. DOI: 10.1007/s11538-006-9131-7. https://doi.org/10.1007/s11538-006-9131-7
* Jilkine A, Marée AFM, Edelstein-Keshet L. *Mathematical model for spatial segregation of the Rho-family GTPases based on inhibitory crosstalk.* Bulletin of Mathematical Biology. 2007;69:1943–1978. DOI: 10.1007/s11538-007-9200-6. https://doi.org/10.1007/s11538-007-9200-6
* Holmes WR, Lin B, Levchenko A, Edelstein-Keshet L. *Modelling Cell Polarization Driven by Synthetic Spatially Graded Rac Activation.* PLoS Computational Biology. 2012;8:e1002366. DOI: 10.1371/journal.pcbi.1002366. https://doi.org/10.1371/journal.pcbi.1002366
* Marée AFM, Grieneisen VA, Edelstein-Keshet L. *How Cells Integrate Complex Stimuli: The Effect of Feedback from Phosphoinositides and Cell Shape on Cell Polarization and Motility.* PLoS Computational Biology. 2012;8:e1002402. DOI: 10.1371/journal.pcbi.1002402. https://doi.org/10.1371/journal.pcbi.1002402
* Spill F, Andasari V, Mak M, Kamm RD, Zaman MH. *Effects of 3D geometries on cellular gradient sensing and polarization.* Physical Biology. 2016;13:036008. DOI: 10.1088/1478-3975/13/3/036008. https://doi.org/10.1088/1478-3975/13/3/036008
* Comer S et al. *On the influence of cell shape on dynamic reaction-diffusion polarization patterns.* PLoS ONE. 2021;16:e0248293. DOI: 10.1371/journal.pone.0248293. https://doi.org/10.1371/journal.pone.0248293

