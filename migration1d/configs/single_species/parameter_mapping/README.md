# Wave-pinning parameterization for a coupled membrane–cytoplasm strip

## Geometry

We consider a one-dimensional representation of a cell of length (L), consisting of two adjacent compartments:

* a membrane strip running from (x=0) to (x=L);
* an adjacent cytoplasmic strip over the same axial coordinate.

The active GTPase, (u(x,t)), is membrane associated and diffuses slowly along the membrane. The inactive GTPase, (v(x,t)), is cytoplasmic and diffuses more rapidly along the bulk strip.

For a physically explicit reduction of a bulk–surface cell geometry, define

* (s_m): membrane area represented per unit axial length ([\mu\mathrm{m}]);
* (A_c): cytoplasmic volume represented per unit axial length ([\mu\mathrm{m}^2]).

The effective cytoplasmic depth associated with the modeled membrane is then

[
h = \frac{A_c}{s_m}.
]

If (u) is a membrane surface concentration and (v) is a cytoplasmic volume concentration, the total amount of GTPase represented by the strip is

[
M =
\int_0^L
\left(
s_m u + A_c v
\right),dx.
]

A membrane reaction flux (J(u,v)), with units of amount per membrane area per time, gives

[
\frac{\partial u}{\partial t}
=============================

D_u\frac{\partial^2u}{\partial x^2}
+
J(u,v),
]

[
\frac{\partial v}{\partial t}
=============================

## D_v\frac{\partial^2v}{\partial x^2}

\frac{1}{h}J(u,v).
]

The factor (1/h=s_m/A_c) is important: activation of a given membrane area removes material from the adjacent cytoplasmic volume. The resulting system conserves (M).

This geometry is therefore intermediate between the original one-dimensional wave-pinning models, in which both variables are represented on the same line, and a full bulk–surface model such as Giese et al., in which active GTPase occupies the membrane and inactive GTPase occupies the cell volume.

---

## Common wave-pinning kinetics

The Mori/Giese form of the wave-pinning reaction is

[
J(u,v)
======

v
\left(
k_0+
\gamma\frac{u^n}{K^n+u^n}
\right)
-\delta u.
]

The parameters are:

| Parameter   | Meaning                                                           |
| ----------- | ----------------------------------------------------------------- |
| (D_u)       | diffusion coefficient of active/membrane GTPase                   |
| (D_v)       | diffusion coefficient of inactive/cytoplasmic GTPase              |
| (k_0)       | basal activation rate                                             |
| (\gamma)    | maximal positive-feedback activation rate                         |
| (K)         | active-GTPase level giving half-maximal feedback                  |
| (\delta)    | first-order inactivation rate                                     |
| (n)         | Hill coefficient                                                  |
| (M)         | total conserved GTPase                                            |
| (L)         | cell/domain length                                                |
| (h=A_c/s_m) | effective cytoplasmic volume per membrane area in the strip model |

Mori and Buttenschön use closely related forms. Chiou uses a slightly different but mechanistically equivalent saturating positive-feedback function.

---

## Published parameterizations

### Mori et al. 2008

Mori et al. introduced the canonical one-dimensional wave-pinning model. Active and inactive species are both represented as one-dimensional concentrations on a domain of length (L), so the model does not explicitly distinguish membrane area from cytoplasmic volume.

Their canonical parameterization is approximately

| Parameter        |                             Mori 2008 |
| ---------------- | ------------------------------------: |
| (L)              |                    (10~\mu\mathrm{m}) |
| (D_u)            | (0.1~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| (D_v)            |  (10~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| (k_0)            |               (0.067~\mathrm{s}^{-1}) |
| (\gamma)         |                   (1~\mathrm{s}^{-1}) |
| (\delta)         |                   (1~\mathrm{s}^{-1}) |
| (K)              |  (1), in their 1D concentration units |
| (n)              |                                     2 |
| initial (u)      |                             0.2683312 |
| initial (v)      |                                   2.0 |
| total amount (M) |                  (22.6833) for (L=10) |

The initial state satisfies

[
M=L(u_0+v_0)
=10(0.2683312+2)
=22.6833.
]

The key transport ratio is

[
\frac{D_v}{D_u}=100.
]

Mori explicitly describes the concentration units as effectively "amount per length"; therefore (K), (u), (v), and (M) should not be transferred literally to a membrane–volume geometry. The reaction timescales and diffusion coefficients are more directly portable. ([PubMed Central (PMC)][1])

### Chiou et al. 2018

Chiou et al. study the generic properties of mass-conserved activator–substrate models. Their wave-pinning example uses

[
J(u,v)
======

\frac{a u^2}{1+k u^2}v-bu.
]

Representative parameters include

| Parameter |                                                     Chiou 2018 |
| --------- | -------------------------------------------------------------: |
| (D_u)     |                         (0.01~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| (D_v)     | (1~\mu\mathrm{m}^2,\mathrm{s}^{-1}), or (D_v\rightarrow\infty) |
| (a)       |                            (1~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| (b)       |                                            (1~\mathrm{s}^{-1}) |
| (k)       |            (0.01)–(1~\mu\mathrm{m}^2), depending on experiment |
| (L)       |                               commonly (10)–(20~\mu\mathrm{m}) |
| (M)       |                            examples include 30, 40, 50 and 200 |

Their systematic parameter analysis uses a basal set

[
a=1,\qquad
b=1,\qquad
k=0.01,\qquad
D_u=0.01,\qquad
M=40,\qquad
L=20,
]

and varies individual parameters over approximately (0.2\times) to (3\times) their basal values. With finite cytoplasmic diffusion they use (D_v=1~\mu\mathrm{m}^2,\mathrm{s}^{-1}). ([PLOS][2])

The natural feedback concentration scale for the Chiou kinetics is

[
K_{\mathrm{eff}}=\frac{1}{\sqrt{k}},
]

because

[
1+k u^2
=======

1+\left(\frac{u}{K_{\mathrm{eff}}}\right)^2.
]

The corresponding dimensionless maximal feedback strength is

[
\Gamma_{\mathrm{Chiou}}
=======================

\frac{a}{bk}.
]

For the basal (k=0.01) parameterization,

[
\Gamma_{\mathrm{Chiou}}=100.
]

This large number should not be compared directly with (\gamma/\delta) in Mori without accounting for the different concentration normalization.

### Buttenschön et al. 2020

Buttenschön et al. use a Mori-like GTPase module on a changing one-dimensional domain,

[
J(G,G_i)
========

\left(
\beta+
\gamma\frac{G^n}{G_0^n+G^n}
\right)G_i
----------

IG.
]

Their literature-motivated parameter table gives

| Parameter                    |                 Supported range/value |
| ---------------------------- | ------------------------------------: |
| basal activation (\beta)     |          (0.02)–(0.1~\mathrm{s}^{-1}) |
| feedback activation (\gamma) |              (1)–(25~\mathrm{s}^{-1}) |
| inactivation (I)             |                   (1~\mathrm{s}^{-1}) |
| half-activation (G_0)        |        (1), in 1D concentration units |
| active diffusion (D)         | (0.1~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| inactive diffusion (D_i)     | (100~\mu\mathrm{m}^2,\mathrm{s}^{-1}) |
| total GTPase (G_T)           |         (1)–(100), in 1D amount units |
| Hill coefficient (n)         |                                  3–25 |

Thus their biologically motivated diffusion ratio is

[
\frac{D_i}{D}=1000.
]

Some of their actual numerical experiments use less extreme transport parameters. For example, their domain-size experiment uses

[
L=55~\mu\mathrm{m},\quad
D=0.01,\quad
D_i=1,\quad
\gamma=2,\quad
\beta=0.1,\quad
I=1,\quad
G_T=2.5,\quad
n=4,
]

giving

[
D_i/D=100.
]

A major feature of this model is that (G_T) is an amount. If (G_T) is fixed while (L) changes, the mean GTPase concentration scales approximately as

[
\bar G_T = \frac{G_T}{L},
]

so growth produces dilution and shrinkage produces concentration. ([ResearchGate][3])

### Giese et al. 2015

Giese et al. provide the closest conceptual match to the membrane-plus-cytoplasm strip geometry because they explicitly extend wave-pinning to a bulk–surface system. Active GTPase occupies the membrane and inactive GTPase occupies the cytoplasmic domain. ([PubMed][4])

They use the Mori reaction function

[
J(u,v)
======

v\left(
k_0+
\gamma\frac{u^2}{K^2+u^2}
\right)
-------

\delta u
]

with

[
k_0=0.067,\qquad
\gamma=1,\qquad
\delta=1,\qquad
K=1,\qquad
n=2.
]

Their reported WP initial condition is

[
u_0=0.2,\qquad v_0=1.95.
]

These values represent the stable low-activity homogeneous state used in their simulations. ([researchgate.net][5])

The important conceptual change from Mori is conservation over two geometrically different compartments:

[
M
=

\int_{\Gamma}u,dA
+
\int_{\Omega}v,dV.
]

In the strip geometry this becomes

[
M
=

\int_0^L
(s_m u+A_c v),dx
================

s_m
\int_0^L
(u+h v),dx.
]

Thus the effective cytoplasmic depth (h=A_c/s_m) becomes an additional control parameter even when all biochemical concentrations and rate constants are unchanged.

---

## Mapping the published models onto the strip geometry

There are two possible levels of abstraction.

### Option 1: purely one-dimensional wave-pinning

The simplest implementation treats both variables as linear densities:

[
U(x,t)=\text{active amount per unit length},
]

[
V(x,t)=\text{inactive amount per unit length}.
]

Then

[
\frac{\partial U}{\partial t}
=============================

D_u U_{xx}+F(U,V),
]

[
\frac{\partial V}{\partial t}
=============================

D_v V_{xx}-F(U,V),
]

with

[
M=\int_0^L(U+V),dx.
]

This corresponds most directly to Mori, Chiou and Buttenschön. It is useful for studying wave-pinning itself but does not represent the changing membrane-to-volume ratio of a real cell.

### Option 2: membrane–cytoplasm strip

For cell-size scaling, the preferred representation is

[
u(x,t)=\text{membrane surface concentration},
]

[
v(x,t)=\text{cytoplasmic volume concentration}.
]

Then

[
u_t=D_u u_{xx}+J(u,v),
]

[
v_t=D_v v_{xx}-\frac{1}{h}J(u,v),
]

where

[
h=\frac{A_c}{s_m}.
]

The conserved amount per unit membrane area is

[
\frac{M}{s_mL}
==============

\bar u+h\bar v.
]

This is the one-dimensional strip analogue of the Giese bulk–surface model.

It also makes explicit a quantity that is hidden in standard 1D wave-pinning models:

[
\boxed{h = V/A_{\mathrm{mem}}}.
]

For geometrically similar cells, (h) scales with the characteristic cell size. Consequently, increasing cell size at fixed membrane and cytoplasmic concentrations increases the cytoplasmic reservoir available per unit membrane area.

---

## Useful nondimensionalization

For Mori-type kinetics, choose

[
U=\frac{u}{K},
]

[
\tau=\delta t,
]

and the active-GTPase diffusion length

[
\ell_u=\sqrt{\frac{D_u}{\delta}}.
]

Define

[
X=\frac{x}{\ell_u}.
]

The dimensionless cell length is then

[
\Lambda
=======

L\sqrt{\frac{\delta}{D_u}}.
]

Several useful dimensionless parameter combinations follow.

### Diffusion ratio

[
\rho
====

\frac{D_v}{D_u}.
]

Published wave-pinning implementations strongly support

[
\boxed{\rho\sim10^2-10^3}.
]

Mori and many numerical examples use approximately (100); Buttenschön's biologically motivated table extends this to (1000).

### Basal activation relative to inactivation

[
\kappa
======

\frac{k_0}{\delta}
]

for same-dimensional Mori-like kinetics.

Mori and Giese both use approximately

[
\boxed{\kappa\approx0.067}.
]

Buttenschön gives

[
\frac{\beta}{I}=0.02-0.1.
]

Thus a reasonable supported range is approximately

[
\boxed{\kappa\sim0.02-0.1}.
]

### Feedback relative to inactivation

For Mori/Buttenschön-type kinetics,

[
\Gamma
======

\frac{\gamma}{\delta}.
]

Mori uses

[
\Gamma=1,
]

whereas Buttenschön considers approximately

[
\Gamma=1-25.
]

A conservative starting range is therefore

[
\boxed{\Gamma\sim1-5},
]

with extension to

[
\boxed{\Gamma\sim25}
]

for broader sensitivity analysis.

Chiou demonstrates that substantially stronger effective feedback is also compatible with wave-pinning once the model is normalized differently, so exact values of (\Gamma) should not be compared across kinetic forms without normalizing (K).

### Hill coefficient

Mori and Giese use

[
n=2.
]

Buttenschön considers

[
n=3-25,
]

with (n=4) used in an important size-dependent polarization experiment.

For a general exploration,

[
\boxed{n=2-4}
]

is therefore a natural core range, while larger (n) can be treated as a sharp-switch limit rather than the default.

### Dimensionless domain size

[
\Lambda=L\sqrt{\frac{\delta}{D_u}}.
]

For example, Mori's canonical values give

[
\Lambda
=======

10\sqrt{\frac{1}{0.1}}
\approx31.6.
]

Using (D_u=0.01) instead gives

[
\Lambda=10L
]

when (L) is expressed in microns. Thus published models cover dimensionless domains from tens to several hundred.

Rather than assigning a universal value of (L), it is preferable to regard

[
\boxed{\Lambda=L/\ell_u}
]

as the relevant quantity: it measures how many active-GTPase diffusion lengths fit along the cell.

### Geometry parameter

For the membrane–cytoplasm strip, define a dimensionless cytoplasmic depth

[
H=\frac{h}{\ell_{\mathrm{exchange}}},
]

where (\ell_{\mathrm{exchange}}) is the length scale implicit in the membrane–cytoplasm reaction kinetics. Equivalently, if a specific bulk–surface kinetic formulation is chosen, the corresponding reaction parameters can be used to construct this length explicitly.

More generally, the important geometric quantity is simply

[
\boxed{h=V/A_{\mathrm{mem}}}
]

because total mass per membrane area is

[
m
=

\bar u+h\bar v.
]

This is the parameter through which cell size enters even when membrane and cytoplasmic concentrations remain unchanged.

---

## Supported parameter space for the strip model

The literature does not support a unique parameter set. It does, however, support a reasonably constrained region of parameter space.

### Transport

A useful core range is

[
\boxed{
D_u=0.01-0.1~\mu\mathrm{m}^2,\mathrm{s}^{-1}
}
]

and

[
\boxed{
D_v=1-100~\mu\mathrm{m}^2,\mathrm{s}^{-1}
}
]

with the additional constraint

[
\boxed{
D_v/D_u\approx100-1000.
}
]

The diffusion **ratio** is better supported across models than any individual value.

### Reaction timescale

Published models repeatedly use an inactivation timescale of approximately one second:

[
\boxed{
\delta\sim1~\mathrm{s}^{-1}.
}
]

Basal activation relative to this timescale is approximately

[
\boxed{
k_0/\delta\sim0.02-0.1.
}
]

Positive-feedback activation is typically comparable to or faster than inactivation:

[
\boxed{
\gamma/\delta\sim1-5
}
]

as a core range, with values up to approximately 25 represented in published sensitivity analyses.

### Cooperativity

[
\boxed{
n=2-4
}
]

is a well-supported core range.

### Cell length

Published models commonly use dimensional lengths of order

[
\boxed{
L\sim10-50~\mu\mathrm{m}
}
]

although the correct quantity for comparison is

[
\boxed{
\Lambda=L\sqrt{\delta/D_u}.
}
]

### GTPase abundance

Absolute GTPase abundance is the least transferable parameter because Mori, Chiou and Buttenschön use effective one-dimensional concentrations, whereas Giese separates surface and bulk concentrations.

For the strip model it is preferable to define the conserved abundance as

[
\boxed{
m=\frac{M}{s_mL}
================

\bar u+h\bar v
}
]

and explore abundance relative to the feedback concentration (K):

[
\boxed{
\mu
===

\frac{\bar u+h\bar v}{K}
}
]

after ensuring dimensional consistency between the surface and bulk variables.

The most defensible strategy is therefore not to import Mori's (M=22.68), Chiou's (M=40), or Buttenschön's (G_T=2.5) literally. Instead, each published parameterization should be reproduced in its own native units, its polarization regime identified, and the equivalent dimensionless abundance in the strip model chosen to place the system in the same regime.

---

## Recommended baseline and parameter sweep

For a membrane–cytoplasm strip implementation, a reasonable baseline inspired by the common region of the published models is

| Parameter              |                Suggested baseline |         Supported exploration |
| ---------------------- | --------------------------------: | ----------------------------: |
| (D_u)                  |          (0.05~\mu\mathrm{m}^2/s) |                    (0.01-0.1) |
| (D_v/D_u)              |                               100 |                    (100-1000) |
| (\delta)               |                        (1~s^{-1}) |           approximately fixed |
| (k_0/\delta)           |                             0.067 |                    (0.02-0.1) |
| (\gamma/\delta)        |                                 1 |       (1-5), optionally to 25 |
| (n)                    |                                 2 |                         (2-4) |
| (L)                    |                geometry dependent |   approximately (10-50~\mu m) |
| (h=V/A_{\mathrm{mem}}) |                geometry dependent | determined from cell geometry |
| abundance              | chosen within polarization regime |              systematic sweep |

The most useful sensitivity analysis is likely to vary the following **dimensionless combinations** rather than every dimensional parameter independently:

[
\boxed{
\rho=\frac{D_v}{D_u},
\qquad
\Lambda=L\sqrt{\frac{\delta}{D_u}},
\qquad
\kappa=\frac{k_0}{\delta},
\qquad
\Gamma=\frac{\gamma}{\delta},
\qquad
n,
\qquad
\mu,
\qquad
h\ \text{or its dimensionless equivalent}.
}
]

These combinations separate distinct biological questions:

* (\rho): separation of membrane and cytoplasmic transport timescales;
* (\Lambda): cell length relative to the intrinsic polarity length scale;
* (\kappa) and (\Gamma): strength of basal and feedback activation relative to GTPase inactivation;
* (n): cooperativity of feedback;
* (\mu): total available GTPase relative to the activation threshold;
* (h): cytoplasmic reservoir available per unit membrane.

For cell-size scaling, the last three quantities are especially important. Changing cell dimensions can alter (\Lambda), (h), and potentially (\mu) even if all molecular kinetic parameters remain unchanged.

## Interpretation

The four published models support a common qualitative parameter regime rather than a single canonical parameter set. Active GTPase diffuses slowly, inactive GTPase diffuses approximately (10^2)-(10^3) times faster, GTPase turnover occurs on an approximately second timescale, positive feedback is at least comparable to the inactivation rate, and cooperative activation typically has Hill coefficient (n\approx2-4).

Mori, Chiou and Buttenschön establish the kinetic and dynamical range over which wave-pinning occurs. Giese provides the most appropriate precedent for mapping those kinetics onto a membrane–cytoplasm geometry. In the strip representation, this introduces an explicit geometric control parameter,

[
h=\frac{V}{A_{\mathrm{mem}}},
]

which determines the size of the cytoplasmic reservoir available to each unit of membrane. Consequently, the most useful parameterization for subsequent size-scaling experiments is not simply a set of dimensional reaction constants, but the combination

[
\boxed{
(\rho,\Lambda,\kappa,\Gamma,n,\mu,h)
}
]

or a fully nondimensional equivalent. This representation permits published wave-pinning parameterizations to be compared on a common basis while retaining the membrane-to-cytoplasm geometry that is central to the cell-size problem.

## References

* Mori Y, Jilkine A, Edelstein-Keshet L. *Wave-pinning and cell polarity from a bistable reaction-diffusion system.* Biophysical Journal. 2008;94:3684–3697. DOI: 10.1529/biophysj.107.120824.
* Giese W, Eigel M, Westerheide S, Engwer C, Klipp E. *Influence of cell shape, inhomogeneities and diffusion barriers in cell polarization models.* Physical Biology. 2015;12:066014.
* Chiou J-G et al. *Principles that govern competition or co-existence in Rho-GTPase driven polarization.* PLoS Computational Biology. 2018;14:e1006095. DOI: 10.1371/journal.pcbi.1006095.
* Buttenschön A, Liu Y, Edelstein-Keshet L. *Cell size, mechanical tension, and GTPase signaling in the single cell.* Bulletin of Mathematical Biology. 2020;82:28. DOI: 10.1007/s11538-020-00702-5.

[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2292363/?utm_source=chatgpt.com "Wave-Pinning and Cell Polarity from a Bistable Reaction-Diffusion System - PMC"
[2]: https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1006095&utm_source=chatgpt.com "Principles that govern competition or co-existence in Rho-GTPase driven polarization | PLOS Computational Biology"
[3]: https://www.researchgate.net/publication/335462743_Cell_size_mechanical_tension_and_GTPase_signaling_in_the_Single_Cell?utm_source=chatgpt.com "(PDF) Cell size, mechanical tension, and GTPase signaling in the Single Cell"
[4]: https://pubmed.ncbi.nlm.nih.gov/26599916/?utm_source=chatgpt.com "Influence of cell shape, inhomogeneities and diffusion ..."
[5]: https://www.researchgate.net/publication/264911228_Influence_of_cell_shape_inhomogeneities_and_diffusion_barriers_in_cell_polarization_models?utm_source=chatgpt.com "(PDF) Influence of cell shape, inhomogeneities and diffusion barriers in cell polarization models"

