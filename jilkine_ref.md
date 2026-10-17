Yes — the PDF helps a lot. Jilkine gives essentially the full model equations, parameterizations, stimulation protocols, and the size sweep in the Methods. The main caveat is that **a few quantities needed for exact numerical reproduction are not fully specified in this review**, especially some initial conditions/total mass choices and the constant (E_c) treatment in the Goryachev model. 

## Common framework

For WP, OT, and GOR:

[
\frac{\partial u}{\partial t}
=============================

D_u\frac{\partial^2u}{\partial x^2}+f(u,v)
]

[
\frac{\partial v}{\partial t}
=============================

D_v\frac{\partial^2v}{\partial x^2}+g(u,v)
]

on

[
0\le x\le L
]

with **no-flux boundary conditions**. Default:

[
L=10;\mu{\rm m}.
]

For the mass-conserved systems,

[
g(u,v)=-f(u,v).
]

Jilkine calibrates all models to physical units of seconds and microns. 

---

# 1. Wave-pinning model

Variables:

* (u): active, membrane-associated Rho GTPase
* (v): inactive, rapidly diffusing cytosolic Rho GTPase

Reaction:

[
f(u,v)
======

v\left(
k_0+\frac{c u^2}{K^2+u^2}
\right)-du
]

and

[
g(u,v)=-f(u,v).
]

Parameters used by Jilkine:

| Parameter   |                          Value |
| ----------- | -----------------------------: |
| (D_u)       | (0.1;\mu{\rm m}^2{\rm s}^{-1}) |
| (D_v)       |  (10;\mu{\rm m}^2{\rm s}^{-1}) |
| (k_0)       |           (0.067;{\rm s}^{-1}) |
| (c)         |               (1;{\rm s}^{-1}) |
| (K)         |                 (1;\mu{\rm M}) |
| (d)         |               (1;{\rm s}^{-1}) |
| default (L) |                (10;\mu{\rm m}) |

The system conserves total (u+v). 

### Important missing item

The review **does not state the exact total GTPase abundance / homogeneous initial (u,v)** used for the Figure 2 simulations. That matters substantially for wave-pinning and especially for its size threshold. Jilkine explicitly notes that the critical size for transition to persistent two peaks depends on total material. 

So exact reproduction should recover that from Mori et al. 2008 or any original simulation code rather than infer it from the figure.

---

# 2. Otsuji model

Again (u,v) are active/inactive Rho-family proteins.

Jilkine writes:

[
f(u,v)
======

a_1\left[
v-
\frac{u+v}
{\left(a_2s(u+v)+1\right)^2}
\right]
]

with

[
g(u,v)=-f(u,v).
]

Parameters:

| Parameter |                          Value |
| --------- | -----------------------------: |
| (a_1)     |              (25;{\rm s}^{-1}) |
| (a_2)     |          (0.7;\mu{\rm M}^{-1}) |
| (s)       |                          (1.0) |
| (u+v)     |                 (2;\mu{\rm M}) |
| (D_u)     | (0.1;\mu{\rm m}^2{\rm s}^{-1}) |
| (D_v)     |  (10;\mu{\rm m}^2{\rm s}^{-1}) |
| (L)       |                (10;\mu{\rm m}) |

These are **not exactly the original Otsuji parameters**. Jilkine rescales the original model's time by a factor of 10 so that its diffusion coefficients match the WP model. With this parameterization, spatial modes (n=1) and (n=2) are both unstable. 

That is important: for WP1 we should reproduce **Jilkine's rescaled OT**, not blindly implement the original Otsuji 2007 numbers.

---

# 3. Goryachev reduced model

Variables:

* (u): GTP-bound Cdc42
* (v): GDP-bound Cdc42
* (E_c): cytoplasmic Cdc24–Bem1 complex

Reaction:

[
f(u,v)
======

aE_cu^2v+bE_cuv-cu
]

with

[
g(u,v)=-f(u,v).
]

The complete reduced formulation presented is

[
E_c=
\frac{E_c^0}
{1+\int_S \phi(u),ds}
]

where the integral represents conservation of total Cdc24–Bem1 complex.

However, **Jilkine explicitly says that for simplicity (E_c) was taken to be constant in their simulations.** 

Parameters:

[
D_u=0.0025;\mu{\rm m}^2{\rm s}^{-1}
]

[
D_v=10;\mu{\rm m}^2{\rm s}^{-1}
]

and

[
a=\frac{0.33}{1+\phi}
]

[
b=\frac{0.67}{1+\phi}
]

[
\phi=100
]

[
c=0.01733;{\rm s}^{-1}.
]

Numerically,

[
a\approx0.003267
]

and

[
b\approx0.006634.
]

### Important missing item

Jilkine does **not clearly provide the numerical constant (E_c)** they actually used after deciding to hold it fixed.

Nor does the review give an explicit homogeneous starting concentration for (u,v).

So exact GOR reproduction requires checking Goryachev & Pokhilko 2008 and/or recovering the original implementation.

This is probably the biggest unresolved detail in the four-model set.

---

# 4. LEGI model

Variables:

* (A): local activator
* (I): global inhibitor
* (R): response element
* (S(t,x)): external signal

Equations:

[
\frac{\partial A}{\partial t}
=============================

k_A S(t,x)-k_{-A}A
]

[
\frac{\partial I}{\partial t}
=============================

D_I\frac{\partial^2I}{\partial x^2}
+k_I S(t,x)-k_{-I}I
]

[
\frac{\partial R}{\partial t}
=============================

k_R A(R_T-R)-k_{-R}IR.
]

Parameters:

| Parameter |                           Value |
| --------- | ------------------------------: |
| (k_A)     |                (2,{\rm s}^{-1}) |
| (k_{-A})  |                (2,{\rm s}^{-1}) |
| (k_I)     |                (1,{\rm s}^{-1}) |
| (k_{-I})  |                (1,{\rm s}^{-1}) |
| (k_R)     | (1,\mu{\rm M}^{-1}{\rm s}^{-1}) |
| (k_{-R})  | (1,\mu{\rm M}^{-1}{\rm s}^{-1}) |
| (R_T)     |                  (1,\mu{\rm M}) |
| (D_I)     |   (10,\mu{\rm m}^2{\rm s}^{-1}) |

Jilkine notes that adding

[
D_A=D_R=0.1;\mu{\rm m}^2{\rm s}^{-1}
]

does not materially change results. Their main implementation effectively treats (A) and (R) as local. 

LEGI differs conceptually from the other three because it is not intrinsically pattern-forming: when the signal is removed, polarity disappears.

---

# 5. Common stimulus coupling for WP / OT / GOR

Jilkine modifies the active/inactive equations to:

[
\frac{\partial u}{\partial t}
=============================

D_u u_{xx}
+
f(u,v)
+
k_Sv
]

[
\frac{\partial v}{\partial t}
=============================

## D_v v_{xx}

## f(u,v)

k_Sv.
]

Thus, the external stimulus acts by converting inactive (v) into active (u), preserving total material.

For LEGI, the equivalent (k_S) profile is simply used as the external signal (S(t,x)). 

This is important for a common implementation API: every model can receive the same abstract stimulus profile even though LEGI consumes it differently.

---

# 6. Stimulus protocols

These are probably the most directly reusable part for WP1.

## A. Single transient localized stimulus

Stimulate the leftmost 10% of a default (L=10,\mu m) domain:

[
k_S^{\rm loc}(x,t)=
s(t)(1+\cos \pi x),
\qquad
0\le x\le1,
]

and zero otherwise.

Temporal envelope:

[
s(t)=
\begin{cases}
S/2,&0\le t\le t_1[4pt]
\dfrac{S}{4}
\left[
1+\cos\left(
\pi\frac{t-t_1}{t_2-t_1}
\right)
\right],
&t_1\le t\le t_2[8pt]
0,&t>t_2.
\end{cases}
]

Ranges tested:

[
S=0.01-0.5
]

[
t_1=10-25,{\rm s}
]

[
t_2=15-50,{\rm s}.
]



---

## B. Two competing transient localized stimuli

Use the same stimulus on:

* 0–10% of the domain;
* 90–100% of the domain,

with equal magnitude, reflected about the center.

This is the protocol used for their principal **cell-size comparison**. 

---

## C. Persistent gradient

For (L=10,\mu m):

[
k_S^{\rm grad}
==============

S(10-x)
]

for (0\le x\le10).

They use

[
S=0,0.01,\ldots,0.09.
]



---

# 7. Gradient reversal

First apply:

[
k_S^{\rm grad}=s(t)(10-x)
]

then after reversal:

[
k_S^{\rm grad}=s(t)x.
]

Initial-direction temporal envelope:

[
s(t)=
\begin{cases}
S,&0\le t\le t_1\
S\left(1-\frac{t-t_1}{t_2-t_1}\right),
&t_1\le t\le t_2\
0,&{\rm otherwise}.
\end{cases}
]

For the reversed gradient, equivalent timing is shifted relative to (t_{\rm reversal}).

Ranges:

[
S=0.01-0.05
]

[
t_1=t_{\rm reversal}+10\text{–}20,{\rm s}
]

[
t_2=t_{\rm reversal}+15\text{–}25,{\rm s}.
]



For WP1 this protocol is especially important because Jilkine tests reversal at the reference size but **does not systematically combine reversal with their size sweep**.

---

# 8. Noise protocol

Initial conditions are perturbed around the homogeneous steady state:

[
u(x,0)
======

u_{\rm ss}-A(1-R)
]

where

[
R\sim U(0,1)
]

and

[
A=0.01-0.1.
]

They describe this as nonuniform initial conditions around the homogeneous steady state.

The printed equation is slightly awkward relative to their statement of “zero mean”; I would inspect/reproduce the original implementation rather than trust this literally without verification.

---

# 9. Jilkine cell-size experiment

This is crucial.

They solve on a **unit computational domain**, with space

[
\xi=x/L.
]

Then:

[
D_u^{\rm effective}=\frac{D_u}{L^2},
\qquad
D_v^{\rm effective}=\frac{D_v}{L^2}.
]

They evaluate

[
L=1,3,5,\ldots,39;\mu{\rm m}.
]

They apply the **two-opposing-transient-stimulus protocol** and inspect states at

[
T=2000,{\rm s}
]

(~33 min). 

This is exactly the experiment WP1 should first reproduce rather than redefine.

Their reported outcomes:

* all three pattern-forming models fail to polarize at very small sizes;
* larger GOR and OT domains resolve competing peaks increasingly slowly;
* apparent OT two-peak states at (T=2000) are described as **quasi-stationary**;
* WP shows a transition from one peak to two persistent peaks beyond a critical size;
* that WP critical size depends on total material;
* LEGI essentially mirrors the imposed signal and shows no analogous domain-size dependence.

One slight tension in their wording is worth preserving: they describe WP's large-(L) states as “two persistent peaks,” but their summary then says all models ultimately produce a single front. For our reproduction, **long-time convergence should be checked explicitly rather than assuming either interpretation**. 

---

# 10. Parameter-robustness perturbations they already tested

Jilkine also perturbs selected parameters by **±1.5-fold**.

### GOR

Vary:

* (a)
* (b)
* initial (u)
* (c)

Strongest sensitivity: (c).

### OT

Vary:

* (a_1)
* (a_2)
* (s)
* initial (u)

Strong sensitivity to (a_2) and (s).

### LEGI

Vary:

* (k_I)
* (R_T)
* (k_R)
* (k_A)

Generally relatively insensitive.

### WP

Vary:

* initial total (u)/protein abundance
* (k_0)
* (d)
* (K)

WP is especially sensitive to (d) and (K); a 1.5-fold change could eliminate polarization altogether. 

For our purposes, this provides a useful **minimal robustness panel** to reproduce before constructing broader sweeps.

---

# What is still missing for exact implementation

I would not yet hand only the PDF formulas to a coding agent and claim they are sufficient for bit-for-bit reproduction.

We still need to resolve:

1. **WP total mass / exact homogeneous (u_0,v_0).**
2. **GOR exact constant (E_c)** used in Jilkine's comparison.
3. **GOR exact initial (u,v)** and total Cdc42.
4. Exact baseline stimulus parameters used for each panel of Fig. 2, because Methods gives **ranges**, whereas the plotted figure corresponds to particular choices.
5. Exact (t_{\rm reversal}) used in Fig. 2d.
6. Numerical discretization / solver / grid resolution — not specified here.
7. Possibly the intended implementation of the random-initial-condition equation.

Those are straightforward to recover next from the **four cited source papers and/or any archived model code**.

### Recommended WP1 implementation object

I would structure the common framework around:

```text
model:
    rhs(state, x, params, stimulus)
    default_parameters
    homogeneous_initial_state
    conserved_quantities
    observable

protocol:
    single_local
    dual_local
    persistent_gradient
    gradient_reversal
    noise
    size_sweep

analysis:
    polarization_amplitude
    number_of_peaks
    winning_side
    resolution_time
    reversal_time
    long_time_state
```

Then each of WP / GOR / OT / LEGI becomes a plugin using the exact Jilkine equations.

The **next useful step is therefore not more conceptual work**: it is to chase down those seven missing implementation details from Mori 2008, Otsuji 2007, Goryachev 2008, and Levchenko/Iglesias 2002 so that WP1 has a fully specified reproduction manifest.
