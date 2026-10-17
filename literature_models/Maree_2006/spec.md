Treat the following as the **reconstruction specification for the original Marée et al. 2006 model**, using the 2012 paper only where it explicitly clarifies numerical mechanics that the 2006 paper leaves vague. The biochemical/actin model itself is unusually well specified; the remaining uncertainty is mostly in low-level numerical implementation and moving-boundary bookkeeping.

I use three tags below: **[06]** = explicitly specified in the 2006 paper; **[12]** = explicit in the 2012 continuation and probably reflects the same implementation lineage, but is not guaranteed identical to 2006; **[INF]** = reconstruction choice/inference where neither paper nails the implementation down.

## 1. Model state and domain

Use a 2-D **hexagonal lattice**, \(500\times500\), with toroidal outer boundaries. The characteristic lattice spacing is

$$
\Delta x=100\ {\rm nm}.
$$

One CPM Monte Carlo step corresponds to

$$
\Delta t_{\rm MCS}=0.1\ {\rm s}.
$$

The cell occupies approximately \(3\times10^4\) lattice sites, with target projected area

$$
A_0=3\times10^8\ {\rm nm^2},
$$

and is initially circular with diameter approximately \(20\,\mu{\rm m}\). ([ResearchGate][1])

At each intracellular site track:

| Variables            | Meaning                              |               Units |
| -------------------- | ------------------------------------ | ------------------: |
| \(C,R,\rho\)         | active Cdc42, Rac, Rho               |      \(\mu{\rm M}\) |
| \(C_i,R_i,\rho_i\)   | inactive forms                       |      \(\mu{\rm M}\) |
| \(A\)                | activated Arp2/3                     |      \(\mu{\rm M}\) |
| \(F_m,\;m=0\ldots5\) | F-actin density in six orientations  |         nm\(^{-1}\) |
| \(B_m,\;m=0\ldots5\) | freely extending barbed ends         |         nm\(^{-2}\) |
| \(P_m,\;m=0\ldots5\) | membrane-engaged pushing barbed ends | nm\(^{-1}\) of edge |

The six actin orientations are

$$
\theta_m=\frac{2\pi m}{6},
\qquad
\mathbf V_m=v_0(\cos\theta_m,\sin\theta_m),
$$

so the true \(\sim70^\circ\) Arp2/3 branching angle is approximated as \(60^\circ\). ([ResearchGate][1])

Thus the 2006 system contains **7 biochemical fields + 12 bulk actin fields + 6 edge actin variables**, plus the CPM cell mask. The 2012 paper explicitly describes the actin part as 18 coupled orientation variables \(F,B,P\). ([PLOS][2])

---

# 2. Rho-family signalling system

Define the local reaction terms

$$
f_C=
\frac{I_C}{1+(\rho/\beta_\rho)^n}
\frac{C_i}{C_{\rm tot}}
-d_C C,
$$

$$
f_R=
(I_R+\alpha_C C)
\frac{R_i}{R_{\rm tot}}
-d_R R,
$$

$$
f_\rho=
\frac{I_\rho+\alpha_R R}
{1+(C/\beta_C)^n}
\frac{\rho_i}{\rho_{\rm tot}}
-d_\rho\rho .
$$

Then

$$
\frac{\partial C}{\partial t}
=f_C+D_m\nabla^2C,
$$

$$
\frac{\partial R}{\partial t}
=f_R+D_m\nabla^2R,
$$

$$
\frac{\partial \rho}{\partial t}
=f_\rho+D_m\nabla^2\rho,
$$

and

$$
\frac{\partial C_i}{\partial t}
=-f_C+D_{mc}\nabla^2C_i,
$$

$$
\frac{\partial R_i}{\partial t}
=-f_R+D_{mc}\nabla^2R_i,
$$

$$
\frac{\partial \rho_i}{\partial t}
=-f_\rho+D_{mc}\nabla^2\rho_i.
$$

These are essentially Eqs. 5–8 of the original model. ([ResearchGate][1])

All six species are confined to the cell and obey **no-flux/Neumann conditions at the moving cell boundary**. The active membrane-associated forms diffuse much more slowly than the inactive membrane/cytosolic pool. ([ResearchGate][1])

The parameter set is:

| Parameter          |                        Value |
| ------------------ | ---------------------------: |
| \(C_{\rm tot}\)    |          \(2.4\,\mu{\rm M}\) |
| \(R_{\rm tot}\)    |          \(7.5\,\mu{\rm M}\) |
| \(\rho_{\rm tot}\) |          \(3.1\,\mu{\rm M}\) |
| \(C^\*\)           |            \(1\,\mu{\rm M}\) |
| \(R^\*\)           |            \(3\,\mu{\rm M}\) |
| \(\rho^\*\)        |         \(1.25\,\mu{\rm M}\) |
| \(I_C\)            |  \(3.4\,\mu{\rm M\,s^{-1}}\) |
| \(I_R\)            |  \(0.5\,\mu{\rm M\,s^{-1}}\) |
| \(I_\rho\)         |  \(3.3\,\mu{\rm M\,s^{-1}}\) |
| \(\beta_\rho\)     |         \(1.25\,\mu{\rm M}\) |
| \(\beta_C\)        |            \(1\,\mu{\rm M}\) |
| \(n\)              |                            3 |
| \(\alpha_C\)       |        \(4.5\,{\rm s^{-1}}\) |
| \(\alpha_R\)       |        \(0.3\,{\rm s^{-1}}\) |
| \(d_C=d_R=d_\rho\) |          \(1\,{\rm s^{-1}}\) |
| \(D_m\)            | \(10^5\,{\rm nm^2\,s^{-1}}\) |
| \(D_{mc}\)         | \(10^7\,{\rm nm^2\,s^{-1}}\) |

These are directly from Table 1. ([ResearchGate][1])

### Initialization

**[INF, but strongly implied]**

Uniformly initialize

$$
C=C^\*=1,\qquad
R=R^\*=3,\qquad
\rho=\rho^\*=1.25
$$

and therefore

$$
C_i=C_{\rm tot}-C^\*=1.4,
$$

$$
R_i=R_{\rm tot}-R^\*=4.5,
$$

$$
\rho_i=\rho_{\rm tot}-\rho^\*=1.85
$$

in \(\mu{\rm M}\).

The paper says that all variables are initialized at their “typical levels” and that the initial cell is unpolarized, but does not print this explicit initialization vector. ([ResearchGate][1])

---

# 3. Arp2/3

Let

$$
F_t=\sum_{m=0}^{5}F_m.
$$

The branching/nucleation function is

$$
\eta(A,F_t)
=
\eta_0
\frac{A}
{K_m+A+\ell F_t}.
$$

Arp2/3 evolves as

$$
\frac{\partial A}{\partial t}
=
\mu_C C+\mu_RR
-\eta(A,F_t)F_t
-d_AA
+D_A\nabla^2A.
$$

Again use no-flux boundary conditions. ([ResearchGate][1])

Parameters:

$$
A^\*=2\ \mu{\rm M},
$$

$$
\mu_C=\mu_R=0.16\ {\rm s^{-1}},
$$

$$
d_A=0.1\ {\rm s^{-1}},
$$

$$
D_A=10^6\ {\rm nm^2\,s^{-1}},
$$

$$
\eta_0=60\ \mu{\rm M\,nm\,s^{-1}},
$$

$$
K_m=2\ \mu{\rm M},
$$

$$
\ell=255\ \mu{\rm M\,nm}.
$$

These come from Table 2 and its Appendix derivation. ([ResearchGate][1])

---

# 4. Oriented actin network

For every orientation \(m\),

$$
\frac{\partial F_m}{\partial t}
=
v_0B_m-d_FF_m,
$$

and

$$
\frac{\partial B_m}{\partial t}
=
-\nabla\cdot(B_m\mathbf V_m)
+
\frac{1}{2}k\eta(A,F_t)
\left(F_{m-1}+F_{m+1}\right)
-\kappa(R)B_m,
$$

with orientation indices taken modulo six. ([ResearchGate][1])

The Rac-dependent capping law is

$$
\kappa(R)
=
\kappa_{\max}
-\kappa_{\rm Rac}
\frac{R}{K_R+R}.
$$

Equivalently,

$$
\kappa(R)
=
\kappa_{\rm basal}
+
\frac{\kappa_{\rm Rac}}
{1+R/K_R},
\qquad
\kappa_{\rm basal}=0.7\,{\rm s^{-1}}.
$$

([ResearchGate][1])

Parameters:

| Parameter            |                                            Value |
| -------------------- | -----------------------------------------------: |
| \(v_0\)              |                        \(500\,{\rm nm\,s^{-1}}\) |
| \(d_F\)              |                           \(0.03\,{\rm s^{-1}}\) |
| \(k\)                | \(1.06\times10^{-4}\ {\rm nm^{-2}\,\mu M^{-1}}\) |
| \(\kappa_{\max}\)    |                            \(2.8\,{\rm s^{-1}}\) |
| \(\kappa_{\rm Rac}\) |                            \(2.1\,{\rm s^{-1}}\) |
| \(K_R\)              |                                \(3\,\mu{\rm M}\) |
| \(r\)                |                                             0.14 |
| \(F^\*\)             |                         \(0.278\,{\rm nm^{-1}}\) |
| \(B^\*\)             |              \(1.7\times10^{-5}\,{\rm nm^{-2}}\) |

([ResearchGate][1])

For an isotropic initial condition I would use

$$
F_m=\frac{F^\*}{6},
\qquad
B_m=\frac{B^\*}{6}.
$$

That division by six is **[INF]**: it is the obvious isotropic interpretation, but I have not found a sentence in the 2006 paper stating it.

---

# 5. Pushing barbed ends at the membrane

At the edge, define \(P_m\) as pushing barbed ends per unit membrane length. The paper gives

$$
\frac{\partial P_m}{\partial t}
=
B_m(\mathbf V_m\cdot\hat{\mathbf n})
-r\kappa(R)P_m,
$$

where \(\hat{\mathbf n}\) is the outward membrane normal. ([ResearchGate][1])

Physically, this equation is conversion of freely growing \(B_m\) into stalled/load-bearing \(P_m\) when the filament reaches the membrane.

### Important implementation ambiguity

The printed equation contains \(\mathbf V_m\cdot\hat n\), not an explicitly rectified

$$
[\mathbf V_m\cdot\hat n]_+.
$$

Negative values would imply negative production of \(P_m\), which makes no physical sense. The surrounding text says these are barbed ends **reaching and pushing against the edge**, strongly suggesting that only outward-facing classes participate. ([ResearchGate][1])

I would therefore implement

$$
\frac{dP_m}{dt}
=
B_m[\mathbf V_m\cdot\hat n]_+
-r\kappa(R)P_m
$$

and record this as an **[INF] implementation choice**.

The typical reported total pushing-end density is

$$
P^\*=0.05\ {\rm nm^{-1}}.
$$

([ResearchGate][1])

I would **not** simply initialize every \(P_m=P^\*/6\) around a circular membrane. The natural equilibrium depends upon the local normal. Better is to initialize \(P_m=0\) and allow the isotropic actin network to populate the membrane before stimulation, or analytically initialize the local \(P_m\) from the above ODE. This initialization detail is not specified in 2006.

---

# 6. CPM mechanics

For the one-cell/medium system, the Hamiltonian is effectively

$$
H
=
J_{CM}N_{CM}
+
\lambda(a-A_0)^2,
$$

where \(N_{CM}\) is the number of cell-medium neighbour contacts under the chosen CPM interaction neighbourhood, \(a\) is current area and \(A_0\) target area.

The paper says the interfacial term is evaluated over neighbours **“up to third order.”** ([ResearchGate][1])

Parameters:

$$
J_{CM}=7.5\times10^{-4}\ {\rm nm^{-1}},
$$

$$
\lambda=10^{-9}\ {\rm nm^{-3}},
$$

$$
A_0=3\times10^8\ {\rm nm^2}.
$$

The chosen \(\lambda\) keeps area fluctuations within a few percent. ([ResearchGate][1])

## Copy algorithm

During every MCS:

1. consider lattice sites in random order;
2. each lattice site receives one copy attempt;
3. choose a random neighbouring site;
4. propose copying its cell/medium state into the target;
5. calculate the change in Hamiltonian. ([ResearchGate][1])

The standard acceptance probability used here is

$$
P_{\rm acc}
=
\begin{cases}
1, & \Delta H<-H_b,\\[4pt]
\exp[-(\Delta H+H_b)/T], & \Delta H\ge -H_b .
\end{cases}
$$

with

$$
T=0.008\ {\rm nm^{-1}},
\qquad
H_b=0.046\ {\rm nm^{-1}}.
$$

([ResearchGate][1])

Notice that this is **not** the most common textbook CPM expression \(\min[1,e^{-\Delta H/T}]\); the membrane-yield offset \(H_b\) matters.

---

# 7. Actin/Rho modification of CPM copy energies

For an attempted **extension**,

$$
\Delta H'
=
\Delta H
-\sum_m P_m
+\xi(\rho-\rho_{\rm th}),
$$

whereas for a **retraction**,

$$
\Delta H'
=
\Delta H
+\sum_mP_m
-\xi(\rho-\rho_{\rm th}).
$$

Use \(\Delta H'\), rather than \(\Delta H\), in the Metropolis probability. ([ResearchGate][1])

The logic is:

* pushing barbed ends favour extension;
* pushing ends resist retraction;
* elevated Rho favours retraction.

The Rho threshold is

$$
\rho_{\rm th}=1.25\ \mu{\rm M},
$$

and

$$
\xi=0.06
$$

for the standard simulation. ([ResearchGate][1])

### Two important ambiguities

First, the prose explicitly says the Rho contractile contribution acts **where \(\rho>\rho_{\rm th}\)**. Thus I would implement

$$
\xi[\rho-\rho_{\rm th}]_+
$$

rather than allowing low Rho to generate an opposite “negative contraction.” ([ResearchGate][1])

Second, \(\sum_mP_m\) is described as the ends pushing **towards the empty site** involved in the proposed extension. Precisely how membrane segments/orientations are selected for a particular copy attempt is not specified. ([ResearchGate][1])

The clean reconstruction is:

$$
P_{\rm push}(x,\hat d)
=
\sum_{m:\mathbf V_m\cdot\hat d>0}P_m(x),
$$

where \(\hat d\) is the direction of the attempted boundary displacement. I would initially **not cosine-weight this term**, because Eq. 24 simply sums \(P_m\), but I would make this geometric selection switchable because it is one of the uncertainties that could materially affect results.

### \(\xi\) units

The machine-readable 2006 copy renders the unit as something resembling \(\mu{\rm M}\,{\rm nm^{-1}}\). Dimensionally Eq. 24 requires

$$
[\xi]
=
\mu{\rm M}^{-1}{\rm nm^{-1}},
$$

because \(\xi(\rho-\rho_{\rm th})\) must have the same nm\(^{-1}\) units as \(P\). This is exactly the sort of typography issue for which I would like the native PDF before declaring the specification completely settled.

---

# 8. Moving boundary and actin bookkeeping

This is one of the most important implementation pieces.

The **2006 paper explicitly says** that if a membrane site retracts, its filaments and barbed ends are pushed backward with the membrane and accumulate at adjacent intracellular sites according to their orientations; the barbed ends then encounter the new membrane. ([ResearchGate][1])

The **2012 paper clarifies this considerably**:

* retraction transfers \(F_m/B_m\) inward while preserving their orientations;
* some free \(B_m\) thereby become pushing \(P_m\);
* upon extension, pushing \(P_m\) lose membrane contact and are demoted to free \(B_m\);
* no filaments or barbed ends are created or destroyed merely because the CPM boundary moved. ([PLOS][2])

A sensible explicit implementation is therefore:

### Accepted retraction \(x:\mathrm{cell}\rightarrow\mathrm{medium}\)

For each orientation \(m\):

$$
F_m(x)\rightarrow F_m(x-\mathbf e_m),
$$

$$
B_m(x)\rightarrow B_m(x-\mathbf e_m),
$$

using the closest intracellular site opposite the filament growth direction. If that receiving site is now at the membrane and the filament direction intersects the membrane, convert the appropriate number of \(B_m\) into \(P_m\).

### Accepted extension \(y:\mathrm{medium}\rightarrow\mathrm{cell}\)

Any \(P_m\) on the old membrane segment that has become internal should become free \(B_m\), with the correct area/length units conversion, and thereafter advect outward until it reaches the new edge again.

The 2012 paper explicitly says a correction factor is used for the \(B\,[{\rm nm^{-2}}]\leftrightarrow P\,[{\rm nm^{-1}}]\) conversion, but the HTML equations are image-rendered and I have not recovered the precise factor from accessible text. ([PLOS][2])

**That is a genuine unresolved item.** The obvious dimensional conversion for a lattice layer is proportional to \(\Delta x\), i.e.

$$
P\sim B\,\Delta x,
$$

and inversely \(B\sim P/\Delta x\), but I would not call that the historical implementation without seeing the original equation/code.

---

# 9. Moving-domain bookkeeping for soluble fields

This is barely specified in 2006 but explicit in 2012.

The later implementation does the following after CPM deformation:

* when the cell extends into an extracellular pixel, the new site inherits the concentrations of the intracellular source pixel, giving a locally flat profile;
* concentrations are then globally renormalized so membrane movement itself neither creates nor destroys molecular mass;
* the correction was kept below \(0.1\%\) per timestep;
* intracellular hydrodynamic convection is **not** modelled. ([PLOS][2])

I would reproduce that.

For each scalar field \(q\), after the complete MCS, enforce

$$
\sum_{x\in\Omega_{\rm new}}q_{\rm corrected}(x)A_{\rm pix}
=
\sum_{x\in\Omega_{\rm old}}q_{\rm preCPM}(x)A_{\rm pix}.
$$

**[12 rather than 06]**.

Whether their exact implementation renormalized each field individually or conserved combined active+inactive GTPase pools is not perfectly clear from the prose. I would independently conserve each field's value across the purely geometric CPM operation; reaction terms have already handled biochemical conversion.

---

# 10. Numerical integration

This is the largest unresolved historical detail.

The 2006 paper states all three of the following:

* one MCS \(=0.1\) s;
* “the same timestep” is used for the PDE integration;
* diffusion was integrated with explicit Euler using sufficiently small timesteps to guarantee stability. ([ResearchGate][1])

Those statements cannot literally mean **one 0.1-s explicit diffusion update**. With

$$
D_{mc}=10^7\ {\rm nm^2s^{-1}},
\qquad
\Delta x=100\ {\rm nm},
$$

the characteristic explicit-diffusion stability scale is

$$
\frac{\Delta x^2}{D}
=
10^{-3}\ {\rm s},
$$

and the actual stability limit is smaller by the stencil-dependent factor, i.e. typically several \(10^{-4}\) s. A \(0.1\)-s Euler step would be hundreds of times too large.

So I would use:

**[INF]**

$$
\Delta t_{\rm PDE}\lesssim2\times10^{-4}\ {\rm s}
$$

for an explicit implementation and perform enough PDE substeps to span each 0.1-s MCS. Calculate the precise CFL bound from your chosen hexagonal Laplacian.

Alternatively, if the goal is biological rather than numerical archaeology, the 2012 authors switched explicitly to a **hexagonal ADI implicit solver**, advancing successively along the three lattice axes for one-third timesteps. ([PLOS][2])

For reproducing **2006**, however, I would use explicit substepping first.

### Spatial differencing

The 2006 paper does **not** specify:

* exact discrete hexagonal Laplacian weights;
* the scheme for \(-\nabla\cdot(B_m\mathbf V_m)\);
* whether the barbed-end transport is donor-cell/upwind, centered, etc.

Because \(\mathbf V_m\) points exactly along a hexagonal lattice direction, first-order conservative upwind transport is an extremely natural implementation:

$$
B_m^{n+1}(x)
=
B_m^n(x)
-
\frac{v_0\Delta t}{\Delta x}
\left[B_m(x)-B_m(x-\mathbf e_m)\right]+\cdots.
$$

But that is **[INF]**, not documented historical fact.

---

# 11. Recommended macro-update order

The most defensible reconstruction is:

$$
\boxed{\text{reaction/diffusion/actin evolution}
\rightarrow
\text{edge }P_m\text{ evolution}
\rightarrow
\text{one CPM MCS}
\rightarrow
\text{moving-boundary bookkeeping}}
$$

for every \(0.1\) s macrostep.

The 2012 paper explicitly states that each timestep is **reaction-diffusion first, then one CPM timestep**. ([PLOS][2])

Within the CPM MCS, I would apply the \(F/B/P\) promotion/demotion immediately after each accepted boundary copy, since subsequent copy attempts can involve the newly altered membrane. That fine-grained ordering is not stated explicitly.

---

# 12. Polarization protocol

The original 2006 model does **not spontaneously polarize** from the stated uniform starting state. Instead a spatial gradient in \(I_C\), the basal Cdc42 activation rate, is imposed for

$$
10\ {\rm s}=100\ {\rm MCS},
$$

then completely removed. ([ResearchGate][1])

The cell should then maintain polarity autonomously.

I have **not found the magnitude of the initial 2006 polarization gradient** in the accessible 2006 text. This appears genuinely underspecified.

The 2012 implementation uses a **roughly 15% variation in \(I_C\) across the cell for 10 s**, and describes that as its “standard polarization” protocol. ([PubMed][3])

Therefore, for a working 2006 reproduction I would use

$$
I_C(x)
=
I_C^0[1+g\,s(x)]
$$

with \(s\) varying from \(-1/2\) to \(+1/2\) across the cell and

$$
g=0.15
$$

for 10 s, then set \(g=0\).

Label that **[12-derived]**, not exact 2006.

---

# 13. Validation ladder

This model gives us unusually good internal checks. I would validate it **bottom-up**, and not even run the full migration model until the first few pass.

## A. Hard numerical/conservation checks

### 1. GTPase mass conservation

With a frozen cell and no boundary movement,

$$
M_C=\int(C+C_i)dA,
$$

$$
M_R=\int(R+R_i)dA,
$$

$$
M_\rho=\int(\rho+\rho_i)dA
$$

must remain constant to numerical precision.

Failure here means either reaction signs or no-flux diffusion are wrong.

With CPM movement enabled, these totals should remain invariant apart from your explicitly quantified moving-domain numerical correction.

### 2. Isotropic state must remain isotropic

Initialize the circular cell uniformly with

$$
F_0=\cdots=F_5,\quad
B_0=\cdots=B_5.
$$

Without stimulation, there should be no deterministic directional drift in any orientation class. Across stochastic replicates, mean COM displacement must vanish.

### 3. Resting-cell control

This is one of the paper's clearest whole-model controls.

With the standard parameter set but **no initial \(I_C\) gradient**, the model remains:

* essentially circular;
* unpolarized;
* essentially stationary.

([ResearchGate][1])

If your implementation spontaneously produces a persistent front, the signalling/noise/boundary bookkeeping is wrong.

### 4. Cell area

At steady state, area should remain around

$$
3\times10^8\,{\rm nm^2}
$$

with fluctuations only of order a few percent. The authors explicitly chose \(\lambda\) for this. ([ResearchGate][1])

---

# 14. Signalling-subsystem checks

These are particularly useful because they isolate the chemistry from CPM complications.

### 5. Well-mixed bifurcation behaviour

The authors explicitly analyze the homogeneous GTPase equations.

With \(n=3\), the well-mixed system should **not** show ordinary bistability; with \(n=6\), it does. ([ResearchGate][1])

That is an excellent unit test of the reaction equations before introducing diffusion.

### 6. 1-D spatial polarity

On a fixed 20-\(\mu\)m 1-D domain with no-flux boundaries and explicitly diffusing inactive forms, an imposed initial front/back perturbation should settle to a persistent spatial state with high Cdc42 at one end and high Rho at the other.

The paper finds the spatially polarized regime much broader than the homogeneous bistability region. For \(n=6\), they report persistent spatial bistability approximately over

$$
1.7<I_C<43.4\ \mu{\rm M\,s^{-1}},
$$

versus approximately

$$
2.5<I_C<6.5
$$

for the corresponding well-mixed bistability. ([ResearchGate][1])

The \(n=3\) standard system should also support a stable spatial front/back state even though the well-mixed system is not bistable. This is one of the core phenomena of their model.

If you cannot reproduce this with a **fixed domain**, don't proceed to the CPM.

---

# 15. Actin-subsystem checks

### 7. Basic steady-state arithmetic

From Eq. 13,

$$
v_0B^\*
\approx
d_FF^\*.
$$

Numerically,

$$
500(1.7\times10^{-5})
=
8.5\times10^{-3}\ {\rm s^{-1}}
$$

and

$$
0.03(0.278)
=
8.34\times10^{-3}\ {\rm s^{-1}},
$$

so they agree almost exactly.

The authors explicitly used this relationship to derive \(B^\*\). ([ResearchGate][1])

This is a very useful implementation sanity check.

### 8. Capping rates

At low Rac,

$$
\kappa\rightarrow2.8\ {\rm s^{-1}},
$$

while at high Rac,

$$
\kappa\rightarrow0.7\ {\rm s^{-1}}.
$$

At the membrane this is multiplied by \(r=0.14\), giving a slowest edge capping rate of approximately

$$
0.1\ {\rm s^{-1}}.
$$

The Appendix explicitly states this. ([ResearchGate][1])

---

# 16. Protrusion-mechanics check — probably the best hard CPM validation

The authors deliberately calibrated \(T\) and \(H_b\) so that an isolated straight membrane pushed by actin reproduces the thermal-ratchet force/velocity law

$$
v_{\rm TR}(b)
=
v_0e^{-w/b},
$$

where

$$
w=0.05\ {\rm nm^{-1}}.
$$

([ResearchGate][1])

For the CPM they derive the implicit relation

$$
v=
\frac{\Delta x}{\Delta t}
e^{-H_b/T}
\left[
e^{(1-v/v_0)b/T}
-
e^{-(1-v/v_0)b/T}
\right].
$$

([ResearchGate][1])

With

$$
T=0.008\,{\rm nm^{-1}},
\qquad
H_b=0.046\,{\rm nm^{-1}},
$$

the two curves should nearly coincide for biologically relevant

$$
0.05\lesssim b\lesssim0.25\ {\rm nm^{-1}}.
$$

This is an unusually powerful reproduction check because it tests your **Metropolis rule, \(P\) coupling, time conversion and space conversion simultaneously**.

I would implement this as an automated regression test.

---

# 17. Whole-cell hard/semihard checks

After the 10-s polarization pulse, the cell should:

1. establish **high Cdc42 + Rac at the front**;
2. establish **high Rho at the rear**;
3. show a sharp front/back biochemical interface;
4. create a broad smooth protrusive leading edge;
5. generate rear contraction;
6. settle into a crescent-like keratocyte shape;
7. move persistently at roughly

$$
v_{\rm COM}\approx90\ {\rm nm\,s^{-1}}
=5.4\ \mu{\rm m\,min^{-1}}.
$$

The reported speed is close to constant once established. ([ResearchGate][1])

I would not require 90.000 nm/s. Given the underspecified numerical details, something like 70–110 nm/s with the right architecture would already tell me we were in the right model regime. A factor-of-two miss would make me investigate implementation choices.

---

# 18. Spatial actin architecture checks

The original paper gives surprisingly specific qualitative spatial predictions:

* free barbed ends should peak **very close to the leading edge**;
* F-actin density should peak **a few microns behind the leading edge**;
* F-actin then declines toward the rear;
* in a broad region behind the front, most filaments and barbed ends should point within approximately **\(\pm60^\circ\)** of the direction toward the leading edge. ([ResearchGate][1])

That last point is particularly good for validating your orientation transport/boundary-selection code.

Calculate, for example,

$$
\mathbf p_F(x)
=
\frac{\sum_m F_m(x)\hat{\mathbf e}_m}
{\sum_mF_m(x)}
$$

and similarly \(\mathbf p_B\).

Near the leading edge you should see a substantial nonzero orientation vector pointing forward; the resting cell should instead remain approximately isotropic.

---

# 19. Contractility sweep — excellent phenotype test

The paper sweeps the Rho contraction coefficient \(\xi\).

It reports approximately:

|         \(\xi\) | Expected behaviour                                            |
| --------------: | ------------------------------------------------------------- |
|     \(\le0.01\) | no net translocation; polarized cell tends to rotate/wobble   |
|   \(0.02-0.04\) | slow, irregular migration with sudden turns                   |
|   \(0.05-0.07\) | regular persistent movement and typical keratocyte morphology |
| \(\gtrsim0.08\) | contraction dominates; increasingly triangular/unstable shape |

Standard value is \(0.06\). ([ResearchGate][1])

This is an extremely valuable **soft fingerprint**: reproducing not just one simulation but this sequence of phenotype transitions makes accidental agreement much less likely.

---

# 20. Chemotactic turning check

Once a cell is already moving, apply a shallow \(I_C\) gradient perpendicular to its trajectory.

The original model is remarkably sensitive: it responds to gradients as small as

$$
10^{-6}\,
\mu{\rm M\,s^{-1}\,nm^{-1}},
$$

corresponding to roughly a **1% variation across a 30-\(\mu\)m cell**. ([ResearchGate][1])

For a stronger **10% \(I_C\) gradient**, the paper's example rotates roughly

$$
90^\circ
$$

in about

$$
3\ {\rm min}.
$$

([ResearchGate][1])

This is a good late-stage validation of the chemistry–shape feedback.

---

# 21. What I would automate as the reproduction test suite

If this were my implementation, I would require these checks, in roughly this order:

| Test                            | Pass criterion                                                         |
| ------------------------------- | ---------------------------------------------------------------------- |
| Reaction-only pool conservation | \(C+C_i,\;R+R_i,\;\rho+\rho_i\) conserved                              |
| Isotropy                        | six orientation classes remain equivalent without stimulus             |
| Fixed 1-D chemistry             | persistent Cdc42-front/Rho-back spatial state                          |
| Actin equilibrium               | \(v_0B^\*\approx d_FF^\*\)                                             |
| Thermal-ratchet test            | CPM \(v(b)\) tracks Eq. 25 over \(b=0.05-0.25\) nm\(^{-1}\)            |
| Rest cell                       | round + essentially stationary                                         |
| Polarization                    | 10-s pulse produces persistent Cdc42/Rac-front, Rho-back               |
| Cell area                       | within a few percent of \(3\times10^8\) nm²                            |
| Migration speed                 | order \(90\) nm/s                                                      |
| Actin organization              | \(B\) edge peak; \(F\) peak a few \(\mu\)m behind; forward orientation |
| \(\xi\) sweep                   | stationary/irregular → keratocyte → triangular sequence                |
| Turning                         | shallow \(I_C\) gradients turn an already polarized cell               |

The **thermal-ratchet test + 1-D signalling test + \(\xi\) phenotype sweep** are probably the three most diagnostic things beyond simply looking at a movie.

---

# 22. Remaining uncertainty ledger

These are the pieces that prevent me from calling the above a pixel-for-pixel reproduction specification.

| Issue                                                          | Status                                                        |                              Importance |
| -------------------------------------------------------------- | ------------------------------------------------------------- | --------------------------------------: |
| Exact hexagonal diffusion stencil                              | not given in 2006                                             |                                **high** |
| Explicit-Euler diffusion substep                               | contradictory/underspecified wording                          |                                **high** |
| \(B_m\) advection discretization                               | not specified                                                 |                                **high** |
| CPM “third-order neighbours” exact coordinate set              | phrase given, exact stencil not enumerated                    |                             medium-high |
| Copy-proposal neighbour set                                    | “random neighbour”, exact neighbourhood not stated            |                                  medium |
| Connectivity constraint                                        | not mentioned                                                 |                                  medium |
| Initial \(F_m/B_m\) distribution                               | isotropic division by six inferred                            |                              low-medium |
| Initial \(P_m\)                                                | not specified                                                 |                                  medium |
| Exact magnitude of initial 2006 polarization gradient          | not found                                                     |                  medium; 2012 gives 15% |
| \(P\leftrightarrow B\) conversion factor during membrane moves | later paper says factor exists; equation inaccessible in HTML |                                **high** |
| Exact receiving pixel for retracted oriented actin             | qualitative rule given                                        |                                  medium |
| Which \(P_m\) contribute to an individual CPM copy             | prose but no explicit geometric algorithm                     |                                **high** |
| Positive-part handling in \(P_m\) production                   | physically implied, not printed                               |                                  medium |
| Positive-part handling of Rho contraction                      | prose strongly implies it                                     |                                  medium |
| Typography/units of \(\xi\)                                    | text extraction ambiguous                                     | low mathematically, but worth verifying |
| Soluble-field mass redistribution in 2006                      | clarified only in 2012                                        |                                  medium |

Those first, third, tenth and twelfth are the ones I would investigate most aggressively if a first reconstruction did not reproduce the paper.

---

# 23. Sources actually accessed

The primary source is **Marée AFM, Jilkine A, Dawes A, Grieneisen VA, Edelstein-Keshet L. “Polarization and Movement of Keratocytes: A Multiscale Modelling Approach.” Bulletin of Mathematical Biology 68 (2006), 1169–1211. PMID 16794915.** ([PubMed][4])

I was able to access essentially the **complete text, equations, tables, appendix and figure captions** through an author-uploaded ResearchGate rendering. That is where almost all equations and values above came from. The direct “download full-text PDF” route failed, however, so I have **not visually inspected a pristine/native PDF**. That matters particularly for superscripts, Table 3's \(\xi\) units, and equations or diagrams whose formatting may be lost in text extraction.

The main secondary source is **Marée, Grieneisen & Edelstein-Keshet 2012, PLoS Computational Biology 8:e1002402**, which is fully open-access. It explicitly states that it builds on the earlier 2-D model and gives much more detail on moving-domain concentration conservation, \(B/P\) promotion/demotion, reaction-diffusion/CPM update order, and their later ADI solver. ([PubMed][3]) Its Supporting Information is simulation videos rather than a hidden methods/code supplement. ([PLOS][2])

I also accessed the ResearchGate rendering of the **Marée, Grieneisen & Hogeweg 2007 CPM chapter**, but its extracted body text is badly corrupted, so I have **not used it to fill any implementation gaps**.

Finally, the old Utrecht `theory.bio.uu.nl/stan/keratocyte` companion site is referenced by the 2006 paper as the location of the simulation movies, but it is no longer retrievable through my browser. The paper itself describes that resource as **movies**, not source code or a methodological supplement. ([ResearchGate][1])


[1]: https://www.researchgate.net/publication/6989407_Polarization_and_Movement_of_Keratocytes_A_Multiscale_Modelling_Approach "(PDF) Polarization and Movement of Keratocytes: A Multiscale Modelling Approach"
[2]: https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002402 "How Cells Integrate Complex Stimuli: The Effect of Feedback from Phosphoinositides and Cell Shape on Cell Polarization and Motility | PLOS Computational Biology"
[3]: https://pubmed.ncbi.nlm.nih.gov/22396633/ "How cells integrate complex stimuli: the effect of feedback from phosphoinositides and cell shape on cell polarization and motility - PubMed"
[4]: https://pubmed.ncbi.nlm.nih.gov/16794915/ "pubmed.ncbi.nlm.nih.gov"
