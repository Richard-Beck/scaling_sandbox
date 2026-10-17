ao et al. (2020) mechanochemical motility model — reproduction specification

Source: K. Tao et al., Tuning Cell Motility via Cell Tension with a Mechanochemical Cell Migration Model, Biophysical Journal 118 (2020), main article + supplied Supporting Information.

Reproducibility status

The papers specify the governing equations, almost all biological/mechanical parameters, the initial cell geometry and chemical state, periodic boundary conditions, and a semi-implicit Fourier-spectral numerical method. They do not specify enough to make an exact numerical reproduction unique. Missing or ambiguous items are: phase-field width ε; grid resolution; time step / exact temporal update; stimulus amplitude k_samp used for the migration figures; how the Gaussian “local noise” is spatially sampled/correlated and whether/when it is redrawn; the exact coordinate convention/update used in Eq. 11; and the operational procedure used to obtain f_max and hence f_c. Table S1 also lists σ=0.005–0.015 pN/μm for Fig. 4, whereas Fig. 4 explicitly includes σ=0.020 pN/μm.

Accordingly, the specification below is sufficient to reconstruct the reported model, but reproducing the published trajectories to numerical precision requires choosing/inferencing the omitted numerical/stimulus details.

1. State, geometry and initialization

Computational domain: a 40 × 40 μm periodic 2-D square, Ω.

Initial cell: disk Ω₀ of radius R = 5 μm at (x₀,y₀) (the paper does not state the coordinates; centering it in Ω is the natural choice). Prescribed area:

[
A_0=78.3110\ \mu\mathrm m^2.
]

Species:

u: membrane Rac-GTP (active Rac)

v: cytosolic Rac-GDP (inactive Rac)

f: polymerized F-actin

m: myosin

φ: phase field (φ≈1 inside, φ≈0 outside)

Initial phase field:

[
\phi(x,y)=\frac12\left[\tanh!\left(\frac{R-\sqrt{(x-x_0)^2+(y-y_0)^2}}{\varepsilon}\right)+1\right].
]

Interface localization function:

[
B(\phi)=3\phi^2(1-\phi)^2,
\qquad
B'(\phi)=6\phi(1-\phi)(1-2\phi).
]

Homogeneous chemical initial conditions:

[
u=2,\qquad v=6,\qquad f=0,\qquad m=2.
]

The Methods print μm² for the u,v,m initial values, but the parameter dimensions (K₁,K₃,K₄ in μm⁻²) imply that these are areal concentrations in μm⁻²; this appears to be a sign/typesetting error in the article.

2. Reaction terms and global tension

Define Rac activation rate

[
H(u,f;k_s)=b+\frac{c_1u^2}{u^2+K_1^2}+\frac{c_2f^2}{f^2+K_2^2}+k_s.
]

The diffuse-interface approximation to perimeter is

[
L_\phi=\int_\Omega\left[\frac{\varepsilon}{2}|\nabla\phi|^2+\frac{B(\phi)}{\varepsilon}\right]dx,dy.
]

The global cell tension used to inhibit F-actin is

[
mt(f,\phi)=\sigma\left(1+\int_{\Omega_0}f,dx,dy\right)L_\phi.
]

(The paper leaves implicit how Ω₀ is numerically extracted from the diffuse field when evaluating the F-actin integral.)

3. Governing diffuse-domain equations actually solved

Rac-GTP

[
\frac{\partial[B(\phi)u]}{\partial t}
=D_u\nabla!\cdot[B(\phi)\nabla u]
+B(\phi)H(u,f;k_s)v-B(\phi)ru.
\tag{6}
]

Rac-GDP

[
\frac{\partial(\phi v)}{\partial t}
=D_v\nabla!\cdot(\phi\nabla v)
-|\nabla\phi|H(u,f;k_s)v+|\nabla\phi|ru.
\tag{7}
]

F-actin

[
\frac{\partial(\phi f)}{\partial t}
=D_f\nabla!\cdot(\phi\nabla f)
+|\nabla\phi|\left[
\frac{c_3u^2}{u^2+K_3^2}
\frac{K_F}{mt(f,\phi)+K_F}
\right]
+\phi\frac{c_4K_4^2}{m^2+K_4^2}
-\phi d_ff.
\tag{8}
]

Myosin

[
\frac{\partial(\phi m)}{\partial t}
=D_m\nabla!\cdot(\phi\nabla m)
+\phi\frac{c_5K_5^2}{f^2+K_5^2}m
-\phi d_mm.
\tag{9}
]

Cell shape / motion

The membrane is overdamped. Protrusion occurs only where f > f_c, with

[
f_c=0.8f_{\max},
]

where the paper defines f_max as the maximum steady-state F-actin concentration but does not give a more precise calibration protocol.

The phase field evolves as

\left[\alpha I(f>f_c)f-\beta m\right]|\nabla\phi|
+\sigma\left(1+\int_{\Omega_0}f,dx,dy\right)
\left(\nabla^2\phi-\frac{B'(\phi)}{\varepsilon^2}\right)
-M\left(\int_\Omega\phi,dx,dy-A_0\right)|\nabla\phi|.
\tag{10}
]

The three terms are respectively protrusion/retraction, tension/curvature, and area conservation. This is obtained from overdamped balance of tension, area-pressure, protrusion/myosin and friction forces.

4. Stimuli

Stimuli modify only Rac-GDP → Rac-GTP conversion, i.e. enter as +k_s in H above.

Persistent external gradient

[
k_s^{\mathrm{grad}}=k_{s\mathrm{amp}},[D(t)\pm x],\qquad 0<x<D(t).
\tag{11}
]

The gradient is horizontal and is described as linearly decreasing from a maximum on one side of the cell to zero on the other; ± selects direction. D(t) tracks cell extent during deformation/migration. The text calls it the “maximal radius” but also states D(0)=2R (the initial diameter), so the exact coordinate convention is ambiguous.

Persistent random local stimulus

[
k_s^{\mathrm{loc}}=k_{s\mathrm{amp}}R_0,
\qquad R_0\sim N(0,1).
\tag{12}
]

The paper calls this local random internal noise and states that it is persistent. It does not specify the spatial sampling/correlation scale or whether the random values are ever refreshed. External and local stimuli may also coexist.

5. Parameters

Parameter

Value

Unit

D_u

0.1

μm² s⁻¹

D_v

50

μm² s⁻¹

b

0.001

s⁻¹

c₁

1

s⁻¹

K₁

6

μm⁻²

c₂

0.2

s⁻¹

K₂

3

μm⁻¹

r

1

s⁻¹

c₃

6

μm⁻² s⁻¹

K₃

3

μm⁻²

K_F

1

pN μm⁻¹

c₄

1

μm⁻² s⁻¹

K₄

3

μm⁻²

d_f

1

s⁻¹

D_f

0.8

μm² s⁻¹

D_m

0.3

μm² s⁻¹

c₅

1

s⁻¹

K₅

1

μm⁻¹

d_m

1

s⁻¹

β

1

pN μm⁻¹

τ

2.62

pN s μm⁻²

M

10

pN μm⁻³

A₀

78.3110

μm²

Two main mechanical parameter regimes are used:

Figures

σ (pN μm⁻¹)

α (pN μm⁻¹)

Fig. 1b–d, Fig. 3a–c, Fig. 5b

0.2–0.8

7

Fig. 2a,c, Fig. 4a–d, Fig. 5a

0.005–0.015*

0.2

* Fig. 4 itself explicitly evaluates σ=0.020 pN/μm, despite Table S1 listing only 0.005–0.015 for this figure group.

The authors note that effective tension after the F-actin contribution is roughly 10–25 pN/μm for σ=0.005–0.015 and 70–140 pN/μm for σ=0.2–0.8.

6. Numerical prescription and measured outputs

Boundary condition: periodic on the 40 × 40 μm computational domain.

Solver description: semi-implicit Fourier-spectral method for Eqs. 6–10.

Unreported: ε, grid size, Δt, exact semi-implicit split, convergence tolerances/dealiasing.

Reported observables:

instantaneous migration speed: based on successive cell-centroid positions;

Fig. 4 average speed: total distance traveled / total simulation time;

migration classes: stationary <0.1 μm/s; low 0.1–0.3 μm/s; high >0.3 μm/s;

area ratio: A/(area of equal-perimeter circle) = 4πA/P²;

aspect ratio: maximum radius / minimum radius;

persistent time: duration for which an elongated cell retains persistent movement (no numerical aspect-ratio/direction threshold is specified).

Suggested reproduction-validation suite

Because several numerical/stimulus details are omitted, these tests should be treated as model-faithfulness tests, not bitwise reproduction tests. Choose the missing numerical settings once, demonstrate grid/time-step convergence, and do not retune them between conditions. Likewise choose one k_samp per stimulus protocol and hold it fixed across the compared tensions.

V1. Deterministic persistent-gradient migration (highest-value trajectory test)

Setup: external gradient only; parameter regime α=0.2, β=1; compare σ=0.005 and σ=0.020 pN/μm; identical k_samp; run at least 160 s.

Published targets:

average speed falls from about 0.080 μm/s at σ=0.005 to 0.044 μm/s at σ=0.020;

Fig. 4 snapshots label displacement/position approximately:

σ=0.005: 0.73, 4.99, 9.83, 15.21 μm at 10, 60, 110, 160 s;

σ=0.020: 0.39, 2.77, 5.46, 8.23 μm at the same times;

low tension produces a visibly elongated cell; high tension remains much rounder.

Suggested pass: both reported average speeds within ±10%, and the low/high speed ratio within ±15% of 0.080/0.044 ≈ 1.82. If using the same position origin as the figure, require the 160-s values within ±15%. Failure of the direction or monotonic tension effect is a hard fail.

V2. Random-noise deformation statistics

Setup: persistent local random stimulus only; α=0.2; compare σ=0.005 and 0.015 pN/μm; 60-min simulations. Use at least 20 independent realizations (paper reports five for the quoted deformation statistics).

Published targets:

maximum perimeter: 35.5 ± 1.4 μm at low tension vs 31.5 ± 0.1 μm at high tension;

area ratio: 0.79 ± 0.06 vs 0.99 ± 0.01;

low tension develops multiple F-actin clusters/protrusions; high tension stays close to circular.

Suggested pass: ensemble means within 10% of all four published means, with the correct ordering (P_low>P_high, area-ratio_low < area-ratio_high). The high-tension area ratio should be ≥0.95.

V3. Random-noise migration-mode distribution

Setup: persistent local random stimulus only; α=7; compare σ=0.2 and 0.8 pN/μm; classify runs by mean instantaneous speed using the authors’ thresholds. Use ≥50 realizations for a stable reproduction check.

Published targets (10-run samples):

σ=0.2: 8/10 high-migration, 2/10 low-migration, essentially no stationary runs;

σ=0.8: 7/10 stationary, about 3/10 migrating.

Suggested pass: at low tension, ≥70% high-migration and ≤10% stationary; at high tension, ≥60% stationary. A reversal of these distributions is a hard fail.

V4. Tension suppression of polarized F-actin

Setup: persistent external gradient; α=7; identical stimulus at σ=0.2 and 0.8 pN/μm; evaluate the steady polarized state.

Published target: maximum F-actin falls from about 0.7 μm⁻¹ at σ=0.2 to about 0.3 μm⁻¹ at σ=0.8.

Suggested pass: maxima within ±15% of those values and a high/low ratio of 0.35–0.55. The F-actin cap must occur on the stimulated side in both conditions.

Critical implementation decisions that must be recorded in any reproduction

For a genuinely auditable reproduction, explicitly report the chosen values/rules for: ε; grid; Δt; the semi-implicit splitting; k_samp; the spatial representation and persistence/redrawing rule of R₀; the Eq. 11 coordinate convention and how D(t) is computed; the numerical definition of Ω₀ inside the diffuse field; the procedure for determining f_max/f_c; and centroid/perimeter/radius extraction from φ. These are not fully specified in the supplied article or Supporting Information.