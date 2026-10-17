# Methods for the `migration1d` experiments

## Purpose and scope

This document consolidates the methodology implemented in
`migration1d`. It describes the completed single-species and three-species
calibration, perturbation, cell-size-scaling, and moving-interval experiments.
It is intentionally limited to experimental design, mathematical definitions,
numerical procedures, and operational outcome definitions. It does not propose
an interpretation or a figure plan.

The completed experiments fall into four methodological groups:

1. fixed-strip calibration of a one-GTPase wave-pinning model and a reduced
   Cdc42-Rac-Rho model;
2. fixed-strip tests of simultaneous polarity sites and sequential polarity
   reversal;
3. concentration-preserving isometric cell-volume scaling (the WGD null);
4. moving-interval sweeps over reference length, endpoint mechanics, cue
   protocol, and tension feedback.

The workbook `EXPERIMENT_SUMMARIES.xlsx` contains descriptive summaries of the
completed outputs. `FULL_DATA_CATALOG.md` describes the condition-level and
trajectory-level data available if those summaries are insufficient.

## Geometry and state variables

### One-dimensional strip

The cell is represented by a one-dimensional interval with endpoints
`a(t)` and `b(t)`, length

```text
L(t) = b(t) - a(t),
```

and an affine material grid between the endpoints. Fixed-strip experiments set
the mechanics activation time to `None`, so `a` and `b` do not move. Moving
experiments activate mechanics from `t=0`.

The strip represents a membrane-associated compartment and an adjacent
cytoplasmic or combined-inactive compartment. For the single-species model,
`h = A_c/s_m` is the effective cytoplasmic volume per represented membrane
area. The conserved amount is therefore

```text
M = integral_a^b [u(x,t) + h v(x,t)] dx.
```

For the reduced three-species model, each active/inactive pair is represented
in the same one-dimensional concentration convention and has its own conserved
amount:

```text
M_C = integral (C + C_i) dx,
M_R = integral (R + R_i) dx,
M_rho = integral (rho + rho_i) dx.
```

The six-field model does not explicitly resolve inactive membrane-bound and
inactive cytosolic pools. Consequently it has no explicit surface-to-volume
reservoir parameter corresponding to `h`; such a representation would require
a nine-field membrane/exchange model.

## Single-species chemistry

### Variables and reaction

The model contains active, slowly diffusing membrane GTPase `u(x,t)` and
inactive, rapidly diffusing cytoplasmic GTPase `v(x,t)`:

```text
partial_t u = D_u partial_xx u + J(u,v),
partial_t v = D_v partial_xx v - J(u,v)/h,

J(u,v) = delta v [kappa + Gamma u^n/(K^n + u^n) + S(x,t)] - delta u.
```

Here `S(x,t)` is an externally imposed activation cue. The primary controls are

```text
rho    = D_v/D_u,
Lambda = L sqrt(delta/D_u),
kappa  = k_0/delta,
Gamma  = gamma/delta,
mu     = mean(u + h v)/K,
h      = cytoplasmic depth,
n      = Hill coefficient.
```

Unless otherwise noted, `K=1`, `delta=1 s^-1`, `rho=100`, `kappa=0.067`,
`Gamma=1`, and `n=2`. For a requested physical/model length `L`, the active
diffusivity is inferred from `Lambda` as

```text
D_u = delta (L/Lambda)^2,
D_v = rho D_u.
```

The initial homogeneous state is the lowest stable root on the conservation
line `u + h v = mu K`. When two active regions are supplied as an initial
condition, the same amount is removed from `v` locally (divided by `h`) so the
represented material is conserved pointwise.

### Parameter lineage

The reaction form and baseline dimensionless regime are based on the
wave-pinning lineage of Mori et al. (2008), with membrane/cytoplasm geometry
motivated by later bulk-surface formulations and parameter ranges compared
against Giese et al. (2015), Chiou et al. (2018), and Buttenschön et al. (2020).
Published one-dimensional total amounts are not transferred literally into the
strip geometry. Instead, the implementation uses the dimensionless controls
above and explicitly represents the reservoir factor `h`.

## Three-species chemistry

### Six-field Holmes reduction

The model contains active and combined-inactive pools for Cdc42, Rac, and Rho:

```text
(C, C_i), (R, R_i), (rho, rho_i).
```

With the inactive-accessibility factor equal to one, the implemented active
reaction rates are

```text
dC/dt = I_C/[1 + (rho/a_1)^n] * C_i/C_tot - C,

dR/dt = [I_R1 + alpha C + S_R(x,t)] * R_i/R_tot - R,

drho/dt = I_rho/[1 + (R/a_2)^n] * rho_i/rho_tot - rho.
```

The inactive rates are the negatives of the corresponding active rates. The
mapped totals and fixed controls are

| Quantity | Value |
|---|---:|
| `C_tot` | 2.4 |
| `R_tot` | 7.5 |
| `rho_tot` | 3.1 |
| `I_C` | 2.95 s^-1 |
| `I_rho` | 6.6 s^-1 |
| `a_1` | 1.25 |
| `a_2` | 1.0 |
| Cdc42-to-Rac coupling `alpha` | 0.65 s^-1 in standardized screens |
| active diffusivity `D_m` | 0.1 um^2 s^-1 |
| inactive diffusivity `D_i` | 50 um^2 s^-1 |
| diffusion ratio | 500 |

The Rac basal control is reported as

```text
kappa_R = I_R1/(R_tot d_R),
```

with `d_R=1 s^-1`. The standardized experiments vary `kappa_R`, `Lambda`, and
the inhibitory Hill exponent `n`. Physical length is

```text
L = Lambda sqrt(D_m/d_R) = Lambda sqrt(0.1) um.
```

Each parameter condition is initialized at a numerically relaxed homogeneous
equilibrium. A localized stimulus enters additively in the Rac activation
rate. This reduced network follows Holmes et al. (2012): Rho inhibits Cdc42,
Cdc42 promotes Rac, and Rac inhibits Rho. The broader model lineage also
includes Marée et al. (2006, 2012), Jilkine et al. (2007), Spill et al. (2016),
and Comer et al. (2021).

## External cues and initial perturbations

### Gaussian cue

All localized activation cues use relative material position
`xi=(x-a)/(b-a)` and the form

```text
S(x,t) = A exp[-0.5 ((xi-xi_0)/w)^2]
```

during the specified time interval and zero otherwise. Thus cue centers and
widths scale with the current interval rather than remaining fixed in laboratory
coordinates.

The common fixed-strip single cue has amplitude 1, center `xi_0=0.1`, width
`w=0.05`, and acts from 10 to 30 s. Sequential reversal adds an equal cue at
`xi_0=0.9` from 150 to 170 s. Moving-interval reversal instead applies the
second cue from 250 to 270 s and follows the trajectory to 500 s.

### Simultaneous dual cues

The moving-interval `dual_equal` protocol applies equal left and right cues at
relative positions 0.1 and 0.9 from 10 to 30 s. `dual_biased` multiplies the
left amplitude by `1+0.01` and the right amplitude by `1-0.01`.

### Two seeded domains

The fixed-strip competition tests do not apply time-limited cues. They initialize
two Gaussian active enrichments centered at relative positions 0.2 and 0.8,
with width 0.04 and nominal amplitude 2. Left and right amplitudes are
`A(1+epsilon)` and `A(1-epsilon)`, for
`epsilon = 0, 0.001, 0.01, 0.05`. The transferred material is removed from the
corresponding inactive pool.

## Numerical discretization and integration

### Conservative finite volumes

The interval is discretized into nodal control volumes. Endpoint control
volumes have half the adjacent grid spacing; interior control volumes have half
the sum of the adjacent spacings. The integrated state stores chemical amounts
rather than concentrations. Concentrations are reconstructed by dividing by the
current control-volume weights.

For concentration `c` and diffusivity `D`, the interface flux is

```text
F_(i+1/2) = -D (c_(i+1)-c_i)/(x_(i+1)-x_i).
```

Flux differences update control-volume amounts. Flux is zero at both endpoints.
Reaction rates are multiplied by the current control-volume weights. This form
accounts for changing control-volume sizes without adding a separate advection
term and directly conserves represented material up to solver error.

### Time integration

The method of lines is integrated with SciPy `solve_ivp`, using the BDF method,
relative tolerance `2e-5`, absolute tolerance `1e-8`, and no random initial
noise in production experiments. Fixed-strip calibration and WGD experiments
normally use 81 nodes, output every 2 s, and maximum internal step 0.5 s.
Moving-interval experiments output every 1 s, use maximum step 0.25 s, and set

```text
nodes = max(81, ceil(L_0/0.25) + 1).
```

The three-species simulations enable a reaction-diffusion Jacobian sparsity
pattern for BDF. Chemical rows couple local reactions, nearest-neighbor
diffusion, endpoint-derived mechanical quantities, and—when present—the local
spatial-tension field. Endpoint rows remain dense because endpoint force
mappings may inspect all fields.

### Conservation diagnostics

For every trajectory, the maximum relative change from the initial conserved
amount is reported. The moving-interval merge rejects a campaign if a species
conservation error exceeds `1e-7` or if the interval falls below 0.25 um.
WGD extended merges use a conservation threshold of `1e-6`. The workbook
reports observed errors separately from biological outcomes.

## Moving-interval mechanics

### Kelvin-Voigt endpoint model

Only the two endpoints are mechanical degrees of freedom. Interior nodes follow
the affine interpolation of endpoint motion; the model is not an interior
active-gel or continuum force-balance model.

Let `V_l` and `V_r` be endpoint velocities, `F_l` and `F_r` the applied
laboratory-frame endpoint forces, `k` the spring stiffness, `eta` the dashpot
viscosity, `L_0` the reference length, and `gamma_l=gamma_r=1` the endpoint
frictions. The internal force is

```text
T = k(L-L_0) + eta(V_r-V_l),
```

and force balance is

```text
gamma_l V_l = F_l + T,
gamma_r V_r = F_r - T.
```

These two linear equations are solved directly at each right-hand-side
evaluation. The spring contribution is `k(L-L_0)` and the dashpot contribution
is `eta(V_r-V_l)`. The reported elastic tension is the positive part of the
spring contribution; compression and dashpot force are excluded from that
quantity.

### Mechanical cases

For reference length `L_0` in micrometers,

```text
soft:  k = 1/L_0,    stiff: k = 100/L_0,
eta = 20/L_0,
gamma_l = gamma_r = 1.
```

This is a 100-fold stiffness contrast. Scaling `k` and `eta` inversely with
reference length normalizes the nominal mechanical response across the length
sweep.

### Active endpoint force

For the single-species model, the active-GTPase endpoint value is passed through
a Hill function with exponent 4 and half level `K=1`. The homogeneous resting
Hill value is subtracted, the result is multiplied by gain 0.08, and the force
is directed outward at each endpoint.

For the three-species model, endpoint force is outward Rac activity minus inward
Rho activity. Both use gain 0.08 and Hill exponent 4; the half levels are the
relaxed homogeneous Rac and Rho values. Laboratory-frame signs are assigned so
Rac at the leading edge and Rho at the trailing edge can contribute to motion
in the same direction.

## Tension feedback

Tension feedback adds a nonnegative deactivation rate to the single active
GTPase or specifically to active Rac. Deactivated material is transferred to
the corresponding inactive pool, preserving that GTPase total.

### Global feedback

The global mode uses the positive elastic spring force uniformly:

```text
r_added(t) = alpha_t max[0, k(L-L_0)].
```

Compression and dashpot force do not contribute. The tested strengths are
`alpha_t = 0, 5, 20`.

### Spatial feedback

The spatial mode integrates an excess-tension signal `theta(x,t)` with no-flux
boundaries:

```text
partial_t theta = D_T partial_xx theta - theta/tau_T + beta P(x,t),
r_added(x,t) = alpha_t max(theta,0).
```

`P` is active GTPase in the single-species model and active Rac in the
three-species model. The implemented parameters are

| Quantity | Value |
|---|---:|
| `D_T` | 10 |
| `tau_T` | 10 s |
| `beta` | 0.01 |
| baseline tension | 0 |
| communication length `sqrt(D_T tau_T)` | 10 um |
| `alpha_t` | 0, 5, 20 |

The spatial signal is phenomenological and has the units implied by these
couplings; it is not an independently calibrated material stress field.

## Operational outcome definitions

### Common active-field classification

The primary active field is `u` for the single-species model and active Rac for
the matched three-species endpoint analyses. For a profile `p(x)`, normalized
polarization amplitude is

```text
P = [max(p)-min(p)]/mean(p).
```

A profile qualifies as polarized when `P >= 0.5`. A high domain is a connected
interval satisfying

```text
p > min(p) + 0.5 [max(p)-min(p)].
```

High domains are counted only after the amplitude criterion is satisfied. A
small numerical ripple in an effectively homogeneous profile therefore counts
as zero domains.

Late stationarity is assessed over the final 50 s. Consecutive saved profiles
are compared by root-mean-square change divided by mean absolute late activity.
The largest normalized change is used. A value above 0.01 is nonstationary.

Terminal labels are assigned as follows:

| Label | Operational definition |
|---|---|
| `homogeneous_low` | stationary, `P<0.5`, final mean active level `<K` |
| `homogeneous_high` | stationary, `P<0.5`, final mean active level `>=K` |
| `single_polarized` | stationary, `P>=0.5`, one high domain |
| `dual_polarized` | stationary, `P>=0.5`, two high domains |
| `multi_peak` | stationary, `P>=0.5`, more than two high domains |
| `transient_coarsening` | nonstationary and late domain count decreases |
| `oscillatory_dynamic` | nonstationary without a late domain-count decrease |

`time_to_polarization` is the first time at which `P>=0.5` remains true for
10 s. “Terminal” and “persistent” refer to the finite observation window, not
formal asymptotic stability.

### Endpoint contrast and reversal

Endpoint contrast is

```text
C_end = (p_right-p_left)/(|p_right|+|p_left|).
```

A fixed-strip reversal run is eligible only if, immediately before the right
cue, it has exactly one polarized domain and `C_end <= -0.5`. An eligible run is
classified as:

| Label | Endpoint definition |
|---|---|
| `reversed` | one terminal domain and `C_end >= 0.5` |
| `retained_original_pole` | one terminal domain and `C_end <= -0.5` |
| `lost_polarity` | zero terminal domains |
| `multiple_polarized_domains` | at least two terminal domains |
| `ambiguous` | eligible but none of the preceding rules apply |
| `not_initially_left_polarized` | eligibility criterion not met |

Reversal time is measured from the beginning of the opposing cue to the first
single right-biased domain that persists for 10 s. Ineligible cases have no
reversal time and must not be counted as failed reversals.

### Three-species composite annotations

For each active Cdc42, Rac, and Rho profile, range contrast is
`(max-min)/(max+min)`. The spatial first moment uses a coordinate from -1 to 1
and is normalized by profile integral. Coherent polarity requires:

- Cdc42 and Rac range contrasts at least 0.2;
- Cdc42-Rac correlation at least 0.5;
- negative Cdc42-Rho correlation;
- nonzero Cdc42 moment (`|moment|>=0.01`);
- Cdc42 and Rac moments with the same sign and Rho moment with the opposite
  sign.

The homogeneous equilibrium is also linearized. The zero mode respects the
three conserved totals. Nonzero no-flux modes `m=1,...,64` add
`-D(m pi/L)^2` to the reaction Jacobian. The largest real growth rate and its
fastest mode are retained. These composite and stability quantities are
secondary annotations; the matched endpoint classifications use active Rac and
the common domain rules.

### Moving-cell polarity and motion metrics

The single-species signed polarity is endpoint contrast. The three-species
signed polarity is the mean of the Cdc42 and Rac spatial moments minus the Rho
moment:

```text
polarity = (moment_C + moment_R - moment_rho)/3.
```

Cell center and velocity are

```text
x_center = (a+b)/2,
V_center = (V_l+V_r)/2.
```

The output records net center displacement, absolute speed, signed velocity,
and polarity-projected speed `V_center sign(polarity)`. Interval means use:

| Interval | Time span |
|---|---:|
| cue | 10-30 s |
| early | 30-100 s |
| established | 100-250 s |
| reversal cue | 250-270 s |
| reorientation | 270-350 s |
| late | 350-500 s |

These intervals are calculated for every protocol; intervals named for reversal
still exist for non-reversal protocols and should be interpreted only as time
windows. Additional outputs include final/minimum/maximum length, maximum
absolute strain, elastic-tension quantiles, fraction of saved times under
positive elastic tension, tension-signal quantiles, added-deactivation-rate
quantiles, and time-integrated mean added-deactivation exposure.

The directional atlas integrates positive and negative center velocity over
time to obtain cumulative rightward and leftward travel. Direction is `left` or
`right` if that direction contributes at least 80% of total travel; otherwise it
is `both`. This classification is a descriptive reduction of the complete
velocity trajectory.

## Completed experimental designs

### Single-species local bistability

At `mu=2.2683312`, `Lambda=31.6227766017`, `h=1`, and the common fixed
kinetics, the local equation `J(u,v)=0` is solved over 301 inactive
concentrations from 0 to 3. Roots are located by bracketing on a dense active
grid and Brent root finding. Stability is assigned from the sign of the local
derivative with respect to active concentration.

### Single-species single-cue regime map

The fixed interval has model length 10. The full factorial grid contains 17
values each of `mu` (0.5 to 5) and `Lambda` (10 to 100), at `h=1,2,4`, for 867
trajectories. Each receives the common left Gaussian cue and is observed for
300 s.

### Single-species two-domain map

The grid contains 13 values each of `mu` (0.5 to 5) and `Lambda` (10 to 100),
at `h=1,2,4` and four `epsilon` values, for 2,028 trajectories. Domain counts
are recorded at 0, 100, and 300 s.

### Single-species reversal map

The reversal grid matches the 17 by 17 by 3 calibration grid (867 trajectories).
The left cue acts from 10 to 30 s and the right cue from 150 to 170 s; outcomes
are evaluated at 300 s.

### Three-species tier-1 calibration

The grid contains 13 `Lambda` values from 31.6227766017 to 94.8683298051,
13 `kappa_R` values from 0 to 0.08, and `n=2,3,4`, for 507 trajectories. All
other kinetics, totals, and diffusivities are fixed. A left active-Rac cue acts
from 10 to 30 s and trajectories are observed for 300 s. The output includes
the common active-Rac classification, composite Cdc42-Rac-Rho annotations, and
linear stability quantities.

### Isometric WGD-null scaling

Each WGD-null experiment evaluates matched volume multipliers
`V=0.5,1,2`. The linear multiplier is

```text
s = V^(1/3).
```

Physical length and dimensionless length scale as `L -> sL` and
`Lambda -> s Lambda`. Reaction kinetics, diffusivities, concentrations, cue
amplitude and timing, and cue fractional position and width are fixed.

For the single-species model, `h -> s h`. The resting baseline concentrations
`u_0` and `v_0` are preserved, so the scaled abundance is

```text
mu_scaled = (u_0 + s h v_0)/K.
```

The endpoint screen evaluates the 867 baseline `(mu,Lambda,h)` conditions at
three volumes: 2,601 simulations. The extended two-domain experiment uses the
13 by 13 by 3 grid, four asymmetries, and three volumes: 6,084 simulations. The
extended reversal experiment uses the 17 by 17 by 3 grid and three volumes:
2,601 simulations.

For the three-species model, concentrations are preserved directly and no `h`
is present. The endpoint screen evaluates 507 baseline
`(Lambda,kappa_R,n)` conditions at three volumes: 1,521 simulations. The
extended two-domain and reversal experiments have the same respective counts
as the single-species extensions: 6,084 and 2,601 simulations.

WGD comparisons are paired by their baseline molecular condition. Reported
transitions compare the exact 1x outcome with its matched 0.5x or 2x outcome;
aggregate prevalence by volume does not preserve these pairings.

### Moving-interval focused-length sweeps

Four matched campaigns combine:

- chemistry: single species or three species;
- tension mode: global or spatial;
- reference length: 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, or 50 um;
- mechanics: soft or stiff;
- protocol: single, dual equal, dual 1% biased, or reversal;
- tension-feedback strength: 0, 5, or 20.

Each campaign therefore contains 240 trajectories and the four campaigns
contain 960 trajectories. Each is followed for 500 s. The single-species
chemistry fixes `mu=2.1875`, `h=1`, `D_u=0.05`, `rho=100`, `kappa=0.067`,
`Gamma=1`, and `n=2`. Physical length is varied at fixed concentrations and
fixed `h`; this is not isometric volume scaling. The three-species chemistry
uses `kappa_R=0.04` and `n=3` with the mapped fixed controls above.

## Methodological boundaries

- These are numerical regime maps and finite-duration perturbation experiments,
  not formal bifurcation or infinite-time stability analyses.
- Parameter-grid counts describe the sampled grid and are not estimates of
  biological probabilities.
- The WGD null isolates geometric scaling while preserving local concentrations;
  it is not fitted to a particular WGD cell line.
- The moving mechanics are endpoint Kelvin-Voigt mechanics with affine interior
  motion, not a distributed active material model.
- Tension feedback is phenomenological and its tested strengths are model
  controls rather than fitted molecular rate constants.
- The six-field Holmes reduction combines inactive pools and therefore cannot
  represent size-dependent membrane/cytosol partitioning explicitly.
- Outcome thresholds are operational analysis definitions. Sensitivity to those
  thresholds was not supplied as a separate completed campaign.

## Primary implementation sources

The consolidated description above was checked against the implemented source
modules for the moving finite-volume solver, force laws, tension feedback,
single-species chemistry and observables, Holmes chemistry and observables,
WGD transformations, and focused-length output calculations. Exact source and
data paths are listed in `FULL_DATA_CATALOG.md`.

## Literature cited for model provenance

- Mori Y, Jilkine A, Edelstein-Keshet L. Wave-pinning and cell polarity from a
  bistable reaction-diffusion system. *Biophysical Journal*. 2008;94:3684-3697.
  doi:10.1529/biophysj.107.120824.
- Giese W, Eigel M, Westerheide S, Engwer C, Klipp E. Influence of cell shape,
  inhomogeneities and diffusion barriers in cell polarization models.
  *Physical Biology*. 2015;12:066014.
- Chiou J-G et al. Principles that govern competition or co-existence in
  Rho-GTPase driven polarization. *PLoS Computational Biology*.
  2018;14:e1006095. doi:10.1371/journal.pcbi.1006095.
- Buttenschön A, Liu Y, Edelstein-Keshet L. Cell size, mechanical tension, and
  GTPase signaling in the single cell. *Bulletin of Mathematical Biology*.
  2020;82:28. doi:10.1007/s11538-020-00702-5.
- Marée AFM et al. Polarization and movement of keratocytes: a multiscale
  modelling approach. *Bulletin of Mathematical Biology*. 2006;68:1169-1211.
  doi:10.1007/s11538-006-9131-7.
- Jilkine A, Marée AFM, Edelstein-Keshet L. Mathematical model for spatial
  segregation of the Rho-family GTPases based on inhibitory crosstalk.
  *Bulletin of Mathematical Biology*. 2007;69:1943-1978.
  doi:10.1007/s11538-007-9200-6.
- Holmes WR, Lin B, Levchenko A, Edelstein-Keshet L. Modelling cell polarization
  driven by synthetic spatially graded Rac activation. *PLoS Computational
  Biology*. 2012;8:e1002366. doi:10.1371/journal.pcbi.1002366.
- Marée AFM, Grieneisen VA, Edelstein-Keshet L. How cells integrate complex
  stimuli. *PLoS Computational Biology*. 2012;8:e1002402.
  doi:10.1371/journal.pcbi.1002402.
- Spill F et al. Effects of 3D geometries on cellular gradient sensing and
  polarization. *Physical Biology*. 2016;13:036008.
  doi:10.1088/1478-3975/13/3/036008.
- Comer S et al. On the influence of cell shape on dynamic reaction-diffusion
  polarization patterns. *PLoS ONE*. 2021;16:e0248293.
  doi:10.1371/journal.pone.0248293.
