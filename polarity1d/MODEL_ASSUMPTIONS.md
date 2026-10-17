# Model provenance and explicit assumptions

## Frozen reproduction scope

The Jilkine implementations and existing qualification assets are frozen as a
documented equation-level, qualitative reconstruction. Exact figure matching
is not a project objective: the source does not specify a reproducible
symmetry-breaking mechanism and omits other details that materially affect the
reported branches. In particular, a winner selected from a perfectly symmetric
cue by floating-point or solver perturbations is a numerical branch selection,
not a biological observation. Production probability estimates must use
declared cue asymmetry or seeded initial-condition perturbations.

The optimized De Belly finite-volume port is accepted as the Phase 1 screening
implementation. Its documented discretization differences remain in scope as
numerical assumptions, but a new native-Gridap trajectory cross-check is not
required before the production campaign.

## Jilkine repertoire

Equations, dimensional parameters, and stimulus definitions follow
`../jilkine_ref.md`. All three two-field models impose no-flux boundaries and
apply stimulus as inactive→active conversion, so total material is conserved.
LEGI consumes the same scalar profile as its signal.

### Wave-pinning

Defaults are `Du=0.1`, `Dv=10`, `k0=0.067`, feedback `c=1`, `K=1`, and
`d=1`. The review omits the Figure 2 total pool. We use the homogeneous fixed
point already present in this repository's Mori report, `u0=0.2683 µM` and
`v0=2 µM`. Total concentration is held fixed as length changes, so total
material grows with domain length.

### Otsuji

The implementation uses Jilkine's time-rescaled values (`a1=25`, `a2=0.7`,
`s=1`, total `u+v=2`) rather than the original Otsuji clock. The default
homogeneous state is calculated analytically from `f=0`:

```text
v0 = total / (a2*s*total + 1)^2
u0 = total - v0
```

### Goryachev

Jilkine states that cytoplasmic Cdc24–Bem1 (`Ec`) was held constant but does
not give its value or the initial Cdc42 abundance. The documented defaults are
therefore reproduction assumptions: `Ec=1`, `u0=1`, and
`v0=c/[Ec*(a*u0+b)]`, which puts the system exactly on its positive homogeneous
fixed point. Both `Ec` and the initial concentrations are constructor
parameters and should be treated as sensitivity dimensions before biological
claims are made.

### LEGI

Activator and response diffusion default to zero, inhibitor diffusion to
`10 µm²/s`, matching the review's main local-A/local-R implementation. Setting
`da=dr=0.1` permits the stated diffusion control. The response is not a
mass-conserved intrinsic pattern former and should be interpreted as sensing,
not persistent biochemical polarity.

## Stimulation and noise

`LocalizedPulse` normalizes the cosine cap to a specified maximum at the cell
edge; at 10% width this is algebraically the Jilkine profile. `LinearGradient`
uses a physical slope, so `slope=0.01` at `L=10` gives the published maximum
rate `0.1/s`. `ReversingGradient` and `AlternatingGradient` implement the WP1
extensions.

The current qualification uses the literal symmetric two-cue protocol at
`L=10 µm`: equal reflected stimuli, homogeneous model initial states, and no
imposed perturbation or cue asymmetry. Its selected in-range stimulus values are
`S=0.2`, `t1=20 s`, and `t2=40 s`. Any eventual left/right selection is therefore
floating-point/solver-seeded rather than prescribed by the protocol. This is a
diagnostic of the deterministic implementation, not evidence for a biological
source of symmetry breaking. The same protocol is reused without retuning for
the `L=1, 3, 17, 25, 39 µm` profiles at `T=2000 s`; space is plotted as `x/L`,
while model observables retain their raw scales.

The third qualification row deliberately restores the earlier model-specific,
figure-guided assumptions for comparison at the same lengths and `T=2000 s`:

- WP: `u0=0.49`, `v0=1.8801794752073646`, transient `S=0.04`, and 2% cue
  asymmetry.
- GOR: equal transient `S=1e-4` and conservative initial noise `A=0.1`, seed 11.
- OT: transient `S=0.04` and 2% cue asymmetry.
- LEGI: `DA=0.1`, `DR=0`, `R0=0.1`, and persistent equal `S=2e-4`.

This row is explicitly a qualitative figure-guided reconstruction, not the
shared-protocol reproduction represented by rows one and two.

The review's printed noise equation conflicts with its zero-mean description.
`noisy_initial_state` therefore uses seeded zero-mean uniform noise and moves
material between the active/inactive pair pointwise. This is conservative and
explicitly avoids negative pools. For a model whose observable homogeneous
state is zero (the current LEGI default), signed zero-mean noise would create
negative concentrations; the helper instead shifts that realization to be
one-sided and nonnegative. This exception is recorded rather than hidden.

## 2026 de Belly 1D mechanochemical model

The source of truth is the retained Julia code at
`../literature_models/debelly_2026_mechanochemical/src/debelly/one_d/`. The
membrane-cortex mechanics descend from the 2023 work, while the Rac-Rho-MCA
mechanochemical polarity circuit and runnable 1D implementation are attributed
to the 2026 paper.
The optimized model retains:

- Rac and Rho mutual antagonism;
- bound and unbound MCA with published on/off rates;
- membrane displacement and cortical velocity force balances;
- Rho-driven active stress;
- tension→Rho and MCA depletion→Rac feedback;
- Rac-dependent leading-edge protrusion;
- independent `rac_to_rho_coupling` and `rho_to_rac_coupling` multipliers.

The older arm multipliers remain API compatibility adapters, but they do not
feed the common-protocol campaign. That campaign directly changes only
`rho_mechanical_feedback` (β=0), only `rac_mechanical_feedback` (α=0), or
`protrusion_speed`, both feedbacks, and `active_stress` (local inhibition).

For legacy experiments, the arm multipliers behave as follows:

- `rac_to_rho_coupling` scales Rac-driven protrusion and tension-dependent Rho;
- `rho_to_rac_coupling` scales Rho active stress and MCA-dependent Rac.

Setting one to zero removes that complete long-range arm while leaving basal
Rac/Rho antagonism intact. Setting both to zero is the local biochemical
control. `opto_target="rac"` (default) or `"rho"` maps the common scalar
protocol to either released basal-activation input.

This is a deforming mechanochemical model, but not a whole-cell migration
model. The `published_one_sided` assay mode reproduces the released geometry:
the right boundary is anchored while Rac sets a protrusive velocity at the
left boundary. It is retained for native Julia/Gridap trajectory validation.
The optional `center_boundary_motion=False` setting in the symmetric assay
retains independent Rac/tension switches at both endpoints but skips both
mean-displacement and mean-cortical-velocity subtraction. Each end retains its
own raw displacement and protrusion velocity. An inactive end starting at zero
therefore remains at zero displacement without being explicitly anchored.
The diagnostic L=21 run changes only this setting relative to the retained
Figure 2 task, uses its original 0.5-s timestep, and has no substrate clutch or
added translation. Regenerate it with `scripts/run_debelly_L21_uncentered.py`,
then `scripts/build_debelly_L21_deformed_kymographs.py --uncentered`.
Its left endpoint reaches -1.359 µm and its right endpoint remains at 21 µm.
The Jilkine production campaign instead uses the explicitly labelled
`reflection_symmetric_extension`: the same tension-inhibited Rac protrusion
law is applied with opposite outward normals at the left and right endpoints,
both boundary displacements use the same relaxation rule, and their mean
displacement and velocity are subtracted. **Boundary audit (2026-10-08):** this
in-loop operation changes local strain/tension and is not a harmless removal of
rigid translation. The independently tested correction disables it and retains
each endpoint's own Rac/tension command; see
[the resolution audit](manuscript/data/runs/debelly_boundary_resolution_20261008/README.md).
An inactive endpoint remains constrained in that corrected prescribed-boundary
assay. No stress-free passive end or Rho-driven boundary contraction is added. Membrane displacement, cortical velocity, and tension evolve
internally; Rho-generated active stress drives cortical flow, and bound/unbound
MCA is transported by the membrane/cortical velocities. Rac and Rho evolve by
reaction and diffusion; MCA is the mechanically advected species. The reversal
and tracking assays therefore measure redistribution of intracellular polarity
on a deforming reference cell, not turning of a migrating cell.

Differences from Gridap are explicit numerical choices: cell-centred finite
volumes replace P1 finite elements; transport is conservative upwind instead
of unstabilized Galerkin advection; the released staggered MCA exchange order
is retained; the prescribed cortical and membrane velocities supply the MCA
flux at the moving left boundary; and boundary-lost MCA is restored as a
spatially uniform source split by the resting exchange fractions. This is the
finite-volume counterpart of the released global MCA-restoration term and
preserves the requested total abundance exactly. Production assays use the
timestep-independent homogeneous algebraic Rac-Rho fixed point. Direct
trajectory validation alone selects `initialization_mode="julia_warmup"` to
reproduce the released 251-step finite warm-up. It also selects
`input_sampling_mode="julia_element_right_edge"`, because the Julia
common-protocol adapter assigns one forcing value per element at `j*dx`, while
production finite volumes evaluate physical inputs at their cell centres.
Julia's exported tension is a backward difference of interior displacement
DOFs labelled at the right node; cross-validation reconstructs this same
diagnostic from Python displacement and separately retains the direct
centred-gradient discrepancy. Native validation failures are failed-audit
results, not waived equivalence.
