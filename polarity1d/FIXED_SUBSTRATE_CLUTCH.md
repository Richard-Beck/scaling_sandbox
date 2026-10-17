# Fixed-substrate clutch extension

This is an exploratory extension of the L=21 micrometre full-model Figure 2
assay, with **coupled** substrate traction and cortical dynamics. Clutch
parameters are illustrative, not fitted or claimed to reproduce a published
clutch experiment. The previous no-substrate model and its figures remain
available. Call `simulate_debelly(..., substrate_clutch=FixedSubstrateClutch(...))`
to enable the new mechanics; omit the argument for the original code path.

## Mechanics and reference frames

The cortex equation uses the original centered cell/reference coordinates:

    eta * v_xx + sigma * Rho_x - chi * MCA_bound * (v - u_t) - traction_ref = 0

Here `v` is cortical velocity in the cell frame; `U` is rigid cell translation
velocity. Cortical velocity relative to the stationary substrate is `v + U`.
Laboratory position is reference position + membrane displacement + cell
translation. Membrane-cortex relative velocity is unchanged by translation.
The original membrane force balance, Rho/Rac circuit, MCA transport, and
symmetric endpoint protrusion law are retained. The internal mechanics remain
the original reference-domain/small-deformation formulation; this is not a
finite-strain or full spreading-cell model.

The cell is assumed freely translating with no other external force. At each
step `U` and the cortical velocity are solved together such that integrated
clutch traction is zero. A nonzero positive/negative local traction pattern
can therefore coexist with zero net traction. The imposed symmetric endpoint
kinematics are inherited from the earlier model; this extension does not
replace them with a new free-edge force law or contact-area model.

## Fixed substrate sites and bond resetting

Adhesion sites are fixed at laboratory positions `(integer + 0.5) * spacing`.
Their states are stored by these persistent integer IDs, not transported with
the cortex or the membrane. Contact is the interval between the deformed cell
endpoints. New sites entering contact start unbound and unstretched. Sites
leaving contact detach immediately; they do not retain extension if the cell
later returns. Initially covered sites have uniform zero-load equilibrium
occupancy and zero extension. Potential site density and all kinetic parameters
are uniform and independent of Rac.

Each substrate site stores bound fraction `n` and extension moment `q=n*e`:

    n_t = kon * (1 - n) - koff * n
    q_t = n * (v + U) - koff * q
    koff = koff0 * exp(abs(k * e) / Fb)
    traction_lab = density * k * q

Newly formed bonds enter with zero extension; detachment removes stored
extension. `e=q/n` is the mean extension of bound bonds, set to zero when the
bound fraction is numerically negligible. Positive traction here denotes
force exerted **on the substrate**; the cortex experiences its negative.
The force scale in the off-rate is the single-bond force `k*e`, not the
population traction. This is a mean-extension closure, not a resolved bond
extension distribution or stochastic individual-bond simulation.

## Numerical coupling

Occupancy and extension moment use backward Euler. The slip-bond off-rate is
evaluated from the previous extension. The new moment is affine in the new
slip velocity, so two banded cortex solves plus a scalar force-balance equation
give cortical velocity and translation without an iterative nonlinear solve.
The off-rate is nonlinear and updated every timestep; the velocity-affine
step is an implicit discretization of elastic bond extension, not replacement
of the clutch law by instantaneous linear drag.

Sites sample piecewise-constant cortical velocity in their current reference
bin. Substrate control-interval weights exactly tile the footprint. The same
weights transfer traction back to the cortex and enforce zero total force.
Traction is divided by reference-bin width for the reference-coordinate force
balance. The plotted physical traction instead divides by physical quadrature
coverage. Thus substrate density is per physical substrate length, while the
mechanical right-hand side is per reference length.

Contact uses the new membrane shape and the previous translation; translation
then advances by `dt*U`. This is a first-order staggered footprint update.
Raw fixed-lattice histories record that step's contact states; their footprint
can consequently differ from the displayed final endpoint positions by one
translation increment. Sites are stored only while in contact; NaN in the
history export means outside the recorded footprint, not a missing solve.

For this new branch, boundary relaxation is scaled as `0.98 ** (dt / 0.5)`
to retain the original Figure 2 relaxation timescale when refining timestep.
The no-clutch branch keeps its legacy per-step behavior. Runs use a 0.125 s
timestep, 0.2 micrometre reference bins, and a 0.025 micrometre substrate lattice.
Local length change is 100 times (deformed bin width / reference bin width - 1).
Original field amplitudes are plotted directly, without a concentration
Jacobian correction.

## Illustrative parameters and protocol

| Parameter | Value |
|---|---:|
| Potential clutch density | 1, 10, or 50 per micrometre |
| Single-clutch stiffness | 1 pN/micrometre |
| Binding rate | 0.1 per second |
| Unloaded detachment rate | 0.05 per second |
| Slip-bond rupture force scale | 0.5 pN |

All trials use the retained L=21 full-model task `00040-8bf0845d4b06`:
left Rac cue from 200 to 360 s with an 80 s withdrawal taper, amplitude 3.6,
width fraction 0.2, and 10,000 s total duration. The cue remains tied to the
cell's reference coordinates as the cell translates. Geometry and cellular
parameters otherwise match the earlier Figure 2 setup.

## Reproduction and outputs

From the repository root:

    mkdir -p polarity1d/results/campaigns/debelly_fixed_substrate_L21/logs
    sbatch polarity1d/workflow/run_debelly_fixed_substrate.slurm
    # After all tasks finish:
    conda run -n scaling_sandbox python polarity1d/scripts/check_debelly_fixed_substrate.py
    conda run -n scaling_sandbox python polarity1d/scripts/build_debelly_fixed_substrate_kymographs.py

The SLURM array includes three production cases and numerical checks. Outputs
are under `results/campaigns/debelly_fixed_substrate_L21`, with one trajectory,
raw fixed-substrate history, and JSON summary per case. `pilot` is an early
development run predating the boundary-timescale correction and is not used
in the final figures or checks.

`results/figures/Figure2_L21_fixed_substrate_clutch_comparison` compares Rac,
signed cell-frame cortical velocity, and traction across densities. Individual
case figures show ten fields plus displacement and length traces. The separate
`Figure2_L21_fixed_substrate_adhesion_history` plots occupancy and extension
moment directly on the fixed substrate lattice. PNG and PDF are generated.
Each field has shared color limits across production densities; signed fields
have symmetric limits about zero. Coordinates are laboratory micrometres,
with the exterior left blank and the moving endpoints drawn in black.

Validation includes fixed-site state persistence, stress-free new binding,
extension loss on turnover, pointwise cortical force balance, zero integrated
substrate force, stationary unstimulated cells, reversal under reflected cues,
MCA conservation, and timestep/substrate/reference-grid refinement. The
reproducible numerical comparison is saved as `numerical_checks.json`.

## Results from this run

At 10,000 s, cell-center displacements are -28.99, -25.89, and -9.84
micrometres for densities 1, 10, and 50, respectively. Negative means leftward,
toward the stimulated end. Peak total length is approximately 22.06–22.07
micrometres. These values describe these chosen parameters and do not establish
a general monotonic relationship between adhesion density and migration.

The 19 targeted model/clutch tests pass. Integrated substrate force residuals
are below 8e-15 pN, and relative MCA mass errors are below 6e-16. At density 10,
halving timestep from 0.125 to 0.0625 s changes the full-trajectory displacement
by at most 0.119 micrometres (0.46% of the production displacement scale).
Peak instantaneous differences in cell extension and cortical velocity are
0.047 micrometres and 0.0023 micrometres/s, respectively. Local occupancy has
larger pointwise differences near contact edges (up to 0.039 in bound fraction).
Over the first 1,000 s, halving substrate spacing changes displacement by less
than 0.0001 micrometres; doubling reference-grid resolution changes it by
0.0012 micrometres. These are numerical sensitivity checks for this exploratory
configuration, not a calibration or validation of the biological clutch law.

## Uncentered endpoint variant

The separate `debelly_fixed_substrate_L21_uncentered` campaign disables
`center_boundary_motion` while keeping the same clutch parameters, stimulus,
and independent Rac/tension-dependent rules at both ends. Both raw endpoint
displacements and cortical endpoint velocities are used without subtracting
their means. Endpoint displacement still relaxes toward its reference value
when protrusion is inactive; these remain prescribed boundary conditions,
not free-stress boundaries.

The clutch solve still adds a spatially uniform reference-frame translation
and solves for zero net substrate traction. Lab position is therefore still
reference position + membrane displacement + added translation. In this mode,
**added translation is not cell-midpoint displacement**: midpoint displacement
also includes half the sum of the two raw endpoint displacements. The figures
plot the actual midpoint, and summaries export both quantities. An end with
zero raw displacement can move in the laboratory through the added translation.
The saved cortical `velocity` is relative to the model reference frame, which
is no longer constrained to track the cell midpoint.

Reproduce this variant with:

    mkdir -p polarity1d/results/campaigns/debelly_fixed_substrate_L21_uncentered/logs
    sbatch polarity1d/workflow/run_debelly_fixed_substrate_uncentered.slurm
    # After the array finishes:
    conda run -n scaling_sandbox python polarity1d/scripts/check_debelly_fixed_substrate.py --uncentered
    conda run -n scaling_sandbox python polarity1d/scripts/build_debelly_fixed_substrate_kymographs.py --uncentered

New figures carry the `_uncentered` filename suffix. Previous centered figures
and trajectories are preserved. The numerical check exports include comparison
of actual midpoint displacement as well as added translation.

The uncentered production runs use a smaller timestep of 0.03125 s, with a
0.015625-s full-duration check at density 10. The preliminary 0.125-s run showed
a 0.416 µm maximum midpoint difference on halving timestep, so it was refined;
its audit is retained as `numerical_checks_dt0p125.json` and is not the final
numerical check. Spatial resolution and all physical parameters remain the
same as in the centered trials.

For the uncentered production runs, final midpoint shifts at 10,000 s are
+12.03, +3.62, and +3.39 micrometres for densities 1, 10, and 50, respectively.
Positive is rightward, away from the left-stimulated end. This differs from
the centered runs; the plots report the result of the stated boundary and
clutch model without prescribing a migration direction.

All 21 targeted tests pass, including reflection symmetry with uncentered
clutch coupling. Halving the production timestep changes the density-10
midpoint trajectory by at most 0.101 µm over 10,000 s (1.89% of peak midpoint
displacement); the final midpoint differs by 0.097 µm. Over the first 1,000 s,
doubling reference-grid resolution changes midpoint position by 0.021 µm, and
halving substrate spacing changes it by less than 0.000035 µm. Integrated
substrate force residuals remain below 1.2e-14 pN in all three production runs.

## Option 2: passive rear with shared membrane/cortex velocity

The `debelly_fixed_substrate_L21_passive_rear` campaign keeps the original
left-end displacement/protrusion prescription for the entire left-cue assay.
The right end instead obeys two conditions, in the same reference frame:

    v_right = (u_right_new - u_right_old) / dt
    K * u_x(right) + eta * v_x(right) + sigma * Rho(right) = S_rest
    S_rest = sigma * initial_resting_Rho

`S_rest` is constant (about 69.26 in model stress units for this setup). It
represents the load balancing the initial active stress; it is not a
zero-total-stress boundary. Neither position nor velocity at the passive end
is prescribed. Its previous 0.98 relaxation-to-reference update is removed;
the active end retains the same physical relaxation timescale. No centering
is applied. The right side is designated passive throughout this assay, not
selected dynamically by instantaneous Rac. The mirrored validation exchanges
the active and passive sides and reverses the cue.

The two passive conditions couple displacement and cortical velocity, so the
mechanical fields and passive endpoint are solved jointly with backward Euler.
Membrane-cortex friction and Rho stress use the previous biochemical state;
new velocities then drive the usual MCA exchange/transport update, followed
by Rac/Rho. This replaces the older staggered mechanical update in this branch
only, while retaining the same interior continuum equations. Conservative
face stresses and the same face-based Rho gradient make the integrated
mechanical force balance checkable without a boundary quadrature mismatch.
The substrate population update uses the previous full footprint in this branch.

Rigid reference-frame translation still enforces zero net substrate traction.
We additionally verify that right total boundary stress minus left total
boundary stress equals integrated substrate traction. With equal resting
loads at the two ends, the boundary forces therefore cancel; this is checked
on the computed active-edge reaction, not silently assumed from its prescribed
velocity. Independent tests evaluate both bulk force balances, the passive
stress condition, the shared edge velocity, and the unchanged active boundary
values. A homogeneous unstimulated cell remains stationary.

This version uses the same three clutch densities and other physical parameters,
with a 0.0625-s production timestep and a 0.03125-s density-10 refinement. The
spatial checks and full duration are unchanged. Reproduce with:

    mkdir -p polarity1d/results/campaigns/debelly_fixed_substrate_L21_passive_rear/logs
    sbatch polarity1d/workflow/run_debelly_fixed_substrate_passive_rear.slurm
    # After all jobs complete:
    conda run -n scaling_sandbox python polarity1d/scripts/check_debelly_fixed_substrate.py --passive-rear
    conda run -n scaling_sandbox python polarity1d/scripts/build_debelly_fixed_substrate_kymographs.py --passive-rear

Figures have a `_passive_rear` suffix. The additional
`Figure2_L21_passive_rear_boundary_audit` figure plots both raw endpoint
displacements, the membrane/cortical/total stresses at the passive edge, and
its membrane and cortical velocities. The latter curves should coincide.
Boundary stress diagnostics use the Rho values at the mechanical solve,
before the subsequent biochemical update.

Production results at 10,000 s give midpoint shifts of +8.58, -1.26, and
-3.48 µm for densities 1, 10, and 50. Positive is rightward. The minimum
lengths are 19.62, 19.46, and 19.06 µm, respectively. Right-end Rac stays
below the protrusion threshold throughout all three runs, consistent with
designating this end passive. The 23 targeted tests pass. Maximum passive
stress residuals are below 8e-12 and global boundary/traction balance errors
below 3e-11 in model stress units; the two passive endpoint velocities agree
to better than 8e-15 µm/s in these runs.
