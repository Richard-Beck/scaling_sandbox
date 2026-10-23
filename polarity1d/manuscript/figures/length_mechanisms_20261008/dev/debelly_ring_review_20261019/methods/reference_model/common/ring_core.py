"""Conservative variable-cell deBelly ring operators, with prescribed point cuts.

This module has no threshold, peak detector, event timing or gap tracking.  Face
``i`` is between left cell ``(i-1) % N`` and right cell ``i`` in one fixed global
orientation.  A key in ``state['cuts']`` opens that face mechanically.  Each cut
stores ``Uleft, Uright, raw_left, raw_right``: displacement controllers on its
left-cell and right-cell flanks.  The left flank is an interval's *right* end;
its positive outward cortical command therefore has the opposite sign to the
right flank's native left-end command.  All vectors remain globally oriented.

Public API::

    core = RingCore(frozen_native_module, parameters, widths)
    new_state, ledger = core.step(state, dt=.5, external=cue, target_mass=6.2)

State fields are ``R, Q, b, c, u, v, T`` and ``cuts``.  Native long field names
are accepted as input aliases; output uses short names.  Inputs are not mutated.
The optional ``cuts=`` argument replaces state['cuts'] for this step.  Opening,
closing, initializing controllers and conservative cell splitting belong to the
caller. ``split_cell`` and ``face_trace`` support that work without event logic.

Spatial scheme
--------------
For cell widths hL,hR, every connected diffusive or mechanical link has series
half-cell conductance ``D / ((hL+hR)/2)``.  Its common linear trace weights the
left/right averages by hR/hL.  Open mechanical faces replace that link by two
Dirichlet half-cell conductances ``2 D/h``.  The cell gradient is the difference
of reconstructed right and left face values divided by its width.  Displacement
uses its prescribed open-face value; Rho uses one-sided linear extrapolation
from the nearest connected cell on its own interval.  An isolated one-cell
interval has zero Rho slope.  These rules are linear-exact on unequal cells and
reduce exactly to native ``_gradient`` on uniform intervals with >=2 cells,
including native Dirichlet ghost gradients for tension.  They do not claim
second-order accuracy at abrupt cell-width changes for curved profiles.

Connected advection uses one globally conservative upwind flux with the same
linearly reconstructed face velocity.  At open faces it retains the *signed*
native own-cell donor convention independently on each flank, including inward
velocities.  Bound diffusion stops there; c,R,Q diffusion remains reciprocal.
MCA stages, global uniform replacement, old-Q cortical forcing and staggered
R-then-Q updates follow the frozen mixed implementation.  The sole deliberate
time-refinement change is native displacement relaxation ** (dt / 0.5).

All matrices are cyclic tridiagonal and solved in O(N) using banded
Sherman--Morrison reduction; no dense N-by-N assembly occurs for N >= 3.
"""
from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

import numpy as np
from scipy.linalg import solve_banded


FIELDS = ('R', 'Q', 'b', 'c', 'u', 'v', 'T')
NATIVE_FIELDS = dict(R='rac', Q='rho', b='mca_bound', c='mca_unbound',
                     u='displacement', v='velocity', T='tension')


def face_trace(values, widths, face: int) -> float:
    """Common zero-resistance half-cell trace, in fixed global ring order."""
    z = np.asarray(values, dtype=float)
    h = np.asarray(widths, dtype=float)
    right = int(face) % len(h)
    left = (right - 1) % len(h)
    return float((h[right] * z[left] + h[left] * z[right])
                 / (h[left] + h[right]))


def split_cell(widths, state: Mapping[str, Any], index: int, fraction: float = .5):
    """Split one control volume without a concentration or material impulse.

    ``fraction`` of the parent's width goes to its globally left child. All
    seven child field averages equal the parent, preserving every integral.
    Existing face indices shift as needed.  No new cut is opened.  Returns
    ``(new_widths, new_state, new_internal_face)``.  The caller can initialize
    both new flank controllers to ``face_trace(new_state['u'], new_widths,
    new_internal_face)`` before opening.  Retaining a split face after healing
    needs no remapping.
    """
    h = np.asarray(widths, dtype=float)
    size = len(h)
    index = int(index)
    if not 0 <= index < size or not 0 < fraction < 1:
        raise ValueError('split needs an existing cell and 0 < fraction < 1')
    if h.ndim != 1 or np.any(h <= 0) or not np.all(np.isfinite(h)):
        raise ValueError('cell widths must be finite and positive')
    left_width = float(h[index] * fraction)
    new_h = np.insert(h, index + 1, h[index] - left_width)
    new_h[index] = left_width
    result = {k: deepcopy(v) for k, v in state.items()
              if k not in set(FIELDS) | set(NATIVE_FIELDS.values()) | {'cuts'}}
    for name in FIELDS:
        values = np.asarray(state[name] if name in state else state[NATIVE_FIELDS[name]],
                            dtype=float)
        if values.shape != (size,):
            raise ValueError(f'{name} must have shape ({size},)')
        result[name] = np.insert(values, index + 1, values[index])
    result['cuts'] = {(int(face) + (int(face) > index)): deepcopy(cut)
                      for face, cut in state.get('cuts', {}).items()}
    return new_h, result, index + 1


class RingCore:
    """Fixed mesh operator; instantiate again after a conservative mesh split."""

    def __init__(self, native, parameters, widths):
        self.native = native
        self.p = parameters
        self.h = np.asarray(widths, dtype=float).copy()
        if self.h.ndim != 1 or len(self.h) < 1 or not np.all(np.isfinite(self.h)) or np.any(self.h <= 0):
            raise ValueError('cell widths must be a nonempty positive finite vector')
        self.N = len(self.h)
        self.C = float(np.sum(self.h))
        self.left = (np.arange(self.N) - 1) % self.N
        self.distance = .5 * (self.h[self.left] + self.h)
        self.threshold = 1.3 * native._homogeneous_rac_rho(parameters, .5)[0]

    def _opened(self, cuts):
        opened = np.zeros(self.N, dtype=bool)
        for face in cuts:
            face = int(face)
            if face < 0 or face >= self.N:
                raise ValueError('cut face index is outside this mesh')
            opened[face] = True
        return opened

    def _solve(self, diagonal, left, right, rhs):
        """Solve row i: d_i x_i + l_i x_(i-1) + r_i x_(i+1)."""
        if self.N < 3:
            # Coincident periodic neighbours require summation for N=1 or 2.
            matrix = np.diag(np.asarray(diagonal, dtype=float).copy())
            for i in range(self.N):
                matrix[i, (i - 1) % self.N] += left[i]
                matrix[i, (i + 1) % self.N] += right[i]
            return np.linalg.solve(matrix, rhs)
        ab = np.zeros((3, self.N))
        ab[1] = diagonal
        ab[0, 1:] = right[:-1]
        ab[2, :-1] = left[1:]
        alpha, beta = float(left[0]), float(right[-1])
        if alpha == 0.0 and beta == 0.0:
            return solve_banded((1, 1), ab, rhs, overwrite_ab=True, check_finite=False)
        gamma = -float(diagonal[0])
        if gamma == 0.0:
            raise FloatingPointError('zero pivot in cyclic ring solve')
        ab[1, 0] -= gamma
        ab[1, -1] -= alpha * beta / gamma
        uv = np.zeros(self.N)
        uv[0], uv[-1] = gamma, beta
        solved = solve_banded((1, 1), ab, np.column_stack([rhs, uv]),
                              overwrite_ab=True, check_finite=False)
        x, z = solved[:, 0], solved[:, 1]
        return x - z * ((x[0] + alpha * x[-1] / gamma)
                        / (1.0 + z[0] + alpha * z[-1] / gamma))

    def _add_links(self, diagonal, left, right, conductance):
        """Add reciprocal physical face fluxes divided by each cell's width."""
        # Conductance[face] joins self.left[face] to face. roll(-1) supplies
        # each cell's globally right face. Weighted sums cancel exactly in
        # the underlying flux form, independent of cell widths.
        to_left = conductance / self.h
        to_right = np.roll(conductance, -1) / self.h
        diagonal += to_left + to_right
        left -= to_left
        right -= to_right

    def elliptic(self, diagonal, coefficient, rhs, cuts, prescribed):
        """Connected FV elasticity/viscosity, independent Dirichlet open flanks."""
        diagonal = np.asarray(diagonal, dtype=float).copy()
        rhs = np.asarray(rhs, dtype=float).copy()
        opened = self._opened(cuts)
        left, right = np.zeros(self.N), np.zeros(self.N)
        conductance = coefficient / self.distance
        self._add_links(diagonal, left, right, np.where(opened, 0., conductance))
        for face in cuts:
            il, ir = (face - 1) % self.N, face
            uleft, uright = prescribed[face]
            al, ar = 2 * coefficient / self.h[il]**2, 2 * coefficient / self.h[ir]**2
            diagonal[il] += al
            diagonal[ir] += ar
            rhs[il] += al * uleft
            rhs[ir] += ar * uright
        return self._solve(diagonal, left, right, rhs)

    def gradient(self, values, cuts, prescribed=None):
        """FV face difference; matches native gradients on uniform open arcs."""
        z = np.asarray(values, dtype=float)
        opened = self._opened(cuts)
        traces = (self.h * z[self.left] + self.h[self.left] * z) / (self.h + self.h[self.left])
        # Face arrays are indexed by the cell immediately to their right.
        left_trace = traces.copy()
        right_trace = np.roll(traces, -1)
        for face in cuts:
            il, ir = (face - 1) % self.N, face
            if prescribed is not None:
                ul, ur = prescribed[face]
            else:
                before, after = (il - 1) % self.N, (ir + 1) % self.N
                slope_l = 0. if opened[il] else (z[il] - z[before]) / self.distance[il]
                slope_r = 0. if opened[after] else (z[after] - z[ir]) / self.distance[after]
                ul = z[il] + .5 * self.h[il] * slope_l
                ur = z[ir] - .5 * self.h[ir] * slope_r
            right_trace[il], left_trace[ir] = ul, ur
        return (right_trace - left_trace) / self.h

    def transport(self, values, velocity, diffusion, dt, cuts, cap_velocity=None,
                  *, shared_diffusion=True, decay=0., source=0.):
        """Implicit conservative diffusion/upwind transport and signed cap flux.

        Returns concentration and diagnostics; cap fluxes use implicit donors.
        ``cap_velocity[face] == (global_left_flank, global_right_flank)``.
        """
        z = np.asarray(values, dtype=float)
        velocity = np.broadcast_to(np.asarray(velocity, dtype=float), (self.N,))
        source = np.broadcast_to(np.asarray(source, dtype=float), (self.N,))
        cap_velocity = {} if cap_velocity is None else cap_velocity
        opened = self._opened(cuts)
        diagonal = np.full(self.N, 1. + dt * decay)
        left, right = np.zeros(self.N), np.zeros(self.N)
        conductance = dt * diffusion / self.distance
        if not shared_diffusion:
            conductance = np.where(opened, 0., conductance)
        self._add_links(diagonal, left, right, conductance)
        speed = (self.h * velocity[self.left] + self.h[self.left] * velocity) / (self.h + self.h[self.left])
        speed = np.where(opened, 0., speed)
        positive, negative = np.maximum(speed, 0.), np.minimum(speed, 0.)
        # Face j, positive speed: donor j-1; negative speed: donor j.
        diagonal += dt * (np.roll(positive, -1) - negative) / self.h
        left -= dt * positive / self.h
        right += dt * np.roll(negative, -1) / self.h
        for face in cuts:
            vl, vr = cap_velocity.get(face, (0., 0.))
            diagonal[(face - 1) % self.N] += dt * vl / self.h[(face - 1) % self.N]
            diagonal[face] -= dt * vr / self.h[face]
        result = self._solve(diagonal, left, right, z + dt * source)
        if not np.all(np.isfinite(result)) or np.min(result) < -1e-10:
            raise FloatingPointError('ring transport produced invalid/negative concentration')
        result = np.maximum(result, 0.)
        caps = {face: (float(cap_velocity.get(face, (0., 0.))[0] * result[(face - 1) % self.N]),
                       float(-cap_velocity.get(face, (0., 0.))[1] * result[face]))
                for face in cuts}
        cap = float(sum(sum(pair) for pair in caps.values()))
        reaction = float(np.dot(self.h, source - decay * result))
        change = float(np.dot(self.h, result - z))
        return result, dict(cap_outflux=cap, cap_flank_outflux=caps,
                            reaction_rate=reaction, amount_change=change,
                            budget_error=change - dt * (reaction - cap))

    def step(self, state: Mapping[str, Any], dt: float, external=0.,
             target_mass: float | None = None, cuts=None):
        """One frozen-native staggered step with caller-prescribed open faces."""
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('dt must be positive and finite')
        old = {}
        for name in FIELDS:
            old[name] = np.asarray(state[name] if name in state else state[NATIVE_FIELDS[name]], dtype=float)
            if old[name].shape != (self.N,) or not np.all(np.isfinite(old[name])):
                raise ValueError(f'{name} must be finite with shape ({self.N},)')
        if min(old[name].min() for name in ('R', 'Q', 'b', 'c')) < 0:
            raise ValueError('chemical concentrations must be nonnegative')
        incoming = state.get('cuts', {}) if cuts is None else cuts
        cuts = {int(face): deepcopy(value) for face, value in incoming.items()}
        self._opened(cuts)
        p, n = self.p, self.native
        R, Q, b, c, u, v, T = (old[name] for name in FIELDS)
        old_mass = float(np.dot(self.h, b + c))
        target_mass = old_mass if target_mass is None else float(target_mass)
        boundary_u, commands, membrane_cap_v, cut_ledger = {}, {}, {}, {}
        relaxation = p.boundary_relaxation ** (dt / .5)
        for face, controller in cuts.items():
            il, ir = (face - 1) % self.N, face
            sensor = face_trace(R, self.h, face)
            switch = float(n.smooth_threshold(sensor, p.protrusion_switch_scale, self.threshold))
            command_l = -p.rac_to_rho_coupling * p.protrusion_speed * switch / (1 + T[il]**2 / p.protrusion_tension_scale)
            command_r = p.rac_to_rho_coupling * p.protrusion_speed * switch / (1 + T[ir]**2 / p.protrusion_tension_scale)
            old_ul, old_ur = float(controller['Uleft']), float(controller['Uright'])
            raw_l = float(controller.get('raw_left', old_ul))
            raw_r = float(controller.get('raw_right', old_ur))
            ul, ur = relaxation * raw_l + dt * command_l, relaxation * raw_r + dt * command_r
            vml, vmr = (ul - old_ul) / dt, (ur - old_ur) / dt
            boundary_u[face], commands[face], membrane_cap_v[face] = (ul, ur), (command_l, command_r), (vml, vmr)
            controller.update(Uleft=ul, Uright=ur, raw_left=ul, raw_right=ur)
            cut_ledger[face] = dict(old_Rac_trace=sensor, smooth_switch=switch,
                old_tension_left=float(T[il]), old_tension_right=float(T[ir]),
                left_command=float(command_l), right_command=float(command_r),
                left_membrane_velocity=float(vml), right_membrane_velocity=float(vmr),
                old_Uleft=old_ul, old_Uright=old_ur, old_raw_left=raw_l, old_raw_right=raw_r,
                Uleft=float(ul), Uright=float(ur), raw_left=float(ul), raw_right=float(ur))
        chi = p.linker_friction * np.maximum(b, 1e-10)
        coef = chi / dt
        un = self.elliptic(coef, p.membrane_stiffness, coef * (u + dt * v), cuts, boundary_u)
        vm = (un - u) / dt
        tn = p.membrane_stiffness * self.gradient(un, cuts, boundary_u)
        bn, bl = self.transport(b, v, p.mca_bound_diffusion, dt, cuts, commands,
            shared_diffusion=False, decay=p.mca_off, source=p.mca_on * c)
        cn, cl = self.transport(c, vm, p.mca_unbound_diffusion, dt, cuts, membrane_cap_v,
            shared_diffusion=True, decay=p.mca_on, source=p.mca_off * bn)
        temporary = float(np.dot(self.h, bn + cn))
        exchange_residual = dt * p.mca_on * float(np.dot(self.h, c - cn))
        cap = bl['cap_outflux'] + cl['cap_outflux']
        deficit = target_mass - temporary
        increment = deficit / self.C
        fraction = p.mca_on / (p.mca_on + p.mca_off)
        # Keep implicit transport-stage cells for native donor diagnostics.
        transport_b, transport_c = bn, cn
        bn = bn + fraction * increment
        cn = cn + (1 - fraction) * increment
        if min(bn.min(), cn.min()) < -1e-10:
            raise FloatingPointError('global MCA replacement produced negative concentration')
        chi = p.linker_friction * np.maximum(bn, 1e-10)
        gQ = p.active_stress * p.rho_to_rac_coupling * self.gradient(Q, cuts)
        vn = self.elliptic(chi, p.cortical_viscosity, chi * vm + gQ, cuts, commands)
        alpha, beta = n._derived_rates(p, bn, tn, np.asarray(external, dtype=float))
        rsource = p.rac_deactivation * alpha / (1 + Q**2)
        rn, rl = self.transport(R, np.zeros(self.N), p.rac_deactivation * p.rac_diffusion,
            dt, cuts, shared_diffusion=True, decay=p.rac_deactivation, source=rsource)
        qsource = p.rho_deactivation * beta / (1 + rn**2)
        qn, ql = self.transport(Q, np.zeros(self.N), p.rho_deactivation * p.rho_diffusion,
            dt, cuts, shared_diffusion=True, decay=p.rho_deactivation, source=qsource)
        new = dict(R=rn, Q=qn, b=bn, c=cn, u=un, v=vn, T=tn, cuts=cuts)
        for face, details in cut_ledger.items():
            il, ir = (face - 1) % self.N, face
            details['bound_cap_flank_outflux'] = list(bl['cap_flank_outflux'][face])
            details['unbound_cap_flank_outflux'] = list(cl['cap_flank_outflux'][face])
            details['bound_implicit_donors'] = [float(transport_b[il]), float(transport_b[ir])]
            details['unbound_implicit_donors'] = [float(transport_c[il]), float(transport_c[ir])]
            details['shared_traces'] = {name: face_trace(new[name], self.h, face) for name in ('R', 'Q', 'c')}
            # Positive diffusive flux goes from the global left cell to right.
            details['unbound_diffusive_flux_left_to_right'] = float(p.mca_unbound_diffusion * (transport_c[il] - transport_c[ir]) / self.distance[face])
            details['Rac_diffusive_flux_left_to_right'] = float(p.rac_deactivation * p.rac_diffusion * (rn[il] - rn[ir]) / self.distance[face])
            details['Rho_diffusive_flux_left_to_right'] = float(p.rho_deactivation * p.rho_diffusion * (qn[il] - qn[ir]) / self.distance[face])
        final_mass = float(np.dot(self.h, bn + cn))
        ledger = dict(mass_before=old_mass, temporary_MCA=temporary, mass_after=final_mass,
            bound_amount_before=float(np.dot(self.h, b)),
            unbound_amount_before=float(np.dot(self.h, c)),
            bound_amount_after_transport=float(np.dot(self.h, transport_b)),
            unbound_amount_after_transport=float(np.dot(self.h, transport_c)),
            bound_amount_after=float(np.dot(self.h, bn)),
            unbound_amount_after=float(np.dot(self.h, cn)),
            bound_net_binding_rate=float(np.dot(self.h, p.mca_on*c-p.mca_off*transport_b)),
            unbound_net_binding_rate=float(np.dot(self.h, p.mca_off*transport_b-p.mca_on*transport_c)),
            bound_uniform_replacement=fraction*deficit,
            unbound_uniform_replacement=(1-fraction)*deficit,
            target_mass=target_mass, global_deficit=deficit,
            uniform_total_concentration_increment=increment, replacement_bound_fraction=fraction,
            bound_cap_outflux=bl['cap_outflux'], unbound_cap_outflux=cl['cap_outflux'],
            signed_cap_loss=dt * cap, sequential_exchange_residual=exchange_residual,
            temporary_pool_budget_error=temporary - old_mass + dt * cap - exchange_residual,
            bound_transport_budget_error=bl['budget_error'], unbound_transport_budget_error=cl['budget_error'],
            mass_error=final_mass - target_mass,
            mass_relative_error=abs(final_mass / target_mass - 1) if target_mass else abs(final_mass),
            rac_reaction_rate=rl['reaction_rate'], rho_reaction_rate=ql['reaction_rate'],
            rac_budget_error=rl['budget_error'], rho_budget_error=ql['budget_error'],
            boundary_relaxation=relaxation, protrusion_switch_midpoint=self.threshold,
            cuts=cut_ledger)
        return new, ledger
