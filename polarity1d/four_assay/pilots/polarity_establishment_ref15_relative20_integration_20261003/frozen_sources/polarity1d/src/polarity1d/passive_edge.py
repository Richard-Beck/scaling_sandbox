"""Coupled membrane/cortex mechanics with a material, stress-loaded passive edge.

The other edge keeps prescribed displacement and cortical velocity. The
passive edge obeys v_b=(u_b-u_b_old)/dt and K*u_x+eta*v_x+sigma*rho=S_rest.
Face stresses and the Rho gradient use the same conservative FV convention.
"""
import numpy as np
from scipy.linalg import solve_banded


def solve_passive_edge(*, old_u, old_boundary, active_displacement, active_velocity,
                       passive_side, friction, stiffness, viscosity, active_stress,
                       rho, resting_rho, dx, dt, memory, drag):
    count=len(old_u)
    if passive_side not in ('left','right'):
        raise ValueError('Specify one passive edge')
    offset=2 if passive_side=='left' else 0
    ui=offset+2*np.arange(count)
    vi=ui+1
    ub=0 if passive_side=='left' else 2*count
    vb=ub+1
    size=2*count+2
    matrix=np.zeros((7,size))
    rhs=np.zeros((size,2))
    def add(rows,cols,values):
        matrix[3+np.asarray(rows)-np.asarray(cols),cols]+=values
    rho_faces=np.r_[rho[0],.5*(rho[:-1]+rho[1:]),rho[-1]]
    active_gradient=active_stress*np.diff(rho_faces)/dx
    ku=stiffness/dx**2
    kv=viscosity/dx**2
    add(ui,ui,friction/dt+2*ku)
    add(ui,vi,-friction)
    add(vi,vi,friction+drag+2*kv)
    add(vi,ui,-friction/dt)
    add(ui[1:],ui[:-1],-ku);add(ui[:-1],ui[1:],-ku)
    add(vi[1:],vi[:-1],-kv);add(vi[:-1],vi[1:],-kv)
    rhs[ui,0]=friction*old_u/dt
    rhs[vi,0]=active_gradient-friction*old_u/dt-memory
    rhs[vi,1]=-drag
    for index,side in ((0,'left'),(-1,'right')):
        add(ui[index],ui[index],ku)
        add(vi[index],vi[index],kv)
        if side==passive_side:
            add(ui[index],ub,-2*ku)
            add(vi[index],vb,-2*kv)
        else:
            rhs[ui[index],0]+=2*ku*active_displacement
            rhs[vi[index],0]+=2*kv*active_velocity
    index=0 if passive_side=='left' else -1
    direction=-1 if passive_side=='left' else 1
    # Scale the stress row to avoid an unnecessarily large matrix row norm.
    add(ub,ub,direction*stiffness)
    add(ub,ui[index],-direction*stiffness)
    add(ub,vb,direction*viscosity)
    add(ub,vi[index],-direction*viscosity)
    rhs[ub,0]=-.5*dx*active_stress*(rho[index]-resting_rho)
    add(vb,vb,dt);add(vb,ub,-1)
    rhs[vb,0]=-old_boundary
    response=solve_banded((3,3),matrix,rhs,check_finite=False)
    denominator=float(np.sum(drag*(1+response[vi,1])))
    if denominator<=1e-15:
        raise FloatingPointError('Passive-edge free-translation system is singular')
    speed=-float(np.sum(memory+drag*response[vi,0]))/denominator
    solution=response[:,0]+speed*response[:,1]
    u,v=solution[ui],solution[vi]
    boundary_u=[active_displacement,active_displacement]
    boundary_v=[active_velocity,active_velocity]
    slot=0 if passive_side=='left' else 1
    boundary_u[slot]=float(solution[ub]);boundary_v[slot]=float(solution[vb])
    membrane_stress=stiffness*np.array([2*(u[0]-boundary_u[0])/dx,2*(boundary_u[1]-u[-1])/dx])
    cortical_stress=viscosity*np.array([2*(v[0]-boundary_v[0])/dx,2*(boundary_v[1]-v[-1])/dx])+active_stress*rho[[0,-1]]
    total_stress=membrane_stress+cortical_stress
    traction=memory+drag*(v+speed)
    residuals={
        'passive_stress_residual':float(total_stress[slot]-active_stress*resting_rho),
        'passive_kinematic_residual':float(boundary_v[slot]-(boundary_u[slot]-old_boundary)/dt),
        'global_force_residual':float(total_stress[1]-total_stress[0]-dx*np.sum(traction)),
        'left_total_boundary_stress':float(total_stress[0]),
        'right_total_boundary_stress':float(total_stress[1]),
        'passive_membrane_stress':float(membrane_stress[slot]),
        'passive_cortical_stress':float(cortical_stress[slot]),
        'passive_membrane_velocity':float((boundary_u[slot]-old_boundary)/dt),
        'passive_cortical_velocity':float(boundary_v[slot]),
    }
    if max(abs(residuals[k]) for k in ('passive_stress_residual','global_force_residual'))>1e-6:
        raise FloatingPointError(f'Passive-edge stress balance failed: {residuals}')
    if abs(residuals['passive_kinematic_residual'])>1e-9:
        raise FloatingPointError('Passive-edge membrane/cortex velocities disagree')
    return u,v,speed,tuple(boundary_u),tuple(boundary_v),residuals


def solve_local_rac_edges(*, old_u, old_boundaries, activated, commanded_velocity,
                          friction, stiffness, viscosity, active_stress, rho,
                          resting_rho, dx, dt, memory, drag):
    """Both ends carry resting stress; local Rac selects cortex BC, in lab frame."""
    count = len(old_u)
    ui = 2 + 2*np.arange(count)
    vi = ui + 1
    matrix = np.zeros((7, 2*count+4))
    rhs = np.zeros(2*count+4)
    def add(rows, cols, values):
        matrix[3+np.asarray(rows)-np.asarray(cols), cols] += values
    ku, kv = stiffness/dx**2, viscosity/dx**2
    faces = np.r_[rho[0], .5*(rho[:-1]+rho[1:]), rho[-1]]
    add(ui, ui, friction/dt+2*ku)
    add(ui, vi, -friction)
    add(vi, vi, friction+drag+2*kv)
    add(vi, ui, -friction/dt)
    add(ui[1:], ui[:-1], -ku)
    add(ui[:-1], ui[1:], -ku)
    add(vi[1:], vi[:-1], -kv)
    add(vi[:-1], vi[1:], -kv)
    rhs[ui] = friction*old_u/dt
    rhs[vi] = active_stress*np.diff(faces)/dx-friction*old_u/dt-memory
    for slot, index, ub, direction in ((0,0,0,-1),(1,-1,2*count+2,1)):
        vb = ub+1
        add(ui[index], ui[index], ku)
        add(vi[index], vi[index], kv)
        add(ui[index], ub, -2*ku)
        add(vi[index], vb, -2*kv)
        add(ub, ub, direction*stiffness)
        add(ub, ui[index], -direction*stiffness)
        add(ub, vb, direction*viscosity)
        add(ub, vi[index], -direction*viscosity)
        rhs[ub] = -.5*dx*active_stress*(rho[index]-resting_rho)
        if activated[slot]:
            add(vb, vb, 1.)
            rhs[vb] = commanded_velocity[slot]
        else:
            add(vb, vb, dt)
            add(vb, ub, -1.)
            rhs[vb] = -old_boundaries[slot]
    solution = solve_banded((3,3), matrix, rhs, check_finite=False)
    u, v = solution[ui], solution[vi]
    boundary_u = solution[[0,2*count+2]]
    boundary_v = solution[[1,2*count+3]]
    membrane = stiffness*np.array([2*(u[0]-boundary_u[0])/dx,2*(boundary_u[1]-u[-1])/dx])
    cortex = viscosity*np.array([2*(v[0]-boundary_v[0])/dx,2*(boundary_v[1]-v[-1])/dx])+active_stress*rho[[0,-1]]
    total = membrane+cortex
    membrane_velocity = (boundary_u-np.asarray(old_boundaries))/dt
    residuals = dict(passive_stress_residual=float(np.max(np.abs(total-active_stress*resting_rho))),
        passive_kinematic_residual=float(max([abs(boundary_v[s]-membrane_velocity[s]) for s in (0,1) if not activated[s]] or [0.])),
        global_force_residual=float(total[1]-total[0]-dx*np.sum(memory+drag*v)),
        left_total_boundary_stress=float(total[0]), right_total_boundary_stress=float(total[1]),
        passive_membrane_stress=float(membrane[1]), passive_cortical_stress=float(cortex[1]),
        passive_membrane_velocity=float(membrane_velocity[1]), passive_cortical_velocity=float(boundary_v[1]))
    return u, v, 0., tuple(boundary_u), tuple(boundary_v), residuals
