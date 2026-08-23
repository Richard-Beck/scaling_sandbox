# Mechanochemical model on an axisymmetric cell surface.
# Simulates coupled Rac/Rho signalling, MCA (Moesin-like cortical actin) dynamics,
# cortical flow, and membrane deformation on a sphere.
#
# Written by Andreu F Gallen (Turlier lab) and Eric Neiva,
# in collaboration with Orion Weiner's lab.
#
# Reference: INSERT DOI
#
# To reproduce pure local-inhibition behaviour, set:
#   coup.vCTE = 0,  coup.α = 0,  coup.β = 0,  σₐ⁰ = 0

include("Plots_RhoRacSinglet.jl")

# ── Parameter structs ────────────────────────────────────────────────────────

@kwdef struct MechanicalParams
    η::Real           # Viscosity
    χ::Real           # Friction coefficient
    χ₀::Real          # Basal friction
    sigmaₐ⁰::Real     # Basal active tension
    sigmaρ⁰::Real     # Rho-dependent active tension coefficient
    S::Real           # Prestress / geometric parameter
    Λ::Real           # Lamé-like parameter
    M::Real           # Lamé-like parameter
    R::Real           # Cell radius [μm]
end

@kwdef struct KineticParams
    koff::Real; kon::Real; M0::Real; D::Real
    dᵃ::Real; dᵇ::Real; λᵇ::Real; λʳᴬ::Real
    Drac::Real; Drho::Real; α₀::Real; β₀::Real; wrac::Real
end

@kwdef struct CouplingParams
    rac0::Real; rho0::Real; ten0::Real; vCTE::Real
    tenth::Real; sig0::Real; MCAbth::Real
    α::Real; β::Real; αopto::Real; βopto::Real
end

@kwdef struct SimControl
    domain::Tuple
    ls::Any           # Level-set object
    n::Int            # Mesh partition
    Δt::Real
    T::Real
    order::Int
    output_frequency::Int
    γᶜ::Real
    τᵈkₒ::Real
end


# ── Helper functions ─────────────────────────────────────────────────────────

function plotting(ylabel_str, data, path, label)
    plot(data)
    xlabel!("ξ [μm]")
    ylabel!(ylabel_str)
    savefig(path * "$ylabel_str" * label * ".png")
end

"""Mass conservation residual for MCA protein."""
function conservation(sMCAu, sMCAb, Minitial)
    return 1.0 * (Minitial - (sMCAu + sMCAb))
end

threshold(x, x₀, xth)  = 0.5 * (tanh ∘ (x / x₀ - xth / x₀) + 1)
threshold2(x, x₀, xth) = 0.5 * (tanh.(x / x₀ .- xth / x₀) .+ 1)

sinθ(x)  = x[2] / norm(x)
cotθ(x)  = x[1] / x[2]
cotθ2(x) = cotθ(x) * cotθ(x)
cotθ3(x) = cotθ(x) * cotθ(x) * cotθ(x)
cscθ2(x) = 1 + cotθ2(x)


# ── Main simulation function ─────────────────────────────────────────────────

"""
Run a single mechanochemical axisymmetric simulation.

# Arguments
- `mech`     : Mechanical parameters
- `kin`      : Kinetic/chemical parameters
- `coup`     : Mechanochemical coupling and optogenetics parameters
- `control`  : Simulation control (mesh, time stepping, output)
- `name`     : Output folder name
- `writesol` : Write VTU output if true
"""
function run_mechanochemical_axisymmetric_vector(
    mech::MechanicalParams,
    kin::KineticParams,
    coup::CouplingParams,
    control::SimControl;
    name::String,
    writesol::Bool = true
)


    # ── MCA weak forms ───────────────────────────────────────────────────────────

    """
    Weak forms for MCA bound and unbound species.
    Returns (mMCA, aMCAb, bMCAb, aMCAu, bMCAu).
    """
    function MCA_bound_unbound_weak_forms(Δt, kon, koff, λᵇ, λ, R2, D, nΓ, dΓ)

        τ = TensorValue(0.0, -1.0, 1.0, 0.0) ⋅ nΓ   # unit tangent vector

        # Bound MCA
        mMCA(Δt, MCA_b, w) = ∫(((MCA_b * w) / Δt) * y)dΓ

        aMCAb(MCA_b, v, w) =
            ∫(((MCA_b * w) / Δt) * y)dΓ +
            ∫((0.0005 * (∇ᵈ(MCA_b, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +   # diffusion for stabilisation
            ∫((koff * (MCA_b * w)) * y)dΓ +
            ∫((w * ((v ⋅ τ) * (∇ᵈ(MCA_b, nΓ) ⋅ τ))) * y)dΓ +
            ∫(w * (MCA_b * ((τ ⋅ ∇ᵈ(v, nΓ)) ⋅ τ)) * y)dΓ

        bMCAb(w, MCA_u, MCAb_old, λ) =
            ∫((kon * (MCA_u * w) + (MCAb_old * w) / Δt) * y)dΓ +
            ∫((λ / (π * R2) * (kon / (kon + koff)) * w) * y)dΓ

        # Unbound MCA
        aMCAu(MCA_u, x, x_old, w) =
            ∫(((MCA_u * w) / Δt) * y)dΓ +
            ∫((D * (∇ᵈ(MCA_u, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
            ∫((kon * (MCA_u * w)) * y)dΓ +
            ∫((w * (((x - x_old) ⋅ τ / Δt) * (∇ᵈ(MCA_u, nΓ) ⋅ τ) +
                MCA_u * ((τ ⋅ ∇ᵈ(x, nΓ) ⋅ τ) - (τ ⋅ ∇ᵈ(x_old, nΓ) ⋅ τ)) / Δt)) * y)dΓ

        bMCAu(w, MCA_b, MCAu_old, λ) =
            ∫((koff * (MCA_b * w)) * y)dΓ +
            ∫(((MCAu_old * w) / Δt) * y)dΓ +
            ∫((λ / (π * R2) * (koff / (kon + koff)) * w) * y)dΓ

        mMCA, aMCAb, bMCAb, aMCAu, bMCAu
    end

    activity::Function = unit_activity_axisymmetric
    redistance_frequency::Int = 1

    # Background mesh
    cells    = (control.n, div(control.n, 2))
    h        = (control.domain[2] - control.domain[1]) / control.n
    bgmodel  = CartesianDiscreteModel(control.domain, cells)
    Ω        = Triangulation(bgmodel)
    degree   = control.order < 3 ? 3 : 2 * control.order
    R2       = mech.R

    # Level-set buffer (updated each time step)
    buffer = Ref{Any}((Ωᶜ=nothing, dΩᶜ=nothing, dΓ=nothing, nΓ=nothing,
                       φ₋=nothing, cp₋=nothing, t=nothing, Vbg=nothing))

    function update_buffer!(i, t, dt, v₋₂, mv₋₂)
        buffer[].t == t && return true

        Ωᶜ = buffer[].Ωᶜ; Vbg = buffer[].Vbg

        if buffer[].Ωᶜ === nothing
            Vbg  = TestFESpace(Ω, ReferenceFE(lagrangian, Float64, control.order))
            _φ₋  = interpolate_everywhere(control.ls.φ, Vbg)
        else
            cp₋₂ = buffer[].cp₋; φ₋₂ = buffer[].φ₋
            __φ  = get_free_dof_values(φ₋₂.φ)
            Ωⱽ   = get_triangulation(Vbg)
            _ϕ₋  = compute_normal_displacement(cp₋₂, φ₋₂, v₋₂, dt, Ωⱽ)
            _φ₋  = FEFunction(Vbg, __φ - _ϕ₋)
        end

        φ₋ = AlgoimCallLevelSetFunction(_φ₋, ∇(_φ₋))
        (i % redistance_frequency == 0) && begin
            _φ₋ = compute_distance_fe_function(bgmodel, Vbg, φ₋, control.order, cppdegree=3)
            φ₋  = AlgoimCallLevelSetFunction(_φ₋, ∇(_φ₋))
        end

        cp₋   = compute_closest_point_projections(Vbg, φ₋, control.order,
                    cppdegree=3, trim=true, limitstol=1.0e-2)
        squad = Quadrature(algoim, φ₋, degree, phase=CUT)
        s_cell_quad, is_c₋ = CellQuadratureAndActiveMask(bgmodel, squad)

        δ₋     = 2.0 * mv₋₂ * dt
        _, is_nᶜ = narrow_band_triangulation(Ω, _φ₋, Vbg, is_c₋, δ₋)
        Ωᶜ, dΓ = TriangulationAndMeasure(Ω, s_cell_quad, is_nᶜ, is_c₋)
        dΩᶜ    = Measure(Ωᶜ, 2 * control.order)
        nΓ     = normal(φ₋, Ω)

        buffer[] = (Ωᶜ=Ωᶜ, dΩᶜ=dΩᶜ, dΓ=dΓ, nΓ=nΓ, cp₋=cp₋, φ₋=φ₋, t=t, Vbg=Vbg)
        return true
    end

    N      = num_dims(bgmodel)
    reffeʷ = ReferenceFE(lagrangian, VectorValue{N,Float64}, control.order - 1)
    reffeᵉ = ReferenceFE(lagrangian, Float64, control.order - 1)

    function update_all!(i, t, dt, disp, val)
        Ωᶜ  = buffer[].Ωᶜ; dΩᶜ = buffer[].dΩᶜ
        dΓ  = buffer[].dΓ;  nΓ  = buffer[].nΓ
        φ   = buffer[].φ₋

        τ = TensorValue(0.0,-1.0, 1.0, 0.0) ⋅ nΓ    # vector tangente 
        Vʷ  = TestFESpace(Ωᶜ, reffeʷ, dirichlet_tags=[5, 8])
        UXʷ = TrialFESpace(Vʷ, [p -> VectorValue(0, x₀),  p -> -xₗ * τ(p)])
        UVʷ = TrialFESpace(Vʷ, [p -> VectorValue(0, v₀),  p -> -vₗ * τ(p)])

        Vᵉ = TestFESpace(Ωᶜ, reffeᵉ)
        Vᴿ = TestFESpace(Ωᶜ, reffeᵉ)
        Vˡ = ConstantFESpace(bgmodel)

        Uʷ = TrialFESpace(Vʷ)
        Uᵉ = TrialFESpace(Vᵉ)
        Uᴿ = TrialFESpace(Vᴿ)
        Uˡ = TrialFESpace(Vˡ)

        Yᵛ = Vʷ
        Xᵛ = MultiFieldFESpace([Uʷ, Uˡ])
        Yʳ = MultiFieldFESpace([Vᵉ, Vˡ])
        Xʳ = MultiFieldFESpace([Uᵉ, Uˡ])

        UXʷ, UVʷ, Vʷ, Xᵛ, Yᵛ, Xʳ, Yʳ, Uᵉ, Vᵉ, Vᴿ, Uᴿ, Ωᶜ, dΩᶜ, dΓ, nΓ, φ
    end

    # Create output directories only when output is requested. This keeps smoke
    # tests and timing benchmarks from spending most of their time rendering.
    pVTU = "./output/" * name * "VTU/"
    pPNG = "./output/" * name
    writesol && mkpath(pVTU)
    writesol && mkpath(pPNG)
    writesol && mkpath(pPNG * "Rac_time/")
    writesol && mkpath(pPNG * "Rho_time/")
    writesol && mkpath(pPNG * "MCAb_time/")
    writesol && mkpath(pPNG * "Rac_time_initial/")
    writesol && mkpath(pPNG * "Rho_time_initial/")

    # Save source files alongside results for reproducibility
    writesol && cp(@__FILE__, pPNG * split(@__FILE__, "/")[end], force=true)
    writesol && cp("./src/WeakForms.jl", pPNG * "WeakForms.jl", force=true)
    writesol && cp("./examples/SurfaceViscousFlows.jl",
       pPNG * "SurfaceViscousFlows.jl", force=true)

    # Time discretisation
    t₀  = 0.0
    u₀  = VectorValue(0.0, 0.0)
    m₀  = 2.0
    nΔt = trunc(Int, control.T / control.Δt + 0.5) + 1
    tol = 1e-8

    # Boundary condition values (updated during time loop)
    x₀ = 0.0; xₗ = 0.0
    v₀ = 0.0; vₗ = 0.0

    update_buffer!(0, t₀, control.Δt, u₀, m₀)
    UXʷ, UVʷ, Vʷ, Xᵛ, Yᵛ, Xʳ, Yʳ, Uᵉ, Vᵉ, Vᴿ, Uᴿ, Ωᶜ, dΩᶜ, dΓ, nΓ, φ =
        update_all!(0, t₀, control.Δt, u₀, m₀)

    τ = TensorValue(0.0, -1.0, 1.0, 0.0) ⋅ nΓ

    γʷ = control.γᶜ / h   # velocity stabilisation parameter
    γᵉ = control.γᶜ / h   # concentration stabilisation parameter

    # Initial conditions
    _υₕ(x) = VectorValue(0.0, 0.0)
    υₕ     = interpolate_everywhere(_υₕ, UVʷ)
    _xₕ(x) = VectorValue(0.0, 0.0)
    xₕ     = interpolate_everywhere(_xₕ, UXʷ)
    xₕ_old = xₕ

    Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
    Tv = Vector{PetscScalar}
    ps = PETScLinearSolver(mykspsetup)

    i = 0
    t = t₀

    # Arc-length function and optogenetic activation profiles
    arclength(x) = R2 * atan(x[2], -x[1])

    α₀opto(x)  = kin.α₀
    β₀opto(x)  = kin.β₀
    α₀opto2(x) = kin.α₀ + coup.αopto * exp(-0.5 * (arclength(x) - π * R2)^2 / (kin.wrac)^2)
    β₀opto2(x) = kin.β₀ + coup.βopto * exp(-0.5 * (arclength(x))^2           / (kin.wrac)^2)

    γ₀  = 0.1 / h
    γ₀R = 0.1 / h
    γ₀M = 0.1 / h

    m₀opto(u, v)  = ∫(u * v)dΓ
    s₀opto(u, v)  = ∫(γ₀  * ((nΓ ⋅ ∇(u)) ⊙ (nΓ ⋅ ∇(v))))dΩᶜ
    s₀R(u, v)     = ∫(γ₀R * ((nΓ ⋅ ∇(u)) ⊙ (nΓ ⋅ ∇(v))))dΩᶜ
    s₀MCA(u, v)   = ∫(γ₀M * ((nΓ ⋅ ∇(u)) ⊙ (nΓ ⋅ ∇(v))))dΩᶜ
    s₀x(υ, μ)     = ∫(γʷ  * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

    λ = 0.0

    A₀opto(u, v)   = m₀opto(u, v) + s₀opto(u, v)
    bα₀opto(v)     = m₀opto(α₀opto,  v)
    bβ₀opto(v)     = m₀opto(β₀opto,  v)
    bα₀opto2(v)    = m₀opto(α₀opto2, v)
    bβ₀opto2(v)    = m₀opto(β₀opto2, v)

    op_α₀ = AffineFEOperator(A₀opto, bα₀opto, Uᴿ, Vᴿ)
    op_β₀ = AffineFEOperator(A₀opto, bβ₀opto, Uᴿ, Vᴿ)
    α₀v   = solve(op_α₀)
    β₀v   = solve(op_β₀)

    # Initial Rac and Rho (coarse equilibration run with larger diffusion)
    Rₕ     = interpolate_everywhere(0.0, Uᴿ)
    ρₕ     = interpolate_everywhere(1.0, Uᴿ)
    Rₕ_old = Rₕ;  ρₕ_old = ρₕ

    a_R, b_R, a_ρ, b_ρ = rac_rho_weak_forms2(
        control.Δt, 200 * kin.dᵃ, 200 * kin.dᵇ, kin.Drac, kin.Drho,
        nΓ, dΓ, coup.α, coup.β)

    # MCA initial conditions (equilibrium bound/unbound split)
    mMCA, aMCAb, bMCAb, aMCAu, bMCAu = MCA_bound_unbound_weak_forms(
        control.Δt, kin.kon, kin.koff, kin.λᵇ, λ, R2, kin.D, nΓ, dΓ)

    uh_MCAb     = interpolate_everywhere(kin.kon  * kin.M0 / (π * R2) / (kin.koff + kin.kon), Uᴿ)
    uh_MCAb_old = uh_MCAb
    uh_MCAu     = interpolate_everywhere(kin.koff * kin.M0 / (π * R2) / (kin.koff + kin.kon), Uᴿ)
    uh_MCAu_old = uh_MCAu
    sum_uh_MCAu = ∑(∫(uh_MCAu)dΓ)
    sum_uh_MCAb = ∑(∫(uh_MCAb)dΓ)
    Minitial    = sum_uh_MCAu + sum_uh_MCAb
    println("Initial total MCA: ", Minitial)
    λ           = conservation(sum_uh_MCAu, sum_uh_MCAb, Minitial)

    I = TensorValue(1.0, 0.0, 0.0, 1.0)

    # Membrane tension bilinear / linear forms
    N(u)      = mech.S * I + mech.Λ * tr(εᶜ(u, nΓ)) * I + mech.M * εᶜ(u, nΓ)
    ∂u(u)     = mech.R * (τ ⋅ (∇ᶜ(u, nΓ)) ⋅ τ)
    mten(u, v)  = ∫((u * v) * y)dΓ
    mten2(u, v) = ∫(((τ ⋅ (N(u) ⋅ τ)) * v) * y)dΓ
    sten(u, v)  = ∫(10 * γ₀ * ((nΓ ⋅ ∇(u)) ⊙ (nΓ ⋅ ∇(v))))dΩᶜ

    Aten(u, v) = mten(u, v) + sten(u, v)
    bten(v)    = mten2(xₕ, v)
    op_ten     = AffineFEOperator(Aten, bten, Uᴿ, Vᴿ)
    ten        = solve(op_ten)

    # Membrane displacement (nonlinear strain bilinear forms)
    aᴹ(M, R, x_old, x, w) =
        ∫((2 * M * (x ⋅ w / 2 +
                    R * R * (∇ᶜ(x, nΓ) ⊙ ∇ᶜ(w, nΓ)) +
                    (cotθ2) * (x ⋅ w))) * sinθ)dΓ +
        ∫((2 * M / R / 2 *
            ((x_old ⋅ x + R * R * (∇ᶜ(x_old, nΓ) ⊙ ∇ᶜ(x, nΓ))) * (R * ∇ᶜ(w, nΓ) ⋅ τ) +
             (cotθ3 * x_old) * (x ⋅ w))) ⋅ τ * sinθ)dΓ

    aᴸ(L, R, x_old, x, w) =
        ∫(L * (R * R * (∇ᶜ(x, nΓ) ⊙ ∇ᶜ(w, nΓ)) +
               (cotθ * (x ⋅ ∇ᶜ(w, nΓ)) + cotθ * (w ⋅ ∇ᶜ(x, nΓ))) ⋅ τ +
               (cotθ2) * (x ⋅ w)) * sinθ)dΓ +
        ∫(0.5 / R * L *
            ((cscθ2) * (x_old ⋅ x) + R * R * ∇ᶜ(x, nΓ) ⊙ ∇ᶜ(x_old, nΓ)) *
            (R * ∇ᶜ(w, nΓ) ⋅ τ + cotθ * w) ⋅ τ * sinθ)dΓ

    bₓ(MCA_b, v, w)     = ∫((mech.χ * (MCA_b) * v ⋅ w) * (mech.R * mech.R) * sinθ)dΓ
    m(MCA_b, Δt, x, w)  = ∫(((mech.χ * MCA_b) * (x ⋅ w) / Δt) * (mech.R * mech.R) * sinθ)dΓ

    mυ(Δt, v, w) = ∫(((v ⋅ w) / Δt) * y)dΓ

    # Weak tangentiality penalty
    bo = 10.0 / ((2 / 40)^2)
    wt(x, w) = ∫(bo * ((x ⋅ nΓ) * (w ⋅ nΓ)))dΓ

    # Equilibrate Rac/Rho at t = 0
    Arac(rac, w) = a_R(rac, w, υₕ) + s₀R(rac, w)
    Brac(w)      = b_R(w, ρₕ, α₀v, Rₕ_old, uh_MCAb, coup.MCAbth, coup.rho0)
    op_rac       = AffineFEOperator(Arac, Brac, Uᴿ, Vᴿ)
    Rₕ           = solve(op_rac);  Rₕ_old = Rₕ

    Arho(rho, w) = a_ρ(rho, w, υₕ) + s₀R(rho, w)
    Brho(w)      = b_ρ(w, Rₕ, β₀v, ρₕ_old, ten, coup.sig0, coup.tenth)
    op_rho       = AffineFEOperator(Arho, Brho, Uᴿ, Vᴿ)
    ρₕ           = solve(op_rho);  ρₕ_old = ρₕ

    # MCA operators at t = 0
    AMCAb(MCA_b, w) = aMCAb(MCA_b, υₕ, w) + s₀MCA(MCA_b, w)
    BMCAb(w)        = bMCAb(w, uh_MCAu, uh_MCAb_old, λ)
    AMCAu(MCA_u, w) = aMCAu(MCA_u, xₕ, xₕ_old, w) + s₀MCA(MCA_u, w)
    BMCAu(w)        = bMCAu(w, uh_MCAb, uh_MCAu_old, λ)

    # Extract quadrature points ordered by arc length
    xΓ         = dΓ.quad.cell_point.values
    xΓ         = lazy_map(Reindex(xΓ), dΓ.quad.cell_point.ptrs)
    alenΓ      = lazy_map(Broadcasting(x -> atan(x[2], -x[1])), xΓ)
    flat_xΓ    = vcat(xΓ...)
    flat_alenΓ = vcat(alenΓ...)
    num_qpoints = length(flat_xΓ)
    perm        = sortperm(flat_alenΓ)
    flat_alenΓ  = R2 * flat_alenΓ[perm]

    # Pre-allocate time-series arrays
    ract  = zeros(nΔt, num_qpoints)
    rhot  = zeros(nΔt, num_qpoints)
    MCAbt = zeros(nΔt, num_qpoints)
    σₐt   = zeros(nΔt, num_qpoints)
    χt    = zeros(nΔt, num_qpoints)
    vt    = zeros(nΔt, num_qpoints)
    vnt   = zeros(nΔt, num_qpoints)
    xt    = zeros(nΔt, num_qpoints)
    tent  = zeros(nΔt, num_qpoints)

    # Record initial state
    vt[1,:]    = vcat(lazy_map(υₕ ⋅ τ,  xΓ)...)[perm]
    vnt[1,:]   = vcat(lazy_map(υₕ ⋅ nΓ, xΓ)...)[perm]
    xt[1,:]    = vcat(lazy_map(xₕ ⋅ τ,  xΓ)...)[perm]
    MCAbt[1,:] = vcat(lazy_map(uh_MCAb,  xΓ)...)[perm]
    tent[1,:]  = vcat(lazy_map(ten,       xΓ)...)[perm]

    # Pre-equilibration loop (coarse diffusion coefficients)
    equilibration_steps = parse(Int, get(ENV, "DEBELLY_EQUILIBRATION_STEPS", "100"))
    for ti in 1:equilibration_steps
        op_rho = AffineFEOperator(Arho, Brho, Uᴿ, Vᴿ)
        ρₕ     = solve(op_rho);  ρₕ_old = ρₕ
        op_rac = AffineFEOperator(Arac, Brac, Uᴿ, Vᴿ)
        Rₕ     = solve(op_rac);  Rₕ_old = Rₕ

        ractt = vcat(lazy_map(Rₕ, xΓ)...)[perm]
        rhott = vcat(lazy_map(ρₕ, xΓ)...)[perm]
        writesol && plotting("rac", ractt, pPNG * "Rac_time_initial/", "$ti")
        writesol && plotting("rho", rhott, pPNG * "Rho_time_initial/", "$ti")
    end

    # Switch to physical diffusion coefficients
    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)
    a_R, b_R, a_ρ, b_ρ = rac_rho_weak_forms2(
        control.Δt, kin.dᵃ, kin.dᵇ, kin.Drac, kin.Drho,
        nΓ, dΓ, coup.α, coup.β)
    #Arac(rac, w) = a_R(rac, w, υₕ) + s₀R(rac, w)

    # ── Main time loop ───────────────────────────────────────────────────────

    while t < control.T + tol

        # Optogenetic activation window: on after step 50, off after step 150
        if i == 50
            op_α₀ = AffineFEOperator(A₀opto, bα₀opto2, Uᴿ, Vᴿ)
            op_β₀ = AffineFEOperator(A₀opto, bβ₀opto2, Uᴿ, Vᴿ)
            α₀v = solve(op_α₀);  β₀v = solve(op_β₀)
        end
        if i == 150
            op_α₀ = AffineFEOperator(A₀opto, bα₀opto, Uᴿ, Vᴿ)
            op_β₀ = AffineFEOperator(A₀opto, bβ₀opto, Uᴿ, Vᴿ)
            α₀v = solve(op_α₀);  β₀v = solve(op_β₀)
        end

        # Update mass conservation Lagrange multiplier
        i1 = ∑(∫(uh_MCAb)dΓ);  i2 = ∑(∫(uh_MCAu)dΓ)
        λ  = conservation(i2, i1, Minitial)

        @info "Time step $i, time $(trunc(t, digits=4)), Δt = $(control.Δt)"

        # Solve velocity
        aᵛ, bᵛ = cortical_flow_problem_mechanochemical_axisymmetric_dimensional(
            xₕ, xₕ_old, control.Δt, mech.η,
            ρₕ, uh_MCAb, dΩᶜ, dΓ, nΓ, γʷ,
            mech.χ, mech.χ₀, activity, mech.sigmaₐ⁰, mech.sigmaρ⁰)
        op  = AffineFEOperator(aᵛ, bᵛ, UVʷ, Yᵛ)
        υₕ  = solve(op)
        υₕtan = to_tangent_vector(υₕ, nΓ)

        # Solve membrane displacement
        aˣ(x, w) = m(uh_MCAb, control.Δt, x, w) +
                   aᴸ(mech.Λ, R2, xₕ, x, w) +
                   aᴹ(mech.M, R2, xₕ, x, w) +
                   s₀x(x, w) + wt(x, w)
        bˣ(w)    = m(uh_MCAb, control.Δt, xₕ, w) + bₓ(uh_MCAb, υₕ, w)
        op_x     = AffineFEOperator(aˣ, bˣ, UXʷ, Vʷ)
        xₕ       = solve(op_x)

        # Update tension
        op_ten = AffineFEOperator(Aten, bten, Uᴿ, Vᴿ)
        ten    = solve(op_ten)
        xₕ_old = xₕ

        i = i + 1
        t = t + control.Δt

        writesol && postprocess_all_with_tangent(
            φ, dΩᶜ.quad.trian, Rₕ, ρₕ, xₕ, υₕ, υₕtan, uh_MCAb, ten;
            i=i, of=control.output_frequency, name=pVTU)

        # Update leading-edge boundary condition (velocity and displacement)
        Rₕaux = vcat(lazy_map(Rₕ, xΓ)...)[perm]
        _ten  = vcat(lazy_map(ten, xΓ)...)[perm]
        _ten  = _ten[end] * threshold2(_ten[end], 0.0001, 0.0)
        vₗ    = coup.vCTE *
                threshold2(Rₕaux[end], coup.rac0, 1.3 * Rₕaux[1]) /
                (1 + _ten * _ten / coup.ten0)
        xₗ   = xₗ * 0.9 + control.Δt * vₗ   # slight spatial relaxation at boundary

        # Rebuild FE spaces on updated geometry
        UXʷ, UVʷ, Vʷ, Xᵛ, Yᵛ, Xʳ, Yʳ, Uᵉ, Vᵉ, Vᴿ, Uᴿ, Ωᶜ, dΩᶜ, dΓ, nΓ, φ =
            update_all!(i, t, control.Δt, υₕ, msₕ)

        # Solve Rac and Rho
        op_rac = AffineFEOperator(Arac, Brac, Uᴿ, Vᴿ)
        Rₕ     = solve(op_rac);  Rₕ_old = Rₕ
        op_rho = AffineFEOperator(Arho, Brho, Uᴿ, Vᴿ)
        ρₕ     = solve(op_rho);  ρₕ_old = ρₕ

        # Solve MCA bound and unbound
        op_MCAb  = AffineFEOperator(AMCAb, BMCAb, Uᴿ, Vᴿ)
        uh_MCAb  = solve(op_MCAb);  uh_MCAb_old = uh_MCAb
        op_MCAu  = AffineFEOperator(AMCAu, BMCAu, Uᴿ, Vᴿ)
        uh_MCAu  = solve(op_MCAu);  uh_MCAu_old = uh_MCAu

        # Record time-series data
        ract[i,:]  = vcat(lazy_map(Rₕ,      xΓ)...)[perm]
        rhot[i,:]  = vcat(lazy_map(ρₕ,      xΓ)...)[perm]
        MCAbt[i,:] = vcat(lazy_map(uh_MCAb,  xΓ)...)[perm]
        σₐt[i,:]   = mech.sigmaₐ⁰ .+ mech.sigmaρ⁰ * rhot[i,:]
        χt[i,:]    = mech.χ₀ .+ mech.χ * MCAbt[i,:]
        vt[i,:]    = vcat(lazy_map(υₕ ⋅ τ,  xΓ)...)[perm]
        vnt[i,:]   = vcat(lazy_map(υₕ ⋅ nΓ, xΓ)...)[perm]
        xt[i,:]    = vcat(lazy_map(xₕ ⋅ τ,  xΓ)...)[perm]
        tent[i,:]  = vcat(lazy_map(ten,       xΓ)...)[perm]

        writesol && plotting("rac",  ract[i,:],  pPNG * "Rac_time/",  "$i")
        writesol && plotting("rho",  rhot[i,:],  pPNG * "Rho_time/",  "$i")
        writesol && plotting("MCAb", MCAbt[i,:], pPNG * "MCAb_time/", "$i")
    end

    writesol && plots_run_singlet(nΔt, vt, vnt, xt * mech.R, MCAbt, ract, rhot, pPNG,
        num_qpoints, π * R2, control.Δt, control.T, flat_alenΓ, σₐt, χt, tent)
    return (; arclength=flat_alenΓ, rac=ract, rho=rhot, mcab=MCAbt,
        tension=tent, tangential_velocity=vt, normal_velocity=vnt,
        displacement=xt)
end
