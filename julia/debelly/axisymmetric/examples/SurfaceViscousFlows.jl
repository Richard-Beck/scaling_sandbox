# Entry point for a mechanochemical simulation on an axisymmetric cell surface.
# Loads the SurfaceBulkViscousFlows module, sets physical parameters, and runs
# one simulation. Output is written to ./output/<name>/.
#
# See MechanochemicalAxisymmetricVector.jl for the simulation driver.

using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC
using SurfaceBulkViscousFlows

# ── Geometry ─────────────────────────────────────────────────────────────────

R      = 6.0   # Cell radius [μm]
domain = (-1.2 * R, 0.9 * R, 0.0, 1.2 * R)

ls = AlgoimCallLevelSetFunction(
    x -> x[1]^2 + x[2]^2 - R^2,
    x -> VectorValue(2.0 * x[1], 2.0 * x[2]))

# ── Physical parameters ───────────────────────────────────────────────────────

ξ = 1.00
k = 40.0

mech = SurfaceBulkViscousFlows.MechanicalParams(
    η        = 10000.0,
    χ        = 2.0,
    χ₀       = 0.0,
    sigmaₐ⁰  = 100.0,
    sigmaρ⁰  = 200.0,
    S        = trunc(sqrt(3) * k * (1 - 1 / ξ),    digits=2),
    Λ        = trunc(sqrt(3) * 0.25 * k / ξ,        digits=2),
    M        = trunc(sqrt(3) * 0.25 * k / ξ,        digits=2),
    R        = R
)

kin = SurfaceBulkViscousFlows.KineticParams(
    koff  = 1.4,  kon  = 5.0,  M0   = 1.0,  D    = 0.3,
    dᵃ    = 0.04, dᵇ   = 0.04, λᵇ   = 0.0,  λʳᴬ  = 0.0,
    Drac  = 0.05, Drho = 0.05, α₀   = 1.0,  β₀   = 1.0,
    wrac  = R
)

coup = SurfaceBulkViscousFlows.CouplingParams(
    rac0   = 0.01,  rho0   = 0.001, ten0   = 0.02,
    vCTE   = -0.3,
    tenth  = 0.025, sig0   = 0.01,  MCAbth = 0.035,
    α      = 1.0,   β      = 1.0,
    αopto  = 2.0,   βopto  = 0.0
)

control = SurfaceBulkViscousFlows.SimControl(
    domain           = domain,
    ls               = ls,
    n                = parse(Int, get(ENV, "DEBELLY_AXIS_N", "40")),
    Δt               = parse(Float64, get(ENV, "DEBELLY_DT", "1")),
    T                = parse(Float64, get(ENV, "DEBELLY_T", "350")),
    order            = 2,
    output_frequency = 1,
    γᶜ               = 1.0,
    τᵈkₒ             = 10.0
)

# ── Derived quantities ────────────────────────────────────────────────────────

len = trunc(sqrt(mech.η * R * π / (mech.χ * R^2 * kin.M0 * π^2)), digits=2)
println("Hydrodynamic length scale = $len μm")

# ── Output path ───────────────────────────────────────────────────────────────

# name="SurfaceViscousFlows/velocity_on/rho-rac 5 better friction T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
#name="SurfaceViscousFlows/conserved/test lowrac T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ rac_t=$rac_total rho_t=$rho_total  x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
name = "T=$(control.T)_dt=$(control.Δt)_ n=$(control.n)/D=$(kin.D)_" *
       "eta=$(mech.η)_vC=$(coup.vCTE) S=$(mech.S)_M=$(mech.M)_a0=$(kin.α₀) b0=$(kin.β₀) wrac=$(kin.wrac) Drac=$(kin.Drac)" *
       "/rac0=$(coup.rac0) koff=$(kin.koff) aopto = $(coup.αopto) da=$(kin.dᵃ) db=$(kin.dᵇ)_sig=$(mech.sigmaρ⁰)_x=$(mech.χ)  L=$len a=$(coup.α) b=$(coup.β) sig0=$(coup.sig0)_tenth=$(coup.tenth)_MCAbth=$(coup.MCAbth)_rho0=$(coup.rho0) ten0=$(coup.ten0)/" 

#mkpath(name)

# ── Run ───────────────────────────────────────────────────────────────────────

elapsed = @elapsed result = GridapPETSc.with() do
    run_mechanochemical_axisymmetric_vector(mech, kin, coup, control;
        name=name, writesol=get(ENV, "DEBELLY_WRITE_OUTPUTS", "1") == "1")
end
rac_final = result.rac[end, :]
rho_final = result.rho[end, :]
println("DEBELLY_SUMMARY protocol=axisymmetric n=$(control.n) steps=$(size(result.rac, 1)-1) elapsed_seconds=$elapsed rac_range=$(maximum(rac_final)-minimum(rac_final)) rho_range=$(maximum(rho_final)-minimum(rho_final))")
