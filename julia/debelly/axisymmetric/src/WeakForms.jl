# Axisymmetric weight y = x[2] (the radial coordinate r in cylindrical symmetry).
iy(x) = VectorValue(0.0, 1.0 / x[2])
y(x)  = x[2]


# ============================================================
#   Rac / Rho reaction-diffusion weak forms
# ============================================================

"""
Basic Rac–Rho system without mechanical coupling.
Returns bilinear forms (a_rac, a_rho) and linear forms (b_rac, b_rho).
"""
function rac_rho_weak_forms2(Δt, dᵃ, dᵇ, Drac, Drho, nΓ, dΓ) 

  m2(Δt, A, w) = ∫(((A * w) / Δt) * y)dΓ
  advection(rac, w, v) = ∫((w * (∇ᵈ(rac, nΓ) ⋅ v + rac * divᶜ(v, nΓ))) * y)dΓ

  a_rac(rac, w, v) =
    (1 / dᵃ) * m2(Δt, rac, w) +
    ∫((Drac * (∇ᵈ(rac, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    ∫((w * rac) * y)dΓ

  b_rac(w, rho, α₀v, rac_old) =
    (1 / dᵃ) * m2(Δt, rac_old, w) +
    ∫((w * (α₀v) / (1 + rho * rho)) * y)dΓ

  a_rho(rho, w, v) =
    (1 / dᵇ) * m2(Δt, rho, w) +
    advection(rho, w, v) +
    ∫((Drho * (∇ᵈ(rho, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    ∫((w * rho) * y)dΓ

  b_rho(w, rac, β₀v, rho_old) =
    (1 / dᵇ) * m2(Δt, rho_old, w) +
    ∫((w * (β₀v) / (1 + rac * rac)) * y)dΓ

  a_rac, b_rac, a_rho, b_rho
end


"""
Rac–Rho system with mechanochemical coupling via threshold functions.
`α` and `β` are coupling strengths; MCA_b and tension fields modulate activation.
"""
function rac_rho_weak_forms2(Δt, dᵃ, dᵇ, Drac, Drho, nΓ, dΓ, α, β)

  threshold(x, x₀, xth) = 0.5 * (tanh ∘ (x / x₀ - xth / x₀) + 1)

  m2(Δt, A, w) = ∫(((A * w) / Δt) * y)dΓ
  advection(rac, w, v) = ∫((w * (∇ᵈ(rac, nΓ) ⋅ v + rac * divᶜ(v, nΓ))) * y)dΓ

  a_rac(rac, w, v) =
    (1 / dᵃ) * m2(Δt, rac, w) +
    ∫((Drac * (∇ᵈ(rac, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    ∫((w * rac) * y)dΓ

  b_rac(w, rho, α₀v, rac_old, MCA_b, MCAbth, rho0) =
    (1 / dᵃ) * m2(Δt, rac_old, w) +
    ∫((w * (α₀v + α * (1.0 - threshold(MCA_b, rho0, MCAbth))) / (1 + rho * rho)) * y)dΓ

  a_rho(rho, w, v) =
    (1 / dᵇ) * m2(Δt, rho, w) +
    ∫((Drho * (∇ᵈ(rho, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    ∫((w * rho) * y)dΓ

  b_rho(w, rac, β₀v, rho_old, ten, sig0, tenth) =
    (1 / dᵇ) * m2(Δt, rho_old, w) +
    ∫((w * (β₀v + β * threshold(ten, sig0, tenth)) / (1 + rac * rac)) * y)dΓ

  a_rac, b_rac, a_rho, b_rho
end


"""
Mass-conserved Rac–Rho system.
Includes auxiliary fields (rac_i, rho_i) to enforce global conservation.
"""
function rac_rho_weak_forms_conserved(Δt, dᵃ, dᵇ, Drac, Drho, nΓ, dΓ)

  m2(Δt, A, w) = ∫(((A * w) / Δt) * y)dΓ
  advection(rac, w, v) = ∫((w * (∇ᵈ(rac, nΓ) ⋅ v + rac * divᶜ(v, nΓ))) * y)dΓ

  a_rac(rac, w, v) =
    (1 / dᵃ) * m2(Δt, rac, w) +
    ∫((Drac * (∇ᵈ(rac, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    ∫((w * rac) * y)dΓ

  b_rac(w, rho, α₀v, rac_i, rac_old) =
    (1 / dᵃ) * m2(Δt, rac_old, w) +
    ∫(w * (rac_i) * (α₀v / (1 + rho * rho)) * y)dΓ

  a_rac_i(rac_i, w) = ∫(((w * (rac_i)) * y))dΓ
  b_rac_i(w, a_t, a_sum) = ∫((w * (a_t - a_sum)) * y)dΓ

  a_rho(rho, w, v) =
    (1 / dᵇ) * m2(Δt, rho, w) +
    ∫((Drho * (∇ᵈ(rho, nΓ) ⋅ ∇ᵈ(w, nΓ))) * y)dΓ +
    advection(rho, w, v) +
    ∫((w * rho) * y)dΓ

  b_rho(w, rac, β₀v, rho_i, rho_old) =
    (1 / dᵇ) * m2(Δt, rho_old, w) +
    ∫((w * (rho_i) * (β₀v / (1 + rac * rac))) * y)dΓ

  a_rho_i(rho_i, w) = ∫(((w * (rho_i)) * y))dΓ
  b_rho_i(w, b_t, b_sum) = ∫((w * (b_t - b_sum)) * y)dΓ

  a_rac, b_rac, a_rho, b_rho, a_rac_i, b_rac_i, a_rho_i, b_rho_i
end


# ============================================================
#   Cortical flow problems (axisymmetric)
# ============================================================

"""
Basic cortical flow with myosin-driven active stress.
ξ₀: background activity field; Pe: Péclet-like activity number; χ: friction.
"""
function cortical_flow_problem_axisymmetricOG(
    eₕ, dΩᶜ, dΓ, nΓ, γ::Float64, Pe::Float64, χ::Float64, ξ₀::Function)

  aʷ(u, v) =
    ∫((εᶜ(u, nΓ) ⊙ εᵈ(v, nΓ) + divᶜ(u, nΓ) ⋅ divᶜ(v, nΓ) +
       2 * (u ⋅ iy) * (v ⋅ iy) + divᶜ(u, nΓ) * (v ⋅ iy) +
       divᶜ(v, nΓ) * (u ⋅ iy)) * y)dΓ

  aᶠ(u, v) = ∫(χ * (u ⋅ v) * y)dΓ

  ξ(e) = 2.0 * e * e / (1.0 + e * e)

  f(μ, e) = ∫(Pe * (-(divᶜ(μ, nΓ) + μ ⋅ iy) * (ξ ∘ (e))) * ξ₀)dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((u ⋅ (ℓ * RB¹)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ((υ, l¹, l²), (μ, ℓ¹, ℓ²)) =
    aʷ(υ, μ) + aᶠ(υ, μ) + sᵘ(υ, μ) +
    r¹(υ, ℓ¹) + r¹(μ, l¹) + r²(υ, ℓ²) + r²(μ, l²)

  bᵛ((μ, ℓ¹, ℓ²)) = f(μ, eₕ)

  aᵛ, bᵛ
end


"""
Cortical flow with myosin turnover: viscosity scales with myosin field eₕ,
and active stress depends on Rho (ρₕ) and actin (R) concentrations.
"""
function cortical_flow_problem_axisymmetric_turnover(
    ρₕ, R, eₕ, dΩᶜ, dΓ, nΓ, γ::Float64, Pe::Float64,
    χᵣ::Float64, χ₀::Float64, ξ₀, σₐ⁰,
    sigmaₐ⁰, sigmaρ⁰, sigmaR⁰)

  aʷ(u, v, e) =
    ∫((e * εᶜ(u, nΓ) ⊙ εᵈ(v, nΓ) + e * divᶜ(u, nΓ) ⋅ divᶜ(v, nΓ) +
       2 * e * (u ⋅ iy) * (v ⋅ iy) + e * divᶜ(u, nΓ) * (v ⋅ iy) +
       e * divᶜ(v, nΓ) * (u ⋅ iy)) * y)dΓ

  χ(R) = χ₀ + χᵣ * R
  aᶠ(u, v, R) = ∫(χ(R) * (u ⋅ v) * y)dΓ

  function sigmaₐ(ρ, R)
    s = sigmaₐ⁰ + sigmaρ⁰ * ρ - sigmaR⁰ * R
    s > 0 ? s : zero(typeof(s))
  end

  f(μ, ρ, R, e) =
    ∫(e * (-(divᶜ(μ, nΓ) + μ ⋅ iy) * (sigmaₐ ∘ (ρ, R)) * σₐ⁰) * ξ₀)dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((RB¹ ⋅ (ℓ * u)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ((υ, l¹, l²), (μ, ℓ¹, ℓ²)) =
    aʷ(υ, μ, eₕ) + aᶠ(υ, μ, R) + sᵘ(υ, μ) +
    r¹(υ, ℓ¹) + r¹(μ, l¹) + r²(υ, ℓ²) + r²(μ, l²)

  bᵛ((μ, ℓ¹, ℓ²)) = f(μ, ρₕ, R, eₕ)

  aᵛ, bᵛ
end


"""
Cortical flow with Rho-dependent active stress. No myosin turnover.
"""
function cortical_flow_problem_axisymmetric(
    ρₕ, R, dΩᶜ, dΓ, nΓ, γ::Float64, Pe::Float64,
    χᵣ::Float64, χ₀::Float64, ξ₀, σₐ⁰, sigmaₐ⁰, sigmaρ⁰, sigmaR⁰)

  aʷ(u, v) =
    ∫((εᶜ(u, nΓ) ⊙ εᵈ(v, nΓ) + divᶜ(u, nΓ) ⋅ divᶜ(v, nΓ) +
       2 * (u ⋅ iy) * (v ⋅ iy) + divᶜ(u, nΓ) * (v ⋅ iy) +
       divᶜ(v, nΓ) * (u ⋅ iy)) * y)dΓ

  χ(R) = χ₀ + χᵣ * R
  aᶠ(u, v, R) = ∫(χ(R) * (u ⋅ v) * y)dΓ

  function sigmaₐ(ρ, R)
    s = sigmaₐ⁰ + sigmaρ⁰ * ρ - sigmaR⁰ * R
    s > 0 ? s : zero(typeof(s))
  end

  f(μ, ρ, R) =
    ∫((-(divᶜ(μ, nΓ) + μ ⋅ iy) * (sigmaₐ ∘ (ρ, R)) * σₐ⁰) * ξ₀)dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((RB¹ ⋅ (ℓ * u)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ((υ, l¹, l²), (μ, ℓ¹, ℓ²)) =
    aʷ(υ, μ) + aᶠ(υ, μ, R) + sᵘ(υ, μ) +
    r¹(υ, ℓ¹) + r¹(μ, l¹) + r²(υ, ℓ²) + r²(μ, l²)

  bᵛ((μ, ℓ¹, ℓ²)) = f(μ, ρₕ, R)

  aᵛ, bᵛ
end


"""
Mechanochemical cortical flow (axisymmetric, non-dimensional).
Tangentiality enforced weakly via a penalty term.
"""
function cortical_flow_problem_mechanochemical_axisymmetric(
    ρₕ, ez, dΩᶜ, dΓ, nΓ, γ::Float64, Pe::Float64,
    χᵣ::Float64, χ₀::Float64, ξ₀, sigmaₐ⁰, sigmaρ⁰)

  aʷ(u, v) =
    ∫((εᶜ(u, nΓ) ⊙ εᵈ(v, nΓ) + divᶜ(u, nΓ) ⋅ divᶜ(v, nΓ) +
       2 * (u ⋅ iy) * (v ⋅ iy) + divᶜ(u, nΓ) * (v ⋅ iy) +
       divᶜ(v, nΓ) * (u ⋅ iy)) * y)dΓ

  χ(R) = χ₀ + χᵣ * R
  aᶠ(u, v, R) = ∫(χ(R) * (u ⋅ v) * y)dΓ

  function sigmaₐ(ρ, R)
    s = sigmaₐ⁰ + sigmaρ⁰ * ρ
    s > 0 ? s : zero(typeof(s))
  end

  f(μ, ρ, R) =
    ∫((-(divᶜ(μ, nΓ) + μ ⋅ iy) * (sigmaₐ ∘ (ρ, R))) * ξ₀)dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  # Weak tangentiality penalty
  η_pen = 10.0 / ((2 / 40)^2)
  k(u, v) = ∫(η_pen * ((u ⋅ nΓ) * (v ⋅ nΓ)))dΓ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((RB¹ ⋅ (ℓ * u)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ(υ, μ) = aʷ(υ, μ) + sᵘ(υ, μ) + k(υ, μ) + aᶠ(υ, μ, ez)
  bᵛ(μ)    = f(μ, ρₕ, ez)

  aᵛ, bᵛ
end


"""
Mechanochemical cortical flow (axisymmetric, dimensional).
Includes membrane inertia through a backward-Euler mass term.
"""
function cortical_flow_problem_mechanochemical_axisymmetric_dimensional(
    xₕ, xₕ_old, Δt, η_visc,
    ρₕ, ez, dΩᶜ, dΓ, nΓ, γ::Float64,
    χᵣ::Float64, χ₀::Float64, ξ₀, sigmaₐ⁰, sigmaρ⁰)

  aʷ(u, v) =
    ∫(η_visc * (εᶜ(u, nΓ) ⊙ εᵈ(v, nΓ) + divᶜ(u, nΓ) ⋅ divᶜ(v, nΓ) +
       2 * (u ⋅ iy) * (v ⋅ iy) + divᶜ(u, nΓ) * (v ⋅ iy) +
       divᶜ(v, nΓ) * (u ⋅ iy)) * y)dΓ

  dotx  = (xₕ - xₕ_old) / Δt
  χ(R)  = χ₀ + χᵣ * R
  aᶠ(u, v, R) = ∫(χ(R) * (u ⋅ v) * y)dΓ

  function sigmaₐ(ρ, R)
    s = sigmaₐ⁰ + sigmaρ⁰ * ρ
    s > 0 ? s : zero(typeof(s))
  end

  f(μ, ρ, R) =
    ∫((-(divᶜ(μ, nΓ) + μ ⋅ iy) * (sigmaₐ ∘ (ρ, R))) * ξ₀)dΓ

  sᵘ(υ, μ) = ∫((η_visc * γ) * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  # Weak tangentiality penalty
  η_pen = 10000.0 / ((2 / 40)^2)
  k(u, v) = ∫(η_pen * ((u ⋅ nΓ) * (v ⋅ nΓ)))dΓ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((RB¹ ⋅ (ℓ * u)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ(υ, μ) = aʷ(υ, μ) + sᵘ(υ, μ) + k(υ, μ) + aᶠ(υ, μ, ez)
  bᵛ(μ)    = f(μ, ρₕ, ez) + aᶠ(dotx, μ, ez)

  aᵛ, bᵛ
end


"""
Cortical flow coupled to a bulk Stokes layer (axisymmetric).
"""
function cortical_flow_problem_axisymmetric(
    ulₕ, plₕ, eₕ, dΩᶜ, dΓ, nΓ,
    γ::Float64, Pe::Float64, μˡ::Float64, R::Float64, ξ₀::Function)

  aʷ(υ, μ) = ∫(2.0 * (εᶜ(υ, nΓ) ⊙ εᵈ(μ, nΓ) + (υ ⋅ iy) * (μ ⋅ iy)) * y)dΓ

  ξ(e) = 2.0 * e * e / (1.0 + e * e)

  f(μ, e) = ∫(Pe * (-(divᶜ(μ, nΓ) + μ ⋅ iy) * (ξ ∘ (e))) * ξ₀)dΓ

  σᵘ(ε, q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ, u, p) = ∫((μ ⋅ ((σᵘ ∘ (ε(u), p)) ⋅ nΓ)) * y)dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  RB¹ = VectorValue(1.0, 0.0)
  r¹(u, ℓ) = ∫((u ⋅ (ℓ * RB¹)) * y)dΓ
  r²(u, ℓ) = ∫((u ⋅ (ℓ * nΓ)) * y)dΓ

  aᵛ((υ, l¹, l²), (μ, ℓ¹, ℓ²)) =
    aʷ(υ, μ) + sᵘ(υ, μ) +
    r¹(υ, ℓ¹) + r¹(μ, l¹) + r²(υ, ℓ²) + r²(μ, l²)

  bᵛ((μ, ℓ¹, ℓ²)) = f(μ, eₕ) - βʳ(μ, ulₕ, plₕ)

  aᵛ, bᵛ
end


# ============================================================
#   Bulk flow problems
# ============================================================

"""Axisymmetric Stokes bulk flow with Nitsche coupling to the cortex."""
function bulk_flow_problem_axisymmetric(
    υ, dΩ, dΓ, nΓ, μˡ::Float64, R::Float64, γ::Float64, h::Float64)

  aᵇ(u, v) = ∫(2.0 * μˡ * R * (ε(u) ⊙ ε(v) + (u ⋅ iy) * (v ⋅ iy)) * y)dΩ
  bᵇ(v, q) = ∫(q * (∇ ⋅ v + v ⋅ iy) * y)dΩ

  σᵘ(ε, q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u, v, p, q) =
    ∫(((γ / h) * (u ⋅ v) -
       u ⋅ ((σᵘ ∘ (ε(v), q)) ⋅ nΓ) -
       v ⋅ ((σᵘ ∘ (ε(u), p)) ⋅ nΓ)) * y)dΓ
  αʳ(v, q, υ) =
    ∫(((γ / h) * (υ ⋅ v) - υ ⋅ ((σᵘ ∘ (ε(v), q)) ⋅ nΓ)) * y)dΓ

  r(p, ℓ) = ∫((p * ℓ) * y)dΩ

  aᵘ((uˡ, pˡ, l), (vˡ, qˡ, ℓ)) =
    aᵇ(uˡ, vˡ) - bᵇ(vˡ, pˡ) - bᵇ(uˡ, qˡ) +
    αˡ(uˡ, vˡ, pˡ, qˡ) + r(pˡ, ℓ) + r(qˡ, l)

  bᵘ((vˡ, qˡ, l)) = αʳ(vˡ, qˡ, υ)

  aᵘ, bᵘ
end


"""3D Stokes bulk flow with Nitsche coupling."""
function bulk_flow_problem_3D(
    υ, dΩ, dΓ, nΓ, μˡ::Float64, R::Float64, γ::Float64, h::Float64)

  aᵇ(u, v) = ∫(2.0 * μˡ * R * (ε(u) ⊙ ε(v)))dΩ
  bᵇ(v, q) = ∫(q * (∇ ⋅ v))dΩ

  σᵘ(ε, q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u, v, p, q) =
    ∫((γ / h) * (u ⋅ v) -
      u ⋅ ((σᵘ ∘ (ε(v), q)) ⋅ nΓ) -
      v ⋅ ((σᵘ ∘ (ε(u), p)) ⋅ nΓ))dΓ
  αʳ(v, q, υ) =
    ∫((γ / h) * (υ ⋅ v) - υ ⋅ ((σᵘ ∘ (ε(v), q)) ⋅ nΓ))dΓ

  r(p, ℓ) = ∫(p * ℓ)dΩ

  aᵘ((uˡ, pˡ, l), (vˡ, qˡ, ℓ)) =
    aᵇ(uˡ, vˡ) - bᵇ(vˡ, pˡ) - bᵇ(uˡ, qˡ) +
    αˡ(uˡ, vˡ, pˡ, qˡ) + r(pˡ, ℓ) + r(qˡ, l)

  bᵘ((vˡ, qˡ, l)) = αʳ(vˡ, qˡ, υ)

  aᵘ, bᵘ
end


# ============================================================
#   Transport / turnover problems
# ============================================================

"""Axisymmetric surface transport with turnover (reaction + advection + diffusion)."""
function transport_problem_axisymmetric(
    u, eₕ, dΓ, dΩᶜ, nΓ, dt::Float64, γ::Float64, τᵈkₒ::Float64)

  m(e, ε)  = ∫((1 / dt) * (e * ε) * y)dΓ
  sᵈ(e, ε) = ∫((∇ᵈ(e, nΓ) ⋅ ∇ᵈ(ε, nΓ)) * y)dΓ
  c(e, ε)  = ∫(((u ⋅ ∇ᵈ(e, nΓ)) * ε + (tr(∇ᵈ(u, nΓ)) + u ⋅ iy) * (e * ε)) * y)dΓ
  r(e, ε)  = ∫(τᵈkₒ * (e * ε) * y)dΓ
  l(ε)     = ∫(τᵈkₒ * ε * y)dΓ
  s(υ, μ)  = ∫(γ * ((nΓ ⋅ ∇(υ)) ⊙ (nΓ ⋅ ∇(μ))))dΩᶜ

  aᵉ(e, ε) = m(e, ε) + c(e, ε) + r(e, ε) + sᵈ(e, ε) + s(e, ε)
  bᵉ(ε)    = m(eₕ, ε) + l(ε)

  aᵉ, bᵉ
end


"""Axisymmetric myosin turnover with diffusion but without advection."""
function turnover_axisymmetric(
    u, eₕ, dΓ, dΩᶜ, nΓ, De, dt::Float64, γ::Float64, τᵈkₒ::Float64)

  m(e, ε)  = ∫((1 / dt) * (e * ε) * y)dΓ
  sᵈ(e, ε) = ∫((∇ᵈ(e, nΓ) ⋅ ∇ᵈ(ε, nΓ)) * y)dΓ
  c(e, ε)  = ∫(((u ⋅ ∇ᵈ(e, nΓ)) * ε + (tr(∇ᵈ(u, nΓ)) + u ⋅ iy) * (e * ε)) * y)dΓ
  r(e, ε)  = ∫(τᵈkₒ * (e * ε) * y)dΓ
  d(e, ε)  = ∫((De * (∇ᵈ(e, nΓ) ⋅ ∇ᵈ(ε, nΓ))) * y)dΓ
  l(ε)     = ∫(τᵈkₒ * ε * y)dΓ

  aᵉ(e, ε) = m(e, ε) + c(e, ε) + r(e, ε) + d(e, ε) + sᵈ(e, ε)
  bᵉ(ε)    = m(eₕ, ε) + l(ε)

  aᵉ, bᵉ
end


"""3D surface transport with turnover."""
function transport_problem_3D(
    u, eₕ, dΓ, dΩᶜ, nΓ, dt::Float64, γ::Float64, τᵈkₒ::Float64)

  m(e, ε)  = ∫((1 / dt) * (e * ε))dΓ
  sᵈ(e, ε) = ∫(∇ᵈ(e, nΓ) ⋅ ∇ᵈ(ε, nΓ))dΓ
  c(e, ε)  = ∫((u ⋅ ∇ᵈ(e, nΓ)) * ε + tr(∇ᵈ(u, nΓ)) * (e * ε))dΓ
  r(e, ε)  = ∫(τᵈkₒ * (e * ε))dΓ
  l(ε)     = ∫(τᵈkₒ * ε)dΓ
  s(υ, μ)  = ∫(γ * ((nΓ ⋅ ∇(υ)) ⊙ (nΓ ⋅ ∇(μ))))dΩᶜ

  aᵉ(e, ε) = m(e, ε) + c(e, ε) + r(e, ε) + sᵈ(e, ε) + s(e, ε)
  bᵉ(ε)    = m(eₕ, ε) + l(ε)

  aᵉ, bᵉ
end


# ============================================================
#   3D cortical flow problem
# ============================================================

"""3D cortical flow with full rigid-body mode elimination (6 translations/rotations + volume)."""
function cortical_flow_problem_3D(
    ulₕ, plₕ, eₕ, dΩᶜ, dΓ, nΓ,
    γ::Float64, Pe::Float64, μˡ::Float64, R::Float64, ξ₀::Function)

  β = 0.001
  aʷ(υ, μ) = ∫(2.0 * (εᶜ(υ, nΓ) ⊙ εᵈ(μ, nΓ)) + β * (υ ⋅ μ))dΓ

  ξ(e) = 2.0 * e * e / (1.0 + e * e)
  f(μ, e) = ∫(Pe * (-(divᶜ(μ, nΓ)) * (ξ ∘ (e))) * ξ₀)dΓ

  σᵘ(ε, q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ, u, p) = ∫(μ ⋅ ((σᵘ ∘ (ε(u), p)) ⋅ nΓ))dΓ

  sᵘ(υ, μ) = ∫(γ * ((nΓ ⋅ ε(υ)) ⊙ (nΓ ⋅ ε(μ))))dΩᶜ

  # Rigid-body modes: 3 translations + 3 rotations + 1 volume constraint
  RB¹ = VectorValue(1.0, 0.0, 0.0)
  RB² = VectorValue(0.0, 1.0, 0.0)
  RB³ = VectorValue(0.0, 0.0, 1.0)
  RB⁴ = x -> VectorValue(0.0, -x[3],  x[2])
  RB⁵ = x -> VectorValue( x[3], 0.0, -x[1])
  RB⁶ = x -> VectorValue(-x[2],  x[1], 0.0)

  r¹(u, ℓ) = ∫((u ⋅ RB¹) * ℓ)dΓ
  r²(u, ℓ) = ∫((u ⋅ RB²) * ℓ)dΓ
  r³(u, ℓ) = ∫((u ⋅ RB³) * ℓ)dΓ
  r⁴(u, ℓ) = ∫((u ⋅ RB⁴) * ℓ)dΓ
  r⁵(u, ℓ) = ∫((u ⋅ RB⁵) * ℓ)dΓ
  r⁶(u, ℓ) = ∫((u ⋅ RB⁶) * ℓ)dΓ
  r⁷(u, ℓ) = ∫((u ⋅ nΓ) * ℓ)dΓ

  aᵛ((υ, l¹, l², l³, l⁴, l⁵, l⁶, l⁷), (μ, ℓ¹, ℓ², ℓ³, ℓ⁴, ℓ⁵, ℓ⁶, ℓ⁷)) =
    aʷ(υ, μ) + sᵘ(υ, μ) +
    r¹(υ, ℓ¹) + r¹(μ, l¹) + r²(υ, ℓ²) + r²(μ, l²) +
    r³(υ, ℓ³) + r³(μ, l³) + r⁴(υ, ℓ⁴) + r⁴(μ, l⁴) +
    r⁵(υ, ℓ⁵) + r⁵(μ, l⁵) + r⁶(υ, ℓ⁶) + r⁶(μ, l⁶) +
    r⁷(υ, ℓ⁷) + r⁷(μ, l⁷)

  bᵛ((μ, ℓ¹, ℓ², ℓ³, ℓ⁴, ℓ⁵, ℓ⁶, ℓ⁷)) = f(μ, eₕ) - βʳ(μ, ulₕ, plₕ)

  aᵛ, bᵛ
end