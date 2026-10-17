using DeBellyAxisymmetric
using Test

@testset "axisymmetric invariants" begin
    baseline = run_axisymmetric(n=40, duration=20.0,
        params=ModelParams(opto_start_step=30, opto_stop_step=40))
    @test maximum(baseline.rac[end, :]) - minimum(baseline.rac[end, :]) < 1e-4
    @test maximum(baseline.rho[end, :]) - minimum(baseline.rho[end, :]) < 1e-4
    @test baseline.mca_relative_error < 5e-4
    @test minimum(baseline.mcab) > 0
    @test minimum(baseline.mcau) > 0
end

@testset "published trajectory checkpoints" begin
    sim = run_axisymmetric(n=80, duration=100.0)
    # Reference extrema sampled from the authors' supplied _sur_75/_sur_100 VTUs.
    @test collect(extrema(sim.rac[75, :])) ≈ [0.6219, 1.3801] atol=3e-3
    @test collect(extrema(sim.rho[75, :])) ≈ [0.7397, 1.0498] atol=3e-3
    @test collect(extrema(sim.rac[100, :])) ≈ [0.4861, 1.6932] atol=4e-3
    @test collect(extrema(sim.rho[100, :])) ≈ [0.6401, 1.3613] atol=4e-3
    @test maximum(hypot.(sim.velocity_x[100, :], sim.velocity_y[100, :])) ≈ 0.0111 atol=1e-3
    @test maximum(hypot.(sim.displacement_x[100, :], sim.displacement_y[100, :])) ≈ 0.0874 atol=1e-2
end
