using DeBellyAxisymmetric

n = parse(Int, get(ENV, "DEBELLY_AXIS_N", "80"))
duration = parse(Float64, get(ENV, "DEBELLY_T", "350"))
saveat = parse(Float64, get(ENV, "DEBELLY_SAVEAT", "1"))
writesol = get(ENV, "DEBELLY_WRITE_OUTPUTS", "1") == "1"

params = ModelParams(
    opto_start_step=round(Int, parse(Float64, get(ENV, "DEBELLY_OPTO_START", "50"))),
    opto_stop_step=round(Int, parse(Float64, get(ENV, "DEBELLY_OPTO_STOP", "150"))),
)

elapsed = @elapsed sim = run_axisymmetric(; n=n, duration=duration,
                                          saveat=saveat, params=params)
if writesol
    output = get(ENV, "DEBELLY_OUTPUT",
                 joinpath(@__DIR__, "..", "output", "axisymmetric.csv"))
    write_csv(output, sim)
    println("Wrote $output")
end

rac_final = @view sim.rac[end, :]
rho_final = @view sim.rho[end, :]
println("DEBELLY_SUMMARY protocol=axisymmetric n=$n steps=$(length(sim.time)) " *
        "elapsed_seconds=$elapsed rac_range=$(maximum(rac_final)-minimum(rac_final)) " *
        "rho_range=$(maximum(rho_final)-minimum(rho_final)) " *
        "mca_relative_error=$(sim.mca_relative_error) " *
        "rac_front=$(rac_final[end]) rac_back=$(rac_final[1])")
