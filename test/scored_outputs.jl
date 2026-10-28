module ScoredOutputsTests
# Opt-in scored output files (src/scored_outputs.jl): configuration, the monthly and wind accumulators, a short coupled
# run that writes all four files, the restart refusal and the unchanged option-off build.
using ReadyESM, Test, Dates, NCDatasets
const R = ReadyESM
const O = R.Oceananigans

fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
include(joinpath(@__DIR__, "analytic_coupled_fixture.jl"))

@testset "scored outputs: configuration" begin
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    @test production.scored_outputs == false
    @test R.ExperimentConfig().scored_outputs == false
    simplified = R.ForcingConfig(; merge(fields(production.forcing),
        (; radiation = :speedy_simplified, solar_distance_correction = false, cloud_scheme = :none))...)
    err = try
        R.validate(with(production; scored_outputs = true, forcing = simplified); check_input_files = false)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("scored_outputs requires RRTMGP", err.msg)
    @test (R.validate(with(production; forcing = simplified); check_input_files = false); true)
end

@testset "scored outputs: accumulators" begin
    M = R.ScoredMonthlyCore
    sigma = [0.25, 0.5, 0.75, 0.95]
    @test M.pressure_temperature([200.0, 250.0, 270.0, 290.0], sigma, 1.0e5, 50000.0) == 250.0
    @test isnothing(M.pressure_temperature([200.0, 250.0, 270.0, 290.0], sigma, 1.0e5, 99000.0))
    @test isnothing(M.pressure_temperature([200.0, 250.0, 270.0, 290.0], sigma, 1.0e5, 20000.0))
    mid = M.pressure_temperature([200.0, 250.0, 270.0, 290.0], sigma, 1.0e5, 60000.0)
    @test 250.0 < mid < 270.0
    start = DateTime(1993, 1, 1)
    a = M.Accumulator(start, start + Hour(2), 3, sigma, fill(0.25, 4), 9.81, 1000.0; step_ms = 1800000)
    T = fill(260.0, 3, 4)
    q = fill(1.0e-3, 3, 4)
    ps = fill(1.0e5, 3)
    M.state!(a, start, T, q, ps)
    @test_throws ErrorException M.state!(a, start, T, q, ps)               # duplicate hour
    for k in 1:4
        M.flux!(a, start + Minute(30k), fill(1.0e-7, 3), zeros(3))
    end
    M.state!(a, start + Hour(1), T, q, ps)
    data = M.finish(a)
    @test data.state_samples == 2 && data.flux_samples == 4
    @test all(data.surface_pressure_Pa .== 1.0e5)
    @test data.precipitation_kg_m2 ≈ fill(4 * 1.0e-7 * 1000.0 * 1800.0, 3)
    @test all(data.temperature_valid[:, 2])                                # 500 hPa inside the column
    @test_throws ErrorException M.finish(M.Accumulator(start, start + Hour(1), 3, sigma, fill(0.25, 4), 9.81, 1000.0))
    W = R.ScoredWindCore
    w = W.Accumulator(start, start + Hour(1), 3, sigma)
    W.state!(w, start, fill(10.0, 3, 4), fill(-5.0, 3, 4), ps)
    wd = W.finish(w)
    @test wd.samples == 1 && all(wd.u[wd.valid] .== 10.0) && all(wd.v[wd.valid] .== -5.0)
end

@testset "scored outputs: coupled run" begin
    released = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
    output = mktempdir()
    base = analytic_coupled_test_config(released; output_dir = output)
    # 3.75 min at T31 gives a 900 s step at T7. Six hours give 24 steps, six hourly states and one six-hourly sample.
    config = with(base; atmosphere_timestep_at_t31_minutes = 3.75, duration_days = 0.25,
        start_date = DateTime(1993, 1, 1), scored_outputs = true)
    R.validate(config; check_input_files = false)
    simulation = R.build_dynamic_esm(config)
    @test haskey(simulation.callbacks, R.SCORED_OUTPUTS_CALLBACK)
    @test_throws ArgumentError R.restore_dynamic_restart_state!(simulation, joinpath(output, "missing.jld2"))
    O.run!(simulation)
    summary = R.save_scored_outputs(simulation, config)
    @test summary["monthly"]["state_samples"] == 6
    @test summary["monthly"]["flux_samples"] == 24
    @test summary["wind"]["state_samples"] == 6
    @test summary["toa"]["samples"] == 6
    @test summary["sixhourly"]["samples"] == 1
    @test summary["passivity_checks"] == 25
    for name in ("monthly.nc", "wind.nc", "toa_monthly.nc", "eddy_6hourly.nc")
        @test isfile(joinpath(output, name))
    end
    NCDataset(joinpath(output, "monthly.nc")) do ds
        @test ds.attrib["start_time_UTC"] == "1993-01-01T00:00:00"
        @test ds.attrib["stop_time_UTC"] == "1993-01-01T06:00:00"
        @test all(isfinite, ds["surface_pressure"][:])
    end
    @test_throws ErrorException R.save_scored_outputs(simulation, config)  # refuses to overwrite
    off = R.build_dynamic_esm(with(config; scored_outputs = false, output_dir = mktempdir()))
    @test !haskey(off.callbacks, R.SCORED_OUTPUTS_CALLBACK)
end
end
