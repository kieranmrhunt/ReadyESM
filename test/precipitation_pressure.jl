module PrecipitationPressureTests
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather

@testset "Precipitation retains its physics-time pressure" begin
    config = R.ExperimentConfig(; name="precipitation_pressure", device=:cpu,
        truncation=7, nlayers=4, duration_days=0.1, mixed_layer_depth_m=30.)
    simulation, _ = R.run_baseline!(config; steps=2)
    vars, model = simulation.variables, simulation.model
    p = vars.parameterizations
    @test all(p.convective_precipitation_surface_pressure.data .> 0)
    @test all(p.large_scale_precipitation_surface_pressure.data .> 0)
    # Apply the real convection scheme to the evolved analytic atmosphere.
    # Isolate its increment from earlier processes: the last baseline step's
    # rain is near cancellation roundoff and cannot discriminate pressure errors.
    fill!(vars.tendencies.grid.humidity.data, 0)
    fill!(vars.tendencies.grid.temperature.data, 0)
    for name in (:rain_rate, :rain_rate_convection, :rain_convection)
        fill!(getproperty(p, name).data, 0)
    end
    for point in eachindex(vars.grid.pressure_prev)
        S.parameterization!(point, vars, model.convection, model)
    end
    sigma = reshape(Float64.(model.geometry.σ_levels_thick), 1, :)
    drying = -vec(sum(Float64.(p.convective_humidity_tendency.data) .* sigma; dims=2))
    rain = 1000 .* Float64.(p.rain_rate_convection.data)
    process_pressure = Float64.(p.convective_precipitation_surface_pressure.data)
    reconstructed = max.(process_pressure ./ model.planet.gravity .* drying, 0)
    @test maximum(rain) > 1e-4
    @test isapprox(rain, reconstructed; rtol=8eps(Float32), atol=0)

    # A later pressure update must not change the retained mass conversion.
    # Deliberately amplify the phase difference as a negative control.
    vars.grid.pressure_prev.data .= 0.8 .* process_pressure
    @test Float64.(p.convective_precipitation_surface_pressure.data) == process_pressure
    wrong = max.(Float64.(vars.grid.pressure_prev.data) ./ model.planet.gravity .* drying, 0)
    # The wrong phase must fail the production validator's absolute threshold
    # as well as a relative check on this substantial precipitation signal.
    @test maximum(abs.(rain .- wrong)) > 3e-6
    @test maximum(abs.(rain .- wrong)) > 0.1maximum(rain)
    println("PRECIPITATION_PRESSURE_RESULT max_rain_kgm2s=$(maximum(rain)) ",
        "max_reconstruction_error=$(maximum(abs.(rain .- reconstructed))) ",
        "wrong_phase_error=$(maximum(abs.(rain .- wrong)))")

    # Both condensation wrappers retain their input pressure. The observer
    # must leave every physical output identical to the unwrapped scheme.
    point = 1
    observer = model.large_scale_condensation
    physical = (vars.tendencies.grid.humidity, vars.tendencies.grid.temperature,
        p.rain_large_scale, p.snow_large_scale, p.rain_rate_large_scale,
        p.snow_rate_large_scale, p.rain_rate, p.snow_rate, p.cloud_top)
    before = map(f -> copy(f.data), physical)
    S.parameterization!(point, vars, observer.condensation, model)
    expected = map(f -> copy(f.data), physical)
    foreach(zip(physical, before)) do (field, saved)
        field.data .= saved
    end
    S.parameterization!(point, vars, observer, model)
    @test all(f.data == saved for (f,saved) in zip(physical, expected))
    @test p.large_scale_precipitation_surface_pressure[point] == vars.grid.pressure_prev[point]
    S.parameterization!(point, vars, R.NetColumnImplicitCondensation(observer.condensation), model)
    @test p.large_scale_precipitation_surface_pressure[point] == vars.grid.pressure_prev[point]
    S.parameterization!(point, vars, model.convection, model)
    @test p.convective_precipitation_surface_pressure[point] == vars.grid.pressure_prev[point]

    mktempdir() do directory
        path = joinpath(directory, "physics_pressure.jld2")
        R.JLD2.jldsave(path; parameterizations=p)
        restored = R.JLD2.load(path, "parameterizations")
        @test restored.convective_precipitation_surface_pressure.data ==
            p.convective_precipitation_surface_pressure.data
        @test restored.large_scale_precipitation_surface_pressure.data ==
            p.large_scale_precipitation_surface_pressure.data
    end
end
end
