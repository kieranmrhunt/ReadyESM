module SampledGlobalDiagnosticsTests
using ReadyESM, Test
const R = ReadyESM

@testset "Only sampled callback records are output" begin
    temperature = R.SpeedyWeather.GlobalSurfaceTemperatureCallback{Float32}(;
        timestep_counter=3,temperature=Float32[280,281,NaN,0,-999,9999,NaN])
    radiation = R.GlobalRadiationBudgetCallback{Float32,Vector{Float32}}(;
        timestep_counter=3,point_weights=Float32[1])
    for name in R._RADIATION_BUDGET_HISTORY_FIELDS
        setproperty!(radiation,name,Float32[NaN,240,NaN,0,-999,9999,NaN])
    end
    model = (;callbacks=Dict{Symbol,Any}(:global_surface_temperature=>temperature,
        :global_radiation_budget=>radiation),time_stepping=(Δt_sec=900.,))
    history = R._sampled_global_callback_diagnostics(model)
    @test history.time_days == [0.,900.,1800.] ./ 86400
    @test isequal(history.global_surface_temperature,[280.,281.,NaN])
    @test all(length(record) == 3 for record in values(history))
    for name in R._RADIATION_BUDGET_HISTORY_FIELDS
        @test isequal(getproperty(history,Symbol("toa_",name)),[NaN,240.,NaN])
    end
    history.global_surface_temperature[1] = -10
    @test temperature.temperature[1] == 280
    @test length(temperature.temperature) == 7
    @test length(radiation.net_downward) == 7
    temperature.timestep_counter = 8
    @test_throws ErrorException R._sampled_global_callback_diagnostics(model)
    temperature.timestep_counter = -1
    @test_throws ErrorException R._sampled_global_callback_diagnostics(model)
    temperature.timestep_counter = 3
    radiation.timestep_counter = 2
    @test_throws ErrorException R._sampled_global_callback_diagnostics(model)
    radiation.timestep_counter = 8
    @test_throws ErrorException R._sampled_global_callback_diagnostics(model)
    radiation.timestep_counter = 3
    radiation.net_downward = Float32[1,2]
    @test_throws ErrorException R._sampled_global_callback_diagnostics(model)
    radiation.net_downward = Float32[NaN,240,NaN]
    temperature.temperature = Float32[280,281,NaN]
    @test isequal(R._sampled_global_callback_diagnostics(model).global_surface_temperature,
        Float64.(temperature.temperature))
    delete!(model.callbacks,:global_radiation_budget)
    @test isempty(R._sampled_global_callback_diagnostics(model).toa_net_downward)
    @test length(R._sampled_global_callback_diagnostics(model).time_days) == 3
    temperature.timestep_counter = 0
    @test all(isempty,values(R._sampled_global_callback_diagnostics(model)))
end
end
