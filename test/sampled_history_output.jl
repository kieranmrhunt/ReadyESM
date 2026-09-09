module SampledHistoryOutputTests
using ReadyESM, Test
const R = ReadyESM

@testset "Recorded history in atmosphere output" begin
    for mode in (:speedy_simplified,:rrtmgp_clear_sky)
        config = R.ExperimentConfig(;name="sampled_history_$mode",device=:cpu,
            truncation=7,nlayers=8,duration_days=0.1,mixed_layer_depth_m=30.,
            forcing=R.ForcingConfig(;radiation=mode))
        runner = mode == :speedy_simplified ? R.run_baseline! : R.run_rrtmgp!
        simulation, complete = runner(config;steps=2)
        @test complete.metadata.active_radiation == String(mode)
        @test length(complete.time_days) == 3
        @test all(isfinite,complete.global_surface_temperature)
        callback = simulation.model.callbacks[:global_surface_temperature]
        @test callback.timestep_counter == 3
        append!(callback.temperature,Float32[0,NaN,-999,9999])
        radiation = get(simulation.model.callbacks,:global_radiation_budget,nothing)
        if !isnothing(radiation)
            for name in R._RADIATION_BUDGET_HISTORY_FIELDS
                append!(getproperty(radiation,name),Float32[0,NaN,-999,9999])
            end
        end
        sampled = R.collect_diagnostics(simulation,config)
        @test isequal(sampled.time_days,complete.time_days)
        @test isequal(sampled.global_surface_temperature,complete.global_surface_temperature)
        for name in R._RADIATION_BUDGET_HISTORY_FIELDS
            @test isequal(getproperty(sampled,Symbol("toa_",name)),
                getproperty(complete,Symbol("toa_",name)))
        end
        @test length(callback.temperature) == 7
        temporary = R.SpeedyWeather.GlobalSurfaceTemperatureCallback(simulation.model.spectral_grid)
        R.SpeedyWeather.initialize!(temporary,simulation.variables,simulation.model)
        @test temporary.timestep_counter == 1
        @test isfinite(temporary.temperature[1])
        @test all(isnan,temporary.temperature[2:end])
        # Genuine recorded invalid values must remain visible.
        callback.temperature[3] = NaN
        invalid = R.collect_diagnostics(simulation,config)
        @test length(invalid.time_days) == 3
        @test isnan(invalid.global_surface_temperature[3])
        mktempdir() do directory
            path = joinpath(directory,"partial.nc")
            R._write_netcdf(path,invalid)
            R.NCDataset(path,"r") do dataset
                @test length(dataset["time"]) == 3
                @test isnan(dataset["global_surface_temperature"][3])
                @test haskey(dataset,"toa_net_downward") == !isnothing(radiation)
                @test dataset["global_surface_temperature"].attrib["long_name"] ==
                    R._GLOBAL_AIR_TEMPERATURE_LONG_NAME
            end
        end
    end
end
end
