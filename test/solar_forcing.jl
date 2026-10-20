using Dates, Statistics

@testset "Seasonal solar forcing" begin
    R = ReadyESM
    @test !R.ForcingConfig().solar_distance_correction
    @test R._solar_distance_factor(DateTime(1993,1,1)) ≈ 1.03505 atol=1e-14
    @test R._solar_distance_factor(DateTime(1993,1,1)+Hour(365*12)) ≈ 0.966608 atol=1e-14
    for year in (1993,1996,2000,2100)
        start,stop = DateTime(year,1,1),DateTime(year+1,1,1)
        @test abs(R._solar_distance_factor(stop-Millisecond(1))-R._solar_distance_factor(stop)) < 1e-9
        values = R._solar_distance_factor.(start:Day(1):stop-Day(1))
        @test mean(values) ≈ 1.000110 atol=1e-14
        @test .966 < minimum(values) < .968
        @test 1.034 < maximum(values) < 1.036
    end
    for T in (Float32,Float64)
        for time in (DateTime(1993,1,1),DateTime(1993,7,1),DateTime(2000,2,29))
            @test R._solar_irradiance(T(1365),time,false) === T(1365)
            @test R._solar_irradiance(T(1365),time,true) isa T
        end
        @test R._solar_irradiance(T(1365),DateTime(1993,1,1),true) > T(1365)
        @test R._solar_irradiance(T(1365),DateTime(1993,7,1),true) < T(1365)
    end
    times = collect(DateTime(1993,1,1):Hour(1):DateTime(1993,1,3))
    @test R._solar_irradiance.(1365f0,times,true) == vcat(
        R._solar_irradiance.(1365f0,times[1:17],true),
        R._solar_irradiance.(1365f0,times[18:end],true))
    for mode in (:rrtmgp_clear_sky,:rrtmgp_all_sky)
        config = R.ExperimentConfig(forcing=R.ForcingConfig(radiation=mode,solar_distance_correction=true))
        @test R.validate(config;check_input_files=false) === config
    end
    unsupported = R.ExperimentConfig(forcing=R.ForcingConfig(solar_distance_correction=true))
    @test_throws ArgumentError R.validate(unsupported;check_input_files=false)
    mktemp() do path,io
        write(io,"forcing:\n  radiation: rrtmgp_clear_sky\n  solar_distance_correction: true\n")
        flush(io)
        @test R.load_config(path;check_input_files=false).forcing.solar_distance_correction
    end
    @test !R.load_config(joinpath(R.PROJECT_ROOT,"config/production.yml");check_input_files=false).forcing.solar_distance_correction
    for solar in (false,true), carbon in (:off,:emissions_driven)
        config = R.ExperimentConfig(forcing=R.ForcingConfig(solar_distance_correction=solar),
            carbon=R.CarbonConfig(mode=carbon))
        text = R._forcing_time_dependence(config)
        @test occursin("solar",text) == solar
        @test occursin("CO2",text) == (carbon==:emissions_driven)
    end
end

@testset "Solar constructor and restart contract" begin
    R = ReadyESM
    grid = R.SpeedyWeather.SpectralGrid(;NF=Float32,trunc=1,nlayers=8)
    rad = R.RRTMGPRadiation(grid)
    @test !rad.solar_distance_correction
    saved = R._rrtmgp_radiation_state(rad)
    @test saved.solar_distance_correction === false
    legacy = (; (key=>value for (key,value) in pairs(saved) if key!=:solar_distance_correction)...)
    @test R._restore_rrtmgp_radiation_state!(rad,legacy) === rad
    @test R._restore_rrtmgp_radiation_state!(rad,saved) === rad
    on = R.RRTMGPRadiation(grid;solar_distance_correction=true)
    @test on.solar_distance_correction
    active = R._rrtmgp_radiation_state(on)
    @test active.solar_distance_correction === true
    @test R._restore_rrtmgp_radiation_state!(on,active) === on
    for (radiation,state) in ((on,legacy),(on,saved),(rad,active))
        before = radiation.call_counter
        @test_throws ArgumentError R._restore_rrtmgp_radiation_state!(radiation,merge(state,(;call_counter=99)))
        @test radiation.call_counter == before
    end
    @test_throws ArgumentError R._validate_solar_restart(rad,(;solar_distance_correction=0))
    config = R.ExperimentConfig(forcing=R.ForcingConfig(radiation=:rrtmgp_clear_sky,solar_distance_correction=true))
    info = R._forcing_provenance(on,config)
    @test info.requested_solar_distance_correction
    @test info.applied_solar_distance_correction
    @test info.solar_distance_scheme == "spencer_present_day_gregorian_year_fraction_v1"
    @test info.forcing_time_dependence == "seasonal_solar_distance_from_model_clock"
end
