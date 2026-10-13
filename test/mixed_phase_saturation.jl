module MixedPhaseSaturationTests
using ReadyESM, Test, Adapt
const R = ReadyESM
const SW = R.SpeedyWeather

with(config; kwargs...) = R.ExperimentConfig(; merge(
    (; (n => getfield(config, n) for n in fieldnames(typeof(config)))...), (; kwargs...))...)

"""Large-scale condensation tendency of one cold column, 90% saturated over liquid."""
function cold_column_tendency(atmosphere)
    spectral_grid = SW.SpectralGrid(trunc = 5, nlayers = 8)
    model = SW.PrimitiveWetModel(spectral_grid; atmosphere)
    vars = SW.initialize!(model).variables
    σ = model.geometry.σ_levels_full
    earth = SW.EarthAtmosphere(spectral_grid)
    vars.grid.pressure_prev[1] = 1f5
    for k in 1:8
        T = 215f0 + 9f0 * (k - 1)
        vars.grid.temperature_prev[1, k] = T
        vars.grid.humidity_prev[1, k] = 0.9f0 * SW.saturation_humidity(T, σ[k] * 1f5, earth)
        vars.tendencies.grid.humidity[1, k] = 0
        vars.tendencies.grid.temperature[1, k] = 0
    end
    vars.parameterizations.cloud_top[1] = 9
    SW.large_scale_condensation!(1, vars, model.large_scale_condensation,
        model.geometry, model.planet, model.atmosphere, model.time_stepping)
    return [vars.tendencies.grid.humidity[1, k] for k in 1:8],
           [vars.grid.temperature_prev[1, k] for k in 1:8]
end

@testset "Mixed-phase saturation option" begin
    # Preserve old configurations; production explicitly selects the D pair.
    @test R.ExperimentConfig().atmosphere_saturation_law == :liquid
    @test R.ExperimentConfig().atmosphere_large_scale_precipitation == :upstream_observed
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml");
        check_input_files = false)
    @test production.atmosphere_saturation_law == :mixed_phase
    @test production.atmosphere_large_scale_precipitation == :net_column_conservative
    legacy = with(production; atmosphere_saturation_law = :liquid,
        atmosphere_large_scale_precipitation = :upstream_observed)
    mktempdir() do directory
        path = joinpath(directory, "config.yml")
        write(path, "atmosphere_saturation_law: mixed_phase\n")
        @test R.load_config(path; check_input_files = false).atmosphere_saturation_law == :mixed_phase
    end
    @test_throws ArgumentError R.validate(with(production; atmosphere_saturation_law = :ice);
        check_input_files = false)
    @test R.validate(with(production; atmosphere_saturation_law = :mixed_phase);
        check_input_files = false) isa R.ExperimentConfig

    grid = SW.SpectralGrid(trunc = 31, nlayers = 27)
    liquid_atmosphere = R._atmosphere_thermodynamics(legacy, grid)
    mixed_atmosphere = R._atmosphere_thermodynamics(production, grid)
    @test R._atmosphere_large_scale_condensation(production, grid) isa R.NetColumnImplicitCondensation
    @test !(R._atmosphere_large_scale_condensation(legacy, grid) isa R.NetColumnImplicitCondensation)
    earth = SW.EarthAtmosphere(grid)
    @test liquid_atmosphere isa SW.EarthAtmosphere{Float32}
    @test all(isequal(getfield(liquid_atmosphere, n), getfield(earth, n)) for n in fieldnames(typeof(earth)))
    @test mixed_atmosphere isa R.MixedPhaseEarthAtmosphere
    @test mixed_atmosphere isa SW.AbstractWetAtmosphere
    @test eltype(mixed_atmosphere) == Float32
    @test isbits(mixed_atmosphere)
    @test propertynames(mixed_atmosphere) == propertynames(earth)
    @test all(isequal(getproperty(mixed_atmosphere, n), getproperty(earth, n)) for n in propertynames(earth))
    @test Adapt.adapt(Array, mixed_atmosphere) === mixed_atmosphere
    @test R._atmosphere_saturation_provenance(legacy).atmosphere_saturation_law == "liquid"
    @test isnan(R._atmosphere_saturation_provenance(legacy).atmosphere_saturation_all_ice_temperature_K)
    mixed_provenance = R._atmosphere_saturation_provenance(production)
    @test mixed_provenance.atmosphere_saturation_law == "mixed_phase"
    @test mixed_provenance.atmosphere_saturation_all_ice_temperature_K == 250.16

    # The law: SpeedyWeather's own value at and above freezing, lower below it,
    # continuous at freezing and close to Murphy and Koop (2005) over ice.
    esat(T, A) = SW.saturation_vapor_pressure(T, A)
    T₀ = earth.temperature_freezing
    temperatures = Float32.(180:0.5:320)
    for T in temperatures
        if T >= T₀
            @test esat(T, mixed_atmosphere) === esat(T, earth)
            @test SW.saturation_humidity(T, 5f4, mixed_atmosphere) === SW.saturation_humidity(T, 5f4, earth)
        else
            @test esat(T, mixed_atmosphere) < esat(T, earth)
        end
    end
    @test all(diff(esat.(temperatures, Ref(mixed_atmosphere))) .> 0)
    @test isapprox(esat(prevfloat(T₀), mixed_atmosphere), esat(T₀, earth); rtol = 1f-5)
    @test 2.1 < esat(210f0, earth) / esat(210f0, mixed_atmosphere) < 2.3
    for T in 190.0:5.0:245.0
        murphy_koop = exp(9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T)
        @test isapprox(esat(Float32(T), mixed_atmosphere), murphy_koop; rtol = 0.03)
    end

    # SpeedyWeather's compiled large-scale condensation sees the active law:
    # cold layers condense only with the mixed law; warm layers are identical.
    liquid_dq, T = cold_column_tendency(SW.EarthAtmosphere(SW.SpectralGrid(trunc = 5, nlayers = 8)))
    mixed_dq, _ = cold_column_tendency(R.MixedPhaseEarthAtmosphere(SW.SpectralGrid(trunc = 5, nlayers = 8)))
    @test all(iszero, liquid_dq)
    @test all(mixed_dq[k] < 0 for k in 1:8 if T[k] < 265)
    @test all(mixed_dq[k] == liquid_dq[k] for k in 1:8 if T[k] >= 273.15)

    # End to end: a small ReadyESM model builds and steps with the mixed law.
    config = R.ExperimentConfig(device = :cpu, truncation = 7, nlayers = 8,
        land_model = :terrarium, terrarium_soil_layers = 3, terrarium_max_layer_thickness_m = 0.5,
        terrarium_timestep_seconds = 60, atmosphere_saturation_law = :mixed_phase,
        atmosphere_large_scale_precipitation = :net_column_conservative,
        forcing = R.ForcingConfig(radiation = :rrtmgp_clear_sky))
    simulation = R.build_rrtmgp(config)
    @test simulation.model.atmosphere isa R.MixedPhaseEarthAtmosphere
    @test simulation.model.large_scale_condensation isa R.NetColumnImplicitCondensation
    SW.run!(simulation; steps = 4)
    @test all(isfinite, Array(simulation.variables.grid.temperature_prev.data))
    @test all(isfinite, Array(simulation.variables.grid.humidity_prev.data))
end
println("MIXED_PHASE_SATURATION_UNIT_PASS")
end
