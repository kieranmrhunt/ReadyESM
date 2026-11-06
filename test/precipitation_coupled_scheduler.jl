module PrecipitationCoupledSchedulerTests
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans
const CR = R.ConservativeRegridding

# SpeedyWeather invokes this after atmosphere physics and Terrarium's land
# timestep, before NumericalEarth regrids the atmosphere onto the ocean.
mutable struct PrecipitationStepObserver <: S.AbstractCallback
    samples::Vector{Any}
end
S.initialize!(::PrecipitationStepObserver, args...) = nothing
S.finalize!(::PrecipitationStepObserver, args...) = nothing

function S.callback!(observer::PrecipitationStepObserver, vars::S.Variables, model::S.AbstractModel)
    p = vars.parameterizations
    rain = vec(Array(p.rain_rate.data))
    convection = vec(Array(p.rain_rate_convection.data))
    large_scale = vec(Array(p.rain_rate_large_scale.data))
    snow = vec(Array(p.snow_rate.data))
    snow_large_scale = vec(Array(p.snow_rate_large_scale.data))
    snow_convection = vec(Array(p.snow_rate_convection.data))
    land_indices = vec(Array(model.land.model.grid.mask.data))
    land_rainfall = vec(Array(O.interior(vars.prognostic.land.terrarium.inputs.rainfall)))
    land_snowfall = vec(Array(O.interior(vars.prognostic.land.terrarium.inputs.snowfall)))
    push!(observer.samples, (; rain, convection, large_scale, snow,
                               snow_convection, snow_large_scale,
                               land_indices, land_rainfall, land_snowfall))
    return nothing
end

fields(x) = (; (name => getfield(x, name) for name in fieldnames(typeof(x)))...)
with(config; kwargs...) = R.ExperimentConfig(; merge(fields(config), (; kwargs...))...)
mass_flux(rate) = Float64(R._FRESHWATER_DENSITY_KG_M3) .* Float64.(rate)

@testset "two coupled steps preserve precipitation before freshwater exchange" begin
    mktempdir() do dir
        released = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml");
                                 check_input_files = false)
        config = with(released;
            device = :cpu, truncation = 7, nlayers = 8,
            ocean_grid = :idealized, ocean_nlongitude = 32,
            ocean_nlatitude = 16, ocean_nlayers = 4,
            ocean_initial_conditions = :analytic,
            ecco_initial_conditions_directory = "",
            atmosphere_initial_conditions = :analytic,
            terrarium_initial_conditions = :hydrostatic,
            terrarium_evapotranspiration = :bare_ground,
            atmosphere_dissipative_heating = :kinetic_energy_v1,
            terrarium_surface_albedo = :constant,
            terrarium_albedo_cycle_path = "",
            terrarium_albedo_cycle_sha256 = "",
            terrarium_vegetation_cycle = :off,
            terrarium_vegetation_cycle_path = "",
            terrarium_vegetation_cycle_sha256 = "",
            era5_pressure_levels_path = "",
            era5_single_levels_path = "",
            era5_land_state_path = "",
            duration_days = 2 / 24,
            analysis_start_day = 0.0,
            output_dir = joinpath(dir, "output"),
            atmosphere_vertical_advection = :centered_full_temperature,
        )
        R.validate(config; check_input_files = false)
        simulation = R.build_dynamic_esm(config; steps = 2)
        atmosphere_observer = PrecipitationStepObserver(Any[])
        S.add!(simulation.model.atmosphere.model.callbacks,
               :precipitation_scheduler_test => atmosphere_observer)
        ocean_samples = Any[]
        function observe_ocean(sim)
            sim.model.clock.iteration in 1:2 || return nothing
            exchanger = sim.model.interfaces.exchanger.atmosphere
            p = sim.model.atmosphere.variables.parameterizations
            component_sum = p.rain_rate_convection.data .+ p.rain_rate_large_scale.data
            expected = deepcopy(exchanger.state.Jʳⁿ)
            CR.regrid!(expected, exchanger.regridder.from_atmosphere, component_sum)
            R._mirror_right_center_fold_partners!(expected)
            actual_mass = vec(Array(O.interior(exchanger.state.Jʳⁿ)))
            expected_mass = vec(Array(O.interior(expected))) .* 1000
            push!(ocean_samples, (; iteration = sim.model.clock.iteration,
                                   actual_mass, expected_mass))
            return nothing
        end
        O.add_callback!(simulation, observe_ocean, O.IterationInterval(1);
                        name = :precipitation_ocean_test)
        O.run!(simulation)
        @test simulation.model.clock.iteration == 2
        @test length(atmosphere_observer.samples) == 2
        @test length(ocean_samples) == 2
        @test any(sample -> any(>(0), sample.convection), atmosphere_observer.samples)
        for (step, sample) in enumerate(atmosphere_observer.samples)
            # The production validator's absolute floor is in kg m⁻² s⁻¹;
            # passing native m s⁻¹ rates would weaken it by a factor of 1000.
            R._validate_precipitation_phase_components(
                mass_flux(sample.rain), mass_flux(sample.convection),
                mass_flux(sample.large_scale), mass_flux(sample.snow),
                mass_flux(sample.snow_convection), mass_flux(sample.snow_large_scale))
            @test sample.land_rainfall == sample.rain[sample.land_indices] .+
                                          sample.snow[sample.land_indices]
            @test sample.land_snowfall == sample.snow[sample.land_indices]
            ocean = ocean_samples[step]
            @test ocean.iteration == step
            @test all(isfinite, ocean.actual_mass)
            @test all(isapprox.(ocean.actual_mass, ocean.expected_mass;
                                atol = 1e-7, rtol = 1e-5))
            println("COUPLED_PRECIP_STEP_PASS step=$step max_ocean_residual=",
                    maximum(abs.(ocean.actual_mass .- ocean.expected_mass)))
        end
        R._release_dynamic_runtime_caches!(simulation)
    end
end
println("COUPLED_PRECIPITATION_TWO_STEP_PASS")
end # module
