import PrecompileTools

function _precompile_mode(raw::AbstractString)
    raw = lowercase(strip(raw))
    if raw in ("0", "false", "no", "off", "none")
        return :off
    elseif raw in ("full", "all", "coupled")
        return :full
    elseif raw in ("", "bounded", "default", "1", "true", "yes", "on")
        return :bounded
    else
        throw(ArgumentError("unknown ReadyESM precompile workload mode $(repr(raw)); expected bounded, off, or full"))
    end
end
_precompile_mode() = _precompile_mode(get(ENV, "READYESM_PRECOMPILE_WORKLOAD_MODE",
                                        get(ENV, "READYESM_PRECOMPILE_WORKLOAD", "bounded")))
const _precompile_workload_mode = _precompile_mode

function _clear_precompile_global_caches!()
    empty!(_BALANCED_LAUNCH_CONTROLLERS)
    empty!(_BS08_DIAGNOSTIC_HISTORIES)
    empty!(_CARBON_CONTROLLERS)
    empty!(_CARBON_DIAGNOSTIC_HISTORIES)
    empty!(_DETERMINISTIC_LEGENDRE_MMAX_CACHE)
    empty!(_DYNAMIC_RESTART_IDENTITY_BY_MODEL)
    _DYNAMIC_RESTART_IDENTITY_OVERRIDE[] = false
    empty!(_DYNAMIC_RESTART_PENDING_STATE)
    empty!(_SEA_ICE_MOMENTUM_FORCE_OBSERVERS)
    empty!(_INTENSIVE_REGRID_CACHE)
    empty!(_TERRARIUM_LAND_INDEX_CACHE)
    return nothing
end

function _precompile_config_with(config::ExperimentConfig; kwargs...)
    base_fields = (; (name => getfield(config, name) for name in fieldnames(ExperimentConfig))...)
    return ExperimentConfig(; merge(base_fields, (; kwargs...))...)
end

PrecompileTools.@setup_workload begin
    mode = _precompile_workload_mode()
    production_config_path = joinpath(PROJECT_ROOT, "config", "production.yml")
    if mode != :off && isfile(production_config_path)
        PrecompileTools.@compile_workload begin
            Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                tmp_root = mktempdir()
                try
                    released = load_config(production_config_path; check_input_files = false)
                    validate(released; check_input_files = false)

                    analytic_prod = _precompile_config_with(
                        released;
                        device = :cpu,
                        truncation = 7,
                        nlayers = 8,
                        ocean_grid = :idealized,
                        ocean_nlongitude = 32,
                        ocean_nlatitude = 16,
                        ocean_nlayers = 4,
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
                        duration_days = 1 / 24,
                        analysis_start_day = 0.0,
                        output_dir = joinpath(tmp_root, "prod"),
                        atmosphere_vertical_advection = :centered_full_temperature,
                    )
                    validate(analytic_prod; check_input_files = false)

                    analytic_off = _precompile_config_with(
                        analytic_prod;
                        atmosphere_dissipative_heating = :off,
                        output_dir = joinpath(tmp_root, "off"),
                    )
                    validate(analytic_off; check_input_files = false)

                    spectral_grid = SpeedyWeather.SpectralGrid(;
                        NF = Float32,
                        trunc = 7,
                        Grid = SpeedyWeather.RingGrids.FullGaussianGrid,
                        nlayers = 8,
                    )
                    hd_on = _atmosphere_horizontal_diffusion(analytic_prod, spectral_grid)
                    _plain_horizontal_diffusion(hd_on)
                    _atmosphere_speed_limit_drag(analytic_prod, spectral_grid)
                    _atmosphere_horizontal_diffusion(analytic_off, spectral_grid)
                    _atmosphere_speed_limit_drag(analytic_off, spectral_grid)

                    _rrtmgp_atmosphere_model(analytic_prod; coupled_ocean = true)
                    _rrtmgp_atmosphere_model(analytic_off; coupled_ocean = true)

                    precompile(Tuple{typeof(build_dynamic_esm), ExperimentConfig})
                    precompile(
                        Tuple{
                            typeof(Core.kwcall),
                            NamedTuple{(:steps,), Tuple{Int}},
                            typeof(build_dynamic_esm),
                            ExperimentConfig,
                        },
                    )
                    precompile(Tuple{typeof(run_dynamic_esm!), ExperimentConfig})
                    precompile(Tuple{typeof(_build_ocean), ExperimentConfig})

                    if mode == :full
                        sim_prod = build_dynamic_esm(analytic_prod; steps = 1)
                        Oceananigans.run!(sim_prod)
                        diag_prod = collect_dynamic_diagnostics(sim_prod, analytic_prod)
                        validate_dynamic_diagnostics(diag_prod)
                        _release_dynamic_runtime_caches!(sim_prod)

                        sim_off = build_dynamic_esm(analytic_off; steps = 1)
                        _release_dynamic_runtime_caches!(sim_off)
                    end
                finally
                    _clear_precompile_global_caches!()
                    rm(tmp_root; force = true, recursive = true)
                end
            end
        end
    end
end
