Base.@kwdef struct ForcingConfig
    co2_ppm::Float64 = 420.0
    aerosol_optical_depth_550nm::Float64 = 0.0
    aerosol_species::String = "sulfate"
    radiation::Symbol = :speedy_simplified
    radiation_every_n_steps::Int = 1
    cloud_scheme::Symbol = :none
    cloud_humidity_search_min_sigma::Float64 = 0.0
    cloud_liquid_water_path_gm2::Float64 = 60.0
    cloud_ice_water_path_gm2::Float64 = 25.0
    cloud_condensate_retention_fraction::Float64 = 0.0
    cloud_condensate_residence_time_hours::Float64 = 6.0
end

Base.@kwdef struct ExperimentConfig
    name::String = "baseline"
    device::Symbol = :cpu
    truncation::Int = 15
    nlayers::Int = 8
    duration_days::Float64 = 10.0
    analysis_start_day::Float64 = 0.0
    start_date::DateTime = DateTime(2000, 1, 1)
    atmosphere_timestep_at_t31_minutes::Float64 = 40.0
    atmosphere_hyperdiffusion_hours::Float64 = 4.0
    atmosphere_divergence_hyperdiffusion_hours::Float64 = 1.0
    # The metric-v2 operator is numerically stable but failed its matched
    # 30-day climate gate.  Keep it available only through an explicit
    # experiment selection; ordinary configurations inherit the qualified
    # SpeedyWeather 0.21.1 zero-transport control.
    atmosphere_vertical_diffusion::Symbol =
        :speedyweather_0_21_1_zero_operator_control
    atmosphere_large_scale_precipitation::Symbol = :upstream_observed
    atmosphere_convection::Symbol = :betts_miller_constant_rh
    atmosphere_convection_upper_relative_humidity::Float64 = 0.5
    atmosphere_convection_lower_relative_humidity::Float64 = 0.7
    atmosphere_convection_transition_top_sigma::Float64 = 0.3
    atmosphere_convection_transition_bottom_sigma::Float64 = 0.7
    atmosphere_speed_limit_ms::Float64 = 80.0
    atmosphere_speed_limit_drag_per_m::Float64 = 4.0e-7
    atmosphere_surface_speed_limit_ms::Float64 = 80.0
    atmosphere_surface_speed_limit_drag_per_m::Float64 = 4.0e-7
    mixed_layer_depth_m::Float64 = 50.0
    ocean_model::Symbol = :slab
    ocean_grid::Symbol = :idealized
    ocean_dynamics::Bool = false
    sea_ice_dynamics::Bool = false
    sea_ice_momentum_substeps::Int = 240
    sea_ice_pressure_formulation::Symbol = :replacement_pressure
    sea_ice_immersed_boundary_drag_coefficient::Float64 = 0.0
    ocean_nlongitude::Int = 32
    ocean_nlatitude::Int = 16
    ocean_nlayers::Int = 4
    ocean_depth_m::Float64 = 1000.0
    ocean_bathymetry_correction::Symbol = :none
    ocean_tripolar_wet_mask::Symbol = :interpolated_mean_elevation
    ocean_shortwave_scheme::Symbol = :two_color_default
    ocean_river_mouth_mixing::Symbol = :none
    ocean_river_mouth_vertical_diffusivity_m2s::Float64 = 1.0e-2
    ocean_river_mouth_horizontal_diffusivity_m2s::Float64 = 4.0e3
    ocean_river_mouth_mixing_depth_m::Float64 = 50.0
    ocean_river_mouth_reference_freshwater_mass_flux_kgm2s::Float64 = 0.0
    ocean_polar_sponge_start_latitude_degrees::Float64 = 70.0
    ocean_polar_sponge_stop_latitude_degrees::Float64 = 75.0
    ocean_polar_sponge_timescale_hours::Float64 = 3.0
    ocean_initial_conditions::Symbol = :analytic
    ecco_year::Int = 1993
    ecco_month::Int = 1
    ecco_initial_conditions_directory::String = ""
    atmosphere_initial_conditions::Symbol = :analytic
    era5_pressure_levels_path::String = ""
    era5_single_levels_path::String = ""
    era5_land_state_path::String = ""
    era5_month::Int = 1
    land_model::Symbol = :speedy_bucket
    terrarium_runoff_routing::Symbol = :equal_area_flux
    terrarium_initial_conditions::Symbol = :hydrostatic
    terrarium_evapotranspiration::Symbol = :bare_ground
    terrarium_max_leaf_conductance_ms::Float64 = 2.0e-3
    terrarium_min_leaf_conductance_ms::Float64 = 5.0e-5
    terrarium_light_half_saturation_wm2::Float64 = 100.0
    terrarium_vpd_scale_pa::Float64 = 2000.0
    terrarium_canopy_extinction_coefficient::Float64 = 0.5
    terrarium_root_a_m1::Float64 = 7.0
    terrarium_root_b_m1::Float64 = 2.0
    terrarium_soil_layers::Int = 8
    terrarium_max_layer_thickness_m::Float64 = 5.0
    terrarium_timestep_seconds::Float64 = 60.0
    output_dir::String = joinpath(PROJECT_ROOT, "artifacts", "baseline")
    forcing::ForcingConfig = ForcingConfig()
end

const _LOCALIZED_ESTUARY_MIXING_MODES = (
    :localized_estuary_v1,
    :localized_estuary_v2,
    :localized_estuary_v3,
    :localized_estuary_v4,
    :localized_estuary_v5,
)

const _ECCO_INITIAL_CONDITION_MODES = (
    :ecco_1993,
    :ecco_1993_full_state,
    :ecco_monthly,
    :ecco_monthly_full_state,
)

@inline _is_localized_estuary_mixing(mode::Symbol) =
    mode in _LOCALIZED_ESTUARY_MIXING_MODES

@inline _is_ecco_initial_conditions(mode::Symbol) =
    mode in _ECCO_INITIAL_CONDITION_MODES

@inline _is_full_state_ecco_initial_conditions(mode::Symbol) =
    mode in (:ecco_1993_full_state, :ecco_monthly_full_state)

@inline function _ecco_initial_condition_date(config::ExperimentConfig)
    return config.ocean_initial_conditions in (:ecco_1993, :ecco_1993_full_state) ?
        DateTime(1993, 1, 1) : DateTime(config.ecco_year, config.ecco_month, 1)
end

@inline _ecco_initial_condition_month_label(config::ExperimentConfig) =
    Dates.format(_ecco_initial_condition_date(config), dateformat"yyyy-mm")

function _mapping_value(mapping::AbstractDict, key::AbstractString, default)
    return get(mapping, key, get(mapping, Symbol(key), default))
end

function _validate_mapping_keys(mapping, allowed_keys, label)
    mapping isa AbstractDict || throw(ArgumentError(
        "$label must be a YAML mapping",
    ))
    allowed = Set(String.(allowed_keys))
    unknown = sort!(String[
        key isa AbstractString || key isa Symbol ? String(key) : repr(key)
        for key in keys(mapping)
        if !(key isa AbstractString || key isa Symbol) || String(key) ∉ allowed
    ])
    isempty(unknown) || throw(ArgumentError(
        "$label contains unknown key(s): $(join(unknown, ", "))",
    ))
    return mapping
end

function _validate_finite_float_fields(value, label)
    for field in fieldnames(typeof(value))
        field_value = getfield(value, field)
        field_value isa AbstractFloat || continue
        isfinite(field_value) || throw(ArgumentError(
            "$label field $field must be finite",
        ))
    end
    return value
end

function _resolve_optional_path(path_raw)
    path = String(path_raw)
    isempty(path) && return path
    return isabspath(path) ? path : joinpath(PROJECT_ROOT, path)
end

function _parse_start_date(value)
    value isa DateTime && return value
    value isa Date && return DateTime(value)
    try
        return DateTime(String(value))
    catch exception
        throw(ArgumentError(
            "start_date must be an ISO-8601 date or date-time: " *
            sprint(showerror, exception),
        ))
    end
end

"""Load an experiment configuration from YAML."""
function load_config(path::AbstractString)
    raw = YAML.load_file(path)
    _validate_mapping_keys(
        raw,
        fieldnames(ExperimentConfig),
        "experiment configuration",
    )
    forcing_raw = _mapping_value(raw, "forcing", Dict{String, Any}())
    _validate_mapping_keys(
        forcing_raw,
        fieldnames(ForcingConfig),
        "forcing configuration",
    )
    forcing = ForcingConfig(
        co2_ppm = Float64(_mapping_value(forcing_raw, "co2_ppm", 420.0)),
        aerosol_optical_depth_550nm = Float64(
            _mapping_value(forcing_raw, "aerosol_optical_depth_550nm", 0.0),
        ),
        aerosol_species = String(_mapping_value(forcing_raw, "aerosol_species", "sulfate")),
        radiation = Symbol(_mapping_value(forcing_raw, "radiation", "speedy_simplified")),
        radiation_every_n_steps = Int(
            _mapping_value(forcing_raw, "radiation_every_n_steps", 1),
        ),
        cloud_scheme = Symbol(_mapping_value(forcing_raw, "cloud_scheme", "none")),
        cloud_humidity_search_min_sigma = Float64(
            _mapping_value(forcing_raw, "cloud_humidity_search_min_sigma", 0.0),
        ),
        cloud_liquid_water_path_gm2 = Float64(
            _mapping_value(forcing_raw, "cloud_liquid_water_path_gm2", 60.0),
        ),
        cloud_ice_water_path_gm2 = Float64(
            _mapping_value(forcing_raw, "cloud_ice_water_path_gm2", 25.0),
        ),
        cloud_condensate_retention_fraction = Float64(
            _mapping_value(
                forcing_raw,
                "cloud_condensate_retention_fraction",
                0.0,
            ),
        ),
        cloud_condensate_residence_time_hours = Float64(
            _mapping_value(
                forcing_raw,
                "cloud_condensate_residence_time_hours",
                6.0,
            ),
        ),
    )

    output_dir_raw = String(
        _mapping_value(raw, "output_dir", joinpath("artifacts", "baseline")),
    )
    output_dir = isabspath(output_dir_raw) ? output_dir_raw : joinpath(PROJECT_ROOT, output_dir_raw)
    atmosphere_initial_conditions = Symbol(
        _mapping_value(raw, "atmosphere_initial_conditions", "analytic"),
    )
    ocean_initial_conditions = Symbol(
        _mapping_value(raw, "ocean_initial_conditions", "analytic"),
    )
    era5_month = Int(_mapping_value(raw, "era5_month", 1))
    ecco_year = Int(_mapping_value(raw, "ecco_year", 1993))
    ecco_month = Int(_mapping_value(raw, "ecco_month", 1))
    default_start_date = if atmosphere_initial_conditions in (
        :era5_monthly,
        :era5_instantaneous,
    )
        DateTime(1993, era5_month, 1)
    elseif ocean_initial_conditions in (:ecco_1993, :ecco_1993_full_state)
        # The pinned ECCO initializer below loads the January 1993 state. Keep
        # the model calendar aligned even for analytic-atmosphere experiments
        # that omit an explicit start_date.
        DateTime(1993, 1, 1)
    elseif ocean_initial_conditions in (:ecco_monthly, :ecco_monthly_full_state)
        DateTime(ecco_year, ecco_month, 1)
    else
        DateTime(2000, 1, 1)
    end

    config = ExperimentConfig(
        name = String(_mapping_value(raw, "name", "baseline")),
        device = Symbol(_mapping_value(raw, "device", "cpu")),
        truncation = Int(_mapping_value(raw, "truncation", 15)),
        nlayers = Int(_mapping_value(raw, "nlayers", 8)),
        duration_days = Float64(_mapping_value(raw, "duration_days", 10.0)),
        analysis_start_day = Float64(
            _mapping_value(raw, "analysis_start_day", 0.0),
        ),
        start_date = _parse_start_date(
            _mapping_value(raw, "start_date", default_start_date),
        ),
        atmosphere_timestep_at_t31_minutes = Float64(
            _mapping_value(raw, "atmosphere_timestep_at_t31_minutes", 40.0),
        ),
        atmosphere_hyperdiffusion_hours = Float64(
            _mapping_value(raw, "atmosphere_hyperdiffusion_hours", 4.0),
        ),
        atmosphere_divergence_hyperdiffusion_hours = Float64(
            _mapping_value(
                raw,
                "atmosphere_divergence_hyperdiffusion_hours",
                1.0,
            ),
        ),
        atmosphere_vertical_diffusion = Symbol(
            _mapping_value(
                raw,
                "atmosphere_vertical_diffusion",
                "metric_bulk_richardson",
            ),
        ),
        atmosphere_large_scale_precipitation = Symbol(
            _mapping_value(
                raw,
                "atmosphere_large_scale_precipitation",
                "upstream_observed",
            ),
        ),
        atmosphere_convection = Symbol(
            _mapping_value(
                raw,
                "atmosphere_convection",
                "betts_miller_constant_rh",
            ),
        ),
        atmosphere_convection_upper_relative_humidity = Float64(
            _mapping_value(
                raw,
                "atmosphere_convection_upper_relative_humidity",
                0.5,
            ),
        ),
        atmosphere_convection_lower_relative_humidity = Float64(
            _mapping_value(
                raw,
                "atmosphere_convection_lower_relative_humidity",
                0.7,
            ),
        ),
        atmosphere_convection_transition_top_sigma = Float64(
            _mapping_value(
                raw,
                "atmosphere_convection_transition_top_sigma",
                0.3,
            ),
        ),
        atmosphere_convection_transition_bottom_sigma = Float64(
            _mapping_value(
                raw,
                "atmosphere_convection_transition_bottom_sigma",
                0.7,
            ),
        ),
        atmosphere_speed_limit_ms = Float64(
            _mapping_value(raw, "atmosphere_speed_limit_ms", 80.0),
        ),
        atmosphere_speed_limit_drag_per_m = Float64(
            _mapping_value(raw, "atmosphere_speed_limit_drag_per_m", 4.0e-7),
        ),
        atmosphere_surface_speed_limit_ms = Float64(
            _mapping_value(raw, "atmosphere_surface_speed_limit_ms", 80.0),
        ),
        atmosphere_surface_speed_limit_drag_per_m = Float64(
            _mapping_value(
                raw,
                "atmosphere_surface_speed_limit_drag_per_m",
                4.0e-7,
            ),
        ),
        mixed_layer_depth_m = Float64(_mapping_value(raw, "mixed_layer_depth_m", 50.0)),
        ocean_model = Symbol(_mapping_value(raw, "ocean_model", "slab")),
        ocean_grid = Symbol(_mapping_value(raw, "ocean_grid", "idealized")),
        ocean_dynamics = Bool(_mapping_value(raw, "ocean_dynamics", false)),
        sea_ice_dynamics = Bool(_mapping_value(raw, "sea_ice_dynamics", false)),
        sea_ice_momentum_substeps = Int(
            _mapping_value(raw, "sea_ice_momentum_substeps", 240),
        ),
        sea_ice_pressure_formulation = Symbol(
            _mapping_value(
                raw,
                "sea_ice_pressure_formulation",
                "replacement_pressure",
            ),
        ),
        sea_ice_immersed_boundary_drag_coefficient = Float64(
            _mapping_value(
                raw,
                "sea_ice_immersed_boundary_drag_coefficient",
                0.0,
            ),
        ),
        ocean_nlongitude = Int(_mapping_value(raw, "ocean_nlongitude", 32)),
        ocean_nlatitude = Int(_mapping_value(raw, "ocean_nlatitude", 16)),
        ocean_nlayers = Int(_mapping_value(raw, "ocean_nlayers", 4)),
        ocean_depth_m = Float64(_mapping_value(raw, "ocean_depth_m", 1000.0)),
        ocean_bathymetry_correction = Symbol(
            _mapping_value(raw, "ocean_bathymetry_correction", "none"),
        ),
        ocean_tripolar_wet_mask = Symbol(
            _mapping_value(
                raw,
                "ocean_tripolar_wet_mask",
                "interpolated_mean_elevation",
            ),
        ),
        ocean_shortwave_scheme = Symbol(
            _mapping_value(raw, "ocean_shortwave_scheme", "two_color_default"),
        ),
        ocean_river_mouth_mixing = Symbol(
            _mapping_value(raw, "ocean_river_mouth_mixing", "none"),
        ),
        ocean_river_mouth_vertical_diffusivity_m2s = Float64(
            _mapping_value(
                raw,
                "ocean_river_mouth_vertical_diffusivity_m2s",
                1.0e-2,
            ),
        ),
        ocean_river_mouth_horizontal_diffusivity_m2s = Float64(
            _mapping_value(
                raw,
                "ocean_river_mouth_horizontal_diffusivity_m2s",
                4.0e3,
            ),
        ),
        ocean_river_mouth_mixing_depth_m = Float64(
            _mapping_value(raw, "ocean_river_mouth_mixing_depth_m", 50.0),
        ),
        ocean_river_mouth_reference_freshwater_mass_flux_kgm2s = Float64(
            _mapping_value(
                raw,
                "ocean_river_mouth_reference_freshwater_mass_flux_kgm2s",
                0.0,
            ),
        ),
        ocean_polar_sponge_start_latitude_degrees = Float64(
            _mapping_value(
                raw,
                "ocean_polar_sponge_start_latitude_degrees",
                70.0,
            ),
        ),
        ocean_polar_sponge_stop_latitude_degrees = Float64(
            _mapping_value(
                raw,
                "ocean_polar_sponge_stop_latitude_degrees",
                75.0,
            ),
        ),
        ocean_polar_sponge_timescale_hours = Float64(
            _mapping_value(raw, "ocean_polar_sponge_timescale_hours", 3.0),
        ),
        ocean_initial_conditions = ocean_initial_conditions,
        ecco_year = ecco_year,
        ecco_month = ecco_month,
        ecco_initial_conditions_directory = _resolve_optional_path(
            _mapping_value(raw, "ecco_initial_conditions_directory", ""),
        ),
        atmosphere_initial_conditions = atmosphere_initial_conditions,
        era5_pressure_levels_path = _resolve_optional_path(
            _mapping_value(raw, "era5_pressure_levels_path", ""),
        ),
        era5_single_levels_path = _resolve_optional_path(
            _mapping_value(raw, "era5_single_levels_path", ""),
        ),
        era5_land_state_path = _resolve_optional_path(
            _mapping_value(raw, "era5_land_state_path", ""),
        ),
        era5_month = era5_month,
        land_model = Symbol(_mapping_value(raw, "land_model", "speedy_bucket")),
        terrarium_runoff_routing = Symbol(_mapping_value(
            raw,
            "terrarium_runoff_routing",
            "equal_area_flux",
        )),
        terrarium_initial_conditions = Symbol(
            _mapping_value(raw, "terrarium_initial_conditions", "hydrostatic"),
        ),
        terrarium_evapotranspiration = Symbol(
            _mapping_value(raw, "terrarium_evapotranspiration", "bare_ground"),
        ),
        terrarium_max_leaf_conductance_ms = Float64(
            _mapping_value(raw, "terrarium_max_leaf_conductance_ms", 2.0e-3),
        ),
        terrarium_min_leaf_conductance_ms = Float64(
            _mapping_value(raw, "terrarium_min_leaf_conductance_ms", 5.0e-5),
        ),
        terrarium_light_half_saturation_wm2 = Float64(
            _mapping_value(raw, "terrarium_light_half_saturation_wm2", 100.0),
        ),
        terrarium_vpd_scale_pa = Float64(
            _mapping_value(raw, "terrarium_vpd_scale_pa", 2000.0),
        ),
        terrarium_canopy_extinction_coefficient = Float64(
            _mapping_value(raw, "terrarium_canopy_extinction_coefficient", 0.5),
        ),
        terrarium_root_a_m1 = Float64(
            _mapping_value(raw, "terrarium_root_a_m1", 7.0),
        ),
        terrarium_root_b_m1 = Float64(
            _mapping_value(raw, "terrarium_root_b_m1", 2.0),
        ),
        terrarium_soil_layers = Int(_mapping_value(raw, "terrarium_soil_layers", 8)),
        terrarium_max_layer_thickness_m = Float64(
            _mapping_value(raw, "terrarium_max_layer_thickness_m", 5.0),
        ),
        terrarium_timestep_seconds = Float64(
            _mapping_value(raw, "terrarium_timestep_seconds", 60.0),
        ),
        output_dir = output_dir,
        forcing = forcing,
    )
    validate(config)
    return config
end

function validate(config::ExperimentConfig)
    _validate_finite_float_fields(config, "experiment configuration")
    _validate_finite_float_fields(config.forcing, "forcing configuration")
    config.device in (:cpu, :gpu) || throw(ArgumentError("device must be cpu or gpu"))
    config.truncation >= 7 || throw(ArgumentError("truncation must be at least 7"))
    config.nlayers >= 2 || throw(ArgumentError("nlayers must be at least 2"))
    config.duration_days > 0 || throw(ArgumentError("duration_days must be positive"))
    0 <= config.analysis_start_day < config.duration_days || throw(
        ArgumentError("analysis_start_day must satisfy 0 <= start < duration_days"),
    )
    config.atmosphere_timestep_at_t31_minutes > 0 ||
        throw(ArgumentError("atmosphere_timestep_at_t31_minutes must be positive"))
    config.atmosphere_hyperdiffusion_hours > 0 || throw(
        ArgumentError("atmosphere_hyperdiffusion_hours must be positive"),
    )
    config.atmosphere_divergence_hyperdiffusion_hours > 0 || throw(
        ArgumentError(
            "atmosphere_divergence_hyperdiffusion_hours must be positive",
        ),
    )
    config.atmosphere_vertical_diffusion in (
        :metric_bulk_richardson,
        :speedyweather_0_21_1_zero_operator_control,
    ) || throw(ArgumentError(
        "atmosphere_vertical_diffusion must be metric_bulk_richardson or " *
        "speedyweather_0_21_1_zero_operator_control",
    ))
    config.atmosphere_large_scale_precipitation in (
        :upstream_observed,
        :net_column_conservative,
    ) || throw(ArgumentError(
        "atmosphere_large_scale_precipitation must be upstream_observed or " *
        "net_column_conservative",
    ))
    config.atmosphere_convection in (
        :betts_miller_constant_rh,
        :betts_miller_sigma_rh_v1,
    ) || throw(ArgumentError(
        "atmosphere_convection must be betts_miller_constant_rh or " *
        "betts_miller_sigma_rh_v1",
    ))
    0 <= config.atmosphere_convection_upper_relative_humidity <= 1 || throw(
        ArgumentError(
            "atmosphere_convection_upper_relative_humidity must lie in [0, 1]",
        ),
    )
    0 <= config.atmosphere_convection_lower_relative_humidity <= 1 || throw(
        ArgumentError(
            "atmosphere_convection_lower_relative_humidity must lie in [0, 1]",
        ),
    )
    0 <= config.atmosphere_convection_transition_top_sigma <
        config.atmosphere_convection_transition_bottom_sigma <= 1 || throw(
        ArgumentError(
            "atmosphere convection transition sigmas must satisfy " *
            "0 <= top < bottom <= 1",
        ),
    )
    config.atmosphere_speed_limit_ms > 0 || throw(
        ArgumentError("atmosphere_speed_limit_ms must be positive"),
    )
    config.atmosphere_speed_limit_drag_per_m > 0 || throw(
        ArgumentError("atmosphere_speed_limit_drag_per_m must be positive"),
    )
    config.atmosphere_surface_speed_limit_ms > 0 || throw(
        ArgumentError("atmosphere_surface_speed_limit_ms must be positive"),
    )
    config.atmosphere_surface_speed_limit_drag_per_m > 0 || throw(
        ArgumentError(
            "atmosphere_surface_speed_limit_drag_per_m must be positive",
        ),
    )
    config.mixed_layer_depth_m > 0 || throw(ArgumentError("mixed_layer_depth_m must be positive"))
    config.ocean_model in (:slab, :oceananigans) ||
        throw(ArgumentError("ocean_model must be slab or oceananigans"))
    config.ocean_grid in (
        :idealized,
        :latitude_longitude_1degree,
        :tripolar_1degree,
    ) || throw(ArgumentError(
        "ocean_grid must be idealized, latitude_longitude_1degree, or tripolar_1degree",
    ))
    config.ocean_nlongitude >= 8 || throw(ArgumentError("ocean_nlongitude must be at least 8"))
    config.ocean_nlatitude >= 4 || throw(ArgumentError("ocean_nlatitude must be at least 4"))
    config.ocean_nlayers >= 1 || throw(ArgumentError("ocean_nlayers must be positive"))
    if config.ocean_grid == :latitude_longitude_1degree
        (config.ocean_nlongitude, config.ocean_nlatitude) == (360, 150) ||
            throw(ArgumentError(
                "latitude_longitude_1degree requires ocean dimensions 360x150",
            ))
    end
    if config.ocean_grid == :tripolar_1degree
        config.ocean_nlongitude % 4 == 0 || throw(
            ArgumentError("tripolar ocean_nlongitude must be divisible by 4"),
        )
        config.ocean_nlatitude >= 8 || throw(
            ArgumentError("tripolar ocean_nlatitude must be at least 8"),
        )
        config.ocean_nlayers >= 4 || throw(
            ArgumentError("tripolar ocean_nlayers must be at least 4"),
        )
    end
    config.ocean_depth_m > 0 || throw(ArgumentError("ocean_depth_m must be positive"))
    config.ocean_bathymetry_correction in (
        :none,
        :persian_gulf_1degree_v1,
        :shallow_culdesacs_1degree_v1,
        :shallow_culdesacs_hormuz_1degree_v2,
        :shallow_culdesacs_polar_boundary_1degree_v3,
        :shallow_culdesacs_global_40m_1degree_v4,
    ) || throw(ArgumentError(
        "ocean_bathymetry_correction must be none, " *
        "persian_gulf_1degree_v1, shallow_culdesacs_1degree_v1, " *
        "shallow_culdesacs_hormuz_1degree_v2, or " *
        "shallow_culdesacs_polar_boundary_1degree_v3, or " *
        "shallow_culdesacs_global_40m_1degree_v4",
    ))
    if config.ocean_bathymetry_correction != :none &&
       config.ocean_grid != :latitude_longitude_1degree
        throw(ArgumentError(
            "the configured ocean bathymetry correction requires " *
            "ocean_grid=latitude_longitude_1degree",
        ))
    end
    config.ocean_tripolar_wet_mask in (
        :interpolated_mean_elevation,
        :etopo_majority_area_polar,
    ) || throw(ArgumentError(
        "ocean_tripolar_wet_mask must be interpolated_mean_elevation or " *
        "etopo_majority_area_polar",
    ))
    if config.ocean_tripolar_wet_mask != :interpolated_mean_elevation
        config.ocean_grid == :tripolar_1degree || throw(ArgumentError(
            "ocean_tripolar_wet_mask=etopo_majority_area_polar requires " *
            "ocean_grid=tripolar_1degree",
        ))
        (config.ocean_nlongitude, config.ocean_nlatitude) == (360, 180) ||
            throw(ArgumentError(
                "ocean_tripolar_wet_mask=etopo_majority_area_polar is " *
                "qualified only for the 360x180 tripolar grid",
            ))
    end
    config.ocean_shortwave_scheme in (
        :two_color_default,
        :depth_aware_coastal,
    ) || throw(ArgumentError(
        "ocean_shortwave_scheme must be two_color_default or depth_aware_coastal",
    ))
    (config.ocean_river_mouth_mixing == :none ||
     _is_localized_estuary_mixing(config.ocean_river_mouth_mixing)) ||
        throw(ArgumentError(
            "ocean_river_mouth_mixing must be none, localized_estuary_v1, " *
            "localized_estuary_v2, localized_estuary_v3, localized_estuary_v4, " *
            "or localized_estuary_v5",
        ))
    config.ocean_river_mouth_vertical_diffusivity_m2s > 0 || throw(
        ArgumentError("ocean river-mouth vertical diffusivity must be positive"),
    )
    config.ocean_river_mouth_horizontal_diffusivity_m2s > 0 || throw(
        ArgumentError("ocean river-mouth horizontal diffusivity must be positive"),
    )
    config.ocean_river_mouth_mixing_depth_m > 0 || throw(
        ArgumentError("ocean river-mouth mixing depth must be positive"),
    )
    reference_freshwater_flux =
        config.ocean_river_mouth_reference_freshwater_mass_flux_kgm2s
    isfinite(reference_freshwater_flux) && reference_freshwater_flux >= 0 ||
        throw(ArgumentError(
            "ocean river-mouth reference freshwater mass flux must be finite " *
            "and nonnegative",
        ))
    if config.ocean_river_mouth_mixing == :localized_estuary_v5
        reference_freshwater_flux > 0 || throw(ArgumentError(
            "localized_estuary_v5 requires a positive ocean river-mouth " *
            "reference freshwater mass flux",
        ))
    else
        iszero(reference_freshwater_flux) || throw(ArgumentError(
            "ocean river-mouth reference freshwater mass flux is only active " *
            "for localized_estuary_v5",
        ))
    end
    if config.ocean_river_mouth_mixing != :none
        config.ocean_grid == :tripolar_1degree || throw(ArgumentError(
            "localized river-mouth mixing currently requires " *
            "ocean_grid=tripolar_1degree",
        ))
        config.ocean_dynamics || throw(ArgumentError(
            "localized river-mouth mixing requires ocean_dynamics=true",
        ))
        config.land_model == :terrarium || throw(ArgumentError(
            "localized river-mouth mixing requires land_model=terrarium",
        ))
    end
    0 <= config.ocean_polar_sponge_start_latitude_degrees <
         config.ocean_polar_sponge_stop_latitude_degrees <= 90 || throw(
        ArgumentError(
            "ocean polar sponge latitudes must satisfy 0 <= start < stop <= 90 degrees",
        ),
    )
    config.ocean_polar_sponge_timescale_hours > 0 || throw(
        ArgumentError("ocean_polar_sponge_timescale_hours must be positive"),
    )
    (config.ocean_initial_conditions == :analytic ||
     _is_ecco_initial_conditions(config.ocean_initial_conditions)) || throw(
        ArgumentError(
            "ocean_initial_conditions must be analytic, ecco_1993, " *
            "ecco_1993_full_state, ecco_monthly, or ecco_monthly_full_state",
        ),
    )
    1992 <= config.ecco_year <= 2017 || throw(ArgumentError(
        "ecco_year must be in 1992:2017",
    ))
    1 <= config.ecco_month <= 12 || throw(ArgumentError(
        "ecco_month must be in 1:12",
    ))
    if _is_ecco_initial_conditions(config.ocean_initial_conditions)
        config.ocean_model == :oceananigans || throw(ArgumentError(
            "ECCO ocean/ice initialization requires ocean_model=oceananigans",
        ))
        if config.ocean_initial_conditions in (:ecco_1993, :ecco_1993_full_state)
            (config.ecco_year, config.ecco_month) == (1993, 1) || throw(
                ArgumentError(
                    "legacy ECCO 1993 initialization requires ecco_year=1993 " *
                    "and ecco_month=1",
                ),
            )
            config.start_date == DateTime(1993, 1, 1) || throw(ArgumentError(
                "ECCO 1993-01 ocean/ice initialization requires " *
                "start_date=1993-01-01T00:00:00",
            ))
        else
            (year(config.start_date), month(config.start_date)) ==
                (config.ecco_year, config.ecco_month) || throw(ArgumentError(
                    "generic ECCO monthly initialization requires the model " *
                    "start_date to fall within the configured source month",
                ))
            isempty(config.ecco_initial_conditions_directory) && throw(
                ArgumentError(
                    "generic ECCO monthly initialization requires " *
                    "ecco_initial_conditions_directory",
                ),
            )
        end

        if !isempty(config.ecco_initial_conditions_directory)
            isdir(config.ecco_initial_conditions_directory) || throw(ArgumentError(
                "ECCO initial-condition directory not found: " *
                config.ecco_initial_conditions_directory,
            ))
            date_label = Dates.format(
                _ecco_initial_condition_date(config),
                dateformat"yyyy_mm",
            )
            variables = _is_full_state_ecco_initial_conditions(
                config.ocean_initial_conditions,
            ) ? ("THETA", "SALT", "EVEL", "NVEL", "SSH", "SIheff", "SIarea") :
                ("THETA", "SALT", "SIheff", "SIarea")
            missing_paths = String[]
            for variable in variables
                path = joinpath(
                    config.ecco_initial_conditions_directory,
                    "$(variable)_$(date_label).nc",
                )
                isfile(path) || push!(missing_paths, path)
            end
            isempty(missing_paths) || throw(ArgumentError(
                "ECCO initial-condition directory lacks required file(s): " *
                join(missing_paths, ", "),
            ))
        end
    elseif !isempty(config.ecco_initial_conditions_directory)
        throw(ArgumentError(
            "ecco_initial_conditions_directory requires an ECCO ocean initializer",
        ))
    end
    config.atmosphere_initial_conditions in (
        :analytic,
        :era5_monthly,
        :era5_instantaneous,
    ) || throw(
        ArgumentError(
            "atmosphere_initial_conditions must be analytic, era5_monthly, " *
            "or era5_instantaneous",
        ),
    )
    1 <= config.era5_month <= 12 || throw(ArgumentError("era5_month must be in 1:12"))
    if config.atmosphere_initial_conditions in (
        :era5_monthly,
        :era5_instantaneous,
    )
        isempty(config.era5_pressure_levels_path) && throw(
            ArgumentError("ERA5 initialization requires era5_pressure_levels_path"),
        )
        isempty(config.era5_single_levels_path) && throw(
            ArgumentError("ERA5 initialization requires era5_single_levels_path"),
        )
        isfile(config.era5_pressure_levels_path) || throw(
            ArgumentError("ERA5 pressure-level file not found: $(config.era5_pressure_levels_path)"),
        )
        isfile(config.era5_single_levels_path) || throw(
            ArgumentError("ERA5 single-level file not found: $(config.era5_single_levels_path)"),
        )
    end
    config.land_model in (:speedy_bucket, :terrarium) ||
        throw(ArgumentError("land_model must be speedy_bucket or terrarium"))
    config.terrarium_runoff_routing in (
        :equal_area_flux,
        :equal_column_fraction,
    ) || throw(ArgumentError(
        "terrarium_runoff_routing must be equal_area_flux or " *
        "equal_column_fraction",
    ))
    config.terrarium_initial_conditions in (
        :hydrostatic,
        :era5_instantaneous,
    ) || throw(ArgumentError(
        "terrarium_initial_conditions must be hydrostatic or era5_instantaneous",
    ))
    if config.terrarium_initial_conditions == :era5_instantaneous
        config.land_model == :terrarium || throw(ArgumentError(
            "terrarium_initial_conditions=era5_instantaneous requires land_model=terrarium",
        ))
        isempty(config.era5_land_state_path) && throw(ArgumentError(
            "ERA5 Terrarium initialization requires era5_land_state_path",
        ))
        isfile(config.era5_land_state_path) || throw(ArgumentError(
            "ERA5 land-state file not found: $(config.era5_land_state_path)",
        ))
    end
    config.terrarium_evapotranspiration in (
        :bare_ground,
        :era5_prescribed_vegetation,
    ) || throw(ArgumentError(
        "terrarium_evapotranspiration must be bare_ground or " *
        "era5_prescribed_vegetation",
    ))
    if config.terrarium_evapotranspiration == :era5_prescribed_vegetation
        config.land_model == :terrarium || throw(ArgumentError(
            "ERA5 prescribed vegetation requires land_model=terrarium",
        ))
        config.terrarium_initial_conditions == :era5_instantaneous || throw(
            ArgumentError(
                "ERA5 prescribed vegetation requires " *
                "terrarium_initial_conditions=era5_instantaneous",
            ),
        )
    end
    0 < config.terrarium_min_leaf_conductance_ms <=
        config.terrarium_max_leaf_conductance_ms || throw(ArgumentError(
            "Terrarium leaf conductances must satisfy 0 < minimum <= maximum",
        ))
    config.terrarium_light_half_saturation_wm2 > 0 || throw(ArgumentError(
        "Terrarium light half-saturation must be positive",
    ))
    config.terrarium_vpd_scale_pa > 0 || throw(ArgumentError(
        "Terrarium VPD scale must be positive",
    ))
    config.terrarium_canopy_extinction_coefficient > 0 || throw(ArgumentError(
        "Terrarium canopy extinction coefficient must be positive",
    ))
    config.terrarium_root_a_m1 > 0 || throw(ArgumentError(
        "Terrarium root a parameter must be positive",
    ))
    config.terrarium_root_b_m1 > 0 || throw(ArgumentError(
        "Terrarium root b parameter must be positive",
    ))
    config.terrarium_soil_layers >= 2 ||
        throw(ArgumentError("terrarium_soil_layers must be at least 2"))
    config.terrarium_max_layer_thickness_m >= 0.05 ||
        throw(ArgumentError("terrarium_max_layer_thickness_m must be at least 0.05"))
    config.terrarium_timestep_seconds > 0 ||
        throw(ArgumentError("terrarium_timestep_seconds must be positive"))
    config.forcing.radiation in (
        :speedy_simplified,
        :rrtmgp_clear_sky,
        :rrtmgp_all_sky,
    ) || throw(ArgumentError(
        "radiation must be speedy_simplified, rrtmgp_clear_sky, or rrtmgp_all_sky",
    ))
    config.forcing.co2_ppm > 0 || throw(ArgumentError("co2_ppm must be positive"))
    config.forcing.aerosol_optical_depth_550nm >= 0 ||
        throw(ArgumentError("aerosol optical depth must be non-negative"))
    config.forcing.aerosol_species == "sulfate" ||
        throw(ArgumentError("the current prescribed-AOD adapter supports aerosol_species=sulfate"))
    rrtmgp_active = config.forcing.radiation in (
        :rrtmgp_clear_sky,
        :rrtmgp_all_sky,
    )
    if !rrtmgp_active && config.forcing.co2_ppm != 420
        throw(ArgumentError(
            "CO2 perturbations require radiation=rrtmgp_clear_sky or rrtmgp_all_sky; " *
            "SpeedyWeather's simplified radiation does not consume co2_ppm",
        ))
    end
    if !rrtmgp_active && config.forcing.aerosol_optical_depth_550nm > 0
        throw(ArgumentError(
            "prescribed aerosol AOD requires radiation=rrtmgp_clear_sky or rrtmgp_all_sky",
        ))
    end
    config.forcing.radiation_every_n_steps > 0 ||
        throw(ArgumentError("radiation_every_n_steps must be positive"))
    config.forcing.cloud_scheme in (
        :none,
        :diagnostic,
        :prognostic_condensate,
    ) || throw(ArgumentError(
        "cloud_scheme must be none, diagnostic, or prognostic_condensate",
    ))
    0 <= config.forcing.cloud_humidity_search_min_sigma <= 1 ||
        throw(ArgumentError("cloud humidity-search minimum sigma must be within [0, 1]"))
    if config.forcing.cloud_scheme == :none &&
       config.forcing.cloud_humidity_search_min_sigma != 0
        throw(ArgumentError(
            "cloud_humidity_search_min_sigma is active only with " *
            "an all-sky cloud scheme",
        ))
    end
    config.forcing.cloud_liquid_water_path_gm2 >= 0 ||
        throw(ArgumentError("cloud liquid water path must be non-negative"))
    config.forcing.cloud_ice_water_path_gm2 >= 0 ||
        throw(ArgumentError("cloud ice water path must be non-negative"))
    0 <= config.forcing.cloud_condensate_retention_fraction <= 1 ||
        throw(ArgumentError(
            "cloud condensate retention fraction must be within [0, 1]",
        ))
    config.forcing.cloud_condensate_residence_time_hours > 0 ||
        throw(ArgumentError(
            "cloud condensate residence time must be positive",
        ))
    if config.forcing.cloud_scheme == :prognostic_condensate
        config.forcing.cloud_condensate_retention_fraction > 0 ||
            throw(ArgumentError(
                "prognostic_condensate requires a positive condensate retention fraction",
            ))
    elseif config.forcing.cloud_condensate_retention_fraction != 0
        throw(ArgumentError(
            "cloud_condensate_retention_fraction is active only with " *
            "cloud_scheme=prognostic_condensate",
        ))
    elseif config.forcing.cloud_condensate_residence_time_hours != 6
        throw(ArgumentError(
            "cloud_condensate_residence_time_hours is active only with " *
            "cloud_scheme=prognostic_condensate",
        ))
    end
    if config.forcing.cloud_scheme != :none && config.forcing.radiation != :rrtmgp_all_sky
        throw(ArgumentError("all-sky cloud schemes require radiation=rrtmgp_all_sky"))
    end
    if config.ocean_grid == :latitude_longitude_1degree && config.ocean_model != :oceananigans
        throw(ArgumentError("ocean_grid=latitude_longitude_1degree requires ocean_model=oceananigans"))
    end
    if config.ocean_grid == :tripolar_1degree && config.ocean_model != :oceananigans
        throw(ArgumentError("ocean_grid=tripolar_1degree requires ocean_model=oceananigans"))
    end
    if config.ocean_dynamics && config.ocean_model != :oceananigans
        throw(ArgumentError("ocean_dynamics=true requires ocean_model=oceananigans"))
    end
    if config.sea_ice_dynamics && !config.ocean_dynamics
        throw(ArgumentError("sea_ice_dynamics=true requires ocean_dynamics=true"))
    end
    config.sea_ice_momentum_substeps >= 1 || throw(ArgumentError(
        "sea_ice_momentum_substeps must be positive",
    ))
    config.sea_ice_pressure_formulation in (
        :replacement_pressure,
        :ice_strength,
    ) || throw(ArgumentError(
        "sea_ice_pressure_formulation must be replacement_pressure or " *
        "ice_strength",
    ))
    config.sea_ice_immersed_boundary_drag_coefficient >= 0 || throw(
        ArgumentError(
            "sea_ice_immersed_boundary_drag_coefficient must be non-negative",
        ),
    )
    return config
end
