include(joinpath(@__DIR__, "c08_finite", "C08FiniteConvection.jl"))
const FiniteC08Core = C08FiniteConvection.C08MultiPlume

const FINITE_MULTIPLUME_STATUS_NAMES = (
    :inactive, :no_stabilising_response, :partial_adjustment, :target_crossing,
    :valid_crossing_with_invalid_midpoint, :admissibility_limited,
    :entrainment_limited, :no_valid_stabilising_response, :rounding_zero,
)
const FINITE_MULTIPLUME_STATUS_CODES =
    Dict(name => i for (i, name) in enumerate(FINITE_MULTIPLUME_STATUS_NAMES))
const FINITE_MULTIPLUME_FIELDS = (
    :finite_c08_status, :finite_c08_evaluations, :finite_c08_invalid_probes,
    :finite_c08_target_cape, :finite_c08_final_cape, :finite_c08_initial_cape,
)

"""
Finite-response three-plume convection for a Float32 CPU atmosphere.

The closure evaluates profiles and exchanges in Float64, then returns Float32
physical tendencies and precipitation. Detrained liquid and ice re-evaporate
locally with their phase energy. The atmosphere carries no convective condensate
memory. GPU execution requires a separate device-workspace implementation.
"""
struct FiniteMultiPlumeConvection{S} <: SpeedyWeather.AbstractConvection
    spectrum::S
    saturation_law::Symbol
end

function _require_finite_multiplume_cpu(grid)
    grid.architecture isa SpeedyWeather.Architectures.AbstractCPU ||
        throw(ArgumentError("finite_multiplume_c08_v1 requires a CPU atmosphere"))
    grid.NF === Float32 ||
        throw(ArgumentError("finite_multiplume_c08_v1 requires Float32 atmosphere storage"))
    return nothing
end

function FiniteMultiPlumeConvection(
    grid::SpeedyWeather.SpectralGrid;
    saturation_law::Symbol = :mixed_phase,
    time_scale_seconds::Real = 7200,
    entrainment_scale::Real = 0.75,
    precipitation_conversion_m_inv::Real = 0.0015,
    condensate_threshold_kg_kg::Real = 0.0005,
)
    _require_finite_multiplume_cpu(grid)
    saturation_law in (:liquid, :mixed_phase) ||
        throw(ArgumentError("C08 saturation law must be liquid or mixed_phase"))
    isfinite(time_scale_seconds) && time_scale_seconds > 0 ||
        throw(ArgumentError("C08 relaxation time must be positive and finite"))
    for (name, value) in (
        (:entrainment_scale, entrainment_scale),
        (:precipitation_conversion_m_inv, precipitation_conversion_m_inv),
        (:condensate_threshold_kg_kg, condensate_threshold_kg_kg),
    )
        isfinite(value) && value >= 0 ||
            throw(ArgumentError("C08 $name must be non-negative and finite"))
    end
    spectrum = FiniteC08Core.typed_plume_spectrum(Float32;
        entrainment_m = entrainment_scale .* (0.3e-3, 1.0e-3, 3.0e-3),
        weights = (1/3, 1/3, 1/3),
        closure_iterations = 8,
        max_layer_entrainment = 2.0,
        time_scale_s = time_scale_seconds,
        precip_conversion_m = precipitation_conversion_m_inv,
        condensate_threshold = condensate_threshold_kg_kg)
    return FiniteMultiPlumeConvection(spectrum, saturation_law)
end

SpeedyWeather.initialize!(
    ::FiniteMultiPlumeConvection, model::SpeedyWeather.AbstractModel,
) = _require_finite_multiplume_cpu(model.spectral_grid)

function _finite_multiplume_metadata(config::ExperimentConfig)
    config.atmosphere_convection in (:finite_multiplume_c08_v1, :finite_multiplume_c07_v1) || return (;)
    return (
        atmosphere_convection_closure = "finite_host_state_CAPE_relaxation",
        atmosphere_convection_solver_precision = "Float64",
        atmosphere_convection_storage_precision = "Float32",
        atmosphere_convection_timescale_seconds = config.atmosphere_convection_timescale_seconds,
        atmosphere_c08_entrainment_scale = config.atmosphere_c08_entrainment_scale,
        atmosphere_c08_precipitation_conversion_m_inv = config.atmosphere_c08_precipitation_conversion_m_inv,
        atmosphere_c08_condensate_threshold_kg_kg = config.atmosphere_c08_condensate_threshold_kg_kg,
        atmosphere_convection_condensate_policy = "local_reevaporation_with_phase_energy",
        atmosphere_convection_gpu_supported = 0,
    )
end

SpeedyWeather.variables(::FiniteMultiPlumeConvection) = (
    SpeedyWeather.ParameterizationVariable(
        :rain_convection,
        SpeedyWeather.Grid2D();
        desc = "Convective rain accumulation",
        units = "m",
    ),
    SpeedyWeather.ParameterizationVariable(
        :rain_rate_convection,
        SpeedyWeather.Grid2D();
        desc = "Convective rain rate",
        units = "m/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :rain_rate,
        SpeedyWeather.Grid2D();
        desc = "Total rain rate",
        units = "m/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :snow_rate,
        SpeedyWeather.Grid2D();
        desc = "Total snow rate",
        units = "m/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :cloud_top,
        SpeedyWeather.Grid2D();
        desc = "Cloud top layer index",
        units = "1",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_humidity_tendency,
        SpeedyWeather.Grid3D();
        desc = "C08 host specific-humidity tendency",
        units = "kg/kg/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_temperature_tendency,
        SpeedyWeather.Grid3D();
        desc = "C08 host temperature tendency",
        units = "K/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_precipitation_surface_pressure,
        SpeedyWeather.Grid2D();
        desc = "Surface pressure used by C08",
        units = "Pa",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_precipitation_mass_correction,
        SpeedyWeather.Grid2D();
        desc = "C08 precipitation mass correction, identically zero",
        units = "m/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_triggered,
        SpeedyWeather.Grid2D();
        desc = "One when the C08 spectrum triggers",
        units = "1",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_active_plumes,
        SpeedyWeather.Grid2D();
        desc = "Number of active plume bins",
        units = "1",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_mass_flux,
        SpeedyWeather.Grid2D();
        desc = "Composite cloud-base mass flux",
        units = "kg/m2/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_weighted_cape,
        SpeedyWeather.Grid2D();
        desc = "Active-bin weighted dilute CAPE",
        units = "J/kg",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_maximum_cin,
        SpeedyWeather.Grid2D();
        desc = "Maximum active-bin dilute CIN",
        units = "J/kg",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_response,
        SpeedyWeather.Grid2D();
        desc = "Finite-step CAPE response",
        units = "J/kg/(kg/m2)",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_limited,
        SpeedyWeather.Grid2D();
        desc = "One when the entrainment limiter is active",
        units = "1",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_degenerate,
        SpeedyWeather.Grid2D();
        desc = "One when mass flux falls below the rounding guard",
        units = "1",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_water_residual,
        SpeedyWeather.Grid2D();
        desc = "C08 column water closure residual",
        units = "kg/m2/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_water_gross_exchange,
        SpeedyWeather.Grid2D();
        desc = "C08 gross column water exchange",
        units = "kg/m2/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_enthalpy_residual,
        SpeedyWeather.Grid2D();
        desc = "C08 frozen moist-enthalpy closure residual",
        units = "W/m2",
    ),
    SpeedyWeather.ParameterizationVariable(
        :c08_enthalpy_gross_exchange,
        SpeedyWeather.Grid2D();
        desc = "C08 gross frozen moist-enthalpy exchange",
        units = "W/m2",
    ),
    (SpeedyWeather.ParameterizationVariable(
        name, SpeedyWeather.Grid2D();
        desc = "Finite C08 closure diagnostic",
        units = occursin("cape", String(name)) ? "J/kg" : "1")
     for name in FINITE_MULTIPLUME_FIELDS)...,
    SpeedyWeather.ParameterizationVariable(
        :snow_rate_convection, SpeedyWeather.Grid2D();
        desc = "Snow increment across convection", units = "m/s"),
)

function _finite_multiplume_thermo(scheme, model, ::Type{NF}) where {NF}
    atmosphere = model.atmosphere
    return FiniteC08Core.typed_thermo(
        NF;
        law = scheme.saturation_law,
        cp = atmosphere.heat_capacity,
        Lv = atmosphere.latent_heat_condensation,
        Lf = atmosphere.latent_heat_fusion,
        g = model.planet.gravity,
    )
end

@inline function _reset_finite_multiplume_diagnostics!(parameterizations, ij, zero_value)
    parameterizations.rain_rate_convection[ij] = zero_value
    parameterizations.convective_precipitation_mass_correction[ij] =
        zero_value
    parameterizations.c08_triggered[ij] = zero_value
    parameterizations.c08_active_plumes[ij] = zero_value
    parameterizations.c08_mass_flux[ij] = zero_value
    parameterizations.c08_weighted_cape[ij] = zero_value
    parameterizations.c08_maximum_cin[ij] = zero_value
    parameterizations.c08_response[ij] = zero_value
    parameterizations.c08_limited[ij] = zero_value
    parameterizations.c08_degenerate[ij] = zero_value
    parameterizations.c08_water_residual[ij] = zero_value
    parameterizations.c08_water_gross_exchange[ij] = zero_value
    parameterizations.c08_enthalpy_residual[ij] = zero_value
    parameterizations.c08_enthalpy_gross_exchange[ij] = zero_value
    return nothing
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    scheme::FiniteMultiPlumeConvection,
    model,
)
    (; geometry, atmosphere, planet, time_stepping) = model
    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    geopotential = vars.grid.geopotential
    temperature_tendency = vars.tendencies.grid.temperature
    humidity_tendency = vars.tendencies.grid.humidity
    parameterizations = vars.parameterizations
    previous_snow_rate = parameterizations.snow_rate[ij]
    sigma = geometry.σ_levels_full
    sigma_thickness = geometry.σ_levels_thick
    number_layers = length(sigma)
    NF = eltype(temperature)
    zero_value = zero(NF)

    _reset_finite_multiplume_diagnostics!(parameterizations, ij, zero_value)
    surface_pressure = NF(vars.grid.pressure_prev[ij])
    parameterizations.convective_precipitation_surface_pressure[ij] =
        surface_pressure

    column_temperature = NF[
        temperature[ij, layer] for layer in 1:number_layers
    ]
    column_vapor = NF[
        max(humidity[ij, layer], zero_value) for layer in 1:number_layers
    ]
    column_ice = zeros(NF, number_layers)
    column_pressure = NF[
        surface_pressure * sigma[layer] for layer in 1:number_layers
    ]
    column_pressure_thickness = NF[
        surface_pressure * sigma_thickness[layer]
        for layer in 1:number_layers
    ]
    column_geopotential = NF[
        geopotential[ij, layer] for layer in 1:number_layers
    ]
    column = FiniteC08Core.Column(
        column_temperature,
        column_vapor,
        column_ice,
        column_pressure,
        column_pressure_thickness,
        column_geopotential,
    )
    thermo = _finite_multiplume_thermo(scheme, model, NF)
    timestep = NF(time_stepping.Δt_sec)
    exchange, diagnostics = C08FiniteConvection.solve(
        column,
        timestep,
        thermo,
        scheme.spectrum,
    )

    layer_mass = column_pressure_thickness ./ thermo.g
    water = FiniteC08Core.ConvectiveDetrainmentInterface.check_water(
        exchange,
        layer_mass,
    )
    enthalpy = FiniteC08Core.ConvectiveDetrainmentInterface.check_enthalpy(
        exchange,
        layer_mass;
        cp = thermo.cp,
        Lv = thermo.Lv,
        Lf = thermo.Lf,
    )
    parameterizations.c08_water_residual[ij] = water.residual
    parameterizations.c08_water_gross_exchange[ij] = water.gross
    parameterizations.c08_enthalpy_residual[ij] = enthalpy.residual
    parameterizations.c08_enthalpy_gross_exchange[ij] = enthalpy.gross

    weighted_cape = zero_value
    maximum_cin = zero_value
    active_plumes = 0
    cloud_top = number_layers + 1
    for plume in eachindex(diagnostics.active)
        diagnostics.active[plume] || continue
        active_plumes += 1
        weighted_cape += scheme.spectrum.weights[plume] *
            diagnostics.cape[plume]
        maximum_cin = max(maximum_cin, diagnostics.cin[plume])
        diagnostics.top[plume] > 0 &&
            (cloud_top = min(cloud_top, diagnostics.top[plume]))
    end
    parameterizations.c08_triggered[ij] = NF(diagnostics.triggered)
    parameterizations.c08_active_plumes[ij] = NF(active_plumes)
    parameterizations.c08_mass_flux[ij] = diagnostics.Mb
    parameterizations.c08_weighted_cape[ij] = weighted_cape
    parameterizations.c08_maximum_cin[ij] = maximum_cin
    parameterizations.c08_response[ij] = diagnostics.R
    parameterizations.finite_c08_status[ij] = FINITE_MULTIPLUME_STATUS_CODES[diagnostics.finite_closure_status]
    parameterizations.finite_c08_evaluations[ij] = diagnostics.finite_closure_evaluations
    parameterizations.finite_c08_invalid_probes[ij] = diagnostics.finite_closure_invalid_probes
    parameterizations.finite_c08_target_cape[ij] = diagnostics.finite_closure_target
    parameterizations.finite_c08_final_cape[ij] = diagnostics.finite_closure_final_cape
    parameterizations.finite_c08_initial_cape[ij] = diagnostics.finite_closure_initial_cape
    parameterizations.c08_limited[ij] = NF(diagnostics.limited)
    parameterizations.c08_degenerate[ij] = NF(diagnostics.degenerate)

    for layer in 1:number_layers
        host_humidity_tendency = exchange.dqv[layer] +
            exchange.dql[layer] + exchange.dqi[layer]
        host_temperature_tendency = exchange.dT[layer] -
            (thermo.Lv * exchange.dql[layer] +
             (thermo.Lv + thermo.Lf) * exchange.dqi[layer]) / thermo.cp
        humidity_tendency[ij, layer] += host_humidity_tendency
        temperature_tendency[ij, layer] += host_temperature_tendency
        parameterizations.convective_humidity_tendency[ij, layer] =
            host_humidity_tendency
        parameterizations.convective_temperature_tendency[ij, layer] =
            host_temperature_tendency
    end

    water_density = NF(atmosphere.water_density)
    rain_rate = exchange.rain / water_density
    snow_rate = exchange.snow / water_density
    parameterizations.rain_rate_convection[ij] = rain_rate
    parameterizations.rain_rate[ij] += rain_rate
    parameterizations.snow_rate[ij] += snow_rate
    parameterizations.snow_rate_convection[ij] =
        parameterizations.snow_rate[ij] - previous_snow_rate
    parameterizations.rain_convection[ij] += rain_rate * timestep
    cloud_top <= number_layers &&
        (parameterizations.cloud_top[ij] = min(
            parameterizations.cloud_top[ij],
            cloud_top,
        ))
    return nothing
end
