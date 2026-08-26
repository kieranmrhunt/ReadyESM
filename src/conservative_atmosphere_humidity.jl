"""
Remove negative grid-space specific humidity without changing the Gaussian-area
mean of any sigma layer.

Spectral transforms can produce small negative humidity values. Simply clipping
them to zero adds water to the grid state used by physics and radiation. This
scheme removes the negatives and reduces the remaining positive values in the
same layer by the exact amount introduced by clipping. No vertical water
transfer is made, and a layer whose pre-clipping global mean is non-positive is
set to zero.
"""
struct LayerConservativeHumidityHoleFilling{W, Q, L} <:
       SpeedyWeather.AbstractHoleFilling
    point_weights::W
    weighted_humidity::Q
    original_layer_mean::L
    positive_layer_mean::L
    layer_scale::L
end

Adapt.@adapt_structure LayerConservativeHumidityHoleFilling

function LayerConservativeHumidityHoleFilling(
    spectral_grid::SpeedyWeather.SpectralGrid,
)
    NF = spectral_grid.NF
    architecture = spectral_grid.architecture
    point_weights = _global_point_weights(spectral_grid)
    weighted_humidity = SpeedyWeather.on_architecture(
        architecture,
        zeros(NF, spectral_grid.npoints, spectral_grid.nlayers),
    )
    original_layer_mean = SpeedyWeather.on_architecture(
        architecture,
        zeros(NF, 1, spectral_grid.nlayers),
    )
    positive_layer_mean = similar(original_layer_mean)
    layer_scale = similar(original_layer_mean)
    return LayerConservativeHumidityHoleFilling(
        point_weights,
        weighted_humidity,
        original_layer_mean,
        positive_layer_mean,
        layer_scale,
    )
end

function SpeedyWeather.hole_filling!(
    humidity,
    scheme::LayerConservativeHumidityHoleFilling,
    ::SpeedyWeather.AbstractModel,
)
    q = humidity.data
    weights = reshape(scheme.point_weights, :, 1)

    scheme.weighted_humidity .= q .* weights
    sum!(scheme.original_layer_mean, scheme.weighted_humidity)

    q .= max.(q, zero(eltype(q)))
    scheme.weighted_humidity .= q .* weights
    sum!(scheme.positive_layer_mean, scheme.weighted_humidity)

    scheme.layer_scale .= max.(scheme.original_layer_mean, zero(eltype(q))) ./
        max.(scheme.positive_layer_mean, eps(eltype(q)))
    q .*= scheme.layer_scale
    return nothing
end

"""
Use the unchanged Betts--Miller temperature and humidity tendencies, but
diagnose surface convective rain from the net column humidity loss.

The upstream diagnostic sums positive layer drying while ignoring compensating
moistening in other layers. That sends more rain to the coupled surface than is
removed from the atmospheric column. This wrapper replaces only the diagnosed
convective rain rate and accumulation; it does not alter the convection
tendencies or cloud-top diagnosis.
"""
struct NetColumnBettsMillerConvection{C} <: SpeedyWeather.AbstractConvection
    convection::C
end

Adapt.@adapt_structure NetColumnBettsMillerConvection

function NetColumnBettsMillerConvection(
    spectral_grid::SpeedyWeather.SpectralGrid;
    kwargs...,
)
    convection = SpeedyWeather.BettsMillerConvection(spectral_grid; kwargs...)
    return NetColumnBettsMillerConvection{typeof(convection)}(convection)
end

"""
Betts--Miller convection with a piecewise-linear reference-relative-humidity
profile in physical sigma coordinates.

The upstream scheme uses one `0.7` value at every adjusted level. This variant
keeps the upstream pseudo-adiabat, convective criteria, thermodynamic profile
adjustment, relaxation time and precipitation diagnostics unchanged; only the
relative humidity used to construct the first-guess reference profile varies
with sigma. Setting the upper and lower values equal exactly recovers the
constant-profile algebra.
"""
struct SigmaProfileBettsMillerConvection{C, NF} <:
       SpeedyWeather.AbstractConvection
    convection::C
    upper_relative_humidity::NF
    lower_relative_humidity::NF
    transition_top_sigma::NF
    transition_bottom_sigma::NF
end

Adapt.@adapt_structure SigmaProfileBettsMillerConvection

function SigmaProfileBettsMillerConvection(
    spectral_grid::SpeedyWeather.SpectralGrid;
    upper_relative_humidity = 0.5,
    lower_relative_humidity = 0.7,
    transition_top_sigma = 0.3,
    transition_bottom_sigma = 0.7,
    kwargs...,
)
    NF = spectral_grid.NF
    return SigmaProfileBettsMillerConvection(
        SpeedyWeather.BettsMillerConvection(spectral_grid; kwargs...),
        NF(upper_relative_humidity),
        NF(lower_relative_humidity),
        NF(transition_top_sigma),
        NF(transition_bottom_sigma),
    )
end

SpeedyWeather.variables(convection::SigmaProfileBettsMillerConvection) =
    SpeedyWeather.variables(convection.convection)

SpeedyWeather.initialize!(
    convection::SigmaProfileBettsMillerConvection,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(convection.convection, model)

@inline function _sigma_profile_reference_humidity(
    convection::SigmaProfileBettsMillerConvection,
    sigma,
)
    fraction = clamp(
        (sigma - convection.transition_top_sigma) /
            (convection.transition_bottom_sigma -
             convection.transition_top_sigma),
        zero(sigma),
        one(sigma),
    )
    return convection.upper_relative_humidity + fraction *
        (convection.lower_relative_humidity -
         convection.upper_relative_humidity)
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    convection::SigmaProfileBettsMillerConvection,
    model,
)
    (; geometry, planet, atmosphere, time_stepping) = model
    sigma = geometry.σ_levels_full
    sigma_half = geometry.σ_levels_half
    sigma_thickness = geometry.σ_levels_thick
    nlayers = length(sigma)
    timestep_seconds = time_stepping.Δt_sec

    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    geopotential = vars.grid.geopotential
    temperature_tendency = vars.tendencies.grid.temperature
    humidity_tendency = vars.tendencies.grid.humidity
    surface_pressure = vars.grid.pressure_prev[ij]
    NF = eltype(temperature)

    water_density = atmosphere.water_density
    gravity = planet.gravity
    latent_heat = atmosphere.latent_heat_condensation
    heat_capacity = atmosphere.heat_capacity

    temperature_reference = vars.scratch.grid.a
    humidity_reference = vars.scratch.grid.b
    level_zero_buoyancy = SpeedyWeather.pseudo_adiabat!(
        ij,
        temperature_reference,
        temperature,
        humidity,
        geopotential,
        surface_pressure,
        sigma,
        atmosphere,
    )

    for k in level_zero_buoyancy:nlayers
        saturation_humidity = SpeedyWeather.saturation_humidity(
            temperature_reference[ij, k],
            surface_pressure * sigma[k],
            atmosphere,
        )
        humidity_reference[ij, k] = saturation_humidity *
            _sigma_profile_reference_humidity(convection, sigma[k])
    end

    humidity_precipitation = zero(NF)
    temperature_precipitation = zero(NF)
    reference_humidity_integral = zero(NF)
    for k in level_zero_buoyancy:nlayers
        humidity_precipitation +=
            (humidity[ij, k] - humidity_reference[ij, k]) *
            sigma_thickness[k]
        temperature_precipitation -=
            (temperature[ij, k] - temperature_reference[ij, k]) *
            sigma_thickness[k]
    end

    deep_convection = humidity_precipitation > 0 &&
        temperature_precipitation > 0
    shallow_convection = humidity_precipitation <= 0 &&
        temperature_precipitation > 0
    !(deep_convection || shallow_convection) && return nothing

    convective_sigma_depth =
        sigma_half[nlayers + 1] - sigma_half[level_zero_buoyancy]
    if deep_convection
        uniform_temperature_adjustment =
            (temperature_precipitation -
             humidity_precipitation * latent_heat / heat_capacity) /
            convective_sigma_depth
        for k in level_zero_buoyancy:nlayers
            temperature_reference[ij, k] -= uniform_temperature_adjustment
        end
    else
        for k in level_zero_buoyancy:nlayers
            reference_humidity_integral -=
                humidity_reference[ij, k] * sigma_thickness[k]
        end
        humidity_scale =
            1 - humidity_precipitation / reference_humidity_integral
        uniform_temperature_adjustment =
            temperature_precipitation / convective_sigma_depth
        for k in level_zero_buoyancy:nlayers
            humidity_reference[ij, k] *= humidity_scale
            temperature_reference[ij, k] -= uniform_temperature_adjustment
        end
    end

    inverse_timescale = inv(convert(
        NF,
        Second(convection.convection.time_scale).value,
    ))
    convective_rain = zero(NF)
    for k in level_zero_buoyancy:nlayers
        temperature_tendency[ij, k] -=
            (temperature[ij, k] - temperature_reference[ij, k]) *
            inverse_timescale
        humidity_adjustment =
            (humidity[ij, k] - humidity_reference[ij, k]) *
            inverse_timescale
        humidity_tendency[ij, k] -= humidity_adjustment
        convective_rain +=
            max(humidity_adjustment * sigma_thickness[k], zero(NF))
    end

    pressure_time_over_gravity_density =
        surface_pressure * timestep_seconds /
        (gravity * water_density) * deep_convection
    convective_rain = max(
        convective_rain * pressure_time_over_gravity_density,
        zero(NF),
    )
    vars.parameterizations.rain_convection[ij] += convective_rain
    convective_rain_rate = convective_rain / timestep_seconds
    vars.parameterizations.rain_rate_convection[ij] = convective_rain_rate
    vars.parameterizations.rain_rate[ij] += convective_rain_rate
    vars.parameterizations.cloud_top[ij] = min(
        vars.parameterizations.cloud_top[ij],
        level_zero_buoyancy,
    )
    return nothing
end

function _atmosphere_convection(
    config::ExperimentConfig,
    spectral_grid::SpeedyWeather.SpectralGrid,
)
    convection = if config.atmosphere_convection ==
                    :betts_miller_constant_rh
        SpeedyWeather.BettsMillerConvection(spectral_grid)
    elseif config.atmosphere_convection == :betts_miller_sigma_rh_v1
        SigmaProfileBettsMillerConvection(
            spectral_grid;
            upper_relative_humidity = config.
                atmosphere_convection_upper_relative_humidity,
            lower_relative_humidity = config.
                atmosphere_convection_lower_relative_humidity,
            transition_top_sigma = config.
                atmosphere_convection_transition_top_sigma,
            transition_bottom_sigma = config.
                atmosphere_convection_transition_bottom_sigma,
        )
    else
        error("unsupported atmosphere convection $(config.atmosphere_convection)")
    end
    return NetColumnBettsMillerConvection(convection)
end

SpeedyWeather.variables(convection::NetColumnBettsMillerConvection) = (
    SpeedyWeather.variables(convection.convection)...,
    SpeedyWeather.ParameterizationVariable(
        :convective_humidity_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental Betts-Miller specific-humidity tendency",
        units = "kg/kg/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_temperature_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental Betts-Miller temperature tendency",
        units = "K/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :convective_precipitation_mass_correction,
        SpeedyWeather.Grid2D();
        desc = "Net-column minus upstream convective precipitation rate",
        units = "m/s",
    ),
)

SpeedyWeather.initialize!(
    convection::NetColumnBettsMillerConvection,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(convection.convection, model)

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    convection::NetColumnBettsMillerConvection,
    model,
)
    humidity_tendency = vars.tendencies.grid.humidity
    temperature_tendency = vars.tendencies.grid.temperature
    convective_humidity_tendency =
        vars.parameterizations.convective_humidity_tendency
    convective_temperature_tendency =
        vars.parameterizations.convective_temperature_tendency
    layer_thickness = model.geometry.σ_levels_thick
    column_tendency_before = zero(eltype(humidity_tendency))
    for k in 1:length(layer_thickness)
        convective_humidity_tendency[ij, k] = humidity_tendency[ij, k]
        convective_temperature_tendency[ij, k] = temperature_tendency[ij, k]
        column_tendency_before +=
            humidity_tendency[ij, k] * layer_thickness[k]
    end
    rain_rate_before = vars.parameterizations.rain_rate[ij]

    SpeedyWeather.parameterization!(
        ij,
        vars,
        convection.convection,
        model,
    )

    column_tendency_after = zero(column_tendency_before)
    for k in 1:length(layer_thickness)
        convective_humidity_tendency[ij, k] =
            humidity_tendency[ij, k] - convective_humidity_tendency[ij, k]
        convective_temperature_tendency[ij, k] =
            temperature_tendency[ij, k] - convective_temperature_tendency[ij, k]
        column_tendency_after +=
            humidity_tendency[ij, k] * layer_thickness[k]
    end
    gross_convective_rain_rate =
        vars.parameterizations.rain_rate[ij] - rain_rate_before
    net_convective_rain_rate = max(
        -(column_tendency_after - column_tendency_before) *
            vars.grid.pressure_prev[ij] /
            (model.planet.gravity * model.atmosphere.water_density),
        zero(gross_convective_rain_rate),
    )
    correction = net_convective_rain_rate - gross_convective_rain_rate

    vars.parameterizations.rain_rate[ij] += correction
    vars.parameterizations.rain_rate_convection[ij] =
        net_convective_rain_rate
    vars.parameterizations.convective_precipitation_mass_correction[ij] =
        correction
    vars.parameterizations.rain_convection[ij] +=
        model.time_stepping.Δt_sec * correction
    return nothing
end

"""
Run Betts--Miller unchanged, then delay a configured fraction of the diagnosed
surface precipitation in a layer-resolved liquid/ice condensate reservoir.

This experimental wrapper is deliberately column-local. It provides the first
conservative prognostic-cloud test while leaving horizontal condensate
transport and condensate enthalpy for a later tracer implementation. Vapour
tendencies are unchanged; surface precipitation is reduced by retained
condensate and increased by sedimentation, so vapour plus cloud water closes
exactly to the original precipitation ledger.
"""
struct PrognosticCondensateConvection{C, S} <:
       SpeedyWeather.AbstractConvection
    convection::C
    condensate::S
end

Adapt.@adapt_structure PrognosticCondensateConvection

SpeedyWeather.variables(convection::PrognosticCondensateConvection) =
    SpeedyWeather.variables(convection.convection)

SpeedyWeather.initialize!(
    convection::PrognosticCondensateConvection,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(convection.convection, model)

@inline function _advance_prognostic_condensate_column!(
    condensate,
    convective_humidity_tendency,
    large_scale_humidity_tendency,
    temperature,
    sigma_thickness,
    surface_pressure,
    gravity,
    water_density,
    timestep_seconds,
    old_convective_rain_rate,
    old_large_scale_rain_rate,
    old_large_scale_snow_rate,
    ij,
)
    NF = eltype(condensate.liquid_path_gm2)
    zero_nf = zero(NF)
    if condensate.initialized[1] == 0
        return (
            convective_rain_rate = old_convective_rain_rate,
            large_scale_rain_rate = old_large_scale_rain_rate,
            large_scale_snow_rate = old_large_scale_snow_rate,
        )
    end

    total_drying_mass_rate = zero_nf
    nlayers = length(sigma_thickness)
    for k in 1:nlayers
        humidity_tendency =
            convective_humidity_tendency[ij, k] +
            large_scale_humidity_tendency[ij, k]
        layer_air_mass = surface_pressure * sigma_thickness[k] / gravity
        total_drying_mass_rate +=
            max(-humidity_tendency, zero_nf) * layer_air_mass
    end

    convective_rain_rate = max(old_convective_rain_rate, zero_nf)
    large_scale_rain_rate = max(old_large_scale_rain_rate, zero_nf)
    large_scale_snow_rate = max(old_large_scale_snow_rate, zero_nf)
    original_precipitation_rate =
        convective_rain_rate + large_scale_rain_rate + large_scale_snow_rate
    retained_mass_rate = total_drying_mass_rate > 0 ?
        condensate.retention_fraction * original_precipitation_rate *
            water_density : zero_nf

    sedimented_liquid_gm2s = zero_nf
    sedimented_ice_gm2s = zero_nf
    for k in 1:nlayers
        humidity_tendency =
            convective_humidity_tendency[ij, k] +
            large_scale_humidity_tendency[ij, k]
        layer_air_mass = surface_pressure * sigma_thickness[k] / gravity
        drying_mass_rate =
            max(-humidity_tendency, zero_nf) * layer_air_mass
        source_mass_rate = total_drying_mass_rate > 0 ?
            retained_mass_rate * drying_mass_rate / total_drying_mass_rate :
            zero_nf
        liquid_fraction = clamp(
            (temperature[ij, k] - NF(253)) / NF(20),
            zero_nf,
            one(NF),
        )
        liquid_source_gm2s = NF(1_000) * source_mass_rate * liquid_fraction
        ice_source_gm2s = NF(1_000) * source_mass_rate * (1 - liquid_fraction)

        old_liquid = max(condensate.liquid_path_gm2[ij, k], zero_nf)
        old_ice = max(condensate.ice_path_gm2[ij, k], zero_nf)
        liquid_sink_gm2s = min(
            old_liquid / condensate.residence_time_seconds,
            old_liquid / timestep_seconds,
        )
        ice_sink_gm2s = min(
            old_ice / condensate.residence_time_seconds,
            old_ice / timestep_seconds,
        )
        condensate.liquid_path_gm2[ij, k] = max(
            old_liquid + timestep_seconds *
                (liquid_source_gm2s - liquid_sink_gm2s),
            zero_nf,
        )
        condensate.ice_path_gm2[ij, k] = max(
            old_ice + timestep_seconds *
                (ice_source_gm2s - ice_sink_gm2s),
            zero_nf,
        )
        sedimented_liquid_gm2s += liquid_sink_gm2s
        sedimented_ice_gm2s += ice_sink_gm2s
    end

    sedimented_mass_rate =
        (sedimented_liquid_gm2s + sedimented_ice_gm2s) / NF(1_000)
    condensate.cumulative_retained_kgm2[ij] +=
        timestep_seconds * retained_mass_rate
    condensate.cumulative_sedimented_kgm2[ij] +=
        timestep_seconds * sedimented_mass_rate

    retained_fraction = original_precipitation_rate > 0 ?
        retained_mass_rate /
            (original_precipitation_rate * water_density) : zero_nf
    immediate_fraction = clamp(1 - retained_fraction, zero_nf, one(NF))
    liquid_sedimentation_rate =
        sedimented_liquid_gm2s / (NF(1_000) * water_density)
    ice_sedimentation_rate =
        sedimented_ice_gm2s / (NF(1_000) * water_density)
    return (
        convective_rain_rate = convective_rain_rate * immediate_fraction,
        large_scale_rain_rate =
            large_scale_rain_rate * immediate_fraction +
            liquid_sedimentation_rate,
        large_scale_snow_rate =
            large_scale_snow_rate * immediate_fraction +
            ice_sedimentation_rate,
    )
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    convection::PrognosticCondensateConvection,
    model,
)
    SpeedyWeather.parameterization!(
        ij,
        vars,
        convection.convection,
        model,
    )

    parameterizations = vars.parameterizations
    old_convective_rain_rate = parameterizations.rain_rate_convection[ij]
    old_large_scale_rain_rate = parameterizations.rain_rate_large_scale[ij]
    old_large_scale_snow_rate = parameterizations.snow_rate_large_scale[ij]
    rates = _advance_prognostic_condensate_column!(
        convection.condensate,
        parameterizations.convective_humidity_tendency,
        parameterizations.large_scale_condensation_humidity_tendency,
        vars.grid.temperature_prev,
        model.geometry.σ_levels_thick,
        vars.grid.pressure_prev[ij],
        model.planet.gravity,
        model.atmosphere.water_density,
        model.time_stepping.Δt_sec,
        old_convective_rain_rate,
        old_large_scale_rain_rate,
        old_large_scale_snow_rate,
        ij,
    )

    parameterizations.rain_rate_convection[ij] = rates.convective_rain_rate
    parameterizations.rain_rate_large_scale[ij] = rates.large_scale_rain_rate
    parameterizations.snow_rate_large_scale[ij] = rates.large_scale_snow_rate
    parameterizations.rain_rate[ij] =
        rates.convective_rain_rate + rates.large_scale_rain_rate
    parameterizations.snow_rate[ij] = rates.large_scale_snow_rate
    timestep_seconds = model.time_stepping.Δt_sec
    parameterizations.rain_convection[ij] += timestep_seconds *
        (rates.convective_rain_rate - old_convective_rain_rate)
    parameterizations.rain_large_scale[ij] += timestep_seconds *
        (rates.large_scale_rain_rate - old_large_scale_rain_rate)
    parameterizations.snow_large_scale[ij] += timestep_seconds *
        (rates.large_scale_snow_rate - old_large_scale_snow_rate)
    return nothing
end

"""
Run SpeedyWeather's upstream implicit condensation without modifying its
thermodynamic tendencies, surface precipitation rates or accumulators, while
retaining the signed layer-by-layer process contributions.

This is the control implementation for diagnosing the known upstream mismatch
between column humidity drying and surface precipitation. The explicitly saved
mass-correction field is zero because no correction is applied.
"""
struct ObservedImplicitCondensation{C} <:
       SpeedyWeather.AbstractCondensation
    condensation::C
end

Adapt.@adapt_structure ObservedImplicitCondensation

function ObservedImplicitCondensation(
    spectral_grid::SpeedyWeather.SpectralGrid;
    kwargs...,
)
    condensation = SpeedyWeather.ImplicitCondensation(
        spectral_grid;
        kwargs...,
    )
    return ObservedImplicitCondensation{typeof(condensation)}(
        condensation,
    )
end

SpeedyWeather.variables(condensation::ObservedImplicitCondensation) = (
    SpeedyWeather.variables(condensation.condensation)...,
    SpeedyWeather.ParameterizationVariable(
        :large_scale_condensation_humidity_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental large-scale-condensation specific-humidity tendency",
        units = "kg/kg/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :large_scale_condensation_temperature_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental large-scale-condensation temperature tendency",
        units = "K/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :large_scale_precipitation_mass_correction,
        SpeedyWeather.Grid2D();
        desc = "Applied correction to upstream large-scale precipitation rate",
        units = "m/s",
    ),
)

SpeedyWeather.initialize!(
    condensation::ObservedImplicitCondensation,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(condensation.condensation, model)

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    condensation::ObservedImplicitCondensation,
    model,
)
    humidity_tendency = vars.tendencies.grid.humidity
    temperature_tendency = vars.tendencies.grid.temperature
    process_humidity_tendency =
        vars.parameterizations.large_scale_condensation_humidity_tendency
    process_temperature_tendency =
        vars.parameterizations.large_scale_condensation_temperature_tendency
    nlayers = length(model.geometry.σ_levels_thick)

    for k in 1:nlayers
        process_humidity_tendency[ij, k] = humidity_tendency[ij, k]
        process_temperature_tendency[ij, k] = temperature_tendency[ij, k]
    end

    SpeedyWeather.parameterization!(
        ij,
        vars,
        condensation.condensation,
        model,
    )

    for k in 1:nlayers
        process_humidity_tendency[ij, k] =
            humidity_tendency[ij, k] - process_humidity_tendency[ij, k]
        process_temperature_tendency[ij, k] =
            temperature_tendency[ij, k] - process_temperature_tendency[ij, k]
    end
    vars.parameterizations.large_scale_precipitation_mass_correction[ij] =
        zero(eltype(humidity_tendency))
    return nothing
end

"""
Use SpeedyWeather's unchanged implicit large-scale-condensation thermodynamic
tendencies, but diagnose surface rain plus snow from the net column humidity
loss and retain the signed layer-by-layer process contributions.

The upstream scheme removes the full requested rain re-evaporation from its
falling flux before applying the implicit relaxation factor to the humidity
tendency. Consequently, only part of the removed rain returns to atmospheric
water and the remainder disappears. Preserve the diagnosed snow flux where
possible and assign the column-conserving remainder to rain. The layer fields
are overwritten on every physics call and do not feed back on the tendencies.
"""
struct NetColumnImplicitCondensation{C} <:
       SpeedyWeather.AbstractCondensation
    condensation::C
end

@inline function _conservative_large_scale_precipitation_rates(
    target_precipitation_rate,
    gross_snow_rate,
)
    target = max(target_precipitation_rate, zero(target_precipitation_rate))
    corrected_snow_rate = min(
        max(gross_snow_rate, zero(gross_snow_rate)),
        target,
    )
    return target - corrected_snow_rate, corrected_snow_rate
end

Adapt.@adapt_structure NetColumnImplicitCondensation

function NetColumnImplicitCondensation(
    spectral_grid::SpeedyWeather.SpectralGrid;
    kwargs...,
)
    condensation = SpeedyWeather.ImplicitCondensation(
        spectral_grid;
        kwargs...,
    )
    return NetColumnImplicitCondensation{typeof(condensation)}(
        condensation,
    )
end

SpeedyWeather.variables(condensation::NetColumnImplicitCondensation) = (
    SpeedyWeather.variables(condensation.condensation)...,
    SpeedyWeather.ParameterizationVariable(
        :large_scale_condensation_humidity_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental large-scale-condensation specific-humidity tendency",
        units = "kg/kg/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :large_scale_condensation_temperature_tendency,
        SpeedyWeather.Grid3D();
        desc = "Incremental large-scale-condensation temperature tendency",
        units = "K/s",
    ),
    SpeedyWeather.ParameterizationVariable(
        :large_scale_precipitation_mass_correction,
        SpeedyWeather.Grid2D();
        desc = "Net-column minus upstream large-scale precipitation rate",
        units = "m/s",
    ),
)

SpeedyWeather.initialize!(
    condensation::NetColumnImplicitCondensation,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(condensation.condensation, model)

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    condensation::NetColumnImplicitCondensation,
    model,
)
    humidity_tendency = vars.tendencies.grid.humidity
    temperature_tendency = vars.tendencies.grid.temperature
    process_humidity_tendency =
        vars.parameterizations.large_scale_condensation_humidity_tendency
    process_temperature_tendency =
        vars.parameterizations.large_scale_condensation_temperature_tendency
    layer_thickness = model.geometry.σ_levels_thick
    nlayers = length(layer_thickness)
    column_tendency_before = zero(eltype(humidity_tendency))

    for k in 1:nlayers
        process_humidity_tendency[ij, k] = humidity_tendency[ij, k]
        process_temperature_tendency[ij, k] = temperature_tendency[ij, k]
        column_tendency_before +=
            humidity_tendency[ij, k] * layer_thickness[k]
    end

    SpeedyWeather.parameterization!(
        ij,
        vars,
        condensation.condensation,
        model,
    )

    column_tendency_after = zero(column_tendency_before)
    for k in 1:nlayers
        process_humidity_tendency[ij, k] =
            humidity_tendency[ij, k] - process_humidity_tendency[ij, k]
        process_temperature_tendency[ij, k] =
            temperature_tendency[ij, k] - process_temperature_tendency[ij, k]
        column_tendency_after +=
            humidity_tendency[ij, k] * layer_thickness[k]
    end

    parameterizations = vars.parameterizations
    target_precipitation_rate = max(
        -(column_tendency_after - column_tendency_before) *
            vars.grid.pressure_prev[ij] /
            (model.planet.gravity * model.atmosphere.water_density),
        zero(column_tendency_after),
    )
    gross_rain_rate = parameterizations.rain_rate_large_scale[ij]
    gross_snow_rate = parameterizations.snow_rate_large_scale[ij]
    corrected_rain_rate, corrected_snow_rate =
        _conservative_large_scale_precipitation_rates(
            target_precipitation_rate,
            gross_snow_rate,
        )
    rain_correction = corrected_rain_rate - gross_rain_rate
    snow_correction = corrected_snow_rate - gross_snow_rate

    parameterizations.rain_rate_large_scale[ij] = corrected_rain_rate
    parameterizations.snow_rate_large_scale[ij] = corrected_snow_rate
    parameterizations.rain_rate[ij] += rain_correction
    parameterizations.snow_rate[ij] += snow_correction
    parameterizations.large_scale_precipitation_mass_correction[ij] =
        rain_correction + snow_correction
    parameterizations.rain_large_scale[ij] +=
        model.time_stepping.Δt_sec * rain_correction
    parameterizations.snow_large_scale[ij] +=
        model.time_stepping.Δt_sec * snow_correction
    return nothing
end
