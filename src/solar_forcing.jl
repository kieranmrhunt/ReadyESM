"""
Present-day inverse-square Earth–Sun distance approximation.

Spencer coefficients as given in Wang et al. (2024), equations 2–3,
https://doi.org/10.5194/essd-16-961-2024. The phase follows the fraction
of the Gregorian year, including 366-day leap years. This calendar
extension is intended for present-day seasonal forcing.
"""
function _solar_distance_factor(time::DateTime)
    start = DateTime(year(time), 1, 1)
    stop = DateTime(year(time) + 1, 1, 1)
    phase = 2π * Dates.value(time - start) / Dates.value(stop - start)
    return 1.000110 + 0.034221*cos(phase) + 0.001280*sin(phase) +
        0.000719*cos(2phase) + 0.000077*sin(2phase)
end

function _solar_irradiance(solar_constant::T, time::DateTime, enabled::Bool) where {T<:AbstractFloat}
    enabled || return solar_constant
    return solar_constant * T(_solar_distance_factor(time))
end

_solar_distance_scheme(enabled::Bool) = enabled ?
    "spencer_present_day_gregorian_year_fraction_v1" : "constant_distance"

function _forcing_time_dependence(config)
    carbon = config.carbon.mode == :emissions_driven
    solar = config.forcing.solar_distance_correction
    carbon && solar && return "CO2_updated_each_coupling_step_and_seasonal_solar_distance_from_model_clock"
    carbon && return "CO2_updated_each_coupling_step_from_emissions_and_air_sea_exchange"
    solar && return "seasonal_solar_distance_from_model_clock"
    return "constant_for_entire_integration"
end

function _validate_solar_restart(radiation, state)
    saved = hasproperty(state, :solar_distance_correction) ?
        state.solar_distance_correction : false
    saved isa Bool || throw(ArgumentError("checkpoint solar-distance option must be Boolean"))
    saved == radiation.solar_distance_correction || throw(ArgumentError(
        "checkpoint and reconstructed radiation use different solar-distance options",
    ))
    return nothing
end
