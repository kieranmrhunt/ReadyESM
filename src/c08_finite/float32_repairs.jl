# The inherited phase interface stores its constants as Float64. These methods
# keep Float32 inputs in Float32 without changing Float64 dispatch.
@inline function EP.ConvectiveDetrainmentInterface.liquid_fraction(
    temperature::Float32,
)
    freezing = Float32(EP.ConvectiveDetrainmentInterface.T_FREEZE)
    all_ice = Float32(EP.ConvectiveDetrainmentInterface.T_ALL_ICE)
    temperature >= freezing && return one(temperature)
    fraction = clamp(
        (temperature - all_ice) / (freezing - all_ice),
        zero(temperature),
        one(temperature),
    )
    return fraction * fraction
end

@inline function EP.ConvectiveDetrainmentInterface.dliquid_fraction_dT(
    temperature::Float32,
)
    freezing = Float32(EP.ConvectiveDetrainmentInterface.T_FREEZE)
    all_ice = Float32(EP.ConvectiveDetrainmentInterface.T_ALL_ICE)
    fraction = (temperature - all_ice) / (freezing - all_ice)
    return zero(temperature) < fraction < one(temperature) ?
        Float32(2) * fraction / (freezing - all_ice) : zero(temperature)
end

# The generic prototype divides by the Float64 literal 1.0e4 during rain
# evaporation and snow sublimation. This Float32 method retains the original
# equations and operation order with a typed pressure scale.
function EP.fallout(
    exchange::ColumnExchange,
    column::Column,
    profile,
    mass,
    mass_flux::Float32,
    effective_timestep::Float32,
    thermo::Thermo{Float32},
    parameters::PlumeParameters{Float32},
)
    number_layers = length(column.T)
    FT = Float32
    vapor_tendency = copy(exchange.dqv)
    temperature_tendency = copy(exchange.dT)
    rain = zero(FT)
    snow = zero(FT)
    pressure_scale = FT(1.0e4)
    for layer in 1:number_layers
        if rain + snow > zero(FT)
            temperature = column.T[layer] +
                effective_timestep * temperature_tendency[layer]
            if parameters.melt_below_freezing &&
                temperature > thermo.T0 && snow > zero(FT)
                melt = min(
                    snow,
                    thermo.cp * mass[layer] * (temperature - thermo.T0) /
                        (thermo.Lf * effective_timestep),
                )
                temperature_tendency[layer] -= thermo.Lf * melt /
                    (thermo.cp * mass[layer])
                temperature -= effective_timestep * thermo.Lf * melt /
                    (thermo.cp * mass[layer])
                rain += melt
                snow -= melt
            end
            if profile.base > 0 && layer > profile.base
                saturation = EP.saturation_humidity(
                    temperature,
                    column.p[layer],
                    thermo,
                )
                relative_humidity = clamp(
                    column.qv[layer] / saturation,
                    zero(FT),
                    one(FT),
                )
                capacity = max(saturation - column.qv[layer], zero(FT)) /
                    (one(FT) + thermo.Lv / thermo.cp *
                        EP.dqs_dT(temperature, column.p[layer], thermo)) *
                    mass[layer] / effective_timestep
                evaporated_rain = min(
                    rain * (one(FT) - exp(
                        -parameters.rain_evaporation *
                            (one(FT) - relative_humidity) *
                            column.dp[layer] / pressure_scale,
                    )),
                    capacity,
                )
                sublimated_snow = max(
                    min(
                        snow * (one(FT) - exp(
                            -parameters.snow_sublimation *
                                (one(FT) - relative_humidity) *
                                column.dp[layer] / pressure_scale,
                        )),
                        capacity - evaporated_rain,
                    ),
                    zero(FT),
                )
                rain -= evaporated_rain
                snow -= sublimated_snow
                vapor_tendency[layer] +=
                    (evaporated_rain + sublimated_snow) / mass[layer]
                temperature_tendency[layer] -=
                    (thermo.Lv * evaporated_rain +
                        (thermo.Lv + thermo.Lf) * sublimated_snow) /
                    (thermo.cp * mass[layer])
            end
        end
        rain += mass_flux * profile.Pl[layer]
        snow += mass_flux * profile.Pi[layer]
    end
    return ColumnExchange(
        vapor_tendency,
        exchange.dql,
        exchange.dqi,
        temperature_tendency,
        rain,
        snow,
    )
end
