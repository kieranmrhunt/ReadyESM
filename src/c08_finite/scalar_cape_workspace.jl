# Scalar CAPE evaluation used by the finite response closure.

@inline _c08_temperature(column, unit, scale, layer) =
    column.T[layer] + scale * unit.dT[layer]
@inline _c08_vapor(column, unit, scale, layer) =
    max(column.qv[layer] + scale * unit.dqv[layer], zero(eltype(column.T)))
@inline _c08_ice(column, unit, scale, layer) =
    max(column.qi[layer] + scale * unit.dqi[layer], zero(eltype(column.T)))

"""Dilute CAPE for a perturbed environment without constructing profile arrays."""
function dilute_cape_perturbed_scalar(column, thermo, parameters, unit, scale)
    n = length(column.T)
    FT = eltype(column.T)

    source = n
    depth = column.dp[n]
    while source > 1 && depth < parameters.source_depth_Pa
        source -= 1
        depth += column.dp[source]
    end
    top_weight = source == n ? one(FT) : clamp(
        (parameters.source_depth_Pa - (depth - column.dp[source])) /
            column.dp[source],
        zero(FT),
        one(FT),
    )

    source_mass = zero(FT)
    for layer in source:n
        weight = layer == source ? top_weight : one(FT)
        source_mass += weight * column.dp[layer] / thermo.g
    end

    mass_in = zero(FT)
    enthalpy_in = zero(FT)
    total_water_in = zero(FT)
    plume_temperature = FT(NaN)
    plume_vapor = zero(FT)
    plume_liquid = zero(FT)
    plume_ice = zero(FT)
    plume_enthalpy = zero(FT)
    plume_total_water = zero(FT)

    for layer in n:-1:source
        weight = layer == source ? top_weight : one(FT)
        environment_temperature = _c08_temperature(column, unit, scale, layer)
        environment_vapor = _c08_vapor(column, unit, scale, layer)
        environment_ice = _c08_ice(column, unit, scale, layer)
        environment_enthalpy = thermo.cp * environment_temperature +
            column.Φ[layer] + thermo.Lv * environment_vapor -
            thermo.Lf * environment_ice
        environment_total_water = environment_vapor + environment_ice
        entrainment = weight * column.dp[layer] / thermo.g / source_mass
        plume_enthalpy = (mass_in * enthalpy_in +
            entrainment * environment_enthalpy) / (mass_in + entrainment)
        plume_total_water = (mass_in * total_water_in +
            entrainment * environment_total_water) / (mass_in + entrainment)
        plume_temperature, plume_vapor, plume_liquid, plume_ice = EP.adjust(
            plume_enthalpy,
            plume_total_water,
            column.p[layer],
            column.Φ[layer],
            thermo,
        )
        mass_in += entrainment
        enthalpy_in = plume_enthalpy
        total_water_in = plume_total_water
    end

    cloud_base = 0
    level_free_convection = 0
    plume_top = 0
    cape = zero(FT)
    cin = zero(FT)
    base_saturation = zero(FT)
    if plume_liquid + plume_ice > zero(FT)
        cloud_base = source
        base_saturation = EP.saturation_humidity(
            plume_temperature,
            column.p[source],
            thermo,
        )
    end

    layer = source
    while layer > 1
        layer -= 1
        dz = (column.Φ[layer] - column.Φ[layer + 1]) / thermo.g
        environment_temperature = _c08_temperature(column, unit, scale, layer)
        environment_vapor = _c08_vapor(column, unit, scale, layer)
        environment_ice = _c08_ice(column, unit, scale, layer)
        environment_enthalpy = thermo.cp * environment_temperature +
            column.Φ[layer] + thermo.Lv * environment_vapor -
            thermo.Lf * environment_ice
        environment_total_water = environment_vapor + environment_ice

        entrainment_rate = if cloud_base == 0
            parameters.subcloud_entrainment_m
        elseif parameters.rh_dependent_entrainment
            environment_saturation = EP.saturation_humidity(
                environment_temperature,
                column.p[layer],
                thermo,
            )
            relative_humidity = clamp(
                environment_vapor / environment_saturation,
                zero(FT),
                one(FT),
            )
            max(
                parameters.entrainment_m * (FT(1.3) - relative_humidity) *
                    (environment_saturation / base_saturation)^3,
                parameters.min_entrainment_m,
            )
        else
            parameters.entrainment_m
        end
        detrainment_rate = cloud_base == 0 ? zero(FT) :
            parameters.detrainment_m
        entrainment = mass_in * expm1(min(
            entrainment_rate * dz,
            parameters.max_layer_entrainment,
        ))
        mixed_mass = mass_in + entrainment
        mixed_enthalpy = (mass_in * plume_enthalpy +
            entrainment * environment_enthalpy) / mixed_mass
        mixed_total_water = (mass_in * plume_total_water +
            entrainment * environment_total_water) / mixed_mass
        temperature, vapor, liquid, ice = EP.adjust(
            mixed_enthalpy,
            mixed_total_water,
            column.p[layer],
            column.Φ[layer],
            thermo,
        )

        physical = isfinite(temperature) &&
            parameters.min_parcel_temperature < temperature < 400 &&
            isfinite(vapor) && isfinite(liquid) && isfinite(ice)
        if !physical
            if level_free_convection > 0
                plume_top = layer + 1
            else
                plume_top = 0
                level_free_convection = 0
            end
            break
        end

        condensate = liquid + ice
        if condensate > parameters.condensate_threshold
            precipitation = (condensate - parameters.condensate_threshold) *
                (one(FT) - exp(-parameters.precip_conversion_m * dz))
            liquid_fraction = liquid / condensate
            liquid_precipitation = precipitation * liquid_fraction
            ice_precipitation = precipitation * (one(FT) - liquid_fraction)
            liquid -= liquid_precipitation
            ice -= ice_precipitation
            mixed_total_water -= precipitation
            mixed_enthalpy += thermo.Lf * ice_precipitation
        end

        detrainment = -mixed_mass * expm1(-detrainment_rate * dz)
        plume_density_temperature = EP.virtual_density_temperature(
            temperature,
            vapor,
            liquid + ice,
            thermo,
        )
        environment_density_temperature = EP.virtual_density_temperature(
            environment_temperature,
            environment_vapor,
            environment_ice,
            thermo,
        )
        buoyancy = thermo.g *
            (plume_density_temperature - environment_density_temperature) /
            environment_density_temperature

        if cloud_base == 0 && condensate > zero(FT)
            cloud_base = layer
            base_saturation = EP.saturation_humidity(
                temperature,
                column.p[layer],
                thermo,
            )
        end
        if buoyancy > zero(FT) && cloud_base > 0 && level_free_convection == 0
            level_free_convection = layer
        end
        if buoyancy <= zero(FT) && level_free_convection > 0
            detrainment = mixed_mass
            plume_top = layer
        elseif buoyancy <= zero(FT)
            cin -= buoyancy * dz
            if cin > parameters.cin_max
                detrainment = mixed_mass
                plume_top = layer
            end
        else
            cape += buoyancy * dz
        end
        if layer == 1 && plume_top == 0
            detrainment = mixed_mass
            plume_top = 1
        end

        plume_temperature = temperature
        plume_vapor = vapor
        plume_liquid = liquid
        plume_ice = ice
        plume_enthalpy = mixed_enthalpy
        plume_total_water = mixed_total_water
        mass_in = mixed_mass - detrainment
        plume_top > 0 && break
    end
    return cape
end
