const _CARBON_STAGE0B_ITERATIONS = 64
const _CARBON_STAGE0B_BRACKET_RESIDUAL_GUARD = 1e-5
const _CARBON_MOLAR_MASS_G_MOL = 12.011
const _CARBON_PG_PER_G = 1e-15

@inline function _carbon_device_constants(temperature_kelvin::FT, salinity::FT) where {FT<:AbstractFloat}
    temperature_100 = temperature_kelvin / FT(100)
    log_temperature = log(temperature_kelvin)
    root_salinity = sqrt(salinity)
    log_k0 = FT(-60.2409) + FT(93.4517) / temperature_100 +
             FT(23.3585) * log(temperature_100) +
             salinity * (
                 FT(0.023517) - FT(0.023656) * temperature_100 +
                 FT(0.0047036) * temperature_100^2
             )
    pk1 = FT(3633.86) / temperature_kelvin - FT(61.2172) +
          FT(9.6777) * log_temperature - FT(0.011555) * salinity +
          FT(0.0001152) * salinity^2
    pk2 = FT(471.78) / temperature_kelvin + FT(25.9290) -
          FT(3.16967) * log_temperature - FT(0.01781) * salinity +
          FT(0.0001122) * salinity^2
    log_kb = (
        FT(-8966.90) - FT(2890.53) * root_salinity - FT(77.942) * salinity +
        FT(1.728) * salinity * root_salinity - FT(0.0996) * salinity^2
    ) / temperature_kelvin + FT(148.0248) + FT(137.1942) * root_salinity +
        FT(1.62142) * salinity - (
            FT(24.4344) + FT(25.085) * root_salinity +
            FT(0.2474) * salinity
        ) * log_temperature + FT(0.053105) * root_salinity * temperature_kelvin
    log_kw = FT(148.9652) - FT(13847.26) / temperature_kelvin -
             FT(23.6521) * log_temperature + (
                 FT(-5.977) + FT(118.67) / temperature_kelvin +
                 FT(1.0495) * log_temperature
             ) * root_salinity - FT(0.01615) * salinity
    log_ten = log(FT(10))
    return (
        K0 = exp(log_k0),
        K1 = exp(-log_ten * pk1),
        K2 = exp(-log_ten * pk2),
        KB = exp(log_kb),
        Kw = exp(log_kw),
        total_boron = FT(4.157e-4) * salinity / FT(35),
    )
end

@inline function _carbon_device_alkalinity_residual(
    hydrogen::FT,
    dic_mol_kg::FT,
    alkalinity_mol_kg::FT,
    constants,
) where {FT<:AbstractFloat}
    carbonate = dic_mol_kg * (
        hydrogen * constants.K1 + FT(2) * constants.K1 * constants.K2
    ) / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    )
    borate = constants.total_boron * constants.KB / (constants.KB + hydrogen)
    water = constants.Kw / hydrogen - hydrogen
    return carbonate + borate + water - alkalinity_mol_kg
end

@inline function _carbon_device_input_valid(
    dic_umol_kg::FT,
    alkalinity_umol_kg::FT,
    temperature_kelvin::FT,
    salinity::FT,
) where {FT<:AbstractFloat}
    return isfinite(dic_umol_kg) && dic_umol_kg > zero(FT) &&
           isfinite(alkalinity_umol_kg) && alkalinity_umol_kg > zero(FT) &&
           isfinite(temperature_kelvin) &&
           FT(271.15) <= temperature_kelvin <= FT(313.15) &&
           isfinite(salinity) && FT(5) <= salinity <= FT(42)
end

@inline function _carbon_device_bracket_valid(
    dic_umol_kg::FT,
    alkalinity_umol_kg::FT,
    temperature_kelvin::FT,
    salinity::FT,
) where {FT<:AbstractFloat}
    constants = _carbon_device_constants(temperature_kelvin, salinity)
    dic = dic_umol_kg * FT(1e-6)
    alkalinity = alkalinity_umol_kg * FT(1e-6)
    low_residual = _carbon_device_alkalinity_residual(
        FT(1e-10),
        dic,
        alkalinity,
        constants,
    )
    high_residual = _carbon_device_alkalinity_residual(
        FT(1e-6),
        dic,
        alkalinity,
        constants,
    )
    guard = FT(_CARBON_STAGE0B_BRACKET_RESIDUAL_GUARD) * alkalinity
    return isfinite(low_residual) && isfinite(high_residual) &&
           low_residual > guard && high_residual < -guard
end

@inline function _carbon_device_pco2_uatm_validated(
    dic_umol_kg::FT,
    alkalinity_umol_kg::FT,
    temperature_kelvin::FT,
    salinity::FT,
) where {FT<:AbstractFloat}
    constants = _carbon_device_constants(temperature_kelvin, salinity)
    dic = dic_umol_kg * FT(1e-6)
    alkalinity = alkalinity_umol_kg * FT(1e-6)
    low = FT(1e-10)
    high = FT(1e-6)

    for _ in 1:_CARBON_STAGE0B_ITERATIONS
        midpoint = sqrt(low * high)
        residual = _carbon_device_alkalinity_residual(
            midpoint,
            dic,
            alkalinity,
            constants,
        )
        if residual > zero(FT)
            low = midpoint
        else
            high = midpoint
        end
    end
    hydrogen = sqrt(low * high)
    co2 = dic * hydrogen^2 / (
        hydrogen^2 + hydrogen * constants.K1 + constants.K1 * constants.K2
    )
    pco2 = co2 / constants.K0 * FT(1e6)
    valid_output = isfinite(pco2) && pco2 > zero(FT)
    return (ifelse(valid_output, pco2, zero(FT)), valid_output)
end

"""Device-safe carbonate solve. It returns `(pco2, valid)` and never throws."""
@inline function _carbon_device_pco2_uatm(
    dic_umol_kg::FT,
    alkalinity_umol_kg::FT,
    temperature_kelvin::FT,
    salinity::FT,
) where {FT<:AbstractFloat}
    _carbon_device_input_valid(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    ) || return (zero(FT), false)
    _carbon_device_bracket_valid(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    ) || return (zero(FT), false)
    return _carbon_device_pco2_uatm_validated(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    )
end

@inline function _carbon_schmidt_number_co2(temperature_celsius::FT) where {FT<:AbstractFloat}
    return FT(2116.8) - FT(136.25) * temperature_celsius +
           FT(4.7353) * temperature_celsius^2 -
           FT(0.092307) * temperature_celsius^3 +
           FT(0.0007555) * temperature_celsius^4
end

@inline function _carbon_transfer_velocity_m_s(
    wind_speed_m_s::FT,
    temperature_celsius::FT,
) where {FT<:AbstractFloat}
    schmidt = _carbon_schmidt_number_co2(temperature_celsius)
    return FT(0.251) * wind_speed_m_s^2 * sqrt(FT(660) / schmidt) *
           FT(0.01) / FT(3600)
end

"""Water-vapour pressure over seawater in atmospheres (Weiss and Price, 1980)."""
@inline function _carbon_water_vapour_pressure_atm(
    temperature_kelvin::FT,
    salinity::FT,
) where {FT<:AbstractFloat}
    log_pressure = FT(24.4543) - FT(67.4509) * FT(100) / temperature_kelvin -
                   FT(4.8489) * log(temperature_kelvin / FT(100)) -
                   FT(0.000544) * salinity
    return exp(log_pressure)
end

@inline function _carbon_flux_forcing_valid(
    temperature_kelvin::FT,
    salinity::FT,
    wind_speed_m_s::FT,
    ice_fraction::FT,
    atmospheric_co2_ppm::FT,
    surface_pressure_pa::FT,
    seawater_density::FT,
) where {FT<:AbstractFloat}
    return isfinite(temperature_kelvin) && isfinite(salinity) &&
           isfinite(wind_speed_m_s) && wind_speed_m_s >= zero(FT) &&
           isfinite(ice_fraction) && zero(FT) <= ice_fraction <= one(FT) &&
           isfinite(atmospheric_co2_ppm) && atmospheric_co2_ppm > zero(FT) &&
           isfinite(surface_pressure_pa) && surface_pressure_pa > zero(FT) &&
           isfinite(seawater_density) && seawater_density > zero(FT)
end

@inline function _carbon_air_sea_flux_from_pco2(
    pco2::FT,
    temperature_kelvin::FT,
    salinity::FT,
    wind_speed_m_s::FT,
    ice_fraction::FT,
    atmospheric_co2_ppm::FT,
    surface_pressure_pa::FT,
    seawater_density::FT,
) where {FT<:AbstractFloat}
    constants = _carbon_device_constants(temperature_kelvin, salinity)
    water_vapour_pressure = _carbon_water_vapour_pressure_atm(
        temperature_kelvin,
        salinity,
    )
    dry_air_pressure = surface_pressure_pa / FT(101325) - water_vapour_pressure
    valid_dry_air = isfinite(dry_air_pressure) && dry_air_pressure > zero(FT)
    valid_dry_air || return (
        flux = zero(FT),
        pco2,
        air_pco2 = zero(FT),
        valid = false,
    )
    air_pco2 = atmospheric_co2_ppm * dry_air_pressure
    transfer_velocity = _carbon_transfer_velocity_m_s(
        wind_speed_m_s,
        temperature_kelvin - FT(273.15),
    )
    flux = transfer_velocity * seawater_density * constants.K0 *
           (one(FT) - ice_fraction) *
           (air_pco2 - pco2) * FT(1e-6)
    valid_flux = isfinite(flux)
    return (
        flux = ifelse(valid_flux, flux, zero(FT)),
        pco2,
        air_pco2,
        valid = valid_flux,
    )
end

"""Air-sea carbon flux in mol C m-2 s-1, positive into the ocean."""
@inline function _carbon_air_sea_flux_mol_m2_s(
    dic_umol_kg::FT,
    alkalinity_umol_kg::FT,
    temperature_kelvin::FT,
    salinity::FT,
    wind_speed_m_s::FT,
    ice_fraction::FT,
    atmospheric_co2_ppm::FT,
    surface_pressure_pa::FT,
    seawater_density::FT,
) where {FT<:AbstractFloat}
    valid_forcing = _carbon_flux_forcing_valid(
        temperature_kelvin,
        salinity,
        wind_speed_m_s,
        ice_fraction,
        atmospheric_co2_ppm,
        surface_pressure_pa,
        seawater_density,
    )
    valid_forcing || return (
        flux = zero(FT),
        pco2 = zero(FT),
        air_pco2 = zero(FT),
        valid = false,
    )
    pco2, valid_chemistry = _carbon_device_pco2_uatm(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    )
    valid_chemistry || return (
        flux = zero(FT),
        pco2 = zero(FT),
        air_pco2 = zero(FT),
        valid = false,
    )
    return _carbon_air_sea_flux_from_pco2(
        pco2,
        temperature_kelvin,
        salinity,
        wind_speed_m_s,
        ice_fraction,
        atmospheric_co2_ppm,
        surface_pressure_pa,
        seawater_density,
    )
end

@inline function _carbon_air_sea_flux_mol_m2_s(
    dic_umol_kg::Float32,
    alkalinity_umol_kg::Float32,
    temperature_kelvin::Float32,
    salinity::Float32,
    wind_speed_m_s::Float32,
    ice_fraction::Float32,
    atmospheric_co2_ppm::Float32,
    surface_pressure_pa::Float32,
    seawater_density::Float32,
)
    valid_forcing = _carbon_flux_forcing_valid(
        temperature_kelvin,
        salinity,
        wind_speed_m_s,
        ice_fraction,
        atmospheric_co2_ppm,
        surface_pressure_pa,
        seawater_density,
    )
    valid_chemistry = _carbon_device_input_valid(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    ) && _carbon_device_bracket_valid(
        dic_umol_kg,
        alkalinity_umol_kg,
        temperature_kelvin,
        salinity,
    )
    valid_forcing && valid_chemistry || return (
        flux = 0f0,
        pco2 = 0f0,
        air_pco2 = 0f0,
        valid = false,
    )
    wide_pco2, valid_wide = _carbon_device_pco2_uatm_validated(
        Float64(dic_umol_kg),
        Float64(alkalinity_umol_kg),
        Float64(temperature_kelvin),
        Float64(salinity),
    )
    valid_wide || return (
        flux = 0f0,
        pco2 = 0f0,
        air_pco2 = 0f0,
        valid = false,
    )
    wide = _carbon_air_sea_flux_from_pco2(
        wide_pco2,
        Float64(temperature_kelvin),
        Float64(salinity),
        Float64(wind_speed_m_s),
        Float64(ice_fraction),
        Float64(atmospheric_co2_ppm),
        Float64(surface_pressure_pa),
        Float64(seawater_density),
    )
    return (
        flux = Float32(wide.flux),
        pco2 = Float32(wide.pco2),
        air_pco2 = Float32(wide.air_pco2),
        valid = wide.valid,
    )
end

"""Apply a bounded surface-cell flux and return its updated DIC and PgC uptake."""
@inline function _carbon_apply_surface_cell_flux(
    dic_umol_kg::FT,
    raw_flux_mol_m2_s::FT,
    time_step_seconds::FT,
    cell_area_m2::FT,
    top_cell_thickness_m::FT,
    seawater_density::FT,
) where {FT<:AbstractFloat}
    valid = isfinite(dic_umol_kg) && dic_umol_kg > zero(FT) &&
            isfinite(raw_flux_mol_m2_s) &&
            isfinite(time_step_seconds) && time_step_seconds > zero(FT) &&
            isfinite(cell_area_m2) && cell_area_m2 > zero(FT) &&
            isfinite(top_cell_thickness_m) && top_cell_thickness_m > zero(FT) &&
            isfinite(seawater_density) && seawater_density > zero(FT)
    valid || return (
        dic_umol_kg = dic_umol_kg,
        flux_mol_m2_s = zero(FT),
        uptake_pgc = zero(FT),
        valid = false,
        limited = false,
    )
    maximum_outgassing = FT(0.9) * dic_umol_kg * FT(1e-6) *
                         seawater_density * top_cell_thickness_m /
                         time_step_seconds
    flux = max(raw_flux_mol_m2_s, -maximum_outgassing)
    updated_dic = dic_umol_kg + flux * time_step_seconds /
                  (seawater_density * top_cell_thickness_m) * FT(1e6)
    uptake_pgc = flux * cell_area_m2 * time_step_seconds *
                 FT(_CARBON_MOLAR_MASS_G_MOL) * FT(_CARBON_PG_PER_G)
    valid_output = isfinite(updated_dic) && updated_dic > zero(FT) &&
                   isfinite(uptake_pgc)
    return (
        dic_umol_kg = ifelse(valid_output, updated_dic, dic_umol_kg),
        flux_mol_m2_s = ifelse(valid_output, flux, zero(FT)),
        uptake_pgc = ifelse(valid_output, uptake_pgc, zero(FT)),
        valid = valid_output,
        limited = valid_output && flux != raw_flux_mol_m2_s,
    )
end

function _carbon_fixed_order_sum(values)
    total = 0.0
    compensation = 0.0
    for value in values
        number = Float64(value)
        isfinite(number) || throw(ArgumentError("carbon reduction contains a nonfinite value"))
        next_total = total + number
        if abs(total) >= abs(number)
            compensation += (total - next_total) + number
        else
            compensation += (number - next_total) + total
        end
        total = next_total
    end
    return total + compensation
end

struct _Stage0bCarbonLedger{FT<:AbstractFloat}
    mode::Symbol
    atmospheric_carbon_pgc::FT
    ocean_carbon_pgc::FT
    cumulative_emissions_pgc::FT
    cumulative_air_sea_uptake_pgc::FT
    cumulative_implied_emissions_pgc::FT
end

function _Stage0bCarbonLedger(
    mode::Symbol,
    atmospheric_carbon_pgc::FT,
    ocean_carbon_pgc::FT,
) where {FT<:AbstractFloat}
    mode in (:emissions_driven, :prescribed) || throw(ArgumentError(
        "carbon ledger mode must be emissions_driven or prescribed",
    ))
    atmospheric_carbon_pgc > zero(FT) || throw(ArgumentError(
        "atmospheric carbon must be positive",
    ))
    ocean_carbon_pgc > zero(FT) || throw(ArgumentError(
        "ocean carbon must be positive",
    ))
    return _Stage0bCarbonLedger(
        mode,
        atmospheric_carbon_pgc,
        ocean_carbon_pgc,
        zero(FT),
        zero(FT),
        zero(FT),
    )
end

function _carbon_ledger_invariant(state::_Stage0bCarbonLedger)
    return state.atmospheric_carbon_pgc + state.ocean_carbon_pgc -
           state.cumulative_emissions_pgc -
           state.cumulative_implied_emissions_pgc
end

function _carbon_ledger_step(
    state::_Stage0bCarbonLedger{FT},
    uptake_pgc::Real,
    emissions_pgc::Real = 0,
) where {FT<:AbstractFloat}
    uptake = FT(uptake_pgc)
    emissions = FT(emissions_pgc)
    isfinite(uptake) || throw(ArgumentError("air-sea uptake must be finite"))
    isfinite(emissions) && emissions >= zero(FT) || throw(ArgumentError(
        "emissions must be finite and nonnegative",
    ))
    if state.mode === :prescribed
        emissions == zero(FT) || throw(ArgumentError(
            "prescribed mode records held-atmosphere supply as implied emissions",
        ))
        updated = _Stage0bCarbonLedger(
            state.mode,
            state.atmospheric_carbon_pgc,
            state.ocean_carbon_pgc + uptake,
            state.cumulative_emissions_pgc,
            state.cumulative_air_sea_uptake_pgc + uptake,
            state.cumulative_implied_emissions_pgc + uptake,
        )
    else
        new_atmosphere = state.atmospheric_carbon_pgc + emissions - uptake
        new_atmosphere > zero(FT) || throw(ArgumentError(
            "air-sea uptake and emissions would exhaust atmospheric carbon",
        ))
        updated = _Stage0bCarbonLedger(
            state.mode,
            new_atmosphere,
            state.ocean_carbon_pgc + uptake,
            state.cumulative_emissions_pgc + emissions,
            state.cumulative_air_sea_uptake_pgc + uptake,
            state.cumulative_implied_emissions_pgc,
        )
    end
    updated.ocean_carbon_pgc > zero(FT) || throw(ArgumentError(
        "air-sea flux would exhaust ocean carbon",
    ))
    return updated
end

@kernel function _carbon_pco2_kernel!(
    output_pco2,
    output_valid,
    dic,
    alkalinity,
    temperature_kelvin,
    salinity,
)
    index = @index(Global)
    if index <= length(output_pco2)
        pco2, valid = _carbon_device_pco2_uatm(
            dic[index],
            alkalinity[index],
            temperature_kelvin[index],
            salinity[index],
        )
        output_pco2[index] = pco2
        output_valid[index] = ifelse(valid, UInt8(1), UInt8(0))
    end
end

@kernel function _carbon_air_sea_flux_kernel!(
    output_flux,
    output_pco2,
    output_air_pco2,
    output_valid,
    dic,
    alkalinity,
    temperature_kelvin,
    salinity,
    wind_speed_m_s,
    ice_fraction,
    atmospheric_co2_ppm,
    surface_pressure_pa,
    seawater_density,
)
    index = @index(Global)
    if index <= length(output_flux)
        result = _carbon_air_sea_flux_mol_m2_s(
            dic[index],
            alkalinity[index],
            temperature_kelvin[index],
            salinity[index],
            wind_speed_m_s[index],
            ice_fraction[index],
            atmospheric_co2_ppm[index],
            surface_pressure_pa[index],
            seawater_density[index],
        )
        output_flux[index] = result.flux
        output_pco2[index] = result.pco2
        output_air_pco2[index] = result.air_pco2
        output_valid[index] = ifelse(result.valid, UInt8(1), UInt8(0))
    end
end

@kernel function _carbon_surface_cell_flux_kernel!(
    output_dic,
    output_flux,
    output_uptake_pgc,
    output_valid,
    output_limited,
    dic,
    raw_flux,
    time_step_seconds,
    cell_area_m2,
    top_cell_thickness_m,
    seawater_density,
)
    index = @index(Global)
    if index <= length(output_dic)
        result = _carbon_apply_surface_cell_flux(
            dic[index],
            raw_flux[index],
            time_step_seconds[index],
            cell_area_m2[index],
            top_cell_thickness_m[index],
            seawater_density[index],
        )
        output_dic[index] = result.dic_umol_kg
        output_flux[index] = result.flux_mol_m2_s
        output_uptake_pgc[index] = result.uptake_pgc
        output_valid[index] = ifelse(result.valid, UInt8(1), UInt8(0))
        output_limited[index] = ifelse(result.limited, UInt8(1), UInt8(0))
    end
end
