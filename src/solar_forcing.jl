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

"""Time dependence of the radiative forcing, as recorded in diagnostics provenance."""
function _forcing_time_dependence(config)
    base = _forcing_time_dependence_without_ozone(config)
    isempty(config.forcing.ozone_path) && return base
    ozone = config.forcing.ozone_fixed_year > 0 ?
        "zonal_mean_ozone_fixed_at_year_$(config.forcing.ozone_fixed_year)" : PRESCRIBED_OZONE_TIME_DEPENDENCE
    return base == "constant_for_entire_integration" ? ozone : base * "_and_" * ozone
end

function _forcing_time_dependence_without_ozone(config)
    carbon = config.carbon.mode == :emissions_driven
    solar = config.forcing.solar_distance_correction
    if !isempty(config.forcing.well_mixed_gas_path)
        gases = config.forcing.well_mixed_gas_fixed_year > 0 ?
            "well_mixed_greenhouse_gases_fixed_at_year_$(config.forcing.well_mixed_gas_fixed_year)" :
            WELL_MIXED_GAS_TIME_DEPENDENCE
        return solar ? gases * "_and_seasonal_solar_distance_from_model_clock" : gases
    end
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

# Total solar irradiance (design: artifacts/development_wave_20261101/tsi_option_claude/DESIGN.md v2). ReadyESM used to
# leave the SpeedyWeather planet at its default, so every earlier release ran at 1365 W/m2. This helper builds the planet
# with the configured value and is passed to every PrimitiveWetModel in src/baseline.jl, since both SpeedyWeather and
# RRTMGP read planet.solar_constant. At the default the planet is identical to SpeedyWeather's own Earth(spectral_grid).
const DEFAULT_SOLAR_CONSTANT_W_M2 = 1365.0

_readyesm_planet(config, spectral_grid) =
    SpeedyWeather.Earth(spectral_grid; solar_constant = spectral_grid.NF(config.forcing.solar_constant_w_m2))

_solar_constant_source(config) = isempty(strip(config.forcing.solar_constant_source)) ?
    "SpeedyWeather Earth default" : config.forcing.solar_constant_source

# The pinned values of data/cmip6_solar_tsi/tsi_values.toml (design v2 amendment 2), held here because ReadyESM has no
# TOML dependency. test/solar_constant_option.jl checks the table's SHA-256 and that each value and SHA-256 below appears
# in it. A non-default solar constant must name one entry, so its value is traceable to a CEDA file and a reduction.
const SOLAR_CONSTANT_TABLE = "data/cmip6_solar_tsi/tsi_values.toml"
const SOLAR_CONSTANT_TABLE_SHA256 = "417197507c93b03f1df4229354e5480baa300f7f219ff1a96a148dbf54b9f5b0"
const SOLAR_CONSTANT_ENTRIES = Dict(
    "picontrol" => (
        value_w_m2 = 1360.7470703125,
        label = "CMIP6 piControl (Matthes et al. 2017)",
        archive_path = "/badc/cmip6/data/input4MIPs/CMIP6/CMIP/SOLARIS-HEPPA/SOLARIS-HEPPA-3-2/atmos/fx/" *
                       "solarforcing-picontrol-fx/gn/v20170103/solarforcing-picontrol-fx_input4MIPs_solar_CMIP_" *
                       "SOLARIS-HEPPA-3-2_gn_18500101-18730128.nc",
        archive_sha256 = "2b53e8bf7aebde3deaf3a9174ae3f44a6ca3714ff369790019760ca55e21fd48",
        year = "",
        reduction = "the single tsi value of the piControl fx file",
    ),
    "historical_1993" => (
        value_w_m2 = 1361.0793891802225,
        label = "CMIP6 historical annual mean 1993 (SOLARIS-HEPPA v3.2), a constant historical-annual-mean sensitivity",
        archive_path = "/badc/cmip6/data/input4MIPs/CMIP6/CMIP/SOLARIS-HEPPA/SOLARIS-HEPPA-3-2/atmos/mon/" *
                       "solarforcing-ref-mon/gn/v20170103/solarforcing-ref-mon_input4MIPs_solar_CMIP_" *
                       "SOLARIS-HEPPA-3-2_gn_185001-229912.nc",
        archive_sha256 = "4b844a3324b74be16faf7272b64ed7853ea2c606c4db841a318366ad223ebcac",
        year = "1993",
        reduction = "day-weighted mean of the twelve 1993 monthly values (365-day year)",
    ),
)
# A typed value agrees with its entry if it matches to the table's quoted precision (three decimals).
const SOLAR_CONSTANT_ENTRY_TOLERANCE_W_M2 = 1.0e-3

function _validate_solar_constant_entry(forcing)
    value, entry = forcing.solar_constant_w_m2, forcing.solar_constant_entry
    if value == DEFAULT_SOLAR_CONSTANT_W_M2
        isempty(entry) || throw(ArgumentError(
            "forcing.solar_constant_entry = \"$entry\" needs its own forcing.solar_constant_w_m2, not the default",
        ))
        # Otherwise the output would record a free source label for the default value (Barry, port read note 2).
        isempty(strip(forcing.solar_constant_source)) || throw(ArgumentError(
            "forcing.solar_constant_source = \"$(forcing.solar_constant_source)\" needs a pinned entry and its value; " *
            "leave it empty for the default",
        ))
        return nothing
    end
    haskey(SOLAR_CONSTANT_ENTRIES, entry) || throw(ArgumentError(
        "a forcing.solar_constant_w_m2 other than $DEFAULT_SOLAR_CONSTANT_W_M2 needs forcing.solar_constant_entry, " *
        "one of $(sort(collect(keys(SOLAR_CONSTANT_ENTRIES))))",
    ))
    pinned = SOLAR_CONSTANT_ENTRIES[entry]
    abs(value - pinned.value_w_m2) <= SOLAR_CONSTANT_ENTRY_TOLERANCE_W_M2 || throw(ArgumentError(
        "forcing.solar_constant_w_m2 = $value does not match entry \"$entry\" ($(pinned.value_w_m2) W/m2)",
    ))
    forcing.solar_constant_source == pinned.label || throw(ArgumentError(
        "forcing.solar_constant_source must be the label of entry \"$entry\": \"$(pinned.label)\"",
    ))
    return nothing
end

"""
Output check: a non-default solar constant must carry the provenance of a pinned entry it agrees with, and the default
must carry none.
"""
function _validate_solar_constant_provenance(metadata)
    value = metadata.solar_constant_w_m2
    if value == DEFAULT_SOLAR_CONSTANT_W_M2
        hasproperty(metadata, :solar_constant_entry) && error("default solar constant output carries an entry")
        metadata.solar_constant_source == "SpeedyWeather Earth default" ||
            error("default solar constant output carries the source \"$(metadata.solar_constant_source)\"")
        return nothing
    end
    hasproperty(metadata, :solar_constant_entry) && haskey(SOLAR_CONSTANT_ENTRIES, metadata.solar_constant_entry) ||
        error("solar constant $value W/m2 has no pinned entry in its provenance")
    pinned = SOLAR_CONSTANT_ENTRIES[metadata.solar_constant_entry]
    abs(value - pinned.value_w_m2) <= SOLAR_CONSTANT_ENTRY_TOLERANCE_W_M2 &&
        metadata.solar_constant_source == pinned.label &&
        metadata.solar_constant_table_sha256 == SOLAR_CONSTANT_TABLE_SHA256 &&
        metadata.solar_constant_archive_sha256 == pinned.archive_sha256 &&
        metadata.solar_constant_archive_path == pinned.archive_path &&
        metadata.solar_constant_year == pinned.year &&
        metadata.solar_constant_reduction == pinned.reduction ||
        error("solar constant provenance does not match entry \"$(metadata.solar_constant_entry)\"")
    return nothing
end

"Provenance of the applied solar constant: the entry's source file, SHA-256, year and reduction, or the default."
function _solar_constant_provenance(config)
    entry = config.forcing.solar_constant_entry
    isempty(entry) && return (; solar_constant_entry = "", solar_constant_table_sha256 = "",
        solar_constant_archive_path = "", solar_constant_archive_sha256 = "", solar_constant_year = "",
        solar_constant_reduction = "SpeedyWeather Earth default, no source file")
    pinned = SOLAR_CONSTANT_ENTRIES[entry]
    return (; solar_constant_entry = entry, solar_constant_table_sha256 = SOLAR_CONSTANT_TABLE_SHA256,
        solar_constant_archive_path = pinned.archive_path, solar_constant_archive_sha256 = pinned.archive_sha256,
        solar_constant_year = pinned.year, solar_constant_reduction = pinned.reduction)
end
