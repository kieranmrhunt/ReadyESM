#=
Time-varying well-mixed greenhouse gases from a sealed annual input file.

Production prescribes constant gases for the whole integration: `forcing.co2_ppm`, optional constant CH4 and N2O, and
no halocarbons. With `forcing.well_mixed_gas_path` set, RRTMGP instead receives CO2, CH4, N2O, CFC-12 and CFC-11 that
follow the model clock. The file holds global annual means (schema `readiesm_well_mixed_gases_v1`) placed at mid-year
(decimal year y + 0.5), and the schedule interpolates linearly between them. There is no extrapolation. The
halocarbon representation is the CMIP option of Meinshausen et al. (2017): actual CFC-12 plus the CFC-11 equivalent
of every other halocarbon, applied to the RRTMGP gas `cfc11`. `forcing.well_mixed_gas_fixed_year > 0` holds that
year's values for the whole run (piControl or a time-slice experiment).
=#

const WELL_MIXED_GAS_SCHEMA = "readiesm_well_mixed_gases_v1"
const WELL_MIXED_GAS_NAMES = (:co2, :ch4, :n2o, :cfc12, :cfc11eq)
const WELL_MIXED_GAS_RRTMGP_NAMES = (co2 = "co2", ch4 = "ch4", n2o = "n2o", cfc12 = "cfc12", cfc11eq = "cfc11")
const WELL_MIXED_GAS_TIME_DEPENDENCE = "well_mixed_greenhouse_gases_interpolated_from_annual_file"

"""Annual global-mean well-mixed gas mole fractions (mol mol-1), one row per year, read and checked from `path`."""
struct WellMixedGasSchedule
    path::String
    sha256::String
    scenario::String
    years::Vector{Int}
    values::Matrix{Float64}
    fixed_year::Int
end

function WellMixedGasSchedule(path::AbstractString; fixed_year::Integer = 0)
    isfile(path) || throw(ArgumentError("well-mixed gas file not found: $path"))
    fixed_year >= 0 || throw(ArgumentError("well_mixed_gas_fixed_year must be zero or a positive year"))
    years, values, scenario = NCDatasets.NCDataset(path) do ds
        get(ds.attrib, "schema", "") == WELL_MIXED_GAS_SCHEMA ||
            throw(ArgumentError("well-mixed gas file $path does not have schema $WELL_MIXED_GAS_SCHEMA"))
        years = Int.(Array(ds["year"]))
        values = Matrix{Float64}(undef, length(years), length(WELL_MIXED_GAS_NAMES))
        for (j, name) in enumerate(WELL_MIXED_GAS_NAMES)
            haskey(ds, String(name)) || throw(ArgumentError("well-mixed gas file $path has no variable $name"))
            variable = ds[String(name)]
            get(variable.attrib, "units", "") == "mol mol-1" ||
                throw(ArgumentError("well-mixed gas $name in $path must be in mol mol-1"))
            get(variable.attrib, "rrtmgp_gas", "") == WELL_MIXED_GAS_RRTMGP_NAMES[name] ||
                throw(ArgumentError("well-mixed gas $name in $path maps to the wrong RRTMGP gas"))
            values[:, j] .= Float64.(Array(variable))
        end
        years, values, String(get(ds.attrib, "scenario", "unknown"))
    end
    length(years) >= 2 && all(diff(years) .== 1) ||
        throw(ArgumentError("well-mixed gas file $path must have at least two consecutive years"))
    all(isfinite, values) && all(>=(0), values) ||
        throw(ArgumentError("well-mixed gas file $path has a negative or nonfinite mole fraction"))
    all(>(0), values[:, 1:3]) || throw(ArgumentError("well-mixed gas file $path has zero CO2, CH4 or N2O"))
    fixed_year == 0 || first(years) <= fixed_year <= last(years) ||
        throw(ArgumentError("well_mixed_gas_fixed_year $fixed_year is outside $(first(years))-$(last(years))"))
    return WellMixedGasSchedule(String(path), bytes2hex(open(SHA.sha256, path)), scenario, years, values,
                                Int(fixed_year))
end

"""Calendar year plus the elapsed fraction of that year (proleptic Gregorian, millisecond resolution)."""
function _decimal_year(time::DateTime)
    year = Dates.year(time)
    start = DateTime(year)
    return year + Dates.value(time - start) / Dates.value(DateTime(year + 1) - start)
end

"""Time span the schedule can interpolate without extrapolation, as decimal years."""
_well_mixed_gas_span(schedule::WellMixedGasSchedule) = (first(schedule.years) + 0.5, last(schedule.years) + 0.5)

"""Gas mole fractions at `time` as a NamedTuple keyed by `WELL_MIXED_GAS_NAMES`."""
function well_mixed_gas_values(schedule::WellMixedGasSchedule, time::DateTime)
    if schedule.fixed_year > 0
        row = schedule.fixed_year - first(schedule.years) + 1
        return NamedTuple{WELL_MIXED_GAS_NAMES}(Tuple(schedule.values[row, :]))
    end
    x = _decimal_year(time)
    lower, upper = _well_mixed_gas_span(schedule)
    lower <= x <= upper || throw(ArgumentError(
        "model time $time (decimal year $x) is outside the well-mixed gas schedule $(lower)-$(upper)"))
    u = x - 0.5
    row = floor(Int, u) - first(schedule.years) + 1
    weight = u - floor(u)
    row == length(schedule.years) && return NamedTuple{WELL_MIXED_GAS_NAMES}(Tuple(schedule.values[row, :]))
    values = ntuple(j -> (1 - weight) * schedule.values[row, j] + weight * schedule.values[row + 1, j],
                    length(WELL_MIXED_GAS_NAMES))
    return NamedTuple{WELL_MIXED_GAS_NAMES}(values)
end

"""The schedule identity recorded in checkpoints and diagnostics provenance."""
_well_mixed_gas_identity(schedule::WellMixedGasSchedule) = (sha256 = schedule.sha256, fixed_year = schedule.fixed_year)
_well_mixed_gas_identity(::Nothing) = nothing

"""The schedule for `config`, or `nothing` for production constant gases."""
_well_mixed_gas_schedule(config) = isempty(config.forcing.well_mixed_gas_path) ? nothing :
    WellMixedGasSchedule(config.forcing.well_mixed_gas_path; fixed_year = config.forcing.well_mixed_gas_fixed_year)

"""Set the radiation's gas fields from the schedule at `time`. A no-op without a schedule."""
_apply_well_mixed_gas_schedule!(radiation, time) = _apply_well_mixed_gas_schedule!(radiation, radiation.gas_schedule, time)
_apply_well_mixed_gas_schedule!(radiation, ::Nothing, time) = nothing

function _apply_well_mixed_gas_schedule!(radiation, schedule::WellMixedGasSchedule, time)
    values = well_mixed_gas_values(schedule, time)
    NF = typeof(radiation.co2_ppm)
    radiation.co2_ppm = NF(values.co2 * 1.0e6)
    radiation.ch4_vmr = NF(values.ch4)
    radiation.n2o_vmr = NF(values.n2o)
    radiation.cfc12_vmr = NF(values.cfc12)
    radiation.cfc11_vmr = NF(values.cfc11eq)
    return values
end

"""Write the optional halocarbons into the solver; the production solver keeps RRTMGP's zero halocarbons."""
function _set_halocarbon_mixing_ratios!(solver, radiation)
    isnothing(radiation.cfc12_vmr) || RRTMGP.set_volume_mixing_ratio!(solver, "cfc12", radiation.cfc12_vmr)
    isnothing(radiation.cfc11_vmr) || RRTMGP.set_volume_mixing_ratio!(solver, "cfc11", radiation.cfc11_vmr)
    return nothing
end

"""
Reject a checkpoint whose well-mixed gas schedule differs from the reconstructed radiation, before any state changes.
Gas values are then recomputed from the restored clock at the next radiation solve. A checkpoint written before this
option existed has no record and restores only into a configuration without a schedule.
"""
function _validate_well_mixed_gas_restart(radiation, state)
    saved = hasproperty(state, :well_mixed_gas_schedule) ? state.well_mixed_gas_schedule : nothing
    current = _well_mixed_gas_identity(radiation.gas_schedule)
    saved == current || throw(ArgumentError(
        "checkpoint and reconstructed radiation use different well-mixed gas schedules " *
        "(checkpoint $(something(saved, "constant gases")), configuration $(something(current, "constant gases")))",
    ))
    return nothing
end

"""
Current scheduled gas values for the restart record, only when a schedule is active, so that a restored radiation
state equals the saved one exactly. The values are recomputed from the clock at the next radiation solve, so they never
affect physics, but an exact fresh-process restart must still reproduce them. Without a schedule the record is absent
and constant-gas checkpoints are unchanged.
"""
_well_mixed_gas_values_record(radiation) = _well_mixed_gas_values_record(radiation, radiation.gas_schedule)
_well_mixed_gas_values_record(radiation, ::Nothing) = nothing
_well_mixed_gas_values_record(radiation, ::WellMixedGasSchedule) = (; co2_ppm = radiation.co2_ppm,
    ch4_vmr = radiation.ch4_vmr, n2o_vmr = radiation.n2o_vmr, cfc12_vmr = radiation.cfc12_vmr,
    cfc11_vmr = radiation.cfc11_vmr)

"""Restore the scheduled gas values saved with a schedule-active checkpoint (after the identity check)."""
function _restore_well_mixed_gas_values!(radiation, state)
    isnothing(radiation.gas_schedule) && return nothing
    hasproperty(state, :well_mixed_gas_values) || throw(ArgumentError(
        "checkpoint uses a well-mixed gas schedule but has no saved gas values"))
    saved = state.well_mixed_gas_values
    NF = typeof(radiation.co2_ppm)
    radiation.co2_ppm = NF(saved.co2_ppm)
    radiation.ch4_vmr = NF(saved.ch4_vmr)
    radiation.n2o_vmr = NF(saved.n2o_vmr)
    radiation.cfc12_vmr = NF(saved.cfc12_vmr)
    radiation.cfc11_vmr = NF(saved.cfc11_vmr)
    return nothing
end

"""Diagnostics provenance of the schedule, written only when it is active."""
function _well_mixed_gas_provenance(radiation, config, rrtmgp_active)
    isempty(config.forcing.well_mixed_gas_path) && return (;)
    schedule = rrtmgp_active ? radiation.gas_schedule : nothing
    isnothing(schedule) && return (well_mixed_gas_forcing = "requested_but_not_active",)
    start_values = well_mixed_gas_values(schedule, config.start_date)
    return (
        well_mixed_gas_forcing = WELL_MIXED_GAS_TIME_DEPENDENCE,
        well_mixed_gas_path = schedule.path,
        well_mixed_gas_sha256 = schedule.sha256,
        well_mixed_gas_scenario = schedule.scenario,
        well_mixed_gas_fixed_year = schedule.fixed_year,
        well_mixed_gas_halocarbons = "cfc12_plus_cfc11eq_meinshausen2017",
        well_mixed_gas_start_co2_ppm = start_values.co2 * 1.0e6,
        well_mixed_gas_start_ch4_ppb = start_values.ch4 * 1.0e9,
        well_mixed_gas_start_n2o_ppb = start_values.n2o * 1.0e9,
        well_mixed_gas_start_cfc12_ppt = start_values.cfc12 * 1.0e12,
        well_mixed_gas_start_cfc11eq_ppt = start_values.cfc11eq * 1.0e12,
        well_mixed_gas_final_co2_ppm = Float64(radiation.co2_ppm),
        well_mixed_gas_final_ch4_ppb = Float64(radiation.ch4_vmr) * 1.0e9,
        well_mixed_gas_final_n2o_ppb = Float64(radiation.n2o_vmr) * 1.0e9,
        well_mixed_gas_final_cfc12_ppt = Float64(radiation.cfc12_vmr) * 1.0e12,
        well_mixed_gas_final_cfc11eq_ppt = Float64(radiation.cfc11_vmr) * 1.0e12,
        requested_co2_ppm_ignored = true,
    )
end
