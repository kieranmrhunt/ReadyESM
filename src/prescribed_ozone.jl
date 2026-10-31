#=
Prescribed zonal-mean ozone from a monthly input file (SSP step 2, opt-in).

Production gives RRTMGP one analytic ozone profile for every column and time:
o3 = 3e-8 + 7.5e-6 exp(-(ln(p / 1200 Pa))^2 / 2.88). With `forcing.ozone_path` set, RRTMGP instead receives the
CMIP ozone of a prepared file (schema `readiesm_zonal_ozone_v1`, written by ssp_forcing_v1/ozone/prepare_ozone.py):
monthly zonal means on the source latitudes and pressure levels. At each radiation update the profiles are
interpolated linearly in time between the monthly midpoints at the model clock, then linearly in latitude to each
model column. The pack loop (CPU) or kernel (GPU) then interpolates linearly in ln p to each layer pressure, with the
lowest or highest source level held outside the source pressure range.

Time follows `_decimal_year` (proleptic Gregorian) and the file's midpoints are decimal years in the source calendar,
so the two calendars can differ by up to about a day. Inside the file's coverage bounds and outside its first and last
midpoints, the first or last month is held. Outside the coverage bounds the run is refused. `ozone_fixed_year > 0`
repeats the twelve months of that year, interpolating cyclically from December to January (piControl or time slices).
=#

const PRESCRIBED_OZONE_SCHEMA = "readiesm_zonal_ozone_v1"
const PRESCRIBED_OZONE_TIME_DEPENDENCE = "zonal_mean_ozone_interpolated_from_monthly_file"

"""Ozone profiles at the current time: one column of source levels per model column, on the model architecture."""
struct OzoneColumns{V, M}
    log_pressure::V  # ln(p / Pa) of the source levels, strictly decreasing (lowest level first)
    values::M        # (source level, model column), mol mol-1
end

Adapt.@adapt_structure OzoneColumns

"""The host-side contents of a prepared ozone file, read and checked."""
struct PrescribedOzoneData
    path::String
    sha256::String
    family::String
    scenario::String
    time::Vector{Float64}       # monthly midpoints, decimal years in the source calendar
    lower::Float64              # coverage start (first month's lower bound)
    upper::Float64              # coverage end (last month's upper bound)
    years::Vector{Int}
    months::Vector{Int}
    pressure::Vector{Float64}   # Pa, strictly decreasing
    latitude::Vector{Float64}   # degrees, strictly increasing
    values::Array{Float64, 3}   # (level, latitude, time)
end

function PrescribedOzoneData(path::AbstractString)
    isfile(path) || throw(ArgumentError("ozone file not found: $path"))
    data = NCDatasets.NCDataset(path) do ds
        get(ds.attrib, "schema", "") == PRESCRIBED_OZONE_SCHEMA ||
            throw(ArgumentError("ozone file $path does not have schema $PRESCRIBED_OZONE_SCHEMA"))
        for name in ("time", "time_bnds", "year", "month", "plev", "lat", "vmro3")
            haskey(ds, name) || throw(ArgumentError("ozone file $path has no variable $name"))
        end
        get(ds["vmro3"].attrib, "units", "") == "mole mole-1" ||
            throw(ArgumentError("ozone in $path must be in mole mole-1"))
        get(ds["plev"].attrib, "units", "") == "Pa" || throw(ArgumentError("ozone pressure in $path must be in Pa"))
        dimnames = NCDatasets.dimnames(ds["vmro3"])
        dimnames == ("lat", "plev", "time") ||
            throw(ArgumentError("ozone vmro3 in $path must have dimensions (lat, plev, time), found $dimnames"))
        # NCDatasets reverses the NetCDF dimension order: (time, plev, lat) on disk is (lat, plev, time) here.
        values = permutedims(Float64.(Array(ds["vmro3"].var[:, :, :])), (2, 1, 3))
        bounds = Float64.(Array(ds["time_bnds"].var[:, :]))
        (time = Float64.(Array(ds["time"].var[:])), lower = bounds[1, 1], upper = bounds[2, end],
         bounds, years = Int.(Array(ds["year"].var[:])), months = Int.(Array(ds["month"].var[:])),
         pressure = Float64.(Array(ds["plev"].var[:])), latitude = Float64.(Array(ds["lat"].var[:])), values,
         family = String(get(ds.attrib, "family", "unknown")), scenario = String(get(ds.attrib, "scenario", "")))
    end
    n = length(data.time)
    n >= 2 && size(data.values, 3) == n && length(data.years) == n && length(data.months) == n ||
        throw(ArgumentError("ozone file $path needs at least two monthly samples with matching year/month"))
    all(isfinite, data.time) && all(isfinite, data.bounds) ||
        throw(ArgumentError("ozone file $path has nonfinite times or time bounds"))
    all(diff(data.time) .> 0) && all(data.bounds[1, 2:end] .== data.bounds[2, 1:(end - 1)]) &&
        all(data.bounds[1, :] .< data.time .< data.bounds[2, :]) ||
        throw(ArgumentError("ozone file $path has non-increasing or non-contiguous time bounds"))
    _validate_ozone_month_labels(path, data.years, data.months, data.time, data.bounds)
    length(data.latitude) >= 2 ||
        throw(ArgumentError("ozone file $path needs at least two latitudes for the latitude interpolation"))
    all(isfinite, data.pressure) && all(>(0), data.pressure) && all(diff(data.pressure) .< 0) ||
        throw(ArgumentError("ozone pressure levels in $path must be positive and strictly decreasing"))
    all(isfinite, data.latitude) && all(diff(data.latitude) .> 0) && all(abs.(data.latitude) .<= 90) ||
        throw(ArgumentError("ozone latitudes in $path must be strictly increasing within ±90"))
    all(isfinite, data.values) && all(>=(0), data.values) ||
        throw(ArgumentError("ozone file $path has a negative or nonfinite mixing ratio"))
    return PrescribedOzoneData(String(path), bytes2hex(open(SHA.sha256, path)), data.family, data.scenario,
                               data.time, data.lower, data.upper, data.years, data.months, data.pressure,
                               data.latitude, data.values)
end

"""
Check the year/month labels against the time axis without inferring calendar months from day counts: labels are
unique and advance by exactly one month, each month's bounds lie inside its labelled year, and the midpoint of
month m lies strictly inside ((m - 1)/12, m/12) of that year. This holds for 360-, 365- and 366-day years (the
midpoint of each month sits at least 0.037 yr inside these limits), and it rejects labels shifted by a month or more.
"""
function _validate_ozone_month_labels(path, years, months, time, bounds)
    all(m -> 1 <= m <= 12, months) || throw(ArgumentError("ozone file $path has month labels outside 1-12"))
    for i in 2:length(years)
        (years[i] * 12 + months[i]) - (years[i - 1] * 12 + months[i - 1]) == 1 ||
            throw(ArgumentError("ozone file $path has non-contiguous or repeated month labels at sample $i"))
    end
    for i in eachindex(years)
        y, m = years[i], months[i]
        f = time[i] - y
        y - 1e-9 <= bounds[1, i] && bounds[2, i] <= y + 1 + 1e-9 && (m - 1) / 12 < f < m / 12 ||
            throw(ArgumentError("ozone file $path labels sample $i as $y-$m but its midpoint is $(time[i]) " *
                                "and its bounds $(bounds[1, i])-$(bounds[2, i])"))
    end
    return nothing
end

"""Indices of the twelve months of `year` (January first), or an error if the file lacks any of them."""
function _ozone_year_indices(data::PrescribedOzoneData, year::Integer)
    indices = [findfirst(i -> data.years[i] == year && data.months[i] == m, eachindex(data.years)) for m in 1:12]
    any(isnothing, indices) && throw(ArgumentError("ozone file $(data.path) lacks some months of fixed year $year"))
    phase = [data.time[i] - year for i in indices]
    all(0 .< phase .< 1) && all(diff(phase) .> 0) ||
        throw(ArgumentError("ozone file $(data.path) has fixed-year phases outside (0, 1) or out of order"))
    return Int.(indices)
end

"""Check that `time` lies inside the coverage of a time-varying file, or that a fixed year is complete."""
function _validate_ozone_time(data::PrescribedOzoneData, fixed_year, start::DateTime, stop::DateTime)
    fixed_year > 0 && return (_ozone_year_indices(data, fixed_year); nothing)
    data.lower <= _decimal_year(start) && _decimal_year(stop) <= data.upper || throw(ArgumentError(
        "the run $start to $stop leaves the ozone file coverage $(data.lower)-$(data.upper) ($(data.path))"))
    return nothing
end

"""Bracketing samples (i, j) and weight w of sample j for `time` (decimal-year interpolation, see header)."""
function _ozone_time_weights(data::PrescribedOzoneData, fixed_year::Int, time::DateTime)
    x = _decimal_year(time)
    if fixed_year > 0
        idx = _ozone_year_indices(data, fixed_year)
        phase = [data.time[i] - fixed_year for i in idx]  # midpoint fractions of the year, increasing
        f = x - floor(x)
        f < phase[1] && return (idx[12], idx[1], (f - (phase[12] - 1)) / (phase[1] - (phase[12] - 1)))
        f >= phase[12] && return (idx[12], idx[1], (f - phase[12]) / (phase[1] + 1 - phase[12]))
        m = searchsortedlast(phase, f)
        return (idx[m], idx[m + 1], (f - phase[m]) / (phase[m + 1] - phase[m]))
    end
    data.lower <= x <= data.upper || throw(ArgumentError(
        "model time $time (decimal year $x) is outside the ozone file coverage $(data.lower)-$(data.upper)"))
    x <= data.time[1] && return (1, 1, 0.0)
    x >= data.time[end] && return (length(data.time), length(data.time), 0.0)
    i = searchsortedlast(data.time, x)
    return (i, i + 1, (x - data.time[i]) / (data.time[i + 1] - data.time[i]))
end

"""Prescribed ozone for one model grid: the file data, per-column latitude weights and the device profiles."""
struct PrescribedOzone{NF, C}
    data::PrescribedOzoneData
    fixed_year::Int
    latitude_index::Vector{Int}       # lower source latitude per model column
    latitude_weight::Vector{Float64}  # weight of the upper source latitude
    host_values::Matrix{NF}           # (source level, model column) staging
    columns::C                        # OzoneColumns on the model architecture
    point_weights::Vector{Float64}    # model area weights, for the column-ozone provenance
end

function PrescribedOzone(data::PrescribedOzoneData, spectral_grid::SpectralGrid; fixed_year::Integer = 0)
    fixed_year >= 0 || throw(ArgumentError("ozone_fixed_year must be zero or a positive year"))
    fixed_year > 0 && _ozone_year_indices(data, fixed_year)
    NF = spectral_grid.NF
    _, latitude = RG.get_londlatds(spectral_grid.grid)
    lat = data.latitude
    index = Vector{Int}(undef, length(latitude))
    weight = Vector{Float64}(undef, length(latitude))
    for (c, φ) in enumerate(Float64.(latitude))
        φ = clamp(φ, lat[1], lat[end])  # source latitudes reach ±90 in the CMIP files
        l = clamp(searchsortedlast(lat, φ), 1, length(lat) - 1)
        index[c] = l
        weight[c] = (φ - lat[l]) / (lat[l + 1] - lat[l])
    end
    host_values = zeros(NF, length(data.pressure), length(latitude))
    architecture = spectral_grid.architecture
    columns = OzoneColumns(SpeedyWeather.on_architecture(architecture, NF.(log.(data.pressure))),
                           SpeedyWeather.on_architecture(architecture, copy(host_values)))
    point_weights = Float64.(Array(_global_point_weights(spectral_grid)))
    return PrescribedOzone{NF, typeof(columns)}(data, Int(fixed_year), index, weight, host_values, columns,
                                                point_weights)
end

PrescribedOzone(path::AbstractString, spectral_grid::SpectralGrid; fixed_year::Integer = 0) =
    PrescribedOzone(PrescribedOzoneData(path), spectral_grid; fixed_year)

"""Fill the device profiles for `time`. A no-op without prescribed ozone."""
_update_prescribed_ozone!(::Nothing, time) = nothing

function _update_prescribed_ozone!(ozone::PrescribedOzone{NF}, time) where {NF}
    i, j, w = _ozone_time_weights(ozone.data, ozone.fixed_year, time)
    v = ozone.data.values
    for c in eachindex(ozone.latitude_index)
        l, a = ozone.latitude_index[c], ozone.latitude_weight[c]
        for k in axes(v, 1)
            # a + w (b - a) form: exact when the two values are equal (constant fields stay constant).
            at_i = v[k, l, i] + a * (v[k, l + 1, i] - v[k, l, i])
            at_j = v[k, l, j] + a * (v[k, l + 1, j] - v[k, l, j])
            ozone.host_values[k, c] = NF(at_i + w * (at_j - at_i))
        end
    end
    copyto!(ozone.columns.values, ozone.host_values)
    return ozone
end

"""Device profiles for the pack loop or kernel, or `nothing` for the analytic production profile."""
_ozone_columns(::Nothing) = nothing
_ozone_columns(ozone::PrescribedOzone) = ozone.columns

"""
Layer ozone at pressure `p` (Pa) for column `ij`. Without prescribed ozone: the analytic production profile, written
exactly as before (constants in the surface-pressure type `typeof(ps)`), so the off path is unchanged.
"""
@inline _layer_ozone(::Nothing, ij, p, ps) =
    typeof(ps)(3.0e-8) + typeof(ps)(7.5e-6) *
    exp(-(log(p / typeof(ps)(1200)))^2 / typeof(ps)(2.88))

@inline function _layer_ozone(ozone::OzoneColumns, ij, p, ps)
    lnp = ozone.log_pressure
    n = length(lnp)
    x = log(p)
    x >= lnp[1] && return ozone.values[1, ij]
    x <= lnp[n] && return ozone.values[n, ij]
    k = 1
    while lnp[k + 1] > x
        k += 1
    end
    w = (lnp[k] - x) / (lnp[k] - lnp[k + 1])
    return ozone.values[k, ij] + w * (ozone.values[k + 1, ij] - ozone.values[k, ij])
end

"""The ozone identity recorded in checkpoints and diagnostics provenance."""
_prescribed_ozone_identity(::Nothing) = nothing
_prescribed_ozone_identity(ozone::PrescribedOzone) = (sha256 = ozone.data.sha256, fixed_year = ozone.fixed_year)

"""The prescribed ozone for `config`, or `nothing` for the analytic production profile."""
_prescribed_ozone(config, spectral_grid) = isempty(config.forcing.ozone_path) ? nothing :
    PrescribedOzone(config.forcing.ozone_path, spectral_grid; fixed_year = config.forcing.ozone_fixed_year)

"""
Reject a checkpoint whose ozone differs from the reconstructed radiation, before any state changes. A checkpoint
written before this option existed has no record and restores only into a configuration with the analytic profile.
"""
function _validate_prescribed_ozone_restart(radiation, state)
    saved = hasproperty(state, :prescribed_ozone) ? state.prescribed_ozone : nothing
    current = _prescribed_ozone_identity(radiation.ozone)
    saved == current || throw(ArgumentError(
        "checkpoint and reconstructed radiation use different ozone " *
        "(checkpoint $(something(saved, "analytic profile")), configuration $(something(current, "analytic profile")))",
    ))
    return nothing
end

"""Column ozone (DU) of each model column from RRTMGP's layer mixing ratios and level pressures (host arrays)."""
function _model_column_ozone_du(vmr_o3, p_lev, gravity, molar_mass_dry_air)
    vmr, p = Float64.(Array(vmr_o3)), Float64.(Array(p_lev))
    nlay = size(vmr, 1)
    dp = p[1:nlay, :] .- p[2:(nlay + 1), :]  # RRTMGP levels run bottom to top
    return vec(sum(vmr .* dp; dims = 1)) ./ gravity ./ molar_mass_dry_air .* 6.02214076e23 ./ 2.6867e20
end

"""Diagnostics provenance of the prescribed ozone, written only when requested."""
function _prescribed_ozone_provenance(radiation, config, rrtmgp_active)
    isempty(config.forcing.ozone_path) && return (;)
    ozone = rrtmgp_active ? radiation.ozone : nothing
    isnothing(ozone) && return (prescribed_ozone = "requested_but_not_active",)
    solver = radiation.solver
    column = _model_column_ozone_du(RRTMGP.volume_mixing_ratio(solver, "o3"), RRTMGP.level_pressure(solver),
                                    9.80665, 0.0289644)
    weights = ozone.point_weights
    return (
        prescribed_ozone = ozone.fixed_year > 0 ?
            "zonal_mean_ozone_fixed_at_year_$(ozone.fixed_year)" : PRESCRIBED_OZONE_TIME_DEPENDENCE,
        ozone_path = ozone.data.path,
        ozone_sha256 = ozone.data.sha256,
        ozone_family = ozone.data.family,
        ozone_scenario = ozone.data.scenario,
        ozone_fixed_year = ozone.fixed_year,
        ozone_final_model_column_du_global_mean = sum(weights .* column) / sum(weights),
    )
end
