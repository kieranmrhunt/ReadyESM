"""Provisional zonal temperatures from one ERA5 analysis for bounded January 1993 experiments."""
struct UpperRadiativeCap{A}
    path::String
    sha256::String
    temperature_5hpa::A
    temperature_1hpa::A
end

"""Experimental 1993 sampled annual cycle, repeated with continuous month-centre interpolation."""
struct MonthlyUpperRadiativeCap{A}
    path::String
    sha256::String
    # field × month × RRTMGP column: two layer means, then 5 hPa, 1 hPa, top faces.
    temperatures::A
end

const UPPER_CAP_SPLIT_5HPA_PA = 500.0
const UPPER_CAP_SPLIT_1HPA_PA = 100.0
# The pinned RRTMGP longwave gas-optics table (rrtmgp-data 1.9) has a
# minimum reference pressure of 1.005183574463 Pa. RRTMGP clips every solver
# pressure below that floor before each flux solve, so use the effective top
# pressure in geometry and provenance from the outset.
const UPPER_CAP_TOP_PA = 1.005183574463
const UPPER_CAP_H2O_VMR = 4.0e-6
const UPPER_CAP_JAN1993_SHA256 = "748910cd885d0d52f57f368373a444bb8aba54bbf05cd17505986dd2cbd2f9c3"
const UPPER_CAP_MONTHLY1993_SHA256 = "5f07536991bd4a6a59536f7e5d86b58af803813c043d154d7a7d0d9184f5ca28"
const UPPER_CAP_MONTHLY_HEADER = "month,latitude_degrees,temperature_5_1hpa_k,temperature_1hpa_1pa_k,temperature_5hpa_boundary_k,temperature_1hpa_boundary_k,temperature_top_1pa_k"
const UPPER_CAP_MONTHLY_METHOD = "era5_1993_eight_samples_per_month_cyclic_month_centres_v1"

function _validate_upper_radiative_cap_options(config, rrtmgp_active, check_input_files)
    path = config.forcing.upper_radiative_cap_path
    isempty(path) && return nothing
    rrtmgp_active || throw(ArgumentError("upper_radiative_cap_path requires RRTMGP radiation"))
    if !isfile(path)
        check_input_files && throw(ArgumentError("upper radiative cap CSV not found: $path"))
        return nothing
    end
    header = open(readline, path)
    if header == UPPER_CAP_MONTHLY_HEADER
        check_input_files || return nothing
        _monthly_upper_cap_data(path)
        bytes2hex(open(SHA.sha256, path)) == UPPER_CAP_MONTHLY1993_SHA256 ||
            throw(ArgumentError("the sampled 1993 monthly upper cap differs from the audited table"))
        return nothing
    end
    header == "latitude_degrees,temperature_5hpa_k,temperature_1hpa_k" ||
        throw(ArgumentError("unrecognized upper radiative cap CSV header"))
    stop = config.start_date + Dates.Millisecond(round(Int64, config.duration_days * 86_400_000))
    DateTime(1993, 1, 1) <= config.start_date && stop <= DateTime(1993, 2, 1) ||
        throw(ArgumentError("the provisional 1993-01-01 upper cap is valid only within January 1993"))
    check_input_files || return nothing
    _upper_cap_zonal_data(path)
    bytes2hex(open(SHA.sha256, path)) == UPPER_CAP_JAN1993_SHA256 ||
        throw(ArgumentError("the provisional upper cap CSV differs from the reviewed ERA5 preparation"))
    return nothing
end

function _monthly_upper_cap_data(path::AbstractString)
    isfile(path) || throw(ArgumentError("monthly upper radiative cap CSV not found: $path"))
    lines = readlines(path)
    length(lines) == 1 + 12 * 181 ||
        throw(ArgumentError("monthly upper cap needs 12 months and 181 latitude rows per month"))
    lines[1] == UPPER_CAP_MONTHLY_HEADER ||
        throw(ArgumentError("unrecognized monthly upper radiative cap CSV header"))
    values = Array{Float64}(undef, 5, 12, 181)
    for month in 1:12, latitude in -90:90
        row = 1 + (month - 1) * 181 + latitude + 91
        cells = split(lines[row], ',')
        length(cells) == 7 || throw(ArgumentError("invalid monthly upper cap row $row"))
        parse(Int, cells[1]) == month && parse(Int, cells[2]) == latitude ||
            throw(ArgumentError("monthly upper cap month/latitude order differs at row $row"))
        values[:, month, latitude + 91] .= parse.(Float64, cells[3:7])
    end
    all(isfinite, values) && all(150 .< values .< 350) ||
        throw(ArgumentError("monthly upper cap temperatures are nonfinite or implausible"))
    return values
end

function _upper_cap_zonal_data(path::AbstractString)
    isfile(path) || throw(ArgumentError("upper radiative cap CSV not found: $path"))
    lines = readlines(path)
    length(lines) == 182 || throw(ArgumentError("upper radiative cap CSV needs 181 latitude rows"))
    lines[1] == "latitude_degrees,temperature_5hpa_k,temperature_1hpa_k" ||
        throw(ArgumentError("unrecognized upper radiative cap CSV header"))
    lat = Vector{Float64}(undef, 181)
    t5 = similar(lat)
    t1 = similar(lat)
    for (i, line) in enumerate(lines[2:end])
        cells = split(line, ',')
        length(cells) == 3 || throw(ArgumentError("invalid upper radiative cap row $i"))
        lat[i], t5[i], t1[i] = parse.(Float64, cells)
    end
    lat == collect(-90.0:90.0) || throw(ArgumentError("upper radiative cap needs the exact 1-degree latitude grid"))
    all(isfinite, t5) && all(isfinite, t1) && all(180 .< t5 .< 330) && all(180 .< t1 .< 330) ||
        throw(ArgumentError("upper radiative cap temperatures are not finite or plausible"))
    return (; lat, t5, t1)
end

@inline function _upper_cap_interpolate_lat(values, latitude)
    x = clamp(Float64(latitude), -90.0, 90.0)
    lo = clamp(floor(Int, x) + 91, 1, 180)
    w = x - (lo - 91)
    return (1 - w) * values[lo] + w * values[lo + 1]
end

function UpperRadiativeCap(path::AbstractString, solver)
    data = _upper_cap_zonal_data(path)
    digest = bytes2hex(open(SHA.sha256, path))
    digest == UPPER_CAP_JAN1993_SHA256 ||
        throw(ArgumentError("the provisional upper cap CSV differs from the reviewed ERA5 preparation"))
    latitude = Array(RRTMGP.latitude(solver))
    NF = eltype(RRTMGP.layer_temperature(solver))
    t5_host = NF[_upper_cap_interpolate_lat(data.t5, lat) for lat in latitude]
    t1_host = NF[_upper_cap_interpolate_lat(data.t1, lat) for lat in latitude]
    t5 = similar(RRTMGP.latitude(solver))
    t1 = similar(t5)
    copyto!(t5, t5_host)
    copyto!(t1, t1_host)
    return UpperRadiativeCap(String(path), digest, t5, t1)
end

function MonthlyUpperRadiativeCap(path::AbstractString, solver)
    data = _monthly_upper_cap_data(path)
    digest = bytes2hex(open(SHA.sha256, path))
    digest == UPPER_CAP_MONTHLY1993_SHA256 ||
        throw(ArgumentError("the sampled 1993 monthly upper cap differs from the audited table"))
    latitude = Array(RRTMGP.latitude(solver))
    NF = eltype(RRTMGP.layer_temperature(solver))
    host = Array{NF}(undef, 5, 12, length(latitude))
    for column in eachindex(latitude), month in 1:12, field in 1:5
        host[field, month, column] = NF(_upper_cap_interpolate_lat(
            view(data, field, month, :), latitude[column],
        ))
    end
    temperatures = similar(RRTMGP.latitude(solver), NF, size(host))
    copyto!(temperatures, host)
    return MonthlyUpperRadiativeCap(String(path), digest, temperatures)
end

function _construct_upper_radiative_cap(path::AbstractString, solver)
    header = open(readline, path)
    return header == UPPER_CAP_MONTHLY_HEADER ?
        MonthlyUpperRadiativeCap(path, solver) : UpperRadiativeCap(path, solver)
end

"""Middle of a Gregorian month, rounded down to the nearest millisecond."""
function _upper_cap_month_centre(year::Int, month::Int)
    start = DateTime(year, month, 1)
    stop = month == 12 ? DateTime(year + 1, 1, 1) : DateTime(year, month + 1, 1)
    return start + Dates.Millisecond(Dates.value(stop - start) ÷ 2)
end

"""Cyclic annual month-centre bracket and interpolation fraction for a model time."""
function _upper_cap_month_bracket(time::DateTime)
    year, month = Dates.year(time), Dates.month(time)
    centre = _upper_cap_month_centre(year, month)
    if time < centre
        left_month = month == 1 ? 12 : month - 1
        left_year = month == 1 ? year - 1 : year
        left = _upper_cap_month_centre(left_year, left_month)
        right = centre
        right_month = month
    else
        left = centre
        left_month = month
        right_month = month == 12 ? 1 : month + 1
        right_year = month == 12 ? year + 1 : year
        right = _upper_cap_month_centre(right_year, right_month)
    end
    fraction = Dates.value(time - left) / Dates.value(right - left)
    0 <= fraction <= 1 || error("upper cap month interpolation left its bracket")
    return left_month, right_month, fraction
end

_upper_cap_identity(::Nothing) = nothing
_upper_cap_identity(cap::UpperRadiativeCap) = (sha256 = cap.sha256, method = "era5_19930101_zonal_start_state_proxy_v1")
_upper_cap_identity(cap::MonthlyUpperRadiativeCap) =
    (sha256 = cap.sha256, method = UPPER_CAP_MONTHLY_METHOD)

function _validate_upper_cap_restart(radiation, state)
    saved = hasproperty(state, :upper_radiative_cap) ? state.upper_radiative_cap : nothing
    current = _upper_cap_identity(radiation.upper_cap)
    saved == current || throw(ArgumentError(
        "checkpoint and reconstructed radiation use different upper radiative cap profiles"))
    return nothing
end

function _upper_cap_provenance(radiation, config, rrtmgp_active)
    isempty(config.forcing.upper_radiative_cap_path) && return (;)
    rrtmgp_active && !isnothing(radiation.upper_cap) ||
        throw(ArgumentError("configured upper radiative cap was not applied"))
    if radiation.upper_cap isa MonthlyUpperRadiativeCap
        return (;
            upper_radiative_cap_method = UPPER_CAP_MONTHLY_METHOD,
            upper_radiative_cap_path = radiation.upper_cap.path,
            upper_radiative_cap_sha256 = radiation.upper_cap.sha256,
            upper_radiative_cap_time_dependence = "1993_sampled_annual_cycle_repeated_all_years_linear_between_gregorian_month_centres",
            upper_radiative_cap_lower_split_pressure_pa = UPPER_CAP_SPLIT_5HPA_PA,
            upper_radiative_cap_upper_split_pressure_pa = UPPER_CAP_SPLIT_1HPA_PA,
            upper_radiative_cap_top_pressure_pa = UPPER_CAP_TOP_PA,
            upper_radiative_cap_temperature_above_1hpa = "era5_pressure_weighted_layer_mean_and_near_1pa_top",
            upper_radiative_cap_heating = "three_radiative_flux_divergences_aggregated_to_dynamic_top_mass_layer",
            upper_radiative_cap_aerosol_scheme = "zero_band_optics_above_5hpa_v1",
        )
    end
    return (;
        upper_radiative_cap_method = "era5_19930101_zonal_start_state_proxy_v1",
        upper_radiative_cap_path = radiation.upper_cap.path,
        upper_radiative_cap_sha256 = radiation.upper_cap.sha256,
        upper_radiative_cap_lower_split_pressure_pa = UPPER_CAP_SPLIT_5HPA_PA,
        upper_radiative_cap_upper_split_pressure_pa = UPPER_CAP_SPLIT_1HPA_PA,
        upper_radiative_cap_top_pressure_pa = UPPER_CAP_TOP_PA,
        upper_radiative_cap_temperature_above_1hpa = "isothermal_at_1hpa_provisional",
        upper_radiative_cap_heating = "three_radiative_flux_divergences_aggregated_to_dynamic_top_mass_layer",
        upper_radiative_cap_aerosol_scheme = "zero_band_optics_above_5hpa_v1",
    )
end

_fill_upper_radiative_cap!(::Nothing, solver, surface_pressure, ozone_columns, nmodel, time) = nothing

function _fill_upper_radiative_cap!(cap::UpperRadiativeCap, solver, surface_pressure, ozone_columns, nmodel, time)
    p_lay = RRTMGP.layer_pressure(solver)
    p_lev = RRTMGP.level_pressure(solver)
    t_lay = RRTMGP.layer_temperature(solver)
    t_lev = RRTMGP.level_temperature(solver)
    size(p_lay, 1) == nmodel + 2 || throw(DimensionMismatch("upper cap needs two extra solver layers"))
    size(p_lev, 1) == nmodel + 3 || throw(DimensionMismatch("upper cap needs three extra solver levels"))
    backend = KernelAbstractions.get_backend(p_lay)
    _fill_upper_radiative_cap_kernel!(backend)(
        p_lay, p_lev, t_lay, t_lev,
        RRTMGP.volume_mixing_ratio(solver, "h2o"),
        RRTMGP.volume_mixing_ratio(solver, "o3"),
        RRTMGP.layer_relative_humidity(solver),
        cap.temperature_5hpa, cap.temperature_1hpa,
        surface_pressure, ozone_columns, nmodel;
        ndrange = size(p_lay, 2),
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

KernelAbstractions.@kernel function _fill_upper_radiative_cap_kernel!(
    p_lay, p_lev, t_lay, t_lev, h2o, o3, rh, t5_columns, t1_columns,
    surface_pressure, ozone_columns, nmodel,
)
    ij = @index(Global, Linear)
    NF = eltype(p_lay)
    p5 = NF(UPPER_CAP_SPLIT_5HPA_PA)
    p1 = NF(UPPER_CAP_SPLIT_1HPA_PA)
    ptop = NF(UPPER_CAP_TOP_PA)
    t5 = t5_columns[ij]
    t1 = t1_columns[ij]
    p_middle = sqrt(p5 * p1)
    p_upper = sqrt(p1 * ptop)
    w = (log(p_middle) - log(p5)) / (log(p1) - log(p5))
    p_lev[nmodel + 1, ij] = p5
    p_lev[nmodel + 2, ij] = p1
    p_lev[nmodel + 3, ij] = ptop
    t_lev[nmodel + 1, ij] = t5
    t_lev[nmodel + 2, ij] = t1
    t_lev[nmodel + 3, ij] = t1
    p_lay[nmodel + 1, ij] = p_middle
    p_lay[nmodel + 2, ij] = p_upper
    t_lay[nmodel + 1, ij] = t5 + w * (t1 - t5)
    t_lay[nmodel + 2, ij] = t1
    h2o[nmodel + 1, ij] = NF(UPPER_CAP_H2O_VMR)
    h2o[nmodel + 2, ij] = NF(UPPER_CAP_H2O_VMR)
    o3[nmodel + 1, ij] = _layer_ozone(ozone_columns, ij, p_middle, surface_pressure[ij])
    o3[nmodel + 2, ij] = _layer_ozone(ozone_columns, ij, p_upper, surface_pressure[ij])
    rh[nmodel + 1, ij] = zero(NF)
    rh[nmodel + 2, ij] = zero(NF)
end

function _fill_upper_radiative_cap!(
    cap::MonthlyUpperRadiativeCap, solver, surface_pressure, ozone_columns, nmodel, time::DateTime,
)
    p_lay = RRTMGP.layer_pressure(solver)
    p_lev = RRTMGP.level_pressure(solver)
    size(p_lay, 1) == nmodel + 2 || throw(DimensionMismatch("monthly upper cap needs two solver layers"))
    size(p_lev, 1) == nmodel + 3 || throw(DimensionMismatch("monthly upper cap needs three solver levels"))
    size(cap.temperatures) == (5, 12, size(p_lay, 2)) ||
        throw(DimensionMismatch("monthly upper cap table differs from solver columns"))
    left, right, fraction = _upper_cap_month_bracket(time)
    backend = KernelAbstractions.get_backend(p_lay)
    _fill_monthly_upper_radiative_cap_kernel!(backend)(
        p_lay, p_lev, RRTMGP.layer_temperature(solver), RRTMGP.level_temperature(solver),
        RRTMGP.volume_mixing_ratio(solver, "h2o"),
        RRTMGP.volume_mixing_ratio(solver, "o3"),
        RRTMGP.layer_relative_humidity(solver), cap.temperatures,
        left, right, fraction, surface_pressure, ozone_columns, nmodel;
        ndrange = size(p_lay, 2),
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

KernelAbstractions.@kernel function _fill_monthly_upper_radiative_cap_kernel!(
    p_lay, p_lev, t_lay, t_lev, h2o, o3, rh, temperatures,
    left, right, fraction, surface_pressure, ozone_columns, nmodel,
)
    ij = @index(Global, Linear)
    NF = eltype(p_lay)
    p5 = NF(UPPER_CAP_SPLIT_5HPA_PA)
    p1 = NF(UPPER_CAP_SPLIT_1HPA_PA)
    ptop = NF(UPPER_CAP_TOP_PA)
    p_middle = sqrt(p5 * p1)
    p_upper = sqrt(p1 * ptop)
    weight = NF(fraction)
    t_lower = (one(NF) - weight) * temperatures[1, left, ij] + weight * temperatures[1, right, ij]
    t_upper = (one(NF) - weight) * temperatures[2, left, ij] + weight * temperatures[2, right, ij]
    t5 = (one(NF) - weight) * temperatures[3, left, ij] + weight * temperatures[3, right, ij]
    t1 = (one(NF) - weight) * temperatures[4, left, ij] + weight * temperatures[4, right, ij]
    ttop = (one(NF) - weight) * temperatures[5, left, ij] + weight * temperatures[5, right, ij]
    p_lev[nmodel + 1, ij] = p5
    p_lev[nmodel + 2, ij] = p1
    p_lev[nmodel + 3, ij] = ptop
    t_lev[nmodel + 1, ij] = t5
    t_lev[nmodel + 2, ij] = t1
    t_lev[nmodel + 3, ij] = ttop
    p_lay[nmodel + 1, ij] = p_middle
    p_lay[nmodel + 2, ij] = p_upper
    t_lay[nmodel + 1, ij] = t_lower
    t_lay[nmodel + 2, ij] = t_upper
    h2o[nmodel + 1, ij] = NF(UPPER_CAP_H2O_VMR)
    h2o[nmodel + 2, ij] = NF(UPPER_CAP_H2O_VMR)
    o3[nmodel + 1, ij] = _layer_ozone(ozone_columns, ij, p_middle, surface_pressure[ij])
    o3[nmodel + 2, ij] = _layer_ozone(ozone_columns, ij, p_upper, surface_pressure[ij])
    rh[nmodel + 1, ij] = zero(NF)
    rh[nmodel + 2, ij] = zero(NF)
end

function _assign_radiative_heating!(heating, solver, ::Nothing, nmodel)
    heating .= RRTMGP.heating_rate(solver)
    return nothing
end

function _assign_radiative_heating!(heating, solver, ::Union{UpperRadiativeCap, MonthlyUpperRadiativeCap}, nmodel)
    size(heating, 1) == nmodel || throw(DimensionMismatch("model-facing heating has wrong layer count"))
    solver_heating = RRTMGP.heating_rate(solver)
    size(solver_heating, 1) == nmodel + 2 ||
        throw(DimensionMismatch("capped solver heating has wrong layer count"))
    @views heating[1:(nmodel - 1), :] .= solver_heating[1:(nmodel - 1), :]
    net = RRTMGP.net_flux(solver)
    p_lev = RRTMGP.level_pressure(solver)
    NF = eltype(heating)
    grav_over_cp = NF(RRTMGP.RP.grav(solver.params) / RRTMGP.RP.cp_d(solver.params))
    # Sum the flux divergence of the model top layer and both radiative cap
    # layers, then divide by the mass of the original top dynamic sigma layer.
    # The dynamic upper face is zero, even though the solver caps it at 1 Pa.
    @views heating[nmodel, :] .= grav_over_cp .* (net[nmodel + 3, :] .- net[nmodel, :]) ./
        (zero(NF) .- p_lev[nmodel, :])
    return nothing
end

_upper_cap_layer_heights(::Nothing, geopotential, orography, gravity, solver, nmodel) =
    _layer_heights(geopotential, orography, gravity, nmodel)

function _upper_cap_layer_heights(
    cap::Union{UpperRadiativeCap, MonthlyUpperRadiativeCap},
    geopotential, orography, gravity, solver, nmodel,
)
    zc_model, zb_model = _layer_heights(geopotential, orography, gravity, nmodel)
    ncol = size(zc_model, 2)
    zc = zeros(nmodel + 2, ncol)
    zb = zeros(nmodel + 3, ncol)
    zc[1:nmodel, :] .= zc_model
    zb[1:nmodel, :] .= zb_model[1:nmodel, :]
    p_lay = Float64.(Array(RRTMGP.layer_pressure(solver)))
    t_lay = Float64.(Array(RRTMGP.layer_temperature(solver)))
    t_lev = Float64.(Array(RRTMGP.level_temperature(solver)))
    gas_constant = Float64(RRTMGP.RP.R_d(solver.params))
    g = Float64(gravity)
    for c in 1:ncol
        p_model = p_lay[nmodel, c]
        p_model > UPPER_CAP_SPLIT_5HPA_PA ||
            throw(ArgumentError("dynamic top-layer centre is at or above the 5 hPa cap split"))
        t_model = t_lay[nmodel, c]
        t5, t1 = t_lev[nmodel + 1, c], t_lev[nmodel + 2, c]
        z5 = zc_model[nmodel, c] + gas_constant / g *
            (t_model + t5) / 2 * log(p_model / UPPER_CAP_SPLIT_5HPA_PA)
        lower_mean = cap isa MonthlyUpperRadiativeCap ? t_lay[nmodel + 1, c] : (t5 + t1) / 2
        z1 = z5 + gas_constant / g * lower_mean *
            log(UPPER_CAP_SPLIT_5HPA_PA / UPPER_CAP_SPLIT_1HPA_PA)
        ztop = z1 + gas_constant / g * t_lay[nmodel + 2, c] *
            log(UPPER_CAP_SPLIT_1HPA_PA / UPPER_CAP_TOP_PA)
        zb[nmodel + 1, c] = z5
        zb[nmodel + 2, c] = z1
        zb[nmodel + 3, c] = ztop
        zc[nmodel + 1, c] = z5 + gas_constant / g *
            (t5 + t_lay[nmodel + 1, c]) / 2 *
            log(UPPER_CAP_SPLIT_5HPA_PA / p_lay[nmodel + 1, c])
        zc[nmodel + 2, c] = z1 + gas_constant / g *
            (t1 + t_lay[nmodel + 2, c]) / 2 *
            log(UPPER_CAP_SPLIT_1HPA_PA / p_lay[nmodel + 2, c])
    end
    all(zb[1:end-1, :] .< zb[2:end, :]) ||
        throw(ArgumentError("upper radiative cap produced nonmonotonic aerosol heights"))
    return zc, zb
end
