# Scored output files (opt-in, `scored_outputs: true`).
#
# The release scorer reads four files besides diagnostics.nc: monthly.nc, wind.nc, toa_monthly.nc and eddy_6hourly.nc.
# Until v0.1.17 these came from observer scripts kept outside the model. This file ports those observers into the model
# so that a plain run writes every scored file. The arithmetic, the sampling times, the file attributes and the variable
# order are those of the observers used for the v0.1.17 release scores:
#   monthly core and capture: artifacts/development_wave_20261019/c02_monthly_capture_v1/{monthly_core,adapter}.jl
#   wind core and capture:    artifacts/development_wave_20261019/c07_monthly_wind_v1/{wind_core,adapter}.jl
#   combined observer:        artifacts/development_wave_20261020/c06_coupled_skill_v1/observer.jl (143c6bc1...)
#   top-of-atmosphere:        artifacts/development_wave_20261022/score_v4_inputs/toa_capture.jl (f2cede6a...)
#   six-hourly eddy fields:   artifacts/development_wave_20261024/full_temperature_screen_v1/eddy_capture.jl (3fbc6f03...)
# Two things differ from the observers, neither of which changes a written byte for a 900 s step:
# - The averaging window starts at config.start_date and lasts config.duration_days. The observers fixed the start at
#   1 January 1993, and the July observer patched that start in its source text.
# - The step counts come from the model timestep. The observers assumed 96 steps per day.

module ScoredMonthlyCore
using Dates

mutable struct Accumulator
    start::DateTime
    stop::DateTime
    step_ms::Int
    sigma::Vector{Float64}
    dsigma::Vector{Float64}
    levels_Pa::Vector{Float64}
    gravity::Float64
    water_density::Float64
    state_count::Int
    flux_count::Int
    water_sum::Vector{Float64}
    ps_sum::Vector{Float64}
    temperature_sum::Matrix{Float64}
    temperature_count::Matrix{Int}
    precipitation_kg_m2::Vector{Float64}
end

function Accumulator(start, stop, npoints, sigma, dsigma, gravity, water_density;
                     step_ms = 900000, levels_Pa = [20000.0, 50000.0, 85000.0])
    span = Dates.value(stop - start)
    span > 0 && span % 3600000 == 0 || error("Window must contain whole hours")
    0 < step_ms <= 3600000 && 3600000 % step_ms == 0 || error("Step must divide one hour")
    npoints > 0 || error("Empty grid")
    length(sigma) == length(dsigma) >= 2 || error("Layer count mismatch")
    all(isfinite, sigma) && all(>(0), sigma) && all(<=(1), sigma) &&
        all(>(0), diff(sigma)) || error("Sigma must increase toward the surface")
    all(isfinite, dsigma) && all(>(0), dsigma) || error("Invalid layer thickness")
    isfinite(gravity) && gravity > 0 && isfinite(water_density) && water_density > 0 || error("Invalid constants")
    all(isfinite, levels_Pa) && all(>(0), levels_Pa) && all(>(0), diff(levels_Pa)) || error("Invalid target levels")
    return Accumulator(start, stop, step_ms, Float64.(sigma), Float64.(dsigma), Float64.(levels_Pa),
        Float64(gravity), Float64(water_density), 0, 0, zeros(npoints), zeros(npoints),
        zeros(npoints, length(levels_Pa)), zeros(Int, npoints, length(levels_Pa)), zeros(npoints))
end

# No extrapolation below the lowest full level or above the highest full level.
function pressure_temperature(temperature, sigma, ps, target)
    target < sigma[1] * ps && return nothing
    target > sigma[end] * ps && return nothing
    x = target / ps
    upper = searchsortedfirst(sigma, x)
    upper > length(sigma) && return nothing
    upper == 1 && return Float64(temperature[1])
    sigma[upper] == x && return Float64(temperature[upper])
    lower = upper - 1
    fraction = log(x / sigma[lower]) / log(sigma[upper] / sigma[lower])
    return (1 - fraction) * Float64(temperature[lower]) + fraction * Float64(temperature[upper])
end

function state!(a::Accumulator, time, temperature, humidity, ps)
    time == a.start + Hour(a.state_count) && time < a.stop || error("Missing, duplicate or out-of-window hourly state")
    n = length(a.ps_sum)
    size(temperature) == size(humidity) == (n, length(a.sigma)) && length(ps) == n || error("State shape mismatch")
    all(isfinite, temperature) && all(isfinite, humidity) && all(isfinite, ps) && all(>(0), ps) || error("Invalid state")
    for i in 1:n
        a.ps_sum[i] += Float64(ps[i])
        a.water_sum[i] += sum(Float64(humidity[i, k]) * a.dsigma[k] for k in eachindex(a.sigma)) * Float64(ps[i]) / a.gravity
        for j in eachindex(a.levels_Pa)
            value = pressure_temperature(view(temperature, i, :), a.sigma, Float64(ps[i]), a.levels_Pa[j])
            isnothing(value) && continue
            a.temperature_sum[i, j] += value
            a.temperature_count[i, j] += 1
        end
    end
    a.state_count += 1
    return nothing
end

function flux!(a::Accumulator, time, rain, snow)
    time == a.start + Millisecond((a.flux_count + 1) * a.step_ms) && time <= a.stop || error("Missing, duplicate or out-of-window flux")
    length(rain) == length(snow) == length(a.ps_sum) || error("Flux shape mismatch")
    all(isfinite, rain) && all(isfinite, snow) || error("Nonfinite precipitation")
    dt = a.step_ms / 1000
    for i in eachindex(rain)
        a.precipitation_kg_m2[i] += (Float64(rain[i]) + Float64(snow[i])) * a.water_density * dt
    end
    a.flux_count += 1
    return nothing
end

function finish(a::Accumulator)
    duration_ms = Dates.value(a.stop - a.start)
    a.state_count == duration_ms ÷ 3600000 || error("Incomplete hourly states")
    a.flux_count == duration_ms ÷ a.step_ms || error("Incomplete precipitation intervals")
    valid = a.temperature_count .== a.state_count
    temperature = fill(NaN, size(a.temperature_sum))
    temperature[valid] .= a.temperature_sum[valid] ./ a.state_count
    seconds = duration_ms / 1000
    return (; temperature_K = temperature, temperature_valid = valid,
        temperature_sample_count = copy(a.temperature_count),
        temperature_sum_K = copy(a.temperature_sum),
        column_water_kg_m2 = a.water_sum ./ a.state_count,
        surface_pressure_Pa = a.ps_sum ./ a.state_count,
        precipitation_kg_m2 = copy(a.precipitation_kg_m2),
        precipitation_kg_m2_s = a.precipitation_kg_m2 ./ seconds,
        precipitation_mm_day = a.precipitation_kg_m2 .* (86400 / seconds),
        state_samples = a.state_count, flux_samples = a.flux_count, duration_seconds = seconds)
end
end

module ScoredWindCore
using Dates
import ..ScoredMonthlyCore

mutable struct Accumulator
    start::DateTime
    stop::DateTime
    sigma::Vector{Float64}
    levels::Vector{Float64}
    samples::Int
    u_sum::Matrix{Float64}
    v_sum::Matrix{Float64}
    count::Matrix{Int}
end

function Accumulator(start, stop, n, sigma; levels = [20000.0, 50000.0, 85000.0])
    duration = Dates.value(stop - start)
    duration > 0 && duration % 3600000 == 0 || error("Invalid window")
    n > 0 && length(sigma) >= 2 && all(isfinite, sigma) && all(>(0), sigma) &&
        all(<=(1), sigma) && all(>(0), diff(sigma)) || error("Invalid sigma grid")
    all(isfinite, levels) && all(>(0), levels) && all(>(0), diff(levels)) || error("Invalid levels")
    return Accumulator(start, stop, Float64.(sigma), Float64.(levels), 0,
        zeros(n, length(levels)), zeros(n, length(levels)), zeros(Int, n, length(levels)))
end

function state!(a, time, u, v, ps)
    time == a.start + Hour(a.samples) && time < a.stop || error("Incomplete wind calendar")
    n = size(a.count, 1)
    size(u) == size(v) == (n, length(a.sigma)) && length(ps) == n || error("Wind shape mismatch")
    all(isfinite, u) && all(isfinite, v) && all(isfinite, ps) && all(>(0), ps) || error("Invalid wind state")
    for i in 1:n, j in eachindex(a.levels)
        x = ScoredMonthlyCore.pressure_temperature(view(u, i, :), a.sigma, Float64(ps[i]), a.levels[j])
        y = ScoredMonthlyCore.pressure_temperature(view(v, i, :), a.sigma, Float64(ps[i]), a.levels[j])
        isnothing(x) == isnothing(y) || error("Wind masks differ")
        isnothing(x) && continue
        a.u_sum[i, j] += x
        a.v_sum[i, j] += y
        a.count[i, j] += 1
    end
    a.samples += 1
    return nothing
end

function finish(a)
    a.samples == Dates.value(a.stop - a.start) ÷ 3600000 || error("Incomplete wind samples")
    valid = a.count .== a.samples
    u = fill(NaN, size(a.count))
    v = copy(u)
    u[valid] .= a.u_sum[valid] ./ a.samples
    v[valid] .= a.v_sum[valid] ./ a.samples
    return (; u, v, valid, count = copy(a.count), u_sum = copy(a.u_sum), v_sum = copy(a.v_sum), samples = a.samples)
end
end

_scored_host(field) = Array(Base.has_offset_axes(field.data) ? parent(field.data) : field.data)

mutable struct ScoredMonthlyCapture <: SpeedyWeather.AbstractCallback
    accumulator::ScoredMonthlyCore.Accumulator
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    weights::Vector{Float64}
end

function _scored_monthly_state!(c, vars)
    ScoredMonthlyCore.state!(c.accumulator, vars.prognostic.clock.time,
        _scored_host(vars.grid.temperature), _scored_host(vars.grid.humidity),
        exp.(Float64.(vec(_scored_host(vars.grid.pressure)))))
end

function SpeedyWeather.initialize!(c::ScoredMonthlyCapture, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel)
    c.accumulator.state_count == c.accumulator.flux_count == 0 || error("Monthly callback reused")
    _scored_monthly_state!(c, vars)
end

function SpeedyWeather.callback!(c::ScoredMonthlyCapture, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel)
    time = vars.prognostic.clock.time
    a = c.accumulator
    ScoredMonthlyCore.flux!(a, time, vec(_scored_host(vars.parameterizations.rain_rate)),
        vec(_scored_host(vars.parameterizations.snow_rate)))
    if time < a.stop && Dates.value(time - a.start) % 3600000 == 0
        _scored_monthly_state!(c, vars)
    end
    return nothing
end
SpeedyWeather.finalize!(c::ScoredMonthlyCapture, args...) = (ScoredMonthlyCore.finish(c.accumulator); nothing)

function _scored_save_monthly(c, path)
    ispath(path) && error("Monthly output exists")
    mkpath(dirname(path))
    a = c.accumulator
    data = ScoredMonthlyCore.finish(a)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "point", length(c.weights))
        NCDatasets.defDim(ds, "pressure", length(a.levels_Pa))
        NCDatasets.defDim(ds, "layer", length(a.sigma))
        ds.attrib["schema"] = "c02_monthly_capture_v1"
        ds.attrib["start_time_UTC"] = string(a.start)
        ds.attrib["stop_time_UTC"] = string(a.stop)
        ds.attrib["state_sampling"] = "Hourly instants, start inclusive and stop exclusive"
        ds.attrib["flux_sampling"] = "Completed-step rates times timestep, start exclusive and stop inclusive"
        ds.attrib["pressure_interpolation"] = "Linear in log pressure, no extrapolation. Complete hourly coverage required."
        ds.attrib["state_samples"] = data.state_samples
        ds.attrib["flux_samples"] = data.flux_samples
        ds.attrib["duration_seconds"] = data.duration_seconds
        ds.attrib["timestep_seconds"] = a.step_ms / 1000
        ds.attrib["gravity_m_s2"] = a.gravity
        ds.attrib["water_density_kg_m3"] = a.water_density
        for (name, values, dims, units) in (
            ("longitude", c.longitude, ("point",), "degrees_east"),
            ("latitude", c.latitude, ("point",), "degrees_north"),
            ("area_weight", c.weights, ("point",), "1"),
            ("pressure", a.levels_Pa, ("pressure",), "Pa"),
            ("sigma_full", a.sigma, ("layer",), "1"),
            ("sigma_thickness", a.dsigma, ("layer",), "1"),
            ("temperature", data.temperature_K, ("point", "pressure"), "K"),
            ("temperature_sum", data.temperature_sum_K, ("point", "pressure"), "K"),
            ("temperature_sample_count", Int32.(data.temperature_sample_count), ("point", "pressure"), "1"),
            ("temperature_valid", Int8.(data.temperature_valid), ("point", "pressure"), "1"),
            ("column_water", data.column_water_kg_m2, ("point",), "kg m-2"),
            ("surface_pressure", data.surface_pressure_Pa, ("point",), "Pa"),
            ("precipitation_amount", data.precipitation_kg_m2, ("point",), "kg m-2"),
            ("precipitation_rate", data.precipitation_kg_m2_s, ("point",), "kg m-2 s-1"),
            ("precipitation_mm_day", data.precipitation_mm_day, ("point",), "mm day-1"))
            v = NCDatasets.defVar(ds, name, eltype(values), dims; attrib = Dict("units" => units))
            v[ntuple(_ -> Colon(), ndims(values))...] = values
        end
    end
    return Dict("state_samples" => data.state_samples, "flux_samples" => data.flux_samples,
        "duration_seconds" => data.duration_seconds, "file" => basename(path),
        "start_time_UTC" => string(a.start), "stop_time_UTC" => string(a.stop))
end

mutable struct ScoredWindCapture <: SpeedyWeather.AbstractCallback
    accumulator::ScoredWindCore.Accumulator
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    weights::Vector{Float64}
end

function _scored_wind_state!(c, vars)
    ScoredWindCore.state!(c.accumulator, vars.prognostic.clock.time, _scored_host(vars.grid.u), _scored_host(vars.grid.v),
        exp.(Float64.(vec(_scored_host(vars.grid.pressure)))))
end
SpeedyWeather.initialize!(c::ScoredWindCapture, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel) =
    _scored_wind_state!(c, vars)
function SpeedyWeather.callback!(c::ScoredWindCapture, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel)
    time = vars.prognostic.clock.time
    a = c.accumulator
    if time < a.stop && Dates.value(time - a.start) % 3600000 == 0
        _scored_wind_state!(c, vars)
    end
    return nothing
end
SpeedyWeather.finalize!(c::ScoredWindCapture, args...) = (ScoredWindCore.finish(c.accumulator); nothing)

function _scored_save_wind(c, path)
    ispath(path) && error("Wind output exists")
    a = c.accumulator
    data = ScoredWindCore.finish(a)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "point", length(c.weights))
        NCDatasets.defDim(ds, "pressure", length(a.levels))
        ds.attrib["schema"] = "monthly_wind_capture_v1"
        ds.attrib["start_time_UTC"] = string(a.start)
        ds.attrib["stop_time_UTC"] = string(a.stop)
        ds.attrib["state_samples"] = data.samples
        ds.attrib["state_sampling"] = "Hourly instants, start inclusive and stop exclusive"
        ds.attrib["pressure_interpolation"] = "Linear in log pressure, no extrapolation. Complete hourly coverage required."
        for (name, values, dims, units) in (
            ("longitude", c.longitude, ("point",), "degrees_east"),
            ("latitude", c.latitude, ("point",), "degrees_north"),
            ("area_weight", c.weights, ("point",), "1"),
            ("pressure", a.levels, ("pressure",), "Pa"),
            ("u", data.u, ("point", "pressure"), "m s-1"),
            ("v", data.v, ("point", "pressure"), "m s-1"),
            ("u_sum", data.u_sum, ("point", "pressure"), "m s-1"),
            ("v_sum", data.v_sum, ("point", "pressure"), "m s-1"),
            ("valid", Int8.(data.valid), ("point", "pressure"), "1"),
            ("sample_count", Int32.(data.count), ("point", "pressure"), "1"))
            variable = NCDatasets.defVar(ds, name, eltype(values), dims; attrib = Dict("units" => units))
            variable[ntuple(_ -> Colon(), ndims(values))...] = values
        end
    end
    return Dict("state_samples" => data.samples, "file" => basename(path),
        "start_time_UTC" => string(a.start), "stop_time_UTC" => string(a.stop))
end

# The monthly and wind captures together, with checks that sampling leaves its inputs unchanged.
mutable struct ScoredSkillObserver <: SpeedyWeather.AbstractCallback
    monthly::ScoredMonthlyCapture
    wind::ScoredWindCapture
    initialized::Bool
    callbacks::Int
    passivity_checks::Int
    initial_hashes::Vector{String}
end

function _scored_input_hashes(vars)
    fields = (vars.grid.temperature, vars.grid.humidity, vars.grid.pressure,
        vars.grid.u, vars.grid.v, vars.parameterizations.rain_rate,
        vars.parameterizations.snow_rate)
    return [bytes2hex(sha256(reinterpret(UInt8, vec(_scored_host(x))))) for x in fields]
end

function _scored_initialize_skill!(c, vars, model)
    c.initialized && error("Repeated coupled skill initialization")
    before = _scored_input_hashes(vars)
    SpeedyWeather.initialize!(c.monthly, vars, model)
    SpeedyWeather.initialize!(c.wind, vars, model)
    before == _scored_input_hashes(vars) || error("Skill initialization changed its inputs")
    c.passivity_checks += 1
    c.initial_hashes = before
    c.initialized = true
    return nothing
end

function SpeedyWeather.callback!(c::ScoredSkillObserver, vars::SpeedyWeather.Variables, model::SpeedyWeather.AbstractModel)
    c.initialized || error("Missing initial coupled skill sample")
    before = _scored_input_hashes(vars)
    SpeedyWeather.callback!(c.monthly, vars, model)
    SpeedyWeather.callback!(c.wind, vars, model)
    before == _scored_input_hashes(vars) || error("Skill callback changed its inputs")
    c.callbacks += 1
    c.passivity_checks += 1
    return nothing
end
SpeedyWeather.finalize!(::ScoredSkillObserver, args...) = nothing

# Gridded top-of-atmosphere radiation: the RRTMGP solver's top-level fluxes at every grid point, the same arrays that
# the global radiation budget averages, sampled at hourly instants after completed steps.
const SCORED_TOA_FIELDS = ("incoming_shortwave", "outgoing_shortwave", "outgoing_longwave",
                           "clear_outgoing_shortwave", "clear_outgoing_longwave")

mutable struct ScoredTOACapture
    sums::Dict{String, Vector{Float64}}
    count::Int
    clear_sky::Bool
end

ScoredTOACapture(npts::Integer) = ScoredTOACapture(Dict(k => zeros(npts) for k in SCORED_TOA_FIELDS), 0, true)

_scored_toa_top(flux) = Float64.(Array(flux[end, :]))

function _scored_sample_toa!(toa::ScoredTOACapture, atmosphere_model)
    solver = atmosphere_model.longwave_radiation.solver
    RT = RRTMGP
    toa.sums["incoming_shortwave"] .+= _scored_toa_top(RT.sw_flux_dn(solver))
    toa.sums["outgoing_shortwave"] .+= _scored_toa_top(RT.sw_flux_up(solver))
    toa.sums["outgoing_longwave"] .+= _scored_toa_top(RT.lw_flux_up(solver))
    if isnothing(solver.clear_flux_lw)
        toa.clear_sky = false
    else
        toa.sums["clear_outgoing_shortwave"] .+= _scored_toa_top(RT.clear_sw_flux_up(solver))
        toa.sums["clear_outgoing_longwave"] .+= _scored_toa_top(RT.clear_lw_flux_up(solver))
    end
    toa.count += 1
    return nothing
end

function _scored_save_toa(toa::ScoredTOACapture, path, longitude, latitude, weights, start, stop)
    toa.count > 0 || error("TOA capture has no samples")
    w = Float64.(weights) ./ sum(weights)
    summary = Dict{String, Any}("samples" => toa.count, "clear_sky" => toa.clear_sky)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "point", length(w))
        ds.attrib["schema"] = "toa_monthly_v1"
        ds.attrib["sampling"] = "Hourly instants after completed steps, RRTMGP top-level fluxes"
        ds.attrib["samples"] = toa.count
        ds.attrib["window_start"] = string(start)
        ds.attrib["window_end"] = string(stop)
        ds.attrib["clear_sky"] = Int(toa.clear_sky)
        for (name, values) in (("longitude", Float64.(longitude)), ("latitude", Float64.(latitude)),
                               ("area_weight", Float64.(weights)))
            NCDatasets.defVar(ds, name, values, ("point",))
        end
        for k in SCORED_TOA_FIELDS
            startswith(k, "clear") && !toa.clear_sky && continue
            mean = toa.sums[k] ./ toa.count
            NCDatasets.defVar(ds, k, mean, ("point",); attrib = ["units" => "W m-2"])
            summary["global_mean_" * k] = sum(w .* mean)
        end
    end
    return summary
end

# Six-hourly pressure-level fields for eddy heat and momentum-flux diagnostics. MSLP uses the reduction in
# SpeedyWeather's MeanSeaLevelPressureOutput.
mutable struct ScoredSixHourlyCapture
    start::DateTime
    stop::DateTime
    sigma::Vector{Float64}
    levels_pa::Vector{Float64}
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    weights::Vector{Float64}
    samples::Vector{Any}
end

function ScoredSixHourlyCapture(start, stop, sigma, longitude, latitude, weights;
        levels_pa = [25000.0, 50000.0, 85000.0])
    Dates.value(stop - start) > 0 && Dates.value(stop - start) % (6 * 3600000) == 0 ||
        error("Capture window must contain whole six-hour periods")
    n = length(weights)
    length(longitude) == length(latitude) == n || error("Point-axis mismatch")
    length(sigma) >= 2 && all(>(0), diff(sigma)) || error("Invalid sigma levels")
    all(>(0), levels_pa) && all(>(0), diff(levels_pa)) || error("Invalid pressure levels")
    return ScoredSixHourlyCapture(start, stop, Float64.(sigma), Float64.(levels_pa),
        Float64.(longitude), Float64.(latitude), Float64.(weights), Any[])
end

function _scored_reduced_mslp_hpa(lnps, orography_m, surface_air_temperature_k, surface_q,
        gravity, gas_constant_dry, virtual_coefficient)
    n = length(lnps)
    length(orography_m) == length(surface_air_temperature_k) == length(surface_q) == n ||
        error("MSLP point-axis mismatch")
    temperature_virtual = Float64.(surface_air_temperature_k) .*
        (1 .+ Float64(virtual_coefficient) .* Float64.(surface_q))
    all(isfinite, temperature_virtual) && all(>(0), temperature_virtual) ||
        error("Invalid surface virtual temperature")
    result = exp.(Float64.(lnps) .+ Float64(gravity) .* Float64.(orography_m) ./
        (Float64(gas_constant_dry) .* temperature_virtual)) ./ 100
    all(isfinite, result) && all(>(0), result) || error("Invalid reduced pressure")
    return Float32.(result)
end

function _scored_sixhour_pressure_fields(u, v, temperature, ps_pa, sigma, levels_pa)
    n, layers = size(temperature)
    size(u) == size(v) == (n, layers) && length(ps_pa) == n && length(sigma) == layers ||
        error("Pressure-field axes differ")
    all(isfinite, u) && all(isfinite, v) && all(isfinite, temperature) &&
        all(isfinite, ps_pa) && all(>(0), ps_pa) || error("Invalid pressure-field input")
    shape = (n, length(levels_pa))
    fields = (u = fill(Float32(NaN), shape), v = fill(Float32(NaN), shape),
        temperature = fill(Float32(NaN), shape), valid = zeros(Int8, shape))
    for i in 1:n, j in eachindex(levels_pa)
        target = levels_pa[j]
        t = ScoredMonthlyCore.pressure_temperature(view(temperature, i, :), sigma, ps_pa[i], target)
        uu = ScoredMonthlyCore.pressure_temperature(view(u, i, :), sigma, ps_pa[i], target)
        vv = ScoredMonthlyCore.pressure_temperature(view(v, i, :), sigma, ps_pa[i], target)
        isnothing(t) == isnothing(uu) == isnothing(vv) || error("Pressure masks differ")
        isnothing(t) && continue
        fields.temperature[i, j] = Float32(t)
        fields.u[i, j] = Float32(uu)
        fields.v[i, j] = Float32(vv)
        fields.valid[i, j] = 1
    end
    return fields
end

function _scored_sample_sixhourly!(capture::ScoredSixHourlyCapture, sim)
    time = sim.model.atmosphere.variables.prognostic.clock.time
    expected = capture.start + Hour(6 * (length(capture.samples) + 1))
    time == expected && time <= capture.stop || error("Missing or duplicate six-hourly sample")
    model = sim.model.atmosphere.model
    vars = sim.model.atmosphere.variables
    u, v = _scored_host(vars.grid.u), _scored_host(vars.grid.v)
    temperature = _scored_host(vars.grid.temperature)
    lnps = vec(_scored_host(vars.grid.pressure))
    ps = exp.(Float64.(lnps))
    pressure = _scored_sixhour_pressure_fields(u, v, temperature, ps,
        capture.sigma, capture.levels_pa)
    surface_t = vec(_scored_host(vars.parameterizations.surface_air_temperature))
    surface_q = view(_scored_host(vars.grid.humidity), :, size(temperature, 2))
    orography = vec(_scored_host(model.orography.orography))
    mslp = _scored_reduced_mslp_hpa(lnps, orography, surface_t, surface_q,
        model.planet.gravity, model.atmosphere.R_dry, model.atmosphere.μ_virt_temp)
    push!(capture.samples, (time_seconds = Float64(sim.model.clock.time),
        u = pressure.u, v = pressure.v, temperature = pressure.temperature,
        valid = pressure.valid, mslp = mslp))
    return nothing
end

function _scored_save_sixhourly(capture::ScoredSixHourlyCapture, path)
    ispath(path) && error("Six-hourly output exists")
    expected = Dates.value(capture.stop - capture.start) ÷ (6 * 3600000)
    length(capture.samples) == expected || error("Incomplete six-hourly coverage")
    n, levels = length(capture.weights), length(capture.levels_pa)
    samples = length(capture.samples)
    wind_u = Array{Float32}(undef, n, levels, samples)
    wind_v = similar(wind_u)
    temperature = similar(wind_u)
    valid = Array{Int8}(undef, n, levels, samples)
    mslp = Array{Float32}(undef, n, samples)
    times = Vector{Float64}(undef, samples)
    for (j, sample) in enumerate(capture.samples)
        wind_u[:, :, j] = sample.u
        wind_v[:, :, j] = sample.v
        temperature[:, :, j] = sample.temperature
        valid[:, :, j] = sample.valid
        mslp[:, j] = sample.mslp
        times[j] = sample.time_seconds
    end
    times == collect(21600.0:21600.0:21600.0 * samples) || error("Six-hour timestamps differ")
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "point", n)
        NCDatasets.defDim(ds, "pressure", levels)
        NCDatasets.defDim(ds, "time", samples)
        ds.attrib["schema"] = "resolution_sixhourly_v1"
        ds.attrib["sampling"] = "Completed six-hour instants, launch exclusive, final day inclusive"
        ds.attrib["pressure_interpolation"] = "Linear in log pressure; invalid below ground"
        ds.attrib["mslp_method"] = "SpeedyWeather MeanSeaLevelPressureOutput reduction using surface virtual temperature"
        for (name, values, dims, units) in (
                ("longitude", capture.longitude, ("point",), "degrees_east"),
                ("latitude", capture.latitude, ("point",), "degrees_north"),
                ("area_weight", capture.weights, ("point",), "1"),
                ("pressure", capture.levels_pa, ("pressure",), "Pa"),
                ("time_seconds", times .+ Dates.value(capture.start - DateTime(1993, 1, 1)) / 1000,
                    ("time",), "seconds since 1993-01-01 00:00:00"),
                ("u", wind_u, ("point", "pressure", "time"), "m s-1"),
                ("v", wind_v, ("point", "pressure", "time"), "m s-1"),
                ("temperature", temperature, ("point", "pressure", "time"), "K"),
                ("valid", valid, ("point", "pressure", "time"), "1"),
                ("mslp", mslp, ("point", "time"), "hPa"))
            variable = NCDatasets.defVar(ds, name, eltype(values), dims;
                attrib = Dict("units" => units))
            variable[ntuple(_ -> Colon(), ndims(values))...] = values
        end
    end
    return Dict("file" => basename(path), "samples" => samples,
        "pressure_levels_pa" => capture.levels_pa,
        "mslp_min_hpa" => minimum(mslp), "mslp_max_hpa" => maximum(mslp))
end

# The coupled-simulation side: one Oceananigans callback that takes the initial skill sample at iteration 0, the
# hourly top-of-atmosphere sample and the six-hourly sample, in that order, as the release harness did.
mutable struct ScoredOutputs
    started::Bool
    skill::ScoredSkillObserver
    toa::ScoredTOACapture
    sixhour::ScoredSixHourlyCapture
    start::DateTime
    stop::DateTime
    steps_per_hour::Int
    steps::Int
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    weights::Vector{Float64}
end

function (s::ScoredOutputs)(sim)
    i = sim.model.clock.iteration
    # The captures average from the start of the window and are not saved in checkpoints. So a run that does not
    # start at iteration 0 (any restart or pickup route) is refused before its first step.
    if !s.started
        i == 0 || throw(ArgumentError(
            "scored_outputs cannot continue a restarted run. Run the full window from its start, or switch scored_outputs off."))
        s.started = true
    end
    if i == 0
        _scored_initialize_skill!(s.skill, sim.model.atmosphere.variables, sim.model.atmosphere.model)
    end
    i > 0 && i % s.steps_per_hour == 0 && _scored_sample_toa!(s.toa, sim.model.atmosphere.model)
    if i > 0 && i % (6 * s.steps_per_hour) == 0
        before = _scored_input_hashes(sim.model.atmosphere.variables)
        _scored_sample_sixhourly!(s.sixhour, sim)
        before == _scored_input_hashes(sim.model.atmosphere.variables) || error("Six-hourly sampling changed its inputs")
    end
    return nothing
end

const SCORED_OUTPUTS_CALLBACK = :readyesm_scored_outputs

"""Install the scored-output captures on a coupled simulation. Called by `build_dynamic_esm` when
`config.scored_outputs` is true."""
function install_scored_outputs!(simulation, config::ExperimentConfig)
    model = simulation.model.atmosphere.model
    step_ms = Dates.value(model.time_stepping.Δt_millisec)
    0 < step_ms <= 3600000 && 3600000 % step_ms == 0 || throw(ArgumentError(
        "scored_outputs requires an atmosphere timestep that divides one hour"))
    start = config.start_date
    duration_ms = round(Int, config.duration_days * 86400000)
    duration_ms % (6 * 3600000) == 0 || throw(ArgumentError(
        "scored_outputs requires duration_days to be a whole number of six-hour periods"))
    stop = start + Millisecond(duration_ms)
    weights = Float64.(Array(_global_point_weights(model.spectral_grid)))
    longitude, latitude = SpeedyWeather.RingGrids.get_londlatds(model.spectral_grid.grid)
    longitude, latitude = Float64.(Array(longitude)), Float64.(Array(latitude))
    monthly = ScoredMonthlyCapture(
        ScoredMonthlyCore.Accumulator(start, stop, length(weights), Array(model.geometry.σ_levels_full),
            Array(model.geometry.σ_levels_thick), model.planet.gravity, model.atmosphere.water_density;
            step_ms = step_ms),
        longitude, latitude, weights)
    wind = ScoredWindCapture(ScoredWindCore.Accumulator(start, stop, length(weights), Array(model.geometry.σ_levels_full)),
        longitude, latitude, weights)
    skill = ScoredSkillObserver(monthly, wind, false, 0, 0, String[])
    SpeedyWeather.add!(model.callbacks, SCORED_OUTPUTS_CALLBACK => skill)
    sixhour = ScoredSixHourlyCapture(start, stop, Array(model.geometry.σ_levels_full), longitude, latitude, weights)
    outputs = ScoredOutputs(false, skill, ScoredTOACapture(length(weights)), sixhour, start, stop,
        3600000 ÷ step_ms, duration_ms ÷ step_ms, longitude, latitude, weights)
    Oceananigans.add_callback!(simulation, outputs, Oceananigans.IterationInterval(1); name = SCORED_OUTPUTS_CALLBACK)
    return outputs
end

"""Write monthly.nc, wind.nc, toa_monthly.nc and eddy_6hourly.nc into `config.output_dir` after a completed run, and
return a summary. Refuses an incomplete run."""
function save_scored_outputs(simulation, config::ExperimentConfig)
    haskey(simulation.callbacks, SCORED_OUTPUTS_CALLBACK) || error("scored outputs were not installed")
    s = simulation.callbacks[SCORED_OUTPUTS_CALLBACK].func
    s isa ScoredOutputs || error("unexpected scored-output callback")
    simulation.model.clock.iteration == s.steps || error("scored outputs need the complete run")
    skill = s.skill
    skill.initialized || error("Missing initial sample")
    skill.callbacks == s.steps || error("Wrong completed-step count")
    skill.passivity_checks == s.steps + 1 || error("Incomplete input passivity checks")
    output = config.output_dir
    for name in ("toa_monthly.nc", "monthly.nc", "wind.nc", "eddy_6hourly.nc")
        ispath(joinpath(output, name)) && error("Scored output exists: $(name)")
    end
    mkpath(output)
    hours = Dates.value(s.stop - s.start) ÷ 3600000
    toa = _scored_save_toa(s.toa, joinpath(output, "toa_monthly.nc"), s.longitude, s.latitude, s.weights,
        s.start, s.stop)
    toa["samples"] == hours || error("incomplete TOA capture")
    monthly = _scored_save_monthly(skill.monthly, joinpath(output, "monthly.nc"))
    wind = _scored_save_wind(skill.wind, joinpath(output, "wind.nc"))
    monthly["state_samples"] == wind["state_samples"] == hours || error("Incomplete hourly state coverage")
    monthly["flux_samples"] == s.steps || error("Incomplete precipitation coverage")
    eddy = _scored_save_sixhourly(s.sixhour, joinpath(output, "eddy_6hourly.nc"))
    return Dict("monthly" => monthly, "wind" => wind, "toa" => toa, "sixhourly" => eddy,
        "callbacks" => skill.callbacks, "passivity_checks" => skill.passivity_checks,
        "initial_input_hashes" => skill.initial_hashes)
end
