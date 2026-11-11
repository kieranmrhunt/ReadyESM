module ReadyESMMonthlyWindows

using Dates, NCDatasets, ReadyESM, SHA, TOML
const R = ReadyESM
const SW = R.SpeedyWeather
const O = R.Oceananigans

const CALLBACK = :partc_calendar_month_outputs
include("cloud_profile.jl")

struct MonthWindow
    start::DateTime
    stop::DateTime
    monthly::R.ScoredMonthlyCapture
    wind::R.ScoredWindCapture
    toa::R.ScoredTOACapture
    cloud::CloudProfileCapture
end

mutable struct MonthObserver <: SW.AbstractCallback
    windows::Vector{MonthWindow}
    run_start::DateTime
    run_stop::DateTime
    step_ms::Int
    steps_per_hour::Int
    total_steps::Int
    started::Bool
    callbacks::Int
    toa_samples::Int
    output_root::String
    records::Dict{String, Any}
end

"""Return only complete calendar months contained in [start, stop]."""
function complete_months(start::DateTime, stop::DateTime)
    start < stop || throw(ArgumentError("run window is empty"))
    first = DateTime(year(start), month(start), 1)
    first < start && (first += Month(1))
    result = Tuple{DateTime, DateTime}[]
    while first + Month(1) <= stop
        last = first + Month(1)
        push!(result, (first, last))
        first = last
    end
    return result
end

function _window_index(windows::Vector{MonthWindow}, time::DateTime)
    for (index, window) in pairs(windows)
        window.start <= time < window.stop && return index
    end
    return nothing
end

function _check_windows(start::DateTime, stop::DateTime, intervals)
    isempty(intervals) && throw(ArgumentError("no complete output window"))
    previous_stop = start
    for (a, b) in intervals
        start <= a < b <= stop || throw(ArgumentError("output window outside run"))
        a >= previous_stop || throw(ArgumentError("overlapping or unordered output windows"))
        Dates.value(b - a) % 3_600_000 == 0 || throw(ArgumentError("output window lacks whole hours"))
        previous_stop = b
    end
    return nothing
end

function _window_label(window::MonthWindow)
    if day(window.start) == 1 && hour(window.start) == 0 &&
       minute(window.start) == 0 && second(window.start) == 0 &&
       window.stop == window.start + Month(1)
        return Dates.format(window.start, dateformat"yyyy-mm")
    end
    # Short explicit windows exist only for the diagnostic gate.
    return Dates.format(window.start, dateformat"yyyy-mm-ddTHHMM") * "_" *
           Dates.format(window.stop, dateformat"yyyy-mm-ddTHHMM")
end

"""Install independent diagnostic-only windows; explicit intervals support short gate fixtures."""
function install_windows!(sim, config::R.ExperimentConfig, intervals; resume::Bool=false)
    config.scored_outputs || throw(ArgumentError("whole-run scored_outputs must remain enabled"))
    haskey(sim.callbacks, R.SCORED_OUTPUTS_CALLBACK) || error("whole-run scored output callback missing")
    haskey(sim.callbacks, CALLBACK) && error("month observer already installed")
    model = sim.model.atmosphere.model
    step_ms = Dates.value(model.time_stepping.Δt_millisec)
    0 < step_ms <= 3_600_000 && 3_600_000 % step_ms == 0 ||
        throw(ArgumentError("timestep must divide one hour"))
    run_start = config.start_date
    duration_ms = round(Int, config.duration_days * 86_400_000)
    run_stop = run_start + Millisecond(duration_ms)
    duration_ms % step_ms == 0 || throw(ArgumentError("run duration is not a whole number of steps"))
    _check_windows(run_start, run_stop, intervals)

    weights = Float64.(Array(R._global_point_weights(model.spectral_grid)))
    longitude, latitude = SW.RingGrids.get_londlatds(model.spectral_grid.grid)
    longitude, latitude = Float64.(Array(longitude)), Float64.(Array(latitude))
    sigma = Array(model.geometry.σ_levels_full)
    dsigma = Array(model.geometry.σ_levels_thick)
    windows = MonthWindow[]
    for (start, stop) in intervals
        (Dates.value(start - run_start) % step_ms == 0 &&
         Dates.value(stop - run_start) % step_ms == 0) ||
            throw(ArgumentError("output boundary is not a model step"))
        monthly = R.ScoredMonthlyCapture(
            R.ScoredMonthlyCore.Accumulator(start, stop, length(weights), sigma, dsigma,
                model.planet.gravity, model.atmosphere.water_density; step_ms),
            longitude, latitude, weights)
        wind = R.ScoredWindCapture(
            R.ScoredWindCore.Accumulator(start, stop, length(weights), sigma),
            longitude, latitude, weights)
        radiation = model.longwave_radiation
        scheme = hasproperty(radiation, :cloud_scheme) ? string(radiation.cloud_scheme) : "none"
        cloud = CloudProfileCapture(start, stop, weights, sigma, dsigma, scheme)
        push!(windows, MonthWindow(start, stop, monthly, wind,
            R.ScoredTOACapture(length(weights)), cloud))
    end
    output_root = joinpath(abspath(config.output_dir), "months")
    if resume
        isdir(output_root) || error("month output directory missing at pickup")
    else
        ispath(output_root) && error("refusing existing month output directory")
        mkpath(output_root)
    end
    observer = MonthObserver(windows, run_start, run_stop, step_ms, 3_600_000 ÷ step_ms,
        duration_ms ÷ step_ms, false, 0, 0, output_root, Dict{String, Any}())
    resume && reconcile_committed_windows!(observer)
    SW.add!(model.callbacks, CALLBACK => observer)
    O.add_callback!(sim, observer, O.IterationInterval(1); name=CALLBACK)
    return observer
end

"""Install one output window per complete calendar month within a continuous run."""
function install_calendar_months!(sim, config::R.ExperimentConfig; resume::Bool=false)
    start = config.start_date
    stop = start + Millisecond(round(Int, config.duration_days * 86_400_000))
    return install_windows!(sim, config, complete_months(start, stop); resume)
end

function _state!(observer::MonthObserver, vars, model, time::DateTime)
    index = _window_index(observer.windows, time)
    isnothing(index) && return nothing
    window = observer.windows[index]
    if window.monthly.accumulator.state_count == 0
        time == window.start || error("first state missed for output window $index")
        SW.initialize!(window.monthly, vars, model)
        SW.initialize!(window.wind, vars, model)
    else
        R._scored_monthly_state!(window.monthly, vars)
        R._scored_wind_state!(window.wind, vars)
    end
    _cloud_sample!(window.cloud, vars, model, time)
    return nothing
end

function SW.callback!(observer::MonthObserver, vars::SW.Variables, model::SW.AbstractModel)
    observer.started || error("month observer was not initialized")
    observer.callbacks += 1
    time = vars.prognostic.clock.time
    time == observer.run_start + Millisecond(observer.callbacks * observer.step_ms) ||
        error("missing or duplicate coupled step in month observer")
    flux_index = _window_index(observer.windows, time - Millisecond(1))
    if !isnothing(flux_index)
        a = observer.windows[flux_index].monthly.accumulator
        R.ScoredMonthlyCore.flux!(a, time,
            vec(R._scored_host(vars.parameterizations.rain_rate)),
            vec(R._scored_host(vars.parameterizations.snow_rate)))
    end
    if time < observer.run_stop && Dates.value(time - observer.run_start) % 3_600_000 == 0
        _state!(observer, vars, model, time)
    end
    return nothing
end

SW.finalize!(::MonthObserver, args...) = nothing

function (observer::MonthObserver)(sim)
    iteration = Int(sim.model.clock.iteration)
    if !observer.started
        iteration == 0 || error("month observer cannot start from a pickup")
        observer.started = true
        sim.model.atmosphere.variables.prognostic.clock.time == observer.run_start ||
            error("wrong month observer start")
        _state!(observer, sim.model.atmosphere.variables,
            sim.model.atmosphere.model, observer.run_start)
    end
    iteration == observer.callbacks || error("month observer step count differs from coupled clock")
    if iteration > 0 && iteration % observer.steps_per_hour == 0
        time = sim.model.atmosphere.variables.prognostic.clock.time
        index = _window_index(observer.windows, time - Millisecond(1))
        if !isnothing(index)
            R._scored_sample_toa!(observer.windows[index].toa, sim.model.atmosphere.model)
            observer.toa_samples += 1
        end
    end
    if iteration > 0
        time = sim.model.atmosphere.variables.prognostic.clock.time
        for (index, window) in pairs(observer.windows)
            window.stop == time && _commit_window!(observer, index)
        end
    end
    return nothing
end

"""Atomically refresh the inventory after a committed window is visible."""
function _write_window_inventory!(observer::MonthObserver)
    records = observer.records
    result = Dict(
        "schema" => "readyesm_monthly_windows_v2",
        "source_commit" => R._dynamic_source_commit(),
        "run_start" => string(observer.run_start),
        "run_stop" => string(observer.run_stop),
        "run_steps" => observer.total_steps,
        "window_hours" => sum(record["hours"] for record in values(records); init=0),
        "window_steps" => sum(record["steps"] for record in values(records); init=0),
        "committed_count" => length(records),
        "windows" => records,
    )
    temporary = joinpath(observer.output_root,
        "." * basename(tempname(observer.output_root)) * ".inventory-incomplete")
    open(temporary, "w") do io
        TOML.print(io, result; sorted=true)
    end
    mv(temporary, joinpath(observer.output_root, "windows.toml"); force=true)
    return result
end

"""Verify sealed month directories and reconcile a crash after rename, before inventory update."""
function reconcile_committed_windows!(observer::MonthObserver)
    root = observer.output_root
    isdir(root) || error("month output directory missing")
    labels = Dict(_window_label(window) => (index, window) for
        (index, window) in pairs(observer.windows))
    length(labels) == length(observer.windows) || error("duplicate month labels")
    actual = Dict{String, Any}()
    for label in readdir(root)
        path = joinpath(root, label)
        label == "windows.toml" && continue
        startswith(label, ".") && continue # abandoned staging directories are invalid months
        haskey(labels, label) && isdir(path) || error("unknown month output entry $label")
        index, window = labels[label]
        receipt_path = joinpath(path, "receipt.toml")
        isfile(receipt_path) || error("sealed month $label has no receipt")
        receipt = TOML.parsefile(receipt_path)
        hours = Dates.value(window.stop - window.start) ÷ 3_600_000
        steps = Dates.value(window.stop - window.start) ÷ observer.step_ms
        expected = Dict(
            "schema" => "readyesm_monthly_window_receipt_v2",
            "source_commit" => R._dynamic_source_commit(),
            "start" => string(window.start), "stop" => string(window.stop),
            "hours" => hours, "steps" => steps,
            "stop_iteration" => Dates.value(window.stop - observer.run_start) ÷ observer.step_ms,
        )
        all(haskey(receipt, key) && isequal(receipt[key], value) for
            (key, value) in expected) || error("sealed month $label identity differs")
        for (file, key) in (("monthly.nc", "monthly_sha256"),
                            ("wind.nc", "wind_sha256"),
                            ("toa_monthly.nc", "toa_sha256"),
                            ("cloud_profile.nc", "cloud_profile_sha256"))
            target = joinpath(path, file)
            isfile(target) && haskey(receipt, key) &&
                bytes2hex(open(sha256, target)) == receipt[key] ||
                error("sealed month $label file $file differs from receipt")
        end
        actual[label] = merge(receipt, Dict(
            "receipt_sha256" => bytes2hex(open(sha256, receipt_path))))
    end
    missing_seen = false
    for window in observer.windows
        label = _window_label(window)
        missing_seen && haskey(actual, label) &&
            error("sealed months are not a prefix of the requested windows")
        missing_seen |= !haskey(actual, label)
    end
    inventory_path = joinpath(root, "windows.toml")
    if isfile(inventory_path)
        inventory = TOML.parsefile(inventory_path)
        expected = Dict("schema" => "readyesm_monthly_windows_v2",
            "source_commit" => R._dynamic_source_commit(),
            "run_start" => string(observer.run_start),
            "run_stop" => string(observer.run_stop),
            "run_steps" => observer.total_steps)
        all(haskey(inventory, key) && isequal(inventory[key], value) for
            (key, value) in expected) || error("month inventory identity differs")
        recorded = get(inventory, "windows", Dict{String, Any}())
        all(haskey(actual, label) && isequal(record, actual[label]) for
            (label, record) in recorded) ||
            error("month inventory conflicts with sealed receipts")
        get(inventory, "committed_count", -1) == length(recorded) ||
            error("month inventory count differs")
        get(inventory, "window_hours", -1) ==
            sum(record["hours"] for record in values(recorded); init=0) &&
        get(inventory, "window_steps", -1) ==
            sum(record["steps"] for record in values(recorded); init=0) ||
            error("month inventory duration differs")
    end
    empty!(observer.records)
    merge!(observer.records, actual)
    # Only write when the aggregate record is absent or lags a verified rename.
    if !isfile(inventory_path) ||
            length(get(TOML.parsefile(inventory_path), "windows", Dict())) != length(actual)
        _write_window_inventory!(observer)
    end
    return length(actual)
end

const RESTART_SCHEMA = "readyesm_monthly_observer_restart_v1"
_restart_require(condition, message) = condition ||
    throw(ArgumentError("month observer checkpoint: " * message))

function _window_static(window::MonthWindow)
    monthly, wind, cloud = window.monthly.accumulator,
        window.wind.accumulator, window.cloud
    return (
        start = window.start, stop = window.stop,
        longitude = copy(window.monthly.longitude),
        latitude = copy(window.monthly.latitude),
        weights = copy(window.monthly.weights),
        monthly_sigma = copy(monthly.sigma), monthly_dsigma = copy(monthly.dsigma),
        monthly_levels_pa = copy(monthly.levels_Pa),
        monthly_gravity = monthly.gravity, monthly_water_density = monthly.water_density,
        wind_sigma = copy(wind.sigma), wind_levels = copy(wind.levels),
        cloud_weights = copy(cloud.weights), cloud_sigma = copy(cloud.sigma),
        cloud_dsigma = copy(cloud.dsigma), cloud_scheme = cloud.scheme,
    )
end

function _observer_restart_state(observer::MonthObserver)
    return (
        schema = RESTART_SCHEMA,
        static = (
            source_commit = R._dynamic_source_commit(),
            run_start = observer.run_start, run_stop = observer.run_stop,
            step_ms = observer.step_ms, steps_per_hour = observer.steps_per_hour,
            total_steps = observer.total_steps,
            windows = [_window_static(window) for window in observer.windows],
        ),
        started = observer.started,
        callbacks = observer.callbacks,
        toa_samples = observer.toa_samples,
        records = deepcopy(observer.records),
        windows = [begin
            m, w, t, c = window.monthly.accumulator, window.wind.accumulator,
                window.toa, window.cloud
            (
                monthly = (
                    state_count = m.state_count, flux_count = m.flux_count,
                    water_sum = copy(m.water_sum), ps_sum = copy(m.ps_sum),
                    temperature_sum = copy(m.temperature_sum),
                    temperature_count = copy(m.temperature_count),
                    precipitation_kg_m2 = copy(m.precipitation_kg_m2),
                ),
                wind = (
                    samples = w.samples, u_sum = copy(w.u_sum),
                    v_sum = copy(w.v_sum), count = copy(w.count),
                ),
                toa = (
                    count = t.count, clear_sky = t.clear_sky,
                    sums = Dict(key => copy(t.sums[key]) for key in R.SCORED_TOA_FIELDS),
                ),
                cloud = (
                    samples = c.samples, fraction_sum = copy(c.fraction_sum),
                    liquid_sum = copy(c.liquid_sum), ice_sum = copy(c.ice_sum),
                    top_sum = copy(c.top_sum),
                ),
            )
        end for window in observer.windows],
    )
end

function _validate_observer_restart_state(observer::MonthObserver, state)
    _restart_require(hasproperty(state, :schema) && state.schema == RESTART_SCHEMA &&
        hasproperty(state, :static) && hasproperty(state, :windows) &&
        hasproperty(state, :records) && hasproperty(state, :started) &&
        hasproperty(state, :callbacks) && hasproperty(state, :toa_samples),
        "missing or unknown payload")
    current = _observer_restart_state(observer).static
    _restart_require(isequal(state.static, current),
        "source, grid, windows or constants differ")
    i = state.callbacks
    _restart_require(i isa Integer && 0 <= i <= observer.total_steps &&
        state.started === true && length(state.windows) == length(observer.windows),
        "iteration or window count differs")
    time = observer.run_start + Millisecond(i * observer.step_ms)
    expected_toa = 0
    completed = String[]
    for (index, window) in pairs(observer.windows)
        payload = state.windows[index]
        _restart_require(hasproperty(payload, :monthly) && hasproperty(payload, :wind) &&
            hasproperty(payload, :toa) && hasproperty(payload, :cloud),
            "window $index payload incomplete")
        m, w, t, c = payload.monthly, payload.wind, payload.toa, payload.cloud
        elapsed = Dates.value(time - window.start)
        hours = Dates.value(window.stop - window.start) ÷ 3_600_000
        steps = Dates.value(window.stop - window.start) ÷ observer.step_ms
        flux_count = clamp(fld(elapsed, observer.step_ms), 0, steps)
        state_count = elapsed < 0 ? 0 : min(hours, fld(elapsed, 3_600_000) + 1)
        toa_count = clamp(fld(elapsed, 3_600_000), 0, hours)
        expected_toa += toa_count
        _restart_require(m.state_count == w.samples == c.samples == state_count &&
            m.flux_count == flux_count && t.count == toa_count &&
            t.clear_sky isa Bool, "window $index sample counts differ from checkpoint time")
        for (src, dst, label) in (
            (m.water_sum, window.monthly.accumulator.water_sum, "water"),
            (m.ps_sum, window.monthly.accumulator.ps_sum, "pressure"),
            (m.temperature_sum, window.monthly.accumulator.temperature_sum, "temperature"),
            (m.temperature_count, window.monthly.accumulator.temperature_count, "temperature count"),
            (m.precipitation_kg_m2, window.monthly.accumulator.precipitation_kg_m2, "precipitation"),
            (w.u_sum, window.wind.accumulator.u_sum, "wind u"),
            (w.v_sum, window.wind.accumulator.v_sum, "wind v"),
            (w.count, window.wind.accumulator.count, "wind count"),
            (c.fraction_sum, window.cloud.fraction_sum, "cloud fraction"),
            (c.liquid_sum, window.cloud.liquid_sum, "cloud liquid"),
            (c.ice_sum, window.cloud.ice_sum, "cloud ice"),
            (c.top_sum, window.cloud.top_sum, "cloud top"),
        )
            _restart_require(size(src) == size(dst) && eltype(src) == eltype(dst) &&
                all(isfinite, src), "window $index $label shape, type or values differ")
        end
        _restart_require(Set(keys(t.sums)) == Set(R.SCORED_TOA_FIELDS),
            "window $index TOA field set differs")
        for key in R.SCORED_TOA_FIELDS
            src, dst = t.sums[key], window.toa.sums[key]
            _restart_require(size(src) == size(dst) && eltype(src) == eltype(dst) &&
                all(isfinite, src), "window $index TOA $key differs")
        end
        window.stop <= time && push!(completed, _window_label(window))
    end
    _restart_require(state.toa_samples == expected_toa,
        "aggregate TOA count differs from checkpoint time")
    _restart_require(Set(keys(state.records)) == Set(completed) &&
        isequal(state.records, observer.records),
        "committed receipt set differs from checkpoint or disk")
    return nothing
end

function _restore_observer_restart_state!(observer::MonthObserver, state)
    _validate_observer_restart_state(observer, state)
    for (window, payload) in zip(observer.windows, state.windows)
        m, w, t, c = window.monthly.accumulator, window.wind.accumulator,
            window.toa, window.cloud
        for (dst, src) in (
            (m.water_sum, payload.monthly.water_sum),
            (m.ps_sum, payload.monthly.ps_sum),
            (m.temperature_sum, payload.monthly.temperature_sum),
            (m.temperature_count, payload.monthly.temperature_count),
            (m.precipitation_kg_m2, payload.monthly.precipitation_kg_m2),
            (w.u_sum, payload.wind.u_sum), (w.v_sum, payload.wind.v_sum),
            (w.count, payload.wind.count),
            (c.fraction_sum, payload.cloud.fraction_sum),
            (c.liquid_sum, payload.cloud.liquid_sum),
            (c.ice_sum, payload.cloud.ice_sum), (c.top_sum, payload.cloud.top_sum),
        )
            copyto!(dst, src)
        end
        m.state_count, m.flux_count = payload.monthly.state_count, payload.monthly.flux_count
        w.samples = payload.wind.samples
        for key in R.SCORED_TOA_FIELDS
            copyto!(t.sums[key], payload.toa.sums[key])
        end
        t.count, t.clear_sky = payload.toa.count, payload.toa.clear_sky
        c.samples = payload.cloud.samples
    end
    observer.started = state.started
    observer.callbacks = state.callbacks
    observer.toa_samples = state.toa_samples
    return observer
end

function O.prognostic_state(callback::O.Callback{P, F, S, CS}) where
    {P, F <: MonthObserver, S, CS}
    return (schedule = O.prognostic_state(callback.schedule),
        func = _observer_restart_state(callback.func))
end

function O.restore_prognostic_state!(callback::O.Callback{P, F, S, CS}, state) where
    {P, F <: MonthObserver, S, CS}
    _restart_require(hasproperty(state, :schedule) && hasproperty(state, :func),
        "month callback payload missing")
    _validate_observer_restart_state(callback.func, state.func)
    O.restore_prognostic_state!(callback.schedule, state.schedule)
    _restore_observer_restart_state!(callback.func, state.func)
    return callback
end

O.restore_prognostic_state!(::O.Callback{P, F, S, CS}, ::Nothing) where
    {P, F <: MonthObserver, S, CS} =
    throw(ArgumentError("month observer checkpoint: callback payload missing"))

"""Restore the coupled model and reject a month-observer clock mismatch before stepping."""
function restore_dynamic_restart_state!(sim, path::AbstractString)
    haskey(sim.callbacks, CALLBACK) || error("month observer callback missing")
    R.restore_dynamic_restart_state!(sim, path)
    observer = sim.callbacks[CALLBACK].func
    observer isa MonthObserver && observer.started &&
        observer.callbacks == Int(sim.model.clock.iteration) ||
        error("month observer count differs from restored coupled clock")
    return sim
end

"""Write one complete month through a staging directory and sealed receipt."""
function _commit_window!(observer::MonthObserver, index::Int)
    window = observer.windows[index]
    label = _window_label(window)
    haskey(observer.records, label) && return false
    target = joinpath(observer.output_root, label)
    ispath(target) && error("refusing existing unrecorded month directory $target")
    hours = Dates.value(window.stop - window.start) ÷ 3_600_000
    steps = Dates.value(window.stop - window.start) ÷ observer.step_ms
    window.monthly.accumulator.state_count == window.wind.accumulator.samples ==
        window.toa.count == hours || error("$label has incomplete hourly samples")
    window.monthly.accumulator.flux_count == steps ||
        error("$label has incomplete flux steps")
    _cloud_finish(window.cloud)
    stage = mktempdir(observer.output_root; prefix=".$label-incomplete-")
    monthly = R._scored_save_monthly(window.monthly, joinpath(stage, "monthly.nc"))
    wind = R._scored_save_wind(window.wind, joinpath(stage, "wind.nc"))
    toa = R._scored_save_toa(window.toa, joinpath(stage, "toa_monthly.nc"),
        window.monthly.longitude, window.monthly.latitude, window.monthly.weights,
        window.start, window.stop)
    cloud = _save_cloud_profile(window.cloud, joinpath(stage, "cloud_profile.nc"))
    monthly["state_samples"] == wind["state_samples"] == toa["samples"] ==
        cloud["state_samples"] == hours || error("$label output sample count differs")
    monthly["flux_samples"] == steps || error("$label output flux count differs")
    record = Dict{String, Any}(
        "schema" => "readyesm_monthly_window_receipt_v2",
        "source_commit" => R._dynamic_source_commit(),
        "start" => string(window.start), "stop" => string(window.stop),
        "hours" => hours, "steps" => steps,
        "stop_iteration" => Dates.value(window.stop - observer.run_start) ÷ observer.step_ms,
        "monthly_sha256" => bytes2hex(open(sha256, joinpath(stage, "monthly.nc"))),
        "wind_sha256" => bytes2hex(open(sha256, joinpath(stage, "wind.nc"))),
        "toa_sha256" => bytes2hex(open(sha256, joinpath(stage, "toa_monthly.nc"))),
        "cloud_profile_sha256" => bytes2hex(open(sha256, joinpath(stage, "cloud_profile.nc"))),
    )
    open(joinpath(stage, "receipt.toml"), "w") do io
        TOML.print(io, record; sorted=true)
    end
    receipt_sha = bytes2hex(open(sha256, joinpath(stage, "receipt.toml")))
    mv(stage, target)
    observer.records[label] = merge(record, Dict("receipt_sha256" => receipt_sha))
    _write_window_inventory!(observer)
    println("V021_MONTH_COMMITTED label=$label receipt_sha256=$receipt_sha")
    flush(stdout)
    return true
end

"""Verify and return all months after the complete run; files were flushed at boundaries."""
function save_windows!(observer::MonthObserver, sim, output::AbstractString)
    observer.started && observer.callbacks == observer.total_steps &&
        sim.model.clock.iteration == observer.total_steps || error("month observer run incomplete")
    joinpath(abspath(output), "months") == observer.output_root ||
        error("month output directory differs from installed observer")
    for index in eachindex(observer.windows)
        _commit_window!(observer, index)
    end
    result = _write_window_inventory!(observer)
    result["committed_count"] == length(observer.windows) ||
        error("not all month windows were committed")
    sum(window.toa.count for window in observer.windows) == observer.toa_samples ==
        result["window_hours"] || error("TOA sample inventory differs")
    if first(observer.windows).start == observer.run_start &&
       last(observer.windows).stop == observer.run_stop
        result["window_steps"] == observer.total_steps ||
            error("calendar windows do not tile the run steps")
    end
    return result
end

"""Detach only a completed, correctly registered observer before restart serialization."""
function detach_completed!(observer::MonthObserver, sim)
    observer.started && observer.callbacks == observer.total_steps &&
        sim.model.clock.iteration == observer.total_steps ||
        error("month observer cannot detach before the complete run")
    model_callbacks = sim.model.atmosphere.model.callbacks
    haskey(sim.callbacks, CALLBACK) &&
        sim.callbacks[CALLBACK].func === observer ||
        error("coupled month observer registration differs")
    haskey(model_callbacks, CALLBACK) && model_callbacks[CALLBACK] === observer ||
        error("host month observer registration differs")
    total_hours = 0
    total_steps = 0
    for (index, window) in pairs(observer.windows)
        hours = Dates.value(window.stop - window.start) ÷ 3_600_000
        steps = Dates.value(window.stop - window.start) ÷ observer.step_ms
        window.monthly.accumulator.state_count == window.wind.accumulator.samples ==
            window.toa.count == hours || error("window $index hourly samples incomplete")
        window.monthly.accumulator.flux_count == steps ||
            error("window $index step samples incomplete")
        total_hours += hours
        total_steps += steps
    end
    total_hours == observer.toa_samples || error("TOA sample inventory differs at detach")
    if first(observer.windows).start == observer.run_start &&
       last(observer.windows).stop == observer.run_stop
        total_steps == observer.total_steps || error("window steps do not tile the run")
    end
    delete!(sim.callbacks, CALLBACK)
    delete!(model_callbacks, CALLBACK)
    !haskey(sim.callbacks, CALLBACK) && !haskey(model_callbacks, CALLBACK) ||
        error("month observer remains registered")
    return (window_count=length(observer.windows), hours=total_hours, steps=total_steps)
end

end # module
