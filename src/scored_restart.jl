# In-file scored-output checkpoint state. Keep this separate from the physical
# restart payload so source checkpoints without scored outputs remain readable.
const SCORED_RESTART_SCHEMA = "readiesm-scored-outputs-restart-v1"

_scored_restart_copy_sample(sample) = (
    time_seconds = sample.time_seconds,
    u = copy(sample.u),
    v = copy(sample.v),
    temperature = copy(sample.temperature),
    valid = copy(sample.valid),
    mslp = copy(sample.mslp),
)

function _scored_restart_state(s::ScoredOutputs)
    monthly = s.skill.monthly.accumulator
    wind = s.skill.wind.accumulator
    return (
        schema = SCORED_RESTART_SCHEMA,
        iteration = s.skill.callbacks,
        static = (
            start = s.start,
            stop = s.stop,
            steps_per_hour = s.steps_per_hour,
            steps = s.steps,
            longitude = copy(s.longitude),
            latitude = copy(s.latitude),
            weights = copy(s.weights),
            monthly_start = monthly.start,
            monthly_stop = monthly.stop,
            monthly_step_ms = monthly.step_ms,
            monthly_sigma = copy(monthly.sigma),
            monthly_dsigma = copy(monthly.dsigma),
            monthly_levels_pa = copy(monthly.levels_Pa),
            monthly_gravity = monthly.gravity,
            monthly_water_density = monthly.water_density,
            wind_start = wind.start,
            wind_stop = wind.stop,
            wind_sigma = copy(wind.sigma),
            wind_levels = copy(wind.levels),
            sixhour_start = s.sixhour.start,
            sixhour_stop = s.sixhour.stop,
            sixhour_sigma = copy(s.sixhour.sigma),
            sixhour_levels_pa = copy(s.sixhour.levels_pa),
            sixhour_longitude = copy(s.sixhour.longitude),
            sixhour_latitude = copy(s.sixhour.latitude),
            sixhour_weights = copy(s.sixhour.weights),
        ),
        started = s.started,
        skill = (
            initialized = s.skill.initialized,
            callbacks = s.skill.callbacks,
            passivity_checks = s.skill.passivity_checks,
            initial_hashes = copy(s.skill.initial_hashes),
        ),
        monthly = (
            state_count = monthly.state_count,
            flux_count = monthly.flux_count,
            water_sum = copy(monthly.water_sum),
            ps_sum = copy(monthly.ps_sum),
            temperature_sum = copy(monthly.temperature_sum),
            temperature_count = copy(monthly.temperature_count),
            precipitation_kg_m2 = copy(monthly.precipitation_kg_m2),
        ),
        wind = (
            samples = wind.samples,
            u_sum = copy(wind.u_sum),
            v_sum = copy(wind.v_sum),
            count = copy(wind.count),
        ),
        toa = (
            sums = Dict(k => copy(s.toa.sums[k]) for k in SCORED_TOA_FIELDS),
            count = s.toa.count,
            clear_sky = s.toa.clear_sky,
        ),
        sixhour_samples = [_scored_restart_copy_sample(x) for x in s.sixhour.samples],
    )
end

_scored_require(condition, message) = condition || throw(ArgumentError("scored checkpoint: " * message))

function _validate_scored_restart_state(s::ScoredOutputs, state, iteration)
    _scored_require(hasproperty(state, :schema) && state.schema == SCORED_RESTART_SCHEMA,
        "missing or unknown payload schema")
    _scored_require(hasproperty(state, :static) && hasproperty(state, :skill) &&
        hasproperty(state, :monthly) && hasproperty(state, :wind) &&
        hasproperty(state, :toa) && hasproperty(state, :sixhour_samples) &&
        hasproperty(state, :started) && hasproperty(state, :iteration),
        "incomplete payload")
    a, m, w, t = state.static, state.monthly, state.wind, state.toa
    current = _scored_restart_state(s).static
    _scored_require(keys(a) == keys(current) && all(isequal(getproperty(a, k), getproperty(current, k))
        for k in keys(current)), "capture window, grid or constants differ")
    i = state.iteration
    _scored_require(i isa Integer && 0 <= i <= s.steps && iteration == i,
        "iteration differs from model clock")
    _scored_require(state.started === true && state.skill.initialized === true &&
        state.skill.callbacks == i && state.skill.passivity_checks == i + 1,
        "skill callback counters differ")
    _scored_require(length(state.skill.initial_hashes) == 7 &&
        all(x -> x isa AbstractString && occursin(r"^[0-9a-f]{64}$", x), state.skill.initial_hashes),
        "initial input hashes are incomplete")
    hours = Dates.value(s.stop - s.start) ÷ 3600000
    states = 1 + min(fld(i, s.steps_per_hour), hours - 1)
    _scored_require(m.state_count == states && m.flux_count == i && w.samples == states &&
        t.count == fld(i, s.steps_per_hour) && t.clear_sky isa Bool &&
        length(state.sixhour_samples) == fld(i, 6 * s.steps_per_hour),
        "sample counts differ from model iteration")
    for (src, dst, label) in (
        (m.water_sum, s.skill.monthly.accumulator.water_sum, "water sum"),
        (m.ps_sum, s.skill.monthly.accumulator.ps_sum, "pressure sum"),
        (m.temperature_sum, s.skill.monthly.accumulator.temperature_sum, "temperature sum"),
        (m.temperature_count, s.skill.monthly.accumulator.temperature_count, "temperature count"),
        (m.precipitation_kg_m2, s.skill.monthly.accumulator.precipitation_kg_m2, "precipitation sum"),
        (w.u_sum, s.skill.wind.accumulator.u_sum, "wind u sum"),
        (w.v_sum, s.skill.wind.accumulator.v_sum, "wind v sum"),
        (w.count, s.skill.wind.accumulator.count, "wind count"),
    )
        _scored_require(size(src) == size(dst) && eltype(src) == eltype(dst), "$label shape or type differs")
        _scored_require(all(isfinite, src), "$label has nonfinite values")
    end
    _scored_require(all(>=(0), m.temperature_count) && all(>=(0), w.count),
        "negative sample count")
    _scored_require(Set(keys(t.sums)) == Set(SCORED_TOA_FIELDS), "TOA field set differs")
    for k in SCORED_TOA_FIELDS
        values = t.sums[k]
        _scored_require(size(values) == size(s.toa.sums[k]) && eltype(values) == Float64 &&
            all(isfinite, values), "TOA $k values differ in shape or are nonfinite")
    end
    n, levels = length(s.weights), length(s.sixhour.levels_pa)
    for (j, sample) in enumerate(state.sixhour_samples)
        _scored_require(sample.time_seconds == 21600.0 * j &&
            size(sample.u) == size(sample.v) == size(sample.temperature) == size(sample.valid) == (n, levels) &&
            length(sample.mslp) == n && all(x -> x == 0 || x == 1, sample.valid) &&
            all(isfinite, sample.mslp), "six-hour sample $j differs")
    end
    return nothing
end

function _restore_scored_restart_state!(s::ScoredOutputs, state)
    _validate_scored_restart_state(s, state, state.iteration)
    monthly, wind = s.skill.monthly.accumulator, s.skill.wind.accumulator
    m, w = state.monthly, state.wind
    for (dst, src) in (
        (monthly.water_sum, m.water_sum), (monthly.ps_sum, m.ps_sum),
        (monthly.temperature_sum, m.temperature_sum),
        (monthly.temperature_count, m.temperature_count),
        (monthly.precipitation_kg_m2, m.precipitation_kg_m2),
        (wind.u_sum, w.u_sum), (wind.v_sum, w.v_sum), (wind.count, w.count),
    )
        copyto!(dst, src)
    end
    monthly.state_count, monthly.flux_count = m.state_count, m.flux_count
    wind.samples = w.samples
    s.skill.initialized = state.skill.initialized
    s.skill.callbacks = state.skill.callbacks
    s.skill.passivity_checks = state.skill.passivity_checks
    empty!(s.skill.initial_hashes)
    append!(s.skill.initial_hashes, state.skill.initial_hashes)
    for k in SCORED_TOA_FIELDS
        copyto!(s.toa.sums[k], state.toa.sums[k])
    end
    s.toa.count, s.toa.clear_sky = state.toa.count, state.toa.clear_sky
    empty!(s.sixhour.samples)
    append!(s.sixhour.samples, [_scored_restart_copy_sample(x) for x in state.sixhour_samples])
    s.started = state.started
    return s
end

function Oceananigans.prognostic_state(
    callback::Oceananigans.Callback{P, F, S, CS},
) where {P, F <: ScoredOutputs, S, CS}
    return (schedule = Oceananigans.prognostic_state(callback.schedule),
            func = _scored_restart_state(callback.func))
end

function Oceananigans.restore_prognostic_state!(
    callback::Oceananigans.Callback{P, F, S, CS}, state,
) where {P, F <: ScoredOutputs, S, CS}
    _scored_require(hasproperty(state, :func), "missing scored callback payload")
    _validate_scored_restart_state(callback.func, state.func, state.func.iteration)
    Oceananigans.restore_prognostic_state!(callback.schedule, state.schedule)
    _restore_scored_restart_state!(callback.func, state.func)
    return callback
end

Oceananigans.restore_prognostic_state!(
    ::Oceananigans.Callback{P, F, S, CS}, ::Nothing,
) where {P, F <: ScoredOutputs, S, CS} =
    throw(ArgumentError("scored checkpoint: missing scored callback payload"))

function _validate_scored_restart_compatibility(simulation, state)
    receiver = haskey(simulation.callbacks, SCORED_OUTPUTS_CALLBACK)
    source = hasproperty(state, :callbacks) && !isnothing(state.callbacks) &&
        hasproperty(state.callbacks, SCORED_OUTPUTS_CALLBACK)
    _scored_require(receiver == source, "scored mode differs between checkpoint and receiver")
    receiver || return nothing
    callback = simulation.callbacks[SCORED_OUTPUTS_CALLBACK]
    _scored_require(callback.func isa ScoredOutputs &&
        hasproperty(state.callbacks.readyesm_scored_outputs, :func),
        "missing scored callback payload")
    _scored_require(callback.func.skill ===
        simulation.model.atmosphere.model.callbacks[SCORED_OUTPUTS_CALLBACK],
        "constructed atmosphere skill callback is not shared")
    clock = state.model.clock
    i = clock.iteration
    _validate_scored_restart_state(callback.func, state.callbacks.readyesm_scored_outputs.func, i)
    step_seconds = 3600 / callback.func.steps_per_hour
    _scored_require(clock.time isa Real && abs(Float64(clock.time) - i * step_seconds) <= 1e-6,
        "model clock time differs from scored iteration")
    return nothing
end

function Oceananigans.restore_prognostic_state!(
    simulation::Oceananigans.Simulation{<:_ReadyESMRadiativeSpeedyEarthSystem}, state,
)
    _validate_scored_restart_compatibility(simulation, state)
    return invoke(Oceananigans.restore_prognostic_state!,
        Tuple{Oceananigans.Simulation, Any}, simulation, state)
end

Oceananigans.restore_prognostic_state!(
    ::Oceananigans.Simulation{<:_ReadyESMRadiativeSpeedyEarthSystem}, ::Nothing,
) = nothing
