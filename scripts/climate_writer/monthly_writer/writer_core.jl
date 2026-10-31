# Bounded streaming output for long coupled runs: calendar and accumulator core (Claude, prototype).
# Pure Julia (Dates only). The I/O layer (writer_io.jl) writes the events that `record!` returns.
#
# Sampling contract (Codex review, 2026-10-20): a sample is the state after a coupled step, at model time t.
# It belongs to the period (start, end] that contains t, so a sample exactly at a period end belongs to the
# ending period, and that period is flushed at its final step. The launch state (iteration 0) is never in a
# period mean. It is written separately as the initial snapshot.
# Quadrature: the right-endpoint rectangle rule with equal weights. Each sample represents (t - Δt, t]. The
# mean over a complete period is (1/N) Σ x(tₙ), with N = (period length) / Δt, so months and years are weighted
# by their actual number of days, including leap days.
module WriterCore

using Dates

export Calendar, State, record!, finalize!, snapshot, restore_state, sample_time, period_keys

const BEFORE = Millisecond(1)

struct Calendar
    origin::DateTime
    step_seconds::Int
    function Calendar(origin::DateTime, step_seconds::Integer)
        step_seconds > 0 && 86_400 % step_seconds == 0 ||
            throw(ArgumentError("the coupled step must divide one day"))
        Time(origin) == Time(0) || throw(ArgumentError("the calendar origin must be at 00:00"))
        return new(origin, Int(step_seconds))
    end
end

sample_time(c::Calendar, iteration::Integer) = c.origin + Second(Int64(iteration) * c.step_seconds)

# (start, end]: the period that contains t is the period that contains t minus one millisecond.
day_key(t::DateTime) = Date(t - BEFORE)
month_key(t::DateTime) = (u = t - BEFORE; (year(u), month(u)))
year_key(t::DateTime) = year(t - BEFORE)
period_keys(t::DateTime) = (day = day_key(t), month = month_key(t), year = year_key(t))

bounds(d::Date) = (DateTime(d), DateTime(d) + Day(1))
bounds(m::Tuple{Int, Int}) = (DateTime(m[1], m[2], 1), DateTime(m[1], m[2], 1) + Month(1))
bounds(y::Int) = (DateTime(y, 1, 1), DateTime(y + 1, 1, 1))

expected_samples(c::Calendar, key) = (b = bounds(key); Int(Dates.value(b[2] - b[1]) ÷ (1000 * c.step_seconds)))

mutable struct Period
    key::Any
    first_iteration::Int
    last_iteration::Int
    count::Int
    sums::Dict{String, Array{Float64}}
end
Period() = Period(nothing, 0, 0, 0, Dict{String, Array{Float64}}())

function reset!(p::Period)
    p.key = nothing
    p.first_iteration = p.last_iteration = p.count = 0
    empty!(p.sums)
    return p
end

function add!(p::Period, key, iteration::Int, values::AbstractDict{String})
    if p.key === nothing
        p.key = key
        p.first_iteration = iteration
        for (name, value) in values
            p.sums[name] = zeros(Float64, size(value))
        end
    end
    p.key == key || error("period $(p.key) still open at a sample of period $key")
    Set(keys(p.sums)) == Set(keys(values)) || error("sampled fields changed within period $key")
    for (name, value) in values
        sum = p.sums[name]
        size(sum) == size(value) || throw(DimensionMismatch("field $name changed shape in period $key"))
        all(isfinite, value) || error("non-finite value in field $name at iteration $iteration")
        sum .+= value
    end
    p.count += 1
    p.last_iteration = iteration
    return p
end

means(p::Period) = Dict(name => sum ./ p.count for (name, sum) in p.sums)

mutable struct State
    calendar::Calendar
    last_iteration::Int
    month::Period
    day::Period
    year::Period
    daily_rows::Vector{Any}
end

"""Writer state at a launch or restart iteration. The next sample must be iteration + 1."""
State(calendar::Calendar; iteration::Integer = 0) =
    State(calendar, Int(iteration), Period(), Period(), Period(), Any[])

function event(s::State, kind::Symbol, p::Period; final::Bool = false)
    first_time = sample_time(s.calendar, p.first_iteration)
    last_time = sample_time(s.calendar, p.last_iteration)
    start, stop = bounds(p.key)
    expected = expected_samples(s.calendar, p.key)
    complete = p.count == expected && first_time == start + Second(s.calendar.step_seconds) && last_time == stop
    return (; kind, key = p.key, start, stop, first_iteration = p.first_iteration,
            last_iteration = p.last_iteration, first_time, last_time, count = p.count, expected, complete,
            means = means(p), daily = kind == :month ? copy(s.daily_rows) : Any[], final)
end

function close_day!(s::State)
    p = s.day
    push!(s.daily_rows, (; date = p.key, first_iteration = p.first_iteration, last_iteration = p.last_iteration,
                          count = p.count, expected = expected_samples(s.calendar, p.key),
                          means = Dict(name => only(v) for (name, v) in means(p))))
    reset!(p)
    return nothing
end

"""
Add the state after coupled step `iteration`. `fields` are the monthly-mean arrays, `scalars` the daily global
values (0-dimensional or one-element arrays), `annual` the annual-mean arrays. Returns the completed periods.
"""
function record!(s::State, iteration::Integer, fields::AbstractDict{String}, scalars::AbstractDict{String},
                 annual::AbstractDict{String})
    iteration == s.last_iteration + 1 ||
        error("sample at iteration $iteration does not follow the writer clock $(s.last_iteration)")
    all(v -> length(v) == 1, values(scalars)) || throw(ArgumentError("daily scalars must have one value"))
    t = sample_time(s.calendar, iteration)
    keys_t = period_keys(t)
    add!(s.month, keys_t.month, Int(iteration), fields)
    add!(s.day, keys_t.day, Int(iteration), scalars)
    add!(s.year, keys_t.year, Int(iteration), annual)
    s.last_iteration = Int(iteration)
    events = Any[]
    t == bounds(s.day.key)[2] && close_day!(s)
    if t == bounds(s.month.key)[2]
        push!(events, event(s, :month, s.month))
        reset!(s.month)
        empty!(s.daily_rows)
    end
    if t == bounds(s.year.key)[2]
        push!(events, event(s, :year, s.year))
        reset!(s.year)
    end
    return events
end

"""Close the open periods at the end of a run (not at a segment boundary). Their files are marked incomplete."""
function finalize!(s::State)
    events = Any[]
    s.day.key === nothing || close_day!(s)
    if s.month.key !== nothing
        push!(events, event(s, :month, s.month; final = true))
        reset!(s.month)
        empty!(s.daily_rows)
    end
    if s.year.key !== nothing
        push!(events, event(s, :year, s.year; final = true))
        reset!(s.year)
    end
    return events
end

period_data(p::Period) = (; key = p.key, first_iteration = p.first_iteration, last_iteration = p.last_iteration,
                          count = p.count, sums = Dict(k => copy(v) for (k, v) in p.sums))
function period_from(d)
    p = Period(d.key, d.first_iteration, d.last_iteration, d.count, Dict{String, Array{Float64}}())
    for (k, v) in d.sums
        p.sums[k] = copy(v)
    end
    return p
end

"""Plain-data copy of the writer state for a checkpoint."""
snapshot(s::State) = (; schema = "monthly_writer_state_v1", origin = s.calendar.origin,
                      step_seconds = s.calendar.step_seconds, last_iteration = s.last_iteration,
                      month = period_data(s.month), day = period_data(s.day), year = period_data(s.year),
                      daily_rows = deepcopy(s.daily_rows))

"""Rebuild the state from a snapshot. The calendar must match and the model must be at the writer clock."""
function restore_state(data, calendar::Calendar, model_iteration::Integer)
    data.schema == "monthly_writer_state_v1" || error("unknown writer state schema $(data.schema)")
    (data.origin, data.step_seconds) == (calendar.origin, calendar.step_seconds) ||
        error("writer calendar differs from the model calendar")
    data.last_iteration == model_iteration ||
        error("writer clock $(data.last_iteration) differs from the restored model iteration $model_iteration")
    return State(calendar, data.last_iteration, period_from(data.month), period_from(data.day),
                 period_from(data.year), deepcopy(data.daily_rows))
end

end
