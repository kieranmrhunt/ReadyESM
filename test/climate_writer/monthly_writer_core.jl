# Synthetic tests for the streaming writer calendar and accumulator.
using Test, Dates, Serialization
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_core.jl"))
using .WriterCore
const W = WriterCore

# Deterministic synthetic model: Float32 fields as a function of the iteration.
field(n) = Float32[sin(0.001f0 * n + i + 3j) for i in 1:3, j in 1:2]
scalar(n) = fill(Float32(cos(0.002 * n)), 1)
annual(n) = Float32[n % 7 + 0.25f0 * k for k in 1:4]
inputs(n) = (Dict("f" => field(n)), Dict("s" => scalar(n)), Dict("a" => annual(n)))

function run(cal, first, last; state = W.State(cal; iteration = first - 1))
    events = Any[]
    for n in first:last
        append!(events, W.record!(state, n, inputs(n)...))
    end
    return state, events
end

# Independent reference: explicit (start, stop] interval tests on the sample times, sequential Float64 sums.
function reference(cal, first, last, start, stop, f)
    total, count = nothing, 0
    for n in first:last
        t = W.sample_time(cal, n)
        start < t <= stop || continue
        x = Float64.(f(n))
        total = total === nothing ? x : total .+ x
        count += 1
    end
    return count == 0 ? (nothing, 0) : (total ./ count, count)
end

@testset "writer core" begin
    @test_throws ArgumentError W.Calendar(DateTime(1993, 1, 1), 7)
    @test_throws ArgumentError W.Calendar(DateTime(1993, 1, 1, 6), 900)

    # Six-hourly samples across a leap February and two year ends.
    cal = W.Calendar(DateTime(1995, 11, 1), 21_600)
    last = Int(Dates.value(DateTime(1997, 3, 1) - DateTime(1995, 11, 1)) ÷ (1000 * 21_600))
    state, events = run(cal, 1, last)
    months = filter(e -> e.kind == :month, events)
    years = filter(e -> e.kind == :year, events)
    @test length(months) == 16 && length(years) == 2
    @test [e.key for e in years] == [1995, 1996]
    feb96 = only(filter(e -> e.key == (1996, 2), months))
    @test feb96.count == feb96.expected == 29 * 4 && feb96.complete
    @test only(filter(e -> e.key == (1997, 2), months)).count == 28 * 4
    for e in months
        ref, count = reference(cal, 1, last, e.start, e.stop, field)
        @test e.count == count == e.expected && e.complete
        @test isequal(e.means["f"], ref)                     # same sequential Float64 sums: exact
        @test length(e.daily) == Dates.daysinmonth(e.start)
        for row in e.daily
            r, c = reference(cal, 1, last, DateTime(row.date), DateTime(row.date) + Day(1), scalar)
            @test row.count == c == 4 && row.means["s"] == only(r)
        end
    end
    y95, y96 = years
    @test y95.complete == false && y95.count == 61 * 4         # launched on 1 November
    @test y96.complete && y96.count == y96.expected == 366 * 4
    @test isequal(y96.means["a"], first(reference(cal, 1, last, DateTime(1996), DateTime(1997), annual)))
    @test state.month.key === nothing && state.year.key == 1997 && state.year.count == 59 * 4

    # Sample exactly at a boundary belongs to the ending period; the next sample opens the new one.
    @test W.period_keys(DateTime(1996, 3, 1)).month == (1996, 2)
    @test W.period_keys(DateTime(1996, 3, 1, 0, 15)).month == (1996, 3)
    @test W.period_keys(DateTime(1997, 1, 1)).year == 1996

    # Fifteen-minute samples across February 1996, with checkpoints at arbitrary iterations.
    cal15 = W.Calendar(DateTime(1996, 1, 25), 900)
    n_end = Int(Dates.value(DateTime(1996, 3, 3, 6) - DateTime(1996, 1, 25)) ÷ 900_000)
    _, uninterrupted = run(cal15, 1, n_end)
    @test only(filter(e -> e.key == (1996, 2), uninterrupted)).count == 29 * 96
    for cut in (1, 95, 96, 97, 677, 2_976, n_end - 1)
        s1, e1 = run(cal15, 1, cut)
        path = tempname()
        open(io -> serialize(io, W.snapshot(s1)), path, "w")
        restored = W.restore_state(open(deserialize, path), cal15, cut)
        _, e2 = run(cal15, cut + 1, n_end; state = restored)
        @test isequal(vcat(e1, e2), uninterrupted)
        @test_throws ErrorException W.restore_state(open(deserialize, path), cal15, cut + 1)
        @test_throws ErrorException W.restore_state(open(deserialize, path), W.Calendar(DateTime(1996, 1, 26), 900), cut)
    end
    s, e = run(cal15, 1, 200)
    final = W.finalize!(s)
    jan = only(filter(x -> x.kind == :month, final))
    @test jan.final && !jan.complete && jan.count == 200 && length(jan.daily) == 3
    @test [r.count for r in jan.daily] == [96, 96, 8] && jan.daily[end].expected == 96

    # Clock and content errors.
    s = W.State(cal15)
    @test_throws ErrorException W.record!(s, 2, inputs(2)...)          # gap
    W.record!(s, 1, inputs(1)...)
    @test_throws ErrorException W.record!(s, 1, inputs(1)...)          # duplicate
    bad = inputs(2)
    s = W.State(cal15); W.record!(s, 1, inputs(1)...)
    @test_throws ErrorException W.record!(s, 2, Dict("g" => field(2)), bad[2], bad[3])   # field set changed
    s = W.State(cal15); W.record!(s, 1, inputs(1)...)
    @test_throws ErrorException W.record!(s, 2, Dict("f" => fill(NaN32, 3, 2)), bad[2], bad[3])
    s = W.State(cal15); W.record!(s, 1, inputs(1)...)
    @test_throws DimensionMismatch W.record!(s, 2, Dict("f" => zeros(Float32, 2, 2)), bad[2], bad[3])
    @test_throws ArgumentError W.record!(W.State(cal15), 1, bad[1], Dict("s" => zeros(Float32, 2)), bad[3])
end
println("WRITER_CORE_TESTS_PASS")
