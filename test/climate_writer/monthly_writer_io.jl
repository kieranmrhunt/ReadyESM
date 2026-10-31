# Tests for the streaming writer's NetCDF and checkpoint I/O.
using Test, Dates, NCDatasets
const NC = NCDatasets
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
using .WriterIO
const W = WriterIO.WC

@testset "writer I/O" begin
    dir = mktempdir()
    cal = W.Calendar(DateTime(1996, 2, 27), 21_600)
    s = W.State(cal)
    events = Any[]
    for n in 1:20
        append!(events, W.record!(s, n, Dict("f" => Float32[n 2n; 3n 4n]), Dict("s" => [Float32(n)]),
                                  Dict("a" => Float32[n, -n, 0.5n])))
    end
    feb = only(events)
    @test feb.key == (1996, 2) && feb.count == 12 && !feb.complete          # launched on 27 February
    meta = Dict("source_commit" => "abc", "config_sha256" => "def")
    path = write_event(NC, dir, feb, cal, meta)
    @test basename(path) == "monthly_1996_02.nc"
    NC.NCDataset(path) do ds
        @test ds["f"][:, :] == feb.means["f"]
        @test ds.attrib["sample_count"] == 12 && ds.attrib["expected_sample_count"] == 29 * 4
        @test ds.attrib["complete"] == 0 && ds.attrib["final_partial_flush"] == 0
        @test ds.attrib["period_end"] == "1996-03-01T00:00:00" && ds.attrib["last_sample_time"] == "1996-03-01T00:00:00"
        @test ds["day_of_month"][:] == [27, 28, 29] && ds["day_count"][:] == [4, 4, 4]
        @test ds["daily_s"][:] == [2.5, 6.5, 10.5]
        @test ds.attrib["source_commit"] == "abc"
    end
    @test_throws ErrorException write_event(NC, dir, feb, cal, meta)                 # never overwrite
    touch(joinpath(dir, "annual_1996.nc.tmp"))
    final = W.finalize!(s)
    @test_throws ErrorException write_event(NC, dir, only(filter(e -> e.kind == :year, final)), cal, meta)
    init = write_initial(NC, dir, cal, Dict("f" => Float32[1 2; 3 4]), Dict("s" => [5f0]), meta)
    NC.NCDataset(init) do ds
        @test ds["f"][:, :] == [1 2; 3 4] && ds["scalar_s"][] == 5 && ds.attrib["sample_iteration"] == 0
    end
    s2 = W.State(cal)
    for n in 1:7
        W.record!(s2, n, Dict("f" => Float32[n 2n; 3n 4n]), Dict("s" => [Float32(n)]), Dict("a" => Float32[n, -n, 0.5n]))
    end
    state_path = joinpath(dir, "writer_state.jls")
    save_state(state_path, s2)
    @test_throws ErrorException save_state(state_path, s2)
    r = load_state(state_path, cal, 7)
    @test r.last_iteration == 7 && isequal(W.snapshot(r), W.snapshot(s2))
    @test_throws ErrorException load_state(state_path, cal, 8)
    bytes = read(state_path); bytes[end] ⊻= 0x01; write(state_path, bytes)
    @test_throws ErrorException load_state(state_path, cal, 7)
end
println("WRITER_IO_TESTS_PASS")
