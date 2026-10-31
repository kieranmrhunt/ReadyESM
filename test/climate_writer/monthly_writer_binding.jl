using Test, Dates
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
using .WriterIO
const W = WriterIO.WC

@testset "writer-model checkpoint binding" begin
    mktempdir() do dir
        cal = W.Calendar(DateTime(1993, 1, 1), 3600)
        state = W.State(cal)
        for i in 1:7
            W.record!(state, i, Dict("f" => [Float64(i)]), Dict("s" => [Float64(i)]),
                Dict("a" => [Float64(i)]))
        end
        metadata = Dict("source_commit" => "source-a", "configuration_sha256" => "config-a",
            "writer_sha256" => "writer-a", "weights" => "gaussian")
        checkpoint = joinpath(dir, "model.jld2")
        write(checkpoint, "model checkpoint fixture at iteration 7")
        path = joinpath(dir, "writer.jls")
        @test_throws ErrorException save_bound_state(path, state, checkpoint, metadata, 8)
        @test !ispath(path)
        @test_throws ErrorException save_bound_state(path, state, checkpoint, Dict("source_commit"=>"a"), 7)
        @test save_bound_state(path, state, checkpoint, metadata, 7) == path
        @test isfile(path) && isfile(path * ".sha256")
        @test !ispath(path * ".tmp") && !ispath(path * ".sha256.tmp")
        restored = load_bound_state(path, checkpoint, metadata, cal)
        @test isequal(W.snapshot(restored), W.snapshot(state))
        @test check_model_iteration(restored, 7) === restored
        @test_throws ErrorException check_model_iteration(restored, 8)
        @test_throws ErrorException save_bound_state(path, state, checkpoint, metadata, 7)
        for key in keys(metadata)
            changed = copy(metadata)
            changed[key] *= "changed"
            @test_throws ErrorException load_bound_state(path, checkpoint, changed, cal)
        end
        @test_throws ErrorException load_bound_state(path, checkpoint, metadata,
            W.Calendar(DateTime(1993, 2, 1), 3600))
        @test_throws ErrorException load_bound_state(path, checkpoint, metadata,
            W.Calendar(DateTime(1993, 1, 1), 1800))
        relocated = joinpath(dir, "moved_model.jld2")
        cp(checkpoint, relocated)
        @test load_bound_state(path, relocated, metadata, cal).last_iteration == 7
        write(relocated, "different model, same claimed iteration")
        @test_throws ErrorException load_bound_state(path, relocated, metadata, cal)
        @test_throws ErrorException load_bound_state(path, joinpath(dir, "absent"), metadata, cal)

        old = joinpath(dir, "unbound.jls")
        save_state(old, state)
        @test_throws ErrorException load_bound_state(old, checkpoint, metadata, cal)
        for suffix in (".sha256", ".tmp", ".sha256.tmp")
            stale = joinpath(dir, "stale_$(replace(suffix, '.'=>'_')).jls")
            write(stale * suffix, "partial artifact")
            @test_throws ErrorException save_bound_state(stale, state, checkpoint, metadata, 7)
            @test read(stale * suffix, String) == "partial artifact"
        end
        # Checksum is verified before deserializing modified bytes.
        bytes = read(path)
        bytes[end] ⊻= 0x01
        write(path, bytes)
        @test_throws ErrorException load_bound_state(path, checkpoint, metadata, cal)
    end
end
println("WRITER_BINDING_TESTS_PASS")
