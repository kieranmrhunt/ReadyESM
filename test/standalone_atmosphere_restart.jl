module StandaloneAtmosphereRestartTests
using ReadyESM, Test
const R = ReadyESM
const S = R.SpeedyWeather
const O = R.Oceananigans

function work_arrays(simulation)
    scratch = simulation.variables.scratch.transform_memory
    return (scratch.north, scratch.south, scratch.column.north, scratch.column.south)
end

@testset "Standalone atmosphere restores all transform work arrays" begin
    config = R.ExperimentConfig(;
        truncation=7, nlayers=8, duration_days=1,
        forcing=R.ForcingConfig(radiation=:rrtmgp_clear_sky))
    simulation = R.build_rrtmgp(config)
    arrays = work_arrays(simulation)
    @test all(!isempty, arrays)
    @test simulation.variables.scratch.transform_memory ===
        simulation.model.spectral_transform.scratch_memory
    # Distinct nonzero values expose the two nested arrays skipped by
    # SpeedyWeather's generic copy, even before any model timestep.
    for (index, values) in enumerate(arrays)
        fill!(values, ComplexF32(index, -index))
    end
    expected = map(copy, arrays)
    saved = O.prognostic_state(simulation)
    mktempdir() do directory
        file = joinpath(directory, "atmosphere.jld2")
        R.JLD2.jldsave(file; state=saved)
        for source in (saved, R.JLD2.load(file, "state")), repetition in 1:2
            for (index, values) in enumerate(arrays)
                fill!(values, ComplexF32(100 + repetition + index, index))
                @test values != expected[index]
            end
            @test O.restore_prognostic_state!(simulation, source) === simulation
            for (actual, target, original) in zip(work_arrays(simulation), expected, arrays)
                @test actual === original
                @test isequal(actual, target)
            end
            @test simulation.variables.scratch.transform_memory ===
                simulation.model.spectral_transform.scratch_memory
        end
    end
end
println("STANDALONE_ATMOSPHERE_RESTART_PASS")
end
