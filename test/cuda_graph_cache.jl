# Standalone GPU regression; no production input data required.
using ReadyESM, Test

let SW = ReadyESM.SpeedyWeather, C = ReadyESM.CUDA
    C.functional() || error("cuda_graph_cache.jl requires a working CUDA GPU")
    C.allowscalar(false)
    ST = SW.SpeedyTransforms
    ext = Base.get_extension(ST, :SpeedyTransformsCUDAExt)
    @assert ext !== nothing

    @testset "Fourier graph buffer selection and reuse" begin
        for (trunc, nlayers, Grid) in ((7, 2, SW.RingGrids.OctahedralGaussianGrid),
                                      (31, 27, SW.RingGrids.FullGaussianGrid))
            @testset "T$(trunc)L$(nlayers)" begin
                sg = SW.SpectralGrid(; trunc, nlayers, NF=Float32, Grid,
                                      architecture=SW.GPU())
                transform = ST.SpectralTransform(sg; cuda_graphs=true)
                field = zeros(Float32, sg.grid, nlayers)
                north1, south1 = transform.scratch_memory.north, transform.scratch_memory.south
                north2, south2 = zero(north1), zero(south1)
                reference_north, reference_south = zero(north1), zero(south1)
                fill!(field.data, 1f0)
                ST._fourier_batched!(north1, south1, field, transform)
                C.synchronize()
                cache = ext.get_cache(transform, nlayers)
                @test length(cache.forward_execs) == 1
                @test all(exec -> exec !== nothing, values(cache.forward_execs))
                @test maximum(abs, Array(north1)) > 0

                fill!(field.data, 2f0)
                ST._fourier_batched!(north2, south2, field, transform)
                C.synchronize()
                observed_north, observed_south = Array(north2), Array(south2)
                ext.forward_loop!(cache, reference_north, reference_south, field, transform)
                C.synchronize()
                @test observed_north == Array(reference_north)
                @test observed_south == Array(reference_south)

                ST._fourier_batched!(field, reference_north, reference_south, transform)
                C.synchronize()
                fill!(north2, 0f0)
                fill!(south2, 0f0)
                ST._fourier_batched!(field, north2, south2, transform)
                C.synchronize()
                @test all(iszero, Array(field.data))

                # New view wrappers over the same allocations must reuse their
                # graphs, and replay must keep following changing input values.
                counts = (length(cache.forward_execs), length(cache.inverse_execs))
                for value in 1:4
                    viewed_field = SW.RingGrids.Field(view(field.data, :, :), field.grid)
                    fill!(viewed_field.data, Float32(value))
                    ST._fourier_batched!(north2, south2, viewed_field, transform)
                    C.synchronize()
                    @test maximum(abs, Array(north2)) == value * transform.nlon_max
                    fill!(north2, 0f0)
                    fill!(south2, 0f0)
                    ST._fourier_batched!(viewed_field, north2, south2, transform)
                    C.synchronize()
                    @test all(iszero, Array(viewed_field.data))
                end
                @test (length(cache.forward_execs), length(cache.inverse_execs)) == counts
                ext.clear_fourier_graph_cache!()
            end
        end
    end

    # Leave no caller-owned references to the buffers. A raw CUDA executable
    # alone does not prevent Julia from collecting their owning CuArrays.
    @noinline function capture_temporary_buffers(transform, sg)
        field = zeros(Float32, sg.grid, sg.nlayers)
        north = zero(transform.scratch_memory.north)
        south = zero(transform.scratch_memory.south)
        fill!(field.data, 3f0)
        ST._fourier_batched!(north, south, field, transform)
        C.synchronize()
        return (WeakRef(field.data), WeakRef(north), WeakRef(south))
    end

    function delay_device_cycles(cycles::UInt64)
        started = ReadyESM.CUDA.clock(UInt64)
        while ReadyESM.CUDA.clock(UInt64) - started < cycles
        end
        return nothing
    end

    @testset "Cache clearing completes graphs on another stream" begin
        sg = SW.SpectralGrid(; trunc=7, nlayers=2, NF=Float32, architecture=SW.GPU())
        transform = ST.SpectralTransform(sg; cuda_graphs=true)
        field = zeros(Float32, sg.grid, sg.nlayers)
        north, south = transform.scratch_memory.north, transform.scratch_memory.south
        ST._fourier_batched!(north, south, field, transform)
        C.synchronize()
        ReadyESM.CUDA.@cuda threads=1 delay_device_cycles(UInt64(0))
        C.synchronize()
        other_stream = C.CuStream()
        @test C.isdone(other_stream)
        C.stream!(other_stream) do
            ReadyESM.CUDA.@cuda threads=1 delay_device_cycles(UInt64(2_000_000_000))
            ST._fourier_batched!(north, south, field, transform)
        end
        @test !C.isdone(other_stream)
        ext.clear_fourier_graph_cache!()
        @test C.isdone(other_stream)
        C.synchronize(other_stream)
    end

    @testset "Fourier graphs own captured buffers until cache clearing" begin
        sg = SW.SpectralGrid(; trunc=7, nlayers=2, NF=Float32, architecture=SW.GPU())
        transform = ST.SpectralTransform(sg; cuda_graphs=true)
        weak_buffers = capture_temporary_buffers(transform, sg)
        GC.gc(true)
        @test all(reference -> reference.value !== nothing, weak_buffers)
        ext.clear_fourier_graph_cache!()
        GC.gc(true)
        @test all(reference -> reference.value === nothing, weak_buffers)
    end
end
println("CUDA_GRAPH_BUFFER_REGRESSION_PASS")
