# Standalone GPU regression; production geometry needs no input downloads.
using ReadyESM, SparseArrays, Test

const C = ReadyESM.CUDA
C.functional() || error("csr_regridding_gpu.jl requires a working CUDA GPU")
C.allowscalar(false)

function check_csr_regridding(A, label)
    T = eltype(A)
    m, n = size(A)
    gpu_matrix = C.CUSPARSE.CuSparseMatrixCSR{T}(A)
    destination = C.fill(T(NaN), m)
    @testset "$label" begin
        for signed in (false, true)
            source = T[(signed ? 0 : 1) + cos(j / 31) / 2 for j in 1:n]
            # Independent higher-precision CPU sparse product; its absolute
            # product bounds rounding even for rows with cancellation.
            reference = Float64.(A) * Float64.(source)
            scale = abs.(Float64.(A)) * abs.(Float64.(source))
            gpu_source = C.CuArray(source)
            ReadyESM._deterministic_csr_spmv!(destination, gpu_matrix, gpu_source)
            C.synchronize()
            observed = Array(destination)
            @test all(isfinite, observed)
            @test all(abs.(Float64.(observed) .- reference) .<=
                      16eps(T) .* max.(scale, 1.0))
            for repetition in 2:20
                fill!(destination, T(NaN))
                ReadyESM._deterministic_csr_spmv!(destination, gpu_matrix, gpu_source)
                C.synchronize()
                @test reinterpret(UInt8, Array(destination)) == reinterpret(UInt8, observed)
            end
        end
        @test_throws DimensionMismatch ReadyESM._deterministic_csr_spmv!(
            C.zeros(T, m + 1), gpu_matrix, C.zeros(T, n),
        )
        @test_throws DimensionMismatch ReadyESM._deterministic_csr_spmv!(
            destination, gpu_matrix, C.zeros(T, n + 1),
        )
    end
end

for T in (Float32, Float64)
    # Empty first, interior and final rows; signed weights and uneven row sizes.
    A = sparse([2, 2, 4, 6, 6, 6], [1, 4, 2, 1, 3, 6],
               T[1, -2, 0.25, 3, -1, 0.5], 9, 6)
    check_csr_regridding(A, "Sparse row edge cases $T")
end

let O = ReadyESM.Oceananigans, SW = ReadyESM.SpeedyWeather,
    CR = ReadyESM.ConservativeRegridding
    # Use the same horizontal geometry conversion as the production exchanger.
    # Bathymetry and vertical coordinates do not enter these intersections.
    sg = SW.SpectralGrid(; NF=Float32, trunc=31, nlayers=27,
                        Grid=ReadyESM.RG.FullGaussianGrid, architecture=SW.CPU())
    grid = O.TripolarGrid(O.CPU(), Float32; size=(360, 180, 1),
                         z=(-1f0, 0f0), halo=(5, 5, 1))
    geometry = ReadyESM._float64_regridding_geometry_grid(grid)
    manifold = CR.GOCore.best_manifold(geometry)
    regridder = ReadyESM._float32_regridder(CR.Regridder(manifold, geometry, sg.grid))
    check_csr_regridding(regridder.intersections, "Production geometry forward")
    check_csr_regridding(copy(transpose(regridder.intersections)), "Production geometry reverse")
end
println("CSR_REGRIDDING_GPU_PASS")
