using ReadyESM, Test
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_sampling.jl"))
using .WriterSampling

@testset "writer spatial means" begin
    @test weighted_mean([2.0, NaN, 4.0], [1.0, 0.0, 3.0]) == 3.5
    @test_throws ErrorException weighted_mean([NaN, 4.0], [1.0, 3.0])
    @test_throws ErrorException weighted_mean([1.0], [0.0])
    @test_throws ErrorException weighted_mean([1.0], [-1.0])
    @test_throws ErrorException weighted_mean([1.0], [NaN])
    @test_throws DimensionMismatch weighted_mean([1.0], [1.0, 1.0])
    for trunc in (7, 31, 63, 127)
        sg = ReadyESM.SpeedyWeather.SpectralGrid(trunc=trunc, nlayers=8)
        weights = Float64.(Array(ReadyESM._global_point_weights(sg)))
        _, latitude = ReadyESM.SpeedyWeather.RingGrids.get_londlatds(sg.grid)
        field = sind.(Float64.(latitude)).^2
        corrected = weighted_mean(field, weights)
        legacy = weighted_mean(field, cosd.(Float64.(latitude)))
        @test isapprox(sum(weights), 1; atol=2e-7, rtol=0)
        @test weighted_mean(fill(7.0, length(weights)), weights) ≈ 7
        @test isapprox(corrected, 1/3; atol=2e-7, rtol=0)
        @test abs(legacy - 1/3) > 0.01
        println("WRITER_AREA_CHECK T$trunc correct=$corrected legacy=$legacy exact=$(1/3)")
    end
end
println("WRITER_SAMPLING_TESTS_PASS")

@testset "Unique tripolar surface reductions" begin
    active=ones(4,3)
    original=copy(active)
    mask=unique_surface_mask(active;right_center_folded=true)
    @test active==original
    @test sum(mask)==10
    @test all(mask[3:4,3].==0)
    @test all(mask[:,1:2].==1)
    @test unique_surface_mask(active)==active
    # Last-row storage mirrors two physical cells, not four independent cells.
    values=ones(4,3); values[:,3]=[2,4,4,2]
    @test weighted_mean(values,mask)==14/10
    @test weighted_mean(values,active)==20/12
    values[3:4,3].=NaN
    @test weighted_mean(values,mask)==14/10
    active[1,1]=0
    mask=unique_surface_mask(active;right_center_folded=true)
    @test sum(mask)==9
    @test weighted_mean(values,mask)==13/9
    @test_throws ErrorException unique_surface_mask(fill(2,4,3))
    @test_throws ErrorException unique_surface_mask(ones(3,3);right_center_folded=true)
end
println("WRITER_UNIQUE_SURFACE_TESTS_PASS")
