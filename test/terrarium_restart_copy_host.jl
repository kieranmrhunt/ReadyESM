@testset "Terrarium restart copy from a host view" begin
    OffsetArray = ReadyESM.Oceananigans.Fields.OffsetArray
    backing = reshape(collect(Float64, 1:60), 4, 3, 5)
    backing[2, 1, 2] = -0.0
    backing[3, 2, 3] = nextfloat(0.0)
    source_parent = @view backing[2:3, 1:2, 2:3]
    source = OffsetArray(source_parent, 0:1, -1:0, 5:6)
    destination = OffsetArray(fill(-1.0, 2, 2, 2), axes(source)...)
    ReadyESM._copy_terrarium_restart_array!(destination, source, "host_view")
    @test parent(source) isa SubArray
    @test axes(destination) == axes(source)
    @test reinterpret(UInt64, vec(Array(parent(destination)))) ==
          reinterpret(UInt64, vec(Array(source_parent)))
    @test_throws DimensionMismatch ReadyESM._copy_terrarium_restart_array!(
        OffsetArray(fill(-1.0, 2, 2, 2), 1:2, -1:0, 5:6), source,
        "wrong_axes")
end
