using Test
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "extended_fields", "fields.jl"))
using .LongRunFields

@testset "Long-run field arithmetic and metadata" begin
    swd = Float32[100 120; 200 240; 340 360]
    swu = Float32[10 12; 20 24; 90 100]
    lwd = Float32[300 310; 100 110; 0 0]
    lwu = Float32[390 400; 280 290; 240 250]
    cswu, clwu = swu .- 5, lwu .+ 10
    inputs = (swd, swu, lwd, lwu, cswu, clwu)
    before = deepcopy(inputs)
    fields = radiation_fields(inputs...)
    @test isequal(inputs, before)
    @test fields["toa_net_downward"] == [10., 10.]
    @test fields["toa_clear_net_downward"] == [5., 5.]
    @test fields["surface_net_radiation_downward"] == [0., 18.]
    @test fields["toa_incoming_shortwave"] == [340., 360.]
    @test fields["surface_shortwave_down"] == [100., 120.]
    @test all(name -> haskey(SPECS, name), keys(fields))
    @test all(name -> SPECS[name][2] == "W m-2", keys(fields))
    @test_throws DimensionMismatch radiation_fields(swd[1:2,:], swu, lwd, lwu, cswu, clwu)
    bad = copy(swd); bad[1,1] = NaN
    @test_throws ErrorException radiation_fields(bad, swu, lwd, lwu, cswu, clwu)
    @test_throws ErrorException radiation_fields((zeros(0,2) for _ in 1:6)...)
    ice = ice_fields(fill(3f0, 2, 2), Float32[0 0.5; 1 0.25])
    @test ice["ice_thickness"] == fill(3.,2,2)
    @test ice["ice_volume_per_area"] == [0 1.5; 3 .75]
    @test_throws DimensionMismatch ice_fields(zeros(2,2), zeros(2,3))
    registry = Dict("existing" => (:surface,"m","existing field"))
    register_metadata!(registry)
    @test length(registry) == length(SPECS) + 1
    @test_throws ErrorException register_metadata!(registry)
end
println("LONG_RUN_FIELDS_UNIT_PASS")
