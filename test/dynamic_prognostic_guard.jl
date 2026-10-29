@testset "dynamic per-step prognostic guard" begin
    simulation = (model = (clock = (iteration = 3, time = 12.5),),)

    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "ocean.T",
        [-5.0, 0.0, 40.0],
        simulation,
    ))
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "sea_ice.ℵ",
        [0.0, 0.25, 1.0],
        simulation,
    ))
    @test_throws ErrorException ReadyESM._check_dynamic_prognostic_array!(
        "ocean.T",
        [10.0, 40.1],
        simulation,
    )
    @test_throws ErrorException ReadyESM._check_dynamic_prognostic_array!(
        "sea_ice.ℵ",
        [0.0, 1.01],
        simulation,
    )
    @test_throws ErrorException ReadyESM._check_dynamic_prognostic_array!(
        "atmosphere.prognostic.vorticity",
        ComplexF64[0 + 1im, NaN + 0im],
        simulation,
    )
end
