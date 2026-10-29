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

    # Spectral prognostics use a packed wrapper. Check its stored coefficients
    # rather than iterating the wrapper, which scalar-indexes a GPU parent.
    spectral = ReadyESM.SpeedyWeather.LowerTriangularArray(
        ComplexF32[1 + 0im, 2 + 0im, 3 + 0im], 2, 2,
    )
    @test ReadyESM._dynamic_guard_storage(spectral) === parent(spectral)
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "atmosphere.prognostic.vorticity", spectral, simulation,
    ))
    bad_spectral = ReadyESM.SpeedyWeather.LowerTriangularArray(
        ComplexF32[1 + 0im, NaN + 0im, 3 + 0im], 2, 2,
    )
    spectral_error = try
        ReadyESM._check_dynamic_prognostic_array!(
            "atmosphere.prognostic.vorticity", bad_spectral, simulation,
        )
        ""
    catch caught
        sprint(showerror, caught)
    end
    @test occursin("first invalid index=2", spectral_error)
    @test occursin("index in packed spectral storage", spectral_error)

    grid = ReadyESM.SpeedyWeather.RingGrids.FullGaussianGrid(2)
    ring_field = ReadyESM.SpeedyWeather.RingGrids.Field(grid)
    fill!(ring_field.data, 250.0)
    @test ReadyESM._dynamic_guard_storage(ring_field) === parent(ring_field)
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "atmosphere.grid.temperature", ring_field, simulation,
    ))
    ring_field.data[1] = 401.0
    @test_throws ErrorException ReadyESM._check_dynamic_prognostic_array!(
        "atmosphere.grid.temperature", ring_field, simulation,
    )

    offset = ReadyESM.Oceananigans.Fields.OffsetArray([0.0, 1.0], -1:0)
    @test ReadyESM._dynamic_guard_storage(offset) === parent(offset)
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "sea_ice.ℵ", offset, simulation,
    ))
    interior = @view [999.0, 0.0, 1.0, 999.0][2:3]
    @test ReadyESM._dynamic_guard_storage(interior) === interior
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "sea_ice.ℵ", interior, simulation,
    ))

    O = ReadyESM.Oceananigans
    ocean_grid = O.RectilinearGrid(O.CPU(), Float64;
        size=(4, 4, 2), x=(0, 4), y=(0, 4), z=(-2, 0),
        halo=(3, 3, 3), topology=(O.Periodic, O.Periodic, O.Bounded))
    ocean_field = O.CenterField(ocean_grid)
    checked_interior = O.interior(ocean_field)
    checked_interior .= 0.3
    @test ReadyESM._dynamic_guard_storage(ocean_field) isa SubArray
    @test size(ReadyESM._dynamic_guard_storage(ocean_field)) == size(checked_interior)
    @test isnothing(ReadyESM._check_dynamic_prognostic_array!(
        "ocean.T", ocean_field, simulation,
    ))
    checked_interior[1, 1, 1] = NaN
    field_error = try
        ReadyESM._check_dynamic_prognostic_array!("ocean.T", ocean_field, simulation)
        ""
    catch caught
        sprint(showerror, caught)
    end
    @test occursin("index in checked field interior", field_error)
    checked_interior[1, 1, 1] = 1.5
    bounded_field_error = try
        ReadyESM._check_dynamic_prognostic_array!(
            "atmosphere.prognostic.land.saturation_water_ice", ocean_field, simulation,
        )
        ""
    catch caught
        sprint(showerror, caught)
    end
    @test occursin("index in checked field interior", bounded_field_error)

    # The callback must descend through Terrarium.StateVariables, whose
    # prognostic fields are neither a NamedTuple at this level nor `.data`.
    terrarium = ReadyESM.Terrarium.StateVariables(
        Float64, (),
        (; soil_moisture = [NaN]),
        (; soil_moisture = [0.0]),
        (;), (;), (;), nothing, nothing,
    )
    nested = (model = (
        clock = (iteration = 3, time = 12.5),
        atmosphere = (variables = (
            prognostic = (land = (terrarium = terrarium,),),
            grid = (;),
        ),),
        ocean = (model = nothing,),
        sea_ice = (model = nothing,),
    ),)
    land_error = try
        ReadyESM.DynamicPrognosticGuardCallback()(nested)
        ""
    catch caught
        sprint(showerror, caught)
    end
    @test occursin(
        "atmosphere.prognostic.land.terrarium.prognostic.soil_moisture",
        land_error,
    )
end
