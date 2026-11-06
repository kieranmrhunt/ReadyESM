@testset "Conservative sea-ice transport is area normalized" begin
    source = read(joinpath(
        ReadyESM.PROJECT_ROOT,
        "src",
        "conservative_sea_ice_advection.jl",
    ), String)

    # The prognostic inventories are A, A*h and A*hs. Their three-dimensional
    # Oceananigans fluxes must lose their face Δz before the horizontal
    # divergence is normalized by horizontal area.
    @test count(line -> occursin(
        "Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)", line,
    ), split(source, '\n')) == 4
    @test !occursin("Oceananigans.Operators.Vᶜᶜᶜ", source)
    @test occursin("_sea_ice_flux_2D_x, _ice_volume_flux_x", source)
    @test occursin("_sea_ice_flux_2D_y, _ice_volume_flux_y", source)
    @test !occursin("ClimaSeaIce.horizontal_div_Uc(", source)

    # A closed horizontal flux has zero global inventory tendency only under
    # area weighting. With changing z-star column heights, V-normalization
    # instead leaves a nonzero Az-weighted volume residual.
    horizontal_flux_divergence = [-3.0, 1.0, 2.0]
    horizontal_area = [2.0, 3.0, 5.0]
    zstar_cell_height = [0.75, 1.2, 1.8]
    area_normalized = horizontal_flux_divergence ./ horizontal_area
    volume_normalized = horizontal_flux_divergence ./
                        (horizontal_area .* zstar_cell_height)
    @test sum(horizontal_area .* area_normalized) ≈ 0.0 atol = eps(Float64)
    @test abs(sum(horizontal_area .* volume_normalized)) > 1e-3

    # Call the production divergence on identical horizontal states with two
    # distinct ocean-cell heights. The 2-D transport tendency must not acquire
    # a depth factor from Oceananigans' three-dimensional face areas.
    OA = ReadyESM.Oceananigans
    function production_divergence(depth)
        grid = OA.RectilinearGrid(
            OA.CPU(),
            Float64;
            size = (12, 10, 1),
            x = (0, 12),
            y = (0, 10),
            z = (-depth, 0),
            halo = (5, 5, 1),
            topology = (OA.Periodic, OA.Periodic, OA.Bounded),
        )
        scheme = OA.Advection.materialize_advection(
            OA.WENO(Float64; order = 5), grid,
        )
        u = OA.XFaceField(grid)
        v = OA.YFaceField(grid)
        inventory = OA.CenterField(grid)
        OA.set!(u, 0.7)
        OA.set!(v, -0.3)
        OA.set!(inventory, (x, y, z) -> 0.4 + 0.03sin(x) - 0.02cos(y))
        OA.fill_halo_regions!((u, v, inventory))
        velocities = (; u, v)
        return [ReadyESM._conservative_div_UV(
            i, j, 1, grid, scheme, velocities, inventory,
        ) for i in 1:12, j in 1:10]
    end
    reference = production_divergence(1.0)
    stretched = production_divergence(1_500.0)
    @test maximum(abs, reference) > 1e-6
    @test stretched ≈ reference rtol = 5e-13 atol = 5e-13
end
