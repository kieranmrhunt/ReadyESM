@testset "Conservative sea-ice transport is area normalized" begin
    source = read(joinpath(
        ReadyESM.PROJECT_ROOT,
        "src",
        "conservative_sea_ice_advection.jl",
    ), String)

    # The prognostic inventory is A*h (m), so a two-dimensional horizontal
    # divergence must be divided by horizontal area. This source assertion
    # covers both the ordinary and FluxFormAdvection overloads.
    @test count(==("Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)"),
                split(source, '\n')) == 2
    @test !occursin("Oceananigans.Operators.Vᶜᶜᶜ", source)

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
end
