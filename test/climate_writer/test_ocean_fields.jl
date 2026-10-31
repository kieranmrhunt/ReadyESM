# Unit tests for OceanCMIPFields (Claude, 2026-10-29). Julia stdlib only, no model.
# usage: julia --startup-file=no test_ocean_fields.jl
using Test
include(joinpath(@__DIR__, "..", "..", "scripts", "climate_writer", "cmip", "ocean_fields.jl"))
using .OceanCMIPFields: layer_heat_content, band_heat_content, SPECS, BANDS

@testset "ocean CMIP fields" begin
    # Two columns, four layers (k = 1 bottom). Column 1 is 3000 m deep; column 2 has one land cell at the bottom.
    dzc = [1500.0, 1200.0, 250.0, 50.0]                    # bottom to top thickness of column 1
    nx, ny, nz = 2, 1, 4
    dz = zeros(nx, ny, nz)
    dz[1, 1, :] = dzc
    dz[2, 1, :] = [900.0, 1200.0, 250.0, 50.0]
    wet = trues(nx, ny, nz)
    wet[2, 1, 1] = false
    dz[2, 1, 1] = NaN                                     # land: ignored
    # Depth of each cell's upper face below the surface.
    depth_top = zeros(nx, ny, nz)
    for i in 1:nx
        acc = 0.0
        for k in nz:-1:1
            depth_top[i, 1, k] = acc
            wet[i, 1, k] && (acc += dz[i, 1, k])
        end
    end
    theta = fill(10.0, nx, ny, nz)
    theta[2, 1, 1] = NaN
    rho_cp = 1025.0 * 3991.87
    layer = layer_heat_content(theta, dz, wet, rho_cp)
    @test layer[1, 1, 4] ≈ rho_cp * 10 * 50
    @test layer[2, 1, 1] == 0
    @test sum(layer[1, 1, :]) ≈ rho_cp * 10 * 3000
    bands = band_heat_content(theta, dz, depth_top, wet, rho_cp)
    @test bands["cmip_ocean_heat_content_0_300m"][1, 1] ≈ rho_cp * 10 * 300
    @test bands["cmip_ocean_heat_content_300_700m"][1, 1] ≈ rho_cp * 10 * 400   # 400 m of the 1200 m layer
    @test bands["cmip_ocean_heat_content_700_2000m"][1, 1] ≈ rho_cp * 10 * 1300
    @test bands["cmip_ocean_heat_content_below_2000m"][1, 1] ≈ rho_cp * 10 * 1000
    total = sum(bands["cmip_ocean_heat_content_$b"][2, 1] for (b, _, _) in BANDS)
    @test total ≈ sum(layer[2, 1, :])                     # bands partition the column exactly
    @test bands["cmip_ocean_heat_content_below_2000m"][2, 1] == 0     # wet column is 1500 m deep
    @test bands["cmip_ocean_heat_content_700_2000m"][2, 1] ≈ rho_cp * 10 * 800
    # A vertically varying temperature splits the crossing layer by thickness fraction.
    theta2 = copy(theta)
    theta2[1, 1, 2] = 4.0                                 # the 1200 m layer from 300 to 1500 m
    b2 = band_heat_content(theta2, dz, depth_top, wet, rho_cp)
    @test b2["cmip_ocean_heat_content_0_300m"][1, 1] ≈ rho_cp * 10 * 300
    @test b2["cmip_ocean_heat_content_300_700m"][1, 1] ≈ rho_cp * 4 * 400
    @test b2["cmip_ocean_heat_content_700_2000m"][1, 1] ≈ rho_cp * (4 * 800 + 10 * 500)
    @test_throws ErrorException band_heat_content(fill(NaN, nx, ny, nz), dz, depth_top, wet, rho_cp)
    @test_throws ErrorException layer_heat_content(theta, -dz, wet, rho_cp)
    @test_throws DimensionMismatch band_heat_content(theta[:, :, 1:3], dz, depth_top, wet, rho_cp)
    @test length(SPECS) == 4 + length(BANDS)
    @test all(startswith(k, "cmip_") for k in keys(SPECS))
end
println("OCEAN_CMIP_FIELDS_TESTS_PASS")
