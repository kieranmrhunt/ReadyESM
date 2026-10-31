# Unit tests for OceanOverturning (Claude, 2026-10-29). Julia stdlib only, no model.
# usage: julia --startup-file=no test_overturning.jl
using Test
include(joinpath(@__DIR__, "..", "..", "scripts", "climate_writer", "cmip", "overturning.jl"))
using .OceanOverturning

@testset "ocean overturning" begin
    nx, ny, nz = 8, 10, 3
    lat_rows = [-45.0 + 10 * (j - 0.5) for j in 1:ny]          # -40 ... 50
    latitude = [lat_rows[j] for i in 1:nx, j in 1:ny]
    wet_s = trues(nx, ny)
    wet = trues(nx, ny, nz)
    dzk = [100.0, 50.0, 10.0]                                  # bottom to top
    dz = [dzk[k] for i in 1:nx, j in 1:ny, k in 1:nz]
    dx_y = fill(1.0e5, nx, ny)
    dy_x = fill(2.0e5, nx, ny)
    theta = fill(10.0, nx, ny, nz)
    rho0, rho_cp = 1025.0, 1025.0 * 3991.87
    basins = Basins([falses(nx, ny), falses(nx, ny), BitMatrix(wet_s)])

    @testset "uniform northward flow, full-width line" begin
        lines, fold = build_lines(basins, latitude, [0.0])
        @test fold == 50.0
        global_line = only(filter(l -> l.basin == 3, lines))
        @test length(global_line.yfaces) == nx && isempty(global_line.xfaces)
        u = zeros(nx, ny, nz); v = fill(0.1, nx, ny, nz)
        volume, heat = line_transports([global_line], u, v, dz, theta, wet, dy_x, dx_y, rho_cp)
        @test volume[1, :] ≈ nx .* 0.1 .* 1.0e5 .* dzk
        @test heat[1] ≈ rho_cp * 10 * sum(volume)
        psi = streamfunction(volume, rho0)
        @test psi[1, nz + 1] == 0                              # zero at the surface
        @test psi[1, 1] ≈ rho0 * sum(volume)                   # net transport at the bottom
        @test psi[1, 3] ≈ rho0 * nx * 0.1 * 1.0e5 * 10         # below the top level only
    end
    @testset "closed overturning cell" begin
        lines, _ = build_lines(basins, latitude, [0.0])
        line = only(filter(l -> l.basin == 3, lines))
        v = zeros(nx, ny, nz)
        v[:, :, 3] .= 0.5                                      # 10 m layer northward
        v[:, :, 1] .= -0.05                                    # 100 m layer southward: same volume
        volume, _ = line_transports([line], zeros(nx, ny, nz), v, dz, theta, wet, dy_x, dx_y, rho_cp)
        psi = streamfunction(volume, rho0)
        @test isapprox(psi[1, 1], 0; atol=1e-6 * maximum(abs, psi))
        @test psi[1, 3] ≈ rho0 * nx * 0.5 * 1.0e5 * 10
        @test argmax(psi[1, :]) in (2, 3)
    end
    @testset "land cells and partial lines" begin
        wet2 = copy(wet)
        wet2[1:2, :, 1] .= false                               # bottom level is land in two columns
        lines, _ = build_lines(basins, latitude, [0.0])
        line = only(filter(l -> l.basin == 3, lines))
        v = fill(0.1, nx, ny, nz)
        volume, _ = line_transports([line], zeros(nx, ny, nz), v, dz, theta, wet2, dy_x, dx_y, rho_cp)
        @test volume[1, 1] ≈ (nx - 2) * 0.1 * 1.0e5 * 100
        @test volume[1, 3] ≈ nx * 0.1 * 1.0e5 * 10
    end
    @testset "x faces when latitude varies along x" begin
        # Latitude increases with i in each row, so a line crosses x faces. Eastward flow then crosses northward.
        lat2 = [-30.0 + 5 * i + 10 * j for i in 1:nx, j in 1:ny]
        b2 = Basins([falses(nx, ny), falses(nx, ny), BitMatrix(wet_s)])
        lines, _ = build_lines(b2, lat2, [42.0])
        line = only(filter(l -> l.basin == 3, lines))
        @test !isempty(line.xfaces) && !isempty(line.yfaces)
        u = fill(0.2, nx, ny, nz)
        volume, _ = line_transports([line], u, zeros(nx, ny, nz), dz, theta, wet, dy_x, dx_y, rho_cp)
        nplus = count(==(1), line.xsign); nminus = count(==(-1), line.xsign)
        @test volume[1, 3] ≈ (nplus - nminus) * 0.2 * 2.0e5 * 10
        @test nplus > 0
        # A uniform flow in a periodic channel has zero divergence, so the net transport through any closed
        # staircase around a south set equals the flux through its boundary: compare two equivalent lines.
        v = fill(0.1, nx, ny, nz)
        vol_uv, _ = line_transports([line], u, v, dz, theta, wet, dy_x, dx_y, rho_cp)
        vol_v, _ = line_transports([line], zeros(nx, ny, nz), v, dz, theta, wet, dy_x, dx_y, rho_cp)
        @test vol_uv ≈ volume .+ vol_v
    end
    @testset "basins by flood fill" begin
        wet_b = trues(nx, ny)
        wet_b[4, 3:ny] .= false                                # land wall between columns 1-3 and 5-7
        wet_b[8, 3:ny] .= false                                # second wall (periodic wrap)
        b = build_basins(wet_b, latitude; atlantic_seed=(2, 8), pacific_seed=(6, 8), south_cut=-25.0)
        atlantic, pacific, glob = b.masks
        @test all(atlantic[1:3, 3:ny]) && !any(atlantic[5:8, :])
        @test all(pacific[5:7, 3:ny]) && !any(pacific[1:4, :])
        @test !any(atlantic[:, 1:2]) && !any(pacific[:, 1:2])  # rows at -40 and -30 are south of the cut
        @test glob == BitMatrix(wet_b)
        leaky = copy(wet_b); leaky[4, 6] = true                # open a strait north of the cut
        @test_throws ErrorException build_basins(leaky, latitude; atlantic_seed=(2, 8), pacific_seed=(6, 8),
                                                 south_cut=-25.0)
        @test_throws ErrorException build_basins(wet_b, latitude; atlantic_seed=(4, 8), pacific_seed=(6, 8))
        # Basin lines use only basin cells.
        lines, _ = build_lines(b, latitude, [5.0])
        atl = only(filter(l -> l.basin == 1, lines))
        @test length(atl.yfaces) == 3
    end
    @testset "fold limit" begin
        lat3 = copy(latitude); lat3[:, ny] .= 70.0
        lines, fold = build_lines(basins, lat3, [0.0, 60.0, 70.0, 80.0])
        @test fold == 70.0
        @test sort(unique(l.latitude for l in lines)) == [0.0, 60.0]
    end
end
println("OCEAN_OVERTURNING_TESTS_PASS")
