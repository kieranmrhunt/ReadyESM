# Grid tests of sample_overturning.jl without a model (Claude, 2026-10-29).
# 1. Latitude-longitude test grid without land (the compact CPU smoke grid): 1D node axes, Ny + 1 v faces, and the
#    two basins connect, so only global lines are written.
# 2. Tripolar grid without land: the basins connect, which must stay an error on this grid type.
# 3. The production surface mask from a writer file: the basins separate, with their areas and the fold limit.
# Set READYESM_CMIP_WRITER_FIXTURE to a production writer file to run the real-grid basin check.
using Test, NCDatasets
import Oceananigans
const O = Oceananigans
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
using .WriterIO
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "cmip", "sample_overturning.jl"))
const R = (; Oceananigans = O)

fake_earth(grid) = (; ocean = (; model = (; grid, tracers = (; T = O.CenterField(grid)))),
                    interfaces = (; ocean_properties = (; reference_density = 1025.0, heat_capacity = 3991.87)))

@testset "overturning sampler on test grids" begin
    @testset "latitude-longitude grid without land" begin
        grid = O.LatitudeLongitudeGrid(O.CPU(), Float32; size = (32, 16, 4), halo = (5, 5, 5),
                                       longitude = (0, 360), latitude = (-80, 80), z = (-1000, 0))
        nx, ny, nz = size(grid)
        u, v = O.XFaceField(grid), O.YFaceField(grid)
        @test size(O.interior(v)) == (nx, ny + 1, nz)
        theta = fill(10.0, nx, ny, nz)
        dz = fill(250.0, nx, ny, nz)
        wet3d = trues(nx, ny, nz)
        varr = fill(0.1, nx, ny + 1, nz)
        out = cmip_overturning_samples(R, fake_earth(grid), theta, dz, wet3d, zeros(nx, ny, nz), varr)
        cache = _CMIP_OVERTURNING_CACHE[]
        @test occursin("global lines only", cache.basin_mode)
        @test occursin("global lines only", WriterIO.WriterGeometry.SPECS["cmip_msftmz"][3])
        faces = out["cmip_overturning_face_count"]
        @test all(faces[1:2, :] .== 0)
        inside = cache.latitudes .> -75.0                            # a line needs a row centre to its south
        @test all(faces[3, inside] .== nx) && all(faces[3, .!inside] .== 0)
        psi = out["cmip_msftmz"]
        @test size(psi) == (3, nz + 1, length(cache.latitudes))
        @test all(psi[:, end, :] .== 0)                               # zero at the surface
        dx = Array(O.interior(O.Field(O.AbstractOperations.KernelFunctionOperation{O.Center, O.Face, O.Center}(
            O.Operators.Δxᶜᶠᶜ, grid))))[:, :, end]
        l = findfirst(==(0.0), cache.latitudes)
        jface = only(unique(j for (i, j) in cache.lines[cache.order[3, l]].yfaces))
        @test psi[3, 1, l] ≈ 1025.0 * 0.1 * sum(dx[:, jface]) * 1000 rtol = 1e-6
        @test maximum(cache.latitudes) < 75.0                         # below the top row of centres
    end
    @testset "tripolar grid without land stays an error" begin
        grid = O.TripolarGrid(O.CPU(), Float32; size = (32, 16, 4), z = (-1000, 0), halo = (5, 5, 4))
        nx, ny, nz = size(grid)
        _CMIP_OVERTURNING_CACHE[] = nothing
        @test_throws ErrorException cmip_overturning_samples(R, fake_earth(grid), fill(10.0, nx, ny, nz),
            fill(250.0, nx, ny, nz), trues(nx, ny, nz), zeros(nx, ny, nz), zeros(nx, ny, nz))
    end
    production_fixture = get(ENV, "READYESM_CMIP_WRITER_FIXTURE", "")
    if !isempty(production_fixture)
      @testset "production surface mask" begin
        NCDataset(production_fixture) do ds
            lon = Float64.(ds["ocean_longitude"][:, :])
            lat = Float64.(ds["ocean_latitude"][:, :])
            area = Float64.(ds["ocean_cell_area"][:, :])
            wet = ds["ocean_active_mask"][:, :] .== 1
            atl = _nearest_wet(lon, lat, wet, _CMIP_ATLANTIC_SEED)
            pac = _nearest_wet(lon, lat, wet, _CMIP_PACIFIC_SEED)
            basins = OceanOverturning.build_basins(BitMatrix(wet), lat; atlantic_seed = atl, pacific_seed = pac)
            lines, fold = OceanOverturning.build_lines(basins, lat, _CMIP_OVERTURNING_TARGETS)
            areas = [sum(area[m]) / 1e12 for m in basins.masks]
            println("PRODUCTION_BASINS seeds=", (atl, pac), " areas_Mkm2=", round.(areas; digits = 2),
                    " fold_limit=", fold, " lines=", length(lines))
            @test 60 < areas[1] < 100                                  # Atlantic 34 S to 65 N
            @test 120 < areas[2] < 200                                 # Indian plus Pacific 34 S to 65 N
            @test 340 < areas[3] < 380                                 # all active ocean
            @test !any(basins.masks[1] .& basins.masks[2])
            l = findfirst(x -> x.basin == 1 && x.latitude == 26.5, lines)
            @test l !== nothing && length(lines[l].yfaces) + length(lines[l].xfaces) > 40
        end
      end
    end
end
println("SAMPLE_OVERTURNING_GRID_TESTS_DONE")
