using ReadyESM
using Test

function test_tripolar_runoff_delivery(architecture)
    O = ReadyESM.Oceananigans
    @testset "Runoff delivery uses unique tripolar ocean cells" begin
        for fold in (O.RightCenterFolded, O.RightFaceFolded), immersed in (false, true)
            plain_grid = O.TripolarGrid(O.CPU(), Float32;
                size = (16, 12, 2), halo = (3, 3, 3), z = (-100f0, 0f0),
                fold_topology = fold)
            cpu_grid = immersed ? O.ImmersedBoundaryGrid(plain_grid,
                O.GridFittedBottom((longitude, latitude) -> latitude < 0 ? 0f0 : -100f0)) : plain_grid
            wet = ReadyESM._wet_ocean_cells(cpu_grid)
            @test !isempty(wet.wet_i)
            @test all(n -> !(fold == O.RightCenterFolded &&
                             wet.wet_j[n] == 12 && wet.wet_i[n] > 8), eachindex(wet.wet_i))
            if !immersed
                @test length(wet.wet_i) == 16 * 12 - (fold == O.RightCenterFolded ? 8 : 0)
            end
            # Put a source at the fold, where a duplicate storage slot is
            # geographically indistinguishable from its physical partner.
            longitude = [wet.λc[2, 12], wet.λc[6, 12]]
            latitude = [wet.φc[2, 12], wet.φc[6, 12]]
            area = [1e10, 2e10]
            runoff = reshape(Float32[1e-7, 4e-7], 2, 1, 1)
            for weighting in (:equal_area_flux, :equal_column_fraction)
                routing = ReadyESM._build_runoff_routing(longitude, latitude, area,
                    cpu_grid; receivers_per_source = 4, routing_weighting = weighting)
                @test !any(eachindex(routing.target_i)) do target
                    fold == O.RightCenterFolded && routing.target_j[target] == 12 &&
                        routing.target_i[target] > 8
                end
                on_device(value) = O.on_architecture(architecture, value)
                grid = on_device(cpu_grid)
                land = ReadyESM.TerrariumRunoffLand(on_device(copy(runoff)),
                    on_device(routing.contribution_column), on_device(routing.contribution_weight),
                    on_device(routing.target_i), on_device(routing.target_j),
                    on_device(routing.offsets), 2, length(routing.target_i),
                    routing.receivers_per_source, routing.routing_weighting, sum(area), nothing)
                flux = O.Field{O.Center, O.Center, Nothing}(grid)
                ReadyESM._scatter_terrarium_runoff_state!(flux, land)
                O.Architectures.synchronize(architecture)
                source = 1000 * sum(area .* Float64.(vec(runoff)))
                delivered = ReadyESM._budget_scalar(ReadyESM._budget_integral(flux))
                @test delivered ≈ source rtol = 2e-6
                before = Array(O.interior(flux))
                O.fill_halo_regions!(flux)
                O.Architectures.synchronize(architecture)
                @test Array(O.interior(flux)) == before
                @test ReadyESM._budget_scalar(ReadyESM._budget_integral(flux)) ≈ source rtol = 2e-6

                # Exercise the actual dynamic gate, including its copied fold
                # storage. Diagnostics must use the same physical domain for
                # the runoff rate, binary mask, and receiver footprint.
                mixing = (
                    dynamically_gated = true,
                    vertical_diffusivity = O.CenterField(grid),
                    horizontal_diffusivity = O.CenterField(grid),
                    eligible_vertical_diffusivity = O.CenterField(grid),
                    eligible_horizontal_diffusivity = O.CenterField(grid),
                    active_receiver_mask = O.Field{O.Center, O.Center, Nothing}(grid),
                    reference_freshwater_mass_flux_kgm2s = 0f0,
                )
                fill!(mixing.eligible_vertical_diffusivity, 1f-4)
                fill!(mixing.eligible_horizontal_diffusivity, 1f0)
                gated = ReadyESM.TerrariumRunoffLand(land.surface_runoff,
                    land.contribution_column, land.contribution_weight,
                    land.target_i, land.target_j, land.offsets, land.source_columns,
                    land.target_cells, land.receivers_per_source,
                    land.routing_weighting, land.represented_land_area, mixing)
                ReadyESM._scatter_terrarium_runoff_state!(flux, gated)
                O.Architectures.synchronize(architecture)
                live_mask = Array(O.interior(mixing.active_receiver_mask))
                reported_mask = ReadyESM._runoff_surface_grid_matrix(mixing.active_receiver_mask)
                reported_flux = ReadyESM._runoff_surface_grid_matrix(flux)
                @test reported_mask == Float64.(ReadyESM._river_mouth_mixing_active.(reported_flux, 0f0))
                @test count(==(1), reported_mask) == length(routing.target_i)
                @test Array(O.interior(mixing.active_receiver_mask)) == live_mask
                @test ReadyESM._budget_scalar(ReadyESM._budget_integral(flux)) ≈ source rtol = 2e-6
                if fold == O.RightCenterFolded
                    @test all(iszero, reported_mask[9:16, 12])
                    @test reported_mask[1:8, 12] == live_mask[1:8, 12, 1]
                else
                    @test reported_mask == Float64.(live_mask[:, :, 1])
                end
                fill!(gated.surface_runoff, 0)
                ReadyESM._scatter_terrarium_runoff_state!(flux, gated)
                O.Architectures.synchronize(architecture)
                @test all(iszero, Array(O.interior(flux)))
                @test all(iszero, Array(O.interior(mixing.active_receiver_mask)))
                @test all(iszero, Array(O.interior(mixing.vertical_diffusivity)))
                @test all(iszero, Array(O.interior(mixing.horizontal_diffusivity)))
            end
        end
    end
end

test_tripolar_runoff_delivery(ReadyESM.Oceananigans.CPU())
if get(ENV, "READYESM_RUNOFF_DELIVERY_TEST_GPU", "false") == "true"
    @eval import CUDA
    CUDA.functional() || error("runoff-delivery GPU test requires CUDA")
    CUDA.allowscalar(false)
    test_tripolar_runoff_delivery(ReadyESM.Oceananigans.GPU())
end
