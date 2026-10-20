@testset "Ordered ocean z-star updates" begin
    O = ReadyESM.Oceananigans
    HFS = O.Models.HydrostaticFreeSurfaceModels
    names = (:σᶜᶜⁿ, :σᶠᶜⁿ, :σᶜᶠⁿ, :σᶠᶠⁿ)
    @test pkgversion(O) == v"0.110.15"

    function make_grid(F, nx, ny)
        O.LatitudeLongitudeGrid(O.CPU(), F;
            size=(nx, ny, 5), halo=(5, 5, 5),
            longitude=(0, 360), latitude=(-80, 80),
            z=O.MutableVerticalDiscretization(collect(range(-F(200), zero(F); length=6))))
    end

    for F in (Float32, Float64), (nx, ny) in ((8, 8), (32, 24)),
        pattern in (:constant, :varying)
        grid = make_grid(F, nx, ny)
        reference = make_grid(F, nx, ny)
        model = O.HydrostaticFreeSurfaceModel(grid; tracers=(:T,))
        coordinate = HFS.ZStarCoordinate()
        ii, jj = -3:nx+4, -3:ny+4
        height = copy(parent(grid.z.ηⁿ))
        ax = axes(grid.z.ηⁿ)
        for j in ax[2], i in ax[1]
            height[i-first(ax[1])+1, j-first(ax[2])+1, 1] =
                pattern == :constant ? F(1) : F(0.5 + 0.2sin(0.3i) + 0.1cos(0.7j))
        end
        @test size(parent(model.free_surface.displacement)) == size(height)
        copyto!(parent(model.free_surface.displacement), height)
        copyto!(parent(reference.z.ηⁿ), height)
        for name in names
            fill!(parent(getproperty(reference.z, name)), F(0.9))
        end
        for j in jj, i in ii
            HFS.update_grid_scaling!(reference.z, i, j, reference)
        end
        old = copy(height)
        for j in jj, i in ii
            old[i-first(ax[1])+1, j-first(ax[2])+1, 1] = -F(0.75)
        end
        for entry in (:ab2, :rk, :reconcile)
            previous = nothing
            for repeat in 1:2
                copyto!(parent(grid.z.ηⁿ), old)
                for name in names
                    fill!(parent(getproperty(grid.z, name)), F(0.9))
                end
                fill!(parent(grid.z.σᶜᶜ⁻), F(-7))
                if entry == :ab2
                    method = which(HFS.ab2_step_grid!,
                        Tuple{typeof(grid), typeof(model), typeof(coordinate), Float64, Float64})
                    HFS.ab2_step_grid!(grid, model, coordinate, 1., 0.1)
                elseif entry == :rk
                    method = which(HFS.rk_substep_grid!,
                        Tuple{typeof(grid), typeof(model), typeof(coordinate), Float64})
                    HFS.rk_substep_grid!(grid, model, coordinate, 1.)
                else
                    method = which(HFS.reconcile_vertical_coordinate!,
                        Tuple{typeof(coordinate), typeof(model), typeof(grid)})
                    HFS.reconcile_vertical_coordinate!(coordinate, model, grid)
                end
                @test method.module === ReadyESM.OrderedOceanGrid
                @test parent(grid.z.ηⁿ) == height
                actual = Dict(name => copy(parent(getproperty(grid.z, name))) for name in names)
                for name in names
                    @test isapprox(actual[name], parent(getproperty(reference.z, name));
                                   rtol=0, atol=8eps(F))
                end
                if entry == :reconcile
                    @test parent(grid.z.σᶜᶜ⁻) == parent(grid.z.σᶜᶜⁿ)
                else
                    @test all(==(F(0.9)), parent(grid.z.σᶜᶜ⁻))
                end
                previous === nothing || @test isequal(actual, previous)
                previous = actual
            end
        end
        # Preserve the package's generic method for other model types.
        generic = which(HFS.rk_substep_grid!,
            Tuple{typeof(grid), NamedTuple{(),Tuple{}}, typeof(coordinate), Float64})
        @test generic.module === HFS
    end
end
