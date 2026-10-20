"""Ordered surface-height and z-star scale updates for Oceananigans 0.110.15.

The upstream fused kernel writes a cell height and reads neighbouring heights
in the same launch. Separate launches make the complete new height available
before any staggered interpolation. Keep the original launch domain, scaling
formula and previous-scale assignments.
"""
module OrderedOceanGrid

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!
using Oceananigans.Models: surface_kernel_parameters
using Oceananigans.ImmersedBoundaries: MutableGridOfSomeKind

const HFS = Oceananigans.Models.HydrostaticFreeSurfaceModels
pkgversion(Oceananigans) == v"0.110.15" ||
    error("Recheck the ordered z-star workaround before changing Oceananigans version")

@kernel function copy_height!(source, grid)
    i, j = @index(Global, NTuple)
    @inbounds grid.z.ηⁿ[i, j, 1] = source[i, j, grid.Nz + 1]
end

@kernel function scale_only!(grid)
    i, j = @index(Global, NTuple)
    HFS.update_grid_scaling!(grid.z, i, j, grid)
end

function ordered_update!(grid, model)
    params = surface_kernel_parameters(grid)
    arch = architecture(grid)
    launch!(arch, grid, params, copy_height!, model.free_surface.displacement, grid)
    launch!(arch, grid, params, scale_only!, grid)
    return nothing
end

function HFS.ab2_step_grid!(grid::MutableGridOfSomeKind,
                          model::HFS.HydrostaticFreeSurfaceModel,
                          ztype::HFS.ZStarCoordinate, Δt, χ)
    parent(grid.z.σᶜᶜ⁻) .= parent(grid.z.σᶜᶜⁿ)
    ordered_update!(grid, model)
    return nothing
end

function HFS.rk_substep_grid!(grid::MutableGridOfSomeKind,
                            model::HFS.HydrostaticFreeSurfaceModel,
                            ztype::HFS.ZStarCoordinate, Δt)
    parent(grid.z.σᶜᶜ⁻) .= parent(grid.z.σᶜᶜⁿ)
    ordered_update!(grid, model)
    return nothing
end

function HFS.reconcile_vertical_coordinate!(::HFS.ZStarCoordinate,
                                           model::HFS.HydrostaticFreeSurfaceModel,
                                           grid::MutableGridOfSomeKind)
    ordered_update!(grid, model)
    parent(grid.z.σᶜᶜ⁻) .= parent(grid.z.σᶜᶜⁿ)
    return nothing
end

end
