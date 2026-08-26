const _SEA_ICE_INACTIVE_SURFACE_MASK_SCHEME =
    "full-underlying-grid inactive-cell zeroing after initialization and every sea-ice state update"

"""
    _mask_inactive_sea_ice_surface_cells!(sea_ice)

Set sea-ice concentration, conditional thickness, and snow thickness to zero in
columns whose upper ocean cell is inactive.  Gridded ECCO metadata contain
values over land, and `Oceananigans.set!` interpolates those values into the
surface fields even when the ocean uses an immersed bathymetry.

This uses an explicit surface kernel rather than Oceananigans'
`mask_immersed_field!` specialization for z-reduced fields.  In the pinned
Oceananigans version that specialization does not clear a fully immersed test
column, whereas `inactive_cell(i, j, Nz, grid)` is the same authoritative mask
used by the active-ocean dynamics.  Halos are refreshed after the clearing
kernel so tripolar-fold and periodic stencils cannot retain pre-mask ECCO ice.
"""
@kernel function _mask_inactive_sea_ice_surface_cells_kernel!(
    ice_thickness,
    concentration,
    snow_thickness,
    grid,
)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)
    inactive = Oceananigans.Grids.inactive_cell(i, j, k, grid)
    if inactive
        @inbounds begin
            ice_thickness[i, j, 1] = zero(eltype(ice_thickness))
            concentration[i, j, 1] = zero(eltype(concentration))
            if !isnothing(snow_thickness)
                snow_thickness[i, j, 1] = zero(eltype(snow_thickness))
            end
        end
    end
end

# Oceananigans' generic horizontal immersed-field mask launches over the
# immersed grid. When that grid has an active-cell map, `launch!` deliberately
# skips inactive cells -- exactly the cells the mask is intended to clear.
# ClimaSeaIce calls this method after every sea-ice step for its z-reduced
# prognostic fields. This deliberately narrow compatibility method preserves
# Oceananigans' own masking kernel and location semantics while changing only
# the launch domain to the complete underlying horizontal grid.
function Oceananigans.ImmersedBoundaries.mask_immersed_field_xy!(
    field::Oceananigans.ImmersedBoundaries.OnlyZReducedField,
    grid::Oceananigans.ImmersedBoundaries.AGFBIBG,
    loc,
    value,
    k,
)
    architecture = Oceananigans.Architectures.architecture(field)
    instantiated_location =
        Oceananigans.ImmersedBoundaries.instantiate.(loc)
    return Oceananigans.Utils.launch!(
        architecture,
        grid.underlying_grid,
        :xy,
        Oceananigans.ImmersedBoundaries._mask_immersed_field_xy!,
        field,
        instantiated_location,
        grid,
        value,
        k,
    )
end

function _mask_inactive_sea_ice_surface_cells!(sea_ice)
    model = hasproperty(sea_ice, :model) ? sea_ice.model : sea_ice
    grid = model.grid
    architecture = Oceananigans.Architectures.architecture(grid)
    # `launch!` on an immersed grid with an active-cells map deliberately
    # omits the very cells this kernel must clear.  Launch over the complete
    # underlying horizontal grid while passing the immersed grid as the mask.
    launch_grid = grid isa Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ?
        grid.underlying_grid : grid
    Oceananigans.Architectures.synchronize(architecture)
    Oceananigans.Utils.launch!(
        architecture,
        launch_grid,
        :xy,
        _mask_inactive_sea_ice_surface_cells_kernel!,
        model.ice_thickness,
        model.ice_concentration,
        model.snow_thickness,
        grid,
    )
    Oceananigans.Architectures.synchronize(architecture)
    Oceananigans.BoundaryConditions.fill_halo_regions!((
        model.ice_thickness,
        model.ice_concentration,
        model.snow_thickness,
    ))
    Oceananigans.Architectures.synchronize(architecture)
    return nothing
end

@kernel function _ocean_surface_active_mask_kernel!(active, grid)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)
    inactive = Oceananigans.Grids.inactive_cell(i, j, k, grid)
    @inbounds active[i, j, 1] = ifelse(inactive, zero(eltype(active)), one(eltype(active)))
end

"""Return a host `0/1` mask for dynamically active upper-ocean cells."""
function _ocean_surface_active_mask(grid)
    return ones(Float64, size(grid, 1), size(grid, 2))
end

function _ocean_surface_active_mask(
    grid::Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid,
)
    architecture = Oceananigans.Architectures.architecture(grid)
    active = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid.underlying_grid)
    Oceananigans.Utils.launch!(
        architecture,
        grid.underlying_grid,
        :xy,
        _ocean_surface_active_mask_kernel!,
        active,
        grid,
    )
    Oceananigans.Architectures.synchronize(architecture)
    return reshape(
        Float64.(Array(Oceananigans.interior(active))),
        size(grid, 1),
        size(grid, 2),
    )
end
