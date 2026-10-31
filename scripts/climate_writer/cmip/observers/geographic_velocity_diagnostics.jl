const _SURFACE_GEOGRAPHIC_VELOCITY_PROVENANCE =
    "native_face_to_center_average_then_geographic_rotation_v1"

# Preserve the physical face values and their boundary conventions, but give
# halo filling private data and private boundary-condition caches. A tripolar
# halo fill can also rewrite redundant interior faces at the fold.
function _surface_velocity_buffer(source, k)
    LX, LY, LZ = Oceananigans.location(source)
    indices = (:, :, LZ === Nothing ? Colon() : (k:k))
    source_bcs = source.boundary_conditions
    west, east, south, north = deepcopy((
        source_bcs.west, source_bcs.east, source_bcs.south, source_bcs.north,
    ))
    bcs = Oceananigans.FieldBoundaryConditions(
        source.grid, (LX(), LY(), LZ()), indices;
        west, east, south, north, bottom=nothing, top=nothing, immersed=nothing,
    )
    buffer = Oceananigans.Field{LX, LY, LZ}(
        source.grid; indices, boundary_conditions=bcs,
    )
    parent_k = k - first(axes(source.data, 3)) + 1
    copyto!(parent(buffer), @view parent(source)[:, :, parent_k:parent_k])
    return buffer
end

"""
    _geographic_surface_velocity(u, v; k=...)

Return `(eastward, northward)` host matrices at horizontal cell centers from
native staggered velocity Fields. By default sample the top ocean layer, or
`k=1` for vertically reduced sea-ice velocity Fields. On latitude-longitude
and rectilinear grids the rotation is the identity. On a tripolar grid use
Oceananigans' two-component geographic rotation after collocation.

Only private surface buffers have their halos filled. Neither source data
(including halos/fold storage), boundary caches, nor restart histories change.
This helper can also be used by streaming writers before averaging vectors.
"""
function _geographic_surface_velocity(
    u, v;
    k=Oceananigans.location(u)[3] === Nothing ? 1 : size(u.grid, 3),
)
    u.grid === v.grid || throw(ArgumentError("velocity Fields must share a grid"))
    Oceananigans.location(u)[1:2] == (Oceananigans.Face, Oceananigans.Center) &&
        Oceananigans.location(v)[1:2] == (Oceananigans.Center, Oceananigans.Face) ||
        throw(ArgumentError("expected native x-face u and y-face v Fields"))
    Oceananigans.location(u)[3] === Oceananigans.location(v)[3] ||
        throw(ArgumentError("velocity Fields must share their vertical location"))
    grid = u.grid
    buffers = (_surface_velocity_buffer(u, k), _surface_velocity_buffer(v, k))
    Oceananigans.fill_halo_regions!(buffers)
    eastward = Oceananigans.Field{Oceananigans.Center, Oceananigans.Center, Nothing}(grid)
    northward = Oceananigans.Field{Oceananigans.Center, Oceananigans.Center, Nothing}(grid)
    architecture = Oceananigans.architecture(grid)
    Oceananigans.Utils.launch!(
        architecture, grid, :xy, _geographic_surface_velocity_kernel!,
        eastward, northward, buffers..., grid, k,
    )
    Oceananigans.Architectures.synchronize(architecture)
    return (
        eastward=Float64.(Array(Oceananigans.interior(eastward))[:, :, 1]),
        northward=Float64.(Array(Oceananigans.interior(northward))[:, :, 1]),
    )
end

@kernel function _geographic_surface_velocity_kernel!(eastward, northward, u, v, grid, k)
    i, j = @index(Global, NTuple)
    uc = Oceananigans.Operators.ℑxᶜᵃᵃ(i, j, k, grid, u)
    vc = Oceananigans.Operators.ℑyᵃᶜᵃ(i, j, k, grid, v)
    # Use the two-component overload: these are horizontal velocities.
    ue, vn = Oceananigans.Operators.extrinsic_vector(i, j, k, grid, uc, vc)
    @inbounds eastward[i, j, 1] = ue
    @inbounds northward[i, j, 1] = vn
end

function _write_surface_velocity_diagnostics!(dataset, diagnostics, ocean_dimensions)
    specifications = (
        (:ocean_surface_zonal_velocity, "eastward_sea_water_velocity", "eastward ocean surface velocity"),
        (:ocean_surface_meridional_velocity, "northward_sea_water_velocity", "northward ocean surface velocity"),
        (:sea_ice_zonal_velocity, "eastward_sea_ice_velocity", "eastward sea-ice velocity"),
        (:sea_ice_meridional_velocity, "northward_sea_ice_velocity", "northward sea-ice velocity"),
    )
    for (name, standard_name, long_name) in specifications
        variable = defVar(dataset, String(name), Float64, ocean_dimensions)
        variable.attrib["units"] = "m/s"
        variable.attrib["standard_name"] = standard_name
        variable.attrib["long_name"] = long_name
        variable.attrib["coordinates"] = "ocean_longitude ocean_latitude"
        variable.attrib["grid_location"] = "horizontal_cell_center"
        variable.attrib["vector_basis"] = "geographic_east_north"
        variable.attrib["interpolation"] = _SURFACE_GEOGRAPHIC_VELOCITY_PROVENANCE
        variable[:, :] = getproperty(diagnostics, name)
    end
    return nothing
end

# These maxima retain their native staggered meaning. Keep the old names as
# explicit aliases because existing qualification readers consume them.
function _write_native_velocity_extrema!(dataset, diagnostics)
    for (axis, legacy_component) in (("x", "zonal"), ("y", "meridional"))
        canonical = "ocean_max_abs_native_$(axis)_velocity"
        legacy = "ocean_max_abs_$(legacy_component)_velocity"
        value = getproperty(diagnostics, Symbol(canonical))
        dataset.attrib[canonical] = value
        dataset.attrib[legacy] = value
        dataset.attrib[legacy * "_deprecated_alias_for"] = canonical
        values = getproperty(diagnostics, Symbol(canonical * "_series"))
        for name in (canonical * "_series", legacy * "_series")
            variable = defVar(dataset, name, Float64, ("coupled_budget_time",))
            variable.attrib["units"] = "m/s"
            variable.attrib["long_name"] = "maximum absolute native $(axis)-face velocity component"
            variable.attrib["vector_basis"] = "native_grid_axes_on_staggered_faces"
            if name == legacy * "_series"
                variable.attrib["deprecated_alias_for"] = canonical * "_series"
            end
            variable[:] = values
        end
    end
    return nothing
end
