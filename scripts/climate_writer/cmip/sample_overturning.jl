# Included by the streamed runner after sample_ocean.jl (Claude, 2026-10-29). Reads only.
# Monthly basin overturning streamfunction (msftmz) and northward heat transport (hfbasin), DESIGN.md phase 2b.
include("overturning.jl")

const _CMIP_OVERTURNING_TARGETS = collect(-80.0:0.5:90.0)       # includes 26.5 N
const _CMIP_ATLANTIC_SEED = (320.0, 30.0)                         # 40 W, 30 N
const _CMIP_PACIFIC_SEED = (200.0, 30.0)                          # 160 W, 30 N
const _CMIP_OVERTURNING_CACHE = Ref{Any}(nothing)

function _nearest_wet(longitude, latitude, wet, target)
    best, index = Inf, nothing
    for I in CartesianIndices(wet)
        wet[I] || continue
        dlon = mod(longitude[I] - target[1] + 180, 360) - 180
        d = (dlon * cosd(target[2]))^2 + (latitude[I] - target[2])^2
        d < best && ((best, index) = (d, Tuple(I)))
    end
    index === nothing && error("no wet cell for basin seed $(target)")
    return index
end

function _cmip_overturning_cache(R, ocean, wet3d)
    O = R.Oceananigans
    cache = _CMIP_OVERTURNING_CACHE[]
    cache !== nothing && cache.grid === ocean.grid && return cache
    grid = ocean.grid
    KFO = O.AbstractOperations.KernelFunctionOperation
    C, F = O.Center, O.Face
    dy_x_field = O.Field(KFO{F, C, C}(O.Operators.Δyᶠᶜᶜ, grid))
    dx_y_field = O.Field(KFO{C, F, C}(O.Operators.Δxᶜᶠᶜ, grid))
    O.compute!(dy_x_field)
    O.compute!(dx_y_field)
    dy_x = Array(O.interior(dy_x_field))[:, :, end]
    dx_y = Array(O.interior(dx_y_field))[:, :, end]
    x, y, _ = O.nodes(ocean.tracers.T)
    lon_nodes, lat_nodes = Array(x), Array(y)
    # Latitude-longitude test grids give separate 1D axes. The tripolar grid gives 2D arrays.
    longitude, latitude = ndims(lon_nodes) == 1 ?
        ([lon_nodes[i] for i in eachindex(lon_nodes), j in eachindex(lat_nodes)],
         [lat_nodes[j] for i in eachindex(lon_nodes), j in eachindex(lat_nodes)]) : (lon_nodes, lat_nodes)
    nx, ny = size(wet3d, 1), size(wet3d, 2)
    size(longitude) == size(latitude) == (nx, ny) || throw(DimensionMismatch("tracer coordinates"))
    surface = wet3d[:, :, end]
    underlying = grid isa O.ImmersedBoundaries.ImmersedBoundaryGrid ? grid.underlying_grid : grid
    tripolar = underlying isa O.Grids.OrthogonalSphericalShellGrid
    basins, basin_mode = try
        (OceanOverturning.build_basins(BitMatrix(surface), latitude;
            atlantic_seed=_nearest_wet(longitude, latitude, surface, _CMIP_ATLANTIC_SEED),
            pacific_seed=_nearest_wet(longitude, latitude, surface, _CMIP_PACIFIC_SEED)),
         "Atlantic-Arctic and Indo-Pacific masks by flood fill north of 34 S")
    catch err
        # Only a test grid without continents may connect the two basins. On the production grid it is an error.
        (err isa ErrorException && occursin("basins connect", err.msg) && !tripolar) || rethrow()
        @warn "CMIP overturning: Atlantic and Indo-Pacific connect on this non-tripolar grid; writing global lines only"
        (OceanOverturning.Basins([falses(nx, ny), falses(nx, ny), BitMatrix(surface)]),
         "global lines only, because the two basins connect on this test grid; the basin lines are empty")
    end
    lines, fold_limit = OceanOverturning.build_lines(basins, latitude, _CMIP_OVERTURNING_TARGETS)
    nb = length(basins.masks)
    latitudes = sort(unique(l.latitude for l in lines))
    length(lines) == nb * length(latitudes) || error("each basin needs the same set of lines")
    order = [findfirst(l -> l.basin == b && l.latitude == phi, lines) for b in 1:nb, phi in latitudes]
    nz = size(wet3d, 3)
    specs = WriterIO.WriterGeometry.SPECS
    dims = WriterIO.WriterGeometry.CUSTOM_DIMS
    basin_dim = ("cmip_basin", Float64.(1:nb), "1")
    interface_dim = ("cmip_ocean_interface", Float64.(1:(nz + 1)), "1")
    latitude_dim = ("cmip_overturning_latitude", Float64.(latitudes), "degrees_north")
    if !haskey(specs, "cmip_msftmz")
        specs["cmip_msftmz"] = (:custom, "kg s-1",
            "basin overturning mass streamfunction by the staircase method, summed from the surface; basins " *
            join(OceanOverturning.BASIN_NAMES, ", ") * "; interface 1 is the bottom of level 1, $(nz + 1) the surface; " *
            basin_mode)
        specs["cmip_hfbasin"] = (:custom, "W",
            "basin northward ocean heat transport, centred face temperature (not the model's WENO flux), reference 0 C")
        specs["cmip_overturning_face_count"] = (:custom, "1",
            "number of faces in each staircase line; zero where the basin does not reach that latitude (mask there)")
        dims["cmip_msftmz"] = [basin_dim, interface_dim, latitude_dim]
        dims["cmip_hfbasin"] = [basin_dim, latitude_dim]
        dims["cmip_overturning_face_count"] = [basin_dim, latitude_dim]
    end
    cache = (; grid, lines, order, dy_x, dx_y, basins, basin_mode, fold_limit, latitudes)
    _CMIP_OVERTURNING_CACHE[] = cache
    return cache
end

"""Monthly msftmz (basin x interface x latitude) and hfbasin (basin x latitude). Arguments, including the face
velocities `u` and `v`, come from the ocean sampler's per-step host arrays so nothing is copied twice."""
function cmip_overturning_samples(R, earth, theta, dz, wet3d, u, v)
    ocean = earth.ocean.model
    cache = _cmip_overturning_cache(R, ocean, wet3d)
    size(u, 1) == size(theta, 1) && size(u, 3) == size(theta, 3) && size(v, 2) >= size(theta, 2) ||
        throw(DimensionMismatch("velocity shapes $(size(u)), $(size(v)) against tracers $(size(theta))"))
    properties = earth.interfaces.ocean_properties
    rho0 = Float64(properties.reference_density)
    rho_cp = rho0 * Float64(properties.heat_capacity)
    volume, heat = OceanOverturning.line_transports(cache.lines, u, v, dz, theta, wet3d, cache.dy_x, cache.dx_y, rho_cp)
    psi = OceanOverturning.streamfunction(volume, rho0)
    nb, nl = size(cache.order)
    nz1 = size(psi, 2)
    msftmz = [psi[cache.order[b, l], k] for b in 1:nb, k in 1:nz1, l in 1:nl]
    hfbasin = [heat[cache.order[b, l]] for b in 1:nb, l in 1:nl]
    faces = [Float64(length(cache.lines[cache.order[b, l]].xfaces) + length(cache.lines[cache.order[b, l]].yfaces))
             for b in 1:nb, l in 1:nl]
    return Dict{String,Any}("cmip_msftmz" => msftmz, "cmip_hfbasin" => hfbasin, "cmip_overturning_face_count" => faces)
end
