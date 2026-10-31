# Included by the streamed runner after sample.jl (Claude, 2026-10-29). Reads only: no model state is changed.
# Returns (monthly, annual): sea-surface height, mixed-layer thickness and banded heat content every month, and the
# per-layer heat content in the annual files beside the existing 3D temperature and salinity. See DESIGN.md.
include("ocean_fields.jl")
OceanCMIPFields.register_metadata!(WriterIO.WriterGeometry.SPECS)

# The runner binds only ReadyESM, so reach Oceananigans through it.
const _CMIP_OC = ReadyESM.Oceananigans

# Kernel functions for Oceananigans KernelFunctionOperation: the height of each cell's upper face (current z-star
# position) and a static wet mask (1 wet, 0 immersed).
@inline _cmip_face_height(i, j, k, grid) =
    _CMIP_OC.Grids.znode(i, j, k + 1, grid, _CMIP_OC.Center(), _CMIP_OC.Center(), _CMIP_OC.Face())
@inline _cmip_wet(i, j, k, grid) = ifelse(
    _CMIP_OC.Grids.inactive_node(i, j, k, grid, _CMIP_OC.Center(), _CMIP_OC.Center(), _CMIP_OC.Center()),
    zero(eltype(grid)), one(eltype(grid)))

const _CMIP_OCEAN_CACHE = Ref{Any}(nothing)

# Native velocities have their own index dimensions, because v on the tripolar grid may have one more row than the
# tracer grid. They are written as Float32 to halve the annual-file size.
function _register_native_velocity!(usize, vsize)
    specs, dims, storage = WriterIO.WriterGeometry.SPECS, WriterIO.WriterGeometry.CUSTOM_DIMS, WriterIO.WriterGeometry.STORAGE
    haskey(specs, "cmip_native_u") && return nothing
    index(name, n) = (name, Float64.(1:n), "1")
    specs["cmip_native_u"] = (:custom, "m s-1",
        "native-grid x-face velocity (grid x direction, not eastward); rotate with the O-015 observer for uo")
    specs["cmip_native_v"] = (:custom, "m s-1",
        "native-grid y-face velocity (grid y direction, not northward); rotate with the O-015 observer for vo")
    dims["cmip_native_u"] = [index("cmip_xface_i", usize[1]), index("cmip_xface_j", usize[2]), index("cmip_level", usize[3])]
    dims["cmip_native_v"] = [index("cmip_yface_i", vsize[1]), index("cmip_yface_j", vsize[2]), index("cmip_level", vsize[3])]
    storage["cmip_native_u"] = Float32
    storage["cmip_native_v"] = Float32
    return nothing
end

function _cmip_ocean_cache(R, ocean)
    O, NE = R.Oceananigans, R.NumericalEarth
    cache = _CMIP_OCEAN_CACHE[]
    cache !== nothing && cache.grid === ocean.grid && return cache
    grid = ocean.grid
    KFO = O.AbstractOperations.KernelFunctionOperation
    C = O.Center
    formulation = ocean.buoyancy.formulation
    # CMIP mlotst: sigma-t difference of 0.125 kg m-3, as a buoyancy difference with the model's own g and rho_0.
    criterion = formulation.gravitational_acceleration * 0.125 / formulation.equation_of_state.reference_density
    mld = NE.Diagnostics.MixedLayerDepthField(formulation, grid, ocean.tracers; difference_criterion=criterion)
    dz = O.Field(KFO{C, C, C}(O.Operators.Δzᶜᶜᶜ, grid))
    top = O.Field(KFO{C, C, C}(_cmip_face_height, grid))
    wet_field = O.Field(KFO{C, C, C}(_cmip_wet, grid))
    O.compute!(wet_field)
    wet = Array(O.interior(wet_field)) .> 0.5
    any(wet) || error("ocean grid has no wet cells")
    cache = (; grid, mld, dz, top, wet, criterion)
    _CMIP_OCEAN_CACHE[] = cache
    return cache
end

"""Ocean CMIP fields. Pass `temperature` (host interior array) to reuse the runner's existing copy. With
`overturning`, include sample_overturning.jl first."""
function cmip_ocean_samples(R, earth; temperature=nothing, overturning=true)
    O = R.Oceananigans
    ocean = earth.ocean.model
    cache = _cmip_ocean_cache(R, ocean)
    O.compute!(cache.mld)
    O.compute!(cache.dz)
    O.compute!(cache.top)
    eta = Array(O.interior(ocean.free_surface.displacement))[:, :, 1]
    theta = temperature === nothing ? Array(O.interior(ocean.tracers.T)) : temperature
    dz = Array(O.interior(cache.dz))
    face_height = Array(O.interior(cache.top))
    size(theta) == size(dz) == size(face_height) == size(cache.wet) ||
        throw(DimensionMismatch("ocean sampler shapes"))
    size(eta) == size(theta)[1:2] || throw(DimensionMismatch("free-surface shape"))
    depth_top = reshape(eta, size(eta)..., 1) .- face_height
    properties = earth.interfaces.ocean_properties
    rho_cp = Float64(properties.reference_density) * Float64(properties.heat_capacity)
    mld = Array(O.interior(cache.mld))[:, :, 1]
    monthly = Dict{String,Any}("cmip_zos" => Float64.(eta), "cmip_mlotst" => Float64.(mld))
    merge!(monthly, OceanCMIPFields.band_heat_content(theta, dz, depth_top, cache.wet, rho_cp))
    annual = Dict{String,Any}(
        "cmip_ocean_layer_heat_content" => OceanCMIPFields.layer_heat_content(theta, dz, cache.wet, rho_cp),
        "cmip_ocean_wet" => Float64.(cache.wet))
    # Native face velocities for offline uo/vo: centring and rotation are linear and the grid is fixed, so the
    # O-015 observer applied to these annual means gives the annual mean geographic velocity exactly.
    u = Array(O.interior(ocean.velocities.u))
    v = Array(O.interior(ocean.velocities.v))
    _register_native_velocity!(size(u), size(v))
    annual["cmip_native_u"] = Float64.(u)
    annual["cmip_native_v"] = Float64.(v)
    # Overturning and heat transport (sample_overturning.jl), from the same per-step host arrays.
    overturning && merge!(monthly, cmip_overturning_samples(R, earth, theta, dz, cache.wet, u, v))
    return monthly, annual
end
