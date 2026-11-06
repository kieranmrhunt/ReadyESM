# ECCO4's 0.5-degree EVEL/NVEL products are geographic vectors collocated at
# longitude/latitude cell centers. NumericalEarth 0.6.0 labels them as faces.
# Keep this correction local to ReadyESM; ECCO2 and other native readers retain
# their own location contracts.
const ECCO_CENTERED_VELOCITY_SCHEMA = "readiesm-ecco4-centered-velocity-v1"
const ECCO_VELOCITY_INPAINTING = "downward-then-finite-neighbor-mean-v1"
const ECCO_VELOCITY_INPAINTING_ITERATIONS = 5

struct ECCO4CenteredVelocity{M}
    metadata::M
end

function _validate_ecco4_velocity_metadata(metadata)
    metadata.dataset isa NumericalEarth.ECCO4Monthly || throw(ArgumentError(
        "centered velocity adapter requires ECCO4Monthly",
    ))
    metadata.name in (:u_velocity, :v_velocity) || throw(ArgumentError(
        "centered velocity adapter accepts only EVEL/NVEL",
    ))
    isnothing(metadata.region) || metadata.region isa NumericalEarth.BoundingBox ||
        throw(ArgumentError("centered ECCO velocity supports global or BoundingBox data"))
    return nothing
end

function _ecco4_velocity_cache_path(metadata)
    _validate_ecco4_velocity_metadata(metadata)
    stem = splitext(metadata.filename)[1]
    # A regional field must never share a cache with the global launch field.
    region = isnothing(metadata.region) ? "" :
        "_" * bytes2hex(sha256(repr(metadata.region)))[1:16]
    return joinpath(metadata.dir, stem * region * "_centered_velocity_v1.jld2")
end

function _ecco4_velocity_cache_identity(metadata; halo = (3, 3, 3))
    _validate_ecco4_velocity_metadata(metadata)
    DW = NumericalEarth.DataWrangling
    path = DW.metadata_path(metadata)
    isfile(path) || error("ECCO velocity source is missing: $path")
    # Check the file, rather than trusting the upstream metadata location trait.
    NCDataset(path, "r") do ds
        longitude, latitude = DW.read_file_coords(metadata)
        ds["longitude"][:] == longitude || error("ECCO velocity longitude is not cell centered")
        ds["latitude"][:] == latitude || error("ECCO velocity latitude is not cell centered")
        size(ds[DW.dataset_variable_name(metadata)]) == (720, 360, 50, 1) ||
            error("unexpected ECCO4 velocity dimensions")
        zf = DW.z_interfaces(metadata.dataset)
        z = ds["Z"][:]
        # Retain NumericalEarth's existing vertical discretization. This repair
        # changes horizontal collocation; it does not reconstruct vertical faces.
        length(z) == 50 && all(isfinite, z) && all(diff(z) .< 0) &&
            all(depth -> first(zf) < depth < 0, z) ||
            error("ECCO velocity depth coordinates are invalid")
    end
    source_hash = open(path) do io
        bytes2hex(sha256(io))
    end
    return Dict(
        "schema" => ECCO_CENTERED_VELOCITY_SCHEMA,
        "source_sha256" => source_hash,
        "source_file" => basename(path),
        "variable" => string(metadata.name),
        "date" => string(metadata.dates),
        "region" => repr(metadata.region),
        "halo" => repr(halo),
        "location" => "Center,Center,Center",
        "number_type" => string(eltype(metadata)),
        "numericalearth_version" => string(pkgversion(NumericalEarth)),
        "oceananigans_version" => string(pkgversion(Oceananigans)),
        "inpainting" => ECCO_VELOCITY_INPAINTING,
        "inpainting_maxiter" => string(ECCO_VELOCITY_INPAINTING_ITERATIONS),
        "mask" => "NumericalEarth ECCO4 default: zero, nonfinite, outside (-1e5,1e5)",
    )
end

function _populate_ecco4_centered_velocity!(field, metadata, data)
    Oceananigans.location(field) == (Oceananigans.Center, Oceananigans.Center, Oceananigans.Center) ||
        throw(ArgumentError("ECCO4 velocity source must be a CenterField"))
    NumericalEarth.DataWrangling.set_metadata_field!(field, data, metadata)
    Oceananigans.BoundaryConditions.fill_halo_regions!(field)
    return field
end

"""
Inpaint a CPU velocity field, preserving valid source values and downward
continuation. Use the upstream five horizontal passes and zero fallback, but
accept zero/canceling donor means (O-083) and retain downward donors (O-084).
Return counts so a finite result cannot hide extensive fallback to zero.
"""
function _inpaint_ecco4_velocity!(field, mask; maxiter = ECCO_VELOCITY_INPAINTING_ITERATIONS)
    Oceananigans.Architectures.architecture(field) isa Oceananigans.CPU ||
        throw(ArgumentError("ECCO velocity cache generation is a CPU operation"))
    maxiter >= 0 || throw(ArgumentError("inpainting iteration limit must be nonnegative"))
    values = Oceananigans.interior(field)
    missing = Array(Oceananigans.interior(mask))
    reference = copy(values)
    originally_missing = copy(missing)
    initial_missing = count(missing)
    values[missing] .= eltype(values)(NaN)
    nx, ny, nz = size(values)
    downward_filled = 0
    @inbounds for j in 1:ny, i in 1:nx, k in nz-1:-1:1
        if missing[i, j, k] && isfinite(values[i, j, k+1])
            values[i, j, k] = values[i, j, k+1]
            missing[i, j, k] = false
            downward_filled += 1
        end
    end
    next_values = copy(values)
    horizontal_filled = 0
    zero_mean_filled = 0
    passes = 0
    for _ in 1:maxiter
        Oceananigans.BoundaryConditions.fill_halo_regions!(field)
        gained = 0
        @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
            missing[i, j, k] || continue
            donors = 0
            total = zero(eltype(values))
            for neighbor in (field[i-1, j, k], field[i+1, j, k],
                             field[i, j-1, k], field[i, j+1, k])
                if isfinite(neighbor)
                    total += neighbor
                    donors += 1
                end
            end
            if donors > 0
                next_values[i, j, k] = total / donors
                gained += 1
                zero_mean_filled += iszero(total)
            end
        end
        gained == 0 && break
        copyto!(values, next_values)
        missing .= .!isfinite.(values)
        horizontal_filled += gained
        passes += 1
    end
    fallback_zeros = count(missing)
    values[missing] .= zero(eltype(values))
    Oceananigans.BoundaryConditions.fill_halo_regions!(field)
    all(isfinite, parent(field)) || error("nonfinite centered ECCO velocity after inpainting")
    valid_max_change = zero(eltype(values))
    @inbounds for n in eachindex(values, reference, originally_missing)
        if !originally_missing[n]
            valid_max_change = max(valid_max_change, abs(values[n] - reference[n]))
        end
    end
    iszero(valid_max_change) || error("ECCO velocity inpainting changed valid source data")
    return Dict(
        "source_valid_cells" => length(values) - initial_missing,
        "source_masked_cells" => initial_missing,
        "source_valid_max_absolute_change" => valid_max_change,
        "downward_filled_cells" => downward_filled,
        "horizontal_filled_cells" => horizontal_filled,
        "zero_mean_filled_cells" => zero_mean_filled,
        "zero_fallback_cells" => fallback_zeros,
        "horizontal_passes" => passes,
        "minimum" => minimum(values),
        "maximum" => maximum(values),
    )
end

function _write_ecco4_velocity_cache(path, field, identity, statistics)
    isfile(path) && error("refusing to replace an existing ECCO velocity cache: $path")
    all(isfinite, parent(field)) || error("refusing to cache nonfinite ECCO velocity")
    temporary, io = mktemp(dirname(path))
    close(io)
    try
        JLD2.jldopen(temporary, "w") do file
            file["identity"] = identity
            file["statistics"] = statistics
            file["data"] = Array(parent(field))
        end
        mv(temporary, path; force = false)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

function _read_ecco4_velocity_cache!(field, path, identity)
    JLD2.jldopen(path, "r") do file
        haskey(file, "identity") && file["identity"] == identity ||
            error("ECCO centered velocity cache identity mismatch: $path")
        data = file["data"]
        size(data) == size(parent(field)) || error("ECCO velocity cache shape mismatch: $path")
        eltype(data) == eltype(field) || error("ECCO velocity cache precision mismatch: $path")
        all(isfinite, data) || error("ECCO velocity cache contains nonfinite values: $path")
        all(value -> abs(value) < 1e5, data) || error("ECCO velocity cache exceeds source bounds: $path")
        copyto!(parent(field), data)
    end
    return field
end

function _ecco4_centered_velocity_field(metadata; require_cache = true, halo = (3, 3, 3))
    identity = _ecco4_velocity_cache_identity(metadata; halo)
    DW = NumericalEarth.DataWrangling
    field = Oceananigans.CenterField(DW.native_grid(metadata, Oceananigans.CPU(); halo))
    cache = _ecco4_velocity_cache_path(metadata)
    if isfile(cache)
        return _read_ecco4_velocity_cache!(field, cache, identity)
    end
    require_cache && error(
        "centered ECCO velocity cache is missing: $cache; run " *
        "scripts/prepare_centered_ecco_velocity.jl in an isolated input directory",
    )
    data = DW.retrieve_data(metadata)
    _populate_ecco4_centered_velocity!(field, metadata, data)
    mask = DW.compute_mask(metadata, field)
    statistics = _inpaint_ecco4_velocity!(field, mask)
    _write_ecco4_velocity_cache(cache, field, identity, statistics)
    return field
end

# The owned wrapper makes the existing Oceananigans vector initializer use the
# corrected source geometry. It still interpolates both components to target
# centers, rotates the vector, and then interpolates to the target's faces.
function Oceananigans.Fields.set!(target::Oceananigans.Field, source::ECCO4CenteredVelocity)
    Oceananigans.Architectures.architecture(target) isa Oceananigans.CPU ||
        throw(ArgumentError("centered ECCO velocity must use the CPU vector initializer"))
    native = _ecco4_centered_velocity_field(source.metadata)
    target.grid.Lz <= native.grid.Lz * 1.01 || throw(ArgumentError(
        "target ocean is deeper than the ECCO velocity source",
    ))
    NumericalEarth.DataWrangling.interpolate_physical!(target, native)
    return target
end

@kernel function _stagger_centered_ocean_velocities!(u, v, grid, centered_u, centered_v)
    i, j, k = @index(Global, NTuple)
    # These are CenterFields, so only the component's normal axis changes
    # location. The pinned upstream ℑxy operators introduce a transverse shift.
    @inbounds begin
        u[i, j, k] = Oceananigans.Operators.ℑxᶠᵃᵃ(i, j, k, grid, centered_u)
        v[i, j, k] = Oceananigans.Operators.ℑyᵃᶠᵃ(i, j, k, grid, centered_v)
    end
end

function _native_centered_vector_buffer(field)
    BC = Oceananigans.BoundaryConditions
    bcs = field.boundary_conditions
    bcs.north isa BC.ZBC || return field
    # Geographic components are scalar under the fold before rotation. Native
    # horizontal components change sign afterwards. Rebind the private storage
    # without changing any live model field or its boundary conditions.
    north = BC.BoundaryCondition(bcs.north.classification, -1)
    signed_bcs = BC.FieldBoundaryConditions(
        bcs.west, bcs.east, bcs.south, north, bcs.bottom, bcs.top, bcs.immersed,
    )
    return Oceananigans.Field(Oceananigans.instantiated_location(field), field.grid,
        field.data, signed_bcs, field.indices)
end

function _set_centered_extrinsic_velocities!(velocities, grid, eastward, northward)
    centered_u, centered_v = Oceananigans.CenterField(grid), Oceananigans.CenterField(grid)
    Oceananigans.set!(centered_u, eastward)
    Oceananigans.set!(centered_v, northward)
    arch = Oceananigans.Architectures.architecture(grid)
    Oceananigans.Utils.launch!(arch, grid, :xyz,
        Oceananigans.Models.HydrostaticFreeSurfaceModels._rotate_velocities!,
        centered_u, centered_v, grid)
    centered_u = _native_centered_vector_buffer(centered_u)
    centered_v = _native_centered_vector_buffer(centered_v)
    Oceananigans.BoundaryConditions.fill_halo_regions!(centered_u)
    Oceananigans.BoundaryConditions.fill_halo_regions!(centered_v)
    Oceananigans.Utils.launch!(arch, grid, :xyz, _stagger_centered_ocean_velocities!,
        velocities.u, velocities.v, grid, centered_u, centered_v)
    return nothing
end
