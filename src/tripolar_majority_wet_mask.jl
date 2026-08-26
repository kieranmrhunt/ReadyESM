const _TRIPOLAR_MAJORITY_WET_MINIMUM_ABSOLUTE_LATITUDE = 45.0
const _TRIPOLAR_MAJORITY_WET_FRACTION_THRESHOLD = 0.5
const _TRIPOLAR_MAJORITY_WET_MINIMUM_ETOPO_SAMPLES = 16
const _TRIPOLAR_MAJORITY_WET_SCHEME =
    "ETOPO2022_1arcmin_area_majority_poleward_of_45deg_minimum_16_samples"
const _TRIPOLAR_MAJORITY_WET_ETOPO_ELEVATION_SHA256 =
    "094f9ecee776698d68e49416e848aff2600d8f53fb26c5fec3226ce273c1079b"

_tripolar_wrap_to_reference(value, reference) =
    reference + mod(value - reference + 180, 360) - 180

function _tripolar_cell_polygon(grid, longitude, i, j)
    reference = _tripolar_wrap_to_reference(longitude[i, j], 0.0)
    return (
        (_tripolar_wrap_to_reference(Float64(grid.λᶠᶠᵃ[i, j]), reference),
         Float64(grid.φᶠᶠᵃ[i, j])),
        (_tripolar_wrap_to_reference(Float64(grid.λᶠᶠᵃ[i + 1, j]), reference),
         Float64(grid.φᶠᶠᵃ[i + 1, j])),
        (_tripolar_wrap_to_reference(Float64(grid.λᶠᶠᵃ[i + 1, j + 1]), reference),
         Float64(grid.φᶠᶠᵃ[i + 1, j + 1])),
        (_tripolar_wrap_to_reference(Float64(grid.λᶠᶠᵃ[i, j + 1]), reference),
         Float64(grid.φᶠᶠᵃ[i, j + 1])),
    )
end

@inline function _tripolar_point_in_polygon(x, y, polygon)
    inside = false
    previous = 4
    for current in 1:4
        x1, y1 = polygon[current]
        x2, y2 = polygon[previous]
        if (y1 > y) != (y2 > y)
            x_crossing = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            x < x_crossing && (inside = !inside)
        end
        previous = current
    end
    return inside
end

function _tripolar_etopo_fraction_record(
    grid,
    longitude,
    i,
    j,
    file_longitude,
    file_latitude,
    latitude_weights,
    elevation,
)
    polygon = _tripolar_cell_polygon(grid, longitude, i, j)
    longitude_minimum, longitude_maximum = extrema(first.(polygon))
    latitude_minimum, latitude_maximum = extrema(last.(polygon))
    latitude_first = max(
        1,
        searchsortedfirst(file_latitude, latitude_minimum),
    )
    latitude_last = min(
        length(file_latitude),
        searchsortedlast(file_latitude, latitude_maximum),
    )
    sample_count = 0
    wet_sample_count = 0
    total_weight = 0.0
    wet_weight = 0.0
    wet_depth_weighted_sum = 0.0

    for longitude_shift in (-360.0, 0.0, 360.0)
        source_minimum = longitude_minimum - longitude_shift
        source_maximum = longitude_maximum - longitude_shift
        longitude_first = max(
            1,
            searchsortedfirst(file_longitude, source_minimum),
        )
        longitude_last = min(
            length(file_longitude),
            searchsortedlast(file_longitude, source_maximum),
        )
        longitude_first > longitude_last && continue
        for jj in latitude_first:latitude_last
            φ = file_latitude[jj]
            weight = latitude_weights[jj]
            for ii in longitude_first:longitude_last
                λ = file_longitude[ii] + longitude_shift
                _tripolar_point_in_polygon(λ, φ, polygon) || continue
                sample_count += 1
                total_weight += weight
                height = elevation[ii, jj]
                if height < 0
                    wet_sample_count += 1
                    wet_weight += weight
                    wet_depth_weighted_sum -= weight * height
                end
            end
        end
    end

    if sample_count == 0
        return (
            fraction = NaN,
            mean_wet_depth = NaN,
            sample_count = 0,
            wet_sample_count = 0,
        )
    end
    total_weight > 0 || error(
        "zero ETOPO area weight inside tripolar cell ($i,$j)",
    )
    return (
        fraction = wet_weight / total_weight,
        mean_wet_depth = wet_weight > 0 ?
            wet_depth_weighted_sum / wet_weight : 0.0,
        sample_count,
        wet_sample_count,
    )
end

function _compute_tripolar_majority_wet_mask(grid)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, = size(cpu_grid)
    (Nx, Ny) == (360, 180) || throw(ArgumentError(
        "the ETOPO majority-area mask is qualified only for the 360x180 tripolar grid",
    ))
    longitude = Float64.(Array(cpu_grid.λᶜᶜᵃ[1:Nx, 1:Ny]))
    latitude = Float64.(Array(cpu_grid.φᶜᶜᵃ[1:Nx, 1:Ny]))

    metadata = NumericalEarth.DataWrangling.Metadatum(
        :bottom_height;
        dataset = NumericalEarth.DataWrangling.ETOPO.ETOPO2022(),
    )
    etopo_path = NumericalEarth.DataWrangling.metadata_path(metadata)
    etopo_variable = NumericalEarth.DataWrangling.dataset_variable_name(metadata)
    dataset = NCDataset(etopo_path, "r")
    file_longitude = Float64.(Array(dataset["lon"][:]))
    file_latitude = Float64.(Array(dataset["lat"][:]))
    elevation = Array{Float32}(dataset[etopo_variable][:, :])
    close(dataset)
    all(isfinite, elevation) || error("ETOPO source contains non-finite elevations")
    etopo_elevation_sha256 = bytes2hex(SHA.sha256(
        reinterpret(UInt8, vec(elevation)),
    ))
    etopo_elevation_sha256 ==
        _TRIPOLAR_MAJORITY_WET_ETOPO_ELEVATION_SHA256 || error(
            "ETOPO decoded elevation hash does not match the qualified source",
        )
    latitude_weights = cosd.(file_latitude)

    canonical = trues(Nx, Ny)
    canonical[Nx÷2+1:Nx, Ny] .= false
    scan_mask = canonical .&
                (abs.(latitude) .>=
                 _TRIPOLAR_MAJORITY_WET_MINIMUM_ABSOLUTE_LATITUDE)
    scan_indices = findall(scan_mask)
    wet_fraction = fill(Float32(NaN), Nx, Ny)
    mean_wet_depth = fill(Float32(NaN), Nx, Ny)
    sample_count = zeros(Int32, Nx, Ny)
    wet_sample_count = zeros(Int32, Nx, Ny)

    scan_seconds = @elapsed Threads.@threads :dynamic for n in eachindex(scan_indices)
        index = scan_indices[n]
        i, j = Tuple(index)
        record = _tripolar_etopo_fraction_record(
            cpu_grid,
            longitude,
            i,
            j,
            file_longitude,
            file_latitude,
            latitude_weights,
            elevation,
        )
        wet_fraction[index] = Float32(record.fraction)
        mean_wet_depth[index] = Float32(record.mean_wet_depth)
        sample_count[index] = Int32(record.sample_count)
        wet_sample_count[index] = Int32(record.wet_sample_count)
    end

    for i in Nx÷2+1:Nx
        source_i = Nx - i + 1
        wet_fraction[i, Ny] = wet_fraction[source_i, Ny]
        mean_wet_depth[i, Ny] = mean_wet_depth[source_i, Ny]
        sample_count[i, Ny] = sample_count[source_i, Ny]
        wet_sample_count[i, Ny] = wet_sample_count[source_i, Ny]
    end

    reliable_mask = scan_mask .&
                    (sample_count .>=
                     _TRIPOLAR_MAJORITY_WET_MINIMUM_ETOPO_SAMPLES) .&
                    isfinite.(wet_fraction)
    majority_active = falses(Nx, Ny)
    majority_active[reliable_mask] .=
        wet_fraction[reliable_mask] .>=
        _TRIPOLAR_MAJORITY_WET_FRACTION_THRESHOLD
    for i in Nx÷2+1:Nx
        source_i = Nx - i + 1
        reliable_mask[i, Ny] = reliable_mask[source_i, Ny]
        majority_active[i, Ny] = majority_active[source_i, Ny]
    end

    # These are not tunable hotspot overrides. They pin the complete ETOPO
    # algorithm to the independently audited causal cells and prevent a future
    # coordinate/provenance regression from silently changing the experiment.
    isapprox(wet_fraction[257, 176], 0.376554607; atol = 2e-6, rtol = 0) ||
        error("ETOPO majority mask does not reproduce the coastal dead-end fraction")
    isapprox(wet_fraction[108, 180], 1.0; atol = 2e-6, rtol = 0) ||
        error("ETOPO majority mask does not reproduce the fold-cell fraction")
    isapprox(wet_fraction[109, 180], 0.641486845; atol = 2e-6, rtol = 0) ||
        error("ETOPO majority mask does not reproduce the eastern blocker fraction")
    isapprox(wet_fraction[253, 179], 0.606734788; atol = 2e-6, rtol = 0) ||
        error("ETOPO majority mask does not reproduce the northern blocker fraction")

    elevation = nothing
    GC.gc()
    return (;
        longitude,
        latitude,
        wet_fraction,
        mean_wet_depth,
        sample_count,
        wet_sample_count,
        scan_mask,
        reliable_mask,
        majority_active,
        etopo_path,
        etopo_elevation_sha256,
        scan_seconds,
    )
end

function _apply_tripolar_majority_wet_mask!(bottom_height, mask)
    host_bottom = Array(Oceananigans.interior(
        Oceananigans.on_architecture(Oceananigans.CPU(), bottom_height),
        :,
        :,
        1,
    ))
    size(host_bottom) == size(mask.majority_active) || error(
        "majority wet mask and tripolar bathymetry sizes differ",
    )
    current_active = host_bottom .< 0
    wet_to_land = mask.reliable_mask .& current_active .& .!mask.majority_active
    land_to_wet = mask.reliable_mask .& .!current_active .& mask.majority_active
    host_bottom[wet_to_land] .= zero(eltype(host_bottom))
    host_bottom[land_to_wet] .= -max.(
        convert.(eltype(host_bottom), mask.mean_wet_depth[land_to_wet]),
        convert(eltype(host_bottom), 10),
    )
    Oceananigans.set!(bottom_height, host_bottom)
    Oceananigans.fill_halo_regions!(bottom_height)
    Oceananigans.Architectures.synchronize(
        Oceananigans.architecture(bottom_height.grid),
    )
    return (;
        current_active,
        adjusted_bottom = host_bottom,
        wet_to_land,
        land_to_wet,
        changed_cells = count(wet_to_land) + count(land_to_wet),
        wet_to_land_cells = count(wet_to_land),
        land_to_wet_cells = count(land_to_wet),
        reliable_cells = count(mask.reliable_mask),
        undersampled_canonical_cells = count(
            mask.scan_mask .& .!mask.reliable_mask,
        ),
    )
end

function _save_tripolar_majority_wet_provenance(
    output_dir,
    mask,
    adjustment,
    final_bottom,
    final_cleanup,
)
    mkpath(output_dir)
    path = joinpath(output_dir, "tripolar_majority_wet_mask.jld2")
    JLD2.jldsave(
        path;
        scheme = _TRIPOLAR_MAJORITY_WET_SCHEME,
        minimum_absolute_latitude =
            _TRIPOLAR_MAJORITY_WET_MINIMUM_ABSOLUTE_LATITUDE,
        wet_fraction_threshold = _TRIPOLAR_MAJORITY_WET_FRACTION_THRESHOLD,
        minimum_etopo_samples = _TRIPOLAR_MAJORITY_WET_MINIMUM_ETOPO_SAMPLES,
        etopo_path = mask.etopo_path,
        etopo_elevation_sha256 = mask.etopo_elevation_sha256,
        longitude = mask.longitude,
        latitude = mask.latitude,
        wet_fraction = mask.wet_fraction,
        mean_wet_depth = mask.mean_wet_depth,
        sample_count = mask.sample_count,
        wet_sample_count = mask.wet_sample_count,
        scan_mask = mask.scan_mask,
        reliable_mask = mask.reliable_mask,
        majority_active = mask.majority_active,
        pre_adjustment_active = adjustment.current_active,
        wet_to_land = adjustment.wet_to_land,
        land_to_wet = adjustment.land_to_wet,
        adjusted_bottom = adjustment.adjusted_bottom,
        final_bottom,
        final_cleanup,
        scan_seconds = mask.scan_seconds,
    )
    return path
end
