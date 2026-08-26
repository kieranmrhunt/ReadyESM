const _SEA_ICE_MECHANICS_MAP_NAMES = (
    :area_equivalent_thickness,
    :conditional_thickness,
    :concentration,
    :volume_tendency,
    :concentration_tendency,
    :thermodynamic_volume_tendency,
    :cumulative_dynamic_volume_change,
    :cumulative_thermodynamic_volume_change,
    :velocity_divergence,
    :pressure,
    :relaxation_parameter,
    :strain_invariant,
    :bulk_viscosity,
    :sigma11,
    :sigma22,
    :sigma12,
    :ice_speed,
    :solver_velocity_change,
)

@inline function _sea_ice_metric_divergence(i, j, k, grid, u, v)
    k_top = size(grid, 3)
    return ClimaSeaIce.Rheologies.ϵ̇D(i, j, k_top, grid, u, v)
end

@inline function _sea_ice_center_speed(i, j, k, grid, u, v)
    u_center = Oceananigans.Operators.ℑxᶜᵃᵃ(i, j, k, grid, u)
    v_center = Oceananigans.Operators.ℑyᵃᶜᵃ(i, j, k, grid, v)
    return hypot(u_center, v_center)
end

@inline function _sea_ice_center_velocity_change(
    i,
    j,
    k,
    grid,
    u,
    v,
    u_previous,
    v_previous,
)
    du = Oceananigans.Operators.ℑxᶜᵃᵃ(i, j, k, grid, u) -
         Oceananigans.Operators.ℑxᶜᵃᵃ(i, j, k, grid, u_previous)
    dv = Oceananigans.Operators.ℑyᵃᶜᵃ(i, j, k, grid, v) -
         Oceananigans.Operators.ℑyᵃᶜᵃ(i, j, k, grid, v_previous)
    return hypot(du, dv)
end

@inline _sea_ice_center_sigma12(i, j, k, grid, sigma12) =
    Oceananigans.Operators.ℑxyᶜᶜᵃ(i, j, k, grid, sigma12)

@inline _sea_ice_cell_area(i, j, k, grid) =
    Oceananigans.Operators.Azᶜᶜᶠ(i, j, size(grid, 3), grid)

@inline function _sea_ice_thermodynamic_volume_tendency(
    i,
    j,
    k,
    grid,
    thermodynamic_ice_mass_flux,
    sea_ice_density,
)
    @inbounds mass_flux = thermodynamic_ice_mass_flux[i, j, 1]
    @inbounds density = sea_ice_density[i, j, 1]
    return mass_flux / density
end

function _sea_ice_mechanics_operation(function_, grid, arguments)
    operation = Oceananigans.KernelFunctionOperation{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(function_, grid, arguments)
    return Oceananigans.Field(operation)
end

function _sea_ice_mechanics_matrix(field)
    values = Array(Oceananigans.interior(field))
    ndims(values) == 3 && (values = values[:, :, end])
    return Float32.(values)
end

"""
Daily, non-invasive maps for diagnosing sea-ice transport and EVP response.

Velocity and corner-stress stencils are evaluated on private buffers whose
halos are filled immediately before sampling. This prevents the observer from
changing live halos and thereby changing the model trajectory it is measuring.
"""
Base.@kwdef mutable struct SeaIceMechanicsDiagnosticsCallback{S, B, O, M, C, A} <:
                           Function
    sample_every_n_steps::Int
    momentum_substeps::Int = 0
    pressure_formulation::String = "replacement_pressure"
    immersed_boundary_drag_coefficient::Float64 = 0.0
    ocean_tripolar_wet_mask::String = "interpolated_mean_elevation"
    source_fields::S
    buffers::B
    operations::O
    active_surface_mask::M
    coordinates::C
    cumulative_dynamic_volume_change::A
    cumulative_thermodynamic_volume_change::A
    last_time::Float64 = NaN
    time_days::Vector{Float64} = Float64[]
    peak_i::Vector{Int64} = Int64[]
    peak_j::Vector{Int64} = Int64[]
    peak_area_equivalent_thickness::Vector{Float64} = Float64[]
    maps::Dict{Symbol, Vector{Matrix{Float32}}} = Dict(
        name => Matrix{Float32}[] for name in _SEA_ICE_MECHANICS_MAP_NAMES
    )
end

function SeaIceMechanicsDiagnosticsCallback(
    earth,
    Δt;
    sample_every_n_steps = max(1, round(Int, 86_400 / Float64(Δt))),
    pressure_formulation = "replacement_pressure",
    immersed_boundary_drag_coefficient = 0.0,
    ocean_tripolar_wet_mask = "interpolated_mean_elevation",
)
    sample_every_n_steps > 0 ||
        throw(ArgumentError("sample_every_n_steps must be positive"))
    sea_ice = earth.sea_ice.model
    ocean = earth.ocean.model
    grid = sea_ice.grid
    auxiliaries = sea_ice.dynamics.auxiliaries.fields
    source_fields = (
        h = sea_ice.ice_thickness,
        concentration = sea_ice.ice_concentration,
        volume_tendency = sea_ice.timestepper.Gⁿ.h,
        concentration_tendency = sea_ice.timestepper.Gⁿ.ℵ,
        thermodynamic_ice_mass_flux =
            sea_ice.mass_fluxes.thermodynamics.ice,
        sea_ice_density = sea_ice.sea_ice_density,
        pressure = auxiliaries.P,
        relaxation_parameter = auxiliaries.α,
        strain_invariant = auxiliaries.Δ,
        bulk_viscosity = auxiliaries.ζᶜᶜᶜ,
        sigma11 = auxiliaries.σ₁₁,
        sigma22 = auxiliaries.σ₂₂,
        u = sea_ice.velocities.u,
        v = sea_ice.velocities.v,
        u_previous = auxiliaries.uⁿ,
        v_previous = auxiliaries.vⁿ,
        sigma12 = auxiliaries.σ₁₂,
    )
    buffers = (
        u = similar(source_fields.u),
        v = similar(source_fields.v),
        u_previous = similar(source_fields.u_previous),
        v_previous = similar(source_fields.v_previous),
        sigma12 = similar(source_fields.sigma12),
    )
    operations = (
        velocity_divergence = _sea_ice_mechanics_operation(
            _sea_ice_metric_divergence,
            grid,
            (buffers.u, buffers.v),
        ),
        sigma12 = _sea_ice_mechanics_operation(
            _sea_ice_center_sigma12,
            grid,
            (buffers.sigma12,),
        ),
        ice_speed = _sea_ice_mechanics_operation(
            _sea_ice_center_speed,
            grid,
            (buffers.u, buffers.v),
        ),
        solver_velocity_change = _sea_ice_mechanics_operation(
            _sea_ice_center_velocity_change,
            grid,
            (
                buffers.u,
                buffers.v,
                buffers.u_previous,
                buffers.v_previous,
            ),
        ),
        thermodynamic_volume_tendency = _sea_ice_mechanics_operation(
            _sea_ice_thermodynamic_volume_tendency,
            grid,
            (
                source_fields.thermodynamic_ice_mass_flux,
                source_fields.sea_ice_density,
            ),
        ),
    )
    active_surface_mask = Float32.(_ocean_surface_active_mask(ocean.grid))
    longitude, latitude, _ = Oceananigans.nodes(source_fields.h)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    Nx, Ny = size(active_surface_mask)
    if ndims(longitude) == 1 && ndims(latitude) == 1
        longitude = repeat(reshape(longitude, :, 1), 1, Ny)
        latitude = repeat(reshape(latitude, 1, :), Nx, 1)
    end
    size(longitude) == (Nx, Ny) || throw(DimensionMismatch(
        "sea-ice longitude coordinates do not match the active surface mask",
    ))
    size(latitude) == (Nx, Ny) || throw(DimensionMismatch(
        "sea-ice latitude coordinates do not match the active surface mask",
    ))
    cell_area_operation = _sea_ice_mechanics_operation(
        _sea_ice_cell_area,
        grid,
        (),
    )
    Oceananigans.compute!(cell_area_operation)
    architecture = Oceananigans.Architectures.architecture(source_fields.h)
    Oceananigans.Architectures.synchronize(architecture)
    cell_area = Float64.(_sea_ice_mechanics_matrix(cell_area_operation))
    size(cell_area) == (Nx, Ny) || throw(DimensionMismatch(
        "sea-ice cell areas do not match the active surface mask",
    ))
    all(isfinite, cell_area) || error("sea-ice cell areas must be finite")
    active = active_surface_mask .== 1
    all(>(0), cell_area[active]) ||
        error("active sea-ice surface cells must have positive area")
    coordinates = (; longitude, latitude, cell_area)
    cumulative_dynamic_volume_change = similar(source_fields.h)
    cumulative_thermodynamic_volume_change = similar(source_fields.h)
    fill!(parent(cumulative_dynamic_volume_change), 0)
    fill!(parent(cumulative_thermodynamic_volume_change), 0)
    return SeaIceMechanicsDiagnosticsCallback(
        ; sample_every_n_steps,
        momentum_substeps = sea_ice.dynamics.solver.substeps,
        pressure_formulation = String(pressure_formulation),
        immersed_boundary_drag_coefficient = Float64(
            immersed_boundary_drag_coefficient,
        ),
        ocean_tripolar_wet_mask = String(ocean_tripolar_wet_mask),
        source_fields,
        buffers,
        operations,
        active_surface_mask,
        coordinates,
        cumulative_dynamic_volume_change,
        cumulative_thermodynamic_volume_change,
    )
end

function _accumulate_sea_ice_volume_changes!(
    cumulative_dynamic,
    cumulative_thermodynamic,
    dynamic_tendency,
    thermodynamic_tendency,
    active_surface_mask,
    Δt,
)
    size(cumulative_dynamic) == size(cumulative_thermodynamic) ==
        size(dynamic_tendency) == size(thermodynamic_tendency) ==
        size(active_surface_mask) || throw(DimensionMismatch(
        "sea-ice mechanics accumulation maps differ in size",
    ))
    Δt >= 0 || error("sea-ice mechanics accumulation time moved backwards")
    for index in eachindex(cumulative_dynamic)
        if active_surface_mask[index] == 1
            cumulative_dynamic[index] += Δt * dynamic_tendency[index]
            cumulative_thermodynamic[index] += Δt * thermodynamic_tendency[index]
        end
    end
end

function _accumulate_sea_ice_mechanics!(callback, Δt)
    Oceananigans.compute!(callback.operations.thermodynamic_volume_tendency)
    grid = callback.source_fields.h.grid
    architecture = Oceananigans.Architectures.architecture(grid)
    Oceananigans.Utils.launch!(
        architecture,
        grid,
        :xy,
        _accumulate_sea_ice_volume_changes_kernel!,
        callback.cumulative_dynamic_volume_change,
        callback.cumulative_thermodynamic_volume_change,
        callback.source_fields.volume_tendency,
        callback.operations.thermodynamic_volume_tendency,
        grid,
        Δt,
    )
    return nothing
end

@kernel function _accumulate_sea_ice_volume_changes_kernel!(
    cumulative_dynamic,
    cumulative_thermodynamic,
    dynamic_tendency,
    thermodynamic_tendency,
    grid,
    Δt,
)
    i, j = @index(Global, NTuple)
    k = 1
    k_top = size(grid, 3)
    if Oceananigans.Grids.inactive_cell(i, j, k_top, grid)
        @inbounds cumulative_dynamic[i, j, k] = 0
        @inbounds cumulative_thermodynamic[i, j, k] = 0
    else
        @inbounds cumulative_dynamic[i, j, k] +=
            Δt * dynamic_tendency[i, j, k]
        @inbounds cumulative_thermodynamic[i, j, k] +=
            Δt * thermodynamic_tendency[i, j, k]
    end
end

function _prepare_sea_ice_mechanics_buffers!(callback)
    source = callback.source_fields
    buffers = callback.buffers
    for name in (:u, :v, :u_previous, :v_previous, :sigma12)
        destination = getproperty(buffers, name)
        origin = getproperty(source, name)
        copyto!(parent(destination), parent(origin))
    end
    Oceananigans.fill_halo_regions!(
        (
            buffers.u,
            buffers.v,
            buffers.u_previous,
            buffers.v_previous,
            buffers.sigma12,
        ),
    )
    architecture = Oceananigans.Architectures.architecture(source.h)
    Oceananigans.Architectures.synchronize(architecture)
    for operation in values(callback.operations)
        Oceananigans.compute!(operation)
    end
    Oceananigans.Architectures.synchronize(architecture)
    return nothing
end

function _sample_sea_ice_mechanics!(callback, simulation)
    _prepare_sea_ice_mechanics_buffers!(callback)
    source = callback.source_fields
    h = _sea_ice_mechanics_matrix(source.h)
    concentration = _sea_ice_mechanics_matrix(source.concentration)
    area_equivalent = h .* concentration
    active = callback.active_surface_mask .== 1
    size(area_equivalent) == size(active) || throw(DimensionMismatch(
        "sea-ice mechanics map and active surface mask differ",
    ))
    any(active) || error("sea-ice mechanics observer found no active surface cells")
    peak_values = ifelse.(active, area_equivalent, Float32(-Inf))
    peak = argmax(peak_values)

    push!(callback.time_days, Float64(simulation.model.clock.time) / 86_400)
    push!(callback.peak_i, Int64(peak[1]))
    push!(callback.peak_j, Int64(peak[2]))
    push!(
        callback.peak_area_equivalent_thickness,
        Float64(area_equivalent[peak]),
    )

    sample = (
        area_equivalent_thickness = area_equivalent,
        conditional_thickness = h,
        concentration,
        volume_tendency = _sea_ice_mechanics_matrix(source.volume_tendency),
        concentration_tendency = _sea_ice_mechanics_matrix(
            source.concentration_tendency,
        ),
        thermodynamic_volume_tendency = _sea_ice_mechanics_matrix(
            callback.operations.thermodynamic_volume_tendency,
        ),
        cumulative_dynamic_volume_change = Float32.(
            _sea_ice_mechanics_matrix(callback.cumulative_dynamic_volume_change),
        ),
        cumulative_thermodynamic_volume_change = Float32.(
            _sea_ice_mechanics_matrix(
                callback.cumulative_thermodynamic_volume_change,
            ),
        ),
        velocity_divergence = _sea_ice_mechanics_matrix(
            callback.operations.velocity_divergence,
        ),
        pressure = _sea_ice_mechanics_matrix(source.pressure),
        relaxation_parameter = _sea_ice_mechanics_matrix(
            source.relaxation_parameter,
        ),
        strain_invariant = _sea_ice_mechanics_matrix(source.strain_invariant),
        bulk_viscosity = _sea_ice_mechanics_matrix(source.bulk_viscosity),
        sigma11 = _sea_ice_mechanics_matrix(source.sigma11),
        sigma22 = _sea_ice_mechanics_matrix(source.sigma22),
        sigma12 = _sea_ice_mechanics_matrix(callback.operations.sigma12),
        ice_speed = _sea_ice_mechanics_matrix(callback.operations.ice_speed),
        solver_velocity_change = _sea_ice_mechanics_matrix(
            callback.operations.solver_velocity_change,
        ),
    )
    for name in _SEA_ICE_MECHANICS_MAP_NAMES
        push!(callback.maps[name], getproperty(sample, name))
    end
    return nothing
end

function Oceananigans.initialize!(
    callback::SeaIceMechanicsDiagnosticsCallback,
    simulation,
)
    empty!(callback.time_days)
    empty!(callback.peak_i)
    empty!(callback.peak_j)
    empty!(callback.peak_area_equivalent_thickness)
    for snapshots in values(callback.maps)
        empty!(snapshots)
    end
    fill!(parent(callback.cumulative_dynamic_volume_change), 0)
    fill!(parent(callback.cumulative_thermodynamic_volume_change), 0)
    callback.last_time = Float64(simulation.model.clock.time)
    _sample_sea_ice_mechanics!(callback, simulation)
    return nothing
end

function (callback::SeaIceMechanicsDiagnosticsCallback)(simulation)
    current_time = Float64(simulation.model.clock.time)
    if !isfinite(callback.last_time)
        callback.last_time = current_time
    end
    Δt = current_time - callback.last_time
    _accumulate_sea_ice_mechanics!(callback, Δt)
    callback.last_time = current_time
    iteration = simulation.model.clock.iteration
    at_sample = iteration % callback.sample_every_n_steps == 0
    at_end = iteration >= simulation.stop_iteration ||
             simulation.model.clock.time >= simulation.stop_time
    current_day = current_time / 86_400
    if (at_sample || at_end) &&
       (isempty(callback.time_days) || callback.time_days[end] != current_day)
        _sample_sea_ice_mechanics!(callback, simulation)
    end
    return nothing
end

function save_sea_ice_mechanics_diagnostics(callback, path::AbstractString)
    nsamples = length(callback.time_days)
    nsamples > 0 || error("sea-ice mechanics diagnostics contain no samples")
    all(length(callback.maps[name]) == nsamples for name in _SEA_ICE_MECHANICS_MAP_NAMES) ||
        error("sea-ice mechanics map histories have inconsistent lengths")
    first_map = first(callback.maps[:area_equivalent_thickness])
    Nx, Ny = size(first_map)
    mkpath(dirname(path))
    temporary_path = path * ".part"
    isfile(temporary_path) && rm(temporary_path)
    NCDataset(temporary_path, "c") do dataset
        defDim(dataset, "ocean_x", Nx)
        defDim(dataset, "ocean_y", Ny)
        defDim(dataset, "time", nsamples)
        time = defVar(dataset, "time", Float64, ("time",))
        time.attrib["units"] = "days"
        time.attrib["long_name"] = "elapsed model time"
        time[:] = callback.time_days
        active = defVar(
            dataset,
            "ocean_surface_active_mask",
            Int8,
            ("ocean_x", "ocean_y"),
        )
        active[:, :] = Int8.(callback.active_surface_mask)
        longitude = defVar(
            dataset,
            "ocean_longitude",
            Float64,
            ("ocean_x", "ocean_y"),
        )
        latitude = defVar(
            dataset,
            "ocean_latitude",
            Float64,
            ("ocean_x", "ocean_y"),
        )
        longitude.attrib["units"] = "degrees_east"
        latitude.attrib["units"] = "degrees_north"
        longitude[:, :] = callback.coordinates.longitude
        latitude[:, :] = callback.coordinates.latitude
        cell_area = defVar(
            dataset,
            "ocean_cell_area",
            Float64,
            ("ocean_x", "ocean_y"),
        )
        cell_area.attrib["units"] = "m2"
        cell_area.attrib["long_name"] = "native horizontal cell area"
        cell_area[:, :] = callback.coordinates.cell_area
        units = Dict(
            :area_equivalent_thickness => "m",
            :conditional_thickness => "m",
            :concentration => "1",
            :volume_tendency => "m s-1",
            :concentration_tendency => "s-1",
            :thermodynamic_volume_tendency => "m s-1",
            :cumulative_dynamic_volume_change => "m",
            :cumulative_thermodynamic_volume_change => "m",
            :velocity_divergence => "s-1",
            :pressure => "N m-1",
            :relaxation_parameter => "1",
            :strain_invariant => "s-1",
            :bulk_viscosity => "N s m-1",
            :sigma11 => "N m-1",
            :sigma22 => "N m-1",
            :sigma12 => "N m-1",
            :ice_speed => "m s-1",
            :solver_velocity_change => "m s-1",
        )
        for name in _SEA_ICE_MECHANICS_MAP_NAMES
            variable = defVar(
                dataset,
                String(name),
                Float32,
                ("ocean_x", "ocean_y", "time"),
            )
            variable.attrib["units"] = units[name]
            for sample in 1:nsamples
                variable[:, :, sample] = callback.maps[name][sample]
            end
        end
        peak_i = defVar(dataset, "peak_i", Int64, ("time",))
        peak_j = defVar(dataset, "peak_j", Int64, ("time",))
        peak_Ah = defVar(
            dataset,
            "peak_area_equivalent_thickness",
            Float64,
            ("time",),
        )
        peak_i[:] = callback.peak_i
        peak_j[:] = callback.peak_j
        peak_Ah.attrib["units"] = "m"
        peak_Ah[:] = callback.peak_area_equivalent_thickness
        dataset.attrib["diagnostic_phase"] =
            "after coupled update; EVP auxiliaries are from the most recent momentum solve"
        dataset.attrib["observer_mutates_live_state"] = Int8(0)
        dataset.attrib["velocity_halo_policy"] =
            "copy live storage, fill private diagnostic halos, then evaluate native metric stencil"
        dataset.attrib["tendency_accumulation"] =
            "every coupled step; current final-stage dynamic and thermodynamic tendencies applied over the just-completed interval"
        dataset.attrib["map_sample_every_n_steps"] = callback.sample_every_n_steps
        dataset.attrib["sea_ice_momentum_substeps"] = callback.momentum_substeps
        dataset.attrib["sea_ice_pressure_formulation"] =
            callback.pressure_formulation
        dataset.attrib["sea_ice_immersed_boundary_drag_coefficient"] =
            callback.immersed_boundary_drag_coefficient
        dataset.attrib["ocean_tripolar_wet_mask"] =
            callback.ocean_tripolar_wet_mask
    end
    mv(temporary_path, path; force = true)
    return path
end
