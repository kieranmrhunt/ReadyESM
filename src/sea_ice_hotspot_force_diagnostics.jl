const _SEA_ICE_HOTSPOT_FORCE_DEFAULT_INDICES = (
    (256, 175), (257, 175), (258, 175),
    (256, 176), (257, 176), (258, 176),
    (256, 177), (257, 177), (258, 177),
    (108, 180),
)

const _SEA_ICE_COMPACT_FORCE_DIAGNOSTIC_NAMES = (
    _SEA_ICE_MOMENTUM_FORCE_HOTSPOT_NAMES...,
    :surface_tilt_acceleration_u,
    :surface_tilt_acceleration_v,
    :sea_surface_height,
)

@inline function _sea_ice_surface_tilt_u_acceleration(
    i, j, k, grid, sea_surface_height, gravitational_acceleration,
)
    return -gravitational_acceleration *
           Oceananigans.Operators.∂xᶠᶜᶜ(i, j, k, grid, sea_surface_height)
end

@inline function _sea_ice_surface_tilt_v_acceleration(
    i, j, k, grid, sea_surface_height, gravitational_acceleration,
)
    return -gravitational_acceleration *
           Oceananigans.Operators.∂yᶜᶠᶜ(i, j, k, grid, sea_surface_height)
end

KernelAbstractions.@kernel function _sample_sea_ice_hotspot_forces_kernel!(
    output,
    hotspot_i,
    hotspot_j,
    grid,
    Δτ,
    rheology,
    fields,
    clock,
    coriolis,
    u_immersed,
    v_immersed,
    top_stress,
    bottom_stress,
    u_forcing,
    v_forcing,
    sea_surface_height,
    gravitational_acceleration,
)
    hotspot = @index(Global, Linear)
    i = @inbounds hotspot_i[hotspot]
    j = @inbounds hotspot_j[hotspot]
    k = 1

    observed_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_observed_u_acceleration, fields, Δτ,
    )
    observed_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_observed_v_acceleration, fields, Δτ,
    )
    coriolis_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_coriolis_u_acceleration, coriolis, fields,
    )
    coriolis_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_coriolis_v_acceleration, coriolis, fields,
    )
    atmosphere_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_top_u_acceleration, top_stress, fields, clock,
    )
    atmosphere_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_top_v_acceleration, top_stress, fields, clock,
    )
    ocean_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_bottom_u_acceleration, bottom_stress, fields, clock,
    )
    ocean_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_bottom_v_acceleration, bottom_stress, fields, clock,
    )
    internal_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_internal_u_acceleration,
        rheology, fields, clock, u_immersed,
    )
    internal_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_internal_v_acceleration,
        rheology, fields, clock, v_immersed,
    )
    numerical_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_numerical_u_acceleration,
        rheology, fields, u_forcing, Δτ,
    )
    numerical_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_numerical_v_acceleration,
        rheology, fields, v_forcing, Δτ,
    )
    endpoint_u = coriolis_u + atmosphere_u + ocean_u + internal_u + numerical_u
    endpoint_v = coriolis_v + atmosphere_v + ocean_v + internal_v + numerical_v
    surface_tilt_u = _sea_ice_centered_u_term(
        i, j, k, grid, _sea_ice_surface_tilt_u_acceleration,
        sea_surface_height, gravitational_acceleration,
    )
    surface_tilt_v = _sea_ice_centered_v_term(
        i, j, k, grid, _sea_ice_surface_tilt_v_acceleration,
        sea_surface_height, gravitational_acceleration,
    )

    @inbounds begin
        output[hotspot, 1] = observed_u
        output[hotspot, 2] = observed_v
        output[hotspot, 3] = coriolis_u
        output[hotspot, 4] = coriolis_v
        output[hotspot, 5] = atmosphere_u
        output[hotspot, 6] = atmosphere_v
        output[hotspot, 7] = ocean_u
        output[hotspot, 8] = ocean_v
        output[hotspot, 9] = internal_u
        output[hotspot, 10] = internal_v
        output[hotspot, 11] = numerical_u
        output[hotspot, 12] = numerical_v
        output[hotspot, 13] = endpoint_u
        output[hotspot, 14] = endpoint_v
        output[hotspot, 15] = observed_u - endpoint_u
        output[hotspot, 16] = observed_v - endpoint_v
        output[hotspot, 17] = fields.h[i, j, 1]
        output[hotspot, 18] = fields.ℵ[i, j, 1]
        output[hotspot, 19] = fields.h[i, j, 1] * fields.ℵ[i, j, 1]
        output[hotspot, 20] = fields.P[i, j, 1]
        output[hotspot, 21] = fields.α[i, j, 1]
        output[hotspot, 22] = fields.Δ[i, j, 1]
        output[hotspot, 23] = fields.ζᶜᶜᶜ[i, j, 1]
        output[hotspot, 24] = _sea_ice_centered_u_term(
            i, j, k, grid, _sea_ice_field_value, fields.u,
        )
        output[hotspot, 25] = _sea_ice_centered_v_term(
            i, j, k, grid, _sea_ice_field_value, fields.v,
        )
        output[hotspot, 26] = _sea_ice_centered_u_term(
            i, j, k, grid, _sea_ice_field_value, bottom_stress.uₑ,
        )
        output[hotspot, 27] = _sea_ice_centered_v_term(
            i, j, k, grid, _sea_ice_field_value, bottom_stress.vₑ,
        )
        output[hotspot, 28] = _sea_ice_centered_u_term(
            i, j, k, grid, _sea_ice_field_value, top_stress.u,
        )
        output[hotspot, 29] = _sea_ice_centered_v_term(
            i, j, k, grid, _sea_ice_field_value, top_stress.v,
        )
        output[hotspot, 30] = surface_tilt_u
        output[hotspot, 31] = surface_tilt_v
        output[hotspot, 32] = sea_surface_height[i, j, 1]
    end
end

Base.@kwdef mutable struct SeaIceHotspotForceDiagnosticsCallback{M, B, I, C, E, G}
    model::M
    device_buffer::B
    hotspot_indices::I
    coordinates::C
    sea_surface_height::E
    gravitational_acceleration::G
    sample_every_n_steps::Int
    time_days::Vector{Float64} = Float64[]
    samples::Vector{Matrix{Float32}} = Matrix{Float32}[]
end

function SeaIceHotspotForceDiagnosticsCallback(
    earth,
    Δt;
    hotspot_indices = _SEA_ICE_HOTSPOT_FORCE_DEFAULT_INDICES,
    sample_every_n_steps = max(1, round(Int, 3_600 / Float64(Δt))),
)
    model = earth.sea_ice.model
    model.dynamics isa ClimaSeaIce.SeaIceDynamics.SeaIceMomentumEquation ||
        error("sea-ice hotspot force diagnostics require an active momentum equation")
    sample_every_n_steps > 0 ||
        throw(ArgumentError("sample_every_n_steps must be positive"))
    indices = Tuple((Int(i), Int(j)) for (i, j) in hotspot_indices)
    isempty(indices) && throw(ArgumentError("at least one force hotspot is required"))
    allunique(indices) || throw(ArgumentError("force hotspot indices must be unique"))
    Nx, Ny, _ = size(model.grid)
    for index in indices
        1 <= index[1] <= Nx || throw(BoundsError((1:Nx, 1:Ny), index))
        1 <= index[2] <= Ny || throw(BoundsError((1:Nx, 1:Ny), index))
    end

    active_surface_mask = _ocean_surface_active_mask(earth.ocean.model.grid)
    longitude, latitude, _ = Oceananigans.nodes(model.ice_thickness)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    if ndims(longitude) == 1 && ndims(latitude) == 1
        longitude = repeat(reshape(longitude, :, 1), 1, Ny)
        latitude = repeat(reshape(latitude, 1, :), Nx, 1)
    end
    coordinates = (
        longitude = [longitude[i, j] for (i, j) in indices],
        latitude = [latitude[i, j] for (i, j) in indices],
        active = Int8[active_surface_mask[i, j] == 1 for (i, j) in indices],
    )
    architecture = Oceananigans.Architectures.architecture(model.grid)
    ocean_model = earth.ocean.model
    sea_surface_height = ocean_model.free_surface.displacement
    gravitational_acceleration =
        ocean_model.buoyancy.formulation.gravitational_acceleration
    device_buffer = Oceananigans.on_architecture(
        architecture,
        zeros(Float32, length(indices), length(_SEA_ICE_COMPACT_FORCE_DIAGNOSTIC_NAMES)),
    )
    return SeaIceHotspotForceDiagnosticsCallback(
        ; model, device_buffer, hotspot_indices = indices, coordinates,
        sea_surface_height, gravitational_acceleration,
        sample_every_n_steps,
    )
end

function _observe_sea_ice_force_callback!(
    callback::SeaIceHotspotForceDiagnosticsCallback,
    model,
    Δτ,
)
    completed_iteration = model.clock.iteration + 1
    if completed_iteration != 1 &&
       completed_iteration % callback.sample_every_n_steps != 0
        return nothing
    end

    dynamics = model.dynamics
    fields = merge(
        dynamics.auxiliaries.fields,
        model.velocities,
        (; h = model.ice_thickness,
           ℵ = model.ice_concentration,
           ρ = model.sea_ice_density),
    )
    indices = callback.hotspot_indices
    arguments = (
        callback.device_buffer,
        Tuple(first(index) for index in indices),
        Tuple(last(index) for index in indices),
        model.grid,
        Float32(Δτ),
        dynamics.rheology,
        fields,
        model.clock,
        dynamics.coriolis,
        model.velocities.u.boundary_conditions.immersed,
        model.velocities.v.boundary_conditions.immersed,
        dynamics.external_momentum_stresses.top,
        dynamics.external_momentum_stresses.bottom,
        model.forcing.u,
        model.forcing.v,
        callback.sea_surface_height,
        callback.gravitational_acceleration,
    )
    architecture = Oceananigans.Architectures.architecture(model.grid)
    device_arguments = Oceananigans.Architectures.convert_to_device(
        architecture,
        arguments,
    )
    backend = KernelAbstractions.get_backend(callback.device_buffer)
    _sample_sea_ice_hotspot_forces_kernel!(backend)(
        device_arguments...;
        ndrange = length(indices),
    )
    KernelAbstractions.synchronize(backend)
    sample = Float32.(Array(callback.device_buffer))
    all(isfinite, sample) || error(
        "non-finite compact sea-ice momentum force diagnostic at iteration " *
        string(completed_iteration),
    )
    push!(callback.samples, sample)
    push!(
        callback.time_days,
        (Float64(model.clock.time) + Float64(Δτ)) / 86_400,
    )
    return nothing
end

function save_sea_ice_hotspot_force_diagnostics(
    callback::SeaIceHotspotForceDiagnosticsCallback,
    path::AbstractString,
)
    nsamples = length(callback.time_days)
    nsamples > 0 || error("sea-ice hotspot force diagnostics contain no samples")
    length(callback.samples) == nsamples ||
        error("sea-ice hotspot force sample histories have inconsistent lengths")
    nhotspots = length(callback.hotspot_indices)
    nvariables = length(_SEA_ICE_COMPACT_FORCE_DIAGNOSTIC_NAMES)
    all(size(sample) == (nhotspots, nvariables) for sample in callback.samples) ||
        error("sea-ice hotspot force sample has an invalid shape")
    values = Array{Float32}(undef, nhotspots, nsamples, nvariables)
    for sample in 1:nsamples
        values[:, sample, :] = callback.samples[sample]
    end

    mkpath(dirname(path))
    temporary_path = path * ".part"
    isfile(temporary_path) && rm(temporary_path)
    NCDataset(temporary_path, "c") do dataset
        defDim(dataset, "hotspot", nhotspots)
        defDim(dataset, "sample", nsamples)
        time = defVar(dataset, "time", Float64, ("sample",))
        time.attrib["units"] = "days"
        time[:] = callback.time_days
        hotspot_i = defVar(dataset, "hotspot_i", Int32, ("hotspot",))
        hotspot_j = defVar(dataset, "hotspot_j", Int32, ("hotspot",))
        longitude = defVar(dataset, "longitude", Float64, ("hotspot",))
        latitude = defVar(dataset, "latitude", Float64, ("hotspot",))
        active = defVar(dataset, "active_ocean_cell", Int8, ("hotspot",))
        hotspot_i[:] = Int32[first(index) for index in callback.hotspot_indices]
        hotspot_j[:] = Int32[last(index) for index in callback.hotspot_indices]
        longitude[:] = callback.coordinates.longitude
        latitude[:] = callback.coordinates.latitude
        active[:] = callback.coordinates.active
        longitude.attrib["units"] = "degrees_east"
        latitude.attrib["units"] = "degrees_north"
        state_units = Dict(
            :conditional_thickness => "m",
            :concentration => "1",
            :area_equivalent_thickness => "m",
            :pressure => "N m-1",
            :relaxation_parameter => "1",
            :strain_invariant => "s-1",
            :bulk_viscosity => "N s m-1",
            :ice_velocity_u => "m s-1",
            :ice_velocity_v => "m s-1",
            :ocean_velocity_u => "m s-1",
            :ocean_velocity_v => "m s-1",
            :atmosphere_stress_u => "N m-2",
            :atmosphere_stress_v => "N m-2",
        )
        for (column, name) in enumerate(_SEA_ICE_COMPACT_FORCE_DIAGNOSTIC_NAMES)
            variable = defVar(
                dataset,
                string(name),
                Float32,
                ("hotspot", "sample"),
            )
            variable.attrib["units"] = if column <= 16 || name in (
                :surface_tilt_acceleration_u,
                :surface_tilt_acceleration_v,
            )
                "m s-2"
            elseif name == :sea_surface_height
                "m"
            else
                state_units[name]
            end
            variable[:, :] = values[:, :, column]
        end
        dataset.attrib["diagnostic_phase"] =
            "after final SplitRK momentum solve and before dynamic/thermodynamic ice update"
        dataset.attrib["force_convention"] =
            "native-grid centered endpoint accelerations; ocean drag includes semi-implicit endpoint velocity"
        dataset.attrib["residual_convention"] =
            "observed final-stage velocity increment minus endpoint force sum; residual includes split-explicit substep time integration"
        dataset.attrib["observer_mutates_live_state"] = Int8(0)
        dataset.attrib["sample_every_n_steps"] = callback.sample_every_n_steps
        dataset.attrib["full_domain_force_fields_constructed"] = Int8(0)
        dataset.attrib["kernel_launches_per_sample"] = Int32(1)
    end
    mv(temporary_path, path; force = true)
    return path
end
