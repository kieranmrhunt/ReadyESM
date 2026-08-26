const _SEA_ICE_MOMENTUM_FORCE_OBSERVERS = IdDict{Any, Any}()

const _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES = (
    :observed_acceleration_u,
    :observed_acceleration_v,
    :coriolis_acceleration_u,
    :coriolis_acceleration_v,
    :atmosphere_stress_acceleration_u,
    :atmosphere_stress_acceleration_v,
    :ocean_drag_acceleration_u,
    :ocean_drag_acceleration_v,
    :internal_stress_acceleration_u,
    :internal_stress_acceleration_v,
    :numerical_forcing_acceleration_u,
    :numerical_forcing_acceleration_v,
    :endpoint_sum_acceleration_u,
    :endpoint_sum_acceleration_v,
    :endpoint_balance_residual_u,
    :endpoint_balance_residual_v,
)

const _SEA_ICE_MOMENTUM_FORCE_STATE_NAMES = (
    :conditional_thickness,
    :concentration,
    :area_equivalent_thickness,
    :pressure,
    :relaxation_parameter,
    :strain_invariant,
    :bulk_viscosity,
    :ice_velocity_u,
    :ice_velocity_v,
    :ocean_velocity_u,
    :ocean_velocity_v,
    :atmosphere_stress_u,
    :atmosphere_stress_v,
)

const _SEA_ICE_MOMENTUM_FORCE_HOTSPOT_NAMES = (
    _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES...,
    _SEA_ICE_MOMENTUM_FORCE_STATE_NAMES...,
)

@inline function _sea_ice_u_mass(i, j, k, grid, fields)
    return Oceananigans.Operators.ℑxᶠᵃᵃ(
        i,
        j,
        size(grid, 3),
        grid,
        ClimaSeaIce.ice_mass,
        fields.h,
        fields.ℵ,
        fields.ρ,
    )
end

@inline function _sea_ice_v_mass(i, j, k, grid, fields)
    return Oceananigans.Operators.ℑyᵃᶠᵃ(
        i,
        j,
        size(grid, 3),
        grid,
        ClimaSeaIce.ice_mass,
        fields.h,
        fields.ℵ,
        fields.ρ,
    )
end

@inline function _sea_ice_u_concentration(i, j, k, grid, fields)
    return Oceananigans.Operators.ℑxᶠᵃᵃ(
        i,
        j,
        size(grid, 3),
        grid,
        fields.ℵ,
    )
end

@inline function _sea_ice_v_concentration(i, j, k, grid, fields)
    return Oceananigans.Operators.ℑyᵃᶠᵃ(
        i,
        j,
        size(grid, 3),
        grid,
        fields.ℵ,
    )
end

@inline function _sea_ice_observed_u_acceleration(i, j, k, grid, fields, Δτ)
    return @inbounds (fields.u[i, j, 1] - fields.uⁿ[i, j, 1]) / Δτ
end

@inline function _sea_ice_observed_v_acceleration(i, j, k, grid, fields, Δτ)
    return @inbounds (fields.v[i, j, 1] - fields.vⁿ[i, j, 1]) / Δτ
end

@inline function _sea_ice_coriolis_u_acceleration(
    i,
    j,
    k,
    grid,
    coriolis,
    fields,
)
    U = (; u = fields.u, v = fields.v)
    return -Oceananigans.Coriolis.x_f_cross_U(
        i,
        j,
        size(grid, 3),
        grid,
        coriolis,
        U,
    )
end

@inline function _sea_ice_coriolis_v_acceleration(
    i,
    j,
    k,
    grid,
    coriolis,
    fields,
)
    U = (; u = fields.u, v = fields.v)
    return -Oceananigans.Coriolis.y_f_cross_U(
        i,
        j,
        size(grid, 3),
        grid,
        coriolis,
        U,
    )
end

@inline function _sea_ice_top_u_acceleration(
    i,
    j,
    k,
    grid,
    top_stress,
    fields,
    clock,
)
    mass = _sea_ice_u_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    concentration = _sea_ice_u_concentration(i, j, k, grid, fields)
    stress = ClimaSeaIce.SeaIceDynamics.x_momentum_stress(
        i,
        j,
        size(grid, 3),
        grid,
        top_stress,
        clock,
        fields,
    )
    return -stress * concentration / mass
end

@inline function _sea_ice_top_v_acceleration(
    i,
    j,
    k,
    grid,
    top_stress,
    fields,
    clock,
)
    mass = _sea_ice_v_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    concentration = _sea_ice_v_concentration(i, j, k, grid, fields)
    stress = ClimaSeaIce.SeaIceDynamics.y_momentum_stress(
        i,
        j,
        size(grid, 3),
        grid,
        top_stress,
        clock,
        fields,
    )
    return -stress * concentration / mass
end

@inline function _sea_ice_bottom_u_acceleration(
    i,
    j,
    k,
    grid,
    bottom_stress,
    fields,
    clock,
)
    mass = _sea_ice_u_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    concentration = _sea_ice_u_concentration(i, j, k, grid, fields)
    stress = ClimaSeaIce.SeaIceDynamics.x_momentum_stress(
        i,
        j,
        size(grid, 3),
        grid,
        bottom_stress,
        clock,
        fields,
    )
    return stress * concentration / mass
end

@inline function _sea_ice_bottom_v_acceleration(
    i,
    j,
    k,
    grid,
    bottom_stress,
    fields,
    clock,
)
    mass = _sea_ice_v_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    concentration = _sea_ice_v_concentration(i, j, k, grid, fields)
    stress = ClimaSeaIce.SeaIceDynamics.y_momentum_stress(
        i,
        j,
        size(grid, 3),
        grid,
        bottom_stress,
        clock,
        fields,
    )
    return stress * concentration / mass
end

@inline function _sea_ice_internal_u_acceleration(
    i,
    j,
    k,
    grid,
    rheology,
    fields,
    clock,
    immersed_boundary_condition,
)
    mass = _sea_ice_u_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    k_top = size(grid, 3)
    divergence = ClimaSeaIce.Rheologies.∂ⱼ_σ₁ⱼ(
        i,
        j,
        k_top,
        grid,
        rheology,
        clock,
        fields,
    ) + ClimaSeaIce.Rheologies.immersed_∂ⱼ_σ₁ⱼ(
        i,
        j,
        k_top,
        grid,
        immersed_boundary_condition,
        rheology,
        clock,
        fields,
    )
    return divergence / mass
end

@inline function _sea_ice_internal_v_acceleration(
    i,
    j,
    k,
    grid,
    rheology,
    fields,
    clock,
    immersed_boundary_condition,
)
    mass = _sea_ice_v_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    k_top = size(grid, 3)
    divergence = ClimaSeaIce.Rheologies.∂ⱼ_σ₂ⱼ(
        i,
        j,
        k_top,
        grid,
        rheology,
        clock,
        fields,
    ) + ClimaSeaIce.Rheologies.immersed_∂ⱼ_σ₂ⱼ(
        i,
        j,
        k_top,
        grid,
        immersed_boundary_condition,
        rheology,
        clock,
        fields,
    )
    return divergence / mass
end

@inline function _sea_ice_numerical_u_acceleration(
    i,
    j,
    k,
    grid,
    rheology,
    fields,
    forcing,
    Δτ,
)
    mass = _sea_ice_u_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    return ClimaSeaIce.Rheologies.sum_of_forcing_u(
        i,
        j,
        size(grid, 3),
        grid,
        rheology,
        forcing,
        fields,
        Δτ,
    )
end


@inline function _sea_ice_numerical_v_acceleration(
    i,
    j,
    k,
    grid,
    rheology,
    fields,
    forcing,
    Δτ,
)
    mass = _sea_ice_v_mass(i, j, k, grid, fields)
    mass <= 0 && return zero(grid)
    return ClimaSeaIce.Rheologies.sum_of_forcing_v(
        i,
        j,
        size(grid, 3),
        grid,
        rheology,
        forcing,
        fields,
        Δτ,
    )
end

@inline function _sea_ice_centered_u_term(i, j, k, grid, term, arguments...)
    return Oceananigans.Operators.ℑxᶜᵃᵃ(i, j, k, grid, term, arguments...)
end

@inline function _sea_ice_centered_v_term(i, j, k, grid, term, arguments...)
    return Oceananigans.Operators.ℑyᵃᶜᵃ(i, j, k, grid, term, arguments...)
end

@inline _sea_ice_field_value(i, j, k, grid, field) = @inbounds field[i, j, 1]

@inline _sea_ice_area_equivalent_thickness(i, j, k, grid, h, concentration) =
    @inbounds h[i, j, 1] * concentration[i, j, 1]

@inline _sea_ice_force_sum(i, j, k, grid, a, b, c, d, e) =
    @inbounds a[i, j, k] + b[i, j, k] + c[i, j, k] + d[i, j, k] + e[i, j, k]

@inline _sea_ice_force_residual(i, j, k, grid, observed, total) =
    @inbounds observed[i, j, k] - total[i, j, k]

function _sea_ice_force_operation(function_, grid, arguments)
    operation = Oceananigans.KernelFunctionOperation{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(function_, grid, arguments)
    # Every force diagnostic is defined over the complete horizontal sea-ice
    # grid. Avoid AbstractOperations' window-index inference here: the
    # arguments deliberately mix 3D centered state, 2D surface fields and
    # ConstantFields whose formal locations contain `Nothing`. Intersecting
    # those locations attempts an undefined `Center <- Nothing` restriction
    # on a vertically resolved production grid, despite no window being
    # requested. Explicit full-domain indices preserve the intended extent.
    return Oceananigans.Field(operation; indices = (:, :, :))
end

function _sea_ice_force_component_operations(model, Δτ)
    dynamics = model.dynamics
    auxiliaries = dynamics.auxiliaries.fields
    fields = merge(
        auxiliaries,
        model.velocities,
        (;
            h = model.ice_thickness,
            ℵ = model.ice_concentration,
            ρ = model.sea_ice_density,
        ),
    )
    grid = model.grid
    clock = model.clock
    coriolis = dynamics.coriolis
    rheology = dynamics.rheology
    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom
    u_immersed = model.velocities.u.boundary_conditions.immersed
    v_immersed = model.velocities.v.boundary_conditions.immersed

    observed_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (_sea_ice_observed_u_acceleration, fields, Δτ),
    )
    observed_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (_sea_ice_observed_v_acceleration, fields, Δτ),
    )
    coriolis_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (_sea_ice_coriolis_u_acceleration, coriolis, fields),
    )
    coriolis_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (_sea_ice_coriolis_v_acceleration, coriolis, fields),
    )
    atmosphere_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (_sea_ice_top_u_acceleration, top_stress, fields, clock),
    )
    atmosphere_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (_sea_ice_top_v_acceleration, top_stress, fields, clock),
    )
    ocean_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (_sea_ice_bottom_u_acceleration, bottom_stress, fields, clock),
    )
    ocean_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (_sea_ice_bottom_v_acceleration, bottom_stress, fields, clock),
    )
    internal_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (
            _sea_ice_internal_u_acceleration,
            rheology,
            fields,
            clock,
            u_immersed,
        ),
    )
    internal_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (
            _sea_ice_internal_v_acceleration,
            rheology,
            fields,
            clock,
            v_immersed,
        ),
    )
    numerical_u = _sea_ice_force_operation(
        _sea_ice_centered_u_term,
        grid,
        (
            _sea_ice_numerical_u_acceleration,
            rheology,
            fields,
            model.forcing.u,
            Δτ,
        ),
    )
    numerical_v = _sea_ice_force_operation(
        _sea_ice_centered_v_term,
        grid,
        (
            _sea_ice_numerical_v_acceleration,
            rheology,
            fields,
            model.forcing.v,
            Δτ,
        ),
    )
    total_u = _sea_ice_force_operation(
        _sea_ice_force_sum,
        grid,
        (coriolis_u, atmosphere_u, ocean_u, internal_u, numerical_u),
    )
    total_v = _sea_ice_force_operation(
        _sea_ice_force_sum,
        grid,
        (coriolis_v, atmosphere_v, ocean_v, internal_v, numerical_v),
    )
    residual_u = _sea_ice_force_operation(
        _sea_ice_force_residual,
        grid,
        (observed_u, total_u),
    )
    residual_v = _sea_ice_force_operation(
        _sea_ice_force_residual,
        grid,
        (observed_v, total_v),
    )

    force_operations = (
        observed_acceleration_u = observed_u,
        observed_acceleration_v = observed_v,
        coriolis_acceleration_u = coriolis_u,
        coriolis_acceleration_v = coriolis_v,
        atmosphere_stress_acceleration_u = atmosphere_u,
        atmosphere_stress_acceleration_v = atmosphere_v,
        ocean_drag_acceleration_u = ocean_u,
        ocean_drag_acceleration_v = ocean_v,
        internal_stress_acceleration_u = internal_u,
        internal_stress_acceleration_v = internal_v,
        numerical_forcing_acceleration_u = numerical_u,
        numerical_forcing_acceleration_v = numerical_v,
        endpoint_sum_acceleration_u = total_u,
        endpoint_sum_acceleration_v = total_v,
        endpoint_balance_residual_u = residual_u,
        endpoint_balance_residual_v = residual_v,
    )
    state_operations = (
        conditional_thickness = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (model.ice_thickness,),
        ),
        concentration = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (model.ice_concentration,),
        ),
        area_equivalent_thickness = _sea_ice_force_operation(
            _sea_ice_area_equivalent_thickness,
            grid,
            (model.ice_thickness, model.ice_concentration),
        ),
        pressure = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (auxiliaries.P,),
        ),
        relaxation_parameter = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (auxiliaries.α,),
        ),
        strain_invariant = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (auxiliaries.Δ,),
        ),
        bulk_viscosity = _sea_ice_force_operation(
            _sea_ice_field_value,
            grid,
            (auxiliaries.ζᶜᶜᶜ,),
        ),
        ice_velocity_u = _sea_ice_force_operation(
            _sea_ice_centered_u_term,
            grid,
            (_sea_ice_field_value, model.velocities.u),
        ),
        ice_velocity_v = _sea_ice_force_operation(
            _sea_ice_centered_v_term,
            grid,
            (_sea_ice_field_value, model.velocities.v),
        ),
        ocean_velocity_u = _sea_ice_force_operation(
            _sea_ice_centered_u_term,
            grid,
            (
                _sea_ice_field_value,
                dynamics.external_momentum_stresses.bottom.uₑ,
            ),
        ),
        ocean_velocity_v = _sea_ice_force_operation(
            _sea_ice_centered_v_term,
            grid,
            (
                _sea_ice_field_value,
                dynamics.external_momentum_stresses.bottom.vₑ,
            ),
        ),
        atmosphere_stress_u = _sea_ice_force_operation(
            _sea_ice_centered_u_term,
            grid,
            (
                _sea_ice_field_value,
                dynamics.external_momentum_stresses.top.u,
            ),
        ),
        atmosphere_stress_v = _sea_ice_force_operation(
            _sea_ice_centered_v_term,
            grid,
            (
                _sea_ice_field_value,
                dynamics.external_momentum_stresses.top.v,
            ),
        ),
    )
    return merge(force_operations, state_operations)
end

Base.@kwdef mutable struct SeaIceMomentumForceDiagnosticsCallback{M, O, C} <:
                           Function
    model::M
    operations::O
    coordinates::C
    sample_every_n_steps::Int
    hotspot_i::Int
    hotspot_j::Int
    time_days::Vector{Float64} = Float64[]
    hotspot::Dict{Symbol, Vector{Float64}} = Dict(
        name => Float64[] for name in _SEA_ICE_MOMENTUM_FORCE_HOTSPOT_NAMES
    )
    map_time_days::Vector{Float64} = Float64[]
    maps::Dict{Symbol, Vector{Matrix{Float32}}} = Dict(
        name => Matrix{Float32}[]
        for name in _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES
    )
end

function SeaIceMomentumForceDiagnosticsCallback(
    earth,
    Δt;
    hotspot = (253, 179),
    sample_every_n_steps = max(1, round(Int, 86_400 / Float64(Δt))),
)
    model = earth.sea_ice.model
    model.dynamics isa ClimaSeaIce.SeaIceDynamics.SeaIceMomentumEquation ||
        error("sea-ice force diagnostics require an active momentum equation")
    sample_every_n_steps > 0 ||
        throw(ArgumentError("sample_every_n_steps must be positive"))
    Nx, Ny, _ = size(model.grid)
    i, j = hotspot
    1 <= i <= Nx || throw(BoundsError((1:Nx, 1:Ny), hotspot))
    1 <= j <= Ny || throw(BoundsError((1:Nx, 1:Ny), hotspot))
    active_surface_mask = _ocean_surface_active_mask(earth.ocean.model.grid)
    active_surface_mask[i, j] == 1 ||
        error("sea-ice force diagnostic hotspot must be an active ocean cell")
    longitude, latitude, _ = Oceananigans.nodes(model.ice_thickness)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    if ndims(longitude) == 1 && ndims(latitude) == 1
        longitude = repeat(reshape(longitude, :, 1), 1, Ny)
        latitude = repeat(reshape(latitude, 1, :), Nx, 1)
    end
    operations = _sea_ice_force_component_operations(model, Float32(Δt))
    return SeaIceMomentumForceDiagnosticsCallback(
        ;
        model,
        operations,
        coordinates = (; longitude, latitude, active_surface_mask),
        sample_every_n_steps,
        hotspot_i = i,
        hotspot_j = j,
    )
end

function install_sea_ice_momentum_force_diagnostics!(callback)
    model = callback.model
    haskey(_SEA_ICE_MOMENTUM_FORCE_OBSERVERS, model) &&
        error("a sea-ice momentum force observer is already installed")
    _SEA_ICE_MOMENTUM_FORCE_OBSERVERS[model] = callback
    return callback
end

function remove_sea_ice_momentum_force_diagnostics!(callback)
    model = callback.model
    get(_SEA_ICE_MOMENTUM_FORCE_OBSERVERS, model, nothing) === callback &&
        delete!(_SEA_ICE_MOMENTUM_FORCE_OBSERVERS, model)
    return nothing
end

function _sea_ice_force_hotspot_value(field, i, j)
    value = Array(Oceananigans.interior(field, i:i, j:j, 1:1))
    return Float64(only(value))
end

function _observe_sea_ice_momentum_endpoint!(model, Δτ)
    callback = get(_SEA_ICE_MOMENTUM_FORCE_OBSERVERS, model, nothing)
    isnothing(callback) && return nothing
    model.clock.stage == model.timestepper.Nstages || return nothing
    isapprox(Float64(Δτ), Float64(model.clock.last_stage_Δt); rtol = 0, atol = 0) ||
        error("sea-ice momentum observer received an inconsistent final-stage timestep")

    return _observe_sea_ice_force_callback!(callback, model, Δτ)
end

function _observe_sea_ice_force_callback!(
    callback::SeaIceMomentumForceDiagnosticsCallback,
    model,
    Δτ,
)
    architecture = Oceananigans.Architectures.architecture(model.grid)
    for operation in values(callback.operations)
        Oceananigans.compute!(operation)
    end
    Oceananigans.Architectures.synchronize(architecture)

    time_days = (Float64(model.clock.time) + Float64(Δτ)) / 86_400
    push!(callback.time_days, time_days)
    for name in _SEA_ICE_MOMENTUM_FORCE_HOTSPOT_NAMES
        value = _sea_ice_force_hotspot_value(
            getproperty(callback.operations, name),
            callback.hotspot_i,
            callback.hotspot_j,
        )
        isfinite(value) || error(
            "non-finite sea-ice momentum force diagnostic $name at " *
            "($(callback.hotspot_i),$(callback.hotspot_j))",
        )
        push!(callback.hotspot[name], value)
    end

    completed_iteration = model.clock.iteration + 1
    if completed_iteration % callback.sample_every_n_steps == 0
        push!(callback.map_time_days, time_days)
        for name in _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES
            push!(
                callback.maps[name],
                _sea_ice_mechanics_matrix(getproperty(callback.operations, name)),
            )
        end
    end
    return nothing
end

function save_sea_ice_momentum_force_diagnostics(
    callback,
    path::AbstractString,
)
    nsteps = length(callback.time_days)
    nsteps > 0 || error("sea-ice momentum force diagnostics contain no samples")
    all(length(callback.hotspot[name]) == nsteps for name in
        _SEA_ICE_MOMENTUM_FORCE_HOTSPOT_NAMES) ||
        error("sea-ice momentum force hotspot histories have inconsistent lengths")
    nmaps = length(callback.map_time_days)
    all(length(callback.maps[name]) == nmaps for name in
        _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES) ||
        error("sea-ice momentum force map histories have inconsistent lengths")

    Nx, Ny = size(callback.coordinates.active_surface_mask)
    mkpath(dirname(path))
    temporary_path = path * ".part"
    isfile(temporary_path) && rm(temporary_path)
    NCDataset(temporary_path, "c") do dataset
        defDim(dataset, "ocean_x", Nx)
        defDim(dataset, "ocean_y", Ny)
        defDim(dataset, "step", nsteps)
        defDim(dataset, "map_time", nmaps)
        step_time = defVar(dataset, "step_time", Float64, ("step",))
        map_time = defVar(dataset, "map_time", Float64, ("map_time",))
        step_time.attrib["units"] = "days"
        map_time.attrib["units"] = "days"
        step_time[:] = callback.time_days
        map_time[:] = callback.map_time_days
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
        active = defVar(
            dataset,
            "ocean_surface_active_mask",
            Int8,
            ("ocean_x", "ocean_y"),
        )
        longitude.attrib["units"] = "degrees_east"
        latitude.attrib["units"] = "degrees_north"
        longitude[:, :] = callback.coordinates.longitude
        latitude[:, :] = callback.coordinates.latitude
        active[:, :] = Int8.(callback.coordinates.active_surface_mask)
        for name in _SEA_ICE_MOMENTUM_FORCE_COMPONENT_NAMES
            hotspot = defVar(
                dataset,
                "hotspot_$(name)",
                Float64,
                ("step",),
            )
            hotspot.attrib["units"] = "m s-2"
            hotspot[:] = callback.hotspot[name]
            maps = defVar(
                dataset,
                "map_$(name)",
                Float32,
                ("ocean_x", "ocean_y", "map_time"),
            )
            maps.attrib["units"] = "m s-2"
            for sample in 1:nmaps
                maps[:, :, sample] = callback.maps[name][sample]
            end
        end
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
        for name in _SEA_ICE_MOMENTUM_FORCE_STATE_NAMES
            hotspot = defVar(
                dataset,
                "hotspot_$(name)",
                Float64,
                ("step",),
            )
            hotspot.attrib["units"] = state_units[name]
            hotspot[:] = callback.hotspot[name]
        end
        dataset.attrib["hotspot_i"] = callback.hotspot_i
        dataset.attrib["hotspot_j"] = callback.hotspot_j
        dataset.attrib["hotspot_longitude_degrees_east"] =
            callback.coordinates.longitude[callback.hotspot_i, callback.hotspot_j]
        dataset.attrib["hotspot_latitude_degrees_north"] =
            callback.coordinates.latitude[callback.hotspot_i, callback.hotspot_j]
        dataset.attrib["diagnostic_phase"] =
            "after final SplitRK momentum solve and before dynamic/thermodynamic ice update"
        dataset.attrib["force_convention"] =
            "native-grid endpoint accelerations from installed ClimaSeaIce operators; ocean drag includes semi-implicit endpoint velocity"
        dataset.attrib["residual_convention"] =
            "observed final-stage velocity increment minus endpoint force sum; nonzero residual includes split-explicit substep time integration"
        dataset.attrib["observer_mutates_live_state"] = Int8(0)
        dataset.attrib["map_sample_every_n_steps"] = callback.sample_every_n_steps
    end
    mv(temporary_path, path; force = true)
    return path
end
