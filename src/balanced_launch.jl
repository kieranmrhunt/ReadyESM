"""
Mutable clock state for a balanced coupled launch.

`applied_fraction` has explicit prognostic-state methods so, once the schedule
is attached to the coupled restart graph, a segmented run cannot apply an
analysis increment twice. The physical component increments and their source
ledgers are attached separately; this object owns only the normalized time
distribution.
"""
mutable struct BalancedLaunchSchedule
    start_time_seconds::Float64
    window_seconds::Float64
    applied_fraction::Float64
end

function BalancedLaunchSchedule(
    ; start_time_seconds = 0.0,
    window_seconds,
    applied_fraction = 0.0,
)
    start = Float64(start_time_seconds)
    window = Float64(window_seconds)
    applied = Float64(applied_fraction)
    isfinite(start) || throw(ArgumentError(
        "balanced-launch start time must be finite",
    ))
    isfinite(window) && window >= 0 || throw(ArgumentError(
        "balanced-launch window must be finite and nonnegative",
    ))
    isfinite(applied) && 0 <= applied <= 1 || throw(ArgumentError(
        "balanced-launch applied fraction must lie in [0, 1]",
    ))
    return BalancedLaunchSchedule(start, window, applied)
end

"""
Return the normalized cumulative analysis increment at `elapsed_seconds`.

For a positive window this is the integral of a raised-cosine tendency,
`1 - cos(2πt / W)`, normalized over the window `W`. The tendency is zero at
both ends, nonnegative everywhere, and integrates to exactly one. A zero-width
window represents the instantaneous launch branch.
"""
function _balanced_launch_cumulative_fraction(
    elapsed_seconds::Real,
    window_seconds::Real,
)
    elapsed = Float64(elapsed_seconds)
    window = Float64(window_seconds)
    isfinite(elapsed) || throw(ArgumentError(
        "balanced-launch elapsed time must be finite",
    ))
    isfinite(window) && window >= 0 || throw(ArgumentError(
        "balanced-launch window must be finite and nonnegative",
    ))
    elapsed < 0 && return 0.0
    window == 0 && return 1.0
    elapsed >= window && return 1.0
    elapsed == 0 && return 0.0
    phase = elapsed / window
    return phase - sinpi(2phase) / (2π)
end

"""
Return and record the fraction of the analysis increment due at `time_seconds`.

The value is computed from the desired cumulative fraction minus the fraction
already applied. This makes irregular timesteps and exact checkpoint/restart
boundaries equivalent to a continuous run. Once the window closes, the final
call applies `1 - applied_fraction`, eliminating cumulative roundoff from the
analysis total.
"""
function _next_balanced_launch_fraction!(
    schedule::BalancedLaunchSchedule,
    time_seconds::Real,
)
    time = Float64(time_seconds)
    isfinite(time) || throw(ArgumentError(
        "balanced-launch model time must be finite",
    ))
    desired = _balanced_launch_cumulative_fraction(
        time - schedule.start_time_seconds,
        schedule.window_seconds,
    )
    tolerance = 64eps(Float64)
    desired + tolerance >= schedule.applied_fraction || error(
        "balanced-launch model time moved behind its checkpointed applied fraction",
    )
    desired <= schedule.applied_fraction && return 0.0
    fraction = desired == 1.0 ?
        1.0 - schedule.applied_fraction :
        desired - schedule.applied_fraction
    schedule.applied_fraction = desired
    return fraction
end

Oceananigans.prognostic_state(schedule::BalancedLaunchSchedule) = (
    start_time_seconds = schedule.start_time_seconds,
    window_seconds = schedule.window_seconds,
    applied_fraction = schedule.applied_fraction,
)

function Oceananigans.restore_prognostic_state!(
    schedule::BalancedLaunchSchedule,
    state,
)
    restored = BalancedLaunchSchedule(
        ; start_time_seconds = state.start_time_seconds,
        window_seconds = state.window_seconds,
        applied_fraction = state.applied_fraction,
    )
    schedule.start_time_seconds = restored.start_time_seconds
    schedule.window_seconds = restored.window_seconds
    schedule.applied_fraction = restored.applied_fraction
    return schedule
end

const _BALANCED_ATMOSPHERE_SPECTRAL_NAMES = (
    :vorticity,
    :divergence,
    :temperature,
    :pressure,
    :humidity,
)

const _BALANCED_ATMOSPHERE_GRID_NAMES = (
    :vorticity,
    :divergence,
    :temperature,
    :pressure,
    :humidity,
    :u,
    :v,
)

"""Native SpeedyWeather analysis increment plus immutable source provenance."""
struct BalancedAtmosphereIncrement{S, G, P, T}
    spectral::S
    grid::G
    provenance::P
    source_totals::T
end

function _balanced_launch_difference(target, reference, description)
    size(target) == size(reference) || throw(DimensionMismatch(
        "$description target size $(size(target)) does not match " *
        "reference size $(size(reference))",
    ))
    eltype(target) == eltype(reference) || throw(ArgumentError(
        "$description target eltype $(eltype(target)) does not match " *
        "reference eltype $(eltype(reference))",
    ))
    difference = similar(target)
    difference.data .= target.data .- reference.data
    return difference
end

"""
Area-weighted atmospheric mass and energy inventory in physical grid space.

The energy diagnostic separates dry sensible enthalpy, vapour latent energy,
kinetic energy and geopotential energy. Their sum is a transparent analysis-
source ledger, not a claim to reproduce every term of SpeedyWeather's discrete
Hamiltonian. Cloud condensate is absent because the selected cloud scheme has
no prognostic condensate reservoir.
"""
function _balanced_atmosphere_inventory(variables, model)
    temperature = Float64.(Array(variables.grid.temperature.data))
    humidity = Float64.(Array(variables.grid.humidity.data))
    zonal_velocity = Float64.(Array(variables.grid.u.data))
    meridional_velocity = Float64.(Array(variables.grid.v.data))
    geopotential = Float64.(Array(variables.grid.geopotential.data))
    log_surface_pressure = Float64.(Array(variables.grid.pressure.data))
    size(temperature) == size(humidity) == size(zonal_velocity) ==
        size(meridional_velocity) == size(geopotential) || throw(
            DimensionMismatch(
                "balanced atmosphere inventory fields do not share one grid",
            ),
        )
    length(log_surface_pressure) == size(temperature, 1) || throw(
        DimensionMismatch(
            "balanced atmosphere pressure does not span the horizontal grid",
        ),
    )
    all(isfinite, temperature) || error(
        "balanced atmosphere temperature inventory is nonfinite",
    )
    all(isfinite, humidity) || error(
        "balanced atmosphere humidity inventory is nonfinite",
    )
    all(isfinite, log_surface_pressure) || error(
        "balanced atmosphere pressure inventory is nonfinite",
    )

    layer_thickness = Float64.(Array(model.geometry.σ_levels_thick))
    size(temperature, 2) == length(layer_thickness) || throw(
        DimensionMismatch(
            "balanced atmosphere inventory and sigma geometry differ",
        ),
    )
    surface_pressure = exp.(log_surface_pressure)
    gravity = Float64(model.planet.gravity)
    heat_capacity = Float64(model.atmosphere.heat_capacity)
    latent_heat = Float64(
        SpeedyWeather.latent_heat_condensation(model.atmosphere),
    )
    layer_mass = surface_pressure ./ gravity .* reshape(
        layer_thickness,
        1,
        :,
    )
    dry_fraction = 1 .- humidity
    dry_mass_by_column = vec(sum(layer_mass .* dry_fraction; dims = 2))
    water_by_column = vec(sum(layer_mass .* humidity; dims = 2))
    sensible_by_column = vec(sum(
        layer_mass .* heat_capacity .* temperature;
        dims = 2,
    ))
    latent_by_column = vec(sum(
        layer_mass .* latent_heat .* humidity;
        dims = 2,
    ))
    kinetic_by_column = vec(sum(
        layer_mass .* 0.5 .* (
            zonal_velocity .^ 2 .+ meridional_velocity .^ 2
        );
        dims = 2,
    ))
    potential_by_column = vec(sum(
        layer_mass .* geopotential;
        dims = 2,
    ))
    point_weights = Float64.(Array(_global_point_weights(model.spectral_grid)))
    length(point_weights) == length(dry_mass_by_column) || throw(
        DimensionMismatch(
            "balanced atmosphere inventory weights do not span the grid",
        ),
    )
    weighted_mean(values) = sum(point_weights .* values)
    dry_mass = weighted_mean(dry_mass_by_column)
    water = weighted_mean(water_by_column)
    sensible = weighted_mean(sensible_by_column)
    latent = weighted_mean(latent_by_column)
    kinetic = weighted_mean(kinetic_by_column)
    potential = weighted_mean(potential_by_column)
    return (
        dry_air_mass_kg_m2_global = dry_mass,
        water_vapor_kg_m2_global = water,
        total_air_mass_kg_m2_global = dry_mass + water,
        dry_sensible_enthalpy_j_m2_global = sensible,
        water_vapor_latent_energy_j_m2_global = latent,
        kinetic_energy_j_m2_global = kinetic,
        geopotential_energy_j_m2_global = potential,
        diagnostic_total_energy_j_m2_global =
            sensible + latent + kinetic + potential,
    )
end

function _balanced_inventory_difference(target, reference)
    keys(target) == keys(reference) || throw(ArgumentError(
        "balanced source inventories do not share the same quantities",
    ))
    names = keys(target)
    return NamedTuple{names}(map(names) do name
        Float64(getproperty(target, name)) -
            Float64(getproperty(reference, name))
    end)
end

"""
Construct a January-style ERA5 analysis increment on an existing atmosphere's
native spectral and RingGrid representations.

The reference is leapfrog slot 2, which is the state used and produced by every
normal SpeedyWeather step. Vorticity and divergence are radius-scaled exactly as
they are in a running simulation before differencing. Grid increments are kept
alongside spectral increments so applying an increment never has to call the
normal transform that rotates current fields into `*_prev` memory.
"""
function _balanced_atmosphere_increment(
    atmosphere::SpeedyWeather.Simulation,
    target_config::ExperimentConfig,
)
    target_config.atmosphere_initial_conditions == :era5_instantaneous ||
        throw(ArgumentError(
            "balanced atmosphere target requires era5_instantaneous",
        ))
    model = atmosphere.model
    variables = atmosphere.variables
    model.geometry.nlayers == target_config.nlayers || throw(ArgumentError(
        "balanced atmosphere target has $(target_config.nlayers) layers but " *
        "the running model has $(model.geometry.nlayers)",
    ))
    model.spectral_grid.trunc == target_config.truncation || throw(ArgumentError(
        "balanced atmosphere target has truncation " *
        "$(target_config.truncation) but the running model has " *
        "$(model.spectral_grid.trunc)",
    ))
    variables.prognostic.scale[] == model.planet.radius || throw(ArgumentError(
        "balanced atmosphere reference must already use dynamical-core " *
        "radius scaling",
    ))

    target = SpeedyWeather.Variables(model)
    analysis = ERA5Atmosphere(target_config)
    SpeedyWeather.initialize!(target, analysis, model)
    SpeedyWeather.scale_prognostic!(target, model.planet.radius)
    SpeedyWeather.SpeedyTransforms.transform!(
        target,
        1,
        model;
        initialize = true,
    )
    SpeedyWeather.temperature_average!(
        target,
        SpeedyWeather.get_step(target.prognostic.temperature, 1),
        model.spectral_transform,
    )
    SpeedyWeather.geopotential!(target, model)

    spectral = NamedTuple{_BALANCED_ATMOSPHERE_SPECTRAL_NAMES}(map(
        _BALANCED_ATMOSPHERE_SPECTRAL_NAMES,
    ) do name
        _balanced_launch_difference(
            SpeedyWeather.get_step(target.prognostic[name], 1),
            SpeedyWeather.get_step(variables.prognostic[name], 2),
            "atmosphere spectral $name",
        )
    end)
    grid = NamedTuple{_BALANCED_ATMOSPHERE_GRID_NAMES}(map(
        _BALANCED_ATMOSPHERE_GRID_NAMES,
    ) do name
        _balanced_launch_difference(
            target.grid[name],
            variables.grid[name],
            "atmosphere grid $name",
        )
    end)
    provenance = (
        source_mode = String(analysis.source_mode),
        source_time = string(analysis.source_time),
        pressure_levels_path = analysis.pressure_levels_path,
        pressure_levels_sha256 = analysis.pressure_levels_sha256,
        single_levels_path = analysis.single_levels_path,
        single_levels_sha256 = analysis.single_levels_sha256,
    )
    source_totals = _balanced_inventory_difference(
        _balanced_atmosphere_inventory(target, model),
        _balanced_atmosphere_inventory(variables, model),
    )
    return BalancedAtmosphereIncrement(
        spectral,
        grid,
        provenance,
        source_totals,
    )
end

@inline function _apply_scaled_balanced_increment!(destination, increment, fraction)
    destination.data .+= fraction .* increment.data
    return destination
end

"""
Apply a fraction of one native atmosphere increment without advancing time.

The same spectral increment is added to both leapfrog slots, preserving their
computational-mode difference. The same physical grid increment is added to
current and previous diagnostic fields. Previous pressure is stored in Pa, so
the log-pressure increment is applied as a multiplicative ratio. Finally the
ordinary conservative humidity repair and the two derived atmosphere fields
used by parameterizations are recomputed.
"""
function _apply_balanced_atmosphere_increment!(
    atmosphere::SpeedyWeather.Simulation,
    increment::BalancedAtmosphereIncrement,
    fraction::Real,
)
    f = Float64(fraction)
    isfinite(f) && 0 <= f <= 1 || throw(ArgumentError(
        "balanced atmosphere increment fraction must lie in [0, 1]",
    ))
    iszero(f) && return nothing

    variables = atmosphere.variables
    model = atmosphere.model
    for name in _BALANCED_ATMOSPHERE_SPECTRAL_NAMES
        prognostic = variables.prognostic[name]
        analysis_increment = increment.spectral[name]
        _apply_scaled_balanced_increment!(
            SpeedyWeather.get_step(prognostic, 1),
            analysis_increment,
            f,
        )
        _apply_scaled_balanced_increment!(
            SpeedyWeather.get_step(prognostic, 2),
            analysis_increment,
            f,
        )
    end
    for name in _BALANCED_ATMOSPHERE_GRID_NAMES
        _apply_scaled_balanced_increment!(
            variables.grid[name],
            increment.grid[name],
            f,
        )
    end
    for name in (:temperature, :humidity, :u, :v)
        _apply_scaled_balanced_increment!(
            variables.grid[Symbol(name, :_prev)],
            increment.grid[name],
            f,
        )
    end
    variables.grid.pressure_prev.data .*= exp.(f .* increment.grid.pressure.data)

    SpeedyWeather.hole_filling!(
        variables.grid.humidity,
        model.hole_filling,
        model,
    )
    SpeedyWeather.hole_filling!(
        variables.grid.humidity_prev,
        model.hole_filling,
        model,
    )
    SpeedyWeather.temperature_average!(
        variables,
        SpeedyWeather.get_step(variables.prognostic.temperature, 2),
        model.spectral_transform,
    )
    SpeedyWeather.geopotential!(variables, model)
    return nothing
end

const _BALANCED_LAND_CONSERVED_NAMES = (
    :internal_energy,
    :saturation_water_ice,
)

const _BALANCED_LAND_PRESCRIBED_NAMES = (
    :vegetation_fraction,
    :leaf_area_index,
    :root_fraction,
)

"""
Native Terrarium analysis increment and its exact imposed-source inventory.

Only fields with a defensible analysis meaning are included. Soil water and
volumetric internal energy are conserved prognostics. ERA5 vegetation cover
and LAI are prescribed boundary data. `surface_excess_water` is deliberately
preserved because ERA5 does not observe ReadyESM's ponded-water reservoir.
Although Terrarium declares `skin_temperature` prognostic, the selected
implicit surface-energy balance overwrites it during every auxiliary update;
it is therefore reconciled diagnostically instead of being presented as a
conserved analysis increment.
"""
struct BalancedLandIncrement{C, B, P, S}
    conserved::C
    prescribed::B
    provenance::P
    source_totals::S
end

function _balanced_launch_interior_difference(target, reference, description)
    target_interior = Oceananigans.interior(target)
    reference_interior = Oceananigans.interior(reference)
    size(target_interior) == size(reference_interior) || throw(DimensionMismatch(
        "$description target size $(size(target_interior)) does not match " *
        "reference size $(size(reference_interior))",
    ))
    eltype(target_interior) == eltype(reference_interior) || throw(ArgumentError(
        "$description target eltype $(eltype(target_interior)) does not match " *
        "reference eltype $(eltype(reference_interior))",
    ))
    difference = similar(target_interior)
    difference .= target_interior .- reference_interior
    return difference
end

function _balanced_land_water_table_elevation(column_grid)
    face_depths = _terrarium_vertical_face_depths(column_grid)
    isempty(face_depths) && error("Terrarium column has no vertical faces")
    target_depth = 5.0
    aligned_depth = face_depths[argmin(abs.(face_depths .- target_depth))]
    aligned_depth >= 0 || error(
        "Terrarium vertical-face depths must be nonnegative",
    )
    return -eltype(column_grid)(aligned_depth)
end

function _balanced_land_target_state(
    atmosphere::SpeedyWeather.Simulation,
    target_config::ExperimentConfig,
)
    target_config.land_model == :terrarium || throw(ArgumentError(
        "balanced land target requires land_model=terrarium",
    ))
    target_config.terrarium_initial_conditions == :era5_instantaneous ||
        throw(ArgumentError(
            "balanced land target requires terrarium_initial_conditions=" *
            "era5_instantaneous",
        ))
    land = atmosphere.model.land
    land isa _ReadyESMAbstractTerrariumLandModel || throw(ArgumentError(
        "balanced land increment requires a running Terrarium land model",
    ))
    target_config.truncation == atmosphere.model.spectral_grid.trunc ||
        throw(ArgumentError(
            "balanced land target truncation does not match the running model",
        ))
    target_config.terrarium_soil_layers ==
        Terrarium.get_field_grid(land.model.grid).Nz || throw(ArgumentError(
            "balanced land target soil-layer count does not match the running model",
        ))
    target_config.terrarium_evapotranspiration ==
        :era5_prescribed_vegetation || throw(ArgumentError(
            "balanced land target requires ERA5-prescribed vegetation",
        ))
    isempty(land.fields) || error(
        "balanced land target allocation does not support aliased prebuilt fields",
    )

    target_analysis = _era5_land_initializers(
        target_config,
        atmosphere.model.spectral_grid,
        land.model.grid,
        _balanced_land_water_table_elevation(land.model.grid),
    )
    reference_state = atmosphere.variables.prognostic.land.terrarium
    target_state = Terrarium.StateVariables(
        land.model;
        clock = Terrarium.Clock(time = reference_state.clock.time),
        boundary_conditions = land.boundary_conditions,
        input_variables = land.input_variables,
    )
    target_integrator = Terrarium.ModelIntegrator(
        target_state.clock,
        land.model,
        Terrarium.InputSources(eltype(reference_state)),
        target_state,
        target_analysis.initializers,
    )
    Terrarium.initialize!(target_integrator)
    target_backend = KernelAbstractions.get_backend(
        Oceananigans.interior(target_state.internal_energy),
    )
    KernelAbstractions.synchronize(target_backend)
    return target_state, target_analysis
end

function _balanced_land_source_totals(atmosphere, conserved)
    model = atmosphere.model
    land = model.land
    layer_thickness, porosity = _terrarium_soil_geometry(
        model,
        model.spectral_grid,
    )
    isnothing(layer_thickness) && error(
        "balanced land source accounting requires Terrarium soil geometry",
    )
    thickness = Float64.(Array(layer_thickness))
    saturation_increment = Float64.(Array(conserved.saturation_water_ice))
    energy_increment = Float64.(Array(conserved.internal_energy))
    size(saturation_increment, 3) == length(thickness) || throw(
        DimensionMismatch("balanced land saturation and soil geometry differ"),
    )
    size(energy_increment) == size(saturation_increment) || throw(
        DimensionMismatch("balanced land water and energy increments differ"),
    )

    point_weights = Float64.(Array(_global_point_weights(model.spectral_grid)))
    land_fraction = clamp.(
        Float64.(Array(model.land_sea_mask.mask.data)),
        0,
        1,
    )
    land_mask = Bool.(Array(land.model.grid.mask.data))
    length(point_weights) == length(land_fraction) == length(land_mask) ||
        throw(DimensionMismatch(
            "balanced land accounting masks do not share the atmosphere grid",
        ))
    column_weights = point_weights[land_mask] .* land_fraction[land_mask]
    land_area_fraction = sum(column_weights)
    land_area_fraction > 0 || error(
        "balanced land source accounting found zero land area",
    )
    column_weights ./= land_area_fraction
    length(column_weights) == size(saturation_increment, 1) || throw(
        DimensionMismatch(
            "balanced land source fields do not match the weighted land columns",
        ),
    )

    thickness_shape = reshape(thickness, 1, 1, :)
    water_by_column = _FRESHWATER_DENSITY_KG_M3 * Float64(porosity) .* vec(sum(
        saturation_increment .* thickness_shape;
        dims = 3,
    ))
    energy_by_column = vec(sum(
        energy_increment .* thickness_shape;
        dims = 3,
    ))
    water_land_mean = sum(column_weights .* water_by_column)
    energy_land_mean = sum(column_weights .* energy_by_column)
    return (
        soil_water_kg_m2_land = water_land_mean,
        soil_internal_energy_j_m2_land = energy_land_mean,
        soil_water_kg_m2_global = land_area_fraction * water_land_mean,
        soil_internal_energy_j_m2_global =
            land_area_fraction * energy_land_mean,
        land_area_fraction,
    )
end

"""
Construct a January-style ERA5 land increment on the running Terrarium grid.

The target internal energy is obtained by initializing a scratch Terrarium
state through the same temperature-to-energy closure as the production model;
it is not approximated with a constant heat capacity. Source totals are area-
weighted inventories of the complete increment and can therefore be scaled by
the restart-safe schedule's applied fraction without weakening physical budget
residuals.
"""
function _balanced_land_increment(
    atmosphere::SpeedyWeather.Simulation,
    target_config::ExperimentConfig,
)
    reference_state = atmosphere.variables.prognostic.land.terrarium
    target_state, analysis = _balanced_land_target_state(
        atmosphere,
        target_config,
    )
    conserved = NamedTuple{_BALANCED_LAND_CONSERVED_NAMES}(map(
        _BALANCED_LAND_CONSERVED_NAMES,
    ) do name
        _balanced_launch_interior_difference(
            getproperty(target_state, name),
            getproperty(reference_state, name),
            "land conserved $name",
        )
    end)
    prescribed_names = Tuple(filter(
        name -> hasproperty(target_state, name) &&
            hasproperty(reference_state, name),
        _BALANCED_LAND_PRESCRIBED_NAMES,
    ))
    prescribed = NamedTuple{prescribed_names}(map(prescribed_names) do name
        _balanced_launch_interior_difference(
            getproperty(target_state, name),
            getproperty(reference_state, name),
            "land prescribed $name",
        )
    end)
    source = analysis.source
    provenance = (
        source_mode = "era5_instantaneous",
        source_time = string(source.source_time),
        source_time_index = source.source_time_index,
        land_state_path = source.path,
        land_state_sha256 = source.sha256,
        soil_water_policy = "ERA5-H-TESSEL-relative-saturation",
        soil_energy_policy =
            "native-Terrarium-temperature-to-internal-energy-closure",
        surface_excess_water_policy = "preserve-unobserved-reservoir",
        skin_temperature_policy =
            "diagnose-with-active-implicit-surface-energy-balance",
        prescribed_fields = join(string.(prescribed_names), ","),
    )
    source_totals = _balanced_land_source_totals(atmosphere, conserved)
    return BalancedLandIncrement(
        conserved,
        prescribed,
        provenance,
        source_totals,
    )
end

@inline function _apply_scaled_balanced_array_increment!(
    destination,
    increment,
    fraction,
)
    destination .+= fraction .* increment
    return destination
end

function _balanced_launch_array_summary(array)
    values = Float64.(Array(array))
    finite = isfinite.(values)
    nonfinite_indices = findall(!, finite)
    finite_indices = findall(finite)
    finite_values = values[finite_indices]
    limits = isempty(finite_values) ? (NaN, NaN) : extrema(finite_values)
    minimum_index = isempty(finite_values) ? nothing :
        finite_indices[argmin(finite_values)]
    maximum_index = isempty(finite_values) ? nothing :
        finite_indices[argmax(finite_values)]
    return (;
        limits,
        minimum_index,
        maximum_index,
        nonfinite_count = length(nonfinite_indices),
        first_nonfinite_index =
            isempty(nonfinite_indices) ? nothing : first(nonfinite_indices),
    )
end

function _balanced_launch_optional_state_summary(state, name::Symbol)
    hasproperty(state, name) || return nothing
    return _balanced_launch_array_summary(
        Oceananigans.interior(getproperty(state, name)),
    )
end

function _balanced_launch_skin_closure_summary(state, land)
    hasproperty(state, :skin_temperature) || return nothing
    hasproperty(state, :ground_heat_flux) || return nothing
    surface_energy_balance = land.model.surface_energy_balance
    hasproperty(surface_energy_balance, :skin_temperature) || return nothing
    skin_scheme = surface_energy_balance.skin_temperature
    hasproperty(skin_scheme, :κₛ) || return nothing

    face_depths = _terrarium_vertical_face_depths(land.model.grid)
    length(face_depths) >= 2 || error(
        "balanced-launch skin closure requires at least one soil layer",
    )
    top_layer_thickness = abs(face_depths[end] - face_depths[end - 1])
    conductivity = Float64(skin_scheme.κₛ)
    conductivity > 0 || error(
        "balanced-launch skin closure requires positive conductivity",
    )

    skin = vec(Float64.(Array(Oceananigans.interior(
        state.skin_temperature,
    ))))
    temperature = Float64.(Array(Oceananigans.interior(state.temperature)))
    ground = vec(@view temperature[:, :, end])
    ground_heat_flux = vec(Float64.(Array(Oceananigans.interior(
        state.ground_heat_flux,
    ))))
    length(skin) == length(ground) == length(ground_heat_flux) || throw(
        DimensionMismatch(
            "balanced-launch skin, ground and heat-flux columns differ",
        ),
    )

    # The implicit surface condition is
    # Ts = Tg - G Δz / (2κ), with every flux positive upwards.
    residual = skin .- ground .+
        ground_heat_flux .* top_layer_thickness ./ (2 * conductivity)
    return _balanced_launch_array_summary(residual)
end

function _trace_balanced_land_application(atmosphere)
    enabled = lowercase(get(
        ENV,
        "READYESM_TRACE_BALANCED_LAND_APPLICATION",
        "false",
    )) in ("1", "true", "yes", "on")
    enabled || return false
    threshold = tryparse(Int, get(
        ENV,
        "READYESM_TRACE_BALANCED_LAND_AFTER_ITERATION",
        "0",
    ))
    isnothing(threshold) && throw(ArgumentError(
        "READYESM_TRACE_BALANCED_LAND_AFTER_ITERATION must be an integer",
    ))
    threshold >= 0 || throw(ArgumentError(
        "READYESM_TRACE_BALANCED_LAND_AFTER_ITERATION must be non-negative",
    ))
    iteration = Int(atmosphere.variables.prognostic.clock.timestep_counter)
    return iteration >= threshold
end

function _trace_balanced_land_application!(atmosphere, phase, fraction)
    _trace_balanced_land_application(atmosphere) || return nothing
    state = atmosphere.variables.prognostic.land.terrarium
    clock = atmosphere.variables.prognostic.clock
    summaries = (
        internal_energy = _balanced_launch_array_summary(
            Oceananigans.interior(state.internal_energy),
        ),
        saturation_water_ice = _balanced_launch_array_summary(
            Oceananigans.interior(state.saturation_water_ice),
        ),
        temperature = _balanced_launch_array_summary(
            Oceananigans.interior(state.temperature),
        ),
        skin_temperature = _balanced_launch_array_summary(
            Oceananigans.interior(state.skin_temperature),
        ),
        surface_net_radiation = _balanced_launch_optional_state_summary(
            state,
            :surface_net_radiation,
        ),
        sensible_heat_flux = _balanced_launch_optional_state_summary(
            state,
            :sensible_heat_flux,
        ),
        latent_heat_flux = _balanced_launch_optional_state_summary(
            state,
            :latent_heat_flux,
        ),
        ground_heat_flux = _balanced_launch_optional_state_summary(
            state,
            :ground_heat_flux,
        ),
        skin_energy_closure_residual =
            _balanced_launch_skin_closure_summary(state, atmosphere.model.land),
        shortwave_down = _balanced_launch_optional_state_summary(
            state,
            :surface_shortwave_down,
        ),
        longwave_down = _balanced_launch_optional_state_summary(
            state,
            :surface_longwave_down,
        ),
        air_temperature = _balanced_launch_optional_state_summary(
            state,
            :air_temperature,
        ),
        specific_humidity = _balanced_launch_optional_state_summary(
            state,
            :specific_humidity,
        ),
        windspeed = _balanced_launch_optional_state_summary(
            state,
            :windspeed,
        ),
        evaporation_ground = _balanced_launch_optional_state_summary(
            state,
            :evaporation_ground,
        ),
        transpiration = _balanced_launch_optional_state_summary(
            state,
            :transpiration,
        ),
        pressure_head = _balanced_launch_array_summary(
            Oceananigans.interior(state.pressure_head),
        ),
        surface_excess_water = _balanced_launch_array_summary(
            Oceananigans.interior(state.surface_excess_water),
        ),
    )
    println(
        "BALANCED_LAND_APPLICATION phase=$phase " *
        "atmosphere_iteration=$(clock.timestep_counter) " *
        "atmosphere_time=$(clock.time) fraction=$(Float64(fraction)) " *
        "summaries=$summaries",
    )
    flush(stdout)
    return summaries
end

function _synchronize_balanced_land_mirrors!(atmosphere)
    variables = atmosphere.variables
    land = atmosphere.model.land
    state = variables.prognostic.land.terrarium
    land_indices = _terrarium_land_indices(land)
    backend = KernelAbstractions.get_backend(land_indices)
    NF = eltype(state)
    _scatter_terrarium_surface!(
        backend,
        variables.prognostic.land.soil_temperature,
        state.skin_temperature,
        land_indices;
        offset = NF(273.15),
    )
    saturation = Oceananigans.interior(state.saturation_water_ice)
    _scatter_terrarium_bottom_kernel!(backend)(
        variables.prognostic.land.soil_moisture,
        saturation,
        land_indices,
        one(NF),
        zero(NF),
        size(saturation, 3);
        ndrange = length(land_indices),
    )
    if haskey(variables.prognostic.land, :sensible_heat_flux)
        _scatter_terrarium_surface!(
            backend,
            variables.prognostic.land.sensible_heat_flux,
            state.sensible_heat_flux,
            land_indices,
        )
    end
    if haskey(variables.prognostic.land, :surface_humidity_flux)
        latent_heat =
            land.model.constants.thermodynamics.latent_heat_vaporization
        _scatter_terrarium_surface!(
            backend,
            variables.prognostic.land.surface_humidity_flux,
            state.latent_heat_flux,
            land_indices;
            multiplier = inv(latent_heat),
        )
    end
    if haskey(variables.parameterizations, :surface_longwave_up)
        _scatter_terrarium_surface!(
            backend,
            variables.parameterizations.surface_longwave_up,
            state.surface_longwave_up,
            land_indices,
        )
    end
    if haskey(variables.parameterizations, :surface_shortwave_up)
        _scatter_terrarium_surface!(
            backend,
            variables.parameterizations.surface_shortwave_up,
            state.surface_shortwave_up,
            land_indices,
        )
    end
    KernelAbstractions.synchronize(backend)
    return nothing
end

function _reconcile_balanced_land_state!(atmosphere; fraction = NaN)
    land = atmosphere.model.land
    state = atmosphere.variables.prognostic.land.terrarium
    backend = KernelAbstractions.get_backend(
        Oceananigans.interior(state.internal_energy),
    )
    # Saturation is the analysed hydrological variable, so use the forward
    # closure (saturation -> pressure head). Internal energy is already the
    # conserved target variable and the same closure diagnoses temperature.
    Terrarium.closure!(state, land.model)
    KernelAbstractions.synchronize(backend)
    _trace_balanced_land_application!(
        atmosphere,
        :after_forward_closure,
        fraction,
    )
    Terrarium.compute_auxiliary!(state, land.model)
    # Terrarium 0.1.6 separates diagnostic reconstruction from boundary
    # updates and no longer provides an aggregate StateVariables halo fill.
    # Follow its normal post-closure ordering so temperature and pressure-head
    # halos (and the corresponding prognostic flux boundaries) are coherent.
    Terrarium.compute_boundary_conditions!(state, land.model)
    KernelAbstractions.synchronize(backend)
    _trace_balanced_land_application!(
        atmosphere,
        :after_auxiliary_reconstruction,
        fraction,
    )
    _synchronize_balanced_land_mirrors!(atmosphere)
    return nothing
end

"""Apply one scheduled fraction of a native land increment without stepping."""
function _apply_balanced_land_increment!(
    atmosphere::SpeedyWeather.Simulation,
    increment::BalancedLandIncrement,
    fraction::Real,
)
    f = Float64(fraction)
    isfinite(f) && 0 <= f <= 1 || throw(ArgumentError(
        "balanced land increment fraction must lie in [0, 1]",
    ))
    iszero(f) && return nothing
    state = atmosphere.variables.prognostic.land.terrarium
    _trace_balanced_land_application!(atmosphere, :before_increment, f)
    for name in _BALANCED_LAND_CONSERVED_NAMES
        _apply_scaled_balanced_array_increment!(
            Oceananigans.interior(getproperty(state, name)),
            increment.conserved[name],
            f,
        )
    end
    for name in keys(increment.prescribed)
        _apply_scaled_balanced_array_increment!(
            Oceananigans.interior(getproperty(state, name)),
            increment.prescribed[name],
            f,
        )
    end
    _trace_balanced_land_application!(atmosphere, :after_additive_increment, f)
    _reconcile_balanced_land_state!(atmosphere; fraction = f)
    return nothing
end

const _BALANCED_OCEAN_FIELD_NAMES = (
    :temperature,
    :salinity,
    :u_velocity,
    :v_velocity,
    :free_surface,
)

"""
Native ECCO ocean and sea-ice analysis increment.

Sea ice is stored as concentration and area-equivalent thickness (`A*h`), not
conditional thickness. This makes every scheduled application conservative in
ice volume even when concentration changes. Ocean currents are the ECCO
geographic components after the same vector-aware rotation used by production
initialisation.
"""
struct BalancedOceanIceIncrement{O, I, P, S}
    ocean::O
    sea_ice::I
    provenance::P
    source_totals::S
end

function _balanced_ocean_fields(ocean)
    return (
        temperature = ocean.tracers.T,
        salinity = ocean.tracers.S,
        u_velocity = ocean.velocities.u,
        v_velocity = ocean.velocities.v,
        free_surface = ocean.free_surface.displacement,
    )
end

function _balanced_metadata_path(metadata)
    path = NumericalEarth.DataWrangling.metadata_path(metadata)
    path isa AbstractVector && return only(path)
    return path
end

function _balanced_ice_maximum_conditional_thickness(model)
    return model.advection isa ConservativeSeaIceAdvection ?
        model.advection.maximum_conditional_thickness : eltype(model.grid)(15)
end

function _mask_balanced_sea_ice_target!(thickness, concentration, grid)
    Nz = size(grid, 3)
    Oceananigans.ImmersedBoundaries.mask_immersed_field_xy!(
        thickness;
        k = Nz,
    )
    Oceananigans.ImmersedBoundaries.mask_immersed_field_xy!(
        concentration;
        k = Nz,
    )
    Oceananigans.fill_halo_regions!((thickness, concentration))
    Oceananigans.Architectures.synchronize(
        Oceananigans.Architectures.architecture(grid),
    )
    return nothing
end

function _balanced_ocean_ice_target_fields(
    earth,
    target_config::ExperimentConfig,
)
    _is_full_state_ecco_initial_conditions(
        target_config.ocean_initial_conditions,
    ) || throw(ArgumentError(
        "balanced ocean target requires full-state ECCO initial conditions",
    ))
    ocean = earth.ocean.model
    sea_ice = earth.sea_ice.model
    size(ocean.grid) == (
        target_config.ocean_nlongitude,
        target_config.ocean_nlatitude,
        target_config.ocean_nlayers,
    ) || throw(ArgumentError(
        "balanced ocean target dimensions do not match the running grid",
    ))

    metadata = _ecco_ocean_ice_metadata(
        ; full_state = true,
        date = _ecco_initial_condition_date(target_config),
        directory = target_config.ecco_initial_conditions_directory,
    )
    foreach(ClimaOcean.download_with_fallback, values(metadata))

    reference = _balanced_ocean_fields(ocean)
    target = NamedTuple{_BALANCED_OCEAN_FIELD_NAMES}(map(
        _BALANCED_OCEAN_FIELD_NAMES,
    ) do name
        similar(reference[name])
    end)
    Oceananigans.set!(target.temperature, metadata.temperature)
    Oceananigans.set!(target.salinity, metadata.salinity)
    _set_extrinsic_ocean_velocity_fields!(
        target.u_velocity,
        target.v_velocity,
        ocean.grid,
        metadata.u_velocity,
        metadata.v_velocity,
    )
    _set_surface_metadatum!(target.free_surface, metadata.free_surface)
    Oceananigans.fill_halo_regions!(Tuple(values(target)))

    # ECCO SIheff is already area-equivalent thickness. Keep an immutable copy
    # before converting the scratch field to ClimaSeaIce's conditional h.
    target_ice_thickness = similar(sea_ice.ice_thickness)
    target_ice_concentration = similar(sea_ice.ice_concentration)
    _set_surface_metadatum!(
        target_ice_thickness,
        metadata.ice_thickness;
        nonfinite_replacement = 0,
    )
    _set_surface_metadatum!(
        target_ice_concentration,
        metadata.ice_concentration;
        nonfinite_replacement = 0,
    )
    _mask_balanced_sea_ice_target!(
        target_ice_thickness,
        target_ice_concentration,
        ocean.grid,
    )
    target_area_equivalent_thickness = copy(
        Oceananigans.interior(target_ice_thickness),
    )
    _initialize_conditional_ice_fields_from_effective_thickness!(
        target_ice_thickness,
        target_ice_concentration,
        ocean.grid,
        _balanced_ice_maximum_conditional_thickness(sea_ice),
    )
    Oceananigans.fill_halo_regions!((
        target_ice_thickness,
        target_ice_concentration,
    ))
    Oceananigans.Architectures.synchronize(
        Oceananigans.Architectures.architecture(ocean.grid),
    )
    return (
        ocean = target,
        sea_ice = (
            concentration = target_ice_concentration,
            conditional_thickness = target_ice_thickness,
            area_equivalent_thickness = target_area_equivalent_thickness,
        ),
        metadata,
    )
end

function _balanced_ocean_cell_volumes(grid)
    operation = _float64_ocean_volume_operation(grid)
    field = Oceananigans.Field(operation)
    Oceananigans.compute!(field)
    Oceananigans.ImmersedBoundaries.mask_immersed_field!(field)
    Oceananigans.Architectures.synchronize(
        Oceananigans.Architectures.architecture(grid),
    )
    return Float64.(Array(Oceananigans.interior(field)))
end

"""
Reconstruct the numerically represented z-star cell volumes for a proposed
free surface without mutating the live grid.

Oceananigans stores the grid metrics in the grid floating-point type. Merely
scaling host `Float64` copies of the current volumes therefore predicts the
continuous `area * Δη` source, but not necessarily the last few ulps of the
volume that the model will actually use. This routine follows the native
`Az * (Δr * σ)` evaluation order and the native `σ = (h + η) / h` arithmetic.
`reference_volumes` supplies the exact three-dimensional wet-cell mask.
"""
function _balanced_zstar_cell_volumes_at_free_surface(
    grid,
    free_surface,
    reference_volumes,
)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    NF = eltype(cpu_grid)
    NF <: AbstractFloat || throw(ArgumentError(
        "balanced-launch ocean grid must use a floating-point element type",
    ))
    Nx, Ny, Nz = size(cpu_grid)
    size(free_surface) == (Nx, Ny) || throw(DimensionMismatch(
        "balanced-launch free surface does not match the ocean grid",
    ))
    size(reference_volumes) == (Nx, Ny, Nz) || throw(DimensionMismatch(
        "balanced-launch reference volumes do not match the ocean grid",
    ))
    target_volumes = zeros(Float64, Nx, Ny, Nz)
    for j in 1:Ny, i in 1:Nx
        any(>(0), @view(reference_volumes[i, j, :])) || continue
        static_depth = Oceananigans.Grids.static_column_depthᶜᶜᵃ(
            i,
            j,
            cpu_grid,
        )
        static_depth > 0 || error(
            "balanced-launch active ocean column has nonpositive static depth",
        )
        target_height = static_depth + NF(free_surface[i, j])
        target_height > 0 || error(
            "balanced-launch target free surface collapses an ocean column",
        )
        stretching = target_height / static_depth
        for k in 1:Nz
            reference_volumes[i, j, k] > 0 || continue
            horizontal_area = Oceananigans.Operators.Azᶜᶜᶜ(
                i,
                j,
                k,
                cpu_grid,
            )
            static_spacing = Oceananigans.Operators.Δrᶜᶜᶜ(
                i,
                j,
                k,
                cpu_grid,
            )
            target_volumes[i, j, k] = Float64(
                horizontal_area * (static_spacing * stretching),
            )
        end
    end
    return target_volumes
end

function _balanced_surface_matrix(field)
    values = Float64.(Array(Oceananigans.interior(field)))
    size(values, 3) == 1 || throw(DimensionMismatch(
        "balanced-launch surface field has more than one vertical level",
    ))
    return reshape(values, size(values, 1), size(values, 2))
end

function _balanced_constant_value(value, description)
    value isa Number && return Float64(value)
    hasproperty(value, :constant) && return Float64(value.constant)
    throw(ArgumentError(
        "balanced-launch $description must be spatially constant",
    ))
end

function _balanced_ocean_ice_inventory_operations(earth)
    ocean = earth.ocean.model
    sea_ice = earth.sea_ice.model
    concentration = sea_ice.ice_concentration
    return (
        ocean_temperature =
            _float64_tracer_integral_operation(ocean.tracers.T),
        ocean_salinity = _float64_tracer_integral_operation(ocean.tracers.S),
        ocean_volume = _float64_ocean_volume_operation(ocean.grid),
        ocean_kinetic_energy_per_density = _budget_integral(
            Oceananigans.Field(0.5 * (
                ocean.velocities.u ^ 2 +
                ocean.velocities.v ^ 2 +
                ocean.velocities.w ^ 2
            )),
        ),
        sea_ice_area = _budget_integral(concentration),
        sea_ice_volume = _budget_integral(Oceananigans.Field(
            sea_ice.ice_thickness * concentration,
        )),
    )
end

"""Absolute conserved ocean/ice inventory on the live z-star geometry."""
function _balanced_ocean_ice_inventory(earth, operations)
    ocean_temperature = _budget_scalar(operations.ocean_temperature)
    ocean_salinity = _budget_scalar(operations.ocean_salinity)
    ocean_volume = _budget_scalar(operations.ocean_volume)
    ocean_kinetic_energy_per_density = _budget_scalar(
        operations.ocean_kinetic_energy_per_density,
    )
    ice_area = _budget_scalar(operations.sea_ice_area)
    ice_volume = _budget_scalar(operations.sea_ice_volume)
    density = Float64(earth.interfaces.ocean_properties.reference_density)
    heat_capacity = Float64(earth.interfaces.ocean_properties.heat_capacity)
    ice_density = _balanced_constant_value(
        earth.sea_ice.model.sea_ice_density,
        "sea-ice density",
    )
    ocean_water_mass = density * ocean_volume
    ocean_heat = density * heat_capacity * ocean_temperature
    ocean_salt = density / 1000 * ocean_salinity
    ocean_kinetic_energy = density * ocean_kinetic_energy_per_density
    ice_mass = ice_density * ice_volume
    return (
        ocean_volume_m3 = ocean_volume,
        ocean_water_mass_kg = ocean_water_mass,
        ocean_sensible_heat_j_relative_to_0c = ocean_heat,
        ocean_salt_mass_kg = ocean_salt,
        ocean_kinetic_energy_j = ocean_kinetic_energy,
        sea_ice_area_m2 = ice_area,
        sea_ice_volume_m3 = ice_volume,
        sea_ice_mass_kg = ice_mass,
        combined_ocean_ice_water_mass_kg = ocean_water_mass + ice_mass,
    )
end


function _balanced_ocean_ice_inventory(earth)
    return _balanced_ocean_ice_inventory(
        earth,
        _balanced_ocean_ice_inventory_operations(earth),
    )
end

"""
Diagnose the complete ocean heat/salt/volume and sea-ice mass source of the
full ECCO increment. The z-star target volume is evaluated analytically from
the target free surface, so changing SSH is included rather than treating the
reference grid-cell volumes as fixed. Ocean kinetic energy is reported
separately by the runtime controller because it must use the reconciled native
staggered velocity state.
"""
function _balanced_ocean_ice_source_totals(earth, target)
    ocean = earth.ocean.model
    sea_ice = earth.sea_ice.model
    reference = _balanced_ocean_fields(ocean)
    volumes = _balanced_ocean_cell_volumes(ocean.grid)
    areas = _ocean_surface_cell_area(ocean.grid)
    active = _ocean_surface_active_mask(ocean.grid)
    size(areas) == size(active) || throw(DimensionMismatch(
        "balanced ocean source area and active masks differ",
    ))

    reference_eta = _balanced_surface_matrix(reference.free_surface)
    target_eta = _balanced_surface_matrix(target.ocean.free_surface)
    column_volume = reshape(
        sum(volumes; dims = 3),
        size(areas, 1),
        size(areas, 2),
    )
    current_height = zeros(Float64, size(areas))
    wet = active .> 0.5
    current_height[wet] .= column_volume[wet] ./ areas[wet]
    all(current_height[wet] .> 0) || error(
        "balanced ocean source accounting found a nonpositive wet-column height",
    )
    volume_scale = ones(Float64, size(areas))
    volume_scale[wet] .= 1 .+
        (target_eta[wet] .- reference_eta[wet]) ./ current_height[wet]
    all(volume_scale[wet] .> 0) || error(
        "balanced ocean target free surface collapses a wet column",
    )
    continuous_target_volumes =
        volumes .* reshape(volume_scale, size(areas)..., 1)
    reconstructed_reference_volumes =
        _balanced_zstar_cell_volumes_at_free_surface(
            ocean.grid,
            reference_eta,
            volumes,
        )
    target_volumes = _balanced_zstar_cell_volumes_at_free_surface(
        ocean.grid,
        target_eta,
        volumes,
    )
    grid_roundoff_scale = eps(eltype(ocean.grid)) * max(
        sum(abs, volumes),
        sum(abs, target_volumes),
        1.0,
    )
    reference_volume_reconstruction_error = sum(abs,
        reconstructed_reference_volumes .- volumes,
    )
    reference_volume_reconstruction_error <= 8grid_roundoff_scale || error(
        "balanced ocean source accounting cannot reconstruct the live z-star volumes",
    )

    reference_temperature = Float64.(Array(
        Oceananigans.interior(reference.temperature),
    ))
    target_temperature = Float64.(Array(
        Oceananigans.interior(target.ocean.temperature),
    ))
    reference_salinity = Float64.(Array(
        Oceananigans.interior(reference.salinity),
    ))
    target_salinity = Float64.(Array(
        Oceananigans.interior(target.ocean.salinity),
    ))
    temperature_content_source = sum(
        target_temperature .* target_volumes .-
        reference_temperature .* volumes,
    )
    salinity_content_source = sum(
        target_salinity .* target_volumes .-
        reference_salinity .* volumes,
    )
    continuous_ocean_volume_source = sum(
        continuous_target_volumes .- volumes,
    )
    ocean_volume_source = sum(target_volumes .- volumes)
    direct_surface_volume_source = sum(
        areas[wet] .* (target_eta[wet] .- reference_eta[wet]),
    )
    abs(continuous_ocean_volume_source - direct_surface_volume_source) <=
        5e-10 * max(abs(direct_surface_volume_source), 1.0) || error(
            "balanced ocean z-star and surface-volume source ledgers disagree",
        )
    ocean_volume_grid_roundoff =
        ocean_volume_source - continuous_ocean_volume_source
    abs(ocean_volume_grid_roundoff) <= 8grid_roundoff_scale || error(
        "balanced ocean target-volume roundoff exceeds the grid-precision envelope",
    )

    reference_concentration = _balanced_surface_matrix(
        sea_ice.ice_concentration,
    )
    reference_thickness = _balanced_surface_matrix(sea_ice.ice_thickness)
    target_concentration = _balanced_surface_matrix(
        target.sea_ice.concentration,
    )
    target_effective = reshape(
        Float64.(Array(target.sea_ice.area_equivalent_thickness)),
        size(areas, 1),
        size(areas, 2),
    )
    reference_effective = reference_concentration .* reference_thickness
    sea_ice_area_source = sum(
        areas[wet] .* (
            target_concentration[wet] .- reference_concentration[wet]
        ),
    )
    sea_ice_volume_source = sum(
        areas[wet] .* (
            target_effective[wet] .- reference_effective[wet]
        ),
    )

    density = Float64(earth.interfaces.ocean_properties.reference_density)
    heat_capacity = Float64(earth.interfaces.ocean_properties.heat_capacity)
    ice_density = _balanced_constant_value(
        sea_ice.sea_ice_density,
        "sea-ice density",
    )
    ocean_heat_source = density * heat_capacity * temperature_content_source
    ocean_salt_mass_source = density / 1000 * salinity_content_source
    ocean_water_mass_source = density * ocean_volume_source
    sea_ice_mass_source = ice_density * sea_ice_volume_source
    return (
        ocean_volume_m3 = ocean_volume_source,
        ocean_water_mass_kg = ocean_water_mass_source,
        ocean_sensible_heat_j_relative_to_0c = ocean_heat_source,
        ocean_salt_mass_kg = ocean_salt_mass_source,
        sea_ice_area_m2 = sea_ice_area_source,
        sea_ice_volume_m3 = sea_ice_volume_source,
        sea_ice_mass_kg = sea_ice_mass_source,
        combined_ocean_ice_water_mass_kg =
            ocean_water_mass_source + sea_ice_mass_source,
        ocean_volume_grid_roundoff_m3 = ocean_volume_grid_roundoff,
        ocean_reference_volume_reconstruction_error_m3 =
            reference_volume_reconstruction_error,
        ocean_kinetic_energy_policy =
            "measure-before-and-after-reconciled-runtime-application",
    )
end

"""Construct an exact January-style ECCO increment on the live coupled grid."""
function _balanced_ocean_ice_increment(
    earth,
    target_config::ExperimentConfig,
)
    target = _balanced_ocean_ice_target_fields(earth, target_config)
    reference = _balanced_ocean_fields(earth.ocean.model)
    ocean_increment = NamedTuple{_BALANCED_OCEAN_FIELD_NAMES}(map(
        _BALANCED_OCEAN_FIELD_NAMES,
    ) do name
        _balanced_launch_interior_difference(
            target.ocean[name],
            reference[name],
            "ocean $name",
        )
    end)
    sea_ice = earth.sea_ice.model
    reference_area_equivalent =
        Oceananigans.interior(sea_ice.ice_concentration) .*
        Oceananigans.interior(sea_ice.ice_thickness)
    ice_increment = (
        concentration = _balanced_launch_interior_difference(
            target.sea_ice.concentration,
            sea_ice.ice_concentration,
            "sea-ice concentration",
        ),
        area_equivalent_thickness =
            target.sea_ice.area_equivalent_thickness .-
            reference_area_equivalent,
    )
    metadata = target.metadata
    path_and_sha(metadata_value) = begin
        path = _balanced_metadata_path(metadata_value)
        (path = path, sha256 = _era5_sha256(path))
    end
    files = NamedTuple{keys(metadata)}(map(path_and_sha, values(metadata)))
    provenance = (
        source_mode = String(target_config.ocean_initial_conditions),
        source_time = string(_ecco_initial_condition_date(target_config)),
        temporal_semantics = "ECCO-V4r4-monthly-mean-labeled-at-month-start",
        files,
        ocean_vector_policy =
            "ECCO-EVEL-NVEL-geographic-to-native-vector-rotation",
        sea_ice_policy =
            "linear-concentration-and-area-equivalent-volume-increment",
        inactive_sea_ice_policy =
            _SEA_ICE_INACTIVE_SURFACE_MASK_SCHEME,
    )
    source_totals = _balanced_ocean_ice_source_totals(earth, target)
    return BalancedOceanIceIncrement(
        ocean_increment,
        ice_increment,
        provenance,
        source_totals,
    )
end

@kernel function _apply_balanced_sea_ice_increment_kernel!(
    ice_thickness,
    concentration,
    concentration_increment,
    effective_thickness_increment,
    fraction,
    maximum_conditional_thickness,
)
    i, j = @index(Global, NTuple)
    k = 1
    @inbounds begin
        A = concentration[i, j, k]
        V = A * ice_thickness[i, j, k]
        target_A = clamp(
            A + fraction * concentration_increment[i, j, k],
            zero(A),
            one(A),
        )
        target_V = max(
            V + fraction * effective_thickness_increment[i, j, k],
            zero(V),
        )
        h, bounded_A = _bounded_conditional_ice_state(
            target_V,
            target_A,
            maximum_conditional_thickness,
        )
        ice_thickness[i, j, k] = h
        concentration[i, j, k] = bounded_A
    end
end

function _apply_balanced_ocean_ice_increment!(
    earth,
    increment::BalancedOceanIceIncrement,
    fraction::Real,
)
    f = Float64(fraction)
    isfinite(f) && 0 <= f <= 1 || throw(ArgumentError(
        "balanced ocean/ice increment fraction must lie in [0, 1]",
    ))
    iszero(f) && return nothing

    ocean = earth.ocean.model
    reference = _balanced_ocean_fields(ocean)
    for name in _BALANCED_OCEAN_FIELD_NAMES
        _apply_scaled_balanced_array_increment!(
            Oceananigans.interior(reference[name]),
            increment.ocean[name],
            f,
        )
    end
    Oceananigans.fill_halo_regions!(Tuple(values(reference)))
    architecture = Oceananigans.Architectures.architecture(ocean.grid)
    Oceananigans.Architectures.synchronize(architecture)
    # Reconcile the free surface, z-star geometry and staggered velocities at
    # explicit completion points. Repeating the reconciliation retains the
    # proven production-initialisation cure for stale tripolar halo geometry.
    for _ in 1:3
        Oceananigans.TimeSteppers.reconcile_state!(ocean)
        Oceananigans.Architectures.synchronize(architecture)
    end
    Oceananigans.TimeSteppers.update_state!(ocean)
    Oceananigans.Architectures.synchronize(architecture)

    sea_ice = earth.sea_ice.model
    launch_grid = ocean.grid isa
        Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ?
        ocean.grid.underlying_grid : ocean.grid
    Oceananigans.Utils.launch!(
        architecture,
        launch_grid,
        :xy,
        _apply_balanced_sea_ice_increment_kernel!,
        sea_ice.ice_thickness,
        sea_ice.ice_concentration,
        increment.sea_ice.concentration,
        increment.sea_ice.area_equivalent_thickness,
        convert(eltype(ocean.grid), f),
        _balanced_ice_maximum_conditional_thickness(sea_ice),
    )
    Oceananigans.Architectures.synchronize(architecture)
    _mask_inactive_sea_ice_surface_cells!(sea_ice)
    Oceananigans.TimeSteppers.update_state!(sea_ice)
    Oceananigans.Architectures.synchronize(architecture)
    return nothing
end

const _BALANCED_ATMOSPHERE_SOURCE_NAMES = (
    :dry_air_mass_kg_m2_global,
    :water_vapor_kg_m2_global,
    :total_air_mass_kg_m2_global,
    :dry_sensible_enthalpy_j_m2_global,
    :water_vapor_latent_energy_j_m2_global,
    :kinetic_energy_j_m2_global,
    :geopotential_energy_j_m2_global,
    :diagnostic_total_energy_j_m2_global,
)

const _BALANCED_LAND_SOURCE_NAMES = (
    :soil_water_kg_m2_land,
    :soil_internal_energy_j_m2_land,
    :soil_water_kg_m2_global,
    :soil_internal_energy_j_m2_global,
)

const _BALANCED_OCEAN_ICE_SOURCE_NAMES = (
    :ocean_volume_m3,
    :ocean_water_mass_kg,
    :ocean_sensible_heat_j_relative_to_0c,
    :ocean_salt_mass_kg,
    :ocean_kinetic_energy_j,
    :sea_ice_area_m2,
    :sea_ice_volume_m3,
    :sea_ice_mass_kg,
    :combined_ocean_ice_water_mass_kg,
)

function _balanced_increment_storage(value, description)
    # RingGrids and similar field wrappers expose their copyable backing array
    # as `.data`. CUDA.CuArray also has a `.data` property, but there it is an
    # opaque GPU allocation handle rather than array storage. Only unwrap a
    # data property that itself implements the array interface.
    data = hasproperty(value, :data) ? getproperty(value, :data) : nothing
    storage = data isa AbstractArray ? data : value
    storage isa AbstractArray || throw(ArgumentError(
        "balanced-launch $description increment has no array storage",
    ))
    isbitstype(eltype(storage)) || throw(ArgumentError(
        "balanced-launch $description increment eltype $(eltype(storage)) " *
        "is not bit-stable",
    ))
    return storage
end

function _balanced_increment_digest_record(description, value)
    storage = _balanced_increment_storage(value, description)
    host = vec(Array(storage))
    digest = bytes2hex(SHA.sha256(reinterpret(UInt8, host)))
    dimensions = isempty(size(storage)) ? "scalar" : join(size(storage), 'x')
    return "$description|$(eltype(storage))|$dimensions|$digest"
end

"""
Return a content digest of every immutable field that defines one launch.

The digest is intentionally computed only for a forecast-referenced launch.
That path cannot be reconstructed from target-data hashes alone: its increment
also depends on the evolved free-forecast endpoint. Hashing each native field
separately bounds host memory use to one field while making restart provenance
sensitive to any changed forecast state.
"""
function _balanced_launch_increment_sha256(
    atmosphere::BalancedAtmosphereIncrement,
    land::BalancedLandIncrement,
    ocean_ice::BalancedOceanIceIncrement,
)
    records = String[]
    for name in _BALANCED_ATMOSPHERE_SPECTRAL_NAMES
        push!(records, _balanced_increment_digest_record(
            "atmosphere.spectral.$name",
            atmosphere.spectral[name],
        ))
    end
    for name in _BALANCED_ATMOSPHERE_GRID_NAMES
        push!(records, _balanced_increment_digest_record(
            "atmosphere.grid.$name",
            atmosphere.grid[name],
        ))
    end
    for name in _BALANCED_LAND_CONSERVED_NAMES
        push!(records, _balanced_increment_digest_record(
            "land.conserved.$name",
            land.conserved[name],
        ))
    end
    for name in keys(land.prescribed)
        push!(records, _balanced_increment_digest_record(
            "land.prescribed.$name",
            land.prescribed[name],
        ))
    end
    for name in _BALANCED_OCEAN_FIELD_NAMES
        push!(records, _balanced_increment_digest_record(
            "ocean.$name",
            ocean_ice.ocean[name],
        ))
    end
    for name in keys(ocean_ice.sea_ice)
        push!(records, _balanced_increment_digest_record(
            "sea_ice.$name",
            ocean_ice.sea_ice[name],
        ))
    end
    return bytes2hex(SHA.sha256(codeunits(join(records, '\n'))))
end

const _BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS = (
    (component = :atmosphere, field = :dry_air_mass_kg_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_dry_air_mass",
     units = "kg/m^2 global"),
    (component = :atmosphere, field = :water_vapor_kg_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_water_vapor",
     units = "kg/m^2 global"),
    (component = :atmosphere, field = :total_air_mass_kg_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_total_air_mass",
     units = "kg/m^2 global"),
    (component = :atmosphere, field = :dry_sensible_enthalpy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_dry_sensible_enthalpy",
     units = "J/m^2 global"),
    (component = :atmosphere, field = :water_vapor_latent_energy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_water_vapor_latent_energy",
     units = "J/m^2 global"),
    (component = :atmosphere, field = :kinetic_energy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_kinetic_energy",
     units = "J/m^2 global"),
    (component = :atmosphere, field = :geopotential_energy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_geopotential_energy",
     units = "J/m^2 global"),
    (component = :atmosphere, field = :diagnostic_total_energy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_atmosphere_diagnostic_total_energy",
     units = "J/m^2 global"),
    (component = :land, field = :soil_water_kg_m2_land,
     netcdf_name = "cumulative_balanced_launch_land_soil_water",
     units = "kg/m^2 of land"),
    (component = :land, field = :soil_internal_energy_j_m2_land,
     netcdf_name = "cumulative_balanced_launch_land_soil_internal_energy",
     units = "J/m^2 of land"),
    (component = :land, field = :soil_water_kg_m2_global,
     netcdf_name = "cumulative_balanced_launch_global_land_soil_water",
     units = "kg/m^2 global"),
    (component = :land, field = :soil_internal_energy_j_m2_global,
     netcdf_name = "cumulative_balanced_launch_global_land_soil_internal_energy",
     units = "J/m^2 global"),
    (component = :ocean_ice, field = :ocean_volume_m3,
     netcdf_name = "cumulative_balanced_launch_ocean_volume",
     units = "m^3"),
    (component = :ocean_ice, field = :ocean_water_mass_kg,
     netcdf_name = "cumulative_balanced_launch_ocean_water_mass",
     units = "kg"),
    (component = :ocean_ice, field = :ocean_sensible_heat_j_relative_to_0c,
     netcdf_name = "cumulative_balanced_launch_ocean_sensible_heat",
     units = "J relative to 0 degree_Celsius"),
    (component = :ocean_ice, field = :ocean_salt_mass_kg,
     netcdf_name = "cumulative_balanced_launch_ocean_salt_mass",
     units = "kg approximate practical-salinity mass"),
    (component = :ocean_ice, field = :ocean_kinetic_energy_j,
     netcdf_name = "cumulative_balanced_launch_ocean_kinetic_energy",
     units = "J"),
    (component = :ocean_ice, field = :sea_ice_area_m2,
     netcdf_name = "cumulative_balanced_launch_sea_ice_area",
     units = "m^2"),
    (component = :ocean_ice, field = :sea_ice_volume_m3,
     netcdf_name = "cumulative_balanced_launch_sea_ice_volume",
     units = "m^3"),
    (component = :ocean_ice, field = :sea_ice_mass_kg,
     netcdf_name = "cumulative_balanced_launch_sea_ice_mass",
     units = "kg"),
    (component = :ocean_ice, field = :combined_ocean_ice_water_mass_kg,
     netcdf_name = "cumulative_balanced_launch_combined_ocean_ice_water_mass",
     units = "kg"),
)

function _balanced_launch_cumulative_source_vector(model)
    diagnostics = balanced_launch_diagnostics(model)
    isnothing(diagnostics) && return zeros(
        Float64,
        length(_BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS),
    )
    return Float64[
        getproperty(
            getproperty(diagnostics, Symbol(spec.component, :_sources)),
            spec.field,
        )
        for spec in _BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS
    ]
end

function _balanced_launch_source_index(component::Symbol, field::Symbol)
    index = findfirst(_BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS) do spec
        spec.component == component && spec.field == field
    end
    isnothing(index) && error(
        "unknown balanced-launch diagnostic source $component.$field",
    )
    return index
end

"""
Restartable owner of one coupled incremental-analysis launch.

The immutable component increments are reconstructed from the declared
background and exact target data whenever a fresh process builds the model.
Only the schedule and the cumulative *actual* imposed sources are restored
from a checkpoint. Provenance equality is mandatory, so a restart cannot pair
checkpoint fractions with a different set of analysis increments.
"""
mutable struct BalancedLaunchController{S, A, L, O, I, P}
    schedule::S
    atmosphere_increment::A
    land_increment::L
    ocean_ice_increment::O
    ocean_ice_inventory_operations::I
    cumulative_atmosphere_sources::Vector{Float64}
    cumulative_land_sources::Vector{Float64}
    cumulative_ocean_ice_sources::Vector{Float64}
    application_count::Int64
    last_application_time_seconds::Float64
    last_application_fraction::Float64
    provenance::P
end

const _BALANCED_LAUNCH_CONTROLLERS = WeakKeyDict{Any, Any}()
const _BALANCED_LAUNCH_CONTROLLERS_LOCK = ReentrantLock()

function _balanced_source_vector(source, names, description)
    values = Float64[getproperty(source, name) for name in names]
    all(isfinite, values) || error(
        "balanced-launch $description source contains a nonfinite value",
    )
    return values
end

function _balanced_source_named_tuple(values, names)
    length(values) == length(names) || throw(DimensionMismatch(
        "balanced-launch source vector and names differ",
    ))
    return NamedTuple{names}(Tuple(Float64.(values)))
end

function _balanced_launch_controller(model)
    return lock(_BALANCED_LAUNCH_CONTROLLERS_LOCK) do
        get(_BALANCED_LAUNCH_CONTROLLERS, model, nothing)
    end
end

function _register_balanced_launch_controller!(model, controller)
    lock(_BALANCED_LAUNCH_CONTROLLERS_LOCK) do
        haskey(_BALANCED_LAUNCH_CONTROLLERS, model) && error(
            "a balanced-launch controller is already installed on this model",
        )
        _BALANCED_LAUNCH_CONTROLLERS[model] = controller
    end
    return controller
end

function _balanced_land_conserved_snapshot(atmosphere)
    state = atmosphere.variables.prognostic.land.terrarium
    return NamedTuple{_BALANCED_LAND_CONSERVED_NAMES}(map(
        _BALANCED_LAND_CONSERVED_NAMES,
    ) do name
        copy(Oceananigans.interior(getproperty(state, name)))
    end)
end

function _balanced_land_conserved_change(atmosphere, before)
    state = atmosphere.variables.prognostic.land.terrarium
    return NamedTuple{_BALANCED_LAND_CONSERVED_NAMES}(map(
        _BALANCED_LAND_CONSERVED_NAMES,
    ) do name
        Oceananigans.interior(getproperty(state, name)) .- before[name]
    end)
end

function _accumulate_balanced_launch_hydrology_sources!(
    atmosphere,
    atmosphere_source,
    land_source,
)
    callbacks = atmosphere.model.callbacks
    haskey(callbacks, :global_atmosphere_diagnostics) || error(
        "balanced launch requires the global atmosphere diagnostics callback",
    )
    diagnostics = callbacks[:global_atmosphere_diagnostics]
    atmosphere_water = Float64(
        atmosphere_source.water_vapor_kg_m2_global,
    )
    land_water = Float64(land_source.soil_water_kg_m2_land)
    isfinite(atmosphere_water) || error(
        "balanced-launch atmosphere water source is nonfinite",
    )
    isfinite(land_water) || error(
        "balanced-launch land water source is nonfinite",
    )
    diagnostics.balanced_launch_atmosphere_water_source += atmosphere_water
    diagnostics.balanced_launch_land_water_source += land_water
    return nothing
end

function _balanced_launch_provenance(
    target_config,
    atmosphere_increment,
    land_increment,
    ocean_ice_increment,
    reference_provenance = nothing,
)
    common = (
        schema = isnothing(reference_provenance) ?
            "readiesm_balanced_launch_v1" :
            "readiesm_forecast_background_balanced_launch_v1",
        target_name = target_config.name,
        target_start_date = string(target_config.start_date),
    )
    reference = if isnothing(reference_provenance)
        NamedTuple()
    else
        reference_provenance isa NamedTuple || throw(ArgumentError(
            "balanced-launch forecast reference provenance must be a NamedTuple",
        ))
        merge(
            reference_provenance,
            (
                increment_sha256 = _balanced_launch_increment_sha256(
                    atmosphere_increment,
                    land_increment,
                    ocean_ice_increment,
                ),
            ),
        )
    end
    components = (
        atmosphere = atmosphere_increment.provenance,
        land = land_increment.provenance,
        ocean_ice = ocean_ice_increment.provenance,
        distribution = "normalized-raised-cosine-incremental-analysis-update",
        source_accounting =
            "immediate-before-after-physical-inventory-at-each-application",
    )
    return isnothing(reference_provenance) ?
        merge(common, components) :
        merge(common, (; reference), components)
end

function _balanced_launch_provenance_fingerprint(provenance)
    hasproperty(provenance, :atmosphere) || return provenance
    atmosphere = provenance.atmosphere
    land = provenance.land
    ocean_ice = provenance.ocean_ice
    ocean_file_hashes = NamedTuple{keys(ocean_ice.files)}(map(
        file -> file.sha256,
        values(ocean_ice.files),
    ))
    fingerprint = (
        schema = provenance.schema,
        target_start_date = provenance.target_start_date,
        atmosphere = (
            source_mode = atmosphere.source_mode,
            source_time = atmosphere.source_time,
            pressure_levels_sha256 = atmosphere.pressure_levels_sha256,
            single_levels_sha256 = atmosphere.single_levels_sha256,
        ),
        land = (
            source_mode = land.source_mode,
            source_time = land.source_time,
            land_state_sha256 = land.land_state_sha256,
            soil_water_policy = land.soil_water_policy,
            soil_energy_policy = land.soil_energy_policy,
            surface_excess_water_policy = land.surface_excess_water_policy,
            skin_temperature_policy = land.skin_temperature_policy,
            prescribed_fields = land.prescribed_fields,
        ),
        ocean_ice = (
            source_mode = ocean_ice.source_mode,
            source_time = ocean_ice.source_time,
            temporal_semantics = ocean_ice.temporal_semantics,
            file_sha256 = ocean_file_hashes,
            ocean_vector_policy = ocean_ice.ocean_vector_policy,
            sea_ice_policy = ocean_ice.sea_ice_policy,
            inactive_sea_ice_policy = ocean_ice.inactive_sea_ice_policy,
        ),
        distribution = provenance.distribution,
        source_accounting = provenance.source_accounting,
    )
    return hasproperty(provenance, :reference) ?
        merge(fingerprint, (reference = provenance.reference,)) :
        fingerprint
end

function _validate_balanced_launch_reference_provenance(model, provenance)
    isnothing(provenance) && return nothing
    provenance isa NamedTuple || throw(ArgumentError(
        "balanced-launch forecast reference provenance must be a NamedTuple",
    ))
    required = (
        :mode,
        :background_name,
        :background_start_date,
        :forecast_end_time_seconds,
        :forecast_end_iteration,
    )
    all(name -> hasproperty(provenance, name), required) || throw(ArgumentError(
        "balanced-launch forecast reference provenance requires " *
        join(string.(required), ", "),
    ))
    provenance.mode == "free-coupled-forecast-endpoint" || throw(ArgumentError(
        "balanced-launch forecast reference mode must be " *
        "free-coupled-forecast-endpoint",
    ))
    end_time = Float64(provenance.forecast_end_time_seconds)
    end_iteration = Int64(provenance.forecast_end_iteration)
    end_time == Float64(model.clock.time) || throw(ArgumentError(
        "balanced-launch reference provenance time $end_time differs from " *
        "the live forecast time $(model.clock.time)",
    ))
    end_iteration == Int64(model.clock.iteration) || throw(ArgumentError(
        "balanced-launch reference provenance iteration $end_iteration differs " *
        "from the live forecast iteration $(model.clock.iteration)",
    ))
    end_time > 0 && end_iteration > 0 || throw(ArgumentError(
        "balanced-launch free-forecast endpoint must be after initialization",
    ))
    return provenance
end

function BalancedLaunchController(
    model,
    target_config::ExperimentConfig;
    start_time_seconds = Float64(model.clock.time),
    window_seconds,
    reference_provenance = nothing,
)
    model isa _ReadyESMRadiativeSpeedyEarthSystem || throw(ArgumentError(
        "balanced launch requires ReadyESM's radiative coupled Earth system",
    ))
    _validate_balanced_launch_reference_provenance(
        model,
        reference_provenance,
    )
    atmosphere_increment = _balanced_atmosphere_increment(
        model.atmosphere,
        target_config,
    )
    land_increment = _balanced_land_increment(
        model.atmosphere,
        target_config,
    )
    ocean_ice_increment = _balanced_ocean_ice_increment(model, target_config)
    ocean_ice_inventory_operations =
        _balanced_ocean_ice_inventory_operations(model)
    schedule = BalancedLaunchSchedule(
        ; start_time_seconds,
        window_seconds,
    )
    provenance = _balanced_launch_provenance(
        target_config,
        atmosphere_increment,
        land_increment,
        ocean_ice_increment,
        reference_provenance,
    )
    return BalancedLaunchController(
        schedule,
        atmosphere_increment,
        land_increment,
        ocean_ice_increment,
        ocean_ice_inventory_operations,
        zeros(Float64, length(_BALANCED_ATMOSPHERE_SOURCE_NAMES)),
        zeros(Float64, length(_BALANCED_LAND_SOURCE_NAMES)),
        zeros(Float64, length(_BALANCED_OCEAN_ICE_SOURCE_NAMES)),
        0,
        NaN,
        0.0,
        provenance,
    )
end

function _apply_balanced_launch_if_due!(model)
    controller = _balanced_launch_controller(model)
    isnothing(controller) && return 0.0
    fraction = _next_balanced_launch_fraction!(
        controller.schedule,
        model.clock.time,
    )
    iszero(fraction) && return 0.0

    atmosphere = model.atmosphere
    atmosphere_before = _balanced_atmosphere_inventory(
        atmosphere.variables,
        atmosphere.model,
    )
    land_before = _balanced_land_conserved_snapshot(atmosphere)
    ocean_ice_before = _balanced_ocean_ice_inventory(
        model,
        controller.ocean_ice_inventory_operations,
    )

    _apply_balanced_atmosphere_increment!(
        atmosphere,
        controller.atmosphere_increment,
        fraction,
    )
    _apply_balanced_land_increment!(
        atmosphere,
        controller.land_increment,
        fraction,
    )
    _apply_balanced_ocean_ice_increment!(
        model,
        controller.ocean_ice_increment,
        fraction,
    )

    atmosphere_source = _balanced_inventory_difference(
        _balanced_atmosphere_inventory(
            atmosphere.variables,
            atmosphere.model,
        ),
        atmosphere_before,
    )
    land_source = _balanced_land_source_totals(
        atmosphere,
        _balanced_land_conserved_change(atmosphere, land_before),
    )
    ocean_ice_source = _balanced_inventory_difference(
        _balanced_ocean_ice_inventory(
            model,
            controller.ocean_ice_inventory_operations,
        ),
        ocean_ice_before,
    )
    _accumulate_balanced_launch_hydrology_sources!(
        atmosphere,
        atmosphere_source,
        land_source,
    )
    controller.cumulative_atmosphere_sources .+= _balanced_source_vector(
        atmosphere_source,
        _BALANCED_ATMOSPHERE_SOURCE_NAMES,
        "atmosphere",
    )
    controller.cumulative_land_sources .+= _balanced_source_vector(
        land_source,
        _BALANCED_LAND_SOURCE_NAMES,
        "land",
    )
    controller.cumulative_ocean_ice_sources .+= _balanced_source_vector(
        ocean_ice_source,
        _BALANCED_OCEAN_ICE_SOURCE_NAMES,
        "ocean/ice",
    )
    controller.application_count += 1
    controller.last_application_time_seconds = Float64(model.clock.time)
    controller.last_application_fraction = fraction
    return fraction
end

function balanced_launch_diagnostics(model)
    controller = _balanced_launch_controller(model)
    isnothing(controller) && return nothing
    return (
        start_time_seconds = controller.schedule.start_time_seconds,
        window_seconds = controller.schedule.window_seconds,
        applied_fraction = controller.schedule.applied_fraction,
        application_count = controller.application_count,
        last_application_time_seconds =
            controller.last_application_time_seconds,
        last_application_fraction = controller.last_application_fraction,
        atmosphere_sources = _balanced_source_named_tuple(
            controller.cumulative_atmosphere_sources,
            _BALANCED_ATMOSPHERE_SOURCE_NAMES,
        ),
        land_sources = _balanced_source_named_tuple(
            controller.cumulative_land_sources,
            _BALANCED_LAND_SOURCE_NAMES,
        ),
        ocean_ice_sources = _balanced_source_named_tuple(
            controller.cumulative_ocean_ice_sources,
            _BALANCED_OCEAN_ICE_SOURCE_NAMES,
        ),
        provenance = controller.provenance,
    )
end

function Oceananigans.prognostic_state(controller::BalancedLaunchController)
    return (
        schedule = Oceananigans.prognostic_state(controller.schedule),
        cumulative_atmosphere_sources =
            copy(controller.cumulative_atmosphere_sources),
        cumulative_land_sources = copy(controller.cumulative_land_sources),
        cumulative_ocean_ice_sources =
            copy(controller.cumulative_ocean_ice_sources),
        application_count = controller.application_count,
        last_application_time_seconds =
            controller.last_application_time_seconds,
        last_application_fraction = controller.last_application_fraction,
        provenance = deepcopy(controller.provenance),
    )
end

function Oceananigans.restore_prognostic_state!(
    controller::BalancedLaunchController,
    state,
)
    _balanced_launch_provenance_fingerprint(state.provenance) ==
        _balanced_launch_provenance_fingerprint(controller.provenance) || error(
        "balanced-launch checkpoint provenance does not match installed increments",
    )
    Oceananigans.restore_prognostic_state!(controller.schedule, state.schedule)
    for (destination, source, description) in (
        (
            controller.cumulative_atmosphere_sources,
            state.cumulative_atmosphere_sources,
            "atmosphere",
        ),
        (
            controller.cumulative_land_sources,
            state.cumulative_land_sources,
            "land",
        ),
        (
            controller.cumulative_ocean_ice_sources,
            state.cumulative_ocean_ice_sources,
            "ocean/ice",
        ),
    )
        length(destination) == length(source) || throw(DimensionMismatch(
            "balanced-launch $description checkpoint source length differs",
        ))
        copyto!(destination, source)
    end
    controller.application_count = Int64(state.application_count)
    controller.last_application_time_seconds =
        Float64(state.last_application_time_seconds)
    controller.last_application_fraction =
        Float64(state.last_application_fraction)
    return controller
end

function _balanced_launch_restart_state(model)
    controller = _balanced_launch_controller(model)
    return isnothing(controller) ? nothing :
        Oceananigans.prognostic_state(controller)
end

function _completed_balanced_launch_controller(state)
    state.schedule.applied_fraction == 1.0 || error(
        "an incomplete balanced-launch checkpoint requires the exact target " *
        "controller to be installed before pickup",
    )
    schedule = BalancedLaunchSchedule(
        ; start_time_seconds = state.schedule.start_time_seconds,
        window_seconds = state.schedule.window_seconds,
        applied_fraction = state.schedule.applied_fraction,
    )
    return BalancedLaunchController(
        schedule,
        nothing,
        nothing,
        nothing,
        nothing,
        Float64.(state.cumulative_atmosphere_sources),
        Float64.(state.cumulative_land_sources),
        Float64.(state.cumulative_ocean_ice_sources),
        Int64(state.application_count),
        Float64(state.last_application_time_seconds),
        Float64(state.last_application_fraction),
        deepcopy(state.provenance),
    )
end

function _restore_balanced_launch_restart_state!(model, state)
    controller = _balanced_launch_controller(model)
    if isnothing(state)
        isnothing(controller) || error(
            "installed balanced launch is absent from the checkpoint",
        )
        return nothing
    end
    if isnothing(controller)
        controller = _completed_balanced_launch_controller(state)
        _register_balanced_launch_controller!(model, controller)
        @info "Restored completed balanced launch without reconstructing target increments" target =
            hasproperty(state.provenance, :target_name) ?
            state.provenance.target_name : "unknown"
    end
    Oceananigans.restore_prognostic_state!(controller, state)
    return controller
end

"""
    prepare_balanced_launch_controller(simulation, target_config;
        window_seconds, start_time_seconds=simulation.model.clock.time,
        reference_provenance=nothing)

Construct but do not attach a balanced-launch controller. Separating
construction from attachment permits the forecast-background workflow to form
an increment at a free-forecast endpoint, restore the exact initial boundary,
and only then attach the immutable increment for replay.
"""
function prepare_balanced_launch_controller(
    simulation::Oceananigans.Simulation,
    target_config::ExperimentConfig;
    window_seconds,
    start_time_seconds = Float64(simulation.model.clock.time),
    reference_provenance = nothing,
)
    return BalancedLaunchController(
        simulation.model,
        target_config;
        start_time_seconds,
        window_seconds,
        reference_provenance,
    )
end

"""Attach a previously prepared controller to its restored replay model."""
function install_balanced_launch!(
    simulation::Oceananigans.Simulation,
    controller::BalancedLaunchController,
)
    controller.schedule.applied_fraction == 0 || error(
        "a prepared balanced-launch controller must be unapplied at installation",
    )
    Float64(simulation.model.clock.time) <=
        controller.schedule.start_time_seconds || error(
        "cannot install a prepared balanced launch after its start time",
    )
    _register_balanced_launch_controller!(simulation.model, controller)
    if iszero(controller.schedule.window_seconds) &&
       Float64(simulation.model.clock.time) >=
       controller.schedule.start_time_seconds
        _apply_balanced_launch_if_due!(simulation.model)
        NumericalEarth.EarthSystemModels.update_state!(simulation.model)
    end
    return controller
end

"""
    install_balanced_launch!(simulation, target_config; window_seconds,
                             start_time_seconds=simulation.model.clock.time)

Attach a conservative atmosphere/land/ocean/ice incremental-analysis launch to
the coupled state graph. A zero-width window applies the complete increment
immediately and reconstructs the ordinary coupled interfaces. Positive windows
are applied after each component sequence and parent-clock tick, immediately
before NumericalEarth rebuilds exchange state and fluxes.
"""
function install_balanced_launch!(
    simulation::Oceananigans.Simulation,
    target_config::ExperimentConfig;
    window_seconds,
    start_time_seconds = Float64(simulation.model.clock.time),
    reference_provenance = nothing,
)
    controller = prepare_balanced_launch_controller(
        simulation,
        target_config;
        window_seconds,
        start_time_seconds,
        reference_provenance,
    )
    return install_balanced_launch!(simulation, controller)
end
