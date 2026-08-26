"""Shortwave placeholder used when RRTMGP handles both spectral ranges globally."""
struct RRTMGPNoShortwave <: SpeedyWeather.AbstractShortwave end

SpeedyWeather.initialize!(::RRTMGPNoShortwave, ::SpeedyWeather.AbstractModel) = nothing
SpeedyWeather.parameterization!(
    ij,
    vars,
    ::RRTMGPNoShortwave,
    model,
) = nothing

"""
Column-local condensate memory used by the experimental prognostic cloud path.

Liquid and ice are stored as grid-box-mean layer water paths in `g m⁻²`.
This makes the precipitation-delay ledger exact without coupling it to changes
in surface pressure. The first experiment intentionally omits condensate
advection; it is a bounded test of prognostic storage, phase and residence time,
not a complete cloud microphysics package.
"""
struct PrognosticCloudCondensateState{NF, A, V}
    liquid_path_gm2::A
    ice_path_gm2::A
    initial_column_path_kgm2::V
    cumulative_retained_kgm2::V
    cumulative_sedimented_kgm2::V
    initialized::V
    retention_fraction::NF
    residence_time_seconds::NF
end

Adapt.@adapt_structure PrognosticCloudCondensateState

function PrognosticCloudCondensateState(
    spectral_grid::SpectralGrid;
    retention_fraction,
    residence_time_hours,
)
    NF = spectral_grid.NF
    architecture = spectral_grid.architecture
    paths = SpeedyWeather.on_architecture(
        architecture,
        zeros(NF, spectral_grid.npoints, spectral_grid.nlayers),
    )
    columns = SpeedyWeather.on_architecture(
        architecture,
        zeros(NF, spectral_grid.npoints),
    )
    initialized = SpeedyWeather.on_architecture(
        architecture,
        zeros(NF, 1),
    )
    ice_paths = similar(paths)
    cumulative_retained = similar(columns)
    cumulative_sedimented = similar(columns)
    fill!(ice_paths, 0)
    fill!(cumulative_retained, 0)
    fill!(cumulative_sedimented, 0)
    return PrognosticCloudCondensateState(
        paths,
        ice_paths,
        columns,
        cumulative_retained,
        cumulative_sedimented,
        initialized,
        NF(retention_fraction),
        NF(3_600 * residence_time_hours),
    )
end

function _prognostic_cloud_condensate_summary(radiation, point_weights)
    if !(radiation isa RRTMGPRadiation) ||
       radiation.cloud_scheme != :prognostic_condensate
        return (
            cloud_condensate_initial_kgm2 = NaN,
            cloud_condensate_endpoint_kgm2 = NaN,
            cloud_condensate_retained_kgm2 = NaN,
            cloud_condensate_sedimented_kgm2 = NaN,
            cloud_condensate_ledger_residual_kgm2 = NaN,
            cloud_condensate_maximum_column_ledger_residual_kgm2 = NaN,
        )
    end
    condensate = radiation.prognostic_cloud_condensate
    weights = Float64.(Array(point_weights))
    initial = Float64.(Array(condensate.initial_column_path_kgm2))
    retained = Float64.(Array(condensate.cumulative_retained_kgm2))
    sedimented = Float64.(Array(condensate.cumulative_sedimented_kgm2))
    endpoint = vec(sum(
        Float64.(Array(
            condensate.liquid_path_gm2 .+ condensate.ice_path_gm2,
        ));
        dims = 2,
    )) ./ 1_000
    residual = endpoint .- initial .- retained .+ sedimented
    return (
        cloud_condensate_initial_kgm2 = sum(weights .* initial),
        cloud_condensate_endpoint_kgm2 = sum(weights .* endpoint),
        cloud_condensate_retained_kgm2 = sum(weights .* retained),
        cloud_condensate_sedimented_kgm2 = sum(weights .* sedimented),
        cloud_condensate_ledger_residual_kgm2 = sum(weights .* residual),
        cloud_condensate_maximum_column_ledger_residual_kgm2 =
            maximum(abs, residual),
    )
end

const _ALL_SKY_CLOUD_SCHEMES = (:diagnostic, :prognostic_condensate)

"""
Clear-sky RRTMGP radiation coupled to SpeedyWeather.

RRTMGP uses `(bottom-to-top layer, column)` arrays while SpeedyWeather uses
`(column, top-to-bottom layer)`. The adapter owns a persistent solver and a
cached heating-rate field, avoiding lookup-table reloads during integration.
"""
mutable struct RRTMGPRadiation{NF, Solver, Heating, Clouds, Condensate} <:
               SpeedyWeather.AbstractLongwave
    solver::Solver
    heating_rate::Heating
    applied_heating_rate::Heating
    co2_ppm::NF
    surface_emissivity::NF
    target_aod_550nm::NF
    solve_every_n_steps::Int
    call_counter::Int
    cloud_scheme::Symbol
    cloud_humidity_search_min_sigma::NF
    clouds::Clouds
    cloud_liquid_water_path_gm2::NF
    cloud_ice_water_path_gm2::NF
    prognostic_cloud_condensate::Condensate
end

const DIAGNOSTIC_CLOUD_HUMIDITY_SEARCH_LEGACY_SCHEME =
    "legacy_layer_index_3_maximum_relative_humidity_v1"
const DIAGNOSTIC_CLOUD_HUMIDITY_SEARCH_SIGMA_SCHEME =
    "speedy_l8_physical_sigma_floor_maximum_relative_humidity_v2"

@inline _diagnostic_cloud_humidity_search_scheme(minimum_sigma) =
    minimum_sigma > 0 ? DIAGNOSTIC_CLOUD_HUMIDITY_SEARCH_SIGMA_SCHEME :
    DIAGNOSTIC_CLOUD_HUMIDITY_SEARCH_LEGACY_SCHEME

"""
GPU-kernel placeholder for RRTMGP's intentionally empty column operation.

ReadyESM launches the complete RRTMGP solve as a global parameterization. When
SpeedyWeather subsequently adapts its fused column-kernel arguments, retaining
the host-side solver would pass managed GPU arrays and other non-isbits state
into that kernel even though the column method is a no-op.
"""
struct RRTMGPColumnNoop <: SpeedyWeather.AbstractLongwave end

# This rule is deliberately limited to CUDA kernel argument conversion. Ordinary
# host/device model construction retains the complete RRTMGP solver used by the
# global parameterization immediately before the fused column kernel.
Adapt.adapt_structure(::CUDA.KernelAdaptor, ::RRTMGPRadiation) = RRTMGPColumnNoop()

const RRTMGPClearSkyRadiation = RRTMGPRadiation

"""Area-weighted global top-of-atmosphere radiation sampled after each model step."""
Base.@kwdef mutable struct GlobalRadiationBudgetCallback{NF, W} <: SpeedyWeather.AbstractCallback
    timestep_counter::Int = 0
    incoming_shortwave::Vector{NF} = NF[]
    outgoing_shortwave::Vector{NF} = NF[]
    outgoing_longwave::Vector{NF} = NF[]
    clear_outgoing_shortwave::Vector{NF} = NF[]
    clear_outgoing_longwave::Vector{NF} = NF[]
    net_downward::Vector{NF} = NF[]
    clear_net_downward::Vector{NF} = NF[]
    point_weights::W
end

function _global_point_weights(spectral_grid::SpectralGrid)
    grid = spectral_grid.grid
    NF = spectral_grid.NF
    weights = zeros(NF, SpeedyWeather.RingGrids.get_npoints(grid))
    quadrature_weights = SpeedyWeather.RingGrids.get_quadrature_weights(grid)
    for (ring_index, ring) in enumerate(SpeedyWeather.RingGrids.eachring(grid))
        weights[ring] .= NF(quadrature_weights[ring_index] / (2 * length(ring)))
    end
    return SpeedyWeather.on_architecture(spectral_grid.architecture, weights)
end

function GlobalRadiationBudgetCallback(spectral_grid::SpectralGrid)
    point_weights = _global_point_weights(spectral_grid)
    NF = spectral_grid.NF
    return GlobalRadiationBudgetCallback{NF, typeof(point_weights)}(; point_weights)
end

_global_grid_mean(values, point_weights) = sum(values .* point_weights)

# A weighted mean of a bounded fraction is mathematically bounded by the same
# interval. Parallel Float32 reduction can nevertheless return the adjacent
# representable value (for example 1.000000119) when many exactly saturated
# columns are accumulated. Keep that reduction noise out of durable physical
# diagnostics without altering the prognostic land state.
@inline function _bounded_fraction_mean(values, point_weights, ::Type{NF}) where {NF}
    mean_value = NF(_global_grid_mean(values, point_weights))
    return clamp(mean_value, zero(NF), one(NF))
end

function _sample_radiation_budget!(callback::GlobalRadiationBudgetCallback, index, model)
    solver = model.longwave_radiation.solver
    weights = callback.point_weights
    incoming_shortwave = _global_grid_mean(RRTMGP.sw_flux_dn(solver)[end, :], weights)
    outgoing_shortwave = _global_grid_mean(RRTMGP.sw_flux_up(solver)[end, :], weights)
    outgoing_longwave = _global_grid_mean(RRTMGP.lw_flux_up(solver)[end, :], weights)
    if isnothing(solver.clear_flux_lw)
        clear_outgoing_shortwave = outgoing_shortwave
        clear_outgoing_longwave = outgoing_longwave
    else
        clear_outgoing_shortwave = _global_grid_mean(
            RRTMGP.clear_sw_flux_up(solver)[end, :],
            weights,
        )
        clear_outgoing_longwave = _global_grid_mean(
            RRTMGP.clear_lw_flux_up(solver)[end, :],
            weights,
        )
    end
    callback.incoming_shortwave[index] = incoming_shortwave
    callback.outgoing_shortwave[index] = outgoing_shortwave
    callback.outgoing_longwave[index] = outgoing_longwave
    callback.clear_outgoing_shortwave[index] = clear_outgoing_shortwave
    callback.clear_outgoing_longwave[index] = clear_outgoing_longwave
    callback.net_downward[index] = incoming_shortwave - outgoing_shortwave - outgoing_longwave
    callback.clear_net_downward[index] =
        incoming_shortwave - clear_outgoing_shortwave - clear_outgoing_longwave
    return nothing
end

function SpeedyWeather.initialize!(
    callback::GlobalRadiationBudgetCallback{NF},
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
) where {NF}
    nsamples = vars.prognostic.clock.n_timesteps + 1
    # Radiation is first evaluated inside the first atmospheric timestep. The
    # initial slot is intentionally NaN rather than a misleading standard-profile
    # flux inherited from RRTMGP construction.
    callback.incoming_shortwave = fill(NF(NaN), nsamples)
    callback.outgoing_shortwave = fill(NF(NaN), nsamples)
    callback.outgoing_longwave = fill(NF(NaN), nsamples)
    callback.clear_outgoing_shortwave = fill(NF(NaN), nsamples)
    callback.clear_outgoing_longwave = fill(NF(NaN), nsamples)
    callback.net_downward = fill(NF(NaN), nsamples)
    callback.clear_net_downward = fill(NF(NaN), nsamples)
    callback.timestep_counter = 1
    return nothing
end

function SpeedyWeather.callback!(
    callback::GlobalRadiationBudgetCallback,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.timestep_counter += 1
    _sample_radiation_budget!(callback, callback.timestep_counter, model)
    return nothing
end

SpeedyWeather.finalize!(::GlobalRadiationBudgetCallback, args...) = nothing

const _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS = (
    :cumulative_convective_humidity_change_profile,
    :cumulative_convective_temperature_change_profile,
    :cumulative_large_scale_humidity_change_profile,
    :cumulative_large_scale_temperature_change_profile,
    :cumulative_radiative_temperature_change_profile,
    :cumulative_surface_humidity_change_profile,
    :cumulative_surface_sensible_temperature_change_profile,
)

const _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS = (
    :weighted_cumulative_convective_humidity_change,
    :weighted_cumulative_convective_temperature_change,
    :weighted_cumulative_large_scale_humidity_change,
    :weighted_cumulative_large_scale_temperature_change,
    :weighted_cumulative_radiative_temperature_change,
    :weighted_cumulative_surface_humidity_change,
    :weighted_cumulative_surface_sensible_temperature_change,
)

"""Daily global atmospheric, process-profile, and hydrological diagnostics."""
Base.@kwdef mutable struct GlobalAtmosphereDiagnosticsCallback{NF, W, Z, A, R} <:
                           SpeedyWeather.AbstractCallback
    timestep_counter::Int = 0
    sample_every_n_steps::Int = 1
    final_timestep::Int = 0
    time_days::Vector{NF} = NF[]
    rainfall_flux::Vector{NF} = NF[]
    convective_rainfall_flux::Vector{NF} = NF[]
    large_scale_rainfall_flux::Vector{NF} = NF[]
    large_scale_snowfall_flux::Vector{NF} = NF[]
    snowfall_flux::Vector{NF} = NF[]
    surface_water_vapor_flux::Vector{NF} = NF[]
    cumulative_surface_water_vapor::Vector{NF} = NF[]
    cumulative_precipitation::Vector{NF} = NF[]
    cumulative_balanced_launch_atmosphere_water_source::Vector{NF} = NF[]
    balanced_launch_atmosphere_water_source::NF = zero(NF)
    water_budget_residual::Vector{NF} = NF[]
    cloud_condensate_water_path::Vector{NF} = NF[]
    total_water_budget_residual::Vector{NF} = NF[]
    surface_air_temperature::Vector{NF} = NF[]
    surface_specific_humidity::Vector{NF} = NF[]
    surface_pressure::Vector{NF} = NF[]
    near_surface_wind_speed::Vector{NF} = NF[]
    column_water_vapor::Vector{NF} = NF[]
    mass_weighted_temperature::Vector{NF} = NF[]
    land_surface_soil_moisture::Vector{NF} = NF[]
    land_total_water_storage::Vector{NF} = NF[]
    land_soil_layer_temperature::Vector{NF} = NF[]
    land_soil_layer_saturation::Vector{NF} = NF[]
    land_rainfall_flux::Vector{NF} = NF[]
    land_snowfall_flux::Vector{NF} = NF[]
    land_evaporation_flux::Vector{NF} = NF[]
    land_ground_evaporation_flux::Vector{NF} = NF[]
    land_transpiration_flux::Vector{NF} = NF[]
    land_surface_runoff_flux::Vector{NF} = NF[]
    land_infiltration_flux::Vector{NF} = NF[]
    cumulative_land_precipitation::Vector{NF} = NF[]
    cumulative_land_evapotranspiration::Vector{NF} = NF[]
    cumulative_land_surface_runoff::Vector{NF} = NF[]
    cumulative_balanced_launch_land_water_source::Vector{NF} = NF[]
    balanced_launch_land_water_source::NF = zero(NF)
    land_water_budget_residual::Vector{NF} = NF[]
    cumulative_convective_humidity_change_profile::Vector{NF} = NF[]
    cumulative_convective_temperature_change_profile::Vector{NF} = NF[]
    cumulative_large_scale_humidity_change_profile::Vector{NF} = NF[]
    cumulative_large_scale_temperature_change_profile::Vector{NF} = NF[]
    cumulative_radiative_temperature_change_profile::Vector{NF} = NF[]
    cumulative_surface_humidity_change_profile::Vector{NF} = NF[]
    cumulative_surface_sensible_temperature_change_profile::Vector{NF} = NF[]
    process_diagnostic_start_timestep::Int = 0
    timestep_seconds::NF
    surface_sigma_thickness::NF
    point_weights::W
    land_point_weights::W
    land_column_weights::W
    column_cumulative_surface_water_vapor::W
    column_cumulative_precipitation::W
    column_cumulative_land_precipitation::W
    column_cumulative_land_evapotranspiration::W
    column_cumulative_land_surface_runoff::W
    soil_layer_thickness::Z
    soil_porosity::NF
    weighted_cumulative_convective_humidity_change::A
    weighted_cumulative_convective_temperature_change::A
    weighted_cumulative_large_scale_humidity_change::A
    weighted_cumulative_large_scale_temperature_change::A
    weighted_cumulative_radiative_temperature_change::A
    weighted_cumulative_surface_humidity_change::A
    weighted_cumulative_surface_sensible_temperature_change::A
    radiation::R
end

KernelAbstractions.@kernel function _accumulate_atmosphere_process_profiles_kernel!(
    weighted_cumulative_convective_humidity_change,
    weighted_cumulative_convective_temperature_change,
    weighted_cumulative_large_scale_humidity_change,
    weighted_cumulative_large_scale_temperature_change,
    weighted_cumulative_radiative_temperature_change,
    weighted_cumulative_surface_humidity_change,
    weighted_cumulative_surface_sensible_temperature_change,
    convective_humidity_tendency,
    convective_temperature_tendency,
    large_scale_humidity_tendency,
    large_scale_temperature_tendency,
    radiative_heating_rate,
    surface_humidity_flux,
    surface_sensible_heat_flux,
    surface_pressure,
    point_weights,
    timestep_seconds,
    gravity,
    heat_capacity,
    surface_sigma_thickness,
)
    point, layer = @index(Global, NTuple)
    nlayers = size(convective_humidity_tendency, 2)
    weighted_timestep = point_weights[point] * timestep_seconds

    weighted_cumulative_convective_humidity_change[point, layer] +=
        weighted_timestep * convective_humidity_tendency[point, layer]
    weighted_cumulative_convective_temperature_change[point, layer] +=
        weighted_timestep * convective_temperature_tendency[point, layer]
    weighted_cumulative_large_scale_humidity_change[point, layer] +=
        weighted_timestep * large_scale_humidity_tendency[point, layer]
    weighted_cumulative_large_scale_temperature_change[point, layer] +=
        weighted_timestep * large_scale_temperature_tendency[point, layer]
    # RRTMGP stores layers bottom-to-top, opposite to SpeedyWeather.
    weighted_cumulative_radiative_temperature_change[point, layer] +=
        weighted_timestep * radiative_heating_rate[nlayers - layer + 1, point]

    if layer == nlayers
        inverse_surface_air_mass = gravity /
            (surface_pressure[point] * surface_sigma_thickness)
        weighted_cumulative_surface_humidity_change[point, layer] +=
            weighted_timestep * inverse_surface_air_mass *
            surface_humidity_flux[point]
        weighted_cumulative_surface_sensible_temperature_change[point, layer] +=
            weighted_timestep * inverse_surface_air_mass *
            surface_sensible_heat_flux[point] / heat_capacity
    end
end

function _accumulate_atmosphere_process_profiles!(callback, vars, model)
    convective_humidity_tendency =
        vars.parameterizations.convective_humidity_tendency.data
    backend = KernelAbstractions.get_backend(convective_humidity_tendency)
    _accumulate_atmosphere_process_profiles_kernel!(backend)(
        callback.weighted_cumulative_convective_humidity_change,
        callback.weighted_cumulative_convective_temperature_change,
        callback.weighted_cumulative_large_scale_humidity_change,
        callback.weighted_cumulative_large_scale_temperature_change,
        callback.weighted_cumulative_radiative_temperature_change,
        callback.weighted_cumulative_surface_humidity_change,
        callback.weighted_cumulative_surface_sensible_temperature_change,
        convective_humidity_tendency,
        vars.parameterizations.convective_temperature_tendency.data,
        vars.parameterizations.large_scale_condensation_humidity_tendency.data,
        vars.parameterizations.large_scale_condensation_temperature_tendency.data,
        callback.radiation.applied_heating_rate,
        vars.parameterizations.surface_humidity_flux.data,
        vars.parameterizations.sensible_heat_flux.data,
        vars.grid.pressure_prev.data,
        callback.point_weights,
        model.time_stepping.Δt_sec,
        model.planet.gravity,
        model.atmosphere.heat_capacity,
        callback.surface_sigma_thickness;
        ndrange = size(convective_humidity_tendency),
    )
    return nothing
end

function _append_atmosphere_process_profiles!(callback)
    backend = KernelAbstractions.get_backend(
        callback.weighted_cumulative_convective_humidity_change,
    )
    KernelAbstractions.synchronize(backend)
    for (history_name, accumulator_name) in zip(
        _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS,
        _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS,
    )
        profile = vec(Array(sum(getproperty(callback, accumulator_name); dims = 1)))
        append!(getproperty(callback, history_name), profile)
    end
    return nothing
end

function _homogeneous_soil_porosity_scheme(stratigraphy)
    if hasproperty(stratigraphy, :porosity)
        # Terrarium 0.1.2, pinned by ReadyESM's SpeedyWeather extension.
        return stratigraphy.porosity
    elseif hasproperty(stratigraphy, :horizons)
        # Forward-compatible with Terrarium's newer single-horizon layout.
        horizons = stratigraphy.horizons
        length(horizons) == 1 || error(
            "land-water storage diagnostic requires one homogeneous soil horizon",
        )
        return only(horizons).porosity
    end
    error("unsupported Terrarium stratigraphy for land-water storage")
end

function _terrarium_soil_geometry(model, spectral_grid)
    land = model.land
    if !hasproperty(land, :model) || !(land.model isa Terrarium.AbstractModel)
        return nothing, spectral_grid.NF(NaN)
    end

    soil_model = land.model
    field_grid = Terrarium.get_field_grid(soil_model.grid)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), field_grid)
    NF = spectral_grid.NF
    layer_thickness = NF[
        Oceananigans.Δzᵃᵃᶜ(1, 1, k, cpu_grid) for k in 1:cpu_grid.Nz
    ]
    layer_thickness = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        layer_thickness,
    )

    # ReadyESM currently constructs one homogeneous mineral horizon and the
    # default zero-carbon biogeochemistry. In that exact configuration the
    # total pore-water depth is porosity * sum(saturation * layer thickness).
    stratigraphy = soil_model.soil.strat
    porosity_scheme = _homogeneous_soil_porosity_scheme(stratigraphy)
    hasproperty(porosity_scheme, :mineral_porosity) || error(
        "land-water storage diagnostic requires a constant mineral porosity",
    )
    hasproperty(soil_model.soil.biogeochem, :ρ_soc) &&
        !iszero(soil_model.soil.biogeochem.ρ_soc) && error(
            "land-water storage diagnostic does not yet support organic soil carbon",
        )
    return layer_thickness, NF(porosity_scheme.mineral_porosity)
end

function _column_soil_water_storage(
    saturation,
    surface_excess_water,
    layer_thickness,
    porosity,
)
    size(saturation, 3) == length(layer_thickness) || throw(
        DimensionMismatch("soil saturation and vertical spacing differ"),
    )
    pore_water_depth = vec(sum(
        saturation .* reshape(layer_thickness, 1, 1, :);
        dims = 3,
    )) .* porosity
    return _FRESHWATER_DENSITY_KG_M3 .* (
        pore_water_depth .+ vec(surface_excess_water)
    )
end

function GlobalAtmosphereDiagnosticsCallback(model::SpeedyWeather.AbstractModel)
    spectral_grid = model.spectral_grid
    NF = spectral_grid.NF
    weights = _global_point_weights(spectral_grid)
    # Geometry arrays live on the selected architecture. Cache the one scalar
    # needed by the per-step surface-process kernel while constructing the
    # callback; host-side `[end]` on a CuArray is forbidden at runtime.
    surface_sigma_thickness = NF(last(Array(model.geometry.σ_levels_thick)))
    surface_sigma_thickness > 0 || error(
        "global atmosphere diagnostics found a nonpositive surface sigma thickness",
    )
    host_weights = Array(weights)
    host_land_fraction = clamp.(Array(model.land_sea_mask.mask.data), 0, 1)
    host_land_mask = host_land_fraction .> 0
    land_point_weights = host_weights .* host_land_fraction
    land_area_fraction = sum(land_point_weights)
    land_area_fraction > 0 || error("global atmosphere diagnostics found no land points")
    land_point_weights ./= land_area_fraction
    land_column_weights = (
        host_weights[host_land_mask] .* host_land_fraction[host_land_mask]
    )
    land_column_weights ./= sum(land_column_weights)
    land_point_weights = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        land_point_weights,
    )
    land_column_weights = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        land_column_weights,
    )
    soil_layer_thickness, soil_porosity = _terrarium_soil_geometry(
        model,
        spectral_grid,
    )
    column_cumulative_surface_water_vapor = similar(weights)
    column_cumulative_precipitation = similar(weights)
    column_cumulative_land_precipitation = similar(land_column_weights)
    column_cumulative_land_evapotranspiration = similar(land_column_weights)
    column_cumulative_land_surface_runoff = similar(land_column_weights)
    fill!(column_cumulative_surface_water_vapor, 0)
    fill!(column_cumulative_precipitation, 0)
    fill!(column_cumulative_land_precipitation, 0)
    fill!(column_cumulative_land_evapotranspiration, 0)
    fill!(column_cumulative_land_surface_runoff, 0)
    process_accumulator = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(NF, spectral_grid.npoints, spectral_grid.nlayers),
    )
    process_accumulators = map(
        _ -> similar(process_accumulator),
        _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS,
    )
    for accumulator in process_accumulators
        fill!(accumulator, 0)
    end
    return GlobalAtmosphereDiagnosticsCallback{
        NF,
        typeof(weights),
        typeof(soil_layer_thickness),
        typeof(process_accumulator),
        typeof(model.longwave_radiation),
    }(
        ; point_weights = weights,
        timestep_seconds = NF(model.time_stepping.Δt_sec),
        surface_sigma_thickness,
        land_point_weights,
        land_column_weights,
        column_cumulative_surface_water_vapor,
        column_cumulative_precipitation,
        column_cumulative_land_precipitation,
        column_cumulative_land_evapotranspiration,
        column_cumulative_land_surface_runoff,
        soil_layer_thickness,
        soil_porosity,
        NamedTuple{_ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS}(
            process_accumulators,
        )...,
        radiation = model.longwave_radiation,
    )
end

function _accumulate_atmosphere_water_columns!(
    cumulative_surface_water_vapor,
    cumulative_precipitation,
    surface_water_vapor_flux,
    rain_rate,
    snow_rate,
    timestep_seconds,
)
    cumulative_surface_water_vapor .+=
        timestep_seconds .* surface_water_vapor_flux
    cumulative_precipitation .+=
        timestep_seconds .* _FRESHWATER_DENSITY_KG_M3 .* (rain_rate .+ snow_rate)
    return nothing
end

@inline _atmosphere_water_budget_residual(
    column_water_vapor,
    initial_column_water_vapor,
    cumulative_surface_water_vapor,
    cumulative_precipitation,
) = column_water_vapor - initial_column_water_vapor -
    cumulative_surface_water_vapor + cumulative_precipitation

function _sample_global_atmosphere!(callback, vars, model; precipitation = true)
    NF = eltype(callback.time_days)
    weights = callback.point_weights
    nlayers = model.geometry.nlayers
    Δσ = model.geometry.σ_levels_thick
    gravity = model.planet.gravity

    T = vars.grid.temperature_prev.data
    q = vars.grid.humidity_prev.data
    u = vars.grid.u_prev.data
    v = vars.grid.v_prev.data
    ps = vars.grid.pressure_prev.data
    Ts = vars.parameterizations.surface_air_temperature.data
    rain = vars.parameterizations.rain_rate.data
    convective_rain = vars.parameterizations.rain_rate_convection.data
    large_scale_rain = vars.parameterizations.rain_rate_large_scale.data
    large_scale_snow = vars.parameterizations.snow_rate_large_scale.data
    snow = vars.parameterizations.snow_rate.data
    land_soil_moisture = vars.prognostic.land.soil_moisture.data
    # Diagnose the conditional land flux, not the atmosphere's combined
    # land-plus-ocean grid-cell flux. The latter contaminates land means at
    # fractional coastal cells and cannot close a Terrarium-only ET split.
    land_surface_water_vapor_flux = if haskey(
        vars.prognostic.land,
        :surface_humidity_flux,
    )
        vars.prognostic.land.surface_humidity_flux.data
    elseif haskey(vars.parameterizations.land, :surface_humidity_flux)
        vars.parameterizations.land.surface_humidity_flux.data
    else
        vars.parameterizations.surface_humidity_flux.data
    end
    surface_water_vapor_flux = vars.parameterizations.surface_humidity_flux.data
    has_terrarium = haskey(vars.prognostic.land, :terrarium)

    q_surface = view(q, :, nlayers)
    u_surface = view(u, :, nlayers)
    v_surface = view(v, :, nlayers)
    wind_speed = sqrt.(u_surface .^ 2 .+ v_surface .^ 2)
    column_q = vec(sum(q .* reshape(Δσ, 1, :); dims = 2)) .* ps ./ gravity
    column_T = vec(sum(T .* reshape(Δσ, 1, :); dims = 2))

    push!(callback.time_days, NF(callback.timestep_counter * model.time_stepping.Δt_sec / 86_400))
    _append_atmosphere_process_profiles!(callback)
    push!(callback.rainfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 * _global_grid_mean(rain, weights)) : NF(NaN))
    push!(callback.convective_rainfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 *
                             _global_grid_mean(convective_rain, weights)) : NF(NaN))
    push!(callback.large_scale_rainfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 *
                             _global_grid_mean(large_scale_rain, weights)) : NF(NaN))
    push!(callback.large_scale_snowfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 *
                             _global_grid_mean(large_scale_snow, weights)) : NF(NaN))
    push!(callback.snowfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 * _global_grid_mean(snow, weights)) : NF(NaN))
    push!(callback.surface_water_vapor_flux,
          NF(_global_grid_mean(surface_water_vapor_flux, weights)))
    # `surface_air_temperature` is a parameterization work array and has not
    # been diagnosed when this callback takes its time-zero sample.  Its
    # allocation value is zero, which is neither a physical initial condition
    # nor the lowest atmospheric model-level temperature.  Preserve that
    # distinction explicitly and start this series after the first completed
    # atmospheric step.
    push!(callback.surface_air_temperature,
          precipitation ? NF(_global_grid_mean(Ts, weights)) : NF(NaN))
    push!(callback.surface_specific_humidity, NF(_global_grid_mean(q_surface, weights)))
    push!(callback.surface_pressure, NF(_global_grid_mean(ps, weights)))
    push!(callback.near_surface_wind_speed, NF(_global_grid_mean(wind_speed, weights)))
    push!(callback.column_water_vapor, NF(_global_grid_mean(column_q, weights)))
    cumulative_surface_water_vapor = NF(_global_grid_mean(
        callback.column_cumulative_surface_water_vapor,
        weights,
    ))
    cumulative_precipitation = NF(_global_grid_mean(
        callback.column_cumulative_precipitation,
        weights,
    ))
    push!(callback.cumulative_surface_water_vapor, cumulative_surface_water_vapor)
    push!(callback.cumulative_precipitation, cumulative_precipitation)
    push!(
        callback.cumulative_balanced_launch_atmosphere_water_source,
        callback.balanced_launch_atmosphere_water_source,
    )
    push!(callback.water_budget_residual, NF(_atmosphere_water_budget_residual(
        last(callback.column_water_vapor),
        first(callback.column_water_vapor),
        cumulative_surface_water_vapor,
        cumulative_precipitation,
    ) - callback.balanced_launch_atmosphere_water_source))
    cloud_summary = _prognostic_cloud_condensate_summary(
        callback.radiation,
        weights,
    )
    cloud_storage = NF(cloud_summary.cloud_condensate_endpoint_kgm2)
    initial_cloud_storage = NF(cloud_summary.cloud_condensate_initial_kgm2)
    if !isfinite(cloud_storage)
        cloud_storage = zero(NF)
        initial_cloud_storage = zero(NF)
    end
    push!(callback.cloud_condensate_water_path, cloud_storage)
    push!(
        callback.total_water_budget_residual,
        last(callback.water_budget_residual) +
            cloud_storage - initial_cloud_storage,
    )
    push!(callback.mass_weighted_temperature, NF(_global_grid_mean(column_T, weights)))
    push!(callback.land_surface_soil_moisture,
          NF(_global_grid_mean(land_soil_moisture, callback.land_point_weights)))
    push!(callback.land_rainfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 *
                             _global_grid_mean(rain, callback.land_point_weights)) : NF(NaN))
    push!(callback.land_snowfall_flux,
          precipitation ? NF(_FRESHWATER_DENSITY_KG_M3 *
                             _global_grid_mean(snow, callback.land_point_weights)) : NF(NaN))
    push!(callback.land_evaporation_flux,
          NF(_global_grid_mean(
              land_surface_water_vapor_flux,
              callback.land_point_weights,
          )))
    if has_terrarium
        terrarium_state = vars.prognostic.land.terrarium
        soil_temperature = Oceananigans.interior(terrarium_state.temperature)
        saturation = Oceananigans.interior(terrarium_state.saturation_water_ice)
        n_soil_layers = length(callback.soil_layer_thickness)
        size(soil_temperature, 3) == n_soil_layers || error(
            "Terrarium temperature and diagnostic soil geometry differ",
        )
        size(saturation, 3) == n_soil_layers || error(
            "Terrarium saturation and diagnostic soil geometry differ",
        )
        for layer in 1:n_soil_layers
            push!(
                callback.land_soil_layer_temperature,
                NF(_global_grid_mean(
                    view(soil_temperature, :, 1, layer),
                    callback.land_column_weights,
                )),
            )
            push!(
                callback.land_soil_layer_saturation,
                _bounded_fraction_mean(
                    view(saturation, :, 1, layer),
                    callback.land_column_weights,
                    NF,
                ),
            )
        end
        surface_excess_water = Oceananigans.interior(
            terrarium_state.surface_excess_water,
        )
        column_water_storage = _column_soil_water_storage(
            saturation,
            surface_excess_water,
            callback.soil_layer_thickness,
            callback.soil_porosity,
        )
        land_runoff = vec(Oceananigans.interior(terrarium_state.surface_runoff))
        land_infiltration = vec(Oceananigans.interior(terrarium_state.infiltration))
        if hasproperty(terrarium_state, :ground_water_flux) &&
           hasproperty(terrarium_state, :transpiration_water_flux)
            ground_evaporation = _FRESHWATER_DENSITY_KG_M3 .* vec(
                Oceananigans.interior(terrarium_state.ground_water_flux),
            )
            transpiration = _FRESHWATER_DENSITY_KG_M3 .* vec(
                Oceananigans.interior(terrarium_state.transpiration_water_flux),
            )
            push!(callback.land_ground_evaporation_flux,
                  NF(_global_grid_mean(
                      ground_evaporation,
                      callback.land_column_weights,
                  )))
            push!(callback.land_transpiration_flux,
                  NF(_global_grid_mean(
                      transpiration,
                      callback.land_column_weights,
                  )))
        else
            # The upstream bare-ground scheme has no separate liquid-water
            # diagnostic. Its entire atmospheric moisture flux is ground
            # evaporation and its transpiration contribution is exactly zero.
            push!(callback.land_ground_evaporation_flux,
                  last(callback.land_evaporation_flux))
            push!(callback.land_transpiration_flux, zero(NF))
        end
        push!(callback.land_total_water_storage,
              NF(_global_grid_mean(
                  column_water_storage,
                  callback.land_column_weights,
              )))
        push!(callback.land_surface_runoff_flux,
              NF(_FRESHWATER_DENSITY_KG_M3 *
                 _global_grid_mean(land_runoff, callback.land_column_weights)))
        push!(callback.land_infiltration_flux,
              NF(_FRESHWATER_DENSITY_KG_M3 *
                 _global_grid_mean(land_infiltration, callback.land_column_weights)))
        cumulative_land_precipitation = NF(_global_grid_mean(
            callback.column_cumulative_land_precipitation,
            callback.land_column_weights,
        ))
        cumulative_land_evapotranspiration = NF(_global_grid_mean(
            callback.column_cumulative_land_evapotranspiration,
            callback.land_column_weights,
        ))
        cumulative_land_surface_runoff = NF(_global_grid_mean(
            callback.column_cumulative_land_surface_runoff,
            callback.land_column_weights,
        ))
        push!(
            callback.cumulative_land_precipitation,
            cumulative_land_precipitation,
        )
        push!(
            callback.cumulative_land_evapotranspiration,
            cumulative_land_evapotranspiration,
        )
        push!(
            callback.cumulative_land_surface_runoff,
            cumulative_land_surface_runoff,
        )
        push!(
            callback.cumulative_balanced_launch_land_water_source,
            callback.balanced_launch_land_water_source,
        )
        push!(callback.land_water_budget_residual, NF(
            last(callback.land_total_water_storage) -
            first(callback.land_total_water_storage) -
            cumulative_land_precipitation +
            cumulative_land_evapotranspiration +
            cumulative_land_surface_runoff -
            callback.balanced_launch_land_water_source
        ))
    else
        push!(callback.land_ground_evaporation_flux,
              last(callback.land_evaporation_flux))
        push!(callback.land_transpiration_flux, zero(NF))
        push!(callback.land_total_water_storage, NF(NaN))
        push!(callback.land_surface_runoff_flux, NF(NaN))
        push!(callback.land_infiltration_flux, NF(NaN))
        push!(callback.cumulative_land_precipitation, NF(NaN))
        push!(callback.cumulative_land_evapotranspiration, NF(NaN))
        push!(callback.cumulative_land_surface_runoff, NF(NaN))
        push!(callback.cumulative_balanced_launch_land_water_source, NF(NaN))
        push!(callback.land_water_budget_residual, NF(NaN))
    end
    return nothing
end

function SpeedyWeather.initialize!(
    callback::GlobalAtmosphereDiagnosticsCallback{NF},
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
) where {NF}
    _initialize_terrarium_temperature_from_atmosphere!(vars, model)
    # SpeedyWeather may adjust Δt while initializing a simulation so an
    # integer number of steps lands exactly on its output cadence. Refresh the
    # constructor value here; legacy-restart timestamp classification must use
    # the timestep that actually generated the retained histories.
    callback.timestep_seconds = NF(model.time_stepping.Δt_sec)
    callback.timestep_counter = 0
    callback.final_timestep = vars.prognostic.clock.n_timesteps
    callback.sample_every_n_steps = max(
        1,
        round(Int, 86_400 / model.time_stepping.Δt_sec),
    )
    for name in (
        :time_days,
        :rainfall_flux,
        :convective_rainfall_flux,
        :large_scale_rainfall_flux,
        :large_scale_snowfall_flux,
        :snowfall_flux,
        :surface_water_vapor_flux,
        :cumulative_surface_water_vapor,
        :cumulative_precipitation,
        :cumulative_balanced_launch_atmosphere_water_source,
        :water_budget_residual,
        :surface_air_temperature,
        :surface_specific_humidity,
        :surface_pressure,
        :near_surface_wind_speed,
        :column_water_vapor,
        :mass_weighted_temperature,
        :land_surface_soil_moisture,
        :land_total_water_storage,
        :land_soil_layer_temperature,
        :land_soil_layer_saturation,
        :land_rainfall_flux,
        :land_snowfall_flux,
        :land_evaporation_flux,
        :land_ground_evaporation_flux,
        :land_transpiration_flux,
        :land_surface_runoff_flux,
        :land_infiltration_flux,
        :cumulative_land_precipitation,
        :cumulative_land_evapotranspiration,
        :cumulative_land_surface_runoff,
        :cumulative_balanced_launch_land_water_source,
        :land_water_budget_residual,
        _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS...,
    )
        empty!(getproperty(callback, name))
    end
    fill!(callback.column_cumulative_surface_water_vapor, 0)
    fill!(callback.column_cumulative_precipitation, 0)
    fill!(callback.column_cumulative_land_precipitation, 0)
    fill!(callback.column_cumulative_land_evapotranspiration, 0)
    fill!(callback.column_cumulative_land_surface_runoff, 0)
    for name in _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS
        fill!(getproperty(callback, name), 0)
    end
    callback.process_diagnostic_start_timestep = 0
    callback.balanced_launch_atmosphere_water_source = zero(NF)
    callback.balanced_launch_land_water_source = zero(NF)
    _sample_global_atmosphere!(callback, vars, model; precipitation = false)
    return nothing
end

function SpeedyWeather.callback!(
    callback::GlobalAtmosphereDiagnosticsCallback,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.timestep_counter += 1
    # A checkpoint written before this observer may fall between its retained
    # daily samples. Such a restart begins at the next retained diagnostic
    # boundary, where an honest zero can be represented, rather than assigning
    # a partial interval to an earlier timestamp.
    callback.timestep_counter > callback.process_diagnostic_start_timestep &&
        _accumulate_atmosphere_process_profiles!(callback, vars, model)
    _accumulate_atmosphere_water_columns!(
        callback.column_cumulative_surface_water_vapor,
        callback.column_cumulative_precipitation,
        vars.parameterizations.surface_humidity_flux.data,
        vars.parameterizations.rain_rate.data,
        vars.parameterizations.snow_rate.data,
        model.time_stepping.Δt_sec,
    )
    sample = callback.timestep_counter % callback.sample_every_n_steps == 0 ||
             callback.timestep_counter == callback.final_timestep
    sample && _sample_global_atmosphere!(callback, vars, model)
    return nothing
end

SpeedyWeather.finalize!(::GlobalAtmosphereDiagnosticsCallback, args...) = nothing


# SpeedyWeather's built-in global-temperature callback scalar-indexes its
# architecture-resident `temp_average` vector. These model-specific methods are
# strictly more specific than the upstream AbstractModel methods and bulk-copy
# the very small vertical-mean vector only on GPU.
function _lowest_layer_global_temperature(vars, model)
    nlayers = model.geometry.nlayers
    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        return vars.grid.temp_average[nlayers]
    end
    return Array(vars.grid.temp_average)[nlayers]
end

function SpeedyWeather.initialize!(
    callback::SpeedyWeather.GlobalSurfaceTemperatureCallback{NF},
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.PrimitiveEquation,
) where {NF}
    callback.temperature = Vector{NF}(undef, vars.prognostic.clock.n_timesteps + 1)
    callback.temperature[1] = _lowest_layer_global_temperature(vars, model)
    callback.timestep_counter = 1
    return nothing
end

function SpeedyWeather.callback!(
    callback::SpeedyWeather.GlobalSurfaceTemperatureCallback,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.PrimitiveEquation,
)
    callback.timestep_counter += 1
    callback.temperature[callback.timestep_counter] =
        _lowest_layer_global_temperature(vars, model)
    return nothing
end

function _enable_diagnostic_clouds(base_solver)
    grid_params = base_solver.grid_params
    NF = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    nlay = grid_params.nlay
    ncol = grid_params.ncol
    # Retain the parallel clear-sky solve so coupled production diagnostics can
    # separate cloud radiative effects from gas/surface contributions.  This is
    # diagnostic only: the all-sky flux remains the applied flux.
    method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(false, true)
    lookups = RRTMGP.lookup_tables(grid_params, method)
    cloud_state = RRTMGP.AtmosphericStates.CloudState(
        DA{NF}(fill(NF(10), nlay, ncol)),
        DA{NF}(fill(NF(30), nlay, ncol)),
        DA{NF}(zeros(NF, nlay, ncol)),
        DA{NF}(zeros(NF, nlay, ncol)),
        DA{NF}(zeros(NF, nlay, ncol)),
        DA{NF}(zeros(NF, ncol)),
        DA{NF}(zeros(NF, ncol)),
        DA{Bool}(falses(nlay, ncol)),
        DA{Bool}(falses(nlay, ncol)),
        DeterministicMaxRandomOverlap(DA, ncol),
        2,
    )
    old_state = base_solver.as
    state = RRTMGP.AtmosphericStates.AtmosphericState(
        old_state.lon,
        old_state.lat,
        old_state.layerdata,
        old_state.p_lev,
        old_state.t_lev,
        old_state.t_sfc,
        old_state.vmr,
        cloud_state,
        old_state.aerosol_state,
    )
    return RRTMGP.RRTMGPSolver(
        grid_params,
        method,
        base_solver.params,
        base_solver.lws.bcs,
        base_solver.sws.bcs,
        state;
        lookups,
    )
end

function _enable_sulfate_aerosol(base_solver, target_aod; all_sky = false)
    grid_params = base_solver.grid_params
    NF = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    nlay = grid_params.nlay
    ncol = grid_params.ncol
    method = all_sky ?
        RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(true, true) :
        RRTMGP.ClearSkyRadiation(true)
    lookups = RRTMGP.lookup_tables(grid_params, method)
    n_aerosols = length(RRTMGP.aerosol_names())
    n_aerosize = maximum(values(lookups.idx_aerosize_sw))
    aerosol_state = RRTMGP.AtmosphericStates.AerosolState(
        DA{NF}(zeros(NF, ncol)),
        DA{NF}(zeros(NF, ncol)),
        DA{Bool}(falses(nlay, ncol)),
        DA{NF}(zeros(NF, n_aerosize, nlay, ncol)),
        DA{NF}(zeros(NF, n_aerosols, nlay, ncol)),
    )
    old_state = base_solver.as
    state = RRTMGP.AtmosphericStates.AtmosphericState(
        old_state.lon,
        old_state.lat,
        old_state.layerdata,
        old_state.p_lev,
        old_state.t_lev,
        old_state.t_sfc,
        old_state.vmr,
        old_state.cloud_state,
        aerosol_state,
    )
    solver = RRTMGP.RRTMGPSolver(
        grid_params,
        method,
        base_solver.params,
        base_solver.lws.bcs,
        base_solver.sws.bcs,
        state;
        lookups,
    )

    # Prescribe a lower-stratospheric sulfate layer and calibrate its column
    # mass separately in every column to the requested 550 nm AOD.
    sulfate_mass = RRTMGP.aerosol_column_mass_density(solver, "sulfate")
    p_lay = RRTMGP.layer_pressure(solver)
    relative_humidity = RRTMGP.layer_relative_humidity(solver)
    relative_humidity .= NF(0.3)
    trial_column_mass = NF(1.0e-5)
    sulfate_mass .= @. exp(
        -(log(max(p_lay, one(NF)) / NF(15_000)))^2 / NF(0.8),
    )
    sulfate_mass .*= trial_column_mass ./ sum(sulfate_mass; dims = 1)
    _reset_solver_cloud_sampler!(solver, 0)
    RRTMGP.update_fluxes!(solver)
    realized = RRTMGP.aod_sw_extinction(solver)
    sulfate_mass .*= reshape(NF(target_aod) ./ max.(realized, eps(NF)), 1, ncol)
    _reset_solver_cloud_sampler!(solver, 0)
    RRTMGP.update_fluxes!(solver)
    return solver
end

function RRTMGPRadiation(
    spectral_grid::SpectralGrid;
    co2_ppm = 420,
    aerosol_optical_depth_550nm = 0,
    surface_emissivity = 0.98,
    solve_every_n_steps = 1,
    cloud_scheme = :none,
    cloud_humidity_search_min_sigma = 0,
    cloud_liquid_water_path_gm2 = 60,
    cloud_ice_water_path_gm2 = 25,
    prognostic_cloud_condensate = nothing,
)
    NF = spectral_grid.NF
    nlay = spectral_grid.nlayers
    ncol = spectral_grid.npoints
    _, latitude = RG.get_londlatds(spectral_grid.grid)
    profile = RRTMGP.standard_atmosphere(NF; nlay, ncol)
    profile.lat .= NF.(latitude)
    profile.well_mixed_vmr["co2"] = NF(co2_ppm * 1.0e-6)
    context = if spectral_grid.architecture isa SpeedyWeather.GPU
        ClimaComms.context(ClimaComms.CUDADevice())
    else
        ClimaComms.context(ClimaComms.CPUSingleThreaded())
    end
    output = RRTMGP.solve(
        profile;
        context,
        method = RRTMGP.ClearSkyRadiation(false),
        surface_emissivity = NF(surface_emissivity),
    )
    cloud_scheme in (:none, _ALL_SKY_CLOUD_SCHEMES...) ||
        throw(ArgumentError(
            "cloud_scheme must be none, diagnostic, or prognostic_condensate",
        ))
    0 <= cloud_humidity_search_min_sigma <= 1 || throw(ArgumentError(
        "cloud_humidity_search_min_sigma must be within [0, 1]",
    ))
    if cloud_scheme == :prognostic_condensate
        isnothing(prognostic_cloud_condensate) && throw(ArgumentError(
            "prognostic_condensate requires a shared condensate state",
        ))
    elseif !isnothing(prognostic_cloud_condensate)
        throw(ArgumentError(
            "a prognostic condensate state requires cloud_scheme=prognostic_condensate",
        ))
    end
    all_sky = cloud_scheme in _ALL_SKY_CLOUD_SCHEMES
    solver = all_sky ? _enable_diagnostic_clouds(output.solver) : output.solver
    if aerosol_optical_depth_550nm > 0
        solver = _enable_sulfate_aerosol(
            solver,
            NF(aerosol_optical_depth_550nm);
            all_sky,
        )
    end
    # Enabling diagnostic clouds constructs a fresh all-sky solver whose flux
    # presentation and combined net-flux buffers are intentionally allocated
    # with `undef`. Populate that solver before deriving the initial heating
    # field. The first live SpeedyWeather parameterization replaces this
    # profile with model state, but reading the fresh buffer here is still an
    # initialization error (and is reported by CUDA Compute Sanitizer).
    _reset_solver_cloud_sampler!(solver, 0)
    RRTMGP.update_fluxes!(solver)
    heating_rate = copy(RRTMGP.heating_rate(solver))
    applied_heating_rate = copy(heating_rate)
    cloud_diagnostic = all_sky ?
        SpeedyWeather.DiagnosticClouds(spectral_grid) : SpeedyWeather.NoClouds(spectral_grid)
    return RRTMGPRadiation(
        solver,
        heating_rate,
        applied_heating_rate,
        NF(co2_ppm),
        NF(surface_emissivity),
        NF(aerosol_optical_depth_550nm),
        Int(solve_every_n_steps),
        0,
        cloud_scheme,
        NF(cloud_humidity_search_min_sigma),
        cloud_diagnostic,
        NF(cloud_liquid_water_path_gm2),
        NF(cloud_ice_water_path_gm2),
        prognostic_cloud_condensate,
    )
end

SpeedyWeather.initialize!(radiation::RRTMGPRadiation, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.initialize!(radiation.clouds, model)
SpeedyWeather.parameterization!(
    ij,
    vars,
    ::RRTMGPRadiation,
    model,
) = nothing
SpeedyWeather.parameterization!(
    ij,
    vars,
    ::RRTMGPColumnNoop,
    model,
) = nothing

@inline function _surface_temperature(ij, vars, model)
    land_fraction = model.land_sea_mask.mask[ij]
    sst = vars.prognostic.ocean.sea_surface_temperature[ij]
    lst = _land_surface_temperature(vars.prognostic.land.soil_temperature, ij)
    # Match the emitted flux of two fractional surfaces under one RRTMGP column.
    return ((1 - land_fraction) * sst^4 + land_fraction * lst^4)^(1 / 4)
end

# GPU kernel argument adaptation exposes SpeedyWeather's one-dimensional land
# temperature storage directly as a device vector. Dispatch on dimensionality
# before the generic RingGrids wrapper path so device code never performs the
# dynamic `field.data` property lookup rejected by GPUCompiler.
@inline _land_surface_temperature(field::AbstractVector, ij) = field[ij]
@inline _land_surface_temperature(field::AbstractMatrix, ij) = field[ij, 1]
@inline _land_surface_temperature(field, ij) =
    ndims(field.data) == 1 ? field[ij] : field[ij, 1]

@inline _rain_rate_mm_per_day(rain_rate_m_per_s) =
    rain_rate_m_per_s * 86_400 * 1_000

@inline function _speedy_static_stability_gradient(
    surface_dry_static_energy,
    above_dry_static_energy,
    surface_geopotential,
    above_geopotential,
)
    geopotential_difference = above_geopotential - surface_geopotential
    iszero(geopotential_difference) && return zero(surface_dry_static_energy)
    return (above_dry_static_energy - surface_dry_static_energy) /
           geopotential_difference
end

@inline function _diagnostic_cloud_humidity_search_start(
    sigma_levels,
    nlayers,
    minimum_sigma,
)
    layer = 3
    final_layer = nlayers - 2
    while layer <= final_layer && sigma_levels[layer] < minimum_sigma
        layer += 1
    end
    return layer
end

@inline function _diagnose_rrtmgp_cloud_properties(
    ij,
    vars,
    clouds,
    model,
    cloud_humidity_search_min_sigma,
)
    temp = vars.grid.temperature_prev
    humid = vars.grid.humidity_prev
    geopotential = vars.grid.geopotential
    NF = eltype(temp)
    nlayers = size(temp, 2)

    surface_pressure = vars.grid.pressure_prev[ij]
    sigma_levels = model.geometry.σ_levels_full
    land_fraction = model.land_sea_mask.mask[ij]
    heat_capacity = model.atmosphere.heat_capacity

    # SpeedyWeather's diagnosed precipitation is a water-depth rate in m s⁻¹.
    # Its v0.21.1 DiagnosticClouds implementation divides by 1000 while claiming
    # to convert to mm day⁻¹, suppressing this term by 10⁶.  Keep the same
    # SPEEDY cloud law but perform the dimensional conversion explicitly here.
    precipitation_mm_per_day = min(
        clouds.precipitation_max,
        _rain_rate_mm_per_day(vars.parameterizations.rain_rate[ij]),
    )
    precipitation_term = clouds.precipitation_weight * sqrt(precipitation_mm_per_day)

    cloud_top_precipitation = vars.parameterizations.cloud_top[ij]
    humidity_excess::NF = 0
    cloud_top_humidity = nlayers + 1

    # Match SPEEDY.f90's cloud search rather than selecting the highest layer
    # that merely passes an RH threshold.  Seed the candidate from the layer
    # immediately above the PBL, then retain the maximum RH excess in free-
    # tropospheric layers 3:L-2. On the original eight-level SPEEDY coordinate,
    # layer 3 is sigma=0.3125. A physical minimum-sigma option preserves that
    # lower boundary when increasing vertical resolution instead of admitting
    # humidity-selected cloud merely because the same integer index moved into
    # the stratosphere. The top two layers and the PBL remain excluded.
    layer_above_pbl = max(1, nlayers - 1)
    humidity_k = humid[ij, layer_above_pbl]
    saturation = SpeedyWeather.saturation_humidity(
        temp[ij, layer_above_pbl],
        sigma_levels[layer_above_pbl] * surface_pressure,
        model.atmosphere,
    )
    if saturation > 0
        relative_humidity = humidity_k / saturation
        if relative_humidity > clouds.relative_humidity_threshold_min
            humidity_excess =
                relative_humidity - clouds.relative_humidity_threshold_min
            cloud_top_humidity = layer_above_pbl
        end
    end

    if nlayers >= 5
        humidity_search_start = _diagnostic_cloud_humidity_search_start(
            sigma_levels,
            nlayers,
            cloud_humidity_search_min_sigma,
        )
        for k in humidity_search_start:(nlayers - 2)
            humidity_k = humid[ij, k]
            saturation = SpeedyWeather.saturation_humidity(
                temp[ij, k],
                sigma_levels[k] * surface_pressure,
                model.atmosphere,
            )
            saturation > 0 || continue
            relative_humidity = humidity_k / saturation
            excess = relative_humidity - clouds.relative_humidity_threshold_min
            if excess > humidity_excess &&
                    humidity_k > clouds.specific_humidity_threshold_min
                humidity_excess = excess
                cloud_top_humidity = k
            end
        end
    end

    humidity_term = min(
        1,
        humidity_excess /
        (clouds.relative_humidity_threshold_max -
         clouds.relative_humidity_threshold_min),
    )^2

    cloud_cover = min(1, precipitation_term + humidity_term)
    cloud_top = min(cloud_top_humidity, cloud_top_precipitation)

    stratocumulus_cover::NF = 0
    if clouds.use_stratocumulus
        surface = nlayers
        above = max(1, nlayers - 1)
        surface_dry_static_energy =
            heat_capacity * temp[ij, surface] + geopotential[ij, surface]
        above_dry_static_energy =
            heat_capacity * temp[ij, above] + geopotential[ij, above]
        static_stability_gradient = _speedy_static_stability_gradient(
            surface_dry_static_energy,
            above_dry_static_energy,
            geopotential[ij, surface],
            geopotential[ij, above],
        )
        static_stability = clamp(
            (static_stability_gradient - clouds.stratocumulus_stability_min) /
            (clouds.stratocumulus_stability_max - clouds.stratocumulus_stability_min),
            0,
            1,
        )
        ocean_cover = static_stability * max(
            clouds.stratocumulus_cover_max -
            clouds.stratocumulus_cloud_factor * cloud_cover,
            0,
        )
        saturation_surface = SpeedyWeather.saturation_humidity(
            temp[ij, surface],
            sigma_levels[surface] * surface_pressure,
            model.atmosphere,
        )
        relative_humidity_surface = humid[ij, surface] / saturation_surface

        # SPEEDY.f90 retains at least 0.15 land stratocumulus cover before
        # multiplying by near-surface RH.  This constant is absent from the
        # installed SpeedyWeather port but is part of the reference scheme.
        land_cover = max(ocean_cover, NF(0.15)) * relative_humidity_surface
        stratocumulus_cover =
            (1 - land_fraction) * ocean_cover + land_fraction * land_cover
    end

    return (; cloud_cover, cloud_top, stratocumulus_cover)
end

@inline function _set_diagnostic_cloud_column!(
    cloud_fraction,
    liquid_path,
    ice_path,
    layer_temperature,
    sigma_half,
    properties,
    liquid_water_path,
    ice_water_path,
    nlay,
    ij,
)
    # In SPEEDY, ordinary cloud extends from its diagnosed top down to the
    # interface above the PBL; cloud_top does not identify a zero-thickness
    # single-layer cloud.  Preserve the configured column-integrated path by
    # distributing it through that depth in proportion to layer mass, with a
    # temperature-dependent liquid/ice partition in each layer.
    if properties.cloud_top < nlay && properties.cloud_cover > 0
        cloud_top = clamp(Int(properties.cloud_top), 1, nlay - 1)
        sigma_depth = sigma_half[nlay] - sigma_half[cloud_top]
        if sigma_depth > 0
            for k in cloud_top:(nlay - 1)
                r = nlay - k + 1
                layer_mass_fraction =
                    (sigma_half[k + 1] - sigma_half[k]) / sigma_depth
                liquid_fraction =
                    clamp((layer_temperature[r, ij] - 253) / 20, 0, 1)
                cloud_fraction[r, ij] = clamp(properties.cloud_cover, 0, 1)
                liquid_path[r, ij] =
                    liquid_water_path * liquid_fraction * layer_mass_fraction
                ice_path[r, ij] =
                    ice_water_path * (1 - liquid_fraction) * layer_mass_fraction
            end
        end
    end
    if properties.stratocumulus_cover > 0
        cloud_fraction[1, ij] = max(
            cloud_fraction[1, ij],
            clamp(properties.stratocumulus_cover, 0, 1),
        )
        liquid_path[1, ij] = max(liquid_path[1, ij], liquid_water_path)
    end
    return nothing
end

function _populate_diagnostic_clouds!(radiation::RRTMGPRadiation, vars, model)
    solver = radiation.solver
    nlay = model.geometry.nlayers
    cloud_fraction = RRTMGP.cloud_fraction(solver)
    liquid_path = RRTMGP.cloud_liquid_water_path(solver)
    ice_path = RRTMGP.cloud_ice_water_path(solver)
    fill!(cloud_fraction, 0)
    fill!(liquid_path, 0)
    fill!(ice_path, 0)
    layer_temperature = RRTMGP.layer_temperature(solver)
    sigma_half = model.geometry.σ_levels_half

    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        for ij in axes(cloud_fraction, 2)
            properties = _diagnose_rrtmgp_cloud_properties(
                ij,
                vars,
                radiation.clouds,
                model,
                radiation.cloud_humidity_search_min_sigma,
            )
            _set_diagnostic_cloud_column!(
                cloud_fraction,
                liquid_path,
                ice_path,
                layer_temperature,
                sigma_half,
                properties,
                radiation.cloud_liquid_water_path_gm2,
                radiation.cloud_ice_water_path_gm2,
                nlay,
                ij,
            )
        end
        return nothing
    end

    backend = KernelAbstractions.get_backend(cloud_fraction)
    _diagnostic_cloud_kernel!(backend)(
        cloud_fraction,
        liquid_path,
        ice_path,
        layer_temperature,
        sigma_half,
        vars,
        radiation.clouds,
        model,
        radiation.cloud_humidity_search_min_sigma,
        radiation.cloud_liquid_water_path_gm2,
        radiation.cloud_ice_water_path_gm2,
        nlay;
        ndrange = size(cloud_fraction, 2),
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

function _update_diagnostic_clouds!(radiation::RRTMGPRadiation, vars, model)
    radiation.cloud_scheme == :diagnostic || return nothing
    return _populate_diagnostic_clouds!(radiation, vars, model)
end

KernelAbstractions.@kernel function _diagnostic_cloud_kernel!(
    cloud_fraction,
    liquid_path,
    ice_path,
    layer_temperature,
    sigma_half,
    vars,
    clouds,
    model,
    cloud_humidity_search_min_sigma,
    liquid_water_path,
    ice_water_path,
    nlay,
)
    ij = @index(Global, Linear)
    properties = _diagnose_rrtmgp_cloud_properties(
        ij,
        vars,
        clouds,
        model,
        cloud_humidity_search_min_sigma,
    )
    _set_diagnostic_cloud_column!(
        cloud_fraction,
        liquid_path,
        ice_path,
        layer_temperature,
        sigma_half,
        properties,
        liquid_water_path,
        ice_water_path,
        nlay,
        ij,
    )
end

const PROGNOSTIC_CONDENSATE_FRACTION_REFERENCE_PATH_GM2 = 100f0

@inline function _set_prognostic_cloud_optics_column!(
    cloud_fraction,
    liquid_path,
    ice_path,
    stored_liquid_path,
    stored_ice_path,
    initial_column_path,
    initialized,
    nlay,
    ij,
)
    first_update = initialized[1] == 0
    initial_path_gm2 = zero(eltype(stored_liquid_path))
    for k in 1:nlay
        r = nlay - k + 1
        diagnosed_fraction = clamp(cloud_fraction[r, ij], 0, 1)
        if first_update
            # RRTMGP paths are in-cloud paths. Store their expected grid-box
            # mass so the first prognostic solve exactly retains the native
            # diagnostic optics rather than launching from a clear atmosphere.
            stored_liquid_path[ij, k] =
                diagnosed_fraction * liquid_path[r, ij]
            stored_ice_path[ij, k] = diagnosed_fraction * ice_path[r, ij]
            initial_path_gm2 +=
                stored_liquid_path[ij, k] + stored_ice_path[ij, k]
        else
            expected_liquid = max(stored_liquid_path[ij, k], 0)
            expected_ice = max(stored_ice_path[ij, k], 0)
            expected_total = expected_liquid + expected_ice
            condensate_fraction = expected_total > 0 ? clamp(
                sqrt(
                    expected_total /
                    eltype(stored_liquid_path)(
                        PROGNOSTIC_CONDENSATE_FRACTION_REFERENCE_PATH_GM2,
                    ),
                ),
                eltype(stored_liquid_path)(0.01),
                one(expected_total),
            ) : zero(expected_total)
            active_fraction = max(diagnosed_fraction, condensate_fraction)
            cloud_fraction[r, ij] = active_fraction
            if active_fraction > 0
                liquid_path[r, ij] = expected_liquid / active_fraction
                ice_path[r, ij] = expected_ice / active_fraction
            else
                liquid_path[r, ij] = 0
                ice_path[r, ij] = 0
            end
        end
    end
    if first_update
        initial_column_path[ij] = initial_path_gm2 / 1_000
    end
    return nothing
end

KernelAbstractions.@kernel function _prognostic_cloud_optics_kernel!(
    cloud_fraction,
    liquid_path,
    ice_path,
    stored_liquid_path,
    stored_ice_path,
    initial_column_path,
    initialized,
    nlay,
)
    ij = @index(Global, Linear)
    _set_prognostic_cloud_optics_column!(
        cloud_fraction,
        liquid_path,
        ice_path,
        stored_liquid_path,
        stored_ice_path,
        initial_column_path,
        initialized,
        nlay,
        ij,
    )
end

function _update_prognostic_clouds!(radiation::RRTMGPRadiation, vars, model)
    radiation.cloud_scheme == :prognostic_condensate || return nothing
    state = radiation.prognostic_cloud_condensate
    isnothing(state) && error(
        "prognostic cloud scheme has no condensate state",
    )

    # Reuse the existing humidity/precipitation law only for cloud occurrence.
    # The liquid and ice optical paths below come from the evolving reservoir.
    _populate_diagnostic_clouds!(radiation, vars, model)
    cloud_fraction = RRTMGP.cloud_fraction(radiation.solver)
    liquid_path = RRTMGP.cloud_liquid_water_path(radiation.solver)
    ice_path = RRTMGP.cloud_ice_water_path(radiation.solver)
    nlay = model.geometry.nlayers

    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        for ij in axes(state.liquid_path_gm2, 1)
            _set_prognostic_cloud_optics_column!(
                cloud_fraction,
                liquid_path,
                ice_path,
                state.liquid_path_gm2,
                state.ice_path_gm2,
                state.initial_column_path_kgm2,
                state.initialized,
                nlay,
                ij,
            )
        end
    else
        backend = KernelAbstractions.get_backend(cloud_fraction)
        _prognostic_cloud_optics_kernel!(backend)(
            cloud_fraction,
            liquid_path,
            ice_path,
            state.liquid_path_gm2,
            state.ice_path_gm2,
            state.initial_column_path_kgm2,
            state.initialized,
            nlay;
            ndrange = size(cloud_fraction, 2),
        )
        KernelAbstractions.synchronize(backend)
    end
    fill!(state.initialized, 1)
    return nothing
end

function _update_rrtmgp_clouds!(radiation::RRTMGPRadiation, vars, model)
    if radiation.cloud_scheme == :diagnostic
        return _update_diagnostic_clouds!(radiation, vars, model)
    elseif radiation.cloud_scheme == :prognostic_condensate
        return _update_prognostic_clouds!(radiation, vars, model)
    end
    return nothing
end

function _update_rrtmgp_state!(radiation::RRTMGPRadiation, vars, model)
    solver = radiation.solver
    nlay = model.geometry.nlayers
    nlev = nlay + 1
    ncol = model.geometry.npoints
    σ_full = model.geometry.σ_levels_full
    σ_half = model.geometry.σ_levels_half
    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    surface_pressure = vars.grid.pressure_prev

    p_lay = RRTMGP.layer_pressure(solver)
    p_lev = RRTMGP.level_pressure(solver)
    t_lay = RRTMGP.layer_temperature(solver)
    t_lev = RRTMGP.level_temperature(solver)
    t_sfc = RRTMGP.surface_temperature(solver)
    vmr_h2o = RRTMGP.volume_mixing_ratio(solver, "h2o")
    vmr_o3 = RRTMGP.volume_mixing_ratio(solver, "o3")
    relative_humidity = RRTMGP.layer_relative_humidity(solver)

    molar_ratio = model.atmosphere.mol_mass_dry_air / model.atmosphere.mol_mass_vapor
    backend = KernelAbstractions.get_backend(p_lay)
    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        for ij in 1:ncol
            ps = surface_pressure[ij]
            ts = _surface_temperature(ij, vars, model)
            t_sfc[ij] = ts
            for k in 1:nlay
                r = nlay - k + 1
                p = max(ps * σ_full[k], one(ps))
                q = clamp(humidity[ij, k], zero(ps), typeof(ps)(0.5))
                p_lay[r, ij] = p
                t_lay[r, ij] = temperature[ij, k]
                vmr_h2o[r, ij] = q / (1 - q) * molar_ratio
                vapor_pressure = vmr_h2o[r, ij] * p / (1 + vmr_h2o[r, ij])
                saturation_pressure = typeof(ps)(611.2) * exp(
                    typeof(ps)(17.67) * (t_lay[r, ij] - typeof(ps)(273.15)) /
                        (t_lay[r, ij] - typeof(ps)(29.65)),
                )
                relative_humidity[r, ij] = clamp(
                    vapor_pressure / saturation_pressure,
                    zero(ps),
                    one(ps),
                )
                vmr_o3[r, ij] =
                    typeof(ps)(3.0e-8) + typeof(ps)(7.5e-6) *
                    exp(-(log(p / typeof(ps)(1200)))^2 / typeof(ps)(2.88))
            end
            for k in 1:nlev
                r = nlev - k + 1
                p_lev[r, ij] = max(ps * σ_half[k], one(ps))
            end
            t_lev[1, ij] = ts
            for r in 2:nlay
                t_lev[r, ij] = (t_lay[r - 1, ij] + t_lay[r, ij]) / 2
            end
            t_lev[nlev, ij] = t_lay[nlay, ij]
        end
    else
        _pack_rrtmgp_state_kernel!(backend)(
            p_lay,
            p_lev,
            t_lay,
            t_lev,
            t_sfc,
            vmr_h2o,
            vmr_o3,
            relative_humidity,
            temperature,
            humidity,
            surface_pressure,
            vars,
            model,
            σ_full,
            σ_half,
            molar_ratio,
            nlay,
            nlev;
            ndrange = ncol,
        )
        KernelAbstractions.synchronize(backend)
    end

    # Solar geometry is a global SpeedyWeather parameterization and has already run.
    RRTMGP.cos_zenith(solver) .= clamp.(vars.parameterizations.cos_zenith.data, 0, 1)
    RRTMGP.toa_flux(solver) .= model.planet.solar_constant
    RRTMGP.surface_emissivity(solver) .= radiation.surface_emissivity
    RRTMGP.set_volume_mixing_ratio!(solver, "co2", radiation.co2_ppm * 1.0e-6)

    direct_albedo = RRTMGP.direct_sw_surface_albedo(solver)
    diffuse_albedo = RRTMGP.diffuse_sw_surface_albedo(solver)
    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        for ij in 1:ncol
            SpeedyWeather.parameterization!(ij, vars, model.albedo, model)
            land_fraction = model.land_sea_mask.mask[ij]
            ocean_albedo = vars.parameterizations.ocean.albedo[ij]
            land_albedo = vars.parameterizations.land.albedo[ij]
            albedo = (1 - land_fraction) * ocean_albedo + land_fraction * land_albedo
            vars.parameterizations.albedo[ij] = albedo
            direct_albedo[:, ij] .= albedo
            diffuse_albedo[:, ij] .= albedo
        end
    else
        _pack_rrtmgp_albedo_kernel!(backend)(
            direct_albedo,
            diffuse_albedo,
            vars,
            model.albedo,
            model;
            ndrange = ncol,
        )
        KernelAbstractions.synchronize(backend)
    end
    _update_rrtmgp_clouds!(radiation, vars, model)
    return nothing
end

KernelAbstractions.@kernel function _pack_rrtmgp_state_kernel!(
    p_lay,
    p_lev,
    t_lay,
    t_lev,
    t_sfc,
    vmr_h2o,
    vmr_o3,
    relative_humidity,
    temperature,
    humidity,
    surface_pressure,
    vars,
    model,
    σ_full,
    σ_half,
    molar_ratio,
    nlay,
    nlev,
)
    ij = @index(Global, Linear)
    ps = surface_pressure[ij]
    ts = _surface_temperature(ij, vars, model)
    t_sfc[ij] = ts
    for k in 1:nlay
        r = nlay - k + 1
        p = max(ps * σ_full[k], one(ps))
        q = clamp(humidity[ij, k], zero(ps), typeof(ps)(0.5))
        p_lay[r, ij] = p
        t_lay[r, ij] = temperature[ij, k]
        vmr_h2o[r, ij] = q / (1 - q) * molar_ratio
        vapor_pressure = vmr_h2o[r, ij] * p / (1 + vmr_h2o[r, ij])
        saturation_pressure = typeof(ps)(611.2) * exp(
            typeof(ps)(17.67) * (t_lay[r, ij] - typeof(ps)(273.15)) /
                (t_lay[r, ij] - typeof(ps)(29.65)),
        )
        relative_humidity[r, ij] = clamp(
            vapor_pressure / saturation_pressure,
            zero(ps),
            one(ps),
        )
        vmr_o3[r, ij] =
            typeof(ps)(3.0e-8) + typeof(ps)(7.5e-6) *
            exp(-(log(p / typeof(ps)(1200)))^2 / typeof(ps)(2.88))
    end
    for k in 1:nlev
        r = nlev - k + 1
        p_lev[r, ij] = max(ps * σ_half[k], one(ps))
    end
    t_lev[1, ij] = ts
    for r in 2:nlay
        t_lev[r, ij] = (t_lay[r - 1, ij] + t_lay[r, ij]) / 2
    end
    t_lev[nlev, ij] = t_lay[nlay, ij]
end

KernelAbstractions.@kernel function _pack_rrtmgp_albedo_kernel!(
    direct_albedo,
    diffuse_albedo,
    vars,
    albedo_scheme,
    model,
)
    ij = @index(Global, Linear)
    # SpeedyWeather passes GPU-compatible parameterizations separately from
    # its device-adapted model for the same reason: non-device model fields are
    # intentionally omitted by adaptation. Pass the scheme explicitly rather
    # than looking up `model.albedo` inside device code.
    SpeedyWeather.parameterization!(ij, vars, albedo_scheme, model)
    land_fraction = model.land_sea_mask.mask[ij]
    ocean_albedo = vars.parameterizations.ocean.albedo[ij]
    land_albedo = vars.parameterizations.land.albedo[ij]
    albedo = (1 - land_fraction) * ocean_albedo + land_fraction * land_albedo
    vars.parameterizations.albedo[ij] = albedo
    for band in axes(direct_albedo, 1)
        direct_albedo[band, ij] = albedo
        diffuse_albedo[band, ij] = albedo
    end
end

function _apply_rrtmgp_fluxes!(radiation::RRTMGPRadiation, vars, model)
    solver = radiation.solver
    nlay = model.geometry.nlayers
    ncol = model.geometry.npoints
    # Preserve the exact field used by this step. The post-physics callback may
    # replace `heating_rate` before diagnostics run, and SpeedyWeather stores
    # callbacks in a Dict whose iteration order is not a scheduling contract.
    radiation.applied_heating_rate .= radiation.heating_rate
    heating = radiation.applied_heating_rate
    tendency = vars.tendencies.grid.temperature
    lw_up = RRTMGP.lw_flux_up(solver)
    lw_down = RRTMGP.lw_flux_dn(solver)
    sw_up = RRTMGP.sw_flux_up(solver)
    sw_down = RRTMGP.sw_flux_dn(solver)
    stefan_boltzmann = model.atmosphere.stefan_boltzmann

    if model.spectral_grid.architecture isa SpeedyWeather.CPU
        for ij in 1:ncol
            for k in 1:nlay
                tendency[ij, k] += heating[nlay - k + 1, ij]
            end
            sst = vars.prognostic.ocean.sea_surface_temperature[ij]
            lst = _land_surface_temperature(vars.prognostic.land.soil_temperature, ij)
            ocean_lw_up = radiation.surface_emissivity * stefan_boltzmann * sst^4
            land_lw_up = radiation.surface_emissivity * stefan_boltzmann * lst^4
            surface_sw_down = sw_down[1, ij]
            ocean_sw_up = surface_sw_down * vars.parameterizations.ocean.albedo[ij]
            land_sw_up = surface_sw_down * vars.parameterizations.land.albedo[ij]

            vars.parameterizations.surface_longwave_down[ij] = lw_down[1, ij]
            vars.parameterizations.surface_longwave_up[ij] = lw_up[1, ij]
            vars.parameterizations.ocean.surface_longwave_up[ij] = ocean_lw_up
            vars.parameterizations.land.surface_longwave_up[ij] = land_lw_up
            vars.parameterizations.surface_shortwave_down[ij] = surface_sw_down
            vars.parameterizations.ocean.surface_shortwave_down[ij] = surface_sw_down
            vars.parameterizations.land.surface_shortwave_down[ij] = surface_sw_down
            vars.parameterizations.surface_shortwave_up[ij] = sw_up[1, ij]
            vars.parameterizations.ocean.surface_shortwave_up[ij] = ocean_sw_up
            vars.parameterizations.land.surface_shortwave_up[ij] = land_sw_up
            vars.parameterizations.outgoing_longwave[ij] = lw_up[end, ij]
            vars.parameterizations.outgoing_shortwave[ij] = sw_up[end, ij]
        end
        return nothing
    end

    backend = KernelAbstractions.get_backend(heating)
    _apply_rrtmgp_fluxes_kernel!(backend)(
        heating,
        tendency,
        lw_up,
        lw_down,
        sw_up,
        sw_down,
        vars,
        radiation.surface_emissivity,
        stefan_boltzmann,
        nlay;
        ndrange = ncol,
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

KernelAbstractions.@kernel function _apply_rrtmgp_fluxes_kernel!(
    heating,
    tendency,
    lw_up,
    lw_down,
    sw_up,
    sw_down,
    vars,
    surface_emissivity,
    stefan_boltzmann,
    nlay,
)
    ij = @index(Global, Linear)
    for k in 1:nlay
        tendency[ij, k] += heating[nlay - k + 1, ij]
    end
    sst = vars.prognostic.ocean.sea_surface_temperature[ij]
    lst = _land_surface_temperature(vars.prognostic.land.soil_temperature, ij)
    ocean_lw_up = surface_emissivity * stefan_boltzmann * sst^4
    land_lw_up = surface_emissivity * stefan_boltzmann * lst^4
    surface_sw_down = sw_down[1, ij]
    ocean_sw_up = surface_sw_down * vars.parameterizations.ocean.albedo[ij]
    land_sw_up = surface_sw_down * vars.parameterizations.land.albedo[ij]

    vars.parameterizations.surface_longwave_down[ij] = lw_down[1, ij]
    vars.parameterizations.surface_longwave_up[ij] = lw_up[1, ij]
    vars.parameterizations.ocean.surface_longwave_up[ij] = ocean_lw_up
    vars.parameterizations.land.surface_longwave_up[ij] = land_lw_up
    vars.parameterizations.surface_shortwave_down[ij] = surface_sw_down
    vars.parameterizations.ocean.surface_shortwave_down[ij] = surface_sw_down
    vars.parameterizations.land.surface_shortwave_down[ij] = surface_sw_down
    vars.parameterizations.surface_shortwave_up[ij] = sw_up[1, ij]
    vars.parameterizations.ocean.surface_shortwave_up[ij] = ocean_sw_up
    vars.parameterizations.land.surface_shortwave_up[ij] = land_sw_up
    vars.parameterizations.outgoing_longwave[ij] = lw_up[end, ij]
    vars.parameterizations.outgoing_shortwave[ij] = sw_up[end, ij]
end

function _update_fluxes!(radiation::RRTMGPRadiation)
    seed = radiation.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES ?
        radiation.call_counter : nothing
    if !isnothing(seed)
        _reset_solver_cloud_sampler!(radiation.solver, seed)
    end
    RRTMGP.update_fluxes!(radiation.solver, seed)
    return nothing
end

function _recalibrate_aod!(radiation::RRTMGPRadiation)
    radiation.target_aod_550nm > 0 || return nothing
    solver = radiation.solver
    sulfate_mass = RRTMGP.aerosol_column_mass_density(solver, "sulfate")
    realized = RRTMGP.aod_sw_extinction(solver)
    sulfate_mass .*= reshape(
        radiation.target_aod_550nm ./ max.(realized, eps(eltype(realized))),
        1,
        size(sulfate_mass, 2),
    )
    # Re-solve with the adjusted burden so the applied fluxes match target AOD.
    _update_fluxes!(radiation)
    return nothing
end

function _synchronize_rrtmgp_postphysics()
    return lowercase(
        get(ENV, "READYESM_SYNC_RRTMGP_POSTPHYSICS", "false"),
    ) in ("1", "true", "yes", "on")
end

function _solve_rrtmgp!(radiation::RRTMGPRadiation, vars, model)
    _update_rrtmgp_state!(radiation, vars, model)
    _update_fluxes!(radiation)
    _recalibrate_aod!(radiation)
    radiation.heating_rate .= RRTMGP.heating_rate(radiation.solver)
    if _synchronize_rrtmgp_postphysics()
        backend = KernelAbstractions.get_backend(radiation.heating_rate)
        KernelAbstractions.synchronize(backend)
    end
    return nothing
end

"""
Refresh RRTMGP after SpeedyWeather has completed its column physics.

SpeedyWeather resets precipitation and cloud-top work arrays immediately before
running global parameterizations, while condensation and convection run later in
the fused column-physics pass.  Solving RRTMGP in its global parameterization
therefore cannot see the cloud diagnostics from that timestep.  This callback
runs after a completed atmospheric step and refreshes the cached radiation every
configured interval; the resulting fluxes and heating are applied explicitly on
the following step.
"""
Base.@kwdef mutable struct RRTMGPPostPhysicsUpdateCallback{R} <:
                           SpeedyWeather.AbstractCallback
    radiation::R
    timestep_counter::Int = 0
end

function SpeedyWeather.initialize!(callback::RRTMGPPostPhysicsUpdateCallback, args...)
    callback.timestep_counter = 0
    return nothing
end

function SpeedyWeather.callback!(
    callback::RRTMGPPostPhysicsUpdateCallback,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.PrimitiveEquation,
)
    callback.timestep_counter += 1
    radiation = callback.radiation
    if callback.timestep_counter % radiation.solve_every_n_steps == 0
        _solve_rrtmgp!(radiation, vars, model)
    end
    return nothing
end

SpeedyWeather.finalize!(::RRTMGPPostPhysicsUpdateCallback, args...) = nothing

function SpeedyWeather.parameterization!(
    vars::SpeedyWeather.Variables,
    radiation::RRTMGPRadiation,
    model::SpeedyWeather.PrimitiveEquation,
)
    radiation.call_counter += 1
    # Supply a valid initial radiative state for the first Euler tendency.
    # Later refreshes occur in RRTMGPPostPhysicsUpdateCallback, where the
    # current step's diagnosed precipitation and cloud top are still present.
    radiation.call_counter == 1 && _solve_rrtmgp!(radiation, vars, model)
    _apply_rrtmgp_fluxes!(radiation, vars, model)
    return nothing
end
