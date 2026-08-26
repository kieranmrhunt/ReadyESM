"""
Exact checkpoint state for the SpeedyWeather component of a dynamic ReadyESM run.

NumericalEarth and Oceananigans already provide structured checkpoint methods for
the coupled clock, ocean, z-star coordinate, and ClimaSeaIce. Their generic
fallback would, however, try to serialize the complete function-bearing
`SpeedyWeather.Simulation`. These methods retain the complete Variables tree
(including leapfrog, diagnostic scratch, and Terrarium state) plus only the
mutable RRTMGP and ReadyESM callback state consumed after a restart.
"""

_restart_host_array(values::AbstractArray) = Array(values)

# Surface-radiation diagnostics are not forcings themselves, but retaining them
# makes a checkpoint a complete reportable state rather than only the minimum
# state required to advance the next timestep.
function _restart_container_state(container)
    isnothing(container) && return nothing
    names = propertynames(container)
    values = map(names) do name
        Oceananigans.prognostic_state(getproperty(container, name))
    end
    return NamedTuple{names}(values)
end

function _restore_restart_container!(container, state)
    isnothing(state) && return container
    for name in keys(state)
        Oceananigans.restore_prognostic_state!(
            getproperty(container, name),
            getproperty(state, name),
        )
    end
    return container
end

function _radiation_interface_flux_state(interface_fluxes)
    isnothing(interface_fluxes) && return nothing
    names = keys(interface_fluxes)
    values = map(names) do name
        _restart_container_state(getproperty(interface_fluxes, name))
    end
    return NamedTuple{names}(values)
end

function _restore_radiation_interface_flux_state!(interface_fluxes, state)
    isnothing(state) && return interface_fluxes
    isnothing(interface_fluxes) && error(
        "checkpoint contains surface-radiation interface fluxes but the reconstructed radiation bridge does not",
    )
    for name in keys(state)
        haskey(interface_fluxes, name) || error(
            "checkpoint surface-radiation interface $name is absent from the reconstructed bridge",
        )
        _restore_restart_container!(
            getproperty(interface_fluxes, name),
            getproperty(state, name),
        )
    end
    return interface_fluxes
end

function Oceananigans.prognostic_state(radiation::AtmosphereDrivenRadiation)
    return (
        interface_fluxes = _radiation_interface_flux_state(
            radiation.interface_fluxes,
        ),
    )
end

function Oceananigans.restore_prognostic_state!(
    radiation::AtmosphereDrivenRadiation,
    state,
)
    _restore_radiation_interface_flux_state!(
        radiation.interface_fluxes,
        state.interface_fluxes,
    )
    return radiation
end

Oceananigans.restore_prognostic_state!(
    radiation::AtmosphereDrivenRadiation,
    ::Nothing,
) = radiation

function _restart_copy_array!(destination::AbstractArray, source, label)
    size(destination) == size(source) || throw(DimensionMismatch(
        "$label checkpoint size $(size(source)) does not match destination $(size(destination))",
    ))
    copyto!(destination, source)
    return destination
end

function _restart_copy_lower_triangular!(destination, source)
    size(destination) == size(source) || throw(DimensionMismatch(
        "SpeedyWeather lower-triangular checkpoint size $(size(source)) " *
        "does not match destination $(size(destination))",
    ))
    eltype(destination) == eltype(source) || throw(ArgumentError(
        "SpeedyWeather lower-triangular checkpoint eltype " *
        "$(eltype(source)) does not match destination $(eltype(destination))",
    ))
    # LowerTriangularArrays.copyto! uses broadcast. A GPU destination then
    # tries to capture the host source in a device kernel, which is invalid.
    # Its flattened storages have identical layout at an exact restart, and
    # CUDA's ordinary array copy has a supported host-to-device path.
    copyto!(destination.data, source.data)
    return destination
end

# SpeedyWeather's Variables copy owns the alias/view policy, so intercept only
# its unsupported concrete cross-architecture leaf. In particular, do not run
# `on_architecture` over the complete Variables tree: ReadyESM attaches a
# Terrarium.StateVariables object, which is intentionally outside
# SpeedyWeather's architecture-conversion method set.
function SpeedyWeather._copy_entry!(
    destination::SpeedyWeather.LowerTriangularArray{T, N, D},
    source::SpeedyWeather.LowerTriangularArray{T, N, S},
) where {T, N, D <: CUDA.CuArray, S <: Array}
    return _restart_copy_lower_triangular!(destination, source)
end

function _restart_copy_ring_field!(destination, source)
    size(destination) == size(source) || throw(DimensionMismatch(
        "SpeedyWeather RingGrids.Field checkpoint size $(size(source)) " *
        "does not match destination $(size(destination))",
    ))
    eltype(destination) == eltype(source) || throw(ArgumentError(
        "SpeedyWeather RingGrids.Field checkpoint eltype $(eltype(source)) " *
        "does not match destination $(eltype(destination))",
    ))
    # RingGrids.Field's generic AbstractArray copy iterates its logical
    # indices. With a GPU-backed destination that reaches scalar setindex! on
    # the CuArray. The Field wrappers have identical flattened storage at an
    # exact restart, so transfer their data arrays through CUDA's supported
    # host-to-device copy path while retaining the reconstructed device grid.
    copyto!(destination.data, source.data)
    return destination
end

function SpeedyWeather._copy_entry!(
    destination::SpeedyWeather.RingGrids.Field{T, N, D, GD},
    source::SpeedyWeather.RingGrids.Field{T, N, S, GS},
) where {
    T,
    N,
    D <: CUDA.CuArray,
    S <: Array,
    GD <: SpeedyWeather.RingGrids.AbstractGrid,
    GS <: SpeedyWeather.RingGrids.AbstractGrid,
}
    return _restart_copy_ring_field!(destination, source)
end

function _restart_copy_speedy_scratch_memory!(destination, source)
    _restart_copy_array!(
        destination.north,
        source.north,
        "SpeedyTransforms northern scratch memory",
    )
    _restart_copy_array!(
        destination.south,
        source.south,
        "SpeedyTransforms southern scratch memory",
    )
    _restart_copy_array!(
        destination.column.north,
        source.column.north,
        "SpeedyTransforms northern column scratch memory",
    )
    _restart_copy_array!(
        destination.column.south,
        source.column.south,
        "SpeedyTransforms southern column scratch memory",
    )
    return destination
end

# SpeedyWeather's generated fallback requires destination and source to have
# the same concrete ScratchMemory type. JLD2 reconstructs the source with
# Arrays while the live atmosphere owns CuArrays, so that fallback cannot
# dispatch. Copy all four work buffers explicitly; this also preserves the
# nested ColumnScratchMemory arrays that the upstream generated fallback does
# not recurse into.
function SpeedyWeather._copy_entry!(
    destination::SpeedyWeather.SpeedyTransforms.ScratchMemory{D3, D1},
    source::SpeedyWeather.SpeedyTransforms.ScratchMemory{S3, S1},
) where {
    D3 <: CUDA.CuArray,
    D1 <: CUDA.CuArray,
    S3 <: Array,
    S1 <: Array,
}
    return _restart_copy_speedy_scratch_memory!(destination, source)
end

function _restart_resize_copy!(destination::Vector, source, label)
    resize!(destination, length(source))
    copyto!(destination, source)
    return destination
end

function _rrtmgp_clear_flux_state(solver)
    isnothing(solver.clear_flux_lw) && return nothing
    return (
        lw_up = _restart_host_array(RRTMGP.clear_lw_flux_up(solver)),
        lw_down = _restart_host_array(RRTMGP.clear_lw_flux_dn(solver)),
        lw_net = _restart_host_array(RRTMGP.clear_lw_flux(solver)),
        sw_up = _restart_host_array(RRTMGP.clear_sw_flux_up(solver)),
        sw_down = _restart_host_array(RRTMGP.clear_sw_flux_dn(solver)),
        sw_net = _restart_host_array(RRTMGP.clear_sw_flux(solver)),
        sw_direct_down = _restart_host_array(
            RRTMGP.clear_sw_direct_flux_dn(solver),
        ),
        net = _restart_host_array(RRTMGP.clear_net_flux(solver)),
    )
end

function _restore_rrtmgp_clear_flux_state!(solver, state)
    isnothing(state) && return nothing
    isnothing(solver.clear_flux_lw) && error(
        "checkpoint contains clear-sky RRTMGP fluxes but the reconstructed solver does not",
    )
    _restart_copy_array!(RRTMGP.clear_lw_flux_up(solver), state.lw_up, "clear LW up")
    _restart_copy_array!(RRTMGP.clear_lw_flux_dn(solver), state.lw_down, "clear LW down")
    _restart_copy_array!(RRTMGP.clear_lw_flux(solver), state.lw_net, "clear LW net")
    _restart_copy_array!(RRTMGP.clear_sw_flux_up(solver), state.sw_up, "clear SW up")
    _restart_copy_array!(RRTMGP.clear_sw_flux_dn(solver), state.sw_down, "clear SW down")
    _restart_copy_array!(RRTMGP.clear_sw_flux(solver), state.sw_net, "clear SW net")
    _restart_copy_array!(
        RRTMGP.clear_sw_direct_flux_dn(solver),
        state.sw_direct_down,
        "clear direct SW down",
    )
    _restart_copy_array!(RRTMGP.clear_net_flux(solver), state.net, "clear net flux")
    return nothing
end

function _rrtmgp_cloud_sampler_state(solver)
    cloud_state = solver.as.cloud_state
    hasproperty(cloud_state, :mask_type) || return nothing
    sampler = cloud_state.mask_type
    sampler isa DeterministicMaxRandomOverlap || return nothing
    return (
        counters = _restart_host_array(sampler.counters),
        seed = _restart_host_array(sampler.seed),
    )
end

function _restore_rrtmgp_cloud_sampler_state!(solver, state)
    isnothing(state) && return nothing
    cloud_state = solver.as.cloud_state
    hasproperty(cloud_state, :mask_type) || error(
        "checkpoint contains deterministic cloud state but the reconstructed solver does not",
    )
    sampler = cloud_state.mask_type
    sampler isa DeterministicMaxRandomOverlap || error(
        "checkpoint requires ReadyESM's deterministic RRTMGP cloud sampler",
    )
    _restart_copy_array!(sampler.counters, state.counters, "cloud sample counters")
    _restart_copy_array!(sampler.seed, state.seed, "cloud sample seed")
    return nothing
end

function _rrtmgp_prognostic_condensate_state(radiation::RRTMGPRadiation)
    condensate = radiation.prognostic_cloud_condensate
    isnothing(condensate) && return nothing
    return (
        liquid_path_gm2 = _restart_host_array(condensate.liquid_path_gm2),
        ice_path_gm2 = _restart_host_array(condensate.ice_path_gm2),
        initial_column_path_kgm2 =
            _restart_host_array(condensate.initial_column_path_kgm2),
        cumulative_retained_kgm2 =
            _restart_host_array(condensate.cumulative_retained_kgm2),
        cumulative_sedimented_kgm2 =
            _restart_host_array(condensate.cumulative_sedimented_kgm2),
        initialized = _restart_host_array(condensate.initialized),
        retention_fraction = condensate.retention_fraction,
        residence_time_seconds = condensate.residence_time_seconds,
    )
end

function _restore_rrtmgp_prognostic_condensate_state!(radiation, state)
    condensate = radiation.prognostic_cloud_condensate
    if isnothing(state)
        isnothing(condensate) || error(
            "checkpoint has no prognostic condensate but reconstructed radiation does",
        )
        return nothing
    end
    isnothing(condensate) && error(
        "checkpoint contains prognostic condensate but reconstructed radiation does not",
    )
    condensate.retention_fraction == state.retention_fraction || error(
        "checkpoint prognostic-condensate retention fraction differs",
    )
    condensate.residence_time_seconds == state.residence_time_seconds || error(
        "checkpoint prognostic-condensate residence time differs",
    )
    for (name, label) in (
        (:liquid_path_gm2, "prognostic cloud liquid path"),
        (:ice_path_gm2, "prognostic cloud ice path"),
        (:initial_column_path_kgm2, "initial prognostic cloud path"),
        (:cumulative_retained_kgm2, "retained prognostic cloud water"),
        (:cumulative_sedimented_kgm2, "sedimented prognostic cloud water"),
        (:initialized, "prognostic cloud initialization flag"),
    )
        _restart_copy_array!(
            getproperty(condensate, name),
            getproperty(state, name),
            label,
        )
    end
    return nothing
end

function _rrtmgp_radiation_state(radiation::RRTMGPRadiation)
    solver = radiation.solver
    sulfate_mass = radiation.target_aod_550nm > 0 ?
        _restart_host_array(
            RRTMGP.aerosol_column_mass_density(solver, "sulfate"),
        ) : nothing
    return (
        call_counter = radiation.call_counter,
        heating_rate = _restart_host_array(radiation.heating_rate),
        applied_heating_rate = _restart_host_array(radiation.applied_heating_rate),
        lw_up = _restart_host_array(RRTMGP.lw_flux_up(solver)),
        lw_down = _restart_host_array(RRTMGP.lw_flux_dn(solver)),
        lw_net = _restart_host_array(RRTMGP.lw_flux_net(solver)),
        sw_up = _restart_host_array(RRTMGP.sw_flux_up(solver)),
        sw_down = _restart_host_array(RRTMGP.sw_flux_dn(solver)),
        sw_net = _restart_host_array(RRTMGP.sw_flux_net(solver)),
        sw_direct_down = _restart_host_array(RRTMGP.sw_direct_flux_dn(solver)),
        net = _restart_host_array(RRTMGP.net_flux(solver)),
        clear = _rrtmgp_clear_flux_state(solver),
        sulfate_mass,
        cloud_sampler = _rrtmgp_cloud_sampler_state(solver),
        prognostic_cloud_condensate =
            _rrtmgp_prognostic_condensate_state(radiation),
    )
end

function _restore_rrtmgp_radiation_state!(radiation::RRTMGPRadiation, state)
    solver = radiation.solver
    radiation.call_counter = state.call_counter
    _restart_copy_array!(radiation.heating_rate, state.heating_rate, "RRTMGP heating rate")
    _restart_copy_array!(
        radiation.applied_heating_rate,
        hasproperty(state, :applied_heating_rate) ?
            state.applied_heating_rate : state.heating_rate,
        "applied RRTMGP heating rate",
    )
    _restart_copy_array!(RRTMGP.lw_flux_up(solver), state.lw_up, "LW up")
    _restart_copy_array!(RRTMGP.lw_flux_dn(solver), state.lw_down, "LW down")
    _restart_copy_array!(RRTMGP.lw_flux_net(solver), state.lw_net, "LW net")
    _restart_copy_array!(RRTMGP.sw_flux_up(solver), state.sw_up, "SW up")
    _restart_copy_array!(RRTMGP.sw_flux_dn(solver), state.sw_down, "SW down")
    _restart_copy_array!(RRTMGP.sw_flux_net(solver), state.sw_net, "SW net")
    _restart_copy_array!(
        RRTMGP.sw_direct_flux_dn(solver),
        state.sw_direct_down,
        "direct SW down",
    )
    _restart_copy_array!(RRTMGP.net_flux(solver), state.net, "net flux")
    _restore_rrtmgp_clear_flux_state!(solver, state.clear)
    if !isnothing(state.sulfate_mass)
        radiation.target_aod_550nm > 0 || error(
            "checkpoint contains sulfate mass but reconstructed radiation has no aerosol forcing",
        )
        _restart_copy_array!(
            RRTMGP.aerosol_column_mass_density(solver, "sulfate"),
            state.sulfate_mass,
            "sulfate column mass",
        )
    end
    _restore_rrtmgp_cloud_sampler_state!(solver, state.cloud_sampler)
    _restore_rrtmgp_prognostic_condensate_state!(
        radiation,
        hasproperty(state, :prognostic_cloud_condensate) ?
            state.prognostic_cloud_condensate : nothing,
    )
    return radiation
end

_component_exchanger_restart_state(::Nothing) = nothing

# Exchanger states may contain lazy diagnostic operations alongside material
# fields. In particular, NumericalEarth's dynamic-ocean exchanger exposes the
# surface tracer diffusivity `κ` as a `KernelFunctionOperation`. It has no
# storage to restore and is evaluated from the restored ocean closure/state on
# demand. Serializing it via Oceananigans' generic fallback makes checkpoint
# pickup attempt `copyto!` into the operation and fail with CanonicalIndexError.
_component_exchanger_restart_entry(
    ::Oceananigans.AbstractOperations.AbstractOperation,
) = nothing
_component_exchanger_restart_entry(value) =
    Oceananigans.prognostic_state(value)

function _component_exchanger_restart_state(exchanger)
    state = exchanger.state
    names = propertynames(state)
    values = map(names) do name
        _component_exchanger_restart_entry(getproperty(state, name))
    end
    return NamedTuple{names}(values)
end

function _restore_component_exchanger_state!(exchanger, state)
    isnothing(state) && return exchanger
    isnothing(exchanger) && error(
        "checkpoint contains a component exchanger absent from the reconstructed model",
    )
    for name in keys(state)
        destination = getproperty(exchanger.state, name)
        destination isa Oceananigans.AbstractOperations.AbstractOperation &&
            continue
        Oceananigans.restore_prognostic_state!(
            destination,
            getproperty(state, name),
        )
    end
    return exchanger
end

_atmosphere_interface_restart_state(::Nothing) = nothing
function _atmosphere_interface_restart_state(interface)
    return (
        fluxes = _restart_container_state(interface.fluxes),
        temperature = Oceananigans.prognostic_state(interface.temperature),
    )
end

function _restore_atmosphere_interface_state!(interface, state)
    isnothing(state) && return interface
    isnothing(interface) && error(
        "checkpoint contains an atmosphere interface absent from the reconstructed model",
    )
    _restore_restart_container!(interface.fluxes, state.fluxes)
    Oceananigans.restore_prognostic_state!(
        interface.temperature,
        state.temperature,
    )
    return interface
end

_sea_ice_ocean_interface_restart_state(::Nothing) = nothing
function _sea_ice_ocean_interface_restart_state(interface)
    return (
        fluxes = _restart_container_state(interface.fluxes),
        temperature = Oceananigans.prognostic_state(interface.temperature),
        salinity = Oceananigans.prognostic_state(interface.salinity),
    )
end

function _restore_sea_ice_ocean_interface_state!(interface, state)
    isnothing(state) && return interface
    isnothing(interface) && error(
        "checkpoint contains a sea-ice/ocean interface absent from the reconstructed model",
    )
    _restore_restart_container!(interface.fluxes, state.fluxes)
    Oceananigans.restore_prognostic_state!(
        interface.temperature,
        state.temperature,
    )
    Oceananigans.restore_prognostic_state!(interface.salinity, state.salinity)
    return interface
end

function _component_interfaces_restart_state(interfaces)
    exchanger = interfaces.exchanger
    return (
        atmosphere_ocean_interface = _atmosphere_interface_restart_state(
            interfaces.atmosphere_ocean_interface,
        ),
        atmosphere_sea_ice_interface = _atmosphere_interface_restart_state(
            interfaces.atmosphere_sea_ice_interface,
        ),
        sea_ice_ocean_interface = _sea_ice_ocean_interface_restart_state(
            interfaces.sea_ice_ocean_interface,
        ),
        atmosphere_land_interface = _atmosphere_interface_restart_state(
            interfaces.atmosphere_land_interface,
        ),
        exchanger = (
            radiation = _component_exchanger_restart_state(exchanger.radiation),
            atmosphere = _component_exchanger_restart_state(exchanger.atmosphere),
            land = _component_exchanger_restart_state(exchanger.land),
            ocean = _component_exchanger_restart_state(exchanger.ocean),
            sea_ice = _component_exchanger_restart_state(exchanger.sea_ice),
        ),
        net_fluxes = Oceananigans.prognostic_state(interfaces.net_fluxes),
    )
end

function _restore_component_interfaces_state!(interfaces, state)
    _restore_atmosphere_interface_state!(
        interfaces.atmosphere_ocean_interface,
        state.atmosphere_ocean_interface,
    )
    _restore_atmosphere_interface_state!(
        interfaces.atmosphere_sea_ice_interface,
        state.atmosphere_sea_ice_interface,
    )
    _restore_sea_ice_ocean_interface_state!(
        interfaces.sea_ice_ocean_interface,
        state.sea_ice_ocean_interface,
    )
    _restore_atmosphere_interface_state!(
        interfaces.atmosphere_land_interface,
        state.atmosphere_land_interface,
    )
    exchanger = interfaces.exchanger
    for name in keys(state.exchanger)
        _restore_component_exchanger_state!(
            getproperty(exchanger, name),
            getproperty(state.exchanger, name),
        )
    end
    Oceananigans.restore_prognostic_state!(
        interfaces.net_fluxes,
        state.net_fluxes,
    )
    return interfaces
end

# Oceananigans deliberately calls `update_state!` once after pickup. For an
# EarthSystemModel that update is not merely halo reconciliation: it overwrites
# the lagged interface and prescribed-atmosphere forcing consumed by the next
# component step. Keep loaded host state here until that mandatory update has
# completed, then reapply it once. The entry is removed before returning, so
# every normal end-of-step update follows NumericalEarth unchanged.
const _DYNAMIC_RESTART_PENDING_STATE = IdDict{Any, Any}()
const _DYNAMIC_RESTART_PENDING_STATE_LOCK = ReentrantLock()
const _DYNAMIC_RESTART_BOUNDARY_CAPTURE_PATH = Ref{Union{Nothing, String}}(
    nothing,
)

function _request_dynamic_restart_boundary_capture!(path)
    _DYNAMIC_RESTART_BOUNDARY_CAPTURE_PATH[] =
        isnothing(path) ? nothing : abspath(path)
    return path
end

function _capture_dynamic_restart_boundary!(model)
    path = _DYNAMIC_RESTART_BOUNDARY_CAPTURE_PATH[]
    isnothing(path) && return nothing
    _DYNAMIC_RESTART_BOUNDARY_CAPTURE_PATH[] = nothing
    directory = dirname(path)
    mkpath(directory)
    temporary = tempname(directory)
    state = (; model = Oceananigans.prognostic_state(model))
    Oceananigans.OutputWriters.jldopen(temporary, "w") do file
        Oceananigans.OutputWriters.serializeproperty!(
            file,
            "simulation",
            state,
        )
    end
    mv(temporary, path; force = true)
    return path
end

function _register_pending_dynamic_restart_state!(model, state)
    lock(_DYNAMIC_RESTART_PENDING_STATE_LOCK) do
        _DYNAMIC_RESTART_PENDING_STATE[model] = state
    end
    return model
end

function _take_pending_dynamic_restart_state!(model)
    return lock(_DYNAMIC_RESTART_PENDING_STATE_LOCK) do
        pop!(_DYNAMIC_RESTART_PENDING_STATE, model, nothing)
    end
end

function _has_pending_dynamic_restart_state(model)
    return lock(_DYNAMIC_RESTART_PENDING_STATE_LOCK) do
        haskey(_DYNAMIC_RESTART_PENDING_STATE, model)
    end
end

function _mark_restored_component_simulation_initialized!(component)
    component isa Oceananigans.Simulation || return component
    component.initialized = true
    return component
end

function _mark_restored_component_simulations_initialized!(model)
    _mark_restored_component_simulation_initialized!(model.ocean)
    _mark_restored_component_simulation_initialized!(model.sea_ice)
    return model
end

_restore_raw_field_storage!(::Nothing, source) = nothing

function _restore_raw_field_storage!(
    destination::Oceananigans.Field,
    source,
)
    isnothing(source) && return destination
    destination_parent = parent(destination)
    source_data = source.data
    size(destination_parent) == size(source_data) || throw(DimensionMismatch(
        "raw restart field size $(size(source_data)) does not match destination " *
        "$(size(destination_parent))",
    ))
    copyto!(destination_parent, source_data)
    return destination
end

# `prognostic_state` represents concrete Oceananigans closure-field structs as
# nested NamedTuples.  At pickup the destination is still the concrete live
# struct, so the Tuple/NamedTuple/Field methods below are not sufficient: the
# traversal otherwise stops at (for example) `CATKEClosureFields` and leaves
# its previous-velocity history at the reconstructed zero state.  Recurse over
# the saved keys while retaining the concrete destination object.
function _restore_raw_field_storage!(destination, source::NamedTuple)
    for name in keys(source)
        hasproperty(destination, name) || throw(ArgumentError(
            "raw restart state contains field $name absent from destination " *
            "$(typeof(destination))",
        ))
        _restore_raw_field_storage!(
            getproperty(destination, name),
            getproperty(source, name),
        )
    end
    return destination
end

# Resolve the intersections between the concrete-struct traversal above and
# the pre-existing leaf/container methods explicitly.
_restore_raw_field_storage!(::Nothing, ::NamedTuple) = nothing

function _restore_raw_field_storage!(
    destination::Oceananigans.Field,
    source::NamedTuple,
)
    destination_parent = parent(destination)
    source_data = source.data
    size(destination_parent) == size(source_data) || throw(DimensionMismatch(
        "raw restart field size $(size(source_data)) does not match destination " *
        "$(size(destination_parent))",
    ))
    copyto!(destination_parent, source_data)
    return destination
end

function _restore_raw_field_storage!(destination::NamedTuple, source)
    isnothing(source) && return destination
    for name in keys(source)
        _restore_raw_field_storage!(
            getproperty(destination, name),
            getproperty(source, name),
        )
    end
    return destination
end


function _restore_raw_field_storage!(
    destination::NamedTuple,
    source::NamedTuple,
)
    keys(destination) == keys(source) || throw(ArgumentError(
        "raw restart NamedTuple keys $(keys(source)) do not match destination " *
        "$(keys(destination))",
    ))
    for name in keys(source)
        _restore_raw_field_storage!(
            getproperty(destination, name),
            getproperty(source, name),
        )
    end
    return destination
end

function _restore_raw_field_storage!(destination::Tuple, source::Tuple)
    length(destination) == length(source) || throw(DimensionMismatch(
        "raw restart tuple length $(length(source)) does not match destination " *
        "$(length(destination))",
    ))
    for index in 1:length(destination)
        _restore_raw_field_storage!(destination[index], source[index])
    end
    return destination
end

_restore_raw_field_storage!(destination, source) = destination

_raw_field_restart_state(::Nothing) = nothing

function _raw_field_restart_state(field::Oceananigans.Field)
    # `Oceananigans.prognostic_state(field)` aliases CPU storage until JLD2
    # serialization. Materialize an independent host copy here so this helper
    # is also safe for in-memory restart tests and for GPU-backed fields.
    return (; data = _restart_host_array(parent(field)))
end

function _raw_field_restart_state(fields::NamedTuple)
    names = keys(fields)
    values = map(names) do name
        _raw_field_restart_state(getproperty(fields, name))
    end
    return NamedTuple{names}(values)
end

_raw_field_restart_state(fields::Tuple) = map(_raw_field_restart_state, fields)

_adaptive_advection_timestep_state(::Nothing) = nothing

function _adaptive_advection_timestep_state(advection::NamedTuple)
    names = keys(advection)
    values = map(names) do name
        _adaptive_advection_timestep_state(getproperty(advection, name))
    end
    return NamedTuple{names}(values)
end

_adaptive_advection_timestep_state(advection::Tuple) =
    map(_adaptive_advection_timestep_state, advection)

function _adaptive_advection_timestep_state(advection)
    discretization = Oceananigans.TimeSteppers.time_discretization(advection)
    discretization isa Oceananigans.AdaptiveVerticallyImplicitDiscretization ||
        return nothing
    return discretization.Δt[]
end

_restore_adaptive_advection_timestep!(advection, ::Nothing) = advection

function _restore_adaptive_advection_timestep!(advection::NamedTuple, state)
    keys(advection) == keys(state) || throw(ArgumentError(
        "ocean advection checkpoint keys $(keys(state)) do not match " *
        "reconstructed keys $(keys(advection))",
    ))
    for name in keys(state)
        _restore_adaptive_advection_timestep!(
            getproperty(advection, name),
            getproperty(state, name),
        )
    end
    return advection
end

function _restore_adaptive_advection_timestep!(advection::Tuple, state::Tuple)
    length(advection) == length(state) || throw(DimensionMismatch(
        "ocean advection checkpoint tuple length $(length(state)) does not " *
        "match reconstructed length $(length(advection))",
    ))
    for index in eachindex(state)
        _restore_adaptive_advection_timestep!(advection[index], state[index])
    end
    return advection
end

function _restore_adaptive_advection_timestep!(advection, state)
    discretization = Oceananigans.TimeSteppers.time_discretization(advection)
    discretization isa Oceananigans.AdaptiveVerticallyImplicitDiscretization ||
        error("checkpoint contains an adaptive-advection timestep for $(typeof(advection))")
    discretization.Δt[] = state
    return advection
end

function _augment_ocean_restart_state(ocean, state)
    ocean isa Oceananigans.Simulation || return state
    timestepper = ocean.model.timestepper
    timestepper isa Oceananigans.TimeSteppers.SplitRungeKuttaTimeStepper ||
        return state
    return merge(
        state,
        (
            ; readyesm_current_tendencies = _raw_field_restart_state(
                timestepper.Gⁿ,
            ),
            readyesm_advection_timestep =
                _adaptive_advection_timestep_state(ocean.model.advection),
        ),
    )
end

function _restore_ocean_current_tendencies!(ocean, state)
    ocean isa Oceananigans.Simulation || return ocean
    hasproperty(state, :readyesm_current_tendencies) || return ocean
    timestepper = ocean.model.timestepper
    timestepper isa Oceananigans.TimeSteppers.SplitRungeKuttaTimeStepper ||
        error("checkpoint contains split-RK ocean tendencies but reconstructed ocean uses $(typeof(timestepper))")
    _restore_raw_field_storage!(
        timestepper.Gⁿ,
        state.readyesm_current_tendencies,
    )
    return ocean
end

function _restore_ocean_advection_timestep!(ocean, state)
    ocean isa Oceananigans.Simulation || return ocean
    hasproperty(state, :readyesm_advection_timestep) || return ocean
    _restore_adaptive_advection_timestep!(
        ocean.model.advection,
        state.readyesm_advection_timestep,
    )
    return ocean
end

function _restore_ocean_closure_storage!(ocean, state)
    ocean isa Oceananigans.Simulation || return ocean
    _restore_raw_field_storage!(
        ocean.model.closure_fields,
        state.model.closure_fields,
    )
    return ocean
end

function _restore_atmosphere_transform_column_storage!(atmosphere, state)
    # SpeedyWeather's Variables.copy! recurses through ScratchMemory's array
    # fields but skips its nested immutable ColumnScratchMemory object. The
    # CPU Legendre implementation reuses these buffers across transforms, so
    # pickup must restore them explicitly to preserve the exact call boundary.
    destination = atmosphere.variables.scratch.transform_memory.column
    source = state.variables.scratch.transform_memory.column
    _restart_copy_array!(
        destination.north,
        source.north,
        "atmosphere transform column north",
    )
    _restart_copy_array!(
        destination.south,
        source.south,
        "atmosphere transform column south",
    )
    return atmosphere
end

function Oceananigans.prognostic_state(
    model::_ReadyESMRadiativeSpeedyEarthSystem,
)
    state = invoke(
        Oceananigans.prognostic_state,
        Tuple{NumericalEarth.EarthSystemModel},
        model,
    )
    ocean = _augment_ocean_restart_state(model.ocean, state.ocean)
    return merge(
        state,
        (
            ; ocean,
            interfaces = _component_interfaces_restart_state(model.interfaces),
            balanced_launch = _balanced_launch_restart_state(model),
        ),
    )
end

function Oceananigans.restore_prognostic_state!(
    model::_ReadyESMRadiativeSpeedyEarthSystem,
    state,
)
    invoke(
        Oceananigans.restore_prognostic_state!,
        Tuple{NumericalEarth.EarthSystemModel, Any},
        model,
        state,
    )
    _restore_ocean_current_tendencies!(model.ocean, state.ocean)
    _restore_ocean_advection_timestep!(model.ocean, state.ocean)
    _restore_ocean_closure_storage!(model.ocean, state.ocean)
    _restore_atmosphere_transform_column_storage!(
        model.atmosphere,
        state.atmosphere,
    )
    _restore_balanced_launch_restart_state!(
        model,
        hasproperty(state, :balanced_launch) ? state.balanced_launch : nothing,
    )
    # Ocean and sea ice are nested Oceananigans `Simulation`s. Upstream's
    # Simulation checkpoint state deliberately omits the lifecycle flag
    # `initialized`; a reconstructed nested simulation therefore defaults to
    # false and would perform an extra initialization/update on its first
    # resumed component step. At any valid coupled checkpoint both components
    # have already crossed initialization, so restore that hidden lifecycle
    # state explicitly.
    _mark_restored_component_simulations_initialized!(model)
    if hasproperty(state, :interfaces) && !isnothing(state.interfaces)
        _register_pending_dynamic_restart_state!(model, state)
    end
    return model
end

Oceananigans.restore_prognostic_state!(
    model::_ReadyESMRadiativeSpeedyEarthSystem,
    ::Nothing,
) = model

function NumericalEarth.EarthSystemModels.update_state!(
    model::_ReadyESMRadiativeSpeedyEarthSystem,
    callbacks = [],
)
    pending_restart = _has_pending_dynamic_restart_state(model)
    pending_restart || _apply_balanced_launch_if_due!(model)
    invoke(
        NumericalEarth.EarthSystemModels.update_state!,
        Tuple{NumericalEarth.EarthSystemModel, Any},
        model,
        callbacks,
    )
    pending = _take_pending_dynamic_restart_state!(model)
    isnothing(pending) && return nothing

    # `update_state!` is mandatory in Oceananigans' pickup initialization, but
    # for this coupled model it is not an idempotent halo-only operation. It
    # also recomputes ocean closure fields and vertical geometry, sea-ice
    # auxiliaries, component states, and lagged coupling buffers. Reapplying
    # only the atmosphere/radiation/interfaces leaves a different state at the
    # start of the first resumed step. Restore the complete checkpoint graph
    # once after that native reconciliation; saved Field state includes halos,
    # so this reproduces the exact end-of-segment state without another update.
    invoke(
        Oceananigans.restore_prognostic_state!,
        Tuple{NumericalEarth.EarthSystemModel, Any},
        model,
        pending,
    )
    _restore_ocean_current_tendencies!(model.ocean, pending.ocean)
    _restore_ocean_advection_timestep!(model.ocean, pending.ocean)
    _restore_ocean_closure_storage!(model.ocean, pending.ocean)
    _restore_atmosphere_transform_column_storage!(
        model.atmosphere,
        pending.atmosphere,
    )
    _mark_restored_component_simulations_initialized!(model)
    _restore_component_interfaces_state!(model.interfaces, pending.interfaces)
    _capture_dynamic_restart_boundary!(model)
    return nothing
end

# Public internal helpers used by the fresh-process split-run verifier. Keeping
# the file payload equal to `prognostic_state(simulation)` means it exercises
# the same state graph as Oceananigans.Checkpointer without serializing model
# construction objects.
function save_dynamic_restart_state(simulation, path::AbstractString)
    directory = dirname(abspath(path))
    mkpath(directory)
    temporary = tempname(directory)
    state = Oceananigans.prognostic_state(simulation)
    Oceananigans.OutputWriters.jldopen(temporary, "w") do file
        Oceananigans.OutputWriters.serializeproperty!(file, "simulation", state)
    end
    mv(temporary, path; force = true)
    return abspath(path)
end

function restore_dynamic_restart_state!(simulation, path::AbstractString)
    state = Oceananigans.OutputWriters.load_checkpoint_state(
        abspath(path);
        base_path = "simulation",
    )
    Oceananigans.restore_prognostic_state!(simulation, state)
    # Perform the same mandatory post-pickup update as Oceananigans.run!. The
    # ReadyESM-specialized update reconciles native state first, then consumes
    # and reapplies the saved lagged interface state exactly once.
    Oceananigans.TimeSteppers.update_state!(simulation.model)
    simulation.initialized = true
    return simulation
end

function _surface_temperature_callback_state(callback)
    count = callback.timestep_counter
    0 <= count <= length(callback.temperature) || error(
        "invalid global surface-temperature callback counter $count",
    )
    return (
        timestep_counter = count,
        temperature = copy(@view callback.temperature[1:count]),
    )
end

function _restore_surface_temperature_callback!(callback, state)
    length(callback.temperature) >= length(state.temperature) || error(
        "reconstructed surface-temperature callback horizon is shorter than checkpoint history",
    )
    callback.timestep_counter = state.timestep_counter
    copyto!(@view(callback.temperature[1:length(state.temperature)]), state.temperature)
    return callback
end

const _RADIATION_BUDGET_HISTORY_FIELDS = (
    :incoming_shortwave,
    :outgoing_shortwave,
    :outgoing_longwave,
    :clear_outgoing_shortwave,
    :clear_outgoing_longwave,
    :net_downward,
    :clear_net_downward,
)

function _radiation_budget_callback_state(callback::GlobalRadiationBudgetCallback)
    count = callback.timestep_counter
    histories = map(_RADIATION_BUDGET_HISTORY_FIELDS) do name
        values = getproperty(callback, name)
        0 <= count <= length(values) || error(
            "invalid radiation-budget callback counter $count for $name",
        )
        copy(@view values[1:count])
    end
    return (
        timestep_counter = count,
        histories = NamedTuple{_RADIATION_BUDGET_HISTORY_FIELDS}(histories),
    )
end

function _restore_radiation_budget_callback!(callback, state)
    callback.timestep_counter = state.timestep_counter
    for name in _RADIATION_BUDGET_HISTORY_FIELDS
        destination = getproperty(callback, name)
        source = getproperty(state.histories, name)
        length(destination) >= length(source) || error(
            "reconstructed radiation-budget horizon for $name is shorter than checkpoint history",
        )
        copyto!(@view(destination[1:length(source)]), source)
    end
    return callback
end

const _ATMOSPHERE_HISTORY_FIELDS = (
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
    :cloud_condensate_water_path,
    :total_water_budget_residual,
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
)

function _atmosphere_diagnostics_callback_state(
    callback::GlobalAtmosphereDiagnosticsCallback,
)
    histories = map(_ATMOSPHERE_HISTORY_FIELDS) do name
        copy(getproperty(callback, name))
    end
    return (
        timestep_counter = callback.timestep_counter,
        sample_every_n_steps = callback.sample_every_n_steps,
        final_timestep = callback.final_timestep,
        timestep_seconds = callback.timestep_seconds,
        surface_sigma_thickness = callback.surface_sigma_thickness,
        histories = NamedTuple{_ATMOSPHERE_HISTORY_FIELDS}(histories),
        layer_process_histories = NamedTuple{
            _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS,
        }(map(
            name -> copy(getproperty(callback, name)),
            _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS,
        )),
        layer_process_accumulators = NamedTuple{
            _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS,
        }(map(
            name -> _restart_host_array(getproperty(callback, name)),
            _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS,
        )),
        process_diagnostic_start_timestep =
            callback.process_diagnostic_start_timestep,
        column_cumulative_surface_water_vapor = _restart_host_array(
            callback.column_cumulative_surface_water_vapor,
        ),
        column_cumulative_precipitation = _restart_host_array(
            callback.column_cumulative_precipitation,
        ),
        column_cumulative_land_precipitation = _restart_host_array(
            callback.column_cumulative_land_precipitation,
        ),
        column_cumulative_land_evapotranspiration = _restart_host_array(
            callback.column_cumulative_land_evapotranspiration,
        ),
        column_cumulative_land_surface_runoff = _restart_host_array(
            callback.column_cumulative_land_surface_runoff,
        ),
        balanced_launch_atmosphere_water_source = Float64(
            callback.balanced_launch_atmosphere_water_source,
        ),
        balanced_launch_land_water_source = Float64(
            callback.balanced_launch_land_water_source,
        ),
    )
end

function _restore_atmosphere_diagnostics_callback!(callback, state)
    callback.timestep_counter = state.timestep_counter
    callback.sample_every_n_steps = state.sample_every_n_steps
    callback.final_timestep = state.final_timestep
    if hasproperty(state, :timestep_seconds)
        state.timestep_seconds == callback.timestep_seconds || error(
            "atmosphere diagnostic timestep changed across restart: " *
            "checkpoint=$(state.timestep_seconds) reconstructed=$(callback.timestep_seconds)",
        )
    end
    if hasproperty(state, :surface_sigma_thickness)
        state.surface_sigma_thickness == callback.surface_sigma_thickness || error(
            "atmosphere diagnostic surface sigma thickness changed across restart: " *
            "checkpoint=$(state.surface_sigma_thickness) " *
            "reconstructed=$(callback.surface_sigma_thickness)",
        )
    end
    for name in _ATMOSPHERE_HISTORY_FIELDS
        source = if hasproperty(state.histories, name)
            getproperty(state.histories, name)
        elseif name in (
            :cumulative_balanced_launch_atmosphere_water_source,
            :cumulative_balanced_launch_land_water_source,
            :cloud_condensate_water_path,
        )
            zeros(
                eltype(getproperty(callback, name)),
                length(state.histories.time_days),
            )
        elseif name == :total_water_budget_residual
            copy(state.histories.water_budget_residual)
        else
            error("checkpoint lacks atmosphere diagnostic history $name")
        end
        _restart_resize_copy!(
            getproperty(callback, name),
            source,
            "atmosphere diagnostic $name",
        )
    end
    if hasproperty(state, :layer_process_histories) &&
       hasproperty(state, :layer_process_accumulators)
        for name in _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS
            _restart_resize_copy!(
                getproperty(callback, name),
                getproperty(state.layer_process_histories, name),
                "atmosphere process-profile diagnostic $name",
            )
        end
        for name in _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS
            _restart_copy_array!(
                getproperty(callback, name),
                getproperty(state.layer_process_accumulators, name),
                "atmosphere process-profile accumulator $name",
            )
        end
        callback.process_diagnostic_start_timestep = hasproperty(
            state,
            :process_diagnostic_start_timestep,
        ) ? state.process_diagnostic_start_timestep : 0
    else
        # Checkpoints written before the cumulative process observer cannot
        # reconstruct past tendencies. If the checkpoint already coincides
        # with a retained diagnostic sample, establish the zero there. For an
        # off-cycle checkpoint, retain all existing samples as missing and
        # begin at the next regular sample; this avoids assigning a zero to an
        # earlier timestamp or a partial process interval to the next one.
        nsamples = length(state.histories.time_days)
        nlayers = size(
            callback.weighted_cumulative_convective_humidity_change,
            2,
        )
        checkpoint_day = state.timestep_counter *
            Float64(callback.timestep_seconds) / 86_400
        # The callback stores `NF(checkpoint_day)` itself. Compare that exact
        # representable value: a many-ULP tolerance would eventually confuse
        # adjacent 15-minute steps once a Float32 day counter reaches a
        # centennial magnitude.
        boundary_is_sampled = nsamples > 0 &&
            last(state.histories.time_days) ==
                eltype(state.histories.time_days)(checkpoint_day)
        next_sample_timestep = boundary_is_sampled ?
            state.timestep_counter :
            (div(state.timestep_counter, state.sample_every_n_steps) + 1) *
                state.sample_every_n_steps
        for name in _ATMOSPHERE_PROCESS_PROFILE_HISTORY_FIELDS
            history = getproperty(callback, name)
            resize!(history, nsamples * nlayers)
            fill!(history, eltype(history)(NaN))
            if boundary_is_sampled
                fill!(@view(history[((nsamples - 1) * nlayers + 1):end]), 0)
            end
        end
        for name in _ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS
            fill!(getproperty(callback, name), 0)
        end
        callback.process_diagnostic_start_timestep = next_sample_timestep
    end
    _restart_copy_array!(
        callback.column_cumulative_surface_water_vapor,
        state.column_cumulative_surface_water_vapor,
        "column cumulative surface water vapor",
    )
    _restart_copy_array!(
        callback.column_cumulative_precipitation,
        state.column_cumulative_precipitation,
        "column cumulative precipitation",
    )
    _restart_copy_array!(
        callback.column_cumulative_land_precipitation,
        state.column_cumulative_land_precipitation,
        "column cumulative land precipitation",
    )
    _restart_copy_array!(
        callback.column_cumulative_land_evapotranspiration,
        state.column_cumulative_land_evapotranspiration,
        "column cumulative land evapotranspiration",
    )
    _restart_copy_array!(
        callback.column_cumulative_land_surface_runoff,
        state.column_cumulative_land_surface_runoff,
        "column cumulative land surface runoff",
    )
    callback.balanced_launch_atmosphere_water_source =
        hasproperty(state, :balanced_launch_atmosphere_water_source) ?
        state.balanced_launch_atmosphere_water_source : 0
    callback.balanced_launch_land_water_source =
        hasproperty(state, :balanced_launch_land_water_source) ?
        state.balanced_launch_land_water_source : 0
    return callback
end

function _speedy_callback_state(callbacks)
    return (
        rrtmgp_postphysics_timestep_counter =
            callbacks[:rrtmgp_postphysics_update].timestep_counter,
        surface_temperature = _surface_temperature_callback_state(
            callbacks[:global_surface_temperature],
        ),
        radiation_budget = _radiation_budget_callback_state(
            callbacks[:global_radiation_budget],
        ),
        atmosphere_diagnostics = _atmosphere_diagnostics_callback_state(
            callbacks[:global_atmosphere_diagnostics],
        ),
    )
end

function _restore_speedy_callback_state!(callbacks, state)
    callbacks[:rrtmgp_postphysics_update].timestep_counter =
        state.rrtmgp_postphysics_timestep_counter
    _restore_surface_temperature_callback!(
        callbacks[:global_surface_temperature],
        state.surface_temperature,
    )
    _restore_radiation_budget_callback!(
        callbacks[:global_radiation_budget],
        state.radiation_budget,
    )
    _restore_atmosphere_diagnostics_callback!(
        callbacks[:global_atmosphere_diagnostics],
        state.atmosphere_diagnostics,
    )
    return callbacks
end

const _SPEEDY_IMPLICIT_RESTART_ARRAY_FIELDS = (
    :temp_profile,
    :R,
    :U,
    :L,
    :W,
    :L0,
    :L1,
    :L2,
    :L3,
    :L4,
    :S,
    :S⁻¹,
)

function _speedy_implicit_restart_state(implicit)
    arrays = NamedTuple{_SPEEDY_IMPLICIT_RESTART_ARRAY_FIELDS}(
        map(_SPEEDY_IMPLICIT_RESTART_ARRAY_FIELDS) do name
            _restart_host_array(getproperty(implicit, name))
        end,
    )
    return merge(
        (
            initialized = implicit.initialized,
            ξ = implicit.ξ[],
        ),
        arrays,
    )
end

function _restore_speedy_implicit_restart_state!(implicit, state)
    implicit.initialized = state.initialized
    implicit.ξ[] = state.ξ
    for name in _SPEEDY_IMPLICIT_RESTART_ARRAY_FIELDS
        _restart_copy_array!(
            getproperty(implicit, name),
            getproperty(state, name),
            "SpeedyWeather implicit $name",
        )
    end
    return implicit
end

function _speedy_time_stepping_restart_state(time_stepping)
    return (; first_step_euler = time_stepping.first_step_euler)
end

function _restore_speedy_time_stepping_restart_state!(time_stepping, state)
    time_stepping.first_step_euler = state.first_step_euler
    return time_stepping
end

function Oceananigans.prognostic_state(simulation::SpeedyWeather.Simulation)
    radiation = simulation.model.longwave_radiation
    radiation isa RRTMGPRadiation || error(
        "ReadyESM dynamic restart requires RRTMGPRadiation",
    )
    variables = deepcopy(
        SpeedyWeather.on_architecture(SpeedyWeather.CPU(), simulation.variables),
    )
    return (
        variables,
        time_stepping = _speedy_time_stepping_restart_state(
            simulation.model.time_stepping,
        ),
        implicit = _speedy_implicit_restart_state(simulation.model.implicit),
        radiation = _rrtmgp_radiation_state(radiation),
        callbacks = _speedy_callback_state(simulation.model.callbacks),
    )
end

function Oceananigans.restore_prognostic_state!(
    simulation::SpeedyWeather.Simulation,
    state,
)
    # JLD2 reconstructs the checkpointed Variables tree on the CPU.
    # SpeedyWeather's alias-aware copy handles ordinary host-to-device arrays.
    # The narrowly typed lower-triangular, RingGrids.Field, and ScratchMemory
    # methods above supply its unsupported wrapper transfers without applying
    # architecture conversion to the attached Terrarium state.
    copy!(simulation.variables, state.variables)
    if hasproperty(state, :time_stepping)
        _restore_speedy_time_stepping_restart_state!(
            simulation.model.time_stepping,
            state.time_stepping,
        )
    end
    if hasproperty(state, :implicit)
        _restore_speedy_implicit_restart_state!(
            simulation.model.implicit,
            state.implicit,
        )
    end
    radiation = simulation.model.longwave_radiation
    radiation isa RRTMGPRadiation || error(
        "ReadyESM dynamic restart requires RRTMGPRadiation",
    )
    _restore_rrtmgp_radiation_state!(radiation, state.radiation)
    _restore_speedy_callback_state!(simulation.model.callbacks, state.callbacks)
    return simulation
end

Oceananigans.restore_prognostic_state!(
    simulation::SpeedyWeather.Simulation,
    ::Nothing,
) = simulation

const _COUPLED_BUDGET_STATIC_FIELDS = (
    :state_integrals,
    :flux_integrals,
    :diagnostic_fields,
)

function _coupled_budget_callback_state(callback::GlobalCoupledBudgetCallback)
    names = Tuple(
        name for name in fieldnames(typeof(callback))
        if name ∉ _COUPLED_BUDGET_STATIC_FIELDS
    )
    values = map(names) do name
        value = getproperty(callback, name)
        value isa AbstractArray ? _restart_host_array(value) : deepcopy(value)
    end
    return NamedTuple{names}(values)
end

function _restore_coupled_budget_callback!(callback, state)
    for name in keys(state)
        destination = getproperty(callback, name)
        source = getproperty(state, name)
        if destination isa Vector
            _restart_resize_copy!(destination, source, "coupled budget $name")
        elseif destination isa AbstractArray
            _restart_copy_array!(destination, source, "coupled budget $name")
        else
            setproperty!(callback, name, source)
        end
    end
    sample_count = length(callback.time_days)
    if !hasproperty(state, :balanced_launch_applied_fraction)
        resize!(callback.balanced_launch_applied_fraction, sample_count)
        fill!(callback.balanced_launch_applied_fraction, 0)
    end
    if !hasproperty(state, :balanced_launch_cumulative_sources)
        resize!(
            callback.balanced_launch_cumulative_sources,
            sample_count * length(_BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS),
        )
        fill!(callback.balanced_launch_cumulative_sources, 0)
    end
    return callback
end

function Oceananigans.prognostic_state(
    callback::Oceananigans.Callback{P, F, S, CS},
) where {P, F <: GlobalCoupledBudgetCallback, S, CS}
    return (
        schedule = Oceananigans.prognostic_state(callback.schedule),
        func = _coupled_budget_callback_state(callback.func),
    )
end

function Oceananigans.restore_prognostic_state!(
    callback::Oceananigans.Callback{P, F, S, CS},
    state,
) where {P, F <: GlobalCoupledBudgetCallback, S, CS}
    Oceananigans.restore_prognostic_state!(callback.schedule, state.schedule)
    _restore_coupled_budget_callback!(callback.func, state.func)
    return callback
end

Oceananigans.restore_prognostic_state!(
    callback::Oceananigans.Callback{P, F, S, CS},
    ::Nothing,
) where {P, F <: GlobalCoupledBudgetCallback, S, CS} = callback

"""
    install_dynamic_checkpointer!(simulation; schedule, dir, prefix,
                                  overwrite_existing=false, cleanup=true)

Install an Oceananigans JLD2 checkpointer whose state methods cover the complete
ReadyESM dynamic stack. Resume with ordinary `Oceananigans.run!(simulation;
pickup=path)` after reconstructing the same configuration and final run horizon.
Oceananigans performs its native pickup `update_state!` before the first resumed
timestep; ReadyESM then reapplies the checkpointed lagged coupling state once.
"""
function install_dynamic_checkpointer!(
    simulation;
    schedule,
    dir,
    prefix = "readiesm_checkpoint",
    overwrite_existing = false,
    cleanup = true,
)
    checkpointer = Oceananigans.Checkpointer(
        simulation.model;
        schedule,
        dir,
        prefix,
        overwrite_existing,
        cleanup,
    )
    simulation.output_writers[:checkpointer] = checkpointer
    return checkpointer
end
