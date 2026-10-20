# BS08 Stage B, routing increment: the runoff enthalpy held in the land exchange follows the native water through
# the consume, the scatter to the ocean receivers and the ocean freshwater heat content (contract [state] and
# [ocean_receiver]).
#
# The option is carried by the exchange storage. With `terrarium_runoff_enthalpy = false` the exchange is the native
# `(; pending_depth)` NamedTuple, so `TerrariumRunoffLand` has its native type and only the native methods run. With
# the option on, the exchange carries the enthalpy storage below, and the methods in this file wrap the native ones
# with `invoke`. Water arrays, routing weights and the native ocean assembly are never changed.

"Exchange storage with the BS08 enthalpy state (option on)."
const _BS08Exchange = NamedTuple{(:pending_depth, :pending_enthalpy, :held_enthalpy_rate, :receiver_enthalpy_flux,
                                  :freshwater_density, :guard)}
const _BS08RunoffLand = TerrariumRunoffLand{<:Any, <:Any, <:Any, <:Any, <:_BS08Exchange}

# `_build_runoff_routing` multiplies every weight by this density (its default); the receiver enthalpy divides by it.
const BS08_ROUTING_FRESHWATER_DENSITY = 1000.0

"""Mutable guard of the ocean receiver correction. `state` is `:idle`, `:assembled` (the native assembly has just
written Jh) or `:applied` (the correction is in Jh). The correction runs only from `:assembled`, so it can never be
added twice to the same native Jh. `iteration` and `time` record the last application."""
mutable struct BS08ReceiverGuard
    state::Symbol
    iteration::Int
    time::Float64
    applications::Int
end
BS08ReceiverGuard() = BS08ReceiverGuard(:idle, -1, NaN, 0)

"Native exchange storage, or the BS08 storage when the Terrarium state owns runoff enthalpy."
function _runoff_exchange_storage(terrarium_state, pending_depth, surface_runoff, ocean)
    _bs08_owns_enthalpy(terrarium_state) || return (; pending_depth)
    grid = ocean.model.grid
    grid isa Oceananigans.ImmersedBoundaries.MutableGridOfSomeKind || error(
        "terrarium_runoff_enthalpy requires the mutable z-star ocean grid; otherwise NumericalEarth ignores Jh for T")
    pending_enthalpy = Oceananigans.interior(terrarium_state.runoff_pending_enthalpy)
    size(pending_enthalpy) == size(pending_depth) || error("runoff enthalpy and water exchange shapes differ")
    held_enthalpy_rate = similar(surface_runoff)
    fill!(held_enthalpy_rate, 0)
    receiver_enthalpy_flux = Oceananigans.Field{Oceananigans.Center, Oceananigans.Center, Nothing}(grid)
    fill!(receiver_enthalpy_flux, 0)
    return (; pending_depth, pending_enthalpy, held_enthalpy_rate, receiver_enthalpy_flux,
            freshwater_density = BS08_ROUTING_FRESHWATER_DENSITY, guard = BS08ReceiverGuard())
end

# ------------------------------------------------------------------------------------------------------------------
# Consume: runs before the native consume and repeats its water arithmetic exactly, so the delivered water depth is
# the one the native kernel will deliver. The pending enthalpy splits at the pending specific enthalpy between the
# delivered water and the remainder (Stage A smaller/larger partner rule, bs08_consume).
# ------------------------------------------------------------------------------------------------------------------
KernelAbstractions.@kernel function _bs08_consume_runoff_enthalpy_kernel!(held_enthalpy_rate, pending_enthalpy,
                                                                          pending_depth, water_rate_prototype, dt)
    column = @index(Global, Linear)
    amount = Float64(@inbounds pending_depth[column])
    enthalpy = Float64(@inbounds pending_enthalpy[column])
    if isfinite(amount) && amount >= 0 && isfinite(enthalpy)
        value = convert(eltype(water_rate_prototype), amount / dt)       # the native rate, bit for bit
        if Float64(value) * dt > amount
            value = prevfloat(value)
        end
        delivered = Float64(value) * dt
        if iszero(amount)
            # no water: no enthalpy leaves; any pending enthalpy stays for the consistency check
            @inbounds held_enthalpy_rate[column] = zero(eltype(held_enthalpy_rate))
        else
            sent, kept = bs08_consume(amount, enthalpy, delivered)
            @inbounds begin
                held_enthalpy_rate[column] = convert(eltype(held_enthalpy_rate), sent / dt)
                pending_enthalpy[column] = kept
            end
        end
    else
        @inbounds held_enthalpy_rate[column] = convert(eltype(held_enthalpy_rate), NaN)
    end
end

function _bs08_consume_runoff_enthalpy!(land::_BS08RunoffLand, Δt)
    dt = Float64(Δt)
    isfinite(dt) && dt > 0 || throw(ArgumentError("runoff exchange requires a finite positive interval"))
    ex = land.runoff_exchange
    backend = KernelAbstractions.get_backend(ex.held_enthalpy_rate)
    _bs08_consume_runoff_enthalpy_kernel!(backend)(ex.held_enthalpy_rate, ex.pending_enthalpy, ex.pending_depth,
        land.surface_runoff, dt; ndrange = length(ex.held_enthalpy_rate))
    KernelAbstractions.synchronize(backend)
    return nothing
end

function Oceananigans.TimeSteppers.time_step!(land::_BS08RunoffLand, Δt)
    _bs08_consume_runoff_enthalpy!(land, Δt)
    return invoke(Oceananigans.TimeSteppers.time_step!, Tuple{TerrariumRunoffLand, Any}, land, Δt)
end

function Oceananigans.Simulations.reset_clock!(land::_BS08RunoffLand)
    invoke(Oceananigans.Simulations.reset_clock!, Tuple{TerrariumRunoffLand}, land)
    ex = land.runoff_exchange
    fill!(ex.pending_enthalpy, 0)
    fill!(ex.held_enthalpy_rate, 0)
    fill!(ex.receiver_enthalpy_flux, 0)
    guard = ex.guard
    guard.state = :idle
    guard.iteration = -1
    guard.time = NaN
    guard.applications = 0
    return nothing
end

# ------------------------------------------------------------------------------------------------------------------
# Scatter: the held enthalpy rate (W m-2 of source land) reaches each receiver with the native contribution weight
# divided by the freshwater density in that weight, so a receiver gets W m-2 of ocean cell.
# ------------------------------------------------------------------------------------------------------------------
KernelAbstractions.@kernel function _bs08_scatter_runoff_enthalpy_kernel!(receiver_enthalpy_flux, held_enthalpy_rate,
        contribution_column, contribution_weight, target_i, target_j, offsets, freshwater_density)
    target = @index(Global, Linear)
    accumulated = zero(eltype(receiver_enthalpy_flux))
    @inbounds for n in offsets[target]:(offsets[target + 1] - 1)
        column = contribution_column[n]
        accumulated += contribution_weight[n] * held_enthalpy_rate[column, 1, 1] / freshwater_density
    end
    @inbounds receiver_enthalpy_flux[target_i[target], target_j[target], 1] = accumulated
end

function _bs08_scatter_runoff_enthalpy!(land::_BS08RunoffLand)
    flux = land.runoff_exchange.receiver_enthalpy_flux
    fill!(flux, 0)
    land.target_cells == 0 && return nothing
    backend = KernelAbstractions.get_backend(parent(flux.data))
    _bs08_scatter_runoff_enthalpy_kernel!(backend)(flux.data, land.runoff_exchange.held_enthalpy_rate,
        land.contribution_column, land.contribution_weight, land.target_i, land.target_j, land.offsets,
        land.runoff_exchange.freshwater_density; ndrange = land.target_cells)
    Oceananigans.fill_halo_regions!(flux)
    Oceananigans.Architectures.synchronize(Oceananigans.architecture(flux.grid))
    return nothing
end

function NumericalEarth.EarthSystemModels.interpolate_state!(exchanger, grid, land::_BS08RunoffLand, coupled_model)
    invoke(NumericalEarth.EarthSystemModels.interpolate_state!, Tuple{Any, Any, TerrariumRunoffLand, Any},
        exchanger, grid, land, coupled_model)
    _bs08_scatter_runoff_enthalpy!(land)
    return nothing
end

# ------------------------------------------------------------------------------------------------------------------
# Ocean receiver: after NumericalEarth's native assembly, on receiver cells only,
#   Jh += H_land / (rho_ocean cp_ocean) - T_oc J_land / rho_ocean
# (bs08_land_heat_content_correction). Only receiver cells are touched, so every other Jh value keeps its bits.
# ------------------------------------------------------------------------------------------------------------------
KernelAbstractions.@kernel function _bs08_receiver_correction_kernel!(heat_content, receiver_enthalpy_flux,
        land_freshwater_flux, ocean_surface_temperature, target_i, target_j, grid, rho_ocean, cp_ocean)
    target = @index(Global, Linear)
    i = @inbounds target_i[target]
    j = @inbounds target_j[target]
    kᴺ = size(grid, 3)
    inactive = Oceananigans.Grids.inactive_node(i, j, kᴺ, grid, Oceananigans.Center(), Oceananigans.Center(),
                                                Oceananigans.Center())
    if !inactive
        @inbounds heat_content[i, j, 1] += bs08_land_heat_content_correction(receiver_enthalpy_flux[i, j, 1],
            land_freshwater_flux[i, j, 1], ocean_surface_temperature[i, j, 1], rho_ocean, cp_ocean)
    end
end

"""Apply the receiver correction to Jh (`heat_content` is the field's data array). Errors before any change unless the
native assembly has just written Jh (guard state `:assembled`)."""
function _bs08_apply_receiver_enthalpy!(heat_content, land::_BS08RunoffLand, land_freshwater_flux,
                                        ocean_surface_temperature, grid, rho_ocean, cp_ocean, clock)
    guard = land.runoff_exchange.guard
    guard.state === :assembled || error(
        "BS08 receiver correction requested without a fresh native assembly (guard state $(guard.state), last " *
        "application at iteration $(guard.iteration), time $(guard.time)); Jh is unchanged")
    if land.target_cells > 0
        backend = KernelAbstractions.get_backend(parent(heat_content))
        _bs08_receiver_correction_kernel!(backend)(heat_content, land.runoff_exchange.receiver_enthalpy_flux.data,
            land_freshwater_flux, ocean_surface_temperature, land.target_i, land.target_j, grid,
            rho_ocean, cp_ocean; ndrange = land.target_cells)
        KernelAbstractions.synchronize(backend)
    end
    guard.state = :applied
    guard.iteration = clock.iteration
    guard.time = Float64(clock.time)
    guard.applications += 1
    return nothing
end

const _BS08RunoffEarthSystem = NumericalEarth.EarthSystemModel{<:Any, <:Any, <:_BS08RunoffLand}

function NumericalEarth.EarthSystemModels.update_net_fluxes!(coupled_model::_BS08RunoffEarthSystem,
                                                            ocean::NumericalEarth.Oceans.OceananigansModelSimulations)
    invoke(NumericalEarth.EarthSystemModels.update_net_fluxes!,
        Tuple{Any, NumericalEarth.Oceans.OceananigansModelSimulations}, coupled_model, ocean)
    land = coupled_model.land
    land.runoff_exchange.guard.state = :assembled
    properties = coupled_model.interfaces.ocean_properties
    _bs08_apply_receiver_enthalpy!(coupled_model.interfaces.net_fluxes.ocean.freshwater_heat_content.data,
        land, coupled_model.interfaces.exchanger.land.state.freshwater_flux.data,
        NumericalEarth.EarthSystemModels.ocean_surface_temperature(ocean), ocean.model.grid,
        properties.reference_density, properties.heat_capacity, coupled_model.clock)
    return nothing
end

# ------------------------------------------------------------------------------------------------------------------
# Checkpoint: schema v2 adds the pending enthalpy and the held enthalpy rate. An option-on restore from a v1
# checkpoint zero-initialises both (the native implicit state: all routed water at 0 C) and records that it did so.
# ------------------------------------------------------------------------------------------------------------------
const BS08_EXCHANGE_LEGACY_ZERO_INITIALIZATIONS = Ref(0)

function Oceananigans.prognostic_state(land::_BS08RunoffLand)
    ex = land.runoff_exchange
    # v3 = v2 plus the host-side BS08 diagnostic history (bs08_diagnostics.jl), or nothing when none is installed.
    return (; schema = "readiesm_native_runoff_exchange_v3",
            runoff_exchange = _native_runoff_exchange_state(land),
            runoff_enthalpy = (pending_enthalpy_jm2 = Array(ex.pending_enthalpy),
                               held_enthalpy_rate_wm2 = Array(ex.held_enthalpy_rate)),
            diagnostic_history = _bs08_history_checkpoint(land))
end

function Oceananigans.restore_prognostic_state!(land::_BS08RunoffLand, state::NamedTuple)
    hasproperty(state, :schema) || error("unsupported native-runoff checkpoint schema")
    ex = land.runoff_exchange
    if state.schema == "readiesm_native_runoff_exchange_v1"
        _restore_native_runoff_exchange!(land, state.runoff_exchange)
        fill!(ex.pending_enthalpy, 0)
        fill!(ex.held_enthalpy_rate, 0)
        legacy_history = _bs08_diagnostic_history(land)
        isnothing(legacy_history) || _reset_bs08_history!(legacy_history)
        BS08_EXCHANGE_LEGACY_ZERO_INITIALIZATIONS[] += 1
        @info "BS08: option-on restore from a runoff exchange without enthalpy; zero-initialised the pending enthalpy and the held enthalpy rate"
        return land
    end
    state.schema in ("readiesm_native_runoff_exchange_v2", "readiesm_native_runoff_exchange_v3") ||
        error("unsupported native-runoff checkpoint schema")
    pairs = ((ex.pending_enthalpy, state.runoff_enthalpy.pending_enthalpy_jm2),
             (ex.held_enthalpy_rate, state.runoff_enthalpy.held_enthalpy_rate_wm2))
    for (destination, source) in pairs                 # validate everything before any array changes
        size(destination) == size(source) || throw(DimensionMismatch("runoff enthalpy checkpoint shape differs"))
        all(isfinite, source) || error("runoff enthalpy checkpoint contains a non-finite value")
    end
    # a checkpoint must not introduce energy without a water carrier
    depth = state.runoff_exchange.pending_depth_m
    size(depth) == size(ex.pending_enthalpy) || throw(DimensionMismatch("runoff exchange checkpoint shape differs"))
    any(c -> iszero(depth[c]) && !iszero(state.runoff_enthalpy.pending_enthalpy_jm2[c]), eachindex(depth)) && error(
        "runoff enthalpy checkpoint has pending enthalpy without pending water")
    # A v2 checkpoint, or a v3 one saved without history, returns an installed history to launch state. The history
    # then refuses its first sample unless the run is at launch (documented launch-only legacy policy).
    history = _bs08_diagnostic_history(land)
    saved_history = state.schema == "readiesm_native_runoff_exchange_v3" ? state.diagnostic_history : nothing
    restore_history = !isnothing(history) && !isnothing(saved_history)
    restore_history && _validate_bs08_history_checkpoint(history, saved_history)
    _restore_native_runoff_exchange!(land, state.runoff_exchange)
    for (destination, source) in pairs
        copyto!(destination, source)
    end
    if restore_history
        _restore_bs08_history!(history, saved_history)
    elseif !isnothing(history)
        _reset_bs08_history!(history)
    end
    return land
end
