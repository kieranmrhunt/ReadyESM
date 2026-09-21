const TERRARIUM_RUNOFF_TIME_INTEGRATION_PROVENANCE =
    "native_step_runoff_amount_held_for_next_ocean_interval_v1"

# This amount belongs to the land-router exchange, not to a diagnostic
# callback. Do not reset it at atmospheric half-steps or auxiliary refreshes.
KernelAbstractions.@kernel function _accumulate_native_runoff_kernel!(pending, runoff, dt)
    column = @index(Global, Linear)
    @inbounds pending[column] += dt * runoff[column]
end

function _accumulate_native_runoff!(state, dt)
    hasproperty(state, :runoff_pending_depth) || return nothing
    pending = Oceananigans.interior(state.runoff_pending_depth)
    runoff = Oceananigans.interior(state.surface_runoff)
    backend = KernelAbstractions.get_backend(pending)
    _accumulate_native_runoff_kernel!(backend)(pending, runoff, dt;
        ndrange=length(pending))
    return nothing
end

KernelAbstractions.@kernel function _consume_native_runoff_kernel!(rate, pending, dt)
    column = @index(Global, Linear)
    amount = Float64(@inbounds pending[column])
    if isfinite(amount) && amount >= 0
        value = convert(eltype(rate), amount / dt)
        # Keep any representational remainder in the exchange store. A rate
        # rounded upward must not export more water than the land produced.
        if Float64(value) * dt > amount
            value = prevfloat(value)
        end
        @inbounds begin
            rate[column] = value
            pending[column] = amount - Float64(value) * dt
        end
    else
        @inbounds rate[column] = convert(eltype(rate), NaN)
    end
end

function _consume_native_runoff!(rate, pending, timestep_seconds)
    dt = Float64(timestep_seconds)
    isfinite(dt) && dt > 0 || throw(ArgumentError(
        "runoff exchange requires a finite positive interval"))
    size(rate) == size(pending) || throw(DimensionMismatch(
        "runoff exchange rate and amount shapes differ"))
    backend = KernelAbstractions.get_backend(rate)
    _consume_native_runoff_kernel!(backend)(rate, pending, dt;ndrange=length(rate))
    KernelAbstractions.synchronize(backend)
    return nothing
end

function _native_runoff_exchange_state(land)
    return (
        pending_depth_m = Array(land.runoff_exchange.pending_depth),
        held_rate_ms = Array(land.surface_runoff),
    )
end

function _restore_native_runoff_exchange!(land, state)
    for (destination, source) in (
        (land.runoff_exchange.pending_depth, state.pending_depth_m),
        (land.surface_runoff, state.held_rate_ms),
    )
        size(destination) == size(source) || throw(DimensionMismatch(
            "runoff exchange checkpoint shape differs"))
        all(value -> isfinite(value) && value >= 0, source) || error(
            "runoff exchange checkpoint contains invalid water")
    end
    copyto!(land.runoff_exchange.pending_depth, state.pending_depth_m)
    copyto!(land.surface_runoff, state.held_rate_ms)
    return nothing
end
