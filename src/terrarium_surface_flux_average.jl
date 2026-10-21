"Scratch flux integrals for one land coupling interval, on the land device."
function _terrarium_surface_flux_accumulator(state)
    sensible = Oceananigans.interior(state.sensible_heat_flux)
    latent = Oceananigans.interior(state.latent_heat_flux)
    size(sensible) == size(latent) || throw(DimensionMismatch("Land surface fluxes"))
    return (; sensible_heat_flux=zero.(sensible), latent_heat_flux=zero.(latent))
end

_reset_terrarium_surface_fluxes!(::Nothing) = nothing
function _reset_terrarium_surface_fluxes!(fluxes)
    for buffer in values(fluxes)
        fill!(buffer, zero(eltype(buffer)))
    end
    return nothing
end

_accumulate_terrarium_surface_fluxes!(::Nothing, state, dt) = nothing
function _accumulate_terrarium_surface_fluxes!(fluxes, state, dt)
    # With Forward Euler and finalize=false these are the fluxes used by
    # the completed substep. Final auxiliary evaluation overwrites them.
    for name in keys(fluxes)
        buffer = getproperty(fluxes, name)
        source = Oceananigans.interior(getproperty(state, name))
        size(buffer) == size(source) || throw(DimensionMismatch("Land surface flux buffer"))
        step = eltype(buffer)(dt)
        buffer .+= step .* source
    end
    return nothing
end

_average_terrarium_surface_fluxes!(::Nothing, duration) = nothing
function _average_terrarium_surface_fluxes!(fluxes, duration)
    isfinite(duration) && duration > 0 || throw(ArgumentError("Positive coupling duration required"))
    for buffer in values(fluxes)
        buffer ./= eltype(buffer)(duration)
    end
    return nothing
end
