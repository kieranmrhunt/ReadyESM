"""
    ConservativeSeaIceAdvection(scheme; maximum_conditional_thickness=15)

Use `scheme` to transport sea-ice concentration and the per-cell ice and snow
volumes `ℵ * h` and `ℵ * hs`. The volume products are materialized in scalar
work fields before every tendency evaluation and reconstructed directly at
faces. ClimaSeaIce 0.5.8 otherwise transports the conditional thicknesses
independently from concentration; reconstructing `ℵ` and `h` independently and
multiplying their face values is not equivalent to reconstructing `ℵ * h`.
Where transport leaves finite volume at vanishing concentration, increase
concentration just enough to keep the conditional thickness below
`maximum_conditional_thickness`. This subgrid consolidation preserves ice and
snow volume exactly; if the area-equivalent volume itself exceeds the bound,
concentration saturates at one and the real excess thickness is retained.
"""
struct ConservativeSeaIceAdvection{A, FT, IV, SV}
    scheme::A
    maximum_conditional_thickness::FT
    ice_volume::IV
    snow_volume::SV
end

const _CONSERVATIVE_SEA_ICE_ADVECTION_PROVENANCE =
    "conservative_A_direct_Ah_direct_Ahs_2D_face_flux_Az_normalized_WENO5_refreshed_halos_hmax_v4"

@kernel function _conditional_ice_from_effective_thickness_kernel!(
    ice_thickness,
    concentration,
    maximum_conditional_thickness,
)
    i, j = @index(Global, NTuple)
    k = 1

    @inbounds begin
        # ECCO SIheff is volume per unit grid area. ClimaSeaIce's `h` is the
        # thickness within the ice-covered fraction and computes mass from
        # `A * h`, so assigning SIheff directly to `h` weights the ECCO volume
        # by concentration a second time. Reconstruct the conditional state
        # while retaining SIheff itself as the represented volume.
        effective_thickness = ice_thickness[i, j, k]
        h, A = _bounded_conditional_ice_state(
            effective_thickness,
            concentration[i, j, k],
            maximum_conditional_thickness,
        )
        ice_thickness[i, j, k] = h
        concentration[i, j, k] = A
    end
end

"""
    _initialize_conditional_ice_from_effective_thickness!(model)

Convert an ECCO effective-thickness field already assigned to
`model.ice_thickness` into ClimaSeaIce's conditional thickness. The conversion
preserves the ECCO area-equivalent ice volume exactly. Extremely small or zero
reported concentration is consolidated only enough to keep the conditional
representation finite; genuine area-equivalent thickness above the configured
representation bound remains visible.
"""
function _initialize_conditional_ice_fields_from_effective_thickness!(
    ice_thickness,
    concentration,
    grid,
    maximum_conditional_thickness,
)
    architecture = Oceananigans.Architectures.architecture(grid)
    launch_grid = grid isa Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ?
        grid.underlying_grid : grid
    hmax = convert(eltype(grid), maximum_conditional_thickness)
    isfinite(hmax) && hmax > 0 || throw(ArgumentError(
        "maximum conditional sea-ice thickness must be finite and positive",
    ))
    Oceananigans.Utils.launch!(
        architecture,
        launch_grid,
        :xy,
        _conditional_ice_from_effective_thickness_kernel!,
        ice_thickness,
        concentration,
        hmax,
    )
    Oceananigans.Architectures.synchronize(architecture)
    return nothing
end

function _initialize_conditional_ice_from_effective_thickness!(model)
    maximum_conditional_thickness = if model.advection isa
                                       ConservativeSeaIceAdvection
        model.advection.maximum_conditional_thickness
    else
        eltype(model.grid)(15)
    end
    return _initialize_conditional_ice_fields_from_effective_thickness!(
        model.ice_thickness,
        model.ice_concentration,
        model.grid,
        maximum_conditional_thickness,
    )
end

function ConservativeSeaIceAdvection(
    scheme;
    maximum_conditional_thickness = 15,
)
    hmax = float(maximum_conditional_thickness)
    isfinite(hmax) && hmax > 0 || throw(ArgumentError(
        "maximum_conditional_thickness must be finite and positive",
    ))
    return ConservativeSeaIceAdvection(scheme, hmax, nothing, nothing)
end

function Oceananigans.Advection.materialize_advection(
    advection::ConservativeSeaIceAdvection,
    grid,
)
    scheme = Oceananigans.Advection.materialize_advection(advection.scheme, grid)
    ice_volume = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    snow_volume = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    return ConservativeSeaIceAdvection(
        scheme,
        eltype(grid)(advection.maximum_conditional_thickness),
        ice_volume,
        snow_volume,
    )
end

@inline function _bounded_conditional_ice_state(
    volume,
    concentration,
    maximum_conditional_thickness,
)
    zero_volume = zero(volume)
    V = max(zero_volume, volume)
    A = clamp(concentration, zero(concentration), one(concentration))
    required_concentration = min(
        one(A),
        V / maximum_conditional_thickness,
    )
    A = max(A, required_concentration)

    # `floatmin` also protects the extremely small subnormal-volume case from
    # a finite-volume divide by zero. The selected value is discarded when
    # V == 0, but GPU `ifelse` evaluates both branches.
    safe_A = max(A, floatmin(typeof(A)))
    h = V / safe_A
    valid_ice = V > zero_volume
    return (
        ifelse(valid_ice, h, zero_volume),
        ifelse(valid_ice, A, zero(A)),
    )
end

@inline function _ice_volume_flux_x(i, j, k, grid, scheme, u, volume)
    return Oceananigans.Advection._advective_tracer_flux_x(
        i, j, k, grid, scheme, u, volume,
    )
end


@inline function _ice_volume_flux_y(i, j, k, grid, scheme, v, volume)
    return Oceananigans.Advection._advective_tracer_flux_y(
        i, j, k, grid, scheme, v, volume,
    )
end


const _ImmersedBoundaryGrid =
    Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid

@inline function _ice_volume_flux_x(i, j, k, grid::_ImmersedBoundaryGrid, scheme, u, volume)
    flux = _ice_volume_flux_x(
        i, j, k, grid.underlying_grid, scheme, u, volume,
    )
    return Oceananigans.Advection.conditional_flux_fcc(
        i, j, k, grid, zero(grid), flux,
    )
end


@inline function _ice_volume_flux_y(i, j, k, grid::_ImmersedBoundaryGrid, scheme, v, volume)
    flux = _ice_volume_flux_y(
        i, j, k, grid.underlying_grid, scheme, v, volume,
    )
    return Oceananigans.Advection.conditional_flux_cfc(
        i, j, k, grid, zero(grid), flux,
    )
end

"""
    _sea_ice_flux_2D_x/y(i, j, k, grid, flux, args...)

Convert Oceananigans' three-dimensional tracer face flux to the horizontal
flux for a two-dimensional sea-ice inventory. `Ax` and `Ay` contain the
instantaneous z-star face height; dividing it out before the horizontal
divergence makes the tendency invariant to that ocean-cell height.
"""
@inline function _sea_ice_flux_2D_x(i, j, k, grid, flux, args...)
    Δz = Oceananigans.Operators.Δzᶠᶜᶜ(i, j, k, grid)
    return ifelse(Δz > 0, flux(i, j, k, grid, args...) / Δz, zero(grid))
end

@inline function _sea_ice_flux_2D_y(i, j, k, grid, flux, args...)
    Δz = Oceananigans.Operators.Δzᶜᶠᶜ(i, j, k, grid)
    return ifelse(Δz > 0, flux(i, j, k, grid, args...) / Δz, zero(grid))
end

@inline _conservative_div_UV(i, j, k, grid, ::Nothing, velocities, volume) =
    zero(eltype(grid))

@inline function _conservative_div_UV(i, j, k, grid, scheme, velocities, volume)
    δx_flux = Oceananigans.Operators.δxᶜᵃᵃ(
        i, j, k, grid, _sea_ice_flux_2D_x, _ice_volume_flux_x,
        scheme, velocities.u, volume,
    )
    δy_flux = Oceananigans.Operators.δyᵃᶜᵃ(
        i, j, k, grid, _sea_ice_flux_2D_y, _ice_volume_flux_y,
        scheme, velocities.v, volume,
    )
    # Sea ice is a two-dimensional inventory: `volume` is area-equivalent
    # thickness (m), so its horizontal flux divergence is normalized by the
    # horizontal cell area. Vᶜᶜᶜ adds the ocean-cell height and would make this
    # tendency vary as z-star changes even when the ice state and fluxes do not.
    return (δx_flux + δy_flux) /
           Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)
end


@inline function _conservative_div_UV(
    i,
    j,
    k,
    grid,
    scheme::Oceananigans.Advection.FluxFormAdvection,
    velocities,
    volume,
)
    δx_flux = Oceananigans.Operators.δxᶜᵃᵃ(
        i, j, k, grid, _sea_ice_flux_2D_x, _ice_volume_flux_x,
        scheme.x, velocities.u, volume,
    )
    δy_flux = Oceananigans.Operators.δyᵃᶜᵃ(
        i, j, k, grid, _sea_ice_flux_2D_y, _ice_volume_flux_y,
        scheme.y, velocities.v, volume,
    )
    return (δx_flux + δy_flux) /
           Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)
end

@inline _conservative_div_Uc(i, j, k, grid, ::Nothing, velocities, concentration) =
    zero(eltype(grid))

@inline function _conservative_div_Uc(i, j, k, grid, scheme, velocities, concentration)
    δx_flux = Oceananigans.Operators.δxᶜᵃᵃ(
        i, j, k, grid, _sea_ice_flux_2D_x, _ice_volume_flux_x,
        scheme, velocities.u, concentration,
    )
    δy_flux = Oceananigans.Operators.δyᵃᶜᵃ(
        i, j, k, grid, _sea_ice_flux_2D_y, _ice_volume_flux_y,
        scheme, velocities.v, concentration,
    )
    return (δx_flux + δy_flux) /
           Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)
end

@inline function _conservative_div_Uc(
    i,
    j,
    k,
    grid,
    scheme::Oceananigans.Advection.FluxFormAdvection,
    velocities,
    concentration,
)
    δx_flux = Oceananigans.Operators.δxᶜᵃᵃ(
        i, j, k, grid, _sea_ice_flux_2D_x, _ice_volume_flux_x,
        scheme.x, velocities.u, concentration,
    )
    δy_flux = Oceananigans.Operators.δyᵃᶜᵃ(
        i, j, k, grid, _sea_ice_flux_2D_y, _ice_volume_flux_y,
        scheme.y, velocities.v, concentration,
    )
    return (δx_flux + δy_flux) /
           Oceananigans.Operators.Azᶜᶜᶜ(i, j, k, grid)
end

@inline function _conservative_snow_volume_tendency!(
    i,
    j,
    k,
    G,
    grid,
    scheme,
    velocities,
    snow_volume,
)
    @inbounds G.hs[i, j, 1] = -_conservative_div_UV(
        i,
        j,
        k,
        grid,
        scheme,
        velocities,
        snow_volume,
    )
    return nothing
end

@inline _conservative_snow_volume_tendency!(
    i,
    j,
    k,
    G,
    grid,
    scheme,
    velocities,
    ::Nothing,
) = nothing

@inline function _store_snow_volume!(
    i,
    j,
    snow_volume,
    active,
    represented_concentration,
    snow_thickness,
)
    @inbounds represented_snow_thickness = ifelse(
        active,
        snow_thickness[i, j, 1],
        zero(snow_thickness[i, j, 1]),
    )
    @inbounds snow_volume[i, j, 1] =
        represented_concentration * represented_snow_thickness
    return nothing
end

@inline _store_snow_volume!(
    i,
    j,
    ::Nothing,
    active,
    represented_concentration,
    ::Nothing,
) = nothing

@kernel function _materialize_conservative_ice_volumes!(
    ice_volume,
    snow_volume,
    grid,
    concentration,
    ice_thickness,
    snow_thickness,
)
    i, j = @index(Global, NTuple)
    surface_k = size(grid, 3)
    active = !Oceananigans.Grids.inactive_cell(i, j, surface_k, grid)
    @inbounds represented_concentration = ifelse(
        active,
        concentration[i, j, 1],
        zero(concentration[i, j, 1]),
    )
    @inbounds represented_ice_thickness = ifelse(
        active,
        ice_thickness[i, j, 1],
        zero(ice_thickness[i, j, 1]),
    )
    @inbounds ice_volume[i, j, 1] =
        represented_concentration * represented_ice_thickness
    _store_snow_volume!(
        i,
        j,
        snow_volume,
        active,
        represented_concentration,
        snow_thickness,
    )
end

@kernel function _compute_conservative_ice_tendencies!(
    G,
    grid,
    velocities,
    scheme,
    ice_volume,
    concentration,
    snow_volume,
)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    # G.h stores the conservative per-cell ice-volume tendency during the
    # dynamic step. It is converted back to conditional thickness below.
    @inbounds begin
        G.h[i, j, 1] = -_conservative_div_UV(
            i,
            j,
            k,
            grid,
            scheme,
            velocities,
            ice_volume,
        )
        G.ℵ[i, j, 1] = -_conservative_div_Uc(
            i,
            j,
            k,
            grid,
            scheme,
            velocities,
            concentration,
        )
    end

    _conservative_snow_volume_tendency!(
        i,
        j,
        k,
        G,
        grid,
        scheme,
        velocities,
        snow_volume,
    )
end

function ClimaSeaIce.compute_tracer_tendencies!(
    model::ClimaSeaIce.SeaIceModel{
        GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP,
        MFX, STF, A,
    },
) where {
    GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP,
    MFX, STF, A <: ConservativeSeaIceAdvection,
}
    grid = model.grid
    arch = Oceananigans.Architectures.architecture(grid)
    ice_volume = model.advection.ice_volume
    snow_volume = isnothing(model.snow_thickness) ?
        nothing : model.advection.snow_volume
    launch_grid = grid isa _ImmersedBoundaryGrid ? grid.underlying_grid : grid

    Oceananigans.Utils.launch!(
        arch,
        launch_grid,
        :xy,
        _materialize_conservative_ice_volumes!,
        ice_volume,
        snow_volume,
        grid,
        model.ice_concentration,
        model.ice_thickness,
        model.snow_thickness,
    )

    fields = isnothing(model.snow_thickness) ?
        (
            model.velocities.u,
            model.velocities.v,
            model.ice_thickness,
            model.ice_concentration,
            ice_volume,
        ) :
        (
            model.velocities.u,
            model.velocities.v,
            model.ice_thickness,
            model.ice_concentration,
            model.snow_thickness,
            ice_volume,
            snow_volume,
        )
    Oceananigans.fill_halo_regions!(fields; async = true)

    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xy,
        _compute_conservative_ice_tendencies!,
        model.timestepper.Gⁿ,
        grid,
        model.velocities,
        model.advection.scheme,
        ice_volume,
        model.ice_concentration,
        snow_volume,
    )
    return nothing
end

@kernel function _update_conservative_ice_state!(
    h,
    concentration,
    h_base,
    concentration_base,
    snow_thickness,
    snow_thickness_base,
    G,
    Δt,
    maximum_conditional_thickness,
)
    i, j = @index(Global, NTuple)
    k = 1

    @inbounds begin
        A₀ = concentration_base[i, j, k]
        V₀ = A₀ * h_base[i, j, k]
        V = max(zero(V₀), V₀ + Δt * G.h[i, j, k])
        A_advected = A₀ + Δt * G.ℵ[i, j, k]
        h_new, A_new = _bounded_conditional_ice_state(
            V,
            A_advected,
            maximum_conditional_thickness,
        )
        concentration[i, j, k] = A_new
        h[i, j, k] = h_new
    end

    _update_conservative_snow_state!(
        i,
        j,
        k,
        concentration,
        snow_thickness,
        snow_thickness_base,
        concentration_base,
        G,
        Δt,
    )
end

@inline _update_conservative_snow_state!(
    i,
    j,
    k,
    concentration,
    ::Nothing,
    ::Nothing,
    concentration_base,
    G,
    Δt,
) = nothing

@inline function _update_conservative_snow_state!(
    i,
    j,
    k,
    concentration,
    snow_thickness,
    snow_thickness_base,
    concentration_base,
    G,
    Δt,
)
    @inbounds begin
        A₀ = concentration_base[i, j, k]
        Vs₀ = A₀ * snow_thickness_base[i, j, k]
        Vs = max(zero(Vs₀), Vs₀ + Δt * G.hs[i, j, k])
        A = concentration[i, j, k]
        snow_thickness[i, j, k] = ifelse(
            A > zero(A) && Vs > zero(Vs),
            Vs / A,
            zero(Vs),
        )
    end
    return nothing
end

@inline _rebase_snow_after_consolidation!(
    i,
    j,
    k,
    ::Nothing,
    snow_volume,
    concentration,
) = nothing

@inline function _rebase_snow_after_consolidation!(
    i,
    j,
    k,
    snow_thickness,
    snow_volume,
    concentration,
)
    safe_A = max(concentration, floatmin(typeof(concentration)))
    snow_thickness[i, j, k] = ifelse(
        concentration > zero(concentration) && snow_volume > zero(snow_volume),
        snow_volume / safe_A,
        zero(snow_volume),
    )
    return nothing
end

@kernel function _consolidate_conditional_ice_state_kernel!(
    h,
    concentration,
    snow_thickness,
    maximum_conditional_thickness,
)
    i, j = @index(Global, NTuple)
    k = 1

    @inbounds begin
        A_old = concentration[i, j, k]
        V = A_old * h[i, j, k]
        Vs = isnothing(snow_thickness) ? zero(V) :
             A_old * snow_thickness[i, j, k]
        h_new, A_new = _bounded_conditional_ice_state(
            V,
            A_old,
            maximum_conditional_thickness,
        )
        h[i, j, k] = h_new
        concentration[i, j, k] = A_new
        _rebase_snow_after_consolidation!(
            i,
            j,
            k,
            snow_thickness,
            Vs,
            A_new,
        )
    end
end

function _consolidate_conditional_ice_state!(model)
    grid = model.grid
    arch = Oceananigans.Architectures.architecture(grid)
    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xy,
        _consolidate_conditional_ice_state_kernel!,
        model.ice_thickness,
        model.ice_concentration,
        model.snow_thickness,
        model.advection.maximum_conditional_thickness,
    )
    return nothing
end

const _SplitRK = Oceananigans.TimeSteppers.SplitRungeKuttaTimeStepper
const _ForwardEuler = ClimaSeaIce.ForwardEulerTimeStepper

function ClimaSeaIce.dynamic_time_step!(
    model::ClimaSeaIce.SeaIceModel{
        GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP,
        MFX, STF, A,
    },
    Δt,
) where {
    GR, TD, SNT, D, TS <: _SplitRK, CL, U, T, IT, IC, SNH, ID, SND, PT,
    CT, SP, MFX, STF, A <: ConservativeSeaIceAdvection,
}
    h = model.ice_thickness
    concentration = model.ice_concentration
    h_base = model.timestepper.Ψ⁻.h
    concentration_base = model.timestepper.Ψ⁻.ℵ
    snow_thickness = model.snow_thickness
    snow_thickness_base = isnothing(snow_thickness) ? nothing : model.timestepper.Ψ⁻.hs
    grid = model.grid
    arch = Oceananigans.Architectures.architecture(grid)

    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xy,
        _update_conservative_ice_state!,
        h,
        concentration,
        h_base,
        concentration_base,
        snow_thickness,
        snow_thickness_base,
        model.timestepper.Gⁿ,
        Δt,
        model.advection.maximum_conditional_thickness,
    )
    return nothing
end

# ClimaSeaIce's SplitRK sequence applies thermodynamics after the dynamic
# update. Thermodynamics may change concentration again and reconstruct h as
# V/A, so consolidate once more before the substep returns to momentum,
# coupling, or the next RK stage. The four upstream calls and their order are
# retained; only the final volume-preserving representation repair is added.
# Two read-only calls of the opt-in sea-ice branch observer
# (sea_ice_branch_observer.jl) bracket the thermodynamic update. With no
# observer installed each returns after testing an empty registry. Installed,
# each copies the final-stage state to host memory and writes no model field.
function Oceananigans.TimeSteppers.rk_substep!(
    model::ClimaSeaIce.SeaIceModel{
        GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP,
        MFX, STF, A,
    },
    Δτ,
    callbacks,
) where {
    GR, TD, SNT, D, TS <: _SplitRK, CL, U, T, IT, IC, SNH, ID, SND, PT,
    CT, SP, MFX, STF, A <: ConservativeSeaIceAdvection,
}
    ClimaSeaIce.compute_tendencies!(model, Δτ)
    ClimaSeaIce.SeaIceDynamics.time_step_momentum!(
        model,
        model.dynamics,
        Δτ,
    )
    _observe_sea_ice_momentum_endpoint!(model, Δτ)
    ClimaSeaIce.dynamic_time_step!(model, Δτ)
    _observe_sea_ice_branch_before_thermodynamics!(model, Δτ)
    ClimaSeaIce.SeaIceThermodynamics.thermodynamic_time_step!(
        model,
        model.ice_thermodynamics,
        model.snow_thermodynamics,
        Δτ,
    )
    _observe_sea_ice_branch_after_thermodynamics!(model, Δτ)
    _consolidate_conditional_ice_state!(model)
    return nothing
end

function ClimaSeaIce.dynamic_time_step!(
    model::ClimaSeaIce.SeaIceModel{
        GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP,
        MFX, STF, A,
    },
    Δt,
) where {
    GR, TD, SNT, D, TS <: _ForwardEuler, CL, U, T, IT, IC, SNH, ID, SND,
    PT, CT, SP, MFX, STF, A <: ConservativeSeaIceAdvection,
}
    h = model.ice_thickness
    concentration = model.ice_concentration
    snow_thickness = model.snow_thickness
    grid = model.grid
    arch = Oceananigans.Architectures.architecture(grid)

    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xy,
        _update_conservative_ice_state!,
        h,
        concentration,
        h,
        concentration,
        snow_thickness,
        snow_thickness,
        model.timestepper.Gⁿ,
        Δt,
        model.advection.maximum_conditional_thickness,
    )
    return nothing
end
