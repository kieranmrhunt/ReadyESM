"""
Build a compact, fully coupled SpeedyWeather–RRTMGP–Oceananigans–ClimaSeaIce
simulation. NumericalEarth (re-exported by ClimaOcean) owns the component
exchange and ConservativeRegridding remaps atmosphere and ocean/ice fields.

The idealized grid remains available for coupling tests. Production experiments
can select a one-degree latitude--longitude or global tripolar ocean with
realistic bathymetry, active momentum/tracer transport, and dynamic/advected sea ice.
Ocean temperature and salinity are in degrees Celsius and practical salinity.
"""

const _TROPICAL_SEA_ICE_MAXIMUM_ABSOLUTE_LATITUDE_DEGREES = 30.0

"""
Two-band penetrating shortwave with bathymetry-aware coastal optical water type.

NumericalEarth's default is Jerlov Type I everywhere, appropriate for clear
open-ocean water but not for shallow, sediment- and chlorophyll-rich shelves.
Because its conservative Beer-law implementation absorbs the beam remaining at
the seabed in the lowest wet cell, using Type I in a 10--40 m one-degree column
can concentrate 7--27% of the deep band there. This scheme uses Paulson &
Simpson's Jerlov Type III coefficients for columns at or shallower than 50 m,
the existing Type I coefficients at or deeper than 200 m, and a smooth linear
blend between. The original conservative bottom absorption is retained, but
the physically inappropriate clear-water residual in shelf cells is removed
without changing bathymetry or losing solar energy from the coupled budget.
"""
struct DepthAwareCoastalRadiation{FT, J, B}
    type_i_fraction::FT
    type_i_first_absorption_coefficient::FT
    type_i_second_absorption_coefficient::FT
    type_iii_fraction::FT
    type_iii_first_absorption_coefficient::FT
    type_iii_second_absorption_coefficient::FT
    coastal_depth_m::FT
    open_ocean_depth_m::FT
    surface_flux::J
    bottom_height::B
end

function Adapt.adapt_structure(to, radiation::DepthAwareCoastalRadiation)
    return DepthAwareCoastalRadiation(
        Adapt.adapt(to, radiation.type_i_fraction),
        Adapt.adapt(to, radiation.type_i_first_absorption_coefficient),
        Adapt.adapt(to, radiation.type_i_second_absorption_coefficient),
        Adapt.adapt(to, radiation.type_iii_fraction),
        Adapt.adapt(to, radiation.type_iii_first_absorption_coefficient),
        Adapt.adapt(to, radiation.type_iii_second_absorption_coefficient),
        Adapt.adapt(to, radiation.coastal_depth_m),
        Adapt.adapt(to, radiation.open_ocean_depth_m),
        Adapt.adapt(to, radiation.surface_flux),
        Adapt.adapt(to, radiation.bottom_height),
    )
end

function DepthAwareCoastalRadiation(grid, bottom_height)
    FT = eltype(grid)
    surface_flux = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    return DepthAwareCoastalRadiation(
        FT(0.58),
        inv(FT(0.35)),
        inv(FT(23)),
        FT(0.78),
        inv(FT(1.4)),
        inv(FT(7.9)),
        FT(50),
        FT(200),
        surface_flux,
        bottom_height,
    )
end

@inline function _coastal_optical_weight(depth_m, coastal_depth_m, open_ocean_depth_m)
    return clamp(
        (open_ocean_depth_m - depth_m) /
        (open_ocean_depth_m - coastal_depth_m),
        zero(depth_m),
        one(depth_m),
    )
end

@inline function _depth_aware_optical_parameters(radiation, i, j)
    bottom = @inbounds radiation.bottom_height[i, j, 1]
    depth_m = max(zero(bottom), -bottom)
    coastal_weight = _coastal_optical_weight(
        depth_m,
        radiation.coastal_depth_m,
        radiation.open_ocean_depth_m,
    )
    open_weight = one(coastal_weight) - coastal_weight
    fraction = open_weight * radiation.type_i_fraction +
               coastal_weight * radiation.type_iii_fraction
    first_absorption_coefficient =
        open_weight * radiation.type_i_first_absorption_coefficient +
        coastal_weight * radiation.type_iii_first_absorption_coefficient
    second_absorption_coefficient =
        open_weight * radiation.type_i_second_absorption_coefficient +
        coastal_weight * radiation.type_iii_second_absorption_coefficient
    return fraction, first_absorption_coefficient, second_absorption_coefficient
end

@inline function _depth_aware_beers_law_radiation(i, j, k, grid, surface_flux, κ)
    center = Oceananigans.Center()
    face = Oceananigans.Face()
    Nz = size(grid, 3)
    z = Oceananigans.Grids.znode(i, j, k, grid, center, center, face)
    η = Oceananigans.Grids.znode(i, j, Nz + 1, grid, center, center, face)
    beam = surface_flux * exp(κ * (z - η))
    return ifelse(
        Oceananigans.Grids.inactive_cell(i, j, k - 1, grid),
        zero(beam),
        beam,
    )
end

@inline function (radiation::DepthAwareCoastalRadiation)(i, j, k, grid, clock, fields)
    surface_flux = @inbounds radiation.surface_flux[i, j, 1]
    fraction, κ₁, κ₂ = _depth_aware_optical_parameters(radiation, i, j)
    dJ₁dz = Oceananigans.Operators.∂zᶜᶜᶜ(
        i, j, k, grid, _depth_aware_beers_law_radiation, surface_flux, κ₁,
    )
    dJ₂dz = Oceananigans.Operators.∂zᶜᶜᶜ(
        i, j, k, grid, _depth_aware_beers_law_radiation, surface_flux, κ₂,
    )
    return fraction * dJ₁dz + (one(fraction) - fraction) * dJ₂dz
end

function NumericalEarth.Oceans.compute_radiative_forcing!(
    radiation::DepthAwareCoastalRadiation,
    downwelling_shortwave_radiation,
    coupled_model,
)
    density = coupled_model.interfaces.ocean_properties.reference_density
    heat_capacity = coupled_model.interfaces.ocean_properties.heat_capacity
    parent(radiation.surface_flux) .=
        -parent(downwelling_shortwave_radiation) ./ (density * heat_capacity)
    return nothing
end

function NumericalEarth.Oceans.get_radiative_forcing(
    radiation::DepthAwareCoastalRadiation,
)
    return radiation
end

@inline function NumericalEarth.Oceans.shortwave_radiative_forcing(
    i,
    j,
    grid,
    radiation::DepthAwareCoastalRadiation,
    transmitted_shortwave,
    ocean_properties,
)
    density = ocean_properties.reference_density
    heat_capacity = ocean_properties.heat_capacity
    @inbounds radiation.surface_flux[i, j, 1] =
        -transmitted_shortwave / (density * heat_capacity)
    return zero(transmitted_shortwave)
end

function _ocean_radiative_forcing(config, grid, bottom_height)
    if config.ocean_shortwave_scheme == :depth_aware_coastal
        return DepthAwareCoastalRadiation(grid, bottom_height)
    end
    return NumericalEarth.Oceans.default_radiative_forcing(grid)
end

function _idealized_ocean_advection(::Type{FT}) where {FT}
    momentum_advection = Oceananigans.WENOVectorInvariant(
        FT;
        order = 5,
        time_discretization = Oceananigans.AdaptiveVerticallyImplicitDiscretization(
            FT;
            cfl = FT(0.5),
        ),
    )
    tracer_advection = Oceananigans.WENO(
        FT;
        order = 5,
        time_discretization = Oceananigans.AdaptiveVerticallyImplicitDiscretization(
            FT;
            cfl = FT(0.5),
        ),
    )
    return momentum_advection, tracer_advection
end

function _production_ocean_advection(::Type{FT}) where {FT}
    # Match ClimaOcean's documented latitude_longitude_ocean_sea_ice example,
    # which is configured for a two-year integration at a 20-minute coupled
    # step: centered vector-invariant momentum and seventh-order WENO tracers.
    # Earlier ReadyESM experiments overrode this with WENO momentum. Even with
    # adaptive vertical implicitness that branch became non-finite in 1--3
    # simulated days, for both analytic and ECCO hydrography.
    momentum_advection = Oceananigans.VectorInvariant(FT)
    # Keep the documented horizontal WENO7 scheme, but use NumericalEarth's
    # current adaptive vertical discretization. It remains explicit below the
    # local vertical-CFL threshold and switches only the vertical tracer flux
    # to the supported implicit treatment during short-lived column spikes.
    tracer_advection = Oceananigans.WENO(
        FT;
        order = 7,
        time_discretization = Oceananigans.AdaptiveVerticallyImplicitDiscretization(
            FT;
            cfl = FT(0.5),
        ),
    )
    return momentum_advection, tracer_advection
end

"""Combine two spatial masks without capturing non-isbits closures on GPU."""
struct _MaximumSpatialMask{M1, M2}
    first::M1
    second::M2
end

@inline function (mask::_MaximumSpatialMask)(x, y, z)
    return max(mask.first(x, y, z), mask.second(x, y, z))
end

function _latitude_longitude_polar_sponge_mask(
    ::Type{FT},
    start_latitude_degrees,
    stop_latitude_degrees,
) where {FT}
    start = FT(start_latitude_degrees)
    stop = FT(stop_latitude_degrees)
    north = Oceananigans.CosineRampMask{:y}(
        start = start,
        stop = stop,
    )
    south = Oceananigans.CosineRampMask{:y}(
        start = -start,
        stop = -stop,
    )
    return _MaximumSpatialMask(north, south)
end

function _production_latitude_longitude_closure(::Type{FT}) where {FT}
    # The package's simplified closure was designed for memory-limited test
    # GPUs. In the coupled year it left coastal salt/freshwater anomalies
    # horizontally trapped and eventually allowed the artificial polar-wall
    # circulation to grow. Use the documented one-degree production closure:
    # CATKE vertical mixing, Gent--McWilliams/isopycnal mixing, biharmonic
    # momentum diffusion, and background vertical diffusivity.
    return ClimaOcean.OceanConfigurations.default_one_degree_closure(
        κ_skew = FT(500),
        κ_symmetric = FT(200),
        biharmonic_timescale = FT(10 * 86_400),
        background_ν = FT(1e-5),
    )
end

"""
Static diffusivity fields used by the opt-in unresolved-estuary closure.

The fields are created with the ocean model, then populated after the
Terrarium runoff routing map is available. This keeps the mixing footprint
identical to the actual runoff receivers rather than applying a global coastal
coefficient. The vertical coefficient is integrated implicitly; the
horizontal coefficient uses Oceananigans' conservative flux-form diffusion.
"""
struct RiverMouthMixingFields{V, H, EV, EH, A, F}
    vertical_diffusivity::V
    horizontal_diffusivity::H
    eligible_vertical_diffusivity::EV
    eligible_horizontal_diffusivity::EH
    active_receiver_mask::A
    dynamically_gated::Bool
    reference_freshwater_mass_flux_kgm2s::F
end

function _river_mouth_mixing_closure(config::ExperimentConfig, grid)
    config.ocean_river_mouth_mixing == :none && return (
        closures = (),
        fields = nothing,
    )
    _is_localized_estuary_mixing(config.ocean_river_mouth_mixing) || error(
        "unsupported river-mouth mixing mode $(config.ocean_river_mouth_mixing)",
    )

    FT = eltype(grid)
    vertical_diffusivity = Oceananigans.CenterField(grid)
    horizontal_diffusivity = Oceananigans.CenterField(grid)
    size(vertical_diffusivity) == size(horizontal_diffusivity) ||
        throw(DimensionMismatch(
            "localized-estuary vertical and horizontal fields must have " *
            "the same size",
        ))
    fill!(vertical_diffusivity, zero(FT))
    fill!(horizontal_diffusivity, zero(FT))

    dynamically_gated = config.ocean_river_mouth_mixing in (
        :localized_estuary_v4,
        :localized_estuary_v5,
    )
    eligible_vertical_diffusivity = dynamically_gated ?
        Oceananigans.CenterField(grid) : nothing
    eligible_horizontal_diffusivity = dynamically_gated ?
        Oceananigans.CenterField(grid) : nothing
    active_receiver_mask = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    isnothing(eligible_vertical_diffusivity) || fill!(
        eligible_vertical_diffusivity,
        zero(FT),
    )
    isnothing(eligible_horizontal_diffusivity) || fill!(
        eligible_horizontal_diffusivity,
        zero(FT),
    )
    fill!(active_receiver_mask, zero(FT))

    tracer_selector = config.ocean_river_mouth_mixing in (
        :localized_estuary_v3,
        :localized_estuary_v4,
        :localized_estuary_v5,
    ) ? LocalizedEstuaryTracerIndex{2}() :
        AllLocalizedEstuaryTracers()
    closure = LocalizedEstuaryDiffusivity(
        vertical_diffusivity,
        horizontal_diffusivity,
        tracer_selector,
    )
    return (
        closures = (closure,),
        fields = RiverMouthMixingFields(
            vertical_diffusivity,
            horizontal_diffusivity,
            eligible_vertical_diffusivity,
            eligible_horizontal_diffusivity,
            active_receiver_mask,
            dynamically_gated,
            FT(config.ocean_river_mouth_reference_freshwater_mass_flux_kgm2s),
        ),
    )
end

@inline function _river_mouth_neighbor_i(i, offset, Nx)
    return mod1(i + offset, Nx)
end

function _initialize_river_mouth_mixing!(
    mixing::RiverMouthMixingFields,
    land::TerrariumRunoffLand,
    grid,
    config::ExperimentConfig,
)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    FT = eltype(grid)
    κz = zeros(FT, Nx, Ny, Nz)
    κh = zeros(FT, Nx, Ny, Nz)
    active_receiver_mask = zeros(FT, Nx, Ny, 1)
    target_i = Int.(Array(land.target_i))
    target_j = Int.(Array(land.target_j))
    length(target_i) == length(target_j) || throw(DimensionMismatch(
        "river-mouth target index vectors have different lengths",
    ))
    isempty(target_i) && error("river-mouth mixing found no runoff receivers")

    center = Oceananigans.Center()
    mixing_depth = FT(config.ocean_river_mouth_mixing_depth_m)
    vertical_value = FT(config.ocean_river_mouth_vertical_diffusivity_m2s)
    horizontal_value = FT(config.ocean_river_mouth_horizontal_diffusivity_m2s)
    cardinal_offsets = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
    horizontal_offsets =
        config.ocean_river_mouth_mixing == :localized_estuary_v1 ?
        cardinal_offsets : ((0, 0),)

    # Vertical mixing acts in the actual routed columns for both versioned
    # schemes. V1 writes the cell-centred horizontal coefficient over the
    # receiver plus one cardinal ring. Face interpolation then gives that
    # field an additional, unintended outer ring of direct support. V2 writes
    # the coefficient only at receivers; conservative face interpolation still
    # exchanges with each adjacent resolved plume cell, without the second ring.
    for (i, j) in zip(target_i, target_j)
        1 <= i <= Nx && 1 <= j <= Ny || error(
            "river-mouth target ($i, $j) lies outside the ocean grid",
        )
        mixing.dynamically_gated || (active_receiver_mask[i, j, 1] = one(FT))
        for k in 1:Nz
            Oceananigans.Grids.inactive_node(
                i,
                j,
                k,
                cpu_grid,
                center,
                center,
                center,
            ) && continue
            z = Oceananigans.Grids.znode(
                i,
                j,
                k,
                cpu_grid,
                center,
                center,
                center,
            )
            z >= -mixing_depth && (κz[i, j, k] = vertical_value)
        end

        for (di, dj) in horizontal_offsets
            ii = _river_mouth_neighbor_i(i, di, Nx)
            jj = j + dj
            1 <= jj <= Ny || continue
            for k in 1:Nz
                Oceananigans.Grids.inactive_node(
                    ii,
                    jj,
                    k,
                    cpu_grid,
                    center,
                    center,
                    center,
                ) && continue
                z = Oceananigans.Grids.znode(
                    ii,
                    jj,
                    k,
                    cpu_grid,
                    center,
                    center,
                    center,
                )
                z >= -mixing_depth && (κh[ii, jj, k] = horizontal_value)
            end
        end
    end

    if mixing.dynamically_gated
        Oceananigans.set!(mixing.eligible_vertical_diffusivity, κz)
        Oceananigans.set!(mixing.eligible_horizontal_diffusivity, κh)
        Oceananigans.fill_halo_regions!(mixing.eligible_vertical_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.eligible_horizontal_diffusivity)
        fill!(mixing.vertical_diffusivity, zero(FT))
        fill!(mixing.horizontal_diffusivity, zero(FT))
        fill!(mixing.active_receiver_mask, zero(FT))
        # Establish valid zero halos before the first coupled interface update.
        # The live coefficient and activity fields are derived from runoff and
        # are therefore intentionally reconstructed, rather than checkpointed.
        Oceananigans.fill_halo_regions!(mixing.vertical_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.horizontal_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.active_receiver_mask)
    else
        Oceananigans.set!(mixing.vertical_diffusivity, κz)
        Oceananigans.set!(mixing.horizontal_diffusivity, κh)
        Oceananigans.fill_halo_regions!(mixing.vertical_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.horizontal_diffusivity)
        Oceananigans.set!(mixing.active_receiver_mask, active_receiver_mask)
        Oceananigans.fill_halo_regions!(mixing.active_receiver_mask)
    end
    return (
        runoff_receiver_cells = length(target_i),
        vertically_eligible_cells = count(!iszero, κz),
        horizontally_eligible_cells = count(!iszero, κh),
        vertically_mixed_cells = mixing.dynamically_gated ? 0 : count(!iszero, κz),
        horizontally_mixed_cells = mixing.dynamically_gated ? 0 : count(!iszero, κh),
        vertical_diffusivity_m2s = Float64(vertical_value),
        horizontal_diffusivity_m2s = Float64(horizontal_value),
        mixing_depth_m = Float64(mixing_depth),
        dynamically_gated = mixing.dynamically_gated,
        reference_freshwater_mass_flux_kgm2s = Float64(
            mixing.reference_freshwater_mass_flux_kgm2s,
        ),
    )
end

_initialize_river_mouth_mixing!(::Nothing, land, grid, config) = nothing

function _river_mouth_surface_footprint(
    land::TerrariumRunoffLand,
    grid,
    mixing_mode::Symbol = :localized_estuary_v1,
)
    (mixing_mode == :none || _is_localized_estuary_mixing(mixing_mode)) ||
        error("unsupported river-mouth footprint mode $mixing_mode")
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    receiver_mask = zeros(Float64, Nx, Ny)
    horizontal_mask = zeros(Float64, Nx, Ny)
    target_i = Int.(Array(land.target_i))
    target_j = Int.(Array(land.target_j))
    length(target_i) == length(target_j) || throw(DimensionMismatch(
        "river-mouth target index vectors have different lengths",
    ))
    center = Oceananigans.Center()
    cardinal_offsets = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
    horizontal_offsets = mixing_mode == :localized_estuary_v1 ?
        cardinal_offsets :
        mixing_mode in (
            :localized_estuary_v2,
            :localized_estuary_v3,
            :localized_estuary_v4,
            :localized_estuary_v5,
        ) ?
        ((0, 0),) : ()

    for (i, j) in zip(target_i, target_j)
        Oceananigans.Grids.inactive_node(
            i,
            j,
            Nz,
            cpu_grid,
            center,
            center,
            center,
        ) || (receiver_mask[i, j] = 1)
        for (di, dj) in horizontal_offsets
            ii = _river_mouth_neighbor_i(i, di, Nx)
            jj = j + dj
            1 <= jj <= Ny || continue
            Oceananigans.Grids.inactive_node(
                ii,
                jj,
                Nz,
                cpu_grid,
                center,
                center,
                center,
            ) || (horizontal_mask[ii, jj] = 1)
        end
    end
    return (; receiver_mask, horizontal_mask)
end

function _river_mouth_mixing_footprint_provenance(mode::Symbol)
    if mode == :none
        return "none"
    elseif mode == :localized_estuary_v1
        return "vertical_at_routed_receivers_horizontal_over_one_cardinal_wet_cell_ring"
    elseif mode == :localized_estuary_v2
        return "vertical_and_horizontal_coefficients_at_routed_receivers_face_exchange_with_cardinal_neighbors"
    elseif mode == :localized_estuary_v3
        return "salinity_only_vertical_and_horizontal_coefficients_at_routed_receivers_face_exchange_with_cardinal_neighbors"
    elseif mode == :localized_estuary_v4
        return "salinity_only_coefficients_gated_each_coupled_step_by_positive_routed_freshwater_flux_at_receiver"
    elseif mode == :localized_estuary_v5
        return "salinity_only_coefficients_smoothly_gated_each_coupled_step_by_routed_freshwater_mass_flux_at_receiver"
    else
        error("unsupported river-mouth mixing mode $mode")
    end
end

function _ocean_mixing_closure_provenance(mode::Symbol)
    if mode == :none
        return "CATKE_GentMcWilliams_isopycnal_biharmonic"
    elseif mode == :localized_estuary_v1
        return "CATKE_GentMcWilliams_isopycnal_biharmonic_localized_estuary"
    elseif mode == :localized_estuary_v2
        return "CATKE_GentMcWilliams_isopycnal_biharmonic_localized_estuary_v2_receiver_centered"
    elseif mode == :localized_estuary_v3
        return "CATKE_GentMcWilliams_isopycnal_biharmonic_localized_estuary_v3_receiver_centered_salinity_only"
    elseif mode == :localized_estuary_v4
        return "CATKE_GentMcWilliams_isopycnal_biharmonic_localized_estuary_v4_active_runoff_salinity_only"
    elseif mode == :localized_estuary_v5
        return "CATKE_GentMcWilliams_isopycnal_biharmonic_localized_estuary_v5_loading_scaled_runoff_salinity_only"
    else
        error("unsupported river-mouth mixing mode $mode")
    end
end

function _ocean_surface_cell_area(grid)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    return [
        Float64(Oceananigans.Operators.Azᶜᶜᶜ(i, j, Nz, cpu_grid))
        for i in 1:Nx, j in 1:Ny
    ]
end

@inline _ocean_build_result(ocean, mixing, ::Val{false}) = ocean
@inline _ocean_build_result(ocean, mixing, ::Val{true}) = (
    ocean = ocean,
    river_mouth_mixing = mixing,
)

"""
Apply a documented coarse-grid marginal-sea correction to an already regridded
latitude--longitude bottom-height array.

All modes change only cells that are already wet and never open positive land
cells. `persian_gulf_1degree_v1` is the original basin-specific diagnostic
candidate. `shallow_culdesacs_1degree_v1` finds every exactly minimum-depth
(10 m) wet terminal with at most one wet cardinal neighbour, deepens the
terminal to 30 m, and also deepens its sole connector to 30 m when that
connector is shallower. `shallow_culdesacs_hormuz_1degree_v2` composes that
global terminal correction with the documented 30 m Persian-Gulf / 50 m
Hormuz-corridor minimum depths. The composition retains the global fix for the
independent Timor/Arafura temperature hotspot while repairing the two shallow
Hormuz bottlenecks left by v1. Coastline and wet-cell connectivity are
preserved. `shallow_culdesacs_polar_boundary_1degree_v3` instead composes v1
with a 30 m minimum depth for already-wet cells in the artificial 70--75 degree
northern boundary band. This targets the independently preserved day-440
fresh and saline extrema without retaining the unrelated Hormuz experiment.
It opens no land, restores no tracer and preserves wet-cell connectivity.
`shallow_culdesacs_global_40m_1degree_v4` composes the v3 targets with a 40 m
minimum depth for every already-wet one-degree cell. The exact v3 day-730 state
showed that eleven of twelve columns at or above 39 degrees Celsius were
shallower than 30 m, while the remaining salinity tails occupied shallow cells
in both the Persian Gulf and artificial Arctic boundary. NumericalEarth deposits
shortwave radiation transmitted to the seabed in the lowest wet cell, so the
global minimum is a resolution-aware treatment of the common mechanism rather
than separate regional tuning. It changes neither land nor wet connectivity.
"""
function _wet_cardinal_neighbours(bottom_height, i, j)
    Nx, Ny = size(bottom_height)
    neighbours = Tuple{Int, Int}[
        (mod1(i - 1, Nx), j),
        (mod1(i + 1, Nx), j),
    ]
    j > 1 && push!(neighbours, (i, j - 1))
    j < Ny && push!(neighbours, (i, j + 1))
    filter!(index -> bottom_height[index...] < 0, neighbours)
    return neighbours
end

function _apply_latitude_longitude_bathymetry_correction!(
    bottom_height::AbstractMatrix,
    longitude,
    latitude,
    correction::Symbol,
)
    correction == :none && return (
        terminal_cells = 0,
        polar_boundary_cells = 0,
        global_minimum_depth_cells = 0,
        changed_cells = 0,
        total_deepening_m = 0.0,
        maximum_deepening_m = 0.0,
    )
    correction in (
        :persian_gulf_1degree_v1,
        :shallow_culdesacs_1degree_v1,
        :shallow_culdesacs_hormuz_1degree_v2,
        :shallow_culdesacs_polar_boundary_1degree_v3,
        :shallow_culdesacs_global_40m_1degree_v4,
    ) || throw(
        ArgumentError("unsupported latitude--longitude bathymetry correction: $correction"),
    )
    size(bottom_height) == (length(longitude), length(latitude)) || throw(
        DimensionMismatch(
            "bottom-height dimensions must match longitude and latitude coordinates",
        ),
    )

    terminal_cells = 0
    polar_boundary_cells = 0
    global_minimum_depth_cells = 0
    target_depths = Dict{Tuple{Int, Int}, Float64}()
    if correction in (
        :shallow_culdesacs_1degree_v1,
        :shallow_culdesacs_hormuz_1degree_v2,
        :shallow_culdesacs_polar_boundary_1degree_v3,
        :shallow_culdesacs_global_40m_1degree_v4,
    )
        for j in eachindex(latitude), i in eachindex(longitude)
            original = bottom_height[i, j]
            -10.01 <= original < 0 || continue
            neighbours = _wet_cardinal_neighbours(bottom_height, i, j)
            length(neighbours) <= 1 || continue
            terminal_cells += 1
            target_depths[(i, j)] = 30.0
            isempty(neighbours) && continue
            connector = only(neighbours)
            connector_height = bottom_height[connector...]
            -30 < connector_height < 0 &&
                (target_depths[connector] = 30.0)
        end
    end
    if correction in (
        :shallow_culdesacs_polar_boundary_1degree_v3,
        :shallow_culdesacs_global_40m_1degree_v4,
    )
        for j in eachindex(latitude), i in eachindex(longitude)
            original = bottom_height[i, j]
            Float64(latitude[j]) >= 70 || continue
            -30 < original < 0 || continue
            polar_boundary_cells += 1
            index = (i, j)
            target_depths[index] = max(
                get(target_depths, index, 0.0),
                30.0,
            )
        end
    end
    if correction == :shallow_culdesacs_global_40m_1degree_v4
        for j in eachindex(latitude), i in eachindex(longitude)
            original = bottom_height[i, j]
            -40 < original < 0 || continue
            global_minimum_depth_cells += 1
            index = (i, j)
            target_depths[index] = max(
                get(target_depths, index, 0.0),
                40.0,
            )
        end
    end
    if correction in (
        :persian_gulf_1degree_v1,
        :shallow_culdesacs_hormuz_1degree_v2,
    )
        for j in eachindex(latitude), i in eachindex(longitude)
            original = bottom_height[i, j]
            original < 0 || continue
            λ = mod(Float64(longitude[i]), 360)
            φ = Float64(latitude[j])

            # Existing Gulf water cells: retain real shelf geometry where
            # deeper, but prevent a one-degree evaporative cell from
            # containing only 10 m of water. The nested corridor captures the
            # resolved Hormuz exit.
            target_depth_m = if 55 <= λ <= 59 && 24 <= φ <= 27
                50.0
            elseif 47 <= λ <= 56 && 24 <= φ <= 30
                30.0
            else
                0.0
            end
            if target_depth_m != 0
                index = (i, j)
                target_depths[index] = max(
                    get(target_depths, index, 0.0),
                    target_depth_m,
                )
            end
        end
    end

    changed_cells = 0
    total_deepening_m = 0.0
    maximum_deepening_m = 0.0
    for ((i, j), target_depth_m) in sort!(collect(target_depths); by = first)
        original = bottom_height[i, j]
        corrected = min(original, -convert(eltype(bottom_height), target_depth_m))
        corrected == original && continue
        deepening_m = Float64(original - corrected)
        bottom_height[i, j] = corrected
        changed_cells += 1
        total_deepening_m += deepening_m
        maximum_deepening_m = max(maximum_deepening_m, deepening_m)
    end
    return (;
        terminal_cells,
        polar_boundary_cells,
        global_minimum_depth_cells,
        changed_cells,
        total_deepening_m,
        maximum_deepening_m,
    )
end

function _apply_latitude_longitude_bathymetry_correction!(
    bottom_height,
    grid,
    correction::Symbol,
)
    correction == :none && return nothing
    longitude = Float64.(collect(Oceananigans.Grids.λnodes(
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )))
    latitude = Float64.(collect(Oceananigans.Grids.φnodes(
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )))
    host_bottom = Array(Oceananigans.interior(bottom_height, :, :, 1))
    statistics = _apply_latitude_longitude_bathymetry_correction!(
        host_bottom,
        longitude,
        latitude,
        correction,
    )
    statistics.changed_cells > 0 || error(
        "configured bathymetry correction $correction changed no wet cells",
    )
    Oceananigans.set!(bottom_height, host_bottom)
    Oceananigans.fill_halo_regions!(bottom_height)
    Oceananigans.Architectures.synchronize(Oceananigans.architecture(grid))
    @info "Applied ocean bathymetry correction" correction statistics
    return nothing
end

function _latitude_longitude_polar_velocity_sponge(
    config::ExperimentConfig,
    ::Type{FT},
) where {FT}
    # A bounded latitude--longitude ocean necessarily has artificial walls at
    # ±75 degrees. Apply a smooth, configurable momentum-only absorbing layer
    # so wind/ice stress cannot accumulate barotropic energy at the wall. The
    # same forcing is applied to dynamic sea ice below: the rejected ocean-only
    # candidate allowed undamped ice to accelerate first and transmit its wall
    # momentum back into the ocean. A narrow five-degree layer subsequently
    # controlled velocity but produced a phase-sensitive free-surface
    # depression. Production therefore uses a wider ramp configured in YAML.
    # Tracers are not restored: salt and heat remain prognostic and are
    # redistributed by the physical one-degree closure above.
    return Oceananigans.Relaxation(
        rate = inv(FT(config.ocean_polar_sponge_timescale_hours * 3_600)),
        mask = _latitude_longitude_polar_sponge_mask(
            FT,
            config.ocean_polar_sponge_start_latitude_degrees,
            config.ocean_polar_sponge_stop_latitude_degrees,
        ),
    )
end

"""ClimaSeaIce-compatible u/v relaxation at the artificial polar walls."""
struct _SeaIcePolarMomentumRelaxation{C, R, M}
    rate::R
    mask::M
end

@inline _sea_ice_velocity_latitude(
    ::Val{:u},
    j,
    grid,
) = Oceananigans.Grids.φnode(j, grid, Oceananigans.Center())

@inline _sea_ice_velocity_latitude(
    ::Val{:v},
    j,
    grid,
) = Oceananigans.Grids.φnode(j, grid, Oceananigans.Face())

@inline function (forcing::_SeaIcePolarMomentumRelaxation{C})(
    i,
    j,
    k,
    grid,
    fields,
) where {C}
    φ = _sea_ice_velocity_latitude(Val(C), j, grid)
    velocity = @inbounds getproperty(fields, C)[i, j, 1]
    return -forcing.rate * forcing.mask(zero(φ), φ, zero(φ)) * velocity
end

function _latitude_longitude_sea_ice_velocity_sponge(
    config::ExperimentConfig,
    ::Type{FT},
    ::Val{C},
) where {FT, C}
    rate = inv(FT(config.ocean_polar_sponge_timescale_hours * 3_600))
    mask = _latitude_longitude_polar_sponge_mask(
        FT,
        config.ocean_polar_sponge_start_latitude_degrees,
        config.ocean_polar_sponge_stop_latitude_degrees,
    )
    return _SeaIcePolarMomentumRelaxation{C, typeof(rate), typeof(mask)}(
        rate,
        mask,
    )
end

"""
Construct NumericalEarth's standard sea-ice simulation while exposing the
`SeaIceModel` momentum-forcing keyword that the current convenience constructor
does not forward. This keeps the upstream thermodynamic and stress setup
unchanged and is used only to damp u/v next to the artificial ±75-degree walls.
"""
function _sea_ice_simulation_with_momentum_forcing(
    grid,
    ocean;
    dynamics,
    advection,
    forcing,
    conductivity,
    internal_heat_flux,
    ice_salinity,
    ice_heat_capacity,
    ice_consolidation_thickness,
    sea_ice_density,
    snow_density,
    boundary_conditions = NamedTuple(),
)
    FT = eltype(grid)
    top_surface_temperature = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    top_heat_boundary_condition =
        ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions.PrescribedTemperature(
            top_surface_temperature.data,
        )
    surface_ocean_salinity = isnothing(ocean) ? zero(FT) :
        NumericalEarth.EarthSystemModels.ocean_surface_salinity(ocean)
    bottom_heat_boundary_condition =
        ClimaSeaIce.SeaIceThermodynamics.IceWaterThermalEquilibrium(
            surface_ocean_salinity,
        )
    ice_thermodynamics = ClimaSeaIce.sea_ice_slab_thermodynamics(
        grid;
        internal_heat_flux,
        top_heat_boundary_condition,
        bottom_heat_boundary_condition,
    )
    snow_thermodynamics =
        NumericalEarth.SeaIces.default_snow_thermodynamics(grid)
    phase_transitions = ClimaSeaIce.PhaseTransitions(
        FT;
        heat_capacity = FT(ice_heat_capacity),
        density = FT(sea_ice_density),
    )
    top_heat_flux = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    bottom_heat_flux = similar(top_heat_flux)
    snowfall = similar(top_heat_flux)

    model = ClimaSeaIce.SeaIceModel(
        grid;
        clock = Oceananigans.Clock(grid),
        ice_salinity,
        advection,
        ice_consolidation_thickness,
        sea_ice_density,
        snow_density,
        phase_transitions,
        ice_thermodynamics,
        snow_thermodynamics,
        snowfall,
        dynamics,
        bottom_heat_flux,
        top_heat_flux,
        forcing,
        boundary_conditions,
    )
    return Oceananigans.Simulation(
        model;
        Δt = 5Oceananigans.Units.minutes,
        stop_time = Inf,
        verbose = false,
    )
end

@inline function _sea_ice_immersed_u_drag(
    i,
    j,
    k,
    grid,
    clock,
    fields,
    coefficient,
)
    return @inbounds -coefficient * fields.u[i, j, k]
end

@inline function _sea_ice_immersed_v_drag(
    i,
    j,
    k,
    grid,
    clock,
    fields,
    coefficient,
)
    return @inbounds -coefficient * fields.v[i, j, k]
end

"""
Construct the lateral immersed-boundary stress used by ClimaSeaIce's own
coastline example while retaining NumericalEarth's vector-aware tripolar fold
conditions. NumericalEarth 0.6.0's convenience constructor does not forward a
`boundary_conditions` keyword, so without this explicit path all eight
land-facing stress fluxes fall through ClimaSeaIce's documented zero/TODO
fallback even though peripheral velocities themselves are set to zero.
"""
function _sea_ice_immersed_velocity_boundary_conditions(grid, coefficient)
    coefficient > 0 || throw(ArgumentError(
        "immersed sea-ice drag coefficient must be positive",
    ))
    u_flux = Oceananigans.FluxBoundaryCondition(
        _sea_ice_immersed_u_drag;
        discrete_form = true,
        parameters = coefficient,
    )
    v_flux = Oceananigans.FluxBoundaryCondition(
        _sea_ice_immersed_v_drag;
        discrete_form = true,
        parameters = coefficient,
    )
    u_immersed = Oceananigans.ImmersedBoundaryCondition(
        top = nothing,
        bottom = nothing,
        west = nothing,
        east = nothing,
        south = u_flux,
        north = u_flux,
    )
    v_immersed = Oceananigans.ImmersedBoundaryCondition(
        top = nothing,
        bottom = nothing,
        south = nothing,
        north = nothing,
        west = v_flux,
        east = v_flux,
    )
    interface = NumericalEarth.EarthSystemModels.InterfaceComputations
    u_base = interface.vector_component_boundary_conditions(
        grid,
        (Oceananigans.Face(), Oceananigans.Center(), nothing),
    )
    v_base = interface.vector_component_boundary_conditions(
        grid,
        (Oceananigans.Center(), Oceananigans.Face(), nothing),
    )
    u = Oceananigans.FieldBoundaryConditions(
        u_base.west,
        u_base.east,
        u_base.south,
        u_base.north,
        u_base.bottom,
        u_base.top,
        u_immersed,
    )
    v = Oceananigans.FieldBoundaryConditions(
        v_base.west,
        v_base.east,
        v_base.south,
        v_base.north,
        v_base.bottom,
        v_base.top,
        v_immersed,
    )
    return (; u, v)
end

"""
Choose the split-explicit substep request for a tripolar grid.

The production 360x180 grid always receives the documented 70 substeps.
Oceananigans' default averaging kernel retains 50 of those and extends the
northern-fold halo to 52 cells.  Compact construction/CUDA gates need fewer
substeps so that this required halo remains smaller than their interior.
NumericalEarth currently allocates its freshwater-flux field on the extended
grid even when `extend_halos=false`, so reducing only compact-gate subcycling
is the reliable package-compatible route.
"""
_tripolar_barotropic_substeps(nlatitude::Integer) = min(70, nlatitude - 3)

"""
    _remove_minor_basins_preserving_alignment!(bottom_height, keep_major_basins)

Remove disconnected ocean basins without moving the bathymetry relative to its
tripolar coordinates. NumericalEarth 0.6.0's periodic extension starts at
`Nx ÷ 2:Nx`, which prepends `Nx ÷ 2 + 1` columns while assigning an offset for
only `Nx ÷ 2`. Copying indices `1:Nx` back therefore shifts all rows below the
north fold one cell east. Build an exact half-domain extension from the field's
interior and call the same connected-component implementation on correctly
aligned offset axes.
"""
function _remove_minor_basins_preserving_alignment!(
    bottom_height,
    keep_major_basins::Integer,
)
    keep_major_basins ≥ 1 || throw(ArgumentError(
        "keep_major_basins must be at least one",
    ))
    grid = bottom_height.grid
    Nx, Ny, = size(grid)
    iseven(Nx) || throw(ArgumentError(
        "periodic basin cleanup requires an even longitudinal grid size",
    ))
    Oceananigans.topology(grid, 1) === Oceananigans.Periodic || throw(
        ArgumentError("alignment-preserving basin cleanup requires periodic longitude"),
    )

    cpu_bottom_height = Oceananigans.on_architecture(
        Oceananigans.CPU(),
        bottom_height,
    )
    source = Array(Oceananigans.interior(cpu_bottom_height, :, :, 1))
    half = Nx ÷ 2
    extended_parent = vcat(
        source[half+1:Nx, :],
        source,
        source[1:half, :],
    )
    extended = NumericalEarth.Bathymetry.OffsetArray(
        extended_parent,
        -half,
        0,
    )
    NumericalEarth.Bathymetry.remove_minor_basins!(
        extended,
        keep_major_basins,
        (Nx, Ny, 1),
    )
    corrected = Array(extended[1:Nx, 1:Ny])

    removed = (source .< 0) .& .!(corrected .< 0)
    unexpected_change = (corrected .!= source) .& .!removed
    any(unexpected_change) && error(
        "minor-basin cleanup changed cells other than removed wet cells",
    )

    Oceananigans.set!(bottom_height, corrected)
    Oceananigans.fill_halo_regions!(bottom_height)
    return (
        removed_cells = count(removed),
        retained_wet_cells = count(<(0), corrected),
    )
end

function _build_ocean(
    config::ExperimentConfig;
    return_auxiliary::Bool = false,
)
    arch = _ocean_architecture(config)
    FT = Float32
    if config.ocean_grid == :tripolar_1degree
        config.ocean_dynamics || throw(
            ArgumentError("the one-degree tripolar ocean requires ocean_dynamics=true"),
        )
        z = Oceananigans.ExponentialDiscretization(
            config.ocean_nlayers,
            -FT(config.ocean_depth_m),
            zero(FT);
            mutable = true,
        )
        grid = Oceananigans.TripolarGrid(
            arch,
            FT;
            size = (
                config.ocean_nlongitude,
                config.ocean_nlatitude,
                config.ocean_nlayers,
            ),
            z,
            halo = (5, 5, 4),
        )
        bottom_height = NumericalEarth.Bathymetry.regrid_bathymetry(
            grid;
            minimum_depth = FT(10),
            interpolation_passes = 10,
            # NumericalEarth 0.6.0's periodic minor-basin extension shifts the
            # complete field one longitude cell before copying it back. Keep
            # the aligned pre-basin result in its cache and apply the corrected
            # cleanup below.
            major_basins = Inf,
        )
        basin_cleanup = _remove_minor_basins_preserving_alignment!(
            bottom_height,
            2,
        )
        @info "Applied alignment-preserving tripolar basin cleanup" basin_cleanup
        if config.ocean_tripolar_wet_mask == :etopo_majority_area_polar
            majority_mask = _compute_tripolar_majority_wet_mask(grid)
            adjustment = _apply_tripolar_majority_wet_mask!(
                bottom_height,
                majority_mask,
            )
            final_cleanup = _remove_minor_basins_preserving_alignment!(
                bottom_height,
                2,
            )
            final_bottom = Array(Oceananigans.interior(
                Oceananigans.on_architecture(
                    Oceananigans.CPU(),
                    bottom_height,
                ),
                :,
                :,
                1,
            ))
            provenance_path = _save_tripolar_majority_wet_provenance(
                config.output_dir,
                majority_mask,
                adjustment,
                final_bottom,
                final_cleanup,
            )
            @info "Applied ETOPO majority-area tripolar wet mask" changed_cells=adjustment.changed_cells wet_to_land_cells=adjustment.wet_to_land_cells land_to_wet_cells=adjustment.land_to_wet_cells reliable_cells=adjustment.reliable_cells undersampled_canonical_cells=adjustment.undersampled_canonical_cells final_cleanup provenance_path
        end
        grid = Oceananigans.ImmersedBoundaryGrid(
            grid,
            Oceananigans.GridFittedBottom(bottom_height);
            active_cells_map = true,
        )
        # Use the same production schemes that completed the exact 730-day
        # latitude--longitude integration. The earlier tripolar development
        # route retained fully explicit WENO5 vector-invariant momentum and
        # tracers and showed trajectory-sensitive failures on 1--100 day
        # horizons. Centered vector-invariant momentum plus horizontal WENO7
        # with adaptive implicit vertical transport removes that unnecessary
        # branch difference while preserving the tripolar fold and free surface.
        momentum_advection, tracer_advection = _production_ocean_advection(FT)
        free_surface = Oceananigans.SplitExplicitFreeSurface(
            grid;
            substeps = _tripolar_barotropic_substeps(config.ocean_nlatitude),
        )
        base_closure = ClimaOcean.OceanConfigurations.default_one_degree_closure(
            κ_skew = FT(500),
            κ_symmetric = FT(200),
            background_ν = FT(1e-5),
        )
        river_mouth_mixing = _river_mouth_mixing_closure(config, grid)
        closure = isempty(river_mouth_mixing.closures) ?
            base_closure :
            (base_closure..., river_mouth_mixing.closures...)
        if _is_localized_estuary_mixing(config.ocean_river_mouth_mixing)
            length(base_closure) == 4 || error(
                "localized-estuary production path expected four base " *
                "Oceananigans closures, found $(length(base_closure))",
            )
            length(closure) == 5 || error(
                "localized-estuary production path must retain the native " *
                "five-closure GPU tuple, found $(length(closure))",
            )
        end
        radiative_forcing = _ocean_radiative_forcing(
            config,
            grid,
            bottom_height,
        )
        ocean = ClimaOcean.ocean_simulation(
            grid;
            momentum_advection,
            tracer_advection,
            free_surface,
            closure,
            radiative_forcing,
            bottom_drag_coefficient = FT(0.003),
        )
        if config.ocean_river_mouth_mixing in (
            :localized_estuary_v3,
            :localized_estuary_v4,
            :localized_estuary_v5,
        )
            tracer_names = propertynames(ocean.model.tracers)
            length(tracer_names) >= 2 && tracer_names[2] == :S || error(
                "salinity-only localized-estuary mixing requires salinity to " *
                "be ocean tracer index 2; " *
                "found $(tracer_names)",
            )
        end
        return _ocean_build_result(
            ocean,
            river_mouth_mixing.fields,
            Val(return_auxiliary),
        )
    end
    if config.ocean_grid == :latitude_longitude_1degree
        config.ocean_dynamics || throw(
            ArgumentError("the one-degree production ocean requires ocean_dynamics=true"),
        )
        momentum_advection, tracer_advection = _production_ocean_advection(FT)
        # ClimaOcean's convenience constructor currently inherits Oceananigans'
        # Float64 global default and does not expose an FT keyword. Construct its
        # documented 360x150 realistic-bathymetry grid explicitly so atmosphere,
        # regridding weights, ocean, and ice all use Float32 on CUDA. Mixed
        # Float64/Float32 sparse matvec falls back from cuSPARSE to scalar indexing.
        z = Oceananigans.ExponentialDiscretization(
            config.ocean_nlayers,
            -FT(config.ocean_depth_m),
            zero(FT);
            mutable = true,
        )
        grid = Oceananigans.LatitudeLongitudeGrid(
            arch,
            FT;
            size = (360, 150, config.ocean_nlayers),
            z,
            halo = (7, 7, 7),
            latitude = (FT(-75), FT(75)),
            longitude = (zero(FT), FT(360)),
        )
        bottom_height = NumericalEarth.Bathymetry.regrid_bathymetry(
            grid;
            minimum_depth = FT(10),
            interpolation_passes = 5,
            major_basins = 3,
        )
        _apply_latitude_longitude_bathymetry_correction!(
            bottom_height,
            grid,
            config.ocean_bathymetry_correction,
        )
        grid = Oceananigans.ImmersedBoundaryGrid(
            grid,
            Oceananigans.GridFittedBottom(bottom_height);
            active_cells_map = true,
        )
        closure = _production_latitude_longitude_closure(FT)
        velocity_sponge = _latitude_longitude_polar_velocity_sponge(config, FT)
        radiative_forcing = _ocean_radiative_forcing(
            config,
            grid,
            bottom_height,
        )
        ocean = ClimaOcean.ocean_simulation(
            grid;
            momentum_advection,
            tracer_advection,
            closure,
            radiative_forcing,
            forcing = (u = velocity_sponge, v = velocity_sponge),
            coriolis = Oceananigans.HydrostaticSphericalCoriolis(
                FT;
                scheme = Oceananigans.EnstrophyConserving(FT),
            ),
            bottom_drag_coefficient = FT(0.003),
        )
        return _ocean_build_result(ocean, nothing, Val(return_auxiliary))
    end

    momentum_advection, tracer_advection = _idealized_ocean_advection(FT)
    halo = config.ocean_dynamics ? (5, 5, 5) : (3, 3, 3)
    grid = Oceananigans.LatitudeLongitudeGrid(
        arch,
        FT;
        size = (
            config.ocean_nlongitude,
            config.ocean_nlatitude,
            config.ocean_nlayers,
        ),
        halo,
        longitude = (0, 360),
        latitude = (-80, 80),
        # Exercise the same z-star freshwater-volume pathway as the production
        # ocean while retaining uniformly spaced levels in the compact gate.
        z = Oceananigans.MutableVerticalDiscretization(collect(range(
            -FT(config.ocean_depth_m),
            zero(FT);
            length = config.ocean_nlayers + 1,
        ))),
    )
    if config.ocean_dynamics
        ocean = ClimaOcean.ocean_simulation(
            grid;
            momentum_advection,
            tracer_advection,
            closure = ClimaOcean.simplified_ocean_closure(FT),
            coriolis = Oceananigans.HydrostaticSphericalCoriolis(
                FT;
                scheme = Oceananigans.EnstrophyConserving(FT),
            ),
            bottom_drag_coefficient = FT(0.003),
        )
        return _ocean_build_result(ocean, nothing, Val(return_auxiliary))
    end
    ocean = ClimaOcean.ocean_simulation(
        grid;
        momentum_advection = nothing,
        tracer_advection = nothing,
        closure = nothing,
        coriolis = nothing,
    )
    return _ocean_build_result(ocean, nothing, Val(return_auxiliary))
end

_sea_ice_pressure_formulation(::Val{:replacement_pressure}) =
    ClimaSeaIce.Rheologies.ReplacementPressure()
_sea_ice_pressure_formulation(::Val{:ice_strength}) =
    ClimaSeaIce.Rheologies.IceStrength()
_sea_ice_pressure_formulation(formulation::Symbol) =
    _sea_ice_pressure_formulation(Val(formulation))

function _build_sea_ice(config::ExperimentConfig, ocean)
    grid = ocean.model.grid
    FT = eltype(grid)
    common_kwargs = (
        conductivity = FT(2),
        internal_heat_flux = ClimaSeaIce.SeaIceThermodynamics.ConductiveFlux(
            FT;
            conductivity = FT(2),
        ),
        ice_salinity = FT(4),
        ice_heat_capacity = FT(2100),
        ice_consolidation_thickness = FT(0.05),
        sea_ice_density = FT(900),
        snow_density = FT(330),
    )
    if config.sea_ice_dynamics
        pressure_formulation = _sea_ice_pressure_formulation(
            config.sea_ice_pressure_formulation,
        )
        solver = ClimaSeaIce.SeaIceDynamics.SplitExplicitSolver(
            grid;
            substeps = config.sea_ice_momentum_substeps,
        )
        dynamics = NumericalEarth.SeaIces.sea_ice_dynamics(
            grid,
            ocean;
            coriolis = ocean.model.coriolis,
            rheology = ClimaSeaIce.Rheologies.ElastoViscoPlasticRheology(
                FT;
                pressure_formulation,
            ),
            solver,
        )
        advection = ConservativeSeaIceAdvection(
            Oceananigans.WENO(FT; order = 5),
        )
        coastal_boundary_conditions =
            config.sea_ice_immersed_boundary_drag_coefficient > 0 ?
            _sea_ice_immersed_velocity_boundary_conditions(
                grid,
                FT(config.sea_ice_immersed_boundary_drag_coefficient),
            ) : nothing
        if config.ocean_grid == :latitude_longitude_1degree
            u_sponge = _latitude_longitude_sea_ice_velocity_sponge(
                config,
                FT,
                Val(:u),
            )
            v_sponge = _latitude_longitude_sea_ice_velocity_sponge(
                config,
                FT,
                Val(:v),
            )
            return _sea_ice_simulation_with_momentum_forcing(
                grid,
                ocean;
                advection,
                dynamics,
                forcing = (u = u_sponge, v = v_sponge),
                boundary_conditions = isnothing(coastal_boundary_conditions) ?
                    NamedTuple() : coastal_boundary_conditions,
                common_kwargs...,
            )
        end
        if !isnothing(coastal_boundary_conditions)
            return _sea_ice_simulation_with_momentum_forcing(
                grid,
                ocean;
                advection,
                dynamics,
                forcing = NamedTuple(),
                boundary_conditions = coastal_boundary_conditions,
                common_kwargs...,
            )
        end
        return ClimaOcean.sea_ice_simulation(
            grid,
            ocean;
            advection,
            dynamics,
            common_kwargs...,
        )
    end
    return ClimaOcean.sea_ice_simulation(
        grid,
        ocean;
        dynamics = nothing,
        advection = nothing,
        common_kwargs...,
    )
end

function _set_analytic_ocean_ice_initial_conditions!(ocean, sea_ice)
    ocean_temperature(λ, φ, z) = max(-1.5, 18 - 0.3 * abs(φ) + 0.002 * z)
    Oceananigans.set!(ocean.model.tracers.T, ocean_temperature)
    Oceananigans.set!(ocean.model.tracers.S, 35)

    ice_thickness(λ, φ) = abs(φ) > 60 ? 1.0 : 0.0
    ice_concentration(λ, φ) = abs(φ) > 60 ? 0.8 : 0.0
    Oceananigans.set!(sea_ice.model; h = ice_thickness, ℵ = ice_concentration)
    return nothing
end

"""
Set geographic eastward/northward ocean velocities on the model's native grid.

ECCO4's `EVEL` and `NVEL` are extrinsic geographic components.  A tripolar
grid's native `u`/`v` axes rotate strongly in the Arctic, so setting the two
staggered fields independently is incorrect.  Oceananigans' vector-aware
initializer first interpolates both components to cell centres, rotates them
together into intrinsic grid axes, and then interpolates them to their native
staggered locations.

Perform that one-off transformation on a CPU copy of the grid and transfer the
completed native fields to the model architecture.  Launching the two generic
tripolar transformation kernels directly on CUDA can leave the driver in
first-use module/JIT work for hours when every production job starts with a
fresh writable depot.  CPU staging preserves Oceananigans' reference
mathematics while keeping initialization outside the production timestep's GPU
specialization barrier.
"""
function _set_extrinsic_ocean_velocity_fields!(
    u_velocity,
    v_velocity,
    grid,
    eastward,
    northward,
)
    cpu_grid = Oceananigans.on_architecture(
        Oceananigans.CPU(),
        grid,
    )
    cpu_velocities = (
        u = Oceananigans.XFaceField(cpu_grid),
        v = Oceananigans.YFaceField(cpu_grid),
    )
    Oceananigans.Models.HydrostaticFreeSurfaceModels.set_from_extrinsic_velocities!(
        cpu_velocities,
        cpu_grid,
        eastward,
        northward,
    )
    Oceananigans.set!(
        u_velocity,
        Array(Oceananigans.interior(cpu_velocities.u)),
    )
    Oceananigans.set!(
        v_velocity,
        Array(Oceananigans.interior(cpu_velocities.v)),
    )
    # Preserve Oceananigans' reference initializer semantics here: it writes
    # the native staggered interiors but deliberately does not refill their
    # halos. On a tripolar grid, an additional generic tuple halo fill also
    # rewrites redundant fold/periodic interior faces. The owning model's
    # subsequent reconciliation performs the correctly ordered halo update.
    Oceananigans.Architectures.synchronize(
        Oceananigans.Architectures.architecture(grid),
    )
    return nothing
end

function _set_extrinsic_ocean_velocities!(ocean_model, eastward, northward)
    return _set_extrinsic_ocean_velocity_fields!(
        ocean_model.velocities.u,
        ocean_model.velocities.v,
        ocean_model.grid,
        eastward,
        northward,
    )
end

"""
    _set_surface_metadatum!(target, metadata)

Interpolate a two-dimensional geophysical `Metadatum` onto a single-layer
surface field.
NumericalEarth's mutable-grid interpolation correctly uses deformed vertical
coordinates for three-dimensional tracers, but currently calls `znode` with a
`Nothing` vertical location for two-dimensional fields. Surface data have no
vertical coordinate to deform, so Oceananigans' ordinary horizontal
interpolator is the correct path here.
"""
function _set_surface_metadatum!(
    target,
    metadata;
    nonfinite_replacement = nothing,
)
    size(target, 3) == 1 || throw(ArgumentError(
        "_set_surface_metadatum! requires a single-layer target field",
    ))
    architecture = Oceananigans.architecture(target.grid)
    native = Oceananigans.Field(metadata, architecture)
    if !isnothing(nonfinite_replacement)
        data = parent(native)
        replacement = convert(eltype(data), nonfinite_replacement)
        data .= ifelse.(isfinite.(data), data, replacement)
    end
    Oceananigans.Fields.interpolate!(target, native)
    Oceananigans.fill_halo_regions!(target)
    return target
end

function _ecco_ocean_ice_metadata(
    ;
    full_state = false,
    date = DateTime(1993, 1, 1),
    directory = "",
)
    dataset = NumericalEarth.ECCO4Monthly()
    metadatum(name) = isempty(directory) ?
        NumericalEarth.Metadatum(name; date, dataset) :
        NumericalEarth.Metadatum(name; date, dataset, dir = directory)
    return (
        temperature = metadatum(:temperature),
        salinity = metadatum(:salinity),
        ice_thickness = metadatum(:sea_ice_thickness),
        ice_concentration = metadatum(:sea_ice_concentration),
        u_velocity = full_state ? metadatum(:u_velocity) : nothing,
        v_velocity = full_state ? metadatum(:v_velocity) : nothing,
        free_surface = full_state ? metadatum(:free_surface) : nothing,
    )
end

function _set_ecco_ocean_ice_initial_conditions!(
    ocean,
    sea_ice;
    full_state = false,
    date = DateTime(1993, 1, 1),
    directory = "",
)
    metadata = _ecco_ocean_ice_metadata(
        ; full_state,
        date,
        directory,
    )

    # For legacy modes ClimaOcean first attempts the native ECCO endpoint and
    # then its public NumericalEarthArtifacts mirror. Generic monthly modes use
    # the caller's already verified local bundle, so compute nodes need no data
    # credentials or network access.
    foreach(
        ClimaOcean.download_with_fallback,
        filter(metadata_value -> !isnothing(metadata_value), values(metadata)),
    )
    Oceananigans.set!(ocean.model.tracers.T, metadata.temperature)
    Oceananigans.set!(ocean.model.tracers.S, metadata.salinity)
    if full_state
        _set_extrinsic_ocean_velocities!(
            ocean.model,
            metadata.u_velocity,
            metadata.v_velocity,
        )
        _set_surface_metadatum!(
            ocean.model.free_surface.displacement,
            metadata.free_surface,
        )
    end
    # ECCO sea-ice payloads are undefined over land. Missing source ice is
    # physically zero, rather than a value to inpaint across coastlines.
    _set_surface_metadatum!(
        sea_ice.model.ice_thickness,
        metadata.ice_thickness;
        nonfinite_replacement = 0,
    )
    _set_surface_metadatum!(
        sea_ice.model.ice_concentration,
        metadata.ice_concentration,
        nonfinite_replacement = 0,
    )
    _initialize_conditional_ice_from_effective_thickness!(sea_ice.model)
    return nothing
end

function _set_ocean_ice_initial_conditions!(ocean, sea_ice, config::ExperimentConfig)
    if _is_ecco_initial_conditions(config.ocean_initial_conditions)
        return _set_ecco_ocean_ice_initial_conditions!(
            ocean,
            sea_ice;
            full_state = _is_full_state_ecco_initial_conditions(
                config.ocean_initial_conditions,
            ),
            date = _ecco_initial_condition_date(config),
            directory = config.ecco_initial_conditions_directory,
        )
    end
    return _set_analytic_ocean_ice_initial_conditions!(ocean, sea_ice)
end

"""
Complete ordered initialization of the mutable z-star ocean state.

Field assignment, halo filling, staggered z-star scaling, and barotropic
velocity reconstruction are separate asynchronous GPU operations. On the
production tripolar grid a single pass can therefore derive geometry before
the preceding displacement-halo update has completed. Three synchronized
normal initialization passes are required: the first completes the assigned
state, the second stabilizes staggered geometry, and the third reconstructs
barotropic and diagnostic fields from that stable geometry.
"""
function _complete_ocean_initialization!(ocean)
    arch = Oceananigans.Architectures.architecture(ocean.model)
    Oceananigans.Architectures.synchronize(arch)
    for _ in 1:3
        Oceananigans.initialize!(ocean)
        Oceananigans.Architectures.synchronize(arch)
    end
    return nothing
end

"""State changes and time-integrated boundary fluxes for conservation diagnosis."""
Base.@kwdef mutable struct GlobalCoupledBudgetCallback{S, F, D} <: Function
    sample_every_n_steps::Int
    state_integrals::S
    flux_integrals::F
    diagnostic_fields::D
    cumulative_fluxes::Vector{Float64}
    previous_fluxes::Vector{Float64}
    last_time::Float64 = NaN
    time_days::Vector{Float64} = Float64[]
    balanced_launch_applied_fraction::Vector{Float64} = Float64[]
    balanced_launch_cumulative_sources::Vector{Float64} = Float64[]
    ocean_temperature_integral::Vector{Float64} = Float64[]
    ocean_salinity_integral::Vector{Float64} = Float64[]
    ocean_volume::Vector{Float64} = Float64[]
    ocean_mean_temperature_profile::Vector{Float64} = Float64[]
    ocean_mean_salinity_profile::Vector{Float64} = Float64[]
    ocean_temperature_minimum_series::Vector{Float64} = Float64[]
    ocean_temperature_maximum_series::Vector{Float64} = Float64[]
    ocean_salinity_minimum_series::Vector{Float64} = Float64[]
    ocean_salinity_maximum_series::Vector{Float64} = Float64[]
    ocean_surface_salinity_minimum_series::Vector{Float64} = Float64[]
    ocean_surface_salinity_minimum_longitude::Vector{Float64} = Float64[]
    ocean_surface_salinity_minimum_latitude::Vector{Float64} = Float64[]
    ocean_surface_salinity_minimum_i::Vector{Int64} = Int64[]
    ocean_surface_salinity_minimum_j::Vector{Int64} = Int64[]
    ocean_runoff_receiver_surface_salinity_minimum_series::Vector{Float64} = Float64[]
    ocean_runoff_receiver_surface_salinity_minimum_longitude::Vector{Float64} = Float64[]
    ocean_runoff_receiver_surface_salinity_minimum_latitude::Vector{Float64} = Float64[]
    ocean_runoff_receiver_surface_salinity_minimum_i::Vector{Int64} = Int64[]
    ocean_runoff_receiver_surface_salinity_minimum_j::Vector{Int64} = Int64[]
    ocean_river_mouth_active_mixing_cells::Vector{Int64} = Int64[]
    ocean_river_mouth_active_mixing_area_m2::Vector{Float64} = Float64[]
    sea_ice_volume::Vector{Float64} = Float64[]
    sea_ice_area::Vector{Float64} = Float64[]
    tropical_sea_ice_area::Vector{Float64} = Float64[]
    # Retain the raw conditional maximum for anomaly visibility, but diagnose
    # physical pack thickness only within the conventional 15% ice extent.
    sea_ice_max_thickness::Vector{Float64} = Float64[]
    sea_ice_extent_max_thickness::Vector{Float64} = Float64[]
    sea_ice_max_area_equivalent_thickness::Vector{Float64} = Float64[]
    sea_ice_marginal_ceiling_cell_count::Vector{Int64} = Int64[]
    sea_ice_marginal_ceiling_max_concentration::Vector{Float64} = Float64[]
    sea_ice_marginal_ceiling_max_area_equivalent_thickness::Vector{Float64} = Float64[]
    ocean_free_surface_minimum::Vector{Float64} = Float64[]
    ocean_free_surface_maximum::Vector{Float64} = Float64[]
    ocean_max_abs_zonal_velocity::Vector{Float64} = Float64[]
    ocean_max_abs_meridional_velocity::Vector{Float64} = Float64[]
    cumulative_ocean_surface_heat::Vector{Float64} = Float64[]
    cumulative_ocean_frazil_heat::Vector{Float64} = Float64[]
    cumulative_ocean_penetrating_shortwave_heat::Vector{Float64} = Float64[]
    cumulative_net_ocean_heat::Vector{Float64} = Float64[]
    cumulative_ocean_freshwater_enthalpy::Vector{Float64} = Float64[]
    cumulative_net_ocean_freshwater_volume::Vector{Float64} = Float64[]
    cumulative_net_ocean_salinity::Vector{Float64} = Float64[]
    cumulative_rain_mass::Vector{Float64} = Float64[]
    cumulative_snow_mass::Vector{Float64} = Float64[]
    cumulative_open_ocean_evaporation_mass::Vector{Float64} = Float64[]
    cumulative_sea_ice_ocean_freshwater_volume::Vector{Float64} = Float64[]
    cumulative_routed_land_runoff_mass::Vector{Float64} = Float64[]
end

function _budget_integral(field)
    field = field isa Oceananigans.Fields.AbstractField ? field : Oceananigans.Field(field)
    return Oceananigans.Field(Oceananigans.Integral(field))
end

function _budget_scalar(field)
    Oceananigans.compute!(field)
    return Float64(only(Array(Oceananigans.interior(field))))
end

_budget_scalar(operation::Oceananigans.AbstractOperations.KernelFunctionOperation) =
    Float64(sum(operation))

_profile_budget_construction() =
    lowercase(get(ENV, "READYESM_PROFILE_BUDGET_CONSTRUCTION", "false")) in
    ("1", "true", "yes", "on")

function _profile_budget_phase(function_, architecture, label)
    _profile_budget_construction() || return function_()
    Oceananigans.Architectures.synchronize(architecture)
    start = time_ns()
    value = function_()
    Oceananigans.Architectures.synchronize(architecture)
    seconds = (time_ns() - start) / 1e9
    println("BUDGET_CONSTRUCTION_PHASE label=$label seconds=$seconds")
    flush(stdout)
    return value
end

@inline function _float64_cell_volume(i, j, k, grid)
    return Float64(Oceananigans.volume(
        i,
        j,
        k,
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    ))
end

@inline function _float64_tracer_content(i, j, k, grid, tracer)
    return Float64(@inbounds tracer[i, j, k]) *
        _float64_cell_volume(i, j, k, grid)
end

function _float64_ocean_volume_operation(grid)
    return Oceananigans.KernelFunctionOperation{
        Oceananigans.Center,
        Oceananigans.Center,
        Oceananigans.Center,
    }(_float64_cell_volume, grid, (), Float64)
end

function _float64_tracer_integral_operation(tracer)
    return Oceananigans.KernelFunctionOperation{
        Oceananigans.Center,
        Oceananigans.Center,
        Oceananigans.Center,
    }(_float64_tracer_content, tracer.grid, (tracer,), Float64)
end

function GlobalCoupledBudgetCallback(
    earth,
    Δt;
    sample_every_n_steps = max(1, round(Int, 86_400 / Float64(Δt))),
)
    ocean = earth.ocean.model
    sea_ice = earth.sea_ice.model
    architecture = Oceananigans.architecture(ocean)
    exchange = earth.interfaces.exchanger.atmosphere.state
    atmosphere_ocean = earth.interfaces.atmosphere_ocean_interface.fluxes
    sea_ice_ocean = earth.interfaces.sea_ice_ocean_interface.fluxes
    net_ocean = earth.interfaces.net_fluxes.ocean
    concentration = sea_ice.ice_concentration
    tropical_mask = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(ocean.grid)
    _profile_budget_phase(architecture, "tropical_mask") do
        Oceananigans.set!(
            tropical_mask,
            (λ, φ) ->
                abs(φ) <=
                _TROPICAL_SEA_ICE_MAXIMUM_ABSOLUTE_LATITUDE_DEGREES ? 1 : 0,
        )
    end

    state_integrals = _profile_budget_phase(
        architecture,
        "state_integral_operations",
    ) do
        (
            # The prognostic fields are Float32, but their global inventories are
            # O(1e25 J) and O(1e19 kg). Reducing in Float32 quantizes physically
            # relevant timestep changes. Accumulate cell content in Float64.
            ocean_temperature = _float64_tracer_integral_operation(ocean.tracers.T),
            ocean_salinity = _float64_tracer_integral_operation(ocean.tracers.S),
            # Follow NumericalEarth's conservation reference exactly: on a z-star
            # grid the conserved ocean volume is the sum of the live cell volumes,
            # not the area integral of the free-surface solver's displacement.
            ocean_volume = _float64_ocean_volume_operation(ocean.grid),
            sea_ice_volume = _budget_integral(
                Oceananigans.Field(sea_ice.ice_thickness * concentration),
            ),
            sea_ice_area = _budget_integral(concentration),
            tropical_sea_ice_area = _budget_integral(
                Oceananigans.Field(tropical_mask * concentration),
            ),
        )
    end
    ocean_reference_density = earth.interfaces.ocean_properties.reference_density
    ocean_heat_capacity = earth.interfaces.ocean_properties.heat_capacity
    # `net_ocean.T` is the ordinary top boundary temperature flux. Frazil is
    # different: NumericalEarth clips sub-freezing ocean cells while computing
    # the *new* interface state, so its diagnosed heat flux applies to the state
    # change of the just-completed step rather than the following one.
    flux_fields = _profile_budget_phase(architecture, "flux_fields") do
        ocean_surface_heat = Oceananigans.Field(
            ocean_reference_density * ocean_heat_capacity * net_ocean.T,
        )
        ocean_frazil_heat = sea_ice_ocean.frazil_heat
        # With NumericalEarth's default TwoColorRadiation, transmitted shortwave
        # is removed from the surface boundary flux and deposited as a volumetric
        # ocean temperature forcing. In the SpeedyWeather-driven interface used
        # here this diagnostic is positive downward (into the ocean), opposite to
        # the positive-upward convention of the ordinary surface heat flux.
        ocean_penetrating_shortwave_heat =
            earth.radiation.interface_fluxes.ocean.downwelling_shortwave
        freshwater_enthalpy = Oceananigans.Field(
            ocean_reference_density * ocean_heat_capacity *
            net_ocean.freshwater_heat_content,
        )
        open_ocean_evaporation = Oceananigans.Field(
            (1 - concentration) * atmosphere_ocean.water_vapor,
        )
        land_exchanger = earth.interfaces.exchanger.land
        routed_land_runoff = if isnothing(land_exchanger)
            Oceananigans.Field(0 * net_ocean.η)
        else
            land_exchanger.state.freshwater_flux
        end
        (;
            ocean_surface_heat,
            ocean_frazil_heat,
            ocean_penetrating_shortwave_heat,
            freshwater_enthalpy,
            open_ocean_evaporation,
            routed_land_runoff,
        )
    end
    flux_integrals = _profile_budget_phase(
        architecture,
        "flux_integral_operations",
    ) do
        (
            ocean_surface_heat = _budget_integral(flux_fields.ocean_surface_heat),
            ocean_frazil_heat = _budget_integral(flux_fields.ocean_frazil_heat),
            ocean_penetrating_shortwave_heat = _budget_integral(
                flux_fields.ocean_penetrating_shortwave_heat,
            ),
            ocean_freshwater_enthalpy = _budget_integral(
                flux_fields.freshwater_enthalpy,
            ),
            net_ocean_freshwater_volume = _budget_integral(net_ocean.η),
            # On the mutable z-star grid, the live S * Jw boundary term cancels
            # salinity carried by grid-volume motion. Total salt content therefore
            # closes against this base (sea-ice) salt flux alone.
            net_ocean_salinity = _budget_integral(net_ocean.S),
            rain_mass = _budget_integral(exchange.Jʳⁿ),
            snow_mass = _budget_integral(exchange.Jˢⁿ),
            open_ocean_evaporation_mass = _budget_integral(
                flux_fields.open_ocean_evaporation,
            ),
            sea_ice_ocean_freshwater_volume = _budget_integral(
                sea_ice_ocean.freshwater,
            ),
            routed_land_runoff_mass = _budget_integral(
                flux_fields.routed_land_runoff,
            ),
        )
    end
    ocean_surface_longitude, ocean_surface_latitude =
        _ocean_surface_coordinate_matrices(ocean.tracers.S)
    ocean_runoff_receiver_surface_cells = if earth.land isa TerrariumRunoffLand
        target_i = Int.(Array(earth.land.target_i))
        target_j = Int.(Array(earth.land.target_j))
        length(target_i) == length(target_j) || throw(DimensionMismatch(
            "runoff receiver index vectors have different lengths",
        ))
        CartesianIndex.(target_i, target_j)
    else
        CartesianIndex{2}[]
    end
    river_mouth_mixing = earth.land isa TerrariumRunoffLand ?
        earth.land.river_mouth_mixing : nothing
    diagnostic_fields = (
        ocean_temperature = ocean.tracers.T,
        ocean_salinity = ocean.tracers.S,
        ocean_mean_temperature_profile = Oceananigans.Field(
            Oceananigans.Average(ocean.tracers.T; dims = (1, 2)),
        ),
        ocean_mean_salinity_profile = Oceananigans.Field(
            Oceananigans.Average(ocean.tracers.S; dims = (1, 2)),
        ),
        sea_ice_thickness = sea_ice.ice_thickness,
        sea_ice_concentration = concentration,
        ocean_free_surface = ocean.free_surface.displacement,
        ocean_zonal_velocity = ocean.velocities.u,
        ocean_meridional_velocity = ocean.velocities.v,
        ocean_surface_active_cells = findall(
            _ocean_surface_active_mask(ocean.grid) .> 0,
        ),
        ocean_surface_longitude,
        ocean_surface_latitude,
        ocean_runoff_receiver_surface_cells,
        ocean_routed_land_runoff = flux_fields.routed_land_runoff,
        ocean_river_mouth_active_mixing_mask =
            isnothing(river_mouth_mixing) ? nothing :
            river_mouth_mixing.active_receiver_mask,
        ocean_river_mouth_dynamically_gated =
            !isnothing(river_mouth_mixing) &&
            river_mouth_mixing.dynamically_gated,
        ocean_river_mouth_reference_freshwater_mass_flux_kgm2s =
            isnothing(river_mouth_mixing) ? 0.0 : Float64(
                river_mouth_mixing.reference_freshwater_mass_flux_kgm2s,
            ),
        ocean_surface_cell_area = _ocean_surface_cell_area(ocean.grid),
    )
    cumulative_fluxes = zeros(Float64, length(flux_integrals))
    previous_fluxes = zeros(Float64, length(flux_integrals))
    return GlobalCoupledBudgetCallback(
        ; sample_every_n_steps,
        state_integrals,
        flux_integrals,
        diagnostic_fields,
        cumulative_fluxes,
        previous_fluxes,
    )
end

@inline _net_ocean_heat(cumulative) = cumulative[1] + cumulative[2] - cumulative[3]

const _SEA_ICE_EXTENT_CONCENTRATION_THRESHOLD = 0.15

function _sea_ice_thickness_extrema(
    thickness,
    concentration;
    extent_threshold = _SEA_ICE_EXTENT_CONCENTRATION_THRESHOLD,
    representation_ceiling = 15.0,
)
    size(thickness) == size(concentration) || throw(
        DimensionMismatch("sea-ice thickness and concentration shapes differ"),
    )
    isempty(thickness) && return (
        raw = 0.0,
        within_extent = 0.0,
        area_equivalent = 0.0,
        marginal_ceiling_cell_count = 0,
        marginal_ceiling_max_concentration = 0.0,
        marginal_ceiling_max_area_equivalent_thickness = 0.0,
    )
    raw = maximum(thickness)
    area_equivalent = maximum(thickness .* concentration)
    within_extent = 0.0
    marginal_ceiling_cell_count = 0
    marginal_ceiling_max_concentration = 0.0
    marginal_ceiling_max_area_equivalent_thickness = 0.0
    for index in eachindex(thickness, concentration)
        if concentration[index] >= extent_threshold
            within_extent = max(within_extent, Float64(thickness[index]))
        end
        at_representation_ceiling =
            thickness[index] >= 0.999 * representation_ceiling &&
            concentration[index] < extent_threshold
        if at_representation_ceiling
            marginal_ceiling_cell_count += 1
            marginal_ceiling_max_concentration = max(
                marginal_ceiling_max_concentration,
                Float64(concentration[index]),
            )
            marginal_ceiling_max_area_equivalent_thickness = max(
                marginal_ceiling_max_area_equivalent_thickness,
                Float64(thickness[index] * concentration[index]),
            )
        end
    end
    return (;
        raw,
        within_extent,
        area_equivalent,
        marginal_ceiling_cell_count,
        marginal_ceiling_max_concentration,
        marginal_ceiling_max_area_equivalent_thickness,
    )
end

function _sample_coupled_budget_state!(callback, simulation)
    state = callback.state_integrals
    diagnostic = callback.diagnostic_fields
    push!(callback.time_days, Float64(simulation.model.clock.time) / 86_400)
    launch = balanced_launch_diagnostics(simulation.model)
    push!(
        callback.balanced_launch_applied_fraction,
        isnothing(launch) ? 0.0 : Float64(launch.applied_fraction),
    )
    append!(
        callback.balanced_launch_cumulative_sources,
        _balanced_launch_cumulative_source_vector(simulation.model),
    )
    push!(callback.ocean_temperature_integral, _budget_scalar(state.ocean_temperature))
    push!(callback.ocean_salinity_integral, _budget_scalar(state.ocean_salinity))
    push!(callback.ocean_volume, _budget_scalar(state.ocean_volume))
    Oceananigans.compute!(diagnostic.ocean_mean_temperature_profile)
    Oceananigans.compute!(diagnostic.ocean_mean_salinity_profile)
    append!(
        callback.ocean_mean_temperature_profile,
        vec(Float64.(Array(Oceananigans.interior(
            diagnostic.ocean_mean_temperature_profile,
        )))),
    )
    append!(
        callback.ocean_mean_salinity_profile,
        vec(Float64.(Array(Oceananigans.interior(
            diagnostic.ocean_mean_salinity_profile,
        )))),
    )
    temperature_extrema = _active_ocean_extrema(diagnostic.ocean_temperature)
    salinity_extrema = _active_ocean_extrema(diagnostic.ocean_salinity)
    push!(callback.ocean_temperature_minimum_series, first(temperature_extrema))
    push!(callback.ocean_temperature_maximum_series, last(temperature_extrema))
    push!(callback.ocean_salinity_minimum_series, first(salinity_extrema))
    push!(callback.ocean_salinity_maximum_series, last(salinity_extrema))
    surface_salinity = Float64.(Array(@view Oceananigans.interior(
        diagnostic.ocean_salinity,
    )[:, :, end]))
    active_surface_cells = diagnostic.ocean_surface_active_cells
    isempty(active_surface_cells) && error("ocean has no active surface cells")
    minimum_surface_cell = active_surface_cells[
        argmin(surface_salinity[active_surface_cells])
    ]
    minimum_i, minimum_j = Tuple(minimum_surface_cell)
    push!(
        callback.ocean_surface_salinity_minimum_series,
        surface_salinity[minimum_surface_cell],
    )
    push!(
        callback.ocean_surface_salinity_minimum_longitude,
        diagnostic.ocean_surface_longitude[minimum_surface_cell],
    )
    push!(
        callback.ocean_surface_salinity_minimum_latitude,
        diagnostic.ocean_surface_latitude[minimum_surface_cell],
    )
    push!(callback.ocean_surface_salinity_minimum_i, minimum_i)
    push!(callback.ocean_surface_salinity_minimum_j, minimum_j)
    runoff_receiver_cells = diagnostic.ocean_runoff_receiver_surface_cells
    if isempty(runoff_receiver_cells)
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_series,
            NaN,
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_longitude,
            NaN,
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_latitude,
            NaN,
        )
        push!(callback.ocean_runoff_receiver_surface_salinity_minimum_i, 0)
        push!(callback.ocean_runoff_receiver_surface_salinity_minimum_j, 0)
    else
        minimum_runoff_receiver = runoff_receiver_cells[
            argmin(surface_salinity[runoff_receiver_cells])
        ]
        runoff_i, runoff_j = Tuple(minimum_runoff_receiver)
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_series,
            surface_salinity[minimum_runoff_receiver],
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_longitude,
            diagnostic.ocean_surface_longitude[minimum_runoff_receiver],
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_latitude,
            diagnostic.ocean_surface_latitude[minimum_runoff_receiver],
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_i,
            runoff_i,
        )
        push!(
            callback.ocean_runoff_receiver_surface_salinity_minimum_j,
            runoff_j,
        )
    end
    active_mixing_field = diagnostic.ocean_river_mouth_active_mixing_mask
    if isnothing(active_mixing_field)
        push!(callback.ocean_river_mouth_active_mixing_cells, 0)
        push!(callback.ocean_river_mouth_active_mixing_area_m2, 0.0)
    else
        active_mixing_mask = Float64.(Array(Oceananigans.interior(
            active_mixing_field,
        )))[:, :, 1]
        all(value -> value == 0 || value == 1, active_mixing_mask) || error(
            "river-mouth active mixing mask is not binary",
        )
        if diagnostic.ocean_river_mouth_dynamically_gated
            routed_runoff = Float64.(Array(Oceananigans.interior(
                diagnostic.ocean_routed_land_runoff,
            )))[:, :, 1]
            all(routed_runoff .>= 0) || error(
                "routed land runoff contains a negative freshwater flux",
            )
            reference_flux = diagnostic.
                ocean_river_mouth_reference_freshwater_mass_flux_kgm2s
            expected_active = Float64.(_river_mouth_mixing_active.(
                routed_runoff,
                reference_flux,
            ))
            if active_mixing_mask != expected_active
                mismatch_count = count(active_mixing_mask .!= expected_active)
                error(
                    "active river-mouth mixing does not match its routed-runoff " *
                    "loading law: mismatches=$mismatch_count, " *
                    "reference_flux_kgm2s=$reference_flux",
                )
            end
        end
        push!(
            callback.ocean_river_mouth_active_mixing_cells,
            count(!iszero, active_mixing_mask),
        )
        push!(
            callback.ocean_river_mouth_active_mixing_area_m2,
            sum(diagnostic.ocean_surface_cell_area[active_mixing_mask .== 1]),
        )
    end
    push!(callback.sea_ice_volume, _budget_scalar(state.sea_ice_volume))
    push!(callback.sea_ice_area, _budget_scalar(state.sea_ice_area))
    push!(callback.tropical_sea_ice_area, _budget_scalar(state.tropical_sea_ice_area))
    sea_ice_thickness = Float64.(Array(Oceananigans.interior(
        diagnostic.sea_ice_thickness,
    )))
    sea_ice_concentration = Float64.(Array(Oceananigans.interior(
        diagnostic.sea_ice_concentration,
    )))
    ocean_free_surface = Float64.(Array(Oceananigans.interior(
        diagnostic.ocean_free_surface,
    )))
    ocean_zonal_velocity = Float64.(Array(Oceananigans.interior(
        diagnostic.ocean_zonal_velocity,
    )))
    ocean_meridional_velocity = Float64.(Array(Oceananigans.interior(
        diagnostic.ocean_meridional_velocity,
    )))
    sea_ice_extrema = _sea_ice_thickness_extrema(
        sea_ice_thickness,
        sea_ice_concentration,
    )
    push!(callback.sea_ice_max_thickness, sea_ice_extrema.raw)
    push!(
        callback.sea_ice_extent_max_thickness,
        sea_ice_extrema.within_extent,
    )
    push!(
        callback.sea_ice_max_area_equivalent_thickness,
        sea_ice_extrema.area_equivalent,
    )
    push!(
        callback.sea_ice_marginal_ceiling_cell_count,
        sea_ice_extrema.marginal_ceiling_cell_count,
    )
    push!(
        callback.sea_ice_marginal_ceiling_max_concentration,
        sea_ice_extrema.marginal_ceiling_max_concentration,
    )
    push!(
        callback.sea_ice_marginal_ceiling_max_area_equivalent_thickness,
        sea_ice_extrema.marginal_ceiling_max_area_equivalent_thickness,
    )
    push!(callback.ocean_free_surface_minimum, minimum(ocean_free_surface))
    push!(callback.ocean_free_surface_maximum, maximum(ocean_free_surface))
    push!(
        callback.ocean_max_abs_zonal_velocity,
        maximum(abs, ocean_zonal_velocity),
    )
    push!(
        callback.ocean_max_abs_meridional_velocity,
        maximum(abs, ocean_meridional_velocity),
    )
    cumulative = callback.cumulative_fluxes
    push!(callback.cumulative_ocean_surface_heat, cumulative[1])
    push!(callback.cumulative_ocean_frazil_heat, cumulative[2])
    push!(callback.cumulative_ocean_penetrating_shortwave_heat, cumulative[3])
    push!(callback.cumulative_net_ocean_heat, _net_ocean_heat(cumulative))
    for (name, value) in zip((
        :cumulative_ocean_freshwater_enthalpy,
        :cumulative_net_ocean_freshwater_volume,
        :cumulative_net_ocean_salinity,
        :cumulative_rain_mass,
        :cumulative_snow_mass,
        :cumulative_open_ocean_evaporation_mass,
        :cumulative_sea_ice_ocean_freshwater_volume,
        :cumulative_routed_land_runoff_mass,
    ), cumulative[4:end])
        push!(getproperty(callback, name), value)
    end
    return nothing
end

function Oceananigans.initialize!(callback::GlobalCoupledBudgetCallback, simulation)
    for name in (
        :time_days,
        :balanced_launch_applied_fraction,
        :balanced_launch_cumulative_sources,
        :ocean_temperature_integral,
        :ocean_salinity_integral,
        :ocean_volume,
        :ocean_mean_temperature_profile,
        :ocean_mean_salinity_profile,
        :ocean_temperature_minimum_series,
        :ocean_temperature_maximum_series,
        :ocean_salinity_minimum_series,
        :ocean_salinity_maximum_series,
        :ocean_surface_salinity_minimum_series,
        :ocean_surface_salinity_minimum_longitude,
        :ocean_surface_salinity_minimum_latitude,
        :ocean_surface_salinity_minimum_i,
        :ocean_surface_salinity_minimum_j,
        :ocean_runoff_receiver_surface_salinity_minimum_series,
        :ocean_runoff_receiver_surface_salinity_minimum_longitude,
        :ocean_runoff_receiver_surface_salinity_minimum_latitude,
        :ocean_runoff_receiver_surface_salinity_minimum_i,
        :ocean_runoff_receiver_surface_salinity_minimum_j,
        :ocean_river_mouth_active_mixing_cells,
        :ocean_river_mouth_active_mixing_area_m2,
        :sea_ice_volume,
        :sea_ice_area,
        :tropical_sea_ice_area,
        :sea_ice_max_thickness,
        :sea_ice_extent_max_thickness,
        :sea_ice_max_area_equivalent_thickness,
        :sea_ice_marginal_ceiling_cell_count,
        :sea_ice_marginal_ceiling_max_concentration,
        :sea_ice_marginal_ceiling_max_area_equivalent_thickness,
        :ocean_free_surface_minimum,
        :ocean_free_surface_maximum,
        :ocean_max_abs_zonal_velocity,
        :ocean_max_abs_meridional_velocity,
        :cumulative_ocean_surface_heat,
        :cumulative_ocean_frazil_heat,
        :cumulative_ocean_penetrating_shortwave_heat,
        :cumulative_net_ocean_heat,
        :cumulative_ocean_freshwater_enthalpy,
        :cumulative_net_ocean_freshwater_volume,
        :cumulative_net_ocean_salinity,
        :cumulative_rain_mass,
        :cumulative_snow_mass,
        :cumulative_open_ocean_evaporation_mass,
        :cumulative_sea_ice_ocean_freshwater_volume,
        :cumulative_routed_land_runoff_mass,
    )
        empty!(getproperty(callback, name))
    end
    fill!(callback.cumulative_fluxes, 0)
    fill!(callback.previous_fluxes, 0)
    callback.last_time = NaN
    return nothing
end

@inline _applied_budget_flux(index, previous, current) =
    index == 2 ? current : previous

function _accumulate_budget_fluxes!(cumulative, previous, current, Δt)
    length(cumulative) == length(previous) == length(current) || throw(
        DimensionMismatch("coupled budget flux vectors differ in length"),
    )
    for index in eachindex(cumulative)
        cumulative[index] += Δt * _applied_budget_flux(
            index,
            previous[index],
            current[index],
        )
    end
    copyto!(previous, current)
    return cumulative
end

function (callback::GlobalCoupledBudgetCallback)(simulation)
    time = Float64(simulation.model.clock.time)
    current_fluxes = [_budget_scalar(field) for field in callback.flux_integrals]
    # Oceananigans fires scheduled callbacks once at initialization, after the
    # coupled model has assembled the time-zero interface state. Seed both the
    # state and rectangle-at-start flux history there; no elapsed interval is
    # accumulated at t=0.
    if !isfinite(callback.last_time)
        copyto!(callback.previous_fluxes, current_fluxes)
        callback.last_time = time
        _sample_coupled_budget_state!(callback, simulation)
        return nothing
    end

    Δt = time - callback.last_time
    Δt >= 0 || error("coupled budget callback time moved backwards")
    # Interface fluxes assembled at the end of step n force the ocean over
    # (n, n+1]. Frazil heat is the exception: it diagnoses and repairs the
    # just-completed state during update_state!, so use its current value.
    _accumulate_budget_fluxes!(
        callback.cumulative_fluxes,
        callback.previous_fluxes,
        current_fluxes,
        Δt,
    )
    callback.last_time = time
    iteration = simulation.model.clock.iteration
    at_sample = iteration % callback.sample_every_n_steps == 0
    at_end = iteration >= simulation.stop_iteration ||
             simulation.model.clock.time >= simulation.stop_time
    if (at_sample || at_end) &&
       (isempty(callback.time_days) || callback.time_days[end] !=
        Float64(simulation.model.clock.time) / 86_400)
        _sample_coupled_budget_state!(callback, simulation)
    end
    return nothing
end

function build_dynamic_esm(config::ExperimentConfig; steps::Union{Nothing, Int} = nothing)
    validate(config)
    config.ocean_model == :oceananigans || throw(
        ArgumentError("build_dynamic_esm requires ocean_model=oceananigans"),
    )

    atmosphere_model = _rrtmgp_atmosphere_model(config; coupled_ocean = true)
    atmosphere = SpeedyWeather.initialize!(
        atmosphere_model;
        time = config.start_date,
    )
    atmosphere_Δt = Float64(atmosphere.model.time_stepping.Δt_sec)
    atmosphere_stop_time = isnothing(steps) ?
        config.duration_days * 86_400 : steps * atmosphere_Δt
    SpeedyWeather.initialize!(
        atmosphere;
        output = false,
        period = SpeedyWeather.Second(round(Int, atmosphere_stop_time)),
    )

    ocean_build = _build_ocean(config; return_auxiliary = true)
    ocean = ocean_build.ocean
    sea_ice = _build_sea_ice(config, ocean)
    _set_ocean_ice_initial_conditions!(ocean, sea_ice, config)
    _mask_inactive_sea_ice_surface_cells!(sea_ice)
    _complete_ocean_initialization!(ocean)
    land = config.land_model == :terrarium ?
        _build_terrarium_runoff_land(
            atmosphere,
            ocean;
            routing_weighting = config.terrarium_runoff_routing,
            river_mouth_mixing = ocean_build.river_mouth_mixing,
        ) : nothing
    river_mouth_mixing = _initialize_river_mouth_mixing!(
        ocean_build.river_mouth_mixing,
        land,
        ocean.model.grid,
        config,
    )
    isnothing(river_mouth_mixing) ||
        @info "Initialized localized river-mouth mixing" river_mouth_mixing

    radiation = AtmosphereDrivenRadiation(eltype(ocean.model.grid))
    earth_model = ClimaOcean.EarthSystemModel(
        radiation,
        atmosphere,
        land,
        sea_ice,
        ocean,
    )
    Δt = eltype(ocean.model.grid)(atmosphere_Δt)
    simulation = if isnothing(steps)
        Oceananigans.Simulation(
            earth_model;
            Δt,
            stop_time = config.duration_days * 86_400,
        )
    else
        steps > 0 || throw(ArgumentError("steps must be positive"))
        Oceananigans.Simulation(earth_model; Δt, stop_iteration = steps)
    end
    Oceananigans.add_callback!(
        simulation,
        ClimaOcean.Progress(),
        Oceananigans.TimeInterval(86_400),
    )
    Oceananigans.add_callback!(
        simulation,
        GlobalCoupledBudgetCallback(
            earth_model,
            Δt;
            sample_every_n_steps = if !isnothing(steps) &&
                                      steps <= round(Int, 86_400 / Float64(Δt))
                1
            else
                max(1, round(Int, 86_400 / Float64(Δt)))
            end,
        ),
        Oceananigans.IterationInterval(1);
        name = :global_coupled_budgets,
    )
    return simulation
end

"""Run the dynamic ocean/ice configuration and return its coupled simulation."""
function run_dynamic_esm!(config::ExperimentConfig; steps::Union{Nothing, Int} = nothing)
    simulation = build_dynamic_esm(config; steps)
    Oceananigans.run!(simulation)
    return simulation
end

"""Release model-owned global caches after a completed multi-case experiment."""
function _release_dynamic_runtime_caches!(simulation)
    atmosphere = simulation.model.atmosphere
    delete!(
        _INTENSIVE_REGRID_CACHE,
        atmosphere.variables.prognostic.clock,
    )
    delete!(
        _DETERMINISTIC_LEGENDRE_MMAX_CACHE,
        atmosphere.model.spectral_transform.kjm_indices,
    )
    land = atmosphere.model.land
    if land isa _ReadyESMAbstractTerrariumLandModel
        delete!(_TERRARIUM_LAND_INDEX_CACHE, land.model.grid.mask.data)
    end
    return nothing
end

_interior_matrix(field) = reshape(
    Float64.(Array(Oceananigans.interior(field))),
    size(field.grid, 1),
    size(field.grid, 2),
)

function _surface_grid_matrix(field)
    values = Float64.(Array(Oceananigans.interior(field)))
    ndims(values) == 3 && (values = values[:, :, end])
    Nx, Ny = size(field.grid, 1), size(field.grid, 2)
    return values[1:Nx, 1:Ny]
end

function _ocean_surface_coordinate_matrices(field)
    longitude, latitude, _ = Oceananigans.nodes(field)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    Nx, Ny = size(field.grid, 1), size(field.grid, 2)
    if size(longitude) == (Nx, Ny) && size(latitude) == (Nx, Ny)
        return longitude, latitude
    end
    length(longitude) == Nx && length(latitude) == Ny || throw(
        DimensionMismatch("ocean surface coordinates do not match the field grid"),
    )
    longitude_matrix = repeat(reshape(longitude, Nx, 1), 1, Ny)
    latitude_matrix = repeat(reshape(latitude, 1, Ny), Nx, 1)
    return longitude_matrix, latitude_matrix
end

function _host_vertical_center_nodes(grid)
    nodes = Oceananigans.Grids.znodes(
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )
    host_nodes = Oceananigans.Architectures.on_architecture(
        Oceananigans.CPU(),
        nodes,
    )
    return Float64.(collect(host_nodes))
end

"""
Return extrema over dynamically active field points.

Oceananigans reductions on an `ImmersedBoundaryGrid` exclude its immersed
periphery.  Reducing the host-side interior array directly does not: dry cells
retain the tracer fill value (zero), which previously made the reported ocean
salinity minimum spuriously equal to 0 psu.
"""
_active_ocean_extrema(field) = Float64.(extrema(field))

"""Return a field mean that excludes immersed ocean points."""
_active_ocean_mean_square(field) = Float64(mean(abs2, field))

function _reshape_profile_history(values, nlayers, nsamples, label)
    expected = nlayers * nsamples
    length(values) == expected || error(
        "$label history has $(length(values)) values, expected $expected",
    )
    return reshape(Float64.(values), nlayers, nsamples)
end

"""Return the configured substep count from an unwrapped ClimaSeaIce model."""
_sea_ice_momentum_substeps(sea_ice) = sea_ice.dynamics.solver.substeps

"""Collect component states and exchanged interface fluxes from a dynamic ESM run."""
function _atmosphere_process_profile_history(callback, name, nlayers)
    values = getproperty(callback, name)
    ntimes = length(callback.time_days)
    length(values) == nlayers * ntimes || error(
        "atmosphere process-profile history $name has $(length(values)) values; " *
        "expected $(nlayers * ntimes)",
    )
    return reshape(Float64.(values), nlayers, ntimes)
end

function collect_dynamic_diagnostics(simulation, config::ExperimentConfig)
    earth = simulation.model
    ocean = earth.ocean.model
    sea_ice = earth.sea_ice.model
    atmosphere = earth.atmosphere
    longitude, latitude, ocean_layer_depth = Oceananigans.nodes(ocean.tracers.T)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    ocean_layer_depth = Float64.(Array(ocean_layer_depth))
    ocean_temperature = Float64.(Array(Oceananigans.interior(ocean.tracers.T)))
    ocean_salinity = Float64.(Array(Oceananigans.interior(ocean.tracers.S)))
    ocean_u = Float64.(Array(Oceananigans.interior(ocean.velocities.u)))
    ocean_v = Float64.(Array(Oceananigans.interior(ocean.velocities.v)))
    ocean_w = Float64.(Array(Oceananigans.interior(ocean.velocities.w)))
    ocean_surface_active_mask = _ocean_surface_active_mask(ocean.grid)
    ocean_cell_area = _ocean_surface_cell_area(ocean.grid)
    runoff_receiver_footprint = if earth.land isa TerrariumRunoffLand
        _river_mouth_surface_footprint(
            earth.land,
            ocean.grid,
            config.ocean_river_mouth_mixing,
        )
    else
        (
            receiver_mask = zeros(Float64, size(ocean_surface_active_mask)),
            horizontal_mask = zeros(Float64, size(ocean_surface_active_mask)),
        )
    end
    river_mouth_footprint = if _is_localized_estuary_mixing(
        config.ocean_river_mouth_mixing,
    )
        earth.land isa TerrariumRunoffLand || error(
            "localized river-mouth diagnostics require Terrarium runoff routing",
        )
        runoff_receiver_footprint
    else
        (
            receiver_mask = zeros(Float64, size(ocean_surface_active_mask)),
            horizontal_mask = zeros(Float64, size(ocean_surface_active_mask)),
        )
    end
    river_mouth_active_mixing_mask = if _is_localized_estuary_mixing(
        config.ocean_river_mouth_mixing,
    )
        earth.land isa TerrariumRunoffLand || error(
            "localized active-mixing diagnostics require Terrarium routing",
        )
        mixing = earth.land.river_mouth_mixing
        isnothing(mixing) && error(
            "localized active-mixing diagnostics lack closure fields",
        )
        active = _surface_grid_matrix(mixing.active_receiver_mask)
        size(active) == size(river_mouth_footprint.receiver_mask) || throw(
            DimensionMismatch(
                "active river-mouth mask does not match the receiver grid",
            ),
        )
        all(value -> value == 0 || value == 1, active) || error(
            "active river-mouth mixing mask is not binary",
        )
        all(active .<= river_mouth_footprint.receiver_mask) || error(
            "active river-mouth mixing extends outside routed receivers",
        )
        active
    else
        zeros(Float64, size(ocean_surface_active_mask))
    end
    ocean_routed_land_runoff_flux = if earth.land isa TerrariumRunoffLand
        _surface_grid_matrix(
            earth.interfaces.exchanger.land.state.freshwater_flux,
        )
    else
        zeros(Float64, size(ocean_surface_active_mask))
    end
    ocean_temperature_extrema = _active_ocean_extrema(ocean.tracers.T)
    ocean_salinity_extrema = _active_ocean_extrema(ocean.tracers.S)
    atmosphere_longitude, atmosphere_latitude = RG.get_londlatds(
        atmosphere.model.spectral_grid.grid,
    )
    atmosphere_land_fraction = Float64.(
        Array(atmosphere.model.land_sea_mask.mask.data),
    )
    land_surface_temperature = _land_surface_field(
        atmosphere.variables.prognostic.land.soil_temperature,
        atmosphere_land_fraction,
    )
    land_surface_moisture = _land_surface_field(
        atmosphere.variables.prognostic.land.soil_moisture,
        atmosphere_land_fraction,
    )
    terrarium_diagnostics = if config.land_model == :terrarium
        terrarium = atmosphere.variables.prognostic.land.terrarium
        soil_temperature = Float64.(
            Array(Oceananigans.interior(terrarium.temperature))[:, 1, :],
        )
        soil_saturation = Float64.(
            Array(Oceananigans.interior(terrarium.saturation_water_ice))[:, 1, :],
        )
        field_grid = Terrarium.get_field_grid(atmosphere.model.land.model.grid)
        soil_layer_depth = _host_vertical_center_nodes(field_grid)
        length(soil_layer_depth) == size(soil_temperature, 2) || error(
            "Terrarium soil depths and state layers differ",
        )
        land_points = findall(atmosphere_land_fraction .> 0)
        length(land_points) == size(soil_temperature, 1) || error(
            "Terrarium land columns and atmosphere land mask differ",
        )
        has_prescribed_vegetation =
            config.terrarium_evapotranspiration == :era5_prescribed_vegetation
        if has_prescribed_vegetation
            for name in (
                :vegetation_fraction,
                :leaf_area_index,
                :root_fraction,
                :root_water_availability,
                :root_uptake_weight_sum,
                :ground_water_flux,
                :transpiration_water_flux,
            )
                hasproperty(terrarium, name) || error(
                    "prescribed-vegetation Terrarium state is missing $name",
                )
            end
        end
        (
            terrarium_land_atmosphere_point = Int64.(land_points),
            terrarium_soil_layer_depth = soil_layer_depth,
            terrarium_soil_temperature = soil_temperature,
            terrarium_soil_saturation = soil_saturation,
            terrarium_skin_temperature = Float64.(vec(Array(
                Oceananigans.interior(terrarium.skin_temperature),
            ))),
            terrarium_surface_excess_water = Float64.(vec(Array(
                Oceananigans.interior(terrarium.surface_excess_water),
            ))),
            terrarium_surface_runoff = Float64.(vec(Array(
                Oceananigans.interior(terrarium.surface_runoff),
            ))),
            terrarium_infiltration = Float64.(vec(Array(
                Oceananigans.interior(terrarium.infiltration),
            ))),
            terrarium_vegetation_fraction = has_prescribed_vegetation ?
                Float64.(vec(Array(Oceananigans.interior(
                    terrarium.vegetation_fraction,
                )))) : Float64[],
            terrarium_leaf_area_index = has_prescribed_vegetation ?
                Float64.(vec(Array(Oceananigans.interior(
                    terrarium.leaf_area_index,
                )))) : Float64[],
            terrarium_root_fraction = has_prescribed_vegetation ?
                Float64.(Array(Oceananigans.interior(
                    terrarium.root_fraction,
                ))[:, 1, :]) : zeros(Float64, 0, 0),
            terrarium_root_water_availability = has_prescribed_vegetation ?
                Float64.(Array(Oceananigans.interior(
                    terrarium.root_water_availability,
                ))[:, 1, :]) : zeros(Float64, 0, 0),
            terrarium_root_uptake_weight_sum = has_prescribed_vegetation ?
                Float64.(vec(Array(Oceananigans.interior(
                    terrarium.root_uptake_weight_sum,
                )))) : Float64[],
            terrarium_ground_evaporation_flux = has_prescribed_vegetation ?
                _FRESHWATER_DENSITY_KG_M3 .* Float64.(vec(Array(
                    Oceananigans.interior(terrarium.ground_water_flux),
                ))) : Float64[],
            terrarium_transpiration_flux = has_prescribed_vegetation ?
                _FRESHWATER_DENSITY_KG_M3 .* Float64.(vec(Array(
                    Oceananigans.interior(terrarium.transpiration_water_flux),
                ))) : Float64[],
        )
    else
        (
            terrarium_land_atmosphere_point = Int64[],
            terrarium_soil_layer_depth = Float64[],
            terrarium_soil_temperature = zeros(Float64, 0, 0),
            terrarium_soil_saturation = zeros(Float64, 0, 0),
            terrarium_skin_temperature = Float64[],
            terrarium_surface_excess_water = Float64[],
            terrarium_surface_runoff = Float64[],
            terrarium_infiltration = Float64[],
            terrarium_vegetation_fraction = Float64[],
            terrarium_leaf_area_index = Float64[],
            terrarium_root_fraction = zeros(Float64, 0, 0),
            terrarium_root_water_availability = zeros(Float64, 0, 0),
            terrarium_root_uptake_weight_sum = Float64[],
            terrarium_ground_evaporation_flux = Float64[],
            terrarium_transpiration_flux = Float64[],
        )
    end
    cloud_diagnostics = _atmosphere_cloud_layer_diagnostics(atmosphere.model)
    atmosphere_radiative_surface_temperature =
        _atmosphere_radiative_surface_temperature(
            atmosphere.variables,
            atmosphere.model,
        )
    ao_fluxes = earth.interfaces.atmosphere_ocean_interface.fluxes
    ai_fluxes = earth.interfaces.atmosphere_sea_ice_interface.fluxes
    ocean_radiation = earth.radiation.interface_fluxes.ocean
    sea_ice_u = if config.sea_ice_dynamics
        _surface_grid_matrix(sea_ice.velocities.u)
    else
        zeros(size(ocean_temperature, 1), size(ocean_temperature, 2))
    end
    sea_ice_v = if config.sea_ice_dynamics
        _surface_grid_matrix(sea_ice.velocities.v)
    else
        zeros(size(ocean_temperature, 1), size(ocean_temperature, 2))
    end
    callback_temperature = atmosphere.model.callbacks[:global_surface_temperature].temperature
    radiation_budget = atmosphere.model.callbacks[:global_radiation_budget]
    atmosphere_budget = atmosphere.model.callbacks[:global_atmosphere_diagnostics]
    cloud_condensate_summary = _prognostic_cloud_condensate_summary(
        atmosphere.model.longwave_radiation,
        atmosphere_budget.point_weights,
    )
    atmosphere_process_profiles = (
        cumulative_global_convective_humidity_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_convective_humidity_change_profile,
                config.nlayers,
            ),
        cumulative_global_convective_temperature_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_convective_temperature_change_profile,
                config.nlayers,
            ),
        cumulative_global_large_scale_humidity_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_large_scale_humidity_change_profile,
                config.nlayers,
            ),
        cumulative_global_large_scale_temperature_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_large_scale_temperature_change_profile,
                config.nlayers,
            ),
        cumulative_global_radiative_temperature_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_radiative_temperature_change_profile,
                config.nlayers,
            ),
        cumulative_global_surface_humidity_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_surface_humidity_change_profile,
                config.nlayers,
            ),
        cumulative_global_surface_sensible_temperature_change_profile =
            _atmosphere_process_profile_history(
                atmosphere_budget,
                :cumulative_surface_sensible_temperature_change_profile,
                config.nlayers,
            ),
    )
    coupled_budget = simulation.callbacks[:global_coupled_budgets].func
    launch_diagnostics = balanced_launch_diagnostics(earth)
    launch_reference = if isnothing(launch_diagnostics) ||
                          !hasproperty(
                              launch_diagnostics.provenance,
                              :reference,
                          )
        nothing
    else
        launch_diagnostics.provenance.reference
    end
    ocean_profile_samples = length(coupled_budget.time_days)
    balanced_launch_source_count =
        length(_BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS)
    length(coupled_budget.balanced_launch_cumulative_sources) ==
        balanced_launch_source_count * ocean_profile_samples || error(
        "balanced-launch source history does not match coupled budget samples",
    )
    balanced_launch_cumulative_sources = reshape(
        copy(coupled_budget.balanced_launch_cumulative_sources),
        balanced_launch_source_count,
        ocean_profile_samples,
    )
    ocean_profile_layers = size(ocean_temperature, 3)
    ocean_mean_temperature_profile = _reshape_profile_history(
        coupled_budget.ocean_mean_temperature_profile,
        ocean_profile_layers,
        ocean_profile_samples,
        "ocean temperature profile",
    )
    ocean_mean_salinity_profile = _reshape_profile_history(
        coupled_budget.ocean_mean_salinity_profile,
        ocean_profile_layers,
        ocean_profile_samples,
        "ocean salinity profile",
    )
    soil_profile_samples = length(atmosphere_budget.time_days)
    soil_profile_layers = length(terrarium_diagnostics.terrarium_soil_layer_depth)
    global_land_soil_layer_temperature = _reshape_profile_history(
        atmosphere_budget.land_soil_layer_temperature,
        soil_profile_layers,
        soil_profile_samples,
        "land soil-temperature profile",
    )
    global_land_soil_layer_saturation = _reshape_profile_history(
        atmosphere_budget.land_soil_layer_saturation,
        soil_profile_layers,
        soil_profile_samples,
        "land soil-saturation profile",
    )
    ocean_reference_density = Float64(earth.interfaces.ocean_properties.reference_density)
    ocean_heat_capacity = Float64(earth.interfaces.ocean_properties.heat_capacity)
    ocean_heat_content = ocean_reference_density * ocean_heat_capacity .*
        coupled_budget.ocean_temperature_integral
    ocean_salt_content = ocean_reference_density / 1000 .*
        coupled_budget.ocean_salinity_integral
    balanced_launch_ocean_heat_source = view(
        balanced_launch_cumulative_sources,
        _balanced_launch_source_index(
            :ocean_ice,
            :ocean_sensible_heat_j_relative_to_0c,
        ),
        :,
    )
    balanced_launch_ocean_volume_source = view(
        balanced_launch_cumulative_sources,
        _balanced_launch_source_index(:ocean_ice, :ocean_volume_m3),
        :,
    )
    balanced_launch_ocean_salt_source = view(
        balanced_launch_cumulative_sources,
        _balanced_launch_source_index(:ocean_ice, :ocean_salt_mass_kg),
        :,
    )
    # NumericalEarth's ocean heat-flux convention is positive upward (out of
    # the ocean).  A mutable z-star grid also admits the enthalpy carried by
    # freshwater volume, which is deliberately stored separately from the
    # temperature boundary flux.  Both terms are required for closure.
    ocean_heat_closure_residual = ocean_heat_content .- first(ocean_heat_content) .+
        coupled_budget.cumulative_net_ocean_heat .-
        coupled_budget.cumulative_ocean_freshwater_enthalpy .-
        balanced_launch_ocean_heat_source
    ocean_heat_state_change = ocean_heat_content .- first(ocean_heat_content)
    ocean_heat_closure_relative = abs.(ocean_heat_closure_residual) ./ max.(
        max.(
            abs.(ocean_heat_state_change),
            abs.(coupled_budget.cumulative_net_ocean_heat),
        ),
        max.(
            max.(
                abs.(coupled_budget.cumulative_ocean_freshwater_enthalpy),
                abs.(balanced_launch_ocean_heat_source),
            ),
            1.0,
        ),
    )
    ocean_freshwater_closure_residual = coupled_budget.ocean_volume .-
        first(coupled_budget.ocean_volume) .-
        coupled_budget.cumulative_net_ocean_freshwater_volume .-
        balanced_launch_ocean_volume_source
    ocean_freshwater_state_change = coupled_budget.ocean_volume .-
        first(coupled_budget.ocean_volume)
    ocean_freshwater_closure_relative = abs.(ocean_freshwater_closure_residual) ./
        max.(
            max.(
                abs.(ocean_freshwater_state_change),
                abs.(coupled_budget.cumulative_net_ocean_freshwater_volume),
            ),
            max.(abs.(balanced_launch_ocean_volume_source), 1.0),
        )
    # Positive salinity boundary flux extracts salt from the ocean.
    ocean_salt_closure_residual = ocean_reference_density / 1000 .* (
        coupled_budget.ocean_salinity_integral .-
        first(coupled_budget.ocean_salinity_integral) .+
        coupled_budget.cumulative_net_ocean_salinity
    ) .- balanced_launch_ocean_salt_source
    # Freshwater carries no salt, so the relevant numerical-conservation
    # measure is drift relative to the initial salt inventory. Normalizing by
    # the near-zero physical salt boundary transport instead makes harmless
    # Float32 roundoff appear order one. This follows NumericalEarth's own
    # tracer-budget test, which checks |Δ∫S dV| / ∫S₀ dV.
    ocean_salt_closure_relative = abs.(ocean_salt_closure_residual) ./
        max(abs(first(ocean_salt_content)), 1.0)
    timestep_days = Float64(atmosphere.model.time_stepping.Δt_sec) / 86_400
    atmosphere_temperature = Float64.(
        Array(atmosphere.variables.grid.temperature_prev.data),
    )
    atmosphere_humidity = Float64.(
        Array(atmosphere.variables.grid.humidity_prev.data),
    )
    atmosphere_convective_humidity_tendency = Float64.(Array(
        atmosphere.variables.parameterizations.convective_humidity_tendency.data,
    ))
    atmosphere_convective_temperature_tendency = Float64.(Array(
        atmosphere.variables.parameterizations.convective_temperature_tendency.data,
    ))
    atmosphere_large_scale_condensation_humidity_tendency = Float64.(Array(
        atmosphere.variables.parameterizations.large_scale_condensation_humidity_tendency.data,
    ))
    atmosphere_large_scale_condensation_temperature_tendency = Float64.(Array(
        atmosphere.variables.parameterizations.large_scale_condensation_temperature_tendency.data,
    ))
    atmosphere_u = Float64.(Array(atmosphere.variables.grid.u_prev.data))
    atmosphere_v = Float64.(Array(atmosphere.variables.grid.v_prev.data))
    atmosphere_surface_pressure = Float64.(
        Array(atmosphere.variables.grid.pressure_prev.data),
    )
    atmosphere_sigma = Float64.(Array(atmosphere.model.geometry.σ_levels_full))
    atmosphere_layer_thickness = Float64.(
        Array(atmosphere.model.geometry.σ_levels_thick),
    )
    atmosphere_pressure = atmosphere_surface_pressure .* reshape(atmosphere_sigma, 1, :)
    atmosphere_column_water_vapor = atmosphere_surface_pressure ./
        Float64(atmosphere.model.planet.gravity) .*
        vec(sum(
            atmosphere_humidity .* reshape(atmosphere_layer_thickness, 1, :);
            dims = 2,
        ))
    atmosphere_rainfall_flux = _FRESHWATER_DENSITY_KG_M3 .* Float64.(
        Array(atmosphere.variables.parameterizations.rain_rate.data),
    )
    atmosphere_convective_rainfall_flux =
        _FRESHWATER_DENSITY_KG_M3 .* Float64.(Array(
            atmosphere.variables.parameterizations.rain_rate_convection.data,
        ))
    atmosphere_convective_precipitation_mass_correction_flux =
        _FRESHWATER_DENSITY_KG_M3 .* Float64.(Array(
            atmosphere.variables.parameterizations.convective_precipitation_mass_correction.data,
        ))
    atmosphere_large_scale_rainfall_flux =
        _FRESHWATER_DENSITY_KG_M3 .* Float64.(Array(
            atmosphere.variables.parameterizations.rain_rate_large_scale.data,
        ))
    atmosphere_large_scale_snowfall_flux =
        _FRESHWATER_DENSITY_KG_M3 .* Float64.(Array(
            atmosphere.variables.parameterizations.snow_rate_large_scale.data,
        ))
    atmosphere_large_scale_precipitation_mass_correction_flux =
        _FRESHWATER_DENSITY_KG_M3 .* Float64.(Array(
            atmosphere.variables.parameterizations.large_scale_precipitation_mass_correction.data,
        ))
    atmosphere_snowfall_flux = _FRESHWATER_DENSITY_KG_M3 .* Float64.(
        Array(atmosphere.variables.parameterizations.snow_rate.data),
    )
    atmosphere_precipitation_cloud_top_layer = Float64.(Array(
        atmosphere.variables.parameterizations.cloud_top.data,
    ))
    atmosphere_vertical_diffusion_cfl_scale =
        haskey(
            atmosphere.variables.parameterizations,
            :vertical_diffusion_cfl_scale,
        ) ? Float64.(Array(
            atmosphere.variables.parameterizations.vertical_diffusion_cfl_scale.data,
        )) : ones(Float64, length(atmosphere_longitude))

    return (
        ocean_longitude = longitude,
        ocean_latitude = latitude,
        ocean_layer_depth,
        ocean_surface_active_mask,
        ocean_cell_area,
        ocean_runoff_receiver_mask = runoff_receiver_footprint.receiver_mask,
        ocean_routed_land_runoff_flux,
        ocean_river_mouth_receiver_mask = river_mouth_footprint.receiver_mask,
        ocean_river_mouth_horizontal_mixing_mask =
            river_mouth_footprint.horizontal_mask,
        ocean_river_mouth_active_mixing_mask =
            river_mouth_active_mixing_mask,
        ocean_surface_temperature = ocean_temperature[:, :, end],
        ocean_surface_salinity = ocean_salinity[:, :, end],
        ocean_temperature_minimum = first(ocean_temperature_extrema),
        ocean_temperature_maximum = last(ocean_temperature_extrema),
        ocean_salinity_minimum = first(ocean_salinity_extrema),
        ocean_salinity_maximum = last(ocean_salinity_extrema),
        ocean_mean_temperature_profile,
        ocean_mean_salinity_profile,
        ocean_surface_zonal_velocity = ocean_u[
            1:size(ocean_temperature, 1),
            1:size(ocean_temperature, 2),
            end,
        ],
        ocean_surface_meridional_velocity = ocean_v[
            1:size(ocean_temperature, 1),
            1:size(ocean_temperature, 2),
            end,
        ],
        ocean_max_abs_zonal_velocity = Float64(maximum(abs, ocean.velocities.u)),
        ocean_max_abs_meridional_velocity = Float64(maximum(abs, ocean.velocities.v)),
        ocean_max_abs_vertical_velocity = Float64(maximum(abs, ocean.velocities.w)),
        ocean_mean_kinetic_energy = 0.5 * (
            _active_ocean_mean_square(ocean.velocities.u) +
            _active_ocean_mean_square(ocean.velocities.v) +
            _active_ocean_mean_square(ocean.velocities.w)
        ),
        sea_ice_concentration = _interior_matrix(sea_ice.ice_concentration),
        sea_ice_thickness = _interior_matrix(sea_ice.ice_thickness),
        sea_ice_zonal_velocity = sea_ice_u,
        sea_ice_meridional_velocity = sea_ice_v,
        atmosphere_ocean_sensible_heat_flux = _interior_matrix(ao_fluxes.sensible_heat),
        atmosphere_ocean_latent_heat_flux = _interior_matrix(ao_fluxes.latent_heat),
        atmosphere_sea_ice_sensible_heat_flux = _interior_matrix(ai_fluxes.sensible_heat),
        atmosphere_sea_ice_latent_heat_flux = _interior_matrix(ai_fluxes.latent_heat),
        ocean_upwelling_longwave = _interior_matrix(ocean_radiation.upwelling_longwave),
        ocean_absorbed_longwave = _interior_matrix(ocean_radiation.downwelling_longwave),
        ocean_absorbed_shortwave = _interior_matrix(ocean_radiation.downwelling_shortwave),
        atmosphere_longitude = Float64.(atmosphere_longitude),
        atmosphere_latitude = Float64.(atmosphere_latitude),
        atmosphere_land_fraction,
        land_surface_temperature,
        land_surface_moisture,
        terrarium_diagnostics...,
        atmosphere_surface_temperature = Float64.(
            Array(atmosphere.variables.parameterizations.surface_air_temperature.data),
        ),
        atmosphere_radiative_surface_temperature,
        atmosphere_sigma,
        atmosphere_sigma_layer_thickness = atmosphere_layer_thickness,
        atmosphere_temperature,
        atmosphere_specific_humidity = atmosphere_humidity,
        atmosphere_convective_humidity_tendency,
        atmosphere_convective_temperature_tendency,
        atmosphere_large_scale_condensation_humidity_tendency,
        atmosphere_large_scale_condensation_temperature_tendency,
        atmosphere_pressure,
        atmosphere_zonal_wind = atmosphere_u,
        atmosphere_meridional_wind = atmosphere_v,
        atmosphere_cloud_fraction = cloud_diagnostics.fraction,
        atmosphere_cloud_liquid_water_path =
            cloud_diagnostics.liquid_water_path,
        atmosphere_cloud_ice_water_path = cloud_diagnostics.ice_water_path,
        atmosphere_surface_pressure,
        atmosphere_surface_specific_humidity = atmosphere_humidity[:, end],
        atmosphere_surface_zonal_wind = atmosphere_u[:, end],
        atmosphere_surface_meridional_wind = atmosphere_v[:, end],
        atmosphere_column_water_vapor,
        atmosphere_rainfall_flux,
        atmosphere_convective_rainfall_flux,
        atmosphere_convective_precipitation_mass_correction_flux,
        atmosphere_large_scale_rainfall_flux,
        atmosphere_large_scale_snowfall_flux,
        atmosphere_large_scale_precipitation_mass_correction_flux,
        atmosphere_snowfall_flux,
        atmosphere_precipitation_cloud_top_layer,
        atmosphere_vertical_diffusion_cfl_scale,
        atmosphere_surface_water_vapor_flux = Float64.(
            Array(atmosphere.variables.parameterizations.surface_humidity_flux.data),
        ),
        atmosphere_surface_sensible_heat_flux = Float64.(
            Array(atmosphere.variables.parameterizations.sensible_heat_flux.data),
        ),
        outgoing_longwave = Float64.(
            Array(atmosphere.variables.parameterizations.outgoing_longwave.data),
        ),
        column_cloud_fraction = cloud_diagnostics.column_fraction,
        time_days = collect((0:(length(callback_temperature) - 1)) .* timestep_days),
        global_surface_temperature = Float64.(callback_temperature),
        toa_incoming_shortwave = Float64.(radiation_budget.incoming_shortwave),
        toa_outgoing_shortwave = Float64.(radiation_budget.outgoing_shortwave),
        toa_outgoing_longwave = Float64.(radiation_budget.outgoing_longwave),
        toa_clear_outgoing_shortwave = Float64.(
            radiation_budget.clear_outgoing_shortwave,
        ),
        toa_clear_outgoing_longwave = Float64.(
            radiation_budget.clear_outgoing_longwave,
        ),
        toa_net_downward = Float64.(radiation_budget.net_downward),
        toa_clear_net_downward = Float64.(radiation_budget.clear_net_downward),
        atmosphere_diagnostic_time_days = Float64.(atmosphere_budget.time_days),
        atmosphere_process_profiles...,
        global_rainfall_flux = Float64.(atmosphere_budget.rainfall_flux),
        global_convective_rainfall_flux = Float64.(
            atmosphere_budget.convective_rainfall_flux,
        ),
        global_large_scale_rainfall_flux = Float64.(
            atmosphere_budget.large_scale_rainfall_flux,
        ),
        global_large_scale_snowfall_flux = Float64.(
            atmosphere_budget.large_scale_snowfall_flux,
        ),
        global_snowfall_flux = Float64.(atmosphere_budget.snowfall_flux),
        global_surface_water_vapor_flux = Float64.(
            atmosphere_budget.surface_water_vapor_flux,
        ),
        cumulative_global_surface_water_vapor = Float64.(
            atmosphere_budget.cumulative_surface_water_vapor,
        ),
        cumulative_global_precipitation = Float64.(
            atmosphere_budget.cumulative_precipitation,
        ),
        cumulative_balanced_launch_atmosphere_water_source = Float64.(
            atmosphere_budget.cumulative_balanced_launch_atmosphere_water_source,
        ),
        global_atmosphere_water_budget_residual = Float64.(
            atmosphere_budget.water_budget_residual,
        ),
        global_cloud_condensate_water_path = Float64.(
            atmosphere_budget.cloud_condensate_water_path,
        ),
        global_total_atmosphere_water_budget_residual = Float64.(
            atmosphere_budget.total_water_budget_residual,
        ),
        global_surface_air_temperature = Float64.(
            atmosphere_budget.surface_air_temperature,
        ),
        global_surface_specific_humidity = Float64.(
            atmosphere_budget.surface_specific_humidity,
        ),
        global_surface_pressure = Float64.(atmosphere_budget.surface_pressure),
        global_near_surface_wind_speed = Float64.(
            atmosphere_budget.near_surface_wind_speed,
        ),
        global_column_water_vapor = Float64.(atmosphere_budget.column_water_vapor),
        global_mass_weighted_atmosphere_temperature = Float64.(
            atmosphere_budget.mass_weighted_temperature,
        ),
        global_land_surface_soil_moisture = Float64.(
            atmosphere_budget.land_surface_soil_moisture,
        ),
        global_land_total_water_storage = Float64.(
            atmosphere_budget.land_total_water_storage,
        ),
        cumulative_global_land_precipitation = Float64.(
            atmosphere_budget.cumulative_land_precipitation,
        ),
        cumulative_global_land_evapotranspiration = Float64.(
            atmosphere_budget.cumulative_land_evapotranspiration,
        ),
        cumulative_global_land_surface_runoff = Float64.(
            atmosphere_budget.cumulative_land_surface_runoff,
        ),
        cumulative_balanced_launch_land_water_source = Float64.(
            atmosphere_budget.cumulative_balanced_launch_land_water_source,
        ),
        global_land_water_budget_residual = Float64.(
            atmosphere_budget.land_water_budget_residual,
        ),
        global_land_soil_layer_temperature,
        global_land_soil_layer_saturation,
        global_land_rainfall_flux = Float64.(atmosphere_budget.land_rainfall_flux),
        global_land_snowfall_flux = Float64.(atmosphere_budget.land_snowfall_flux),
        global_land_evaporation_flux = Float64.(
            atmosphere_budget.land_evaporation_flux,
        ),
        global_land_ground_evaporation_flux = Float64.(
            atmosphere_budget.land_ground_evaporation_flux,
        ),
        global_land_transpiration_flux = Float64.(
            atmosphere_budget.land_transpiration_flux,
        ),
        global_land_surface_runoff_flux = Float64.(
            atmosphere_budget.land_surface_runoff_flux,
        ),
        global_land_infiltration_flux = Float64.(
            atmosphere_budget.land_infiltration_flux,
        ),
        coupled_budget_time_days = Float64.(coupled_budget.time_days),
        balanced_launch_applied_fraction = Float64.(
            coupled_budget.balanced_launch_applied_fraction,
        ),
        balanced_launch_cumulative_sources = Float64.(
            balanced_launch_cumulative_sources,
        ),
        ocean_heat_content = Float64.(ocean_heat_content),
        ocean_salt_content = Float64.(ocean_salt_content),
        ocean_volume = Float64.(coupled_budget.ocean_volume),
        ocean_temperature_minimum_series = Float64.(
            coupled_budget.ocean_temperature_minimum_series,
        ),
        ocean_temperature_maximum_series = Float64.(
            coupled_budget.ocean_temperature_maximum_series,
        ),
        ocean_salinity_minimum_series = Float64.(
            coupled_budget.ocean_salinity_minimum_series,
        ),
        ocean_salinity_maximum_series = Float64.(
            coupled_budget.ocean_salinity_maximum_series,
        ),
        ocean_surface_salinity_minimum_series = Float64.(
            coupled_budget.ocean_surface_salinity_minimum_series,
        ),
        ocean_surface_salinity_minimum_longitude = Float64.(
            coupled_budget.ocean_surface_salinity_minimum_longitude,
        ),
        ocean_surface_salinity_minimum_latitude = Float64.(
            coupled_budget.ocean_surface_salinity_minimum_latitude,
        ),
        ocean_surface_salinity_minimum_i = Int64.(
            coupled_budget.ocean_surface_salinity_minimum_i,
        ),
        ocean_surface_salinity_minimum_j = Int64.(
            coupled_budget.ocean_surface_salinity_minimum_j,
        ),
        ocean_runoff_receiver_surface_salinity_minimum_series = Float64.(
            coupled_budget.ocean_runoff_receiver_surface_salinity_minimum_series,
        ),
        ocean_runoff_receiver_surface_salinity_minimum_longitude = Float64.(
            coupled_budget.ocean_runoff_receiver_surface_salinity_minimum_longitude,
        ),
        ocean_runoff_receiver_surface_salinity_minimum_latitude = Float64.(
            coupled_budget.ocean_runoff_receiver_surface_salinity_minimum_latitude,
        ),
        ocean_runoff_receiver_surface_salinity_minimum_i = Int64.(
            coupled_budget.ocean_runoff_receiver_surface_salinity_minimum_i,
        ),
        ocean_runoff_receiver_surface_salinity_minimum_j = Int64.(
            coupled_budget.ocean_runoff_receiver_surface_salinity_minimum_j,
        ),
        ocean_river_mouth_active_mixing_cells = Int64.(
            coupled_budget.ocean_river_mouth_active_mixing_cells,
        ),
        ocean_river_mouth_active_mixing_area_m2 = Float64.(
            coupled_budget.ocean_river_mouth_active_mixing_area_m2,
        ),
        sea_ice_volume = Float64.(coupled_budget.sea_ice_volume),
        sea_ice_area = Float64.(coupled_budget.sea_ice_area),
        tropical_sea_ice_area = Float64.(coupled_budget.tropical_sea_ice_area),
        sea_ice_max_thickness = Float64.(coupled_budget.sea_ice_max_thickness),
        sea_ice_extent_max_thickness = Float64.(
            coupled_budget.sea_ice_extent_max_thickness,
        ),
        sea_ice_max_area_equivalent_thickness = Float64.(
            coupled_budget.sea_ice_max_area_equivalent_thickness,
        ),
        sea_ice_marginal_ceiling_cell_count = Int64.(
            coupled_budget.sea_ice_marginal_ceiling_cell_count,
        ),
        sea_ice_marginal_ceiling_max_concentration = Float64.(
            coupled_budget.sea_ice_marginal_ceiling_max_concentration,
        ),
        sea_ice_marginal_ceiling_max_area_equivalent_thickness = Float64.(
            coupled_budget.sea_ice_marginal_ceiling_max_area_equivalent_thickness,
        ),
        ocean_free_surface_minimum = Float64.(
            coupled_budget.ocean_free_surface_minimum,
        ),
        ocean_free_surface_maximum = Float64.(
            coupled_budget.ocean_free_surface_maximum,
        ),
        ocean_max_abs_zonal_velocity_series = Float64.(
            coupled_budget.ocean_max_abs_zonal_velocity,
        ),
        ocean_max_abs_meridional_velocity_series = Float64.(
            coupled_budget.ocean_max_abs_meridional_velocity,
        ),
        cumulative_ocean_surface_heat = Float64.(
            coupled_budget.cumulative_ocean_surface_heat,
        ),
        cumulative_ocean_frazil_heat = Float64.(
            coupled_budget.cumulative_ocean_frazil_heat,
        ),
        cumulative_ocean_penetrating_shortwave_heat = Float64.(
            coupled_budget.cumulative_ocean_penetrating_shortwave_heat,
        ),
        cumulative_net_ocean_heat = Float64.(coupled_budget.cumulative_net_ocean_heat),
        cumulative_ocean_freshwater_enthalpy = Float64.(
            coupled_budget.cumulative_ocean_freshwater_enthalpy,
        ),
        cumulative_net_ocean_freshwater_volume = Float64.(
            coupled_budget.cumulative_net_ocean_freshwater_volume,
        ),
        cumulative_net_ocean_salinity = Float64.(
            coupled_budget.cumulative_net_ocean_salinity,
        ),
        cumulative_rain_mass = Float64.(coupled_budget.cumulative_rain_mass),
        cumulative_snow_mass = Float64.(coupled_budget.cumulative_snow_mass),
        cumulative_open_ocean_evaporation_mass = Float64.(
            coupled_budget.cumulative_open_ocean_evaporation_mass,
        ),
        cumulative_sea_ice_ocean_freshwater_volume = Float64.(
            coupled_budget.cumulative_sea_ice_ocean_freshwater_volume,
        ),
        cumulative_routed_land_runoff_mass = Float64.(
            coupled_budget.cumulative_routed_land_runoff_mass,
        ),
        ocean_heat_closure_residual = Float64.(ocean_heat_closure_residual),
        ocean_heat_closure_relative = Float64.(ocean_heat_closure_relative),
        ocean_freshwater_closure_residual = Float64.(
            ocean_freshwater_closure_residual,
        ),
        ocean_freshwater_closure_relative = Float64.(
            ocean_freshwater_closure_relative,
        ),
        ocean_salt_closure_residual = Float64.(ocean_salt_closure_residual),
        ocean_salt_closure_relative = Float64.(ocean_salt_closure_relative),
        metadata = (
            experiment = config.name,
            balanced_launch = isnothing(launch_diagnostics) ?
                "none" : launch_diagnostics.provenance.schema,
            balanced_launch_start_time_seconds =
                isnothing(launch_diagnostics) ? NaN :
                launch_diagnostics.start_time_seconds,
            balanced_launch_window_seconds =
                isnothing(launch_diagnostics) ? NaN :
                launch_diagnostics.window_seconds,
            balanced_launch_final_applied_fraction =
                isnothing(launch_diagnostics) ? 0.0 :
                launch_diagnostics.applied_fraction,
            balanced_launch_target_name = isnothing(launch_diagnostics) ?
                "none" : launch_diagnostics.provenance.target_name,
            balanced_launch_target_start_date =
                isnothing(launch_diagnostics) ? "none" :
                launch_diagnostics.provenance.target_start_date,
            balanced_launch_distribution = isnothing(launch_diagnostics) ?
                "none" : launch_diagnostics.provenance.distribution,
            balanced_launch_source_accounting =
                isnothing(launch_diagnostics) ? "none" :
                launch_diagnostics.provenance.source_accounting,
            balanced_launch_provenance_sha256 =
                isnothing(launch_diagnostics) ? "none" : bytes2hex(SHA.sha256(
                    codeunits(repr(_balanced_launch_provenance_fingerprint(
                        launch_diagnostics.provenance,
                    ))),
                )),
            balanced_launch_reference_mode = isnothing(launch_reference) ?
                "none" : String(launch_reference.mode),
            balanced_launch_reference_background_name =
                isnothing(launch_reference) ? "none" :
                String(launch_reference.background_name),
            balanced_launch_reference_background_start_date =
                isnothing(launch_reference) ? "none" :
                String(launch_reference.background_start_date),
            balanced_launch_reference_forecast_end_time_seconds =
                isnothing(launch_reference) ? NaN :
                Float64(launch_reference.forecast_end_time_seconds),
            balanced_launch_reference_forecast_end_iteration =
                isnothing(launch_reference) ? -1 :
                Int64(launch_reference.forecast_end_iteration),
            balanced_launch_reference_increment_sha256 =
                isnothing(launch_reference) ? "none" :
                String(launch_reference.increment_sha256),
            device = String(config.device),
            truncation = config.truncation,
            nlayers = config.nlayers,
            atmosphere = "SpeedyWeather",
            atmosphere_humidity_nonnegative_scheme =
                "Gaussian-area layer-mean conservative rescaling",
            atmosphere_convective_precipitation_scheme =
                "Betts-Miller net column humidity drying",
            atmosphere_convection = String(config.atmosphere_convection),
            atmosphere_convection_reference_humidity_profile =
                config.atmosphere_convection == :betts_miller_sigma_rh_v1 ?
                "piecewise_linear_upper_to_lower_physical_sigma_v1" :
                "constant_relative_humidity_0.7",
            atmosphere_convection_upper_relative_humidity =
                config.atmosphere_convection == :betts_miller_sigma_rh_v1 ?
                config.atmosphere_convection_upper_relative_humidity : 0.7,
            atmosphere_convection_lower_relative_humidity =
                config.atmosphere_convection == :betts_miller_sigma_rh_v1 ?
                config.atmosphere_convection_lower_relative_humidity : 0.7,
            atmosphere_convection_transition_top_sigma =
                config.atmosphere_convection == :betts_miller_sigma_rh_v1 ?
                config.atmosphere_convection_transition_top_sigma : NaN,
            atmosphere_convection_transition_bottom_sigma =
                config.atmosphere_convection == :betts_miller_sigma_rh_v1 ?
                config.atmosphere_convection_transition_bottom_sigma : NaN,
            atmosphere_large_scale_precipitation_scheme =
                config.atmosphere_large_scale_precipitation ==
                    :upstream_observed ?
                "SpeedyWeather implicit condensation observed without surface-rate correction" :
                "implicit condensation net column humidity drying",
            atmosphere_layer_process_diagnostics =
                "signed_incremental_parameterization_tendencies_v1",
            atmosphere_cumulative_process_profile_diagnostics =
                "Gaussian-area-weighted applied-process integral v1",
            atmosphere_process_profile_start_day =
                atmosphere_budget.process_diagnostic_start_timestep *
                Float64(atmosphere.model.time_stepping.Δt_sec) / 86_400,
            atmosphere_precipitation_mass_correction_diagnostics =
                config.atmosphere_large_scale_precipitation ==
                    :upstream_observed ?
                "convective_net_column_correction_and_zero_large_scale_control_v1" :
                "net_column_minus_upstream_surface_rate_v1",
            atmosphere_cloud_condensate_state =
                config.forcing.cloud_scheme == :prognostic_condensate ?
                "column_local_prognostic_liquid_ice_path_conservative_precipitation_delay_v1" :
                "not_prognostic_diagnostic_cloud_optics_only",
            cloud_condensate_summary...,
            speedy_cuda_graphs = Int(atmosphere.model.spectral_transform.cuda_graphs),
            speedy_cuda_graph_cache_key = atmosphere.model.spectral_transform.cuda_graphs ?
                "device_pointer_upstream_54d500f" : "disabled",
            coupled_gpu_component_handoff = config.device == :gpu &&
                config.land_model == :terrarium ?
                COUPLED_GPU_COMPONENT_HANDOFF_PROVENANCE : "not_applicable",
            speedy_transforms_version = string(
                pkgversion(SpeedyWeather.SpeedyTransforms),
            ),
            speedy_transforms_source = "vendor/SpeedyTransforms",
            atmosphere_to_exchange_sparse_storage = string(nameof(typeof(
                earth.interfaces.exchanger.atmosphere.regridder.from_atmosphere.intersections,
            ))),
            exchange_to_atmosphere_sparse_storage = string(nameof(typeof(
                earth.interfaces.exchanger.atmosphere.regridder.to_atmosphere.intersections,
            ))),
            atmosphere_initial_conditions = String(config.atmosphere_initial_conditions),
            era5_month = config.era5_month,
            _atmosphere_initial_condition_provenance(
                atmosphere.model.initial_conditions,
                config,
            )...,
            start_date = string(config.start_date),
            analysis_start_day = config.analysis_start_day,
            calendar = "proleptic_gregorian",
            radiation = String(config.forcing.radiation),
            atmosphere_radiative_surface_temperature_source =
                atmosphere.model.longwave_radiation isa RRTMGPRadiation ?
                "retained_rrtmgp_solver_state" :
                "diagnosed_fourth_power_surface_mixture",
            co2_ppm = config.forcing.co2_ppm,
            cloud_scheme = String(config.forcing.cloud_scheme),
            _forcing_provenance(atmosphere.model.longwave_radiation, config)...,
            atmosphere_timestep_seconds = atmosphere.model.time_stepping.Δt_sec,
            atmosphere_gravity_ms2 = atmosphere.model.planet.gravity,
            atmosphere_hyperdiffusion_hours = config.atmosphere_hyperdiffusion_hours,
            atmosphere_divergence_hyperdiffusion_hours =
                config.atmosphere_divergence_hyperdiffusion_hours,
            atmosphere_vertical_diffusion =
                _atmosphere_vertical_diffusion_provenance(
                    atmosphere.model.vertical_diffusion,
                ),
            atmosphere_speed_limit_ms = config.atmosphere_speed_limit_ms,
            atmosphere_speed_limit_drag_per_m =
                config.atmosphere_speed_limit_drag_per_m,
            atmosphere_surface_speed_limit_ms =
                config.atmosphere_surface_speed_limit_ms,
            atmosphere_surface_speed_limit_drag_per_m =
                config.atmosphere_surface_speed_limit_drag_per_m,
            land = config.land_model == :terrarium ? "Terrarium" : "SpeedyWeather bucket land",
            coupled_land_radiative_boundary = config.land_model == :terrarium ?
                "shared_constant_albedo_emissivity_v1" : "not_applicable",
            coupled_land_surface_albedo = config.land_model == :terrarium ?
                COUPLED_LAND_SURFACE_ALBEDO : NaN,
            coupled_land_surface_emissivity = config.land_model == :terrarium ?
                COUPLED_SURFACE_EMISSIVITY : NaN,
            terrarium_snowfall_hydrology = config.land_model == :terrarium ?
                "immediate_liquid_water_equivalent_no_snowpack_v1" :
                "not_applicable",
            _land_initial_condition_provenance(config)...,
            land_runoff_routing = earth.land isa TerrariumRunoffLand ?
                "16-nearest-wet-cell conservative spreading" : "none",
            land_runoff_routing_weighting = earth.land isa TerrariumRunoffLand ?
                String(earth.land.routing_weighting) : "none",
            routed_land_columns = earth.land isa TerrariumRunoffLand ?
                earth.land.source_columns : 0,
            routed_ocean_cells = earth.land isa TerrariumRunoffLand ?
                earth.land.target_cells : 0,
            ocean_runoff_receiver_cells = Int(count(
                !iszero,
                runoff_receiver_footprint.receiver_mask,
            )),
            runoff_receivers_per_land_column = earth.land isa TerrariumRunoffLand ?
                earth.land.receivers_per_source : 0,
            ocean = "Oceananigans",
            ocean_grid = String(config.ocean_grid),
            ocean_bathymetry_correction = String(
                config.ocean_bathymetry_correction,
            ),
            ocean_tripolar_wet_mask = String(config.ocean_tripolar_wet_mask),
            ocean_tripolar_wet_mask_scheme =
                config.ocean_tripolar_wet_mask == :etopo_majority_area_polar ?
                _TRIPOLAR_MAJORITY_WET_SCHEME : "none",
            ocean_tripolar_wet_mask_minimum_absolute_latitude_degrees =
                config.ocean_tripolar_wet_mask == :etopo_majority_area_polar ?
                _TRIPOLAR_MAJORITY_WET_MINIMUM_ABSOLUTE_LATITUDE : NaN,
            ocean_tripolar_wet_mask_fraction_threshold =
                config.ocean_tripolar_wet_mask == :etopo_majority_area_polar ?
                _TRIPOLAR_MAJORITY_WET_FRACTION_THRESHOLD : NaN,
            ocean_tripolar_wet_mask_minimum_etopo_samples =
                config.ocean_tripolar_wet_mask == :etopo_majority_area_polar ?
                _TRIPOLAR_MAJORITY_WET_MINIMUM_ETOPO_SAMPLES : 0,
            ocean_shortwave_scheme = String(config.ocean_shortwave_scheme),
            ocean_shortwave_optical_water_types =
                config.ocean_shortwave_scheme == :depth_aware_coastal ?
                "Jerlov_Type_III_at_depth_le_50m_linear_blend_to_Jerlov_Type_I_at_depth_ge_200m" :
                "Jerlov_Type_I",
            ocean_shortwave_type_i_parameters =
                "fraction=0.58 first_scale_m=0.35 second_scale_m=23",
            ocean_shortwave_type_iii_parameters =
                config.ocean_shortwave_scheme == :depth_aware_coastal ?
                "fraction=0.78 first_scale_m=1.4 second_scale_m=7.9" :
                "not_applied",
            ocean_shortwave_bottom_residual_treatment =
                "conservatively_absorbed_in_lowest_wet_cell",
            ocean_river_mouth_mixing = String(
                config.ocean_river_mouth_mixing,
            ),
            ocean_river_mouth_vertical_diffusivity_m2s =
                _is_localized_estuary_mixing(config.ocean_river_mouth_mixing) ?
                config.ocean_river_mouth_vertical_diffusivity_m2s : 0.0,
            ocean_river_mouth_horizontal_diffusivity_m2s =
                _is_localized_estuary_mixing(config.ocean_river_mouth_mixing) ?
                config.ocean_river_mouth_horizontal_diffusivity_m2s : 0.0,
            ocean_river_mouth_mixing_depth_m =
                _is_localized_estuary_mixing(config.ocean_river_mouth_mixing) ?
                config.ocean_river_mouth_mixing_depth_m : 0.0,
            ocean_river_mouth_reference_freshwater_mass_flux_kgm2s =
                config.ocean_river_mouth_mixing == :localized_estuary_v5 ?
                config.ocean_river_mouth_reference_freshwater_mass_flux_kgm2s : 0.0,
            ocean_river_mouth_mixing_footprint =
                _river_mouth_mixing_footprint_provenance(
                    config.ocean_river_mouth_mixing,
                ),
            ocean_river_mouth_mixing_dynamic_gating = if config.ocean_river_mouth_mixing ==
                :localized_estuary_v4
                "positive_routed_freshwater_flux_each_coupled_step"
            elseif config.ocean_river_mouth_mixing == :localized_estuary_v5
                "smooth_clamp_flux_over_reference_minus_0.5_each_coupled_step"
            else
                "none"
            end,
            ocean_river_mouth_dynamic_update_completion =
                config.ocean_river_mouth_mixing in (
                    :localized_estuary_v4,
                    :localized_estuary_v5,
                ) ? DYNAMIC_RIVER_MOUTH_UPDATE_COMPLETION_PROVENANCE :
                "not_applicable",
            ocean_river_mouth_receiver_cells = Int(count(
                !iszero,
                river_mouth_footprint.receiver_mask,
            )),
            ocean_river_mouth_horizontal_footprint_cells = Int(count(
                !iszero,
                river_mouth_footprint.horizontal_mask,
            )),
            ocean_river_mouth_active_mixing_cells = Int(count(
                !iszero,
                river_mouth_active_mixing_mask,
            )),
            ocean_nlongitude = config.ocean_nlongitude,
            ocean_nlatitude = config.ocean_nlatitude,
            ocean_nlayers = config.ocean_nlayers,
            ocean_barotropic_substeps = config.ocean_grid == :tripolar_1degree ?
                _tripolar_barotropic_substeps(config.ocean_nlatitude) : 0,
            ocean_initial_conditions = String(config.ocean_initial_conditions),
            ocean_initial_condition_fields =
                _is_full_state_ecco_initial_conditions(
                    config.ocean_initial_conditions,
                ) ?
                "ECCO4Monthly_$(_ecco_initial_condition_month_label(config))_temperature_salinity_u_velocity_v_velocity_free_surface_sea_ice_thickness_sea_ice_concentration" :
                _is_ecco_initial_conditions(config.ocean_initial_conditions) ?
                "ECCO4Monthly_$(_ecco_initial_condition_month_label(config))_temperature_salinity_sea_ice_thickness_sea_ice_concentration" :
                "analytic_temperature_salinity_sea_ice",
            ocean_initial_condition_source_month =
                _is_ecco_initial_conditions(config.ocean_initial_conditions) ?
                _ecco_initial_condition_month_label(config) : "not_applicable",
            ocean_initial_condition_directory =
                _is_ecco_initial_conditions(config.ocean_initial_conditions) ?
                (isempty(config.ecco_initial_conditions_directory) ?
                 "NumericalEarth_default_cache_or_artifact" :
                 config.ecco_initial_conditions_directory) : "not_applicable",
            ocean_initial_velocity_coordinates =
                _is_full_state_ecco_initial_conditions(
                    config.ocean_initial_conditions,
                ) ?
                "ECCO4_EVEL_NVEL_geographic_east_north_rotated_to_native_grid_axes" :
                "not_initialized_from_observed_velocity",
            sea_ice_initial_thickness_semantics =
                _is_ecco_initial_conditions(config.ocean_initial_conditions) ?
                "ECCO4_SIheff_area_equivalent_converted_to_ClimaSeaIce_conditional_h_preserving_Ah" :
                "analytic_conditional_h",
            ocean_dynamics = Int(config.ocean_dynamics),
            ocean_mixing_closure = config.ocean_grid == :latitude_longitude_1degree ?
                "CATKE_GentMcWilliams_isopycnal_biharmonic" :
                config.ocean_grid == :tripolar_1degree ?
                _ocean_mixing_closure_provenance(
                    config.ocean_river_mouth_mixing,
                ) : "idealized",
            ocean_polar_momentum_sponge =
                config.ocean_grid == :latitude_longitude_1degree ?
                "cosine_ramp_ocean_and_sea_ice" : "none",
            ocean_polar_momentum_sponge_start_latitude_degrees =
                config.ocean_grid == :latitude_longitude_1degree ?
                config.ocean_polar_sponge_start_latitude_degrees : NaN,
            ocean_polar_momentum_sponge_stop_latitude_degrees =
                config.ocean_grid == :latitude_longitude_1degree ?
                config.ocean_polar_sponge_stop_latitude_degrees : NaN,
            ocean_polar_momentum_sponge_timescale_hours =
                config.ocean_grid == :latitude_longitude_1degree ?
                config.ocean_polar_sponge_timescale_hours : NaN,
            sea_ice = "ClimaSeaIce",
            sea_ice_dynamics = Int(config.sea_ice_dynamics),
            sea_ice_momentum_substeps = config.sea_ice_dynamics ?
                _sea_ice_momentum_substeps(sea_ice) : 0,
            sea_ice_pressure_formulation = config.sea_ice_dynamics ?
                String(config.sea_ice_pressure_formulation) : "none",
            sea_ice_immersed_boundary_drag_coefficient =
                config.sea_ice_dynamics ?
                config.sea_ice_immersed_boundary_drag_coefficient : NaN,
            sea_ice_inactive_surface_mask =
                _SEA_ICE_INACTIVE_SURFACE_MASK_SCHEME,
            sea_ice_inactive_surface_mask_value = 0.0,
            sea_ice_advection = config.sea_ice_dynamics ?
                _CONSERVATIVE_SEA_ICE_ADVECTION_PROVENANCE : "none",
            sea_ice_maximum_conditional_thickness_m =
                config.sea_ice_dynamics ? Float64(
                    sea_ice.advection.maximum_conditional_thickness,
                ) : NaN,
            sea_ice_extent_concentration_threshold =
                _SEA_ICE_EXTENT_CONCENTRATION_THRESHOLD,
            tropical_sea_ice_maximum_absolute_latitude_degrees =
                _TROPICAL_SEA_ICE_MAXIMUM_ABSOLUTE_LATITUDE_DEGREES,
            sea_ice_raw_thickness_is_representation_diagnostic = 1,
            sea_ice_physical_acceptance_metric =
                "maximum conditional h within A>=0.15 and maximum A*h",
            surface_radiation_coupling = "live SpeedyWeather/RRTMGP to NumericalEarth",
            coupler = "NumericalEarth/ClimaOcean with ConservativeRegridding",
            coupled_iterations = earth.clock.iteration,
        ),
    )
end

"""Validate finite-state and physical-bound invariants for dynamic-stack output."""
function validate_dynamic_diagnostics(diagnostics)
    finite_fields = (
        :ocean_surface_temperature,
        :ocean_surface_salinity,
        :ocean_river_mouth_receiver_mask,
        :ocean_river_mouth_horizontal_mixing_mask,
        :ocean_surface_zonal_velocity,
        :ocean_surface_meridional_velocity,
        :sea_ice_concentration,
        :sea_ice_thickness,
        :sea_ice_zonal_velocity,
        :sea_ice_meridional_velocity,
        :atmosphere_ocean_sensible_heat_flux,
        :atmosphere_ocean_latent_heat_flux,
        :atmosphere_sea_ice_sensible_heat_flux,
        :atmosphere_sea_ice_latent_heat_flux,
        :ocean_upwelling_longwave,
        :ocean_absorbed_longwave,
        :ocean_absorbed_shortwave,
        :atmosphere_surface_temperature,
        :atmosphere_radiative_surface_temperature,
        :atmosphere_temperature,
        :atmosphere_specific_humidity,
        :atmosphere_convective_humidity_tendency,
        :atmosphere_convective_temperature_tendency,
        :atmosphere_large_scale_condensation_humidity_tendency,
        :atmosphere_large_scale_condensation_temperature_tendency,
        :atmosphere_pressure,
        :atmosphere_zonal_wind,
        :atmosphere_meridional_wind,
        :atmosphere_cloud_fraction,
        :atmosphere_cloud_liquid_water_path,
        :atmosphere_cloud_ice_water_path,
        :atmosphere_surface_pressure,
        :atmosphere_surface_specific_humidity,
        :atmosphere_surface_zonal_wind,
        :atmosphere_surface_meridional_wind,
        :atmosphere_column_water_vapor,
        :atmosphere_rainfall_flux,
        :atmosphere_convective_rainfall_flux,
        :atmosphere_convective_precipitation_mass_correction_flux,
        :atmosphere_large_scale_rainfall_flux,
        :atmosphere_large_scale_snowfall_flux,
        :atmosphere_large_scale_precipitation_mass_correction_flux,
        :atmosphere_snowfall_flux,
        :atmosphere_precipitation_cloud_top_layer,
        :atmosphere_vertical_diffusion_cfl_scale,
        :atmosphere_surface_water_vapor_flux,
        :atmosphere_surface_sensible_heat_flux,
        :outgoing_longwave,
        :column_cloud_fraction,
        :global_surface_temperature,
    )
    for name in finite_fields
        values = getproperty(diagnostics, name)
        if !all(isfinite, values)
            invalid = findall(value -> !isfinite(value), values)
            first_invalid = first(invalid)
            printable_index = first_invalid isa CartesianIndex ?
                Tuple(first_invalid) : first_invalid
            error(
                "dynamic diagnostic $name contains $(length(invalid)) " *
                "non-finite values; first index=$printable_index " *
                "value=$(values[first_invalid])",
            )
        end
    end

    all(0 .<= diagnostics.sea_ice_concentration .<= 1) ||
        error("sea-ice concentration is outside [0, 1]")
    all(diagnostics.sea_ice_thickness .>= 0) ||
        error("sea-ice thickness is negative")
    all(value -> value == 0 || value == 1, diagnostics.ocean_surface_active_mask) ||
        error("ocean surface active mask is not binary")
    size(diagnostics.ocean_cell_area) ==
        size(diagnostics.ocean_surface_active_mask) ||
        error("ocean cell area and surface active mask differ in shape")
    all(isfinite, diagnostics.ocean_cell_area) &&
        all(>(0), diagnostics.ocean_cell_area) ||
        error("ocean cell areas are not finite and positive")
    size(diagnostics.ocean_surface_active_mask) ==
        size(diagnostics.sea_ice_concentration) ||
        error("ocean surface active mask and sea-ice state differ in shape")
    inactive_surface = diagnostics.ocean_surface_active_mask .== 0
    all(iszero, diagnostics.sea_ice_concentration[inactive_surface]) ||
        error("sea-ice concentration is nonzero in inactive ocean columns")
    all(iszero, diagnostics.sea_ice_thickness[inactive_surface]) ||
        error("sea-ice thickness is nonzero in inactive ocean columns")
    all(0 .<= diagnostics.column_cloud_fraction .<= 1) ||
        error("cloud fraction is outside [0, 1]")
    length(diagnostics.atmosphere_sigma_layer_thickness) ==
        length(diagnostics.atmosphere_sigma) || error(
        "atmosphere sigma levels and layer thicknesses differ in length",
    )
    all(>(0), diagnostics.atmosphere_sigma_layer_thickness) && isapprox(
        sum(diagnostics.atmosphere_sigma_layer_thickness),
        1.0;
        atol = 10eps(Float32),
        rtol = 0,
    ) || error("atmosphere sigma-layer thicknesses are invalid")
    cloud_layer_fields = (
        :atmosphere_cloud_fraction,
        :atmosphere_cloud_liquid_water_path,
        :atmosphere_cloud_ice_water_path,
    )
    all(
        size(getproperty(diagnostics, name)) == size(diagnostics.atmosphere_temperature)
        for name in cloud_layer_fields
    ) || error("three-dimensional cloud fields do not match the atmosphere grid")
    layer_process_fields = (
        :atmosphere_convective_humidity_tendency,
        :atmosphere_convective_temperature_tendency,
        :atmosphere_large_scale_condensation_humidity_tendency,
        :atmosphere_large_scale_condensation_temperature_tendency,
    )
    all(
        size(getproperty(diagnostics, name)) == size(diagnostics.atmosphere_temperature)
        for name in layer_process_fields
    ) || error(
        "layer-resolved moist-process fields do not match the atmosphere grid",
    )
    diagnostics.metadata.atmosphere_layer_process_diagnostics ==
        "signed_incremental_parameterization_tendencies_v1" || error(
        "layer-resolved moist-process diagnostics lack exact provenance",
    )
    diagnostics.metadata.atmosphere_cumulative_process_profile_diagnostics ==
        "Gaussian-area-weighted applied-process integral v1" || error(
        "cumulative atmosphere process profiles lack exact provenance",
    )
    process_profile_fields = (
        :cumulative_global_convective_humidity_change_profile,
        :cumulative_global_convective_temperature_change_profile,
        :cumulative_global_large_scale_humidity_change_profile,
        :cumulative_global_large_scale_temperature_change_profile,
        :cumulative_global_radiative_temperature_change_profile,
        :cumulative_global_surface_humidity_change_profile,
        :cumulative_global_surface_sensible_temperature_change_profile,
    )
    expected_process_profile_size = (
        length(diagnostics.atmosphere_sigma),
        length(diagnostics.atmosphere_diagnostic_time_days),
    )
    all(
        size(getproperty(diagnostics, name)) == expected_process_profile_size
        for name in process_profile_fields
    ) || error("cumulative atmosphere process-profile dimensions are inconsistent")
    process_start_day = Float64(
        diagnostics.metadata.atmosphere_process_profile_start_day,
    )
    isfinite(process_start_day) && process_start_day >= 0 || error(
        "atmosphere process-profile start day is invalid",
    )
    process_start_index = findfirst(
        time -> time >= process_start_day - 16eps(Float32),
        diagnostics.atmosphere_diagnostic_time_days,
    )
    isnothing(process_start_index) && error(
        "atmosphere process-profile start lies beyond the diagnostic time axis",
    )
    for name in process_profile_fields
        values = getproperty(diagnostics, name)
        if process_start_index > 1
            all(isnan, @view(values[:, 1:(process_start_index - 1)])) || error(
                "pre-observer atmosphere process profile $name is not missing",
            )
        end
        all(isfinite, @view(values[:, process_start_index:end])) || error(
            "active atmosphere process profile $name is non-finite",
        )
        all(iszero, @view(values[:, process_start_index])) || error(
            "atmosphere process profile $name does not start from zero",
        )
    end
    for name in (
        :cumulative_global_surface_humidity_change_profile,
        :cumulative_global_surface_sensible_temperature_change_profile,
    )
        all(iszero, @view(getproperty(diagnostics, name)[
            1:(end - 1), process_start_index:end,
        ])) || error("surface process profile $name leaks above the lowest layer")
    end
    if diagnostics.metadata.cloud_scheme == "prognostic_condensate"
        diagnostics.metadata.atmosphere_cloud_condensate_state ==
            "column_local_prognostic_liquid_ice_path_conservative_precipitation_delay_v1" ||
            error("prognostic condensate lacks exact provenance")
        abs(diagnostics.metadata.cloud_condensate_ledger_residual_kgm2) <=
            2e-6 || error("global prognostic condensate ledger does not close")
        diagnostics.metadata.cloud_condensate_maximum_column_ledger_residual_kgm2 <=
            2e-5 || error("column prognostic condensate ledger does not close")
        diagnostics.metadata.cloud_condensate_endpoint_kgm2 >= 0 ||
            error("global prognostic condensate storage is negative")
    else
        diagnostics.metadata.atmosphere_cloud_condensate_state ==
            "not_prognostic_diagnostic_cloud_optics_only" || error(
            "diagnostic cloud optics are incorrectly described as prognostic condensate",
        )
    end
    all(0 .<= diagnostics.atmosphere_cloud_fraction .<= 1) ||
        error("three-dimensional cloud fraction is outside [0, 1]")
    all(diagnostics.atmosphere_cloud_liquid_water_path .>= 0) ||
        error("cloud liquid-water path is negative")
    all(diagnostics.atmosphere_cloud_ice_water_path .>= 0) ||
        error("cloud ice-water path is negative")
    vec(maximum(diagnostics.atmosphere_cloud_fraction; dims = 2)) ==
        diagnostics.column_cloud_fraction ||
        error("column and three-dimensional cloud fractions are inconsistent")
    if diagnostics.metadata.ocean_dynamics == 1
        all(isfinite, (
            diagnostics.ocean_max_abs_zonal_velocity,
            diagnostics.ocean_max_abs_meridional_velocity,
            diagnostics.ocean_max_abs_vertical_velocity,
            diagnostics.ocean_mean_kinetic_energy,
        )) || error("ocean circulation summary contains non-finite values")
        diagnostics.ocean_mean_kinetic_energy > 0 ||
            error("ocean dynamics were requested but diagnosed kinetic energy is zero")
    end

    land = diagnostics.atmosphere_land_fraction .> 0
    all(isfinite, diagnostics.land_surface_temperature[land]) ||
        error("land-surface temperature is non-finite on land")
    all(isfinite, diagnostics.land_surface_moisture[land]) ||
        error("land-surface moisture is non-finite on land")
    all(isnan, diagnostics.land_surface_temperature[.!land]) ||
        error("land-surface temperature should be masked over ocean")
    all(isnan, diagnostics.land_surface_moisture[.!land]) ||
        error("land-surface moisture should be masked over ocean")
    all(173.15 .<= diagnostics.land_surface_temperature[land] .<= 353.15) ||
        error("land-surface temperature is outside [-100, 80] degree_Celsius")
    all(0 .<= diagnostics.land_surface_moisture[land] .<= 1) ||
        error("land-surface moisture is outside [0, 1]")

    all(isfinite, (
        diagnostics.ocean_temperature_minimum,
        diagnostics.ocean_temperature_maximum,
        diagnostics.ocean_salinity_minimum,
        diagnostics.ocean_salinity_maximum,
    )) || error("full-depth ocean temperature/salinity extrema are non-finite")
    -5 <= diagnostics.ocean_temperature_minimum <=
        diagnostics.ocean_temperature_maximum <= 40 ||
        error("full-depth ocean temperature is outside [-5, 40] degree_Celsius")
    0 <= diagnostics.ocean_salinity_minimum <=
        diagnostics.ocean_salinity_maximum <= 50 ||
        error("full-depth ocean salinity is outside [0, 50]")
    size(diagnostics.ocean_mean_temperature_profile) == (
        length(diagnostics.ocean_layer_depth),
        length(diagnostics.coupled_budget_time_days),
    ) || error("ocean temperature-profile dimensions are inconsistent")
    size(diagnostics.ocean_mean_salinity_profile) ==
        size(diagnostics.ocean_mean_temperature_profile) ||
        error("ocean temperature and salinity profile dimensions differ")
    all(isfinite, diagnostics.ocean_mean_temperature_profile) ||
        error("ocean mean-temperature profile contains non-finite values")
    all(isfinite, diagnostics.ocean_mean_salinity_profile) ||
        error("ocean mean-salinity profile contains non-finite values")

    if diagnostics.metadata.land == "Terrarium"
        isempty(diagnostics.terrarium_land_atmosphere_point) &&
            error("Terrarium final-state diagnostics are empty")
        terrarium_finite_fields = (
            diagnostics.terrarium_soil_temperature,
            diagnostics.terrarium_soil_saturation,
            diagnostics.terrarium_skin_temperature,
            diagnostics.terrarium_surface_excess_water,
            diagnostics.terrarium_surface_runoff,
            diagnostics.terrarium_infiltration,
        )
        all(all(isfinite, values) for values in terrarium_finite_fields) ||
            error("Terrarium final state contains non-finite values")
        all(-100 .<= diagnostics.terrarium_soil_temperature .<= 80) ||
            error("Terrarium soil temperature is outside [-100, 80] degree_Celsius")
        all(-100 .<= diagnostics.terrarium_skin_temperature .<= 80) ||
            error("Terrarium skin temperature is outside [-100, 80] degree_Celsius")
        all(0 .<= diagnostics.terrarium_soil_saturation .<= 1) ||
            error("Terrarium soil saturation is outside [0, 1]")
        size(diagnostics.global_land_soil_layer_temperature) == (
            length(diagnostics.terrarium_soil_layer_depth),
            length(diagnostics.atmosphere_diagnostic_time_days),
        ) || error("global soil-temperature profile dimensions are inconsistent")
        size(diagnostics.global_land_soil_layer_saturation) ==
            size(diagnostics.global_land_soil_layer_temperature) ||
            error("global soil temperature and saturation profiles differ")
        all(isfinite, diagnostics.global_land_soil_layer_temperature) ||
            error("global soil-temperature profile contains non-finite values")
        all(-100 .<= diagnostics.global_land_soil_layer_temperature .<= 80) ||
            error("global soil-temperature profile is outside [-100, 80] degree_Celsius")
        fraction_reduction_tolerance = 8eps(Float32)
        all(
            -fraction_reduction_tolerance .<=
            diagnostics.global_land_soil_layer_saturation .<=
            1 + fraction_reduction_tolerance,
        ) || error(
            "global soil-saturation profile is outside [0, 1] beyond " *
            "Float32 reduction tolerance",
        )
        # DirectSurfaceRunoff removes ponded water on a one-hour timescale.
        # Five centimetres is therefore already a deliberately generous
        # coarse-grid safety bound; the former porosity/race defect produced
        # 0.257 m in one day and must fail before contaminating coastal ocean
        # salinity.
        all(0 .<= diagnostics.terrarium_surface_excess_water .<= 0.05) ||
            error("Terrarium surface excess water is outside [0, 0.05] m")
        all(0 .<= diagnostics.terrarium_surface_runoff .<= 1e-3) ||
            error("Terrarium surface runoff is outside [0, 1e-3] m/s")
        all(0 .<= diagnostics.terrarium_infiltration .<= 1e-3) ||
            error("Terrarium infiltration is outside [0, 1e-3] m/s")
        if diagnostics.metadata.terrarium_evapotranspiration ==
           "era5_prescribed_vegetation"
            vegetation_fields = (
                diagnostics.terrarium_vegetation_fraction,
                diagnostics.terrarium_leaf_area_index,
                diagnostics.terrarium_root_fraction,
                diagnostics.terrarium_root_water_availability,
                diagnostics.terrarium_root_uptake_weight_sum,
                diagnostics.terrarium_ground_evaporation_flux,
                diagnostics.terrarium_transpiration_flux,
            )
            all(!isempty, vegetation_fields) ||
                error("prescribed-vegetation diagnostics are empty")
            all(all(isfinite, values) for values in vegetation_fields) ||
                error("prescribed-vegetation diagnostics contain non-finite values")
            all(0 .<= diagnostics.terrarium_vegetation_fraction .<= 1) ||
                error("Terrarium vegetation fraction is outside [0, 1]")
            all(0 .<= diagnostics.terrarium_leaf_area_index .<= 10) ||
                error("Terrarium LAI is outside [0, 10]")
            all(0 .<= diagnostics.terrarium_root_fraction .<= 1) ||
                error("Terrarium root fraction is outside [0, 1]")
            all(abs.(sum(diagnostics.terrarium_root_fraction; dims = 2) .- 1) .<=
                5e-6) || error("Terrarium root fractions do not sum to one")
            all(0 .<= diagnostics.terrarium_root_water_availability .<= 1) ||
                error("Terrarium root-water availability is outside [0, 1]")
            all(0 .<= diagnostics.terrarium_root_uptake_weight_sum .<= 1) ||
                error("Terrarium root uptake weights are outside [0, 1]")
            all(diagnostics.terrarium_transpiration_flux .>= 0) ||
                error("Terrarium transpiration is negative")
        end
    end

    budgets = (
        diagnostics.toa_incoming_shortwave,
        diagnostics.toa_outgoing_shortwave,
        diagnostics.toa_outgoing_longwave,
        diagnostics.toa_clear_outgoing_shortwave,
        diagnostics.toa_clear_outgoing_longwave,
        diagnostics.toa_net_downward,
        diagnostics.toa_clear_net_downward,
    )
    all(length(values) == length(diagnostics.time_days) for values in budgets) ||
        error("TOA budget and time dimensions differ")
    all(isnan(values[1]) && all(isfinite, values[2:end]) for values in budgets) ||
        error("TOA budgets must be NaN initially and finite after completed steps")

    minimum(diagnostics.atmosphere_surface_pressure) > 0 ||
        error("surface pressure is not positive")
    all(diagnostics.atmosphere_rainfall_flux .>= 0) ||
        error("rainfall flux is negative")
    all(diagnostics.atmosphere_convective_rainfall_flux .>= 0) ||
        error("convective rainfall flux is negative")
    all(diagnostics.atmosphere_large_scale_rainfall_flux .>= 0) ||
        error("large-scale rainfall flux is negative")
    all(diagnostics.atmosphere_large_scale_snowfall_flux .>= 0) ||
        error("large-scale snowfall flux is negative")
    all(diagnostics.atmosphere_snowfall_flux .>= 0) ||
        error("snowfall flux is negative")
    precipitation_component_tolerance = 1e-10
    all(abs.(
        diagnostics.atmosphere_rainfall_flux .-
        diagnostics.atmosphere_convective_rainfall_flux .-
        diagnostics.atmosphere_large_scale_rainfall_flux
    ) .<= precipitation_component_tolerance) || error(
        "rainfall flux does not equal convective plus large-scale components",
    )
    all(abs.(
        diagnostics.atmosphere_snowfall_flux .-
        diagnostics.atmosphere_large_scale_snowfall_flux
    ) .<= precipitation_component_tolerance) || error(
        "snowfall flux does not equal its large-scale component",
    )
    layer_thickness = reshape(
        diagnostics.atmosphere_sigma_layer_thickness,
        1,
        :,
    )
    isfinite(diagnostics.metadata.atmosphere_gravity_ms2) &&
        diagnostics.metadata.atmosphere_gravity_ms2 > 0 || error(
        "atmosphere gravitational acceleration is invalid",
    )
    surface_column_mass = diagnostics.atmosphere_surface_pressure ./
        diagnostics.metadata.atmosphere_gravity_ms2
    convective_precipitation_from_tendency = max.(
        -surface_column_mass .* vec(sum(
            diagnostics.atmosphere_convective_humidity_tendency .*
            layer_thickness;
            dims = 2,
        )),
        0,
    )
    large_scale_precipitation_from_tendency = max.(
        -surface_column_mass .* vec(sum(
            diagnostics.atmosphere_large_scale_condensation_humidity_tendency .*
            layer_thickness;
            dims = 2,
        )),
        0,
    )
    # The layer fields are Float32 differences of total accumulated
    # tendencies, while the surface rate is formed from the more accurate
    # before/after column sums inside the parameterization kernel. Allow the
    # measured cancellation envelope (3e-9 m/s of water) without admitting
    # the order-1e-8 m/s upstream re-evaporation loss.
    process_precipitation_tolerance = 3e-6
    if diagnostics.metadata.atmosphere_cloud_condensate_state ==
       "not_prognostic_diagnostic_cloud_optics_only"
        all(abs.(
            convective_precipitation_from_tendency .-
            diagnostics.atmosphere_convective_rainfall_flux
        ) .<= process_precipitation_tolerance) || error(
            "convective layer humidity tendency does not reproduce surface rain",
        )
    elseif diagnostics.metadata.atmosphere_cloud_condensate_state !=
           "column_local_prognostic_liquid_ice_path_conservative_precipitation_delay_v1"
        error("unsupported atmosphere cloud-condensate provenance")
    end
    large_scale_scheme =
        diagnostics.metadata.atmosphere_large_scale_precipitation_scheme
    if large_scale_scheme ==
       "implicit condensation net column humidity drying"
        all(abs.(
            large_scale_precipitation_from_tendency .-
            diagnostics.atmosphere_large_scale_rainfall_flux .-
            diagnostics.atmosphere_large_scale_snowfall_flux
        ) .<= process_precipitation_tolerance) || error(
            "conservative large-scale layer humidity tendency does not " *
            "reproduce rain plus snow",
        )
    elseif large_scale_scheme ==
           "SpeedyWeather implicit condensation observed without surface-rate correction"
        all(iszero,
            diagnostics.atmosphere_large_scale_precipitation_mass_correction_flux) ||
            error("upstream-observed large-scale precipitation applied a correction")
    else
        error("unsupported large-scale precipitation provenance $large_scale_scheme")
    end
    cloud_top = diagnostics.atmosphere_precipitation_cloud_top_layer
    all(isinteger, cloud_top) ||
        error("precipitation cloud-top layer is not integer-valued")
    all(1 .<= cloud_top .<= length(diagnostics.atmosphere_sigma) + 1) ||
        error("precipitation cloud-top layer is outside the atmosphere")
    atmosphere_series = (
        diagnostics.global_rainfall_flux,
        diagnostics.global_convective_rainfall_flux,
        diagnostics.global_large_scale_rainfall_flux,
        diagnostics.global_large_scale_snowfall_flux,
        diagnostics.global_snowfall_flux,
        diagnostics.global_surface_water_vapor_flux,
        diagnostics.cumulative_global_surface_water_vapor,
        diagnostics.cumulative_global_precipitation,
        diagnostics.cumulative_balanced_launch_atmosphere_water_source,
        diagnostics.global_atmosphere_water_budget_residual,
        diagnostics.global_cloud_condensate_water_path,
        diagnostics.global_total_atmosphere_water_budget_residual,
        diagnostics.global_surface_air_temperature,
        diagnostics.global_surface_specific_humidity,
        diagnostics.global_surface_pressure,
        diagnostics.global_near_surface_wind_speed,
        diagnostics.global_column_water_vapor,
        diagnostics.global_mass_weighted_atmosphere_temperature,
        diagnostics.global_land_surface_soil_moisture,
        diagnostics.global_land_total_water_storage,
        diagnostics.cumulative_global_land_precipitation,
        diagnostics.cumulative_global_land_evapotranspiration,
        diagnostics.cumulative_global_land_surface_runoff,
        diagnostics.cumulative_balanced_launch_land_water_source,
        diagnostics.global_land_water_budget_residual,
        diagnostics.global_land_rainfall_flux,
        diagnostics.global_land_snowfall_flux,
        diagnostics.global_land_evaporation_flux,
        diagnostics.global_land_ground_evaporation_flux,
        diagnostics.global_land_transpiration_flux,
        diagnostics.global_land_surface_runoff_flux,
        diagnostics.global_land_infiltration_flux,
    )
    all(length(values) == length(diagnostics.atmosphere_diagnostic_time_days)
        for values in atmosphere_series) ||
        error("atmospheric diagnostic series and time dimensions differ")
    precipitation_series = (
        diagnostics.global_rainfall_flux,
        diagnostics.global_convective_rainfall_flux,
        diagnostics.global_large_scale_rainfall_flux,
        diagnostics.global_large_scale_snowfall_flux,
        diagnostics.global_snowfall_flux,
        diagnostics.global_land_rainfall_flux,
        diagnostics.global_land_snowfall_flux,
    )
    all(isnan(values[1]) && all(isfinite, values[2:end])
        for values in precipitation_series) ||
        error("precipitation series must be NaN initially and finite afterwards")
    all(abs.(
        diagnostics.global_rainfall_flux[2:end] .-
        diagnostics.global_convective_rainfall_flux[2:end] .-
        diagnostics.global_large_scale_rainfall_flux[2:end]
    ) .<= precipitation_component_tolerance) || error(
        "global rainfall series does not close across process components",
    )
    all(abs.(
        diagnostics.global_snowfall_flux[2:end] .-
        diagnostics.global_large_scale_snowfall_flux[2:end]
    ) .<= precipitation_component_tolerance) || error(
        "global snowfall series does not close across process components",
    )
    isnan(diagnostics.global_surface_air_temperature[1]) &&
        all(isfinite, diagnostics.global_surface_air_temperature[2:end]) ||
        error("surface-air temperature must be NaN initially and finite afterwards")
    finite_atmosphere_series = (
        diagnostics.global_surface_specific_humidity,
        diagnostics.global_surface_pressure,
        diagnostics.global_near_surface_wind_speed,
        diagnostics.global_column_water_vapor,
        diagnostics.global_mass_weighted_atmosphere_temperature,
        diagnostics.global_land_surface_soil_moisture,
        diagnostics.global_land_evaporation_flux,
        diagnostics.global_land_ground_evaporation_flux,
        diagnostics.global_land_transpiration_flux,
        diagnostics.global_surface_water_vapor_flux,
        diagnostics.cumulative_global_surface_water_vapor,
        diagnostics.cumulative_global_precipitation,
        diagnostics.cumulative_balanced_launch_atmosphere_water_source,
        diagnostics.global_atmosphere_water_budget_residual,
        diagnostics.global_cloud_condensate_water_path,
        diagnostics.global_total_atmosphere_water_budget_residual,
    )
    all(all(isfinite, values) for values in finite_atmosphere_series) ||
        error("global atmospheric diagnostic series contains non-finite values")
    all(diagnostics.cumulative_global_precipitation .>= 0) ||
        error("cumulative global precipitation is negative")
    all(diff(diagnostics.cumulative_global_precipitation) .>= -1e-5) ||
        error("cumulative global precipitation is not monotonic")
    abs(first(diagnostics.cumulative_global_surface_water_vapor)) <= 1e-6 &&
        abs(first(diagnostics.cumulative_global_precipitation)) <= 1e-6 &&
        abs(first(
            diagnostics.cumulative_balanced_launch_atmosphere_water_source,
        )) <= 1e-6 &&
        abs(first(diagnostics.global_atmosphere_water_budget_residual)) <= 1e-6 &&
        abs(first(diagnostics.global_cloud_condensate_water_path)) <= 1e-6 &&
        abs(first(
            diagnostics.global_total_atmosphere_water_budget_residual,
        )) <= 1e-6 ||
        error("atmospheric cumulative water diagnostics do not start at zero")
    all(diagnostics.global_land_transpiration_flux .>= 0) ||
        error("global land transpiration is negative")
    terrarium_series = (
        diagnostics.global_land_total_water_storage,
        diagnostics.cumulative_global_land_precipitation,
        diagnostics.cumulative_global_land_evapotranspiration,
        diagnostics.cumulative_global_land_surface_runoff,
        diagnostics.cumulative_balanced_launch_land_water_source,
        diagnostics.global_land_water_budget_residual,
        diagnostics.global_land_surface_runoff_flux,
        diagnostics.global_land_infiltration_flux,
    )
    if diagnostics.metadata.land == "Terrarium"
        all(all(isfinite, values) for values in terrarium_series) ||
            error("Terrarium diagnostic series contains non-finite values")
        all(diagnostics.global_land_total_water_storage .>= 0) ||
            error("global land total-water storage is negative")
    else
        all(all(isnan, values) for values in terrarium_series) ||
            error("Terrarium-only diagnostics must be NaN without Terrarium")
    end
    all(0 .<= diagnostics.global_land_surface_soil_moisture .<= 1) ||
        error("global land-surface soil saturation is outside [0, 1]")
    if diagnostics.metadata.land == "Terrarium"
        land_nonnegative_cumulative_series = (
            diagnostics.cumulative_global_land_precipitation,
            diagnostics.cumulative_global_land_surface_runoff,
        )
        all(all(values .>= 0) && all(diff(values) .>= -1e-5)
            for values in land_nonnegative_cumulative_series) ||
            error("cumulative Terrarium precipitation/runoff is negative or non-monotonic")
        land_cumulative_series = (
            land_nonnegative_cumulative_series...,
            diagnostics.cumulative_global_land_evapotranspiration,
        )
        all(abs(first(values)) <= 1e-6 for values in land_cumulative_series) &&
            abs(first(
                diagnostics.cumulative_balanced_launch_land_water_source,
            )) <= 1e-6 &&
            abs(first(diagnostics.global_land_water_budget_residual)) <= 1e-6 ||
            error("Terrarium cumulative water diagnostics do not start at zero")
        land_closure_tolerance = 0.25 +
            0.002 * maximum(diagnostics.atmosphere_diagnostic_time_days)
        maximum(abs.(diagnostics.global_land_water_budget_residual)) <=
            land_closure_tolerance ||
            error(
                "Terrarium cumulative water-budget residual exceeds " *
                "$(land_closure_tolerance) kg/m^2 of land",
            )
    end
    # A single-step backend/debug gate can legitimately have no resolvable net
    # hydrological change (for example, no rain/runoff and a change smaller
    # than the diagnostic accumulator precision). Retain the evolution check
    # for every multi-step and production trajectory, where a static land-water
    # series is evidence that the prognostic land state is not being advanced.
    if diagnostics.metadata.land == "Terrarium" &&
       diagnostics.metadata.coupled_iterations > 1
        maximum(abs.(
            diagnostics.global_land_total_water_storage .-
            first(diagnostics.global_land_total_water_storage),
        )) > 0 || error("Terrarium land water storage did not evolve")
    end
    length(diagnostics.balanced_launch_applied_fraction) ==
        length(diagnostics.coupled_budget_time_days) || error(
        "balanced-launch fraction and coupled budget time dimensions differ",
    )
    size(diagnostics.balanced_launch_cumulative_sources) == (
        length(_BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS),
        length(diagnostics.coupled_budget_time_days),
    ) || error("balanced-launch source and coupled budget dimensions differ")
    all(isfinite, diagnostics.balanced_launch_applied_fraction) &&
        all(0 .<= diagnostics.balanced_launch_applied_fraction .<= 1) &&
        all(diff(diagnostics.balanced_launch_applied_fraction) .>= -64eps(Float64)) ||
        error("balanced-launch applied fraction is invalid or non-monotonic")
    all(isfinite, diagnostics.balanced_launch_cumulative_sources) ||
        error("balanced-launch cumulative source history is non-finite")
    all(iszero, view(
        diagnostics.balanced_launch_cumulative_sources,
        :,
        1,
    )) || error("balanced-launch cumulative sources do not start at zero")
    isapprox(
        last(diagnostics.balanced_launch_applied_fraction),
        diagnostics.metadata.balanced_launch_final_applied_fraction;
        atol = 64eps(Float64),
        rtol = 0,
    ) || error("balanced-launch history and metadata fractions differ")
    if diagnostics.metadata.balanced_launch == "none"
        all(iszero, diagnostics.balanced_launch_applied_fraction) &&
            all(iszero, diagnostics.balanced_launch_cumulative_sources) ||
            error("disabled balanced launch has nonzero source diagnostics")
    else
        isfinite(diagnostics.metadata.balanced_launch_start_time_seconds) &&
            isfinite(diagnostics.metadata.balanced_launch_window_seconds) &&
            diagnostics.metadata.balanced_launch_window_seconds > 0 ||
            error("installed balanced launch has invalid timing metadata")
        diagnostics.metadata.balanced_launch_target_name != "none" &&
            diagnostics.metadata.balanced_launch_provenance_sha256 != "none" ||
            error("installed balanced launch lacks target provenance")
    end
    forecast_launch_schema =
        "readiesm_forecast_background_balanced_launch_v1"
    if diagnostics.metadata.balanced_launch == forecast_launch_schema
        diagnostics.metadata.balanced_launch_reference_mode ==
            "free-coupled-forecast-endpoint" || error(
            "forecast-background launch lacks its reference mode",
        )
        diagnostics.metadata.balanced_launch_reference_background_name !=
            "none" || error(
            "forecast-background launch lacks its background identity",
        )
        reference_time = diagnostics.metadata.
            balanced_launch_reference_forecast_end_time_seconds
        reference_iteration = diagnostics.metadata.
            balanced_launch_reference_forecast_end_iteration
        isfinite(reference_time) && reference_time > 0 ||
            error("forecast-background launch has an invalid reference time")
        reference_iteration > 0 ||
            error("forecast-background launch has an invalid reference iteration")
        increment_digest = diagnostics.metadata.
            balanced_launch_reference_increment_sha256
        length(increment_digest) == 64 && all(isxdigit, increment_digest) ||
            error("forecast-background launch has an invalid increment digest")
    else
        diagnostics.metadata.balanced_launch_reference_mode == "none" ||
            error("non-forecast launch declares a forecast reference")
        diagnostics.metadata.balanced_launch_reference_increment_sha256 ==
            "none" || error(
            "non-forecast launch declares a forecast increment digest",
        )
    end
    coupled_budget_series = (
        diagnostics.balanced_launch_applied_fraction,
        diagnostics.ocean_heat_content,
        diagnostics.ocean_salt_content,
        diagnostics.ocean_volume,
        diagnostics.ocean_temperature_minimum_series,
        diagnostics.ocean_temperature_maximum_series,
        diagnostics.ocean_salinity_minimum_series,
        diagnostics.ocean_salinity_maximum_series,
        diagnostics.ocean_surface_salinity_minimum_series,
        diagnostics.ocean_surface_salinity_minimum_longitude,
        diagnostics.ocean_surface_salinity_minimum_latitude,
        diagnostics.ocean_river_mouth_active_mixing_area_m2,
        diagnostics.sea_ice_volume,
        diagnostics.sea_ice_area,
        diagnostics.tropical_sea_ice_area,
        diagnostics.sea_ice_max_thickness,
        diagnostics.sea_ice_extent_max_thickness,
        diagnostics.sea_ice_max_area_equivalent_thickness,
        diagnostics.sea_ice_marginal_ceiling_max_concentration,
        diagnostics.sea_ice_marginal_ceiling_max_area_equivalent_thickness,
        diagnostics.ocean_free_surface_minimum,
        diagnostics.ocean_free_surface_maximum,
        diagnostics.ocean_max_abs_zonal_velocity_series,
        diagnostics.ocean_max_abs_meridional_velocity_series,
        diagnostics.cumulative_ocean_surface_heat,
        diagnostics.cumulative_ocean_frazil_heat,
        diagnostics.cumulative_ocean_penetrating_shortwave_heat,
        diagnostics.cumulative_net_ocean_heat,
        diagnostics.cumulative_ocean_freshwater_enthalpy,
        diagnostics.cumulative_net_ocean_freshwater_volume,
        diagnostics.cumulative_net_ocean_salinity,
        diagnostics.cumulative_rain_mass,
        diagnostics.cumulative_snow_mass,
        diagnostics.cumulative_open_ocean_evaporation_mass,
        diagnostics.cumulative_sea_ice_ocean_freshwater_volume,
        diagnostics.cumulative_routed_land_runoff_mass,
        diagnostics.ocean_heat_closure_residual,
        diagnostics.ocean_heat_closure_relative,
        diagnostics.ocean_freshwater_closure_residual,
        diagnostics.ocean_freshwater_closure_relative,
        diagnostics.ocean_salt_closure_residual,
        diagnostics.ocean_salt_closure_relative,
    )
    all(length(values) == length(diagnostics.coupled_budget_time_days)
        for values in coupled_budget_series) ||
        error("coupled budget series and time dimensions differ")
    length(diagnostics.sea_ice_marginal_ceiling_cell_count) ==
        length(diagnostics.coupled_budget_time_days) ||
        error("marginal sea-ice ceiling counts and time dimensions differ")
    length(diagnostics.ocean_river_mouth_active_mixing_cells) ==
        length(diagnostics.coupled_budget_time_days) || error(
        "active river-mouth counts and coupled-budget times differ",
    )
    all(
        length(values) == length(diagnostics.coupled_budget_time_days)
        for values in (
            diagnostics.ocean_surface_salinity_minimum_i,
            diagnostics.ocean_surface_salinity_minimum_j,
        )
    ) || error("surface-salinity minimum indices and time dimensions differ")
    all(all(isfinite, values) for values in coupled_budget_series) ||
        error("coupled conservation budget contains non-finite values")
    all(
        diagnostics.ocean_temperature_minimum_series .<=
        diagnostics.ocean_temperature_maximum_series,
    ) || error("ocean temperature minimum exceeds maximum")
    all(
        diagnostics.ocean_salinity_minimum_series .<=
        diagnostics.ocean_salinity_maximum_series,
    ) || error("ocean salinity minimum exceeds maximum")
    all(
        diagnostics.ocean_salinity_minimum_series .<=
        diagnostics.ocean_surface_salinity_minimum_series .<=
        diagnostics.ocean_salinity_maximum_series,
    ) || error("surface-salinity minimum is inconsistent with full-depth extrema")
    all(1 .<= diagnostics.ocean_surface_salinity_minimum_i .<=
        size(diagnostics.ocean_surface_active_mask, 1)) &&
        all(1 .<= diagnostics.ocean_surface_salinity_minimum_j .<=
            size(diagnostics.ocean_surface_active_mask, 2)) || error(
        "surface-salinity minimum indices are outside the ocean grid",
    )
    all(-90 .<= diagnostics.ocean_surface_salinity_minimum_latitude .<= 90) ||
        error("surface-salinity minimum latitude is outside [-90, 90]")
    runoff_receiver_series = (
        diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series,
        diagnostics.ocean_runoff_receiver_surface_salinity_minimum_longitude,
        diagnostics.ocean_runoff_receiver_surface_salinity_minimum_latitude,
        diagnostics.ocean_runoff_receiver_surface_salinity_minimum_i,
        diagnostics.ocean_runoff_receiver_surface_salinity_minimum_j,
    )
    all(
        length(values) == length(diagnostics.coupled_budget_time_days)
        for values in runoff_receiver_series
    ) || error("runoff-receiver salinity series and time dimensions differ")
    has_runoff_router = diagnostics.metadata.land == "Terrarium"
    if has_runoff_router
        all(isfinite, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series) &&
            all(isfinite, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_longitude) &&
            all(isfinite, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_latitude) ||
            error("runoff-receiver salinity minimum contains non-finite values")
        all(
            diagnostics.ocean_surface_salinity_minimum_series .<=
            diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series .<=
            diagnostics.ocean_salinity_maximum_series,
        ) || error("runoff-receiver salinity minimum is inconsistent with ocean extrema")
        all(-90 .<=
            diagnostics.ocean_runoff_receiver_surface_salinity_minimum_latitude .<=
            90) || error("runoff-receiver salinity latitude is outside [-90, 90]")
        all(1 .<= diagnostics.ocean_runoff_receiver_surface_salinity_minimum_i .<=
            size(diagnostics.ocean_surface_active_mask, 1)) &&
            all(1 .<= diagnostics.ocean_runoff_receiver_surface_salinity_minimum_j .<=
                size(diagnostics.ocean_surface_active_mask, 2)) || error(
            "runoff-receiver salinity indices are outside the ocean grid",
        )
        last_receiver = CartesianIndex(
            last(diagnostics.ocean_runoff_receiver_surface_salinity_minimum_i),
            last(diagnostics.ocean_runoff_receiver_surface_salinity_minimum_j),
        )
        diagnostics.ocean_runoff_receiver_mask[last_receiver] == 1 || error(
            "last runoff-receiver salinity minimum is outside the receiver mask",
        )
        isapprox(
            diagnostics.ocean_surface_salinity[last_receiver],
            last(diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series);
            atol = 1e-6,
            rtol = 0,
        ) && isapprox(
            minimum(diagnostics.ocean_surface_salinity[
                diagnostics.ocean_runoff_receiver_mask .== 1
            ]),
            last(diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series);
            atol = 1e-6,
            rtol = 0,
        ) || error(
            "last runoff-receiver salinity minimum does not match the final surface field",
        )
    else
        all(isnan, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series) &&
            all(isnan, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_longitude) &&
            all(isnan, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_latitude) &&
            all(iszero, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_i) &&
            all(iszero, diagnostics.ocean_runoff_receiver_surface_salinity_minimum_j) ||
            error("runoff-receiver diagnostics must be empty sentinels without Terrarium")
    end
    issorted(diagnostics.coupled_budget_time_days) ||
        error("coupled budget times are not monotonic")
    maximum(diagnostics.tropical_sea_ice_area) <= 1 ||
        error(
            "material sea ice formed within " *
            "$(_TROPICAL_SEA_ICE_MAXIMUM_ABSOLUTE_LATITUDE_DEGREES) " *
            "degrees of the equator",
        )
    all(diagnostics.sea_ice_max_thickness .>= 0) ||
        error("maximum sea-ice thickness series is negative")
    all(value -> value == 0 || value == 1,
        diagnostics.ocean_runoff_receiver_mask) || error(
        "ocean runoff-receiver mask is not binary",
    )
    size(diagnostics.ocean_runoff_receiver_mask) ==
        size(diagnostics.ocean_surface_active_mask) &&
        all(diagnostics.ocean_runoff_receiver_mask .<=
            diagnostics.ocean_surface_active_mask) || error(
        "ocean runoff-receiver mask extends outside the active surface ocean",
    )
    diagnostics.metadata.ocean_runoff_receiver_cells ==
        count(!iszero, diagnostics.ocean_runoff_receiver_mask) || error(
        "ocean runoff-receiver mask count disagrees with metadata",
    )
    if has_runoff_router
        any(!iszero, diagnostics.ocean_runoff_receiver_mask) || error(
            "Terrarium runoff-receiver mask is empty",
        )
    else
        all(iszero, diagnostics.ocean_runoff_receiver_mask) || error(
            "ocean runoff-receiver mask exists without a land router",
        )
    end
    all(value -> value == 0 || value == 1,
        diagnostics.ocean_river_mouth_receiver_mask) || error(
        "river-mouth receiver mask is not binary",
    )
    all(value -> value == 0 || value == 1,
        diagnostics.ocean_river_mouth_horizontal_mixing_mask) || error(
        "river-mouth horizontal mixing mask is not binary",
    )
    all(value -> value == 0 || value == 1,
        diagnostics.ocean_river_mouth_active_mixing_mask) || error(
        "river-mouth active mixing mask is not binary",
    )
    size(diagnostics.ocean_routed_land_runoff_flux) ==
        size(diagnostics.ocean_runoff_receiver_mask) &&
        all(isfinite, diagnostics.ocean_routed_land_runoff_flux) &&
        all(diagnostics.ocean_routed_land_runoff_flux .>= 0) &&
        all(
            (diagnostics.ocean_routed_land_runoff_flux .> 0) .<=
            (diagnostics.ocean_runoff_receiver_mask .== 1),
        ) || error("routed land runoff flux is invalid or outside receivers")
    all(
        diagnostics.ocean_river_mouth_receiver_mask .<=
        diagnostics.ocean_river_mouth_horizontal_mixing_mask,
    ) || error("river-mouth horizontal footprint excludes a runoff receiver")
    size(diagnostics.ocean_river_mouth_active_mixing_mask) ==
        size(diagnostics.ocean_river_mouth_receiver_mask) && all(
        diagnostics.ocean_river_mouth_active_mixing_mask .<=
        diagnostics.ocean_river_mouth_receiver_mask,
    ) || error("active river-mouth mixing extends outside routed receivers")
    all(0 .<= diagnostics.ocean_river_mouth_active_mixing_cells .<=
        diagnostics.metadata.ocean_river_mouth_receiver_cells) || error(
        "active river-mouth count is outside the routed footprint",
    )
    all(isfinite, diagnostics.ocean_river_mouth_active_mixing_area_m2) &&
        all(0 .<= diagnostics.ocean_river_mouth_active_mixing_area_m2 .<=
            sum(diagnostics.ocean_cell_area)) || error(
        "active river-mouth area is non-finite or outside the ocean area",
    )
    last(diagnostics.ocean_river_mouth_active_mixing_cells) ==
        count(!iszero, diagnostics.ocean_river_mouth_active_mixing_mask) ||
        error("final active river-mouth mask and count history differ")
    isapprox(
        last(diagnostics.ocean_river_mouth_active_mixing_area_m2),
        sum(diagnostics.ocean_cell_area[
            diagnostics.ocean_river_mouth_active_mixing_mask .== 1
        ]);
        atol = 1e-6,
        rtol = 8eps(Float64),
    ) || error("final active river-mouth mask and area history differ")
    river_mouth_mixing_mode = Symbol(
        diagnostics.metadata.ocean_river_mouth_mixing,
    )
    if _is_localized_estuary_mixing(river_mouth_mixing_mode)
        diagnostics.ocean_river_mouth_receiver_mask ==
            diagnostics.ocean_runoff_receiver_mask || error(
            "localized mixing footprint differs from the actual runoff receivers",
        )
        all(
            diagnostics.ocean_river_mouth_horizontal_mixing_mask .<=
            diagnostics.ocean_surface_active_mask,
        ) || error("river-mouth horizontal footprint extends onto dry cells")
        if river_mouth_mixing_mode == :localized_estuary_v1
            count(!iszero,
                diagnostics.ocean_river_mouth_horizontal_mixing_mask) >
                count(!iszero,
                    diagnostics.ocean_river_mouth_receiver_mask) || error(
                "localized_estuary_v1 footprint lacks its wet cardinal ring",
            )
        elseif river_mouth_mixing_mode in (
            :localized_estuary_v2,
            :localized_estuary_v3,
            :localized_estuary_v4,
            :localized_estuary_v5,
        )
            diagnostics.ocean_river_mouth_horizontal_mixing_mask ==
                diagnostics.ocean_river_mouth_receiver_mask || error(
                "$river_mouth_mixing_mode coefficients are not receiver-centred",
            )
        end
        if river_mouth_mixing_mode == :localized_estuary_v4
            diagnostics.metadata.ocean_river_mouth_mixing_footprint ==
                _river_mouth_mixing_footprint_provenance(
                    :localized_estuary_v4,
                ) &&
                diagnostics.metadata.ocean_river_mouth_mixing_dynamic_gating ==
                    "positive_routed_freshwater_flux_each_coupled_step" ||
                error("localized_estuary_v4 lacks dynamic-gating provenance")
            diagnostics.ocean_river_mouth_active_mixing_mask == Float64.(
                diagnostics.ocean_routed_land_runoff_flux .> 0,
            ) || error(
                "localized_estuary_v4 activity does not match final positive runoff",
            )
        elseif river_mouth_mixing_mode == :localized_estuary_v5
            reference_flux = diagnostics.metadata.
                ocean_river_mouth_reference_freshwater_mass_flux_kgm2s
            isfinite(reference_flux) && reference_flux > 0 &&
                diagnostics.metadata.ocean_river_mouth_mixing_footprint ==
                    _river_mouth_mixing_footprint_provenance(
                        :localized_estuary_v5,
                    ) &&
                diagnostics.metadata.ocean_river_mouth_mixing_dynamic_gating ==
                    "smooth_clamp_flux_over_reference_minus_0.5_each_coupled_step" ||
                error("localized_estuary_v5 lacks loading-gate provenance")
            diagnostics.ocean_river_mouth_active_mixing_mask == Float64.(
                diagnostics.ocean_routed_land_runoff_flux .> 0.5reference_flux,
            ) || error(
                "localized_estuary_v5 activity does not match final routed loading",
            )
        else
            diagnostics.metadata.ocean_river_mouth_mixing_dynamic_gating ==
                "none" || error(
                "$river_mouth_mixing_mode incorrectly records dynamic gating",
            )
            diagnostics.ocean_river_mouth_active_mixing_mask ==
                diagnostics.ocean_river_mouth_receiver_mask && all(
                diagnostics.ocean_river_mouth_active_mixing_cells .==
                diagnostics.metadata.ocean_river_mouth_receiver_cells,
            ) || error(
                "$river_mouth_mixing_mode should keep every receiver active",
            )
        end
    else
        all(iszero, diagnostics.ocean_river_mouth_receiver_mask) &&
            all(iszero, diagnostics.ocean_river_mouth_horizontal_mixing_mask) &&
            all(iszero, diagnostics.ocean_river_mouth_active_mixing_mask) &&
            all(iszero, diagnostics.ocean_river_mouth_active_mixing_cells) &&
            all(iszero, diagnostics.ocean_river_mouth_active_mixing_area_m2) ||
            error("disabled localized mixing has a non-empty footprint")
    end
    all(diagnostics.sea_ice_extent_max_thickness .>= 0) ||
        error("extent-qualified maximum sea-ice thickness series is negative")
    all(diagnostics.sea_ice_max_area_equivalent_thickness .>= 0) ||
        error("maximum area-equivalent sea-ice thickness series is negative")
    all(diagnostics.sea_ice_marginal_ceiling_cell_count .>= 0) ||
        error("marginal sea-ice ceiling cell count is negative")
    all(0 .<= diagnostics.sea_ice_marginal_ceiling_max_concentration .<
        _SEA_ICE_EXTENT_CONCENTRATION_THRESHOLD) ||
        error("marginal sea-ice ceiling concentration is outside its declared range")
    all(
        diagnostics.ocean_free_surface_minimum .<=
        diagnostics.ocean_free_surface_maximum,
    ) || error("free-surface minimum exceeds maximum")
    all(diagnostics.ocean_max_abs_zonal_velocity_series .>= 0) &&
        all(diagnostics.ocean_max_abs_meridional_velocity_series .>= 0) ||
        error("absolute ocean-current series is negative")
    return true
end

function _write_dynamic_netcdf(path, diagnostics)
    NCDataset(path, "c") do dataset
        time_origin = replace(diagnostics.metadata.start_date, "T" => " ")
        time_units = "days since $time_origin"
        curvilinear_ocean = ndims(diagnostics.ocean_longitude) == 2
        ocean_dimensions = if curvilinear_ocean
            defDim(dataset, "ocean_x", size(diagnostics.ocean_surface_temperature, 1))
            defDim(dataset, "ocean_y", size(diagnostics.ocean_surface_temperature, 2))
            ("ocean_x", "ocean_y")
        else
            defDim(dataset, "ocean_longitude", length(diagnostics.ocean_longitude))
            defDim(dataset, "ocean_latitude", length(diagnostics.ocean_latitude))
            ("ocean_longitude", "ocean_latitude")
        end
        defDim(dataset, "ocean_layer", length(diagnostics.ocean_layer_depth))
        defDim(dataset, "atmosphere_point", length(diagnostics.atmosphere_longitude))
        defDim(dataset, "atmosphere_layer", length(diagnostics.atmosphere_sigma))
        has_terrarium = !isempty(diagnostics.terrarium_land_atmosphere_point)
        if has_terrarium
            defDim(
                dataset,
                "terrarium_land_column",
                length(diagnostics.terrarium_land_atmosphere_point),
            )
            defDim(
                dataset,
                "terrarium_soil_layer",
                length(diagnostics.terrarium_soil_layer_depth),
            )
        end
        defDim(dataset, "time", length(diagnostics.time_days))
        defDim(
            dataset,
            "atmosphere_diagnostic_time",
            length(diagnostics.atmosphere_diagnostic_time_days),
        )
        defDim(
            dataset,
            "coupled_budget_time",
            length(diagnostics.coupled_budget_time_days),
        )

        ocean_coordinate_dimensions = curvilinear_ocean ?
            ocean_dimensions : (("ocean_longitude",), ("ocean_latitude",))
        ocean_coordinates = curvilinear_ocean ? (
            ("ocean_longitude", diagnostics.ocean_longitude, ocean_dimensions, "degrees_east"),
            ("ocean_latitude", diagnostics.ocean_latitude, ocean_dimensions, "degrees_north"),
        ) : (
            ("ocean_longitude", diagnostics.ocean_longitude, ocean_coordinate_dimensions[1], "degrees_east"),
            ("ocean_latitude", diagnostics.ocean_latitude, ocean_coordinate_dimensions[2], "degrees_north"),
        )
        coordinates = (
            ocean_coordinates...,
            (
                "ocean_layer_depth",
                diagnostics.ocean_layer_depth,
                ("ocean_layer",),
                "m",
            ),
            ("atmosphere_longitude", diagnostics.atmosphere_longitude, ("atmosphere_point",), "degrees_east"),
            ("atmosphere_latitude", diagnostics.atmosphere_latitude, ("atmosphere_point",), "degrees_north"),
            ("atmosphere_sigma", diagnostics.atmosphere_sigma, ("atmosphere_layer",), "1"),
            (
                "atmosphere_sigma_layer_thickness",
                diagnostics.atmosphere_sigma_layer_thickness,
                ("atmosphere_layer",),
                "1",
            ),
            ("time", diagnostics.time_days, ("time",), time_units),
            (
                "atmosphere_diagnostic_time",
                diagnostics.atmosphere_diagnostic_time_days,
                ("atmosphere_diagnostic_time",),
                time_units,
            ),
            (
                "coupled_budget_time",
                diagnostics.coupled_budget_time_days,
                ("coupled_budget_time",),
                time_units,
            ),
        )
        for (name, values, dimensions, units) in coordinates
            variable = defVar(dataset, name, Float64, dimensions)
            variable.attrib["units"] = units
            if name in (
                "time",
                "atmosphere_diagnostic_time",
                "coupled_budget_time",
            )
                variable.attrib["calendar"] = diagnostics.metadata.calendar
            end
            variable[:] = values
        end
        dataset["ocean_layer_depth"].attrib["positive"] = "up"
        if has_terrarium
            land_point = defVar(
                dataset,
                "terrarium_land_atmosphere_point",
                Int64,
                ("terrarium_land_column",),
            )
            land_point.attrib["long_name"] =
                "one-based index into the atmosphere_point coordinate"
            land_point[:] = diagnostics.terrarium_land_atmosphere_point
            soil_depth = defVar(
                dataset,
                "terrarium_soil_layer_depth",
                Float64,
                ("terrarium_soil_layer",),
            )
            soil_depth.attrib["units"] = "m"
            soil_depth.attrib["positive"] = "up"
            soil_depth[:] = diagnostics.terrarium_soil_layer_depth
        end

        ocean_variables = (
            ("ocean_surface_active_mask", diagnostics.ocean_surface_active_mask, "1"),
            ("ocean_cell_area", diagnostics.ocean_cell_area, "m^2"),
            ("ocean_runoff_receiver_mask", diagnostics.ocean_runoff_receiver_mask, "1"),
            ("ocean_routed_land_runoff_flux", diagnostics.ocean_routed_land_runoff_flux, "kg/m^2/s"),
            ("ocean_river_mouth_receiver_mask", diagnostics.ocean_river_mouth_receiver_mask, "1"),
            ("ocean_river_mouth_horizontal_mixing_mask", diagnostics.ocean_river_mouth_horizontal_mixing_mask, "1"),
            ("ocean_river_mouth_active_mixing_mask", diagnostics.ocean_river_mouth_active_mixing_mask, "1"),
            ("ocean_surface_temperature", diagnostics.ocean_surface_temperature, "degree_Celsius"),
            ("ocean_surface_salinity", diagnostics.ocean_surface_salinity, "1e-3"),
            ("ocean_surface_zonal_velocity", diagnostics.ocean_surface_zonal_velocity, "m/s"),
            ("ocean_surface_meridional_velocity", diagnostics.ocean_surface_meridional_velocity, "m/s"),
            ("sea_ice_concentration", diagnostics.sea_ice_concentration, "1"),
            ("sea_ice_thickness", diagnostics.sea_ice_thickness, "m"),
            ("sea_ice_zonal_velocity", diagnostics.sea_ice_zonal_velocity, "m/s"),
            ("sea_ice_meridional_velocity", diagnostics.sea_ice_meridional_velocity, "m/s"),
            ("atmosphere_ocean_sensible_heat_flux", diagnostics.atmosphere_ocean_sensible_heat_flux, "W/m^2"),
            ("atmosphere_ocean_latent_heat_flux", diagnostics.atmosphere_ocean_latent_heat_flux, "W/m^2"),
            ("atmosphere_sea_ice_sensible_heat_flux", diagnostics.atmosphere_sea_ice_sensible_heat_flux, "W/m^2"),
            ("atmosphere_sea_ice_latent_heat_flux", diagnostics.atmosphere_sea_ice_latent_heat_flux, "W/m^2"),
            ("ocean_upwelling_longwave", diagnostics.ocean_upwelling_longwave, "W/m^2"),
            ("ocean_absorbed_longwave", diagnostics.ocean_absorbed_longwave, "W/m^2"),
            ("ocean_absorbed_shortwave", diagnostics.ocean_absorbed_shortwave, "W/m^2"),
        )
        for (name, values, units) in ocean_variables
            variable = defVar(
                dataset,
                name,
                Float64,
                ocean_dimensions,
            )
            variable.attrib["units"] = units
            variable[:, :] = values
        end

        atmosphere_variables = (
            ("atmosphere_land_fraction", diagnostics.atmosphere_land_fraction, "1"),
            ("land_surface_temperature", diagnostics.land_surface_temperature, "K"),
            ("land_surface_moisture", diagnostics.land_surface_moisture, "1"),
            ("atmosphere_surface_temperature", diagnostics.atmosphere_surface_temperature, "K"),
            (
                "atmosphere_radiative_surface_temperature",
                diagnostics.atmosphere_radiative_surface_temperature,
                "K",
            ),
            ("atmosphere_surface_pressure", diagnostics.atmosphere_surface_pressure, "Pa"),
            ("atmosphere_surface_specific_humidity", diagnostics.atmosphere_surface_specific_humidity, "kg/kg"),
            ("atmosphere_surface_zonal_wind", diagnostics.atmosphere_surface_zonal_wind, "m/s"),
            ("atmosphere_surface_meridional_wind", diagnostics.atmosphere_surface_meridional_wind, "m/s"),
            ("atmosphere_column_water_vapor", diagnostics.atmosphere_column_water_vapor, "kg/m^2"),
            ("atmosphere_rainfall_flux", diagnostics.atmosphere_rainfall_flux, "kg/m^2/s"),
            (
                "atmosphere_convective_rainfall_flux",
                diagnostics.atmosphere_convective_rainfall_flux,
                "kg/m^2/s",
            ),
            (
                "atmosphere_convective_precipitation_mass_correction_flux",
                diagnostics.atmosphere_convective_precipitation_mass_correction_flux,
                "kg/m^2/s",
            ),
            (
                "atmosphere_large_scale_rainfall_flux",
                diagnostics.atmosphere_large_scale_rainfall_flux,
                "kg/m^2/s",
            ),
            (
                "atmosphere_large_scale_snowfall_flux",
                diagnostics.atmosphere_large_scale_snowfall_flux,
                "kg/m^2/s",
            ),
            (
                "atmosphere_large_scale_precipitation_mass_correction_flux",
                diagnostics.atmosphere_large_scale_precipitation_mass_correction_flux,
                "kg/m^2/s",
            ),
            ("atmosphere_snowfall_flux", diagnostics.atmosphere_snowfall_flux, "kg/m^2/s"),
            (
                "atmosphere_precipitation_cloud_top_layer",
                diagnostics.atmosphere_precipitation_cloud_top_layer,
                "1",
            ),
            (
                "atmosphere_vertical_diffusion_cfl_scale",
                diagnostics.atmosphere_vertical_diffusion_cfl_scale,
                "1",
            ),
            ("atmosphere_surface_water_vapor_flux", diagnostics.atmosphere_surface_water_vapor_flux, "kg/m^2/s"),
            ("atmosphere_surface_sensible_heat_flux", diagnostics.atmosphere_surface_sensible_heat_flux, "W/m^2"),
            ("outgoing_longwave", diagnostics.outgoing_longwave, "W/m^2"),
            ("column_cloud_fraction", diagnostics.column_cloud_fraction, "1"),
        )
        for (name, values, units) in atmosphere_variables
            variable = defVar(dataset, name, Float64, ("atmosphere_point",))
            variable.attrib["units"] = units
            variable[:] = values
        end
        atmosphere_3d_variables = (
            ("atmosphere_temperature", diagnostics.atmosphere_temperature, "K"),
            ("atmosphere_specific_humidity", diagnostics.atmosphere_specific_humidity, "kg/kg"),
            ("atmosphere_pressure", diagnostics.atmosphere_pressure, "Pa"),
            ("atmosphere_zonal_wind", diagnostics.atmosphere_zonal_wind, "m/s"),
            ("atmosphere_meridional_wind", diagnostics.atmosphere_meridional_wind, "m/s"),
            ("atmosphere_cloud_fraction", diagnostics.atmosphere_cloud_fraction, "1"),
            (
                "atmosphere_cloud_liquid_water_path",
                diagnostics.atmosphere_cloud_liquid_water_path,
                "g/m^2",
            ),
            (
                "atmosphere_cloud_ice_water_path",
                diagnostics.atmosphere_cloud_ice_water_path,
                "g/m^2",
            ),
            (
                "atmosphere_convective_humidity_tendency",
                diagnostics.atmosphere_convective_humidity_tendency,
                "kg/kg/s",
            ),
            (
                "atmosphere_convective_temperature_tendency",
                diagnostics.atmosphere_convective_temperature_tendency,
                "K/s",
            ),
            (
                "atmosphere_large_scale_condensation_humidity_tendency",
                diagnostics.atmosphere_large_scale_condensation_humidity_tendency,
                "kg/kg/s",
            ),
            (
                "atmosphere_large_scale_condensation_temperature_tendency",
                diagnostics.atmosphere_large_scale_condensation_temperature_tendency,
                "K/s",
            ),
        )
        for (name, values, units) in atmosphere_3d_variables
            variable = defVar(
                dataset,
                name,
                Float64,
                ("atmosphere_point", "atmosphere_layer"),
            )
            variable.attrib["units"] = units
            variable[:, :] = values
        end
        atmosphere_process_profiles = (
            (
                "cumulative_global_convective_humidity_change_profile",
                diagnostics.cumulative_global_convective_humidity_change_profile,
                "kg/kg",
            ),
            (
                "cumulative_global_convective_temperature_change_profile",
                diagnostics.cumulative_global_convective_temperature_change_profile,
                "K",
            ),
            (
                "cumulative_global_large_scale_humidity_change_profile",
                diagnostics.cumulative_global_large_scale_humidity_change_profile,
                "kg/kg",
            ),
            (
                "cumulative_global_large_scale_temperature_change_profile",
                diagnostics.cumulative_global_large_scale_temperature_change_profile,
                "K",
            ),
            (
                "cumulative_global_radiative_temperature_change_profile",
                diagnostics.cumulative_global_radiative_temperature_change_profile,
                "K",
            ),
            (
                "cumulative_global_surface_humidity_change_profile",
                diagnostics.cumulative_global_surface_humidity_change_profile,
                "kg/kg",
            ),
            (
                "cumulative_global_surface_sensible_temperature_change_profile",
                diagnostics.cumulative_global_surface_sensible_temperature_change_profile,
                "K",
            ),
        )
        for (name, values, units) in atmosphere_process_profiles
            variable = defVar(
                dataset,
                name,
                Float64,
                ("atmosphere_layer", "atmosphere_diagnostic_time"),
            )
            variable.attrib["units"] = units
            variable.attrib["long_name"] =
                "cumulative Gaussian-area mean applied process contribution"
            variable[:, :] = values
        end
        if has_terrarium
            terrarium_3d_variables = (
                (
                    "terrarium_soil_temperature",
                    diagnostics.terrarium_soil_temperature,
                    "degree_Celsius",
                ),
                (
                    "terrarium_soil_saturation",
                    diagnostics.terrarium_soil_saturation,
                    "1",
                ),
            )
            if !isempty(diagnostics.terrarium_root_fraction)
                terrarium_3d_variables = (
                    terrarium_3d_variables...,
                    (
                        "terrarium_root_fraction",
                        diagnostics.terrarium_root_fraction,
                        "1",
                    ),
                    (
                        "terrarium_root_water_availability",
                        diagnostics.terrarium_root_water_availability,
                        "1",
                    ),
                )
            end
            for (name, values, units) in terrarium_3d_variables
                variable = defVar(
                    dataset,
                    name,
                    Float64,
                    ("terrarium_land_column", "terrarium_soil_layer"),
                )
                variable.attrib["units"] = units
                variable[:, :] = values
            end
            terrarium_surface_variables = (
                (
                    "terrarium_skin_temperature",
                    diagnostics.terrarium_skin_temperature,
                    "degree_Celsius",
                ),
                (
                    "terrarium_surface_excess_water",
                    diagnostics.terrarium_surface_excess_water,
                    "m",
                ),
                (
                    "terrarium_surface_runoff",
                    diagnostics.terrarium_surface_runoff,
                    "m/s",
                ),
                (
                    "terrarium_infiltration",
                    diagnostics.terrarium_infiltration,
                    "m/s",
                ),
            )
            if !isempty(diagnostics.terrarium_vegetation_fraction)
                terrarium_surface_variables = (
                    terrarium_surface_variables...,
                    (
                        "terrarium_vegetation_fraction",
                        diagnostics.terrarium_vegetation_fraction,
                        "1",
                    ),
                    (
                        "terrarium_leaf_area_index",
                        diagnostics.terrarium_leaf_area_index,
                        "1",
                    ),
                    (
                        "terrarium_root_uptake_weight_sum",
                        diagnostics.terrarium_root_uptake_weight_sum,
                        "1",
                    ),
                    (
                        "terrarium_ground_evaporation_flux",
                        diagnostics.terrarium_ground_evaporation_flux,
                        "kg/m^2/s",
                    ),
                    (
                        "terrarium_transpiration_flux",
                        diagnostics.terrarium_transpiration_flux,
                        "kg/m^2/s",
                    ),
                )
            end
            for (name, values, units) in terrarium_surface_variables
                variable = defVar(
                    dataset,
                    name,
                    Float64,
                    ("terrarium_land_column",),
                )
                variable.attrib["units"] = units
                variable[:] = values
            end
            terrarium_profile_series = (
                (
                    "global_land_soil_layer_temperature",
                    diagnostics.global_land_soil_layer_temperature,
                    "degree_Celsius",
                ),
                (
                    "global_land_soil_layer_saturation",
                    diagnostics.global_land_soil_layer_saturation,
                    "1",
                ),
            )
            for (name, values, units) in terrarium_profile_series
                variable = defVar(
                    dataset,
                    name,
                    Float64,
                    ("terrarium_soil_layer", "atmosphere_diagnostic_time"),
                )
                variable.attrib["units"] = units
                variable[:, :] = values
            end
        end
        global_temperature = defVar(dataset, "global_surface_temperature", Float64, ("time",))
        global_temperature.attrib["units"] = "K"
        global_temperature[:] = diagnostics.global_surface_temperature
        budget_variables = (
            ("toa_incoming_shortwave", diagnostics.toa_incoming_shortwave),
            ("toa_outgoing_shortwave", diagnostics.toa_outgoing_shortwave),
            ("toa_outgoing_longwave", diagnostics.toa_outgoing_longwave),
            ("toa_clear_outgoing_shortwave", diagnostics.toa_clear_outgoing_shortwave),
            ("toa_clear_outgoing_longwave", diagnostics.toa_clear_outgoing_longwave),
            ("toa_net_downward", diagnostics.toa_net_downward),
            ("toa_clear_net_downward", diagnostics.toa_clear_net_downward),
        )
        for (name, values) in budget_variables
            variable = defVar(dataset, name, Float64, ("time",))
            variable.attrib["units"] = "W/m^2"
            variable[:] = values
        end
        atmosphere_series = (
            ("global_rainfall_flux", diagnostics.global_rainfall_flux, "kg/m^2/s"),
            (
                "global_convective_rainfall_flux",
                diagnostics.global_convective_rainfall_flux,
                "kg/m^2/s",
            ),
            (
                "global_large_scale_rainfall_flux",
                diagnostics.global_large_scale_rainfall_flux,
                "kg/m^2/s",
            ),
            (
                "global_large_scale_snowfall_flux",
                diagnostics.global_large_scale_snowfall_flux,
                "kg/m^2/s",
            ),
            ("global_snowfall_flux", diagnostics.global_snowfall_flux, "kg/m^2/s"),
            ("global_surface_water_vapor_flux", diagnostics.global_surface_water_vapor_flux, "kg/m^2/s"),
            ("cumulative_global_surface_water_vapor", diagnostics.cumulative_global_surface_water_vapor, "kg/m^2"),
            ("cumulative_global_precipitation", diagnostics.cumulative_global_precipitation, "kg/m^2"),
            ("cumulative_balanced_launch_atmosphere_water_source", diagnostics.cumulative_balanced_launch_atmosphere_water_source, "kg/m^2 global"),
            ("global_atmosphere_water_budget_residual", diagnostics.global_atmosphere_water_budget_residual, "kg/m^2"),
            ("global_cloud_condensate_water_path", diagnostics.global_cloud_condensate_water_path, "kg/m^2"),
            ("global_total_atmosphere_water_budget_residual", diagnostics.global_total_atmosphere_water_budget_residual, "kg/m^2"),
            ("global_surface_air_temperature", diagnostics.global_surface_air_temperature, "K"),
            ("global_surface_specific_humidity", diagnostics.global_surface_specific_humidity, "kg/kg"),
            ("global_surface_pressure", diagnostics.global_surface_pressure, "Pa"),
            ("global_near_surface_wind_speed", diagnostics.global_near_surface_wind_speed, "m/s"),
            ("global_column_water_vapor", diagnostics.global_column_water_vapor, "kg/m^2"),
            ("global_mass_weighted_atmosphere_temperature", diagnostics.global_mass_weighted_atmosphere_temperature, "K"),
            ("global_land_surface_soil_moisture", diagnostics.global_land_surface_soil_moisture, "1"),
            ("global_land_total_water_storage", diagnostics.global_land_total_water_storage, "kg/m^2 of land"),
            ("cumulative_global_land_precipitation", diagnostics.cumulative_global_land_precipitation, "kg/m^2 of land"),
            ("cumulative_global_land_evapotranspiration", diagnostics.cumulative_global_land_evapotranspiration, "kg/m^2 of land"),
            ("cumulative_global_land_surface_runoff", diagnostics.cumulative_global_land_surface_runoff, "kg/m^2 of land"),
            ("cumulative_balanced_launch_land_water_source", diagnostics.cumulative_balanced_launch_land_water_source, "kg/m^2 of land"),
            ("global_land_water_budget_residual", diagnostics.global_land_water_budget_residual, "kg/m^2 of land"),
            ("global_land_rainfall_flux", diagnostics.global_land_rainfall_flux, "kg/m^2/s"),
            ("global_land_snowfall_flux", diagnostics.global_land_snowfall_flux, "kg/m^2/s"),
            ("global_land_evaporation_flux", diagnostics.global_land_evaporation_flux, "kg/m^2/s"),
            ("global_land_ground_evaporation_flux", diagnostics.global_land_ground_evaporation_flux, "kg/m^2/s"),
            ("global_land_transpiration_flux", diagnostics.global_land_transpiration_flux, "kg/m^2/s"),
            ("global_land_surface_runoff_flux", diagnostics.global_land_surface_runoff_flux, "kg/m^2/s"),
            ("global_land_infiltration_flux", diagnostics.global_land_infiltration_flux, "kg/m^2/s"),
        )
        for (name, values, units) in atmosphere_series
            variable = defVar(
                dataset,
                name,
                Float64,
                ("atmosphere_diagnostic_time",),
            )
            variable.attrib["units"] = units
            variable[:] = values
        end
        balanced_fraction = defVar(
            dataset,
            "balanced_launch_applied_fraction",
            Float64,
            ("coupled_budget_time",),
        )
        balanced_fraction.attrib["units"] = "1"
        balanced_fraction[:] = diagnostics.balanced_launch_applied_fraction
        for (source_index, spec) in enumerate(
            _BALANCED_LAUNCH_DIAGNOSTIC_SOURCE_SPECS,
        )
            variable = defVar(
                dataset,
                spec.netcdf_name,
                Float64,
                ("coupled_budget_time",),
            )
            variable.attrib["units"] = spec.units
            variable[:] = view(
                diagnostics.balanced_launch_cumulative_sources,
                source_index,
                :,
            )
        end
        coupled_budget_series = (
            ("ocean_heat_content", diagnostics.ocean_heat_content, "J relative to 0 degree_Celsius"),
            ("ocean_salt_content", diagnostics.ocean_salt_content, "kg approximate practical-salinity mass"),
            ("ocean_volume", diagnostics.ocean_volume, "m^3"),
            ("ocean_temperature_minimum_series", diagnostics.ocean_temperature_minimum_series, "degree_Celsius"),
            ("ocean_temperature_maximum_series", diagnostics.ocean_temperature_maximum_series, "degree_Celsius"),
            ("ocean_salinity_minimum_series", diagnostics.ocean_salinity_minimum_series, "1e-3"),
            ("ocean_salinity_maximum_series", diagnostics.ocean_salinity_maximum_series, "1e-3"),
            ("ocean_surface_salinity_minimum_series", diagnostics.ocean_surface_salinity_minimum_series, "1e-3"),
            ("ocean_surface_salinity_minimum_longitude", diagnostics.ocean_surface_salinity_minimum_longitude, "degrees_east"),
            ("ocean_surface_salinity_minimum_latitude", diagnostics.ocean_surface_salinity_minimum_latitude, "degrees_north"),
            ("ocean_runoff_receiver_surface_salinity_minimum_series", diagnostics.ocean_runoff_receiver_surface_salinity_minimum_series, "1e-3"),
            ("ocean_runoff_receiver_surface_salinity_minimum_longitude", diagnostics.ocean_runoff_receiver_surface_salinity_minimum_longitude, "degrees_east"),
            ("ocean_runoff_receiver_surface_salinity_minimum_latitude", diagnostics.ocean_runoff_receiver_surface_salinity_minimum_latitude, "degrees_north"),
            ("ocean_river_mouth_active_mixing_area_m2", diagnostics.ocean_river_mouth_active_mixing_area_m2, "m^2"),
            ("sea_ice_volume", diagnostics.sea_ice_volume, "m^3"),
            ("sea_ice_area", diagnostics.sea_ice_area, "m^2"),
            ("tropical_sea_ice_area", diagnostics.tropical_sea_ice_area, "m^2"),
            ("sea_ice_max_thickness", diagnostics.sea_ice_max_thickness, "m"),
            ("sea_ice_extent_max_thickness", diagnostics.sea_ice_extent_max_thickness, "m"),
            ("sea_ice_max_area_equivalent_thickness", diagnostics.sea_ice_max_area_equivalent_thickness, "m"),
            ("sea_ice_marginal_ceiling_max_concentration", diagnostics.sea_ice_marginal_ceiling_max_concentration, "1"),
            ("sea_ice_marginal_ceiling_max_area_equivalent_thickness", diagnostics.sea_ice_marginal_ceiling_max_area_equivalent_thickness, "m"),
            ("ocean_free_surface_minimum", diagnostics.ocean_free_surface_minimum, "m"),
            ("ocean_free_surface_maximum", diagnostics.ocean_free_surface_maximum, "m"),
            ("ocean_max_abs_zonal_velocity_series", diagnostics.ocean_max_abs_zonal_velocity_series, "m/s"),
            ("ocean_max_abs_meridional_velocity_series", diagnostics.ocean_max_abs_meridional_velocity_series, "m/s"),
            ("cumulative_ocean_surface_heat", diagnostics.cumulative_ocean_surface_heat, "J positive upward"),
            ("cumulative_ocean_frazil_heat", diagnostics.cumulative_ocean_frazil_heat, "J positive upward"),
            ("cumulative_ocean_penetrating_shortwave_heat", diagnostics.cumulative_ocean_penetrating_shortwave_heat, "J positive downward into ocean"),
            ("cumulative_net_ocean_heat", diagnostics.cumulative_net_ocean_heat, "J positive upward"),
            ("cumulative_ocean_freshwater_enthalpy", diagnostics.cumulative_ocean_freshwater_enthalpy, "J positive into ocean"),
            ("cumulative_net_ocean_freshwater_volume", diagnostics.cumulative_net_ocean_freshwater_volume, "m^3"),
            ("cumulative_net_ocean_salinity", diagnostics.cumulative_net_ocean_salinity, "1e-3 m^3"),
            ("cumulative_rain_mass", diagnostics.cumulative_rain_mass, "kg"),
            ("cumulative_snow_mass", diagnostics.cumulative_snow_mass, "kg"),
            ("cumulative_open_ocean_evaporation_mass", diagnostics.cumulative_open_ocean_evaporation_mass, "kg"),
            ("cumulative_sea_ice_ocean_freshwater_volume", diagnostics.cumulative_sea_ice_ocean_freshwater_volume, "m^3"),
            ("cumulative_routed_land_runoff_mass", diagnostics.cumulative_routed_land_runoff_mass, "kg"),
            ("ocean_heat_closure_residual", diagnostics.ocean_heat_closure_residual, "J"),
            ("ocean_heat_closure_relative", diagnostics.ocean_heat_closure_relative, "1"),
            ("ocean_freshwater_closure_residual", diagnostics.ocean_freshwater_closure_residual, "m^3"),
            ("ocean_freshwater_closure_relative", diagnostics.ocean_freshwater_closure_relative, "1"),
            ("ocean_salt_closure_residual", diagnostics.ocean_salt_closure_residual, "kg approximate practical-salinity mass"),
            ("ocean_salt_closure_relative", diagnostics.ocean_salt_closure_relative, "1"),
        )
        for (name, values, units) in coupled_budget_series
            variable = defVar(dataset, name, Float64, ("coupled_budget_time",))
            variable.attrib["units"] = units
            variable[:] = values
        end
        ocean_profile_series = (
            (
                "ocean_mean_temperature_profile",
                diagnostics.ocean_mean_temperature_profile,
                "degree_Celsius",
            ),
            (
                "ocean_mean_salinity_profile",
                diagnostics.ocean_mean_salinity_profile,
                "1e-3",
            ),
        )
        for (name, values, units) in ocean_profile_series
            variable = defVar(
                dataset,
                name,
                Float64,
                ("ocean_layer", "coupled_budget_time"),
            )
            variable.attrib["units"] = units
            variable[:, :] = values
        end
        marginal_ceiling_count = defVar(
            dataset,
            "sea_ice_marginal_ceiling_cell_count",
            Int64,
            ("coupled_budget_time",),
        )
        marginal_ceiling_count.attrib["long_name"] =
            "marginal cells at the conditional-thickness representation ceiling"
        marginal_ceiling_count[:] = diagnostics.sea_ice_marginal_ceiling_cell_count
        active_mixing_count = defVar(
            dataset,
            "ocean_river_mouth_active_mixing_cells",
            Int64,
            ("coupled_budget_time",),
        )
        active_mixing_count.attrib["long_name"] =
            "routed receiver cells with positive freshwater flux and active localized mixing"
        active_mixing_count[:] = diagnostics.ocean_river_mouth_active_mixing_cells
        for (name, values) in (
            ("ocean_surface_salinity_minimum_i", diagnostics.ocean_surface_salinity_minimum_i),
            ("ocean_surface_salinity_minimum_j", diagnostics.ocean_surface_salinity_minimum_j),
            ("ocean_runoff_receiver_surface_salinity_minimum_i", diagnostics.ocean_runoff_receiver_surface_salinity_minimum_i),
            ("ocean_runoff_receiver_surface_salinity_minimum_j", diagnostics.ocean_runoff_receiver_surface_salinity_minimum_j),
        )
            variable = defVar(
                dataset,
                name,
                Int64,
                ("coupled_budget_time",),
            )
            variable.attrib["long_name"] = startswith(
                name,
                "ocean_runoff_receiver",
            ) ?
                "one-based native-grid index of the daily runoff-receiver surface-salinity minimum" :
                "one-based native-grid index of the daily surface-salinity minimum"
            variable[:] = values
        end

        circulation_scalars = (
            ("ocean_max_abs_zonal_velocity", diagnostics.ocean_max_abs_zonal_velocity),
            ("ocean_max_abs_meridional_velocity", diagnostics.ocean_max_abs_meridional_velocity),
            ("ocean_max_abs_vertical_velocity", diagnostics.ocean_max_abs_vertical_velocity),
            ("ocean_mean_kinetic_energy", diagnostics.ocean_mean_kinetic_energy),
        )
        for (name, value) in circulation_scalars
            dataset.attrib[name] = value
        end

        dataset.attrib["ocean_temperature_minimum"] =
            diagnostics.ocean_temperature_minimum
        dataset.attrib["ocean_temperature_maximum"] =
            diagnostics.ocean_temperature_maximum
        dataset.attrib["ocean_salinity_minimum"] =
            diagnostics.ocean_salinity_minimum
        dataset.attrib["ocean_salinity_maximum"] =
            diagnostics.ocean_salinity_maximum
        dataset.attrib["ocean_full_depth_extrema_domain"] =
            "active ocean cells only; immersed/dry fill values excluded"
        dataset.attrib["ocean_surface_tracer_dry_fill_value"] = 0.0

        for (key, value) in pairs(diagnostics.metadata)
            # NetCDF attributes have no native boolean type.
            dataset.attrib[String(key)] = _netcdf_attribute_value(value)
        end
        dataset.attrib["created_at_utc"] = string(now(UTC))
    end
    return path
end

"""Write dynamic NetCDF diagnostics through a same-directory atomic rename."""
function _write_dynamic_netcdf_atomic(path, diagnostics)
    directory = dirname(path)
    mkpath(directory)
    isdir(path) && error(
        "refusing to replace a directory with dynamic NetCDF diagnostics: $path",
    )
    temporary_path = joinpath(
        directory,
        ".$(basename(path)).tmp-$(getpid())-$(time_ns())",
    )
    try
        _write_dynamic_netcdf(temporary_path, diagnostics)
        mv(temporary_path, path; force = true)
    finally
        isfile(temporary_path) && rm(temporary_path; force = true)
    end
    return path
end

function _dynamic_heatmap!(figure, position, longitude, latitude, values, title; colorrange = nothing)
    axis = Axis(
        _plot_layout(figure, position...);
        title,
        xlabel = "longitude (°E)",
        ylabel = "latitude (°N)",
    )
    plot = if ndims(longitude) == 2
        isnothing(colorrange) ?
            scatter!(axis, vec(longitude), vec(latitude); color = vec(values), marker = :rect, markersize = 3) :
            scatter!(axis, vec(longitude), vec(latitude); color = vec(values), marker = :rect, markersize = 3, colorrange)
    else
        isnothing(colorrange) ?
            heatmap!(axis, longitude, latitude, values) :
            heatmap!(axis, longitude, latitude, values; colorrange)
    end
    Colorbar(_plot_layout(figure, position[1], position[2] + 1), plot)
    return axis
end

function _dynamic_scatter!(figure, position, longitude, latitude, values, title; colorrange = nothing)
    axis = Axis(
        _plot_layout(figure, position...);
        title,
        xlabel = "longitude (°E)",
        ylabel = "latitude (°N)",
    )
    plot = isnothing(colorrange) ?
        scatter!(axis, longitude, latitude; color = values, markersize = 8) :
        scatter!(axis, longitude, latitude; color = values, markersize = 8, colorrange)
    Colorbar(_plot_layout(figure, position[1], position[2] + 1), plot)
    return axis
end

function _write_dynamic_figure(path, diagnostics)
    figure = Figure(size = (1300, 2450))
    budget_axis = Axis(
        _plot_layout(figure, 1, 1:4);
        title = "Global top-of-atmosphere radiation budget",
        xlabel = "model time (days)",
        ylabel = "flux (W m⁻²)",
    )
    lines!(budget_axis, diagnostics.time_days, diagnostics.toa_incoming_shortwave;
           label = "incoming SW", linewidth = 2)
    lines!(budget_axis, diagnostics.time_days, diagnostics.toa_outgoing_shortwave;
           label = "reflected SW", linewidth = 2)
    lines!(budget_axis, diagnostics.time_days, diagnostics.toa_outgoing_longwave;
           label = "outgoing LW", linewidth = 2)
    lines!(budget_axis, diagnostics.time_days, diagnostics.toa_net_downward;
           label = "net downward", linewidth = 2)
    lines!(budget_axis, diagnostics.time_days, diagnostics.toa_clear_net_downward;
           label = "clear-sky net downward", linewidth = 2, linestyle = :dash)
    axislegend(budget_axis; position = :rb, orientation = :horizontal)
    precipitation_axis = Axis(
        _plot_layout(figure, 2, 1:2);
        title = "Global water fluxes",
        xlabel = "model time (days)",
        ylabel = "water equivalent (mm day⁻¹)",
    )
    seconds_per_day = 86_400
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_rainfall_flux;
           label = "global rain", linewidth = 2)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_convective_rainfall_flux;
           label = "convective rain", linewidth = 1, linestyle = :dash)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_large_scale_rainfall_flux;
           label = "large-scale rain", linewidth = 1, linestyle = :dot)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_snowfall_flux;
           label = "global snow", linewidth = 2)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_surface_water_vapor_flux;
           label = "global surface evaporation", linewidth = 2)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_land_evaporation_flux;
           label = "total land ET", linewidth = 2)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_land_ground_evaporation_flux;
           label = "ground evaporation", linewidth = 2, linestyle = :dash)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_land_transpiration_flux;
           label = "transpiration", linewidth = 2, linestyle = :dot)
    lines!(precipitation_axis, diagnostics.atmosphere_diagnostic_time_days,
           seconds_per_day .* diagnostics.global_land_surface_runoff_flux;
           label = "land runoff", linewidth = 2)
    axislegend(precipitation_axis; position = :rt)
    atmosphere_axis = Axis(
        _plot_layout(figure, 2, 3:4);
        title = "Global atmospheric temperature",
        xlabel = "model time (days)",
        ylabel = "temperature (K)",
    )
    lines!(atmosphere_axis, diagnostics.atmosphere_diagnostic_time_days,
           diagnostics.global_surface_air_temperature;
           label = "surface air", linewidth = 2)
    lines!(atmosphere_axis, diagnostics.atmosphere_diagnostic_time_days,
           diagnostics.global_mass_weighted_atmosphere_temperature;
           label = "mass weighted", linewidth = 2)
    axislegend(atmosphere_axis; position = :rb)
    _dynamic_heatmap!(figure, (3, 1), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_surface_temperature, "Ocean surface temperature (°C)")
    _dynamic_heatmap!(figure, (3, 3), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_surface_salinity, "Ocean surface salinity")
    _dynamic_heatmap!(figure, (4, 1), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.sea_ice_concentration, "Sea-ice concentration"; colorrange = (0, 1))
    sea_ice_extent_thickness = ifelse.(
        diagnostics.sea_ice_concentration .>=
        _SEA_ICE_EXTENT_CONCENTRATION_THRESHOLD,
        diagnostics.sea_ice_thickness,
        NaN,
    )
    _dynamic_heatmap!(figure, (4, 3), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      sea_ice_extent_thickness, "Sea-ice thickness within 15% extent (m)")
    _dynamic_scatter!(figure, (5, 1), diagnostics.atmosphere_longitude,
                      diagnostics.atmosphere_latitude, diagnostics.atmosphere_surface_temperature,
                      "Atmosphere surface temperature (K)")
    _dynamic_scatter!(figure, (5, 3), diagnostics.atmosphere_longitude,
                      diagnostics.atmosphere_latitude, diagnostics.outgoing_longwave,
                      "Outgoing longwave (W m⁻²)")
    _dynamic_heatmap!(figure, (6, 1), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.atmosphere_ocean_sensible_heat_flux,
                      "Atmosphere–ocean sensible heat (W m⁻²)")
    _dynamic_scatter!(figure, (6, 3), diagnostics.atmosphere_longitude,
                      diagnostics.atmosphere_latitude, diagnostics.column_cloud_fraction,
                      "Column cloud fraction"; colorrange = (0, 1))
    _dynamic_scatter!(figure, (7, 1), diagnostics.atmosphere_longitude,
                      diagnostics.atmosphere_latitude, diagnostics.land_surface_temperature,
                      "Land-surface temperature (K)")
    _dynamic_scatter!(figure, (7, 3), diagnostics.atmosphere_longitude,
                      diagnostics.atmosphere_latitude, diagnostics.land_surface_moisture,
                      "Land-surface soil moisture"; colorrange = (0, 1))
    _dynamic_heatmap!(figure, (8, 1), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_surface_zonal_velocity,
                      "Ocean surface zonal velocity (m s⁻¹)")
    _dynamic_heatmap!(figure, (8, 3), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_surface_meridional_velocity,
                      "Ocean surface meridional velocity (m s⁻¹)")
    _dynamic_heatmap!(figure, (9, 1), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_absorbed_shortwave,
                      "Ocean absorbed shortwave (W m⁻²)")
    _dynamic_heatmap!(figure, (9, 3), diagnostics.ocean_longitude, diagnostics.ocean_latitude,
                      diagnostics.ocean_absorbed_longwave - diagnostics.ocean_upwelling_longwave,
                      "Ocean net longwave (W m⁻²)")
    save(path, figure; px_per_unit = 1.5)
    return path
end

"""
Atomically save dynamic ESM diagnostics to NetCDF and optionally render a
summary PNG.

Set `render_figure = false` when a runner must preserve the NetCDF evidence
before applying scientific validation. The returned `figure_path` remains the
canonical destination for a subsequent `_write_dynamic_figure` call.
"""
function save_dynamic_diagnostics(
    diagnostics,
    config::ExperimentConfig;
    render_figure::Bool = true,
)
    mkpath(config.output_dir)
    netcdf_path = joinpath(config.output_dir, "diagnostics.nc")
    figure_path = joinpath(config.output_dir, "summary.png")
    _write_dynamic_netcdf_atomic(netcdf_path, diagnostics)
    render_figure && _write_dynamic_figure(figure_path, diagnostics)
    return (; netcdf_path, figure_path)
end
