"""
NumericalEarth land adapter that routes the runoff already evolved by
SpeedyWeather's embedded Terrarium model into the dynamic ocean.

The adapter owns no second land state and takes no timestep. Each Terrarium
column's depth-rate runoff is multiplied by its fractional Gaussian-cell area
to form discharge and spread across the nearest wet ocean cells. The selected
routing weighting controls whether those receivers get equal areal mass flux
(the legacy behavior) or equal fractional column-volume input. Both choices
retain exact global discharge without injecting a whole T31 land column (or a
cluster of inland columns) into one 1-degree ocean cell.
"""
struct TerrariumRunoffLand{R, I, W, M}
    surface_runoff::R
    contribution_column::I
    contribution_weight::W
    target_i::I
    target_j::I
    offsets::I
    source_columns::Int
    target_cells::Int
    receivers_per_source::Int
    routing_weighting::Symbol
    represented_land_area::Float64
    river_mouth_mixing::M
end

const DYNAMIC_RIVER_MOUTH_UPDATE_COMPLETION_PROVENANCE =
    "device_completion_after_runoff_scatter_and_coefficient_halo_fill_v1"

# Preserve the original constructor for non-localized tests and callers. The
# derived mixing fields are attached only by the coupled builder.
TerrariumRunoffLand(
    surface_runoff,
    contribution_column,
    contribution_weight,
    target_i,
    target_j,
    offsets,
    source_columns,
    target_cells,
    receivers_per_source,
    routing_weighting,
    represented_land_area,
) = TerrariumRunoffLand(
    surface_runoff,
    contribution_column,
    contribution_weight,
    target_i,
    target_j,
    offsets,
    source_columns,
    target_cells,
    receivers_per_source,
    routing_weighting,
    represented_land_area,
    nothing,
)

Base.summary(land::TerrariumRunoffLand) =
    "Terrarium runoff router ($(land.source_columns) columns x " *
    "$(land.receivers_per_source) receivers -> $(land.target_cells) ocean cells, " *
    "$(land.routing_weighting))"

Oceananigans.TimeSteppers.time_step!(::TerrariumRunoffLand, Δt) = nothing
Oceananigans.Simulations.reset_clock!(::TerrariumRunoffLand) = nothing
Oceananigans.prognostic_state(::TerrariumRunoffLand) = nothing
Oceananigans.restore_prognostic_state!(land::TerrariumRunoffLand, ::Nothing) = land
NumericalEarth.EarthSystemModels.adopt_clock(land::TerrariumRunoffLand, clock) = land
NumericalEarth.EarthSystemModels.update_net_fluxes!(coupled_model, ::TerrariumRunoffLand) =
    nothing
NumericalEarth.EarthSystemModels.InterfaceComputations.atmosphere_land_interface(
    grid,
    atmosphere,
    ::TerrariumRunoffLand;
    kw...,
) = nothing

function NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
    ::TerrariumRunoffLand,
    grid,
)
    state = (freshwater_flux = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid),)
    return NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
        state,
        nothing,
    )
end

@inline function _wrapped_spherical_distance_squared(λ₁, φ₁, λ₂, φ₂)
    Δλ = mod(λ₂ - λ₁ + 180, 360) - 180
    Δλ *= cosd((φ₁ + φ₂) / 2)
    Δφ = φ₂ - φ₁
    return Δλ^2 + Δφ^2
end

function _wet_ocean_cells(grid; compute_volume = true)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    λc = Float64.(collect(Oceananigans.Grids.λnodes(
        cpu_grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )))
    φc = Float64.(collect(Oceananigans.Grids.φnodes(
        cpu_grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )))
    wet_i = Int[]
    wet_j = Int[]
    wet_area = Float64[]
    wet_volume = Float64[]
    for j in 1:Ny, i in 1:Nx
        inactive = Oceananigans.Grids.inactive_node(
            i,
            j,
            Nz,
            cpu_grid,
            Oceananigans.Center(),
            Oceananigans.Center(),
            Oceananigans.Center(),
        )
        inactive && continue
        area = Float64(Oceananigans.Operators.Azᶜᶜᶜ(
            i,
            j,
            Nz,
            cpu_grid,
        ))
        volume = NaN
        if compute_volume
            volume = 0.0
            for k in 1:Nz
                Oceananigans.Grids.inactive_node(
                    i,
                    j,
                    k,
                    cpu_grid,
                    Oceananigans.Center(),
                    Oceananigans.Center(),
                    Oceananigans.Center(),
                ) && continue
                volume += Float64(Oceananigans.volume(
                    i,
                    j,
                    k,
                    cpu_grid,
                    Oceananigans.Center(),
                    Oceananigans.Center(),
                    Oceananigans.Center(),
                ))
            end
        end
        area > 0 || error("Terrarium runoff routing found non-positive wet-cell area")
        !compute_volume || volume > 0 || error(
            "Terrarium runoff routing found non-positive wet-column volume",
        )
        push!(wet_i, i)
        push!(wet_j, j)
        push!(wet_area, area)
        push!(wet_volume, volume)
    end
    isempty(wet_i) && error("Terrarium runoff routing found no wet ocean cells")
    return (; λc, φc, wet_i, wet_j, wet_area, wet_volume)
end

@inline function _horizontal_coordinate(nodes, i, j, Nx, Ny)
    ndims(nodes) == 2 && return nodes[i, j]
    length(nodes) == Nx && return nodes[i]
    length(nodes) == Ny && return nodes[j]
    throw(DimensionMismatch("horizontal coordinate shape is incompatible with the ocean grid"))
end

function _build_runoff_routing(
    source_longitude,
    source_latitude,
    source_area,
    ocean_grid;
    freshwater_density = 1000,
    receivers_per_source = 16,
    routing_weighting = :equal_area_flux,
)
    length(source_longitude) == length(source_latitude) == length(source_area) ||
        throw(DimensionMismatch("runoff source coordinates and areas must have equal length"))
    receivers_per_source > 0 || throw(
        ArgumentError("receivers_per_source must be positive"),
    )
    routing_weighting in (:equal_area_flux, :equal_column_fraction) || throw(
        ArgumentError(
            "routing_weighting must be equal_area_flux or " *
            "equal_column_fraction",
        ),
    )
    all(isfinite, source_longitude) || throw(ArgumentError(
        "runoff source longitudes must be finite",
    ))
    all(isfinite, source_latitude) || throw(ArgumentError(
        "runoff source latitudes must be finite",
    ))
    all(area -> isfinite(area) && area > 0, source_area) || throw(
        ArgumentError("runoff source areas must be finite and positive"),
    )
    isfinite(freshwater_density) && freshwater_density > 0 || throw(
        ArgumentError("freshwater_density must be finite and positive"),
    )
    wet = _wet_ocean_cells(
        ocean_grid;
        compute_volume = routing_weighting == :equal_column_fraction,
    )
    nreceivers = min(receivers_per_source, length(wet.wet_i))
    contributions = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}()
    distances = Vector{Float64}(undef, length(wet.wet_i))

    for column in eachindex(source_longitude)
        @inbounds for n in eachindex(wet.wet_i)
            i = wet.wet_i[n]
            j = wet.wet_j[n]
            distances[n] = _wrapped_spherical_distance_squared(
                source_longitude[column],
                source_latitude[column],
                _horizontal_coordinate(wet.λc, i, j, size(ocean_grid, 1), size(ocean_grid, 2)),
                _horizontal_coordinate(wet.φc, i, j, size(ocean_grid, 1), size(ocean_grid, 2)),
            )
        end

        nearest = partialsortperm(distances, 1:nreceivers)
        receiving_measure = routing_weighting == :equal_area_flux ?
            sum(@view wet.wet_area[nearest]) :
            sum(@view wet.wet_volume[nearest])
        receiving_measure > 0 || error(
            "Terrarium runoff routing found zero receiving measure",
        )
        for n in nearest
            i = wet.wet_i[n]
            j = wet.wet_j[n]
            # The contribution weight converts runoff depth rate [m s-1]
            # over the source area [m2] to receiving mass flux
            # [kg m-2 s-1]. Equal-area-flux reproduces the legacy constant
            # value. Equal-column-fraction allocates discharge in proportion
            # to wet column volume, so flux / effective depth is identical
            # across the selected cells. The latter matches the mutable
            # z-star volume pathway, which dilutes the complete wet column.
            weight = if routing_weighting == :equal_area_flux
                freshwater_density * source_area[column] / receiving_measure
            else
                freshwater_density * source_area[column] * wet.wet_volume[n] /
                    (receiving_measure * wet.wet_area[n])
            end
            push!(
                get!(contributions, (i, j), Tuple{Int, Float64}[]),
                (column, weight),
            )
        end
    end

    target_i = Int[]
    target_j = Int[]
    offsets = Int[1]
    contribution_column = Int[]
    contribution_weight = Float32[]
    for target in sort!(collect(keys(contributions)))
        push!(target_i, target[1])
        push!(target_j, target[2])
        for (column, weight) in contributions[target]
            push!(contribution_column, column)
            push!(contribution_weight, Float32(weight))
        end
        push!(offsets, length(contribution_column) + 1)
    end

    return (;
        contribution_column,
        contribution_weight,
        target_i,
        target_j,
        offsets,
        receivers_per_source = nreceivers,
        routing_weighting,
    )
end

function _build_terrarium_runoff_land(
    atmosphere,
    ocean;
    routing_weighting = :equal_area_flux,
    river_mouth_mixing = nothing,
)
    terrarium_state = atmosphere.variables.prognostic.land.terrarium
    surface_runoff = terrarium_state.surface_runoff
    spectral_grid = atmosphere.model.spectral_grid
    longitude, latitude = RG.get_londlatds(spectral_grid.grid)
    longitude = Float64.(Array(longitude))
    latitude = Float64.(Array(latitude))
    land_fraction = Float64.(Array(atmosphere.model.land_sea_mask.mask.data))
    land_points = findall(land_fraction .> 0)
    size(Oceananigans.interior(surface_runoff), 1) == length(land_points) ||
        error("Terrarium runoff columns do not match the SpeedyWeather land mask")

    point_weights = Float64.(Array(_global_point_weights(spectral_grid)))
    sphere_area = 4π * Float64(atmosphere.model.planet.radius)^2
    source_area = point_weights[land_points] .* sphere_area .* land_fraction[land_points]
    routing = _build_runoff_routing(
        longitude[land_points],
        latitude[land_points],
        source_area,
        ocean.model.grid,
        ; routing_weighting,
    )
    arch = Oceananigans.architecture(ocean.model.grid)
    to_architecture(values) = Oceananigans.on_architecture(arch, values)
    return TerrariumRunoffLand(
        surface_runoff,
        to_architecture(routing.contribution_column),
        to_architecture(routing.contribution_weight),
        to_architecture(routing.target_i),
        to_architecture(routing.target_j),
        to_architecture(routing.offsets),
        length(land_points),
        length(routing.target_i),
        routing.receivers_per_source,
        routing.routing_weighting,
        sum(source_area),
        river_mouth_mixing,
    )
end

KernelAbstractions.@kernel function _scatter_terrarium_runoff!(
    freshwater_flux,
    surface_runoff,
    contribution_column,
    contribution_weight,
    target_i,
    target_j,
    offsets,
)
    target = @index(Global, Linear)
    accumulated = zero(eltype(freshwater_flux))
    @inbounds for n in offsets[target]:(offsets[target + 1] - 1)
        column = contribution_column[n]
        accumulated += contribution_weight[n] * surface_runoff[column, 1, 1]
    end
    @inbounds freshwater_flux[target_i[target], target_j[target], 1] = accumulated
end

@inline function _river_mouth_mixing_scale(
    accumulated,
    reference_freshwater_mass_flux_kgm2s,
)
    if reference_freshwater_mass_flux_kgm2s >
       zero(reference_freshwater_mass_flux_kgm2s)
        return clamp(
            accumulated / reference_freshwater_mass_flux_kgm2s -
                oftype(accumulated, 0.5),
            zero(accumulated),
            one(accumulated),
        )
    end
    return ifelse(
        accumulated > zero(accumulated),
        one(accumulated),
        zero(accumulated),
    )
end

@inline _river_mouth_mixing_active(
    accumulated,
    reference_freshwater_mass_flux_kgm2s,
) = _river_mouth_mixing_scale(
    accumulated,
    reference_freshwater_mass_flux_kgm2s,
) > zero(accumulated)

KernelAbstractions.@kernel function _scatter_terrarium_runoff_and_gate_mixing!(
    freshwater_flux,
    surface_runoff,
    contribution_column,
    contribution_weight,
    target_i,
    target_j,
    offsets,
    vertical_diffusivity,
    horizontal_diffusivity,
    eligible_vertical_diffusivity,
    eligible_horizontal_diffusivity,
    active_receiver_mask,
    reference_freshwater_mass_flux_kgm2s,
    Nz,
)
    target = @index(Global, Linear)
    accumulated = zero(eltype(freshwater_flux))
    @inbounds for n in offsets[target]:(offsets[target + 1] - 1)
        column = contribution_column[n]
        accumulated += contribution_weight[n] * surface_runoff[column, 1, 1]
    end
    i = @inbounds target_i[target]
    j = @inbounds target_j[target]
    @inbounds freshwater_flux[i, j, 1] = accumulated
    mixing_scale = _river_mouth_mixing_scale(
        accumulated,
        reference_freshwater_mass_flux_kgm2s,
    )
    active = mixing_scale > zero(mixing_scale)
    @inbounds active_receiver_mask[i, j, 1] = ifelse(
        active,
        one(eltype(active_receiver_mask)),
        zero(eltype(active_receiver_mask)),
    )
    @inbounds for k in 1:Nz
        vertical_diffusivity[i, j, k] = ifelse(
            active,
            mixing_scale * eligible_vertical_diffusivity[i, j, k],
            zero(eltype(vertical_diffusivity)),
        )
        horizontal_diffusivity[i, j, k] = ifelse(
            active,
            mixing_scale * eligible_horizontal_diffusivity[i, j, k],
            zero(eltype(horizontal_diffusivity)),
        )
    end
end

function _scatter_terrarium_runoff_state!(freshwater_flux, land::TerrariumRunoffLand)
    fill!(freshwater_flux, 0)
    land.target_cells == 0 && return nothing
    backend = KernelAbstractions.get_backend(parent(freshwater_flux.data))
    mixing = land.river_mouth_mixing
    if !isnothing(mixing) && mixing.dynamically_gated
        isnothing(mixing.eligible_vertical_diffusivity) && error(
            "dynamic river-mouth mixing lacks vertical eligibility storage",
        )
        isnothing(mixing.eligible_horizontal_diffusivity) && error(
            "dynamic river-mouth mixing lacks horizontal eligibility storage",
        )
        isnothing(mixing.active_receiver_mask) && error(
            "dynamic river-mouth mixing lacks its active-receiver mask",
        )
        _scatter_terrarium_runoff_and_gate_mixing!(backend)(
            freshwater_flux.data,
            land.surface_runoff,
            land.contribution_column,
            land.contribution_weight,
            land.target_i,
            land.target_j,
            land.offsets,
            mixing.vertical_diffusivity.data,
            mixing.horizontal_diffusivity.data,
            mixing.eligible_vertical_diffusivity.data,
            mixing.eligible_horizontal_diffusivity.data,
            mixing.active_receiver_mask.data,
            mixing.reference_freshwater_mass_flux_kgm2s,
            size(mixing.vertical_diffusivity, 3);
            ndrange = land.target_cells,
        )
        Oceananigans.fill_halo_regions!(mixing.vertical_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.horizontal_diffusivity)
        Oceananigans.fill_halo_regions!(mixing.active_receiver_mask)

        # The land interpolation hook is called from NumericalEarth's coupled
        # update, while these coefficient fields are consumed by Oceananigans
        # and the following component step.  KernelAbstractions and
        # Oceananigans both launch asynchronously; establish the producer /
        # consumer boundary here rather than allowing the first ClimaSeaIce
        # launch (or a later module load) to become the accidental completion
        # point.  This is required only for the dynamically mutated closure
        # fields; the static and no-mixing paths do not have this dependency.
        Oceananigans.Architectures.synchronize(
            Oceananigans.architecture(mixing.vertical_diffusivity.grid),
        )
    else
        _scatter_terrarium_runoff!(backend)(
            freshwater_flux.data,
            land.surface_runoff,
            land.contribution_column,
            land.contribution_weight,
            land.target_i,
            land.target_j,
            land.offsets;
            ndrange = land.target_cells,
        )
    end
    return nothing
end

function NumericalEarth.EarthSystemModels.interpolate_state!(
    exchanger,
    grid,
    land::TerrariumRunoffLand,
    coupled_model,
)
    freshwater_flux = exchanger.state.freshwater_flux
    return _scatter_terrarium_runoff_state!(freshwater_flux, land)
end
