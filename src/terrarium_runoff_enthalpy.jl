# BS08 Stage B: signed liquid sensible enthalpy carried beside Terrarium's native surface and runoff water.
#
# Reference and units (contract [units]): enthalpy H is relative to liquid water at 0 C (273.15 K), the zero of
# Terrarium's soil internal energy U = T_C C - rho L theta_ice and of the Oceananigans temperature tracer in degC.
# Stores are J m-2, rates W m-2, water depths m and water rates m s-1. Temperatures here are degC.
#
# Every function is a pure scalar of its arguments so the same code runs inside KernelAbstractions kernels and in
# the CPU replay of the GPU gates. The water arguments are the host's own values (host water authority): nothing
# here recomputes a water amount, it only assigns enthalpy to the water the host already moved. Splits compute the
# smaller packet directly and give the larger the subtraction remainder (frozen Stage A rule).

"Liquid water volumetric heat capacity for surface and runoff packets, J m-3 K-1 (1000 kg m-3 x 4181 J kg-1 K-1)."
const BS08_SURFACE_CV = 4.181e6

"""Split a packet of `total_water` carrying `total_enthalpy` into `first_water` and the rest. Returns
(first_enthalpy, second_enthalpy). The smaller packet is computed from the specific enthalpy, the larger by
subtraction, so a tiny remainder keeps its source temperature."""
@inline function bs08_split(total_water, total_enthalpy, first_water)
    z = zero(total_enthalpy)
    iszero(total_water) && return (z, z)
    second_water = total_water - first_water
    specific = total_enthalpy / total_water
    if first_water <= second_water
        first = specific * first_water
        second = iszero(first_water) ? total_enthalpy : total_enthalpy - first
    else
        second = specific * second_water
        first = iszero(second_water) ? total_enthalpy : total_enthalpy - second
    end
    return (first, second)
end

"""Enthalpy rate of the incoming liquid (W m-2). `rainfall` is Terrarium's liquid input (rain + snow-as-liquid,
m s-1), `snowfall` its snow part, `air_temperature` in degC. Rain enters at max(T_air, 0 C); snow-as-liquid at
0 C (0 J), matching native water. Returns (incoming_enthalpy_rate, rain_rate_part)."""
@inline function bs08_incoming_enthalpy_rate(rainfall, snowfall, air_temperature, cv = BS08_SURFACE_CV)
    rain = rainfall - snowfall
    rain_temperature = max(air_temperature, zero(air_temperature))
    return (rain * cv * rain_temperature, rain)
end

"""Assign enthalpy rates to Terrarium's native DirectSurfaceRunoff partition. Host values (all from the same
compute_auxiliary call): store depth `store_water` (m), store enthalpy `store_enthalpy` (J m-2), incoming liquid
`influx` (m s-1) and its enthalpy rate `influx_enthalpy` (W m-2), `drainage` = S/tau_r (m s-1; zero when the store
is empty), `infiltration` and `runoff` (m s-1). Returns a NamedTuple of W m-2 rates:
infiltration, runoff, drainage (removed from the store), store_runoff, and `consistent::Bool` (false when the
host partition is not the native law; the caller must stop before any state change)."""
@inline function bs08_partition_rates(store_water, store_enthalpy, influx, influx_enthalpy, drainage,
                                      infiltration, runoff)
    z = zero(store_enthalpy)
    if store_water > zero(store_water)
        # Case 1 (native): incoming water bypasses to runoff, drainage splits into infiltration and runoff.
        drainage_enthalpy = (store_enthalpy / store_water) * drainage
        infiltration_enthalpy, store_runoff_enthalpy = bs08_split(drainage, drainage_enthalpy, infiltration)
        runoff_enthalpy = influx_enthalpy + store_runoff_enthalpy
        consistent = (infiltration <= drainage) & (runoff == influx + drainage - infiltration)
        return (infiltration = infiltration_enthalpy, runoff = runoff_enthalpy, drainage = drainage_enthalpy,
                store_runoff = store_runoff_enthalpy, consistent = consistent)
    else
        # Case 2 (native): no store, the incoming mixture splits into infiltration and runoff.
        infiltration_enthalpy, runoff_enthalpy = bs08_split(influx, influx_enthalpy, infiltration)
        consistent = iszero(drainage) & iszero(store_enthalpy) & (infiltration <= influx) &
                     (runoff == influx + drainage - infiltration)
        return (infiltration = infiltration_enthalpy, runoff = runoff_enthalpy, drainage = z,
                store_runoff = z, consistent = consistent)
    end
end

"""Enthalpy moved with `water` (m, water + ice volume per area) leaving a soil layer at `temperature` (degC) with
liquid fraction `liquid_fraction`, in Terrarium's convention (J m-2): the liquid part carries `liquid_cv` x T_C,
the ice part carries its in-place energy (`ice_cv` x T_C - `rho_L`). Used for the upward cascade and the downward
deficit fill (debit the source layer, credit the destination layer)."""
@inline function bs08_phase_aware_enthalpy(water, temperature, liquid_fraction, liquid_cv, ice_cv, rho_L)
    liquid = liquid_fraction * water
    ice = water - liquid
    return liquid * liquid_cv * temperature + ice * (ice_cv * temperature - rho_L)
end

"""Top-layer excess pushed to the surface store by the closure. The liquid part leaves the layer with
`liquid_cv` x T_C and enters the store with the same joules. The ice part is melted with the layer's own energy and
enters the store at 0 C: the layer's U is unchanged by it and the store gains 0 J (energy-exact in Terrarium's
convention, contract [classification] exfiltration). Returns (soil_debit, store_credit) in J m-2."""
@inline function bs08_exfiltration_enthalpy(water, temperature, liquid_fraction, liquid_cv)
    liquid = liquid_fraction * water
    transfer = liquid * liquid_cv * temperature
    return (transfer, transfer)
end

"""Consume the pending runoff enthalpy with the native water decision. `pending_water` and `delivered_water` are the
host's amounts (m) before and after its rounded-rate rule; returns (delivered_enthalpy, remaining_enthalpy) in
J m-2 with the smaller packet computed directly."""
@inline function bs08_consume(pending_water, pending_enthalpy, delivered_water)
    return bs08_split(pending_water, pending_enthalpy, delivered_water)
end

"""Ocean freshwater heat-content replacement for the land part, degC m s-1: remove the native land term
T_oc J_land / rho_ocean and add the routed enthalpy H_land / (rho_ocean cp_ocean)."""
@inline function bs08_land_heat_content_correction(land_enthalpy_flux, land_freshwater_flux, ocean_temperature,
                                                   rho_ocean, cp_ocean)
    return land_enthalpy_flux / (rho_ocean * cp_ocean) - ocean_temperature * land_freshwater_flux / rho_ocean
end

# ------------------------------------------------------------------------------------------------------------------
# Host wiring (ExperimentConfig.terrarium_runoff_enthalpy = true only).
#
# The option is carried by state-field presence, like ReadyESM's existing runoff ownership (`runoff_pending_depth`):
# the builder adds the BS08 auxiliaries to Terrarium's `input_variables` only when the option is on, and every hook
# below returns before touching anything when `surface_excess_enthalpy` is absent. Option off therefore allocates
# no field, launches no kernel and takes the native code path unchanged.
# ------------------------------------------------------------------------------------------------------------------

"Number of option-on restores that zero-initialised the BS08 stores from a checkpoint without them."
const BS08_LEGACY_ZERO_INITIALIZATIONS = Ref(0)

const BS08_STATE_KEYS = (
    :surface_excess_enthalpy,
    :runoff_pending_enthalpy,
    :infiltration_enthalpy_flux,
    :runoff_enthalpy_flux,
    :surface_drainage_enthalpy_flux,
    :bs08_partition_consistent,
    :bs08_cumulative_clamp_enthalpy,
)

function _bs08_enthalpy_variables()
    xy = Terrarium.XY()
    return (
        Terrarium.auxiliary(:surface_excess_enthalpy, xy;
            desc = "BS08 surface excess water enthalpy relative to liquid water at 0 C (J m-2)"),
        Terrarium.auxiliary(:runoff_pending_enthalpy, xy;
            desc = "BS08 runoff enthalpy awaiting land-to-ocean exchange (J m-2)"),
        Terrarium.auxiliary(:infiltration_enthalpy_flux, xy;
            desc = "BS08 infiltration enthalpy entering the top soil layer (W m-2, downward)"),
        Terrarium.auxiliary(:runoff_enthalpy_flux, xy;
            desc = "BS08 surface runoff enthalpy (W m-2)"),
        Terrarium.auxiliary(:surface_drainage_enthalpy_flux, xy;
            desc = "BS08 enthalpy drained from the surface store (W m-2)"),
        Terrarium.auxiliary(:bs08_partition_consistent, xy;
            desc = "BS08: 1 where the host surface partition is the native DirectSurfaceRunoff law"),
        Terrarium.auxiliary(:bs08_cumulative_clamp_enthalpy, xy;
            desc = "BS08 enthalpy given to water created by the bottom saturation clamp (J m-2, cumulative, reported)"),
    )
end

_bs08_owns_enthalpy(state) = hasproperty(state, :surface_excess_enthalpy)

"""Test-only stage probe for the g2 hook-local replay. It is `nothing` in every model run, so each call below is one
pointer test on the host. The g2 harness sets it to `(stage, phase, args...) -> ...` to copy the GPU stage inputs
before a BS08 kernel and compare its outputs after it."""
const BS08_STAGE_PROBE = Ref{Any}(nothing)
@inline function _bs08_probe(stage::Symbol, phase::Symbol, args...)
    probe = BS08_STAGE_PROBE[]
    probe === nothing || probe(stage, phase, args...)
    return nothing
end
_bs08_requested(input_variables) =
    any(variable -> Terrarium.varname(variable) == :surface_excess_enthalpy, input_variables)

# Partition: run once per substep after compute_auxiliary! (runoff, infiltration) and the SEB solve (ground heat
# flux), before Terrarium adds the soil top-boundary flux to the internal-energy tendency (soil BCs).
KernelAbstractions.@kernel function _bs08_partition_kernel!(out, grid, fields, canopy_interception, runoff_process)
    i, j = @index(Global, NTuple)
    @inbounds begin
        store_water = fields.surface_excess_water[i, j, 1]
        store_enthalpy = fields.surface_excess_enthalpy[i, j, 1]
        # the host's own incoming liquid, through the same accessor DirectSurfaceRunoff uses
        influx = Terrarium.soil_surface_water_flux(i, j, grid, fields, canopy_interception, nothing)
        influx_enthalpy, _ = bs08_incoming_enthalpy_rate(
            fields.rainfall[i, j, 1], fields.snowfall[i, j, 1], fields.air_temperature[i, j, 1])
        drainage = store_water > zero(store_water) ?
            Terrarium.compute_surface_drainage(runoff_process, store_water) : zero(store_water)
        rates = bs08_partition_rates(store_water, store_enthalpy, influx, influx_enthalpy, drainage,
                                     fields.infiltration[i, j, 1], fields.surface_runoff[i, j, 1])
        out.infiltration_enthalpy_flux[i, j, 1] = rates.infiltration
        out.runoff_enthalpy_flux[i, j, 1] = rates.runoff
        out.surface_drainage_enthalpy_flux[i, j, 1] = rates.drainage
        out.bs08_partition_consistent[i, j, 1] = ifelse(rates.consistent, one(store_water), zero(store_water))
        # Terrarium's top internal-energy boundary is an Oceananigans flux (positive upward, out of the soil).
        # Infiltration enthalpy enters the soil, so it lowers the upward flux.
        out.soil_heat_flux[i, j, 1] = fields.ground_heat_flux[i, j, 1] - rates.infiltration
    end
end

function _bs08_partition!(state, model)
    _bs08_owns_enthalpy(state) || return nothing
    _bs08_probe(:partition, :before, state, model)
    grid = Terrarium.get_grid(model)
    surface_hydrology = model.surface_hydrology
    out = (
        infiltration_enthalpy_flux = state.infiltration_enthalpy_flux,
        runoff_enthalpy_flux = state.runoff_enthalpy_flux,
        surface_drainage_enthalpy_flux = state.surface_drainage_enthalpy_flux,
        bs08_partition_consistent = state.bs08_partition_consistent,
        soil_heat_flux = state.soil_heat_flux,
    )
    state.soil_heat_flux === state.ground_heat_flux && error(
        "BS08 requires a soil heat flux field separate from the ground heat flux")
    fields = (
        surface_excess_water = state.surface_excess_water,
        surface_excess_enthalpy = state.surface_excess_enthalpy,
        rainfall = state.rainfall,
        snowfall = state.snowfall,
        rainfall_ground = state.rainfall_ground,
        air_temperature = state.air_temperature,
        infiltration = state.infiltration,
        surface_runoff = state.surface_runoff,
        ground_heat_flux = state.ground_heat_flux,
    )
    Terrarium.launch!(grid, Terrarium.XY, _bs08_partition_kernel!, out, fields,
        surface_hydrology.canopy_interception, surface_hydrology.surface_runoff)
    minimum(state.bs08_partition_consistent) == 1 || error(
        "BS08: an unclassified surface water source; the host partition is not the native DirectSurfaceRunoff law")
    _bs08_probe(:partition, :after, state, model)
    return nothing
end

# After each substep (Terrarium explicit step + closure done): drain the store enthalpy with the same dt as the water
# and accumulate runoff enthalpy beside the native pending depth.
KernelAbstractions.@kernel function _bs08_after_substep_kernel!(store_enthalpy, pending_enthalpy, drainage_flux,
                                                                runoff_flux, dt)
    column = @index(Global, Linear)
    @inbounds begin
        store_enthalpy[column] -= dt * drainage_flux[column]
        pending_enthalpy[column] += dt * runoff_flux[column]
    end
end

function _bs08_after_substep!(state, dt)
    _bs08_owns_enthalpy(state) || return nothing
    _bs08_probe(:after_substep, :before, state, dt)
    store = Oceananigans.interior(state.surface_excess_enthalpy)
    backend = KernelAbstractions.get_backend(store)
    _bs08_after_substep_kernel!(backend)(
        store,
        Oceananigans.interior(state.runoff_pending_enthalpy),
        Oceananigans.interior(state.surface_drainage_enthalpy_flux),
        Oceananigans.interior(state.runoff_enthalpy_flux),
        dt;
        ndrange = length(store),
    )
    _bs08_probe(:after_substep, :after, state, dt)
    return nothing
end

# Enthalpy-aware replica of Terrarium's adjust_saturation_profile_porosity!: the water arithmetic is copied
# operation for operation (saturation and surface excess water stay bitwise native); each water move also moves
# its phase-aware enthalpy between layer internal energies, and the top excess moves liquid enthalpy to the store.
# Layer temperature and liquid fraction are Terrarium's current auxiliaries (the last energy closure).
KernelAbstractions.@kernel function _bs08_adjust_saturation_kernel!(out, grid, fields, hydrology, strat, bgc,
                                                                   liquid_cv, ice_cv, rho_L)
    i, j = @index(Global, NTuple)
    NF = eltype(out.saturation_water_ice)
    props = Terrarium.get_hydraulic_properties(hydrology)
    sat_min = Terrarium.residual_saturation(props)
    field_grid = Terrarium.get_field_grid(grid)
    N = field_grid.Nz
    sat = out.saturation_water_ice
    surface_excess_water = out.surface_excess_water
    U = out.internal_energy
    T = fields.temperature
    liq = fields.liquid_water_fraction

    @inbounds for k in 1:(N - 1)
        excess_sat = max(sat[i, j, k] - one(NF), zero(NF))
        source_porosity = Terrarium.porosity(i, j, k, grid, fields, strat, bgc)
        destination_porosity = Terrarium.porosity(i, j, k + 1, grid, fields, strat, bgc)
        dz_source = Terrarium.Δzᵃᵃᶜ(i, j, k, field_grid)
        dz_destination = Terrarium.Δzᵃᵃᶜ(i, j, k + 1, field_grid)
        excess_water = excess_sat * source_porosity * dz_source
        sat[i, j, k] -= excess_sat
        sat[i, j, k + 1] += excess_water / (destination_porosity * dz_destination)
        moved = bs08_phase_aware_enthalpy(excess_water, T[i, j, k], liq[i, j, k], liquid_cv, ice_cv, rho_L)
        U[i, j, k] -= moved / dz_source
        U[i, j, k + 1] += moved / dz_destination
    end

    @inbounds for k in N:-1:2
        deficit_sat = max(-sat[i, j, k] + sat_min, zero(NF))
        destination_porosity = Terrarium.porosity(i, j, k, grid, fields, strat, bgc)
        source_porosity = Terrarium.porosity(i, j, k - 1, grid, fields, strat, bgc)
        dz_destination = Terrarium.Δzᵃᵃᶜ(i, j, k, field_grid)
        dz_source = Terrarium.Δzᵃᵃᶜ(i, j, k - 1, field_grid)
        deficit_water = deficit_sat * destination_porosity * dz_destination
        sat[i, j, k] += deficit_sat
        sat[i, j, k - 1] -= deficit_water / (source_porosity * dz_source)
        moved = bs08_phase_aware_enthalpy(deficit_water, T[i, j, k - 1], liq[i, j, k - 1], liquid_cv, ice_cv, rho_L)
        U[i, j, k - 1] -= moved / dz_source
        U[i, j, k] += moved / dz_destination
    end

    @inbounds begin
        excess_sat = max(sat[i, j, N] - one(NF), zero(NF))
        top_porosity = Terrarium.porosity(i, j, N, grid, fields, strat, bgc)
        dz_top = Terrarium.Δzᵃᵃᶜ(i, j, N, field_grid)
        sat[i, j, N] -= excess_sat
        exfiltrated = excess_sat * top_porosity * dz_top
        surface_excess_water[i, j, 1] += exfiltrated
        soil_debit, store_credit = bs08_exfiltration_enthalpy(exfiltrated, T[i, j, N], liq[i, j, N], liquid_cv)
        U[i, j, N] -= soil_debit / dz_top
        out.surface_excess_enthalpy[i, j, 1] += store_credit

        old_bottom = sat[i, j, 1]
        sat[i, j, 1] = max(sat[i, j, 1], sat_min)
        created = (sat[i, j, 1] - old_bottom) * Terrarium.porosity(i, j, 1, grid, fields, strat, bgc) *
            Terrarium.Δzᵃᵃᶜ(i, j, 1, field_grid)
        clamp_enthalpy = created * liquid_cv * T[i, j, 1]
        U[i, j, 1] += clamp_enthalpy / Terrarium.Δzᵃᵃᶜ(i, j, 1, field_grid)
        out.bs08_cumulative_clamp_enthalpy[i, j, 1] += clamp_enthalpy
    end
end

"""Option-on soil hydrology closure: the native sequence of Terrarium's SoilSaturationPressureClosure with the
enthalpy-aware saturation adjustment in place of `adjust_saturation_profile!`."""
function _bs08_hydrology_closure!(state, grid, hydrology, soil, constants)
    _bs08_adjust_saturation!(state, grid, hydrology, soil, constants)
    # the rest of the native SoilSaturationPressureClosure, unchanged
    Terrarium.compute_water_table!(state, grid, hydrology)
    closure = Terrarium.get_closure(hydrology)
    strat = Terrarium.get_stratigraphy(soil)
    bgc = Terrarium.get_biogeochemistry(soil)
    pressure_out = (pressure_head = state.pressure_head,)
    pressure_fields = Terrarium.get_fields(state, hydrology, strat, bgc; except = pressure_out)
    Terrarium.launch!(grid, Terrarium.XYZ, Terrarium.saturation_to_pressure_kernel!,
        pressure_out, pressure_fields, closure, hydrology, strat, bgc)
    return nothing
end

"The enthalpy-aware saturation adjustment alone (one kernel launch), probed for the g2 replay."
function _bs08_adjust_saturation!(state, grid, hydrology, soil, constants)
    _bs08_probe(:adjust_saturation, :before, state, grid, hydrology, soil, constants)
    strat = Terrarium.get_stratigraphy(soil)
    bgc = Terrarium.get_biogeochemistry(soil)
    capacities = soil.energy.thermal_properties.heat_capacities
    rho_L = constants.material.density_water * constants.thermodynamics.latent_heat_fusion
    out = (
        saturation_water_ice = state.saturation_water_ice,
        surface_excess_water = state.surface_excess_water,
        internal_energy = state.internal_energy,
        surface_excess_enthalpy = state.surface_excess_enthalpy,
        bs08_cumulative_clamp_enthalpy = state.bs08_cumulative_clamp_enthalpy,
    )
    fields = merge(Terrarium.get_fields(state, hydrology, strat, bgc),
        (temperature = state.temperature, liquid_water_fraction = state.liquid_water_fraction))
    Terrarium.launch!(grid, Terrarium.XY, _bs08_adjust_saturation_kernel!, out, fields, hydrology, strat, bgc,
        capacities.water, capacities.ice, rho_L)
    _bs08_probe(:adjust_saturation, :after, state, grid, hydrology, soil, constants)
    return nothing
end

"Zero the BS08 stores missing from an option-off checkpoint (the native implicit state); returns the names zeroed."
function _bs08_legacy_zero!(destination::NamedTuple, missing_names)
    for name in missing_names
        name in BS08_STATE_KEYS || error("checkpoint lacks non-BS08 field $name")
        fill!(getproperty(destination, name), 0)
    end
    BS08_LEGACY_ZERO_INITIALIZATIONS[] += 1
    @info "BS08: option-on restore from a checkpoint without enthalpy states; zero-initialised $(collect(missing_names))"
    return missing_names
end
