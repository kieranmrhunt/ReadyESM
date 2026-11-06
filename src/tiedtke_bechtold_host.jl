# Tiedtke-Bechtold convection, stage 1a, part 2: host integration of the pure column (tiedtke_bechtold_column.jl).
#
# Frozen contract: artifacts/development_wave_20261103/tiedtke_bechtold_claude_v1/frozen_v4/FROZEN_CONTRACT_v4.sha256
# (sha256 b07062d259455a207fb45afd1d2a0283387c167f10f8e7a57bb150d36b8c0e76). Specification: DESIGN.md "Integration gate
# after Paul's code audit" items 1-6 and STAGE1_CHOICES_DRAFT.md r15 D43 (outputs, the single unit boundary, the host
# cloud-top rule), T4 (restart), T12 (unit boundary) and T14 (cloud top). The column arithmetic is not changed here.
#
# Two types (DESIGN item 2):
# - `TiedtkeBechtoldConvection`, the scheme. It holds the isbits column parameters, built from the host constants and
#   checked field for field against the constructed model in `initialize!`, and it declares the column work arrays as
#   SpeedyWeather ScratchVariables (DESIGN item 5). It has no `parameterization!` method of its own, so it cannot be
#   used as a model component without the host wrapper.
# - `TiedtkeBechtoldHostConvection`, the dedicated host wrapper. It is not `NetColumnBettsMillerConvection` and it never
#   calls that wrapper's net-humidity rain correction. It is the only unit boundary of D43.
#
# The selection symbol is `atmosphere_convection = tiedtke_bechtold_eqcape_stage1` (DESIGN item 1). The configuration
# refusals are in src/config.jl (`_validate_tiedtke_bechtold_stage1_options`), the factory branch is in
# src/conservative_atmosphere_humidity.jl (`_atmosphere_convection`), and the restart contract is checked in
# src/dynamic_restart.jl before any state is copied.

const TIEDTKE_BECHTOLD_CONVECTION_SYMBOL = :tiedtke_bechtold_eqcape_stage1   # DESIGN item 1
const TIEDTKE_BECHTOLD_STAGE_VERSION = "stage1a"                              # T4 contract value
const CONVECTION_RESTART_CONTRACT_SCHEMA = "readiesm_convection_restart_contract_v1"

# ----------------------------------------------------------------------------------------------------------------
# The scheme and its work arrays (DESIGN item 5)
# ----------------------------------------------------------------------------------------------------------------

"""
    TiedtkeBechtoldConvection(parameters)
    TiedtkeBechtoldConvection(spectral_grid; planet = Earth(spectral_grid), atmosphere = EarthAtmosphere(spectral_grid))

The stage 1a Tiedtke-Bechtold scheme: the isbits `TiedtkeBechtoldParameters` of the pure column. The spectral-grid
constructor builds them with the part-1 helper `tiedtke_bechtold_parameters` from the planet and atmosphere that
`PrimitiveWetModel` constructs by default, with the grid's number format and truncation. `initialize!` then compares
every field with `tiedtke_bechtold_parameters(model)` of the constructed model and refuses any difference, so the
column always uses the live host constants (rule 5 of STAGE1_CHOICES_DRAFT.md).

The scheme declares the column work arrays as ScratchVariables in namespace `:grid` (see
`SpeedyWeather.variables(::TiedtkeBechtoldConvection)`). It has no `parameterization!` method: it runs only inside
`TiedtkeBechtoldHostConvection`.
"""
struct TiedtkeBechtoldConvection{NF <: AbstractFloat} <: SpeedyWeather.AbstractConvection
    parameters::TiedtkeBechtoldParameters{NF}
end

Adapt.@adapt_structure TiedtkeBechtoldConvection

function TiedtkeBechtoldConvection(
    spectral_grid::SpeedyWeather.SpectralGrid;
    planet = SpeedyWeather.Earth(spectral_grid),
    atmosphere = SpeedyWeather.EarthAtmosphere(spectral_grid),
)
    NF = spectral_grid.NF
    return TiedtkeBechtoldConvection(
        tiedtke_bechtold_parameters(NF; planet, atmosphere, truncation = spectral_grid.trunc),
    )
end

"""Refuse a scheme whose parameters differ from those built from the constructed model (live host constants)."""
function _check_tiedtke_bechtold_host_constants(parameters::TiedtkeBechtoldParameters, model)
    live = tiedtke_bechtold_parameters(model)
    typeof(parameters) == typeof(live) || throw(ArgumentError(
        "Tiedtke-Bechtold parameters use $(typeof(parameters)), but the constructed model gives $(typeof(live))",
    ))
    for name in fieldnames(typeof(live))
        getfield(parameters, name) === getfield(live, name) || throw(ArgumentError(
            "Tiedtke-Bechtold parameter $name ($(getfield(parameters, name))) differs from the value built from " *
            "the constructed model ($(getfield(live, name)))",
        ))
    end
    return nothing
end

function SpeedyWeather.initialize!(scheme::TiedtkeBechtoldConvection, model::SpeedyWeather.AbstractModel)
    _check_tiedtke_bechtold_host_constants(scheme.parameters, model)
    return nothing
end

# Work arrays of the column (`tiedtke_bechtold_column_work`), with their units. Layer arrays have one value per full
# level k. Half-level arrays have n+1 values: index j = 2..n+1 (half level j-1/2, the lower half level of layer j-1)
# is stored in a Grid3D at layer j-1, and the top index j = 1 (half level 1/2) in a Grid2D. The column writes every
# entry before it reads it, and after a call the arrays hold the final profiles.
const _TIEDTKE_BECHTOLD_LAYER_WORK = (
    (:entrainment, "kg/m2/s", "updraught entrainment E of each layer (layer-integrated)"),
    (:turbulent_detrainment, "kg/m2/s", "updraught turbulent detrainment D1 of each layer (layer-integrated)"),
    (:organized_detrainment, "kg/m2/s", "updraught organized detrainment D2 of each layer (layer-integrated)"),
    (:generation, "kg/kg/s", "rain generation G"),
    (:rain_evaporation, "kg/kg/s", "sub-cloud rain evaporation e"),
    (:detrained_condensate, "kg/kg/s", "detrained condensate, re-evaporated in the same step"),
)
const _TIEDTKE_BECHTOLD_HALF_LEVEL_WORK = (
    (:temperature_half, "K", "half-level environment temperature"),
    (:humidity_half, "kg/kg", "half-level environment specific humidity"),
    (:geopotential_half, "m2/s2", "half-level geopotential"),
    (:dry_static_energy_half, "J/kg", "half-level dry static energy, stabilised"),
    (:mass_flux_half, "kg/m2/s", "updraught mass flux"),
    (:total_water_flux_half, "kg/m2/s", "net total-water flux F_qt"),
    (:liquid_static_energy_flux_half, "W/m2", "net liquid-static-energy flux F_sl"),
    (:rain_mass_flux_half, "kg/m2/s", "downward rain mass flux P"),
)

"""Scratch name of a work array: `tiedtke_bechtold_<name>` (a Grid3D, layer values or lower half levels)."""
tiedtke_bechtold_scratch_name(name::Symbol) = Symbol("tiedtke_bechtold_", name)
"""Scratch name of the top half level (j = 1) of a half-level work array: `tiedtke_bechtold_<name>_top` (a Grid2D)."""
tiedtke_bechtold_scratch_top_name(name::Symbol) = Symbol("tiedtke_bechtold_", name, "_top")

SpeedyWeather.variables(::TiedtkeBechtoldConvection) = (
    (SpeedyWeather.ScratchVariable(
        tiedtke_bechtold_scratch_name(name), SpeedyWeather.Grid3D();
        desc = "Tiedtke-Bechtold work: $desc", units, namespace = :grid,
    ) for (name, units, desc) in _TIEDTKE_BECHTOLD_LAYER_WORK)...,
    (SpeedyWeather.ScratchVariable(
        tiedtke_bechtold_scratch_name(name), SpeedyWeather.Grid3D();
        desc = "Tiedtke-Bechtold work: $desc at the lower half level of each layer", units, namespace = :grid,
    ) for (name, units, desc) in _TIEDTKE_BECHTOLD_HALF_LEVEL_WORK)...,
    (SpeedyWeather.ScratchVariable(
        tiedtke_bechtold_scratch_top_name(name), SpeedyWeather.Grid2D();
        desc = "Tiedtke-Bechtold work: $desc at the top half level", units, namespace = :grid,
    ) for (name, units, desc) in _TIEDTKE_BECHTOLD_HALF_LEVEL_WORK)...,
)

"""
Column view of a half-level work array (n+1 values, index j holds half level j-1/2): j = 1 is `top[ij]` (a Grid2D)
and j = 2..n+1 is `lower[ij, j-1]` (a Grid3D, the lower half level of layer j-1). Plain array semantics, isbits when
the parents are device arrays, no allocation.
"""
struct _TBHalfLevelColumn{NF, L, U} <: AbstractVector{NF}
    lower::L
    top::U
    ij::Int
    nlayers::Int
end
_TBHalfLevelColumn(lower::L, top::U, ij::Integer, nlayers::Integer) where {L, U} =
    _TBHalfLevelColumn{eltype(lower), L, U}(lower, top, Int(ij), Int(nlayers))
Base.size(x::_TBHalfLevelColumn) = (x.nlayers + 1,)
Base.IndexStyle(::Type{<:_TBHalfLevelColumn}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(x::_TBHalfLevelColumn, j::Int) =
    j == 1 ? x.top[x.ij] : x.lower[x.ij, j - 1]
Base.@propagate_inbounds function Base.setindex!(x::_TBHalfLevelColumn, value, j::Int)
    if j == 1
        x.top[x.ij] = value
    else
        x.lower[x.ij, j - 1] = value
    end
    return x
end

"""Lazy column of pressures p = sigma p_s (Pa) on the host sigma levels (full or half), DESIGN "Column interface v1"."""
struct _TBSigmaPressure{NF, S} <: AbstractVector{NF}
    sigma::S
    surface_pressure::NF
end
Base.size(x::_TBSigmaPressure) = (length(x.sigma),)
Base.IndexStyle(::Type{<:_TBSigmaPressure}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(x::_TBSigmaPressure{NF}, k::Int) where {NF} =
    convert(NF, x.sigma[k]) * x.surface_pressure

"""
    _tiedtke_bechtold_column_call!(dT, dq, ij, vars, scheme, model) -> NamedTuple

Call the pure column on host column `ij`: T, q and geopotential from `vars.grid.temperature_prev`, `humidity_prev` and
`geopotential` (the arrays the Betts-Miller path reads), p_full = sigma_full p_s and p_half = sigma_half p_s from
`vars.grid.pressure_prev` (surface pressure, Pa) and `model.geometry`, dt = `model.time_stepping.Δt_sec` (as the
Betts-Miller path obtains it), the D30 land flag from `model.land_sea_mask.mask` (the land-sea field the other
schemes use) and the latitude from `model.geometry.latds` (diagnostics only, unused by the column). The work arrays
are the scheme's ScratchVariables, sliced per column without allocation. `dT` and `dq` receive the tendencies.
`stage1b` is empty for stage 1a. The stage 1B host passes (shallow parameters, lagged flux, closure), which select the
stage 1B column method (tiedtke_bechtold_stage1b_host.jl).
"""
Base.@propagate_inbounds function _tiedtke_bechtold_column_call!(
    dT, dq, ij, vars, scheme::TiedtkeBechtoldConvection, model, stage1b...,
)
    (; geometry, time_stepping, land_sea_mask) = model
    parameters = scheme.parameters
    nlayers = length(geometry.σ_levels_full)
    scratch = vars.scratch.grid
    surface_pressure = vars.grid.pressure_prev[ij]
    work = (;
        temperature_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_temperature_half, scratch.tiedtke_bechtold_temperature_half_top, ij, nlayers),
        humidity_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_humidity_half, scratch.tiedtke_bechtold_humidity_half_top, ij, nlayers),
        geopotential_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_geopotential_half, scratch.tiedtke_bechtold_geopotential_half_top, ij, nlayers),
        dry_static_energy_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_dry_static_energy_half, scratch.tiedtke_bechtold_dry_static_energy_half_top,
            ij, nlayers),
        mass_flux_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_mass_flux_half, scratch.tiedtke_bechtold_mass_flux_half_top, ij, nlayers),
        total_water_flux_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_total_water_flux_half, scratch.tiedtke_bechtold_total_water_flux_half_top,
            ij, nlayers),
        liquid_static_energy_flux_half = _TBHalfLevelColumn(
            scratch.tiedtke_bechtold_liquid_static_energy_flux_half,
            scratch.tiedtke_bechtold_liquid_static_energy_flux_half_top, ij, nlayers),
        entrainment = view(scratch.tiedtke_bechtold_entrainment, ij, :),
        turbulent_detrainment = view(scratch.tiedtke_bechtold_turbulent_detrainment, ij, :),
        organized_detrainment = view(scratch.tiedtke_bechtold_organized_detrainment, ij, :),
        generation = view(scratch.tiedtke_bechtold_generation, ij, :),
        rain_evaporation = view(scratch.tiedtke_bechtold_rain_evaporation, ij, :),
        detrained_condensate = view(scratch.tiedtke_bechtold_detrained_condensate, ij, :),
    )
    rain_mass_flux_half = _TBHalfLevelColumn(
        scratch.tiedtke_bechtold_rain_mass_flux_half, scratch.tiedtke_bechtold_rain_mass_flux_half_top, ij, nlayers)
    return tiedtke_bechtold_column!(
        dT, dq, rain_mass_flux_half, work,
        view(vars.grid.temperature_prev, ij, :),
        view(vars.grid.humidity_prev, ij, :),
        _TBSigmaPressure(geometry.σ_levels_full, surface_pressure),
        _TBSigmaPressure(geometry.σ_levels_half, surface_pressure),
        view(vars.grid.geopotential, ij, :),
        time_stepping.Δt_sec,
        tb_land_flag(land_sea_mask.mask[ij], parameters),
        geometry.latds[ij],
        parameters,
        stage1b...,
    )
end

# ----------------------------------------------------------------------------------------------------------------
# The host wrapper (DESIGN items 2-4, D43)
# ----------------------------------------------------------------------------------------------------------------

"""
    TiedtkeBechtoldHostConvection(scheme::TiedtkeBechtoldConvection)
    TiedtkeBechtoldHostConvection(spectral_grid; kwargs...)

Host wrapper of the stage 1a Tiedtke-Bechtold column (DESIGN items 2-4, STAGE1_CHOICES_DRAFT.md D43). Per column it
calls the pure column, adds its dT and dq to `vars.tendencies.grid.temperature` and `.humidity`, and converts the
column's own surface rain once, at the single unit boundary:

    rain_rate_convection (m/s) = rain_mass_flux_surface (kg m-2 s-1) / _FRESHWATER_DENSITY_KG_M3,
    rain_convection (m)       += Δt_sec * rain_rate_convection,
    rain_rate (m/s)           += rain_rate_convection           (large-scale rain is already in rain_rate).

Stated equivalence condition (DESIGN item 3, water contract `instantaneous_rain_no_rup_v1`): in each step all net
convective loss of atmospheric water reaches the surface as rain, with no retained condensate, ice or separately
falling precipitation, so sum(dq dp/g) + rain_mass_flux_surface = 0 to rounding. The wrapper only relabels the
column's rain. It never infers rain from the humidity tendency, never repairs a gap, and never writes
`convective_precipitation_mass_correction` (declared only because the coupled diagnostics read it; it stays at its
allocation value zero). A later stage with retained condensate, ice, snow or sedimenting precipitation must bypass
this wrapper (DESIGN item 4).

Cloud top (D43, T14): with no deep cloud the host `cloud_top` is left unchanged; on deep convection it becomes
min(host cloud_top, tb_top_full_level). The process diagnostics `convective_humidity_tendency` and
`convective_temperature_tendency` hold exactly the column's dq and dT, and `convective_precipitation_surface_pressure`
the surface pressure the column used, so `validate_dynamic_diagnostics` checks the column ledger in every coupled run.
"""
struct TiedtkeBechtoldHostConvection{C <: TiedtkeBechtoldConvection} <: SpeedyWeather.AbstractConvection
    convection::C
end

Adapt.@adapt_structure TiedtkeBechtoldHostConvection

TiedtkeBechtoldHostConvection(spectral_grid::SpeedyWeather.SpectralGrid; kwargs...) =
    TiedtkeBechtoldHostConvection(TiedtkeBechtoldConvection(spectral_grid; kwargs...))

SpeedyWeather.initialize!(convection::TiedtkeBechtoldHostConvection, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.initialize!(convection.convection, model)

# TB-named diagnostics (one value per column, from the column result of the current step).
const _TIEDTKE_BECHTOLD_DIAGNOSTICS = (
    (:tiedtke_bechtold_rain_mass_flux_surface, "kg/m2/s",
     "Tiedtke-Bechtold column surface rain mass flux (rain_mass_flux_half at the surface)"),
    (:tiedtke_bechtold_convection_type, "1", "Tiedtke-Bechtold convection type: 0 none, 1 deep"),
    (:tiedtke_bechtold_top_full_level, "1", "Tiedtke-Bechtold top full level (n+1 = none)"),
    (:tiedtke_bechtold_base_half_level, "1", "Tiedtke-Bechtold cloud-base half-level index (0 = none)"),
    (:tiedtke_bechtold_departure_level, "1", "Tiedtke-Bechtold deep departure level (0 = none)"),
    (:tiedtke_bechtold_cloud_base_mass_flux, "kg/m2/s", "Tiedtke-Bechtold closed cloud-base mass flux"),
    (:tiedtke_bechtold_closure_factor, "1", "Tiedtke-Bechtold CAPE closure factor before the limiters"),
    (:tiedtke_bechtold_cfl_factor, "1", "Tiedtke-Bechtold strict CFL limiter factor (1 = inactive)"),
    (:tiedtke_bechtold_positivity_factor, "1", "Tiedtke-Bechtold humidity positivity limiter factor (1 = inactive)"),
    (:tiedtke_bechtold_final_factor, "1", "Tiedtke-Bechtold applied mass-flux factor"),
    (:tiedtke_bechtold_pcape, "J/m3", "Tiedtke-Bechtold unit-flux plume PCAPE"),
    (:tiedtke_bechtold_timescale, "s", "Tiedtke-Bechtold CAPE adjustment time scale"),
    (:tiedtke_bechtold_max_test_cloud_depth, "Pa", "Tiedtke-Bechtold largest test-cloud depth over the departures"),
    (:tiedtke_bechtold_model_top_reached, "1", "One when the Tiedtke-Bechtold plume entered the top layer"),
    (:tiedtke_bechtold_condensation, "kg/m2/s", "Tiedtke-Bechtold column condensation"),
)

SpeedyWeather.variables(convection::TiedtkeBechtoldHostConvection) = (
    SpeedyWeather.variables(convection.convection)...,
    SpeedyWeather.ParameterizationVariable(
        :rain_convection, SpeedyWeather.Grid2D();
        desc = "Convective rain accumulation (Tiedtke-Bechtold column surface rain)", units = "m"),
    SpeedyWeather.ParameterizationVariable(
        :rain_rate_convection, SpeedyWeather.Grid2D();
        desc = "Convective rain rate (Tiedtke-Bechtold column surface rain over water density)", units = "m/s"),
    SpeedyWeather.ParameterizationVariable(
        :rain_rate, SpeedyWeather.Grid2D(); desc = "Rain rate (large-scale + convection)", units = "m/s"),
    SpeedyWeather.ParameterizationVariable(
        :cloud_top, SpeedyWeather.Grid2D(); desc = "Cloud top layer index", units = "1"),
    SpeedyWeather.ParameterizationVariable(
        :convective_humidity_tendency, SpeedyWeather.Grid3D();
        desc = "Tiedtke-Bechtold stage 1a specific-humidity tendency", units = "kg/kg/s"),
    SpeedyWeather.ParameterizationVariable(
        :convective_temperature_tendency, SpeedyWeather.Grid3D();
        desc = "Tiedtke-Bechtold stage 1a temperature tendency", units = "K/s"),
    SpeedyWeather.ParameterizationVariable(
        :convective_precipitation_surface_pressure, SpeedyWeather.Grid2D();
        desc = "Surface pressure used by the Tiedtke-Bechtold column", units = "Pa"),
    SpeedyWeather.ParameterizationVariable(
        :convective_precipitation_mass_correction, SpeedyWeather.Grid2D();
        desc = "Not written by Tiedtke-Bechtold (no rain correction exists): always zero", units = "m/s"),
    (SpeedyWeather.ParameterizationVariable(name, SpeedyWeather.Grid2D(); desc, units)
     for (name, units, desc) in _TIEDTKE_BECHTOLD_DIAGNOSTICS)...,
)

# This runs inside the GPU column kernel, so the message must be one string literal: building it with `*` calls a
# runtime string routine that the GPU compiler refuses (InvalidIRError in GPU gate 65577058).
@noinline _tiedtke_bechtold_cloud_top_error() = throw(ArgumentError(
    "Tiedtke-Bechtold cloud top refused: deep convection must return a full level 1..n and no convection must return n+1 (D43, T14)",
))

"""
    tiedtke_bechtold_host_cloud_top(host_cloud_top, tb_top_full_level, convection_type, nlayers)

D43 host cloud-top rule (T14): a column result with no deep cloud (type 0, top n+1) leaves the host value unchanged;
deep convection (type 1) gives min(host cloud_top, tb_top_full_level). Any other returned index (0, n+2, a deep top of
n+1, a top of 1..n without deep convection, or an unknown type) is refused with an ArgumentError. Never a base, a
half-level index or a zero-based index.
"""
@inline function tiedtke_bechtold_host_cloud_top(host_cloud_top, tb_top_full_level::Integer,
                                                 convection_type::Integer, nlayers::Integer)
    valid = convection_type == 0 ? tb_top_full_level == nlayers + 1 :
            convection_type == 1 ? 1 <= tb_top_full_level <= nlayers : false
    valid || _tiedtke_bechtold_cloud_top_error()
    convection_type == 1 || return host_cloud_top
    return min(host_cloud_top, convert(typeof(host_cloud_top), tb_top_full_level))
end

@inline function _write_tiedtke_bechtold_diagnostics!(p, ij, result)
    NF = typeof(result.rain_mass_flux_surface)
    p.tiedtke_bechtold_rain_mass_flux_surface[ij] = result.rain_mass_flux_surface
    p.tiedtke_bechtold_convection_type[ij] = convert(NF, result.convection_type)
    p.tiedtke_bechtold_top_full_level[ij] = convert(NF, result.tb_top_full_level)
    p.tiedtke_bechtold_base_half_level[ij] = convert(NF, result.cloud_base_half_level)
    p.tiedtke_bechtold_departure_level[ij] = convert(NF, result.departure_level)
    p.tiedtke_bechtold_cloud_base_mass_flux[ij] = result.cloud_base_mass_flux
    p.tiedtke_bechtold_closure_factor[ij] = result.closure_factor
    p.tiedtke_bechtold_cfl_factor[ij] = result.f_CFL
    p.tiedtke_bechtold_positivity_factor[ij] = result.f_pos
    p.tiedtke_bechtold_final_factor[ij] = result.final_factor
    p.tiedtke_bechtold_pcape[ij] = result.PCAPE
    p.tiedtke_bechtold_timescale[ij] = result.tau
    p.tiedtke_bechtold_max_test_cloud_depth[ij] = result.max_test_cloud_depth
    p.tiedtke_bechtold_model_top_reached[ij] = convert(NF, result.model_top_reached)
    p.tiedtke_bechtold_condensation[ij] = result.condensation_column
    return nothing
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    convection::TiedtkeBechtoldHostConvection,
    model,
)
    p = vars.parameterizations
    # The column writes every entry of dT and dq (zeros where it does not act); they are the process diagnostics.
    dT = view(p.convective_temperature_tendency, ij, :)
    dq = view(p.convective_humidity_tendency, ij, :)
    result = _tiedtke_bechtold_column_call!(dT, dq, ij, vars, convection.convection, model)
    temperature_tendency = vars.tendencies.grid.temperature
    humidity_tendency = vars.tendencies.grid.humidity
    nlayers = length(dT)
    for k in 1:nlayers
        temperature_tendency[ij, k] += dT[k]
        humidity_tendency[ij, k] += dq[k]
    end
    # D43: the single unit boundary, with the host constant (never a literal).
    NF = typeof(result.rain_mass_flux_surface)
    rain_rate_convection = result.rain_mass_flux_surface / convert(NF, _FRESHWATER_DENSITY_KG_M3)
    p.rain_rate_convection[ij] = rain_rate_convection
    p.rain_rate[ij] += rain_rate_convection
    p.rain_convection[ij] += model.time_stepping.Δt_sec * rain_rate_convection
    p.convective_precipitation_surface_pressure[ij] = vars.grid.pressure_prev[ij]
    # D43: host cloud top, min rule on deep convection only.
    p.cloud_top[ij] = tiedtke_bechtold_host_cloud_top(
        p.cloud_top[ij], result.tb_top_full_level, result.convection_type, nlayers)
    _write_tiedtke_bechtold_diagnostics!(p, ij, result)
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------
# Restart contract (DESIGN item 5, T4)
# ----------------------------------------------------------------------------------------------------------------

"""
    tiedtke_bechtold_restart_contract() -> NamedTuple

The immutable `convection_restart_contract` record of a Tiedtke-Bechtold stage 1a checkpoint (DESIGN item 5, T4): the
scheme symbol (as a string), the stage version, the water contract, the trigger variant, the parcel-adapter version
and the frozen contract manifest SHA-256, plus a record schema. Plain string fields only, compared by key set and per
key, because JLD2 may reload NamedTuple fields in a different order.
"""
tiedtke_bechtold_restart_contract() = (;
    schema = CONVECTION_RESTART_CONTRACT_SCHEMA,
    scheme = String(TIEDTKE_BECHTOLD_CONVECTION_SYMBOL),
    stage_version = TIEDTKE_BECHTOLD_STAGE_VERSION,
    water_contract = TIEDTKE_BECHTOLD_WATER_CONTRACT,
    trigger_variant = TIEDTKE_BECHTOLD_TRIGGER_VARIANT,
    parcel_adapter = TIEDTKE_BECHTOLD_PARCEL_ADAPTER,
    contract_manifest_sha256 = TIEDTKE_BECHTOLD_CONTRACT_MANIFEST_SHA256,
)

"""The restart contract of a constructed Tiedtke-Bechtold host: stage 1a here, stage 1B in its own file."""
tiedtke_bechtold_restart_contract(::TiedtkeBechtoldHostConvection) = tiedtke_bechtold_restart_contract()
tiedtke_bechtold_restart_contract(convection::PrognosticCondensateConvection) =
    tiedtke_bechtold_restart_contract(convection.convection)

"""Whether a constructed convection component is the Tiedtke-Bechtold host (also under a condensate wrapper)."""
_uses_tiedtke_bechtold_convection(convection) = convection isa TiedtkeBechtoldHostConvection
_uses_tiedtke_bechtold_convection(convection::PrognosticCondensateConvection) =
    _uses_tiedtke_bechtold_convection(convection.convection)

"""
The checkpoint entry: present only with Tiedtke-Bechtold convection, so Betts-Miller (and every other existing scheme)
checkpoints keep exactly their release content.
"""
_convection_restart_state(simulation) =
    _uses_tiedtke_bechtold_convection(simulation.model.convection) ?
        (; convection_restart_contract = tiedtke_bechtold_restart_contract(simulation.model.convection)) : (;)

"""
    _validate_convection_restart(simulation, state)

Refuse a checkpoint whose convection restart contract does not match the reconstructed convection, before any state is
copied. It does not consult the restart-identity override, so `with_dynamic_restart_identity_override` cannot waive it:
- Tiedtke-Bechtold receiver: the checkpoint must carry the contract, with exactly the expected keys and string values.
  A missing contract (a Betts-Miller or other pre-contract checkpoint) is refused.
- Any other receiver (the default Betts-Miller scheme among them): a checkpoint with a contract (a Tiedtke-Bechtold
  checkpoint) is refused. A checkpoint without one is accepted (the legacy rule: no existing scheme writes it).
"""
function _validate_convection_restart(simulation, state)
    convection = simulation.model.convection
    active = _uses_tiedtke_bechtold_convection(convection)
    saved = hasproperty(state, :convection_restart_contract)
    if !active
        saved && throw(ArgumentError(
            "checkpoint carries a convection restart contract (written with Tiedtke-Bechtold convection), but the " *
            "configuration uses $(nameof(typeof(convection))) convection",
        ))
        return nothing
    end
    saved || throw(ArgumentError(
        "the configuration uses Tiedtke-Bechtold convection ($(TIEDTKE_BECHTOLD_CONVECTION_SYMBOL)), but the " *
        "checkpoint has no convection restart contract (it was written with Betts-Miller or another scheme)",
    ))
    record = state.convection_restart_contract
    expected = tiedtke_bechtold_restart_contract(convection)
    names = propertynames(record)
    (length(names) == length(expected) && all(name -> name in names, keys(expected))) || throw(ArgumentError(
        "checkpoint convection restart contract has keys $(Tuple(names)), expected $(keys(expected))",
    ))
    for key in keys(expected)
        value = getproperty(record, key)
        (value isa AbstractString && value == getproperty(expected, key)) || throw(ArgumentError(
            "checkpoint and configuration differ in the convection restart contract field $key " *
            "(checkpoint $(repr(value)), configuration $(repr(getproperty(expected, key))))",
        ))
    end
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------
# Provenance
# ----------------------------------------------------------------------------------------------------------------

"""Diagnostics provenance of the Tiedtke-Bechtold stage 1a convection, empty for every other scheme."""
function _tiedtke_bechtold_metadata(config::ExperimentConfig)
    config.atmosphere_convection == TIEDTKE_BECHTOLD_STAGE1B_CONVECTION_SYMBOL &&
        return _tiedtke_bechtold_stage1b_metadata()                 # tiedtke_bechtold_stage1b_host.jl
    config.atmosphere_convection == TIEDTKE_BECHTOLD_CONVECTION_SYMBOL || return (;)
    return (
        atmosphere_convection_closure = TIEDTKE_BECHTOLD_CLOSURE_NAME,
        atmosphere_convection_stage = TIEDTKE_BECHTOLD_STAGE_VERSION,
        atmosphere_convection_water_contract = TIEDTKE_BECHTOLD_WATER_CONTRACT,
        atmosphere_convection_trigger_variant = TIEDTKE_BECHTOLD_TRIGGER_VARIANT,
        atmosphere_convection_parcel_adapter = TIEDTKE_BECHTOLD_PARCEL_ADAPTER,
        atmosphere_convection_kinetic_energy_version = TIEDTKE_BECHTOLD_KE_VERSION,
        atmosphere_convection_organized_detrainment_version = TIEDTKE_BECHTOLD_ORGDET_VERSION,
        atmosphere_convection_contract_manifest_sha256 = TIEDTKE_BECHTOLD_CONTRACT_MANIFEST_SHA256,
        atmosphere_convection_condensate_policy = "detrained_condensate_reevaporated_in_the_step",
        atmosphere_convection_precipitation_phase = "rain_only_no_convective_snow",
    )
end
