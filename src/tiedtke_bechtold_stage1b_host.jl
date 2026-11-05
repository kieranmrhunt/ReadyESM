# Tiedtke-Bechtold stage 1B, PART B: host realisation of the stage 1B column (tiedtke_bechtold_column.jl, the stage 1B
# method of `tiedtke_bechtold_column!`).
#
# Frozen scientific contract: frozen_v6 (TIEDTKE_BECHTOLD_STAGE1B_CONTRACT_MANIFEST_SHA256), rows D48, D79, D81, D82,
# D83 and gates G1-G6 of DESIGN_STAGE2_v1.md Sect. 2b. Host design: DESIGN_PARTB_v0.md (fbf9cebc, a draft for Codex),
# proposals P1-P9. The column arithmetic is not changed here.
#
# D48 adapter `coupled_prescribed_inputs_v1` (P1). In the coupled build the host applies, in step n, the prescribed
# surface fluxes (1 - f_land) ocean + f_land land of vars.prognostic.{ocean,land}.{sensible_heat_flux,
# surface_humidity_flux} (SpeedyWeather Prescribed{Ocean,Land}{Heat,Humidity}Flux). The land step and the coupler's
# `update_net_fluxes!` write those inputs after the step-(n-1) column loop, reset_variables! does not touch them, and
# the checkpoint holds them (a deep copy of `variables`). So TB reads them with `_pbl_area_weighted`, as the K-profile
# host does: G1 (persistent storage outside vars.parameterizations) and G2 (an n-1 value that no column-loop code
# writes, the same in the CPU and GPU orders) hold without a TB copy. TB only reads them (G4, D74).
#
# D81 n = 0 rule (P3): the per-column prognostic flag `tiedtke_bechtold_lagged_flux_set` is 0 at a fresh start. While
# it is 0 the column receives the unset adapter value and runs the stage 1a path itself. After its column call the
# wrapper sets it to 1 for that column only, so the next reader (TB in step n+1) sees inputs written after step n.
#
# Stage 1a (tiedtke_bechtold_host.jl) is unchanged in behaviour: its wrapper, cloud-top rule and restart record stay.

const TIEDTKE_BECHTOLD_STAGE1B_CONVECTION_SYMBOL = :tiedtke_bechtold_eqcape_stage1b   # P6 (name RE)
const TIEDTKE_BECHTOLD_STAGE1B_STAGE_VERSION = "stage1b"                               # D82
const TIEDTKE_BECHTOLD_STAGE1B_PARCEL_ADAPTER = "coupled_prescribed_inputs_v1"        # D48, D82 (P1)
const TIEDTKE_BECHTOLD_STAGE1B_FLAG_NAMESPACE = :convection

"""
    TiedtkeBechtoldStage1BHostConvection(convection, shallow)
    TiedtkeBechtoldStage1BHostConvection(spectral_grid; kwargs...)

Host wrapper of the stage 1B Tiedtke-Bechtold column with the HOST closure (D71). Per column it reads the D48 lagged
flux from the prescribed coupled inputs (P1), calls the stage 1B column, and then does exactly what the stage 1a
wrapper does with the result: it adds dT and dq to the tendencies, converts the column's own surface rain once at the
single unit boundary (D43, D75) and writes the TB diagnostics. Cloud top follows D79 (a): a shallow column leaves the
host `cloud_top` unchanged (`tiedtke_bechtold_stage1b_host_cloud_top`). `initialize!` refuses a model whose surface
heat or humidity flux components are not the four prescribed schemes (P2).
"""
struct TiedtkeBechtoldStage1BHostConvection{C <: TiedtkeBechtoldConvection, S <: TiedtkeBechtoldShallowParameters} <:
       SpeedyWeather.AbstractConvection
    convection::C
    shallow::S
end

Adapt.@adapt_structure TiedtkeBechtoldStage1BHostConvection

function TiedtkeBechtoldStage1BHostConvection(spectral_grid::SpeedyWeather.SpectralGrid; kwargs...)
    convection = TiedtkeBechtoldConvection(spectral_grid; kwargs...)
    return TiedtkeBechtoldStage1BHostConvection(convection, tiedtke_bechtold_shallow_parameters(spectral_grid.NF))
end

"""Refuse a model whose applied surface fluxes are not the prescribed coupled inputs that TB reads (P2)."""
function _check_tiedtke_bechtold_stage1b_surface_fluxes(model)
    heat, humidity = model.surface_heat_flux, model.surface_humidity_flux
    components = (
        (heat, :ocean, SpeedyWeather.PrescribedOceanHeatFlux), (heat, :land, SpeedyWeather.PrescribedLandHeatFlux),
        (humidity, :ocean, SpeedyWeather.PrescribedOceanHumidityFlux),
        (humidity, :land, SpeedyWeather.PrescribedLandHumidityFlux),
    )
    for (flux, part, expected) in components
        (hasproperty(flux, part) && getproperty(flux, part) isa expected) || throw(ArgumentError(
            "Tiedtke-Bechtold stage 1B reads the prescribed coupled surface fluxes (D48 adapter " *
            "$TIEDTKE_BECHTOLD_STAGE1B_PARCEL_ADAPTER), so the $(nameof(typeof(flux))) $part component must be " *
            "$(nameof(expected)); it is " *
            "$(hasproperty(flux, part) ? nameof(typeof(getproperty(flux, part))) : "absent")",
        ))
    end
    return nothing
end

function SpeedyWeather.initialize!(convection::TiedtkeBechtoldStage1BHostConvection, model::SpeedyWeather.AbstractModel)
    SpeedyWeather.initialize!(convection.convection, model)
    _check_tiedtke_bechtold_stage1b_surface_fluxes(model)
    return nothing
end

# Stage 1B diagnostics (D83), one value per column from the current step's column result.
const _TIEDTKE_BECHTOLD_STAGE1B_DIAGNOSTICS = (
    (:tiedtke_bechtold_trigger_path, "1", "Tiedtke-Bechtold trigger path: 1 stage 1a fallback (D81), 2 stage 1B"),
    (:tiedtke_bechtold_lagged_sensible_heat_flux, "W/m2",
     "Tiedtke-Bechtold D48 lagged area-weighted sensible heat flux read this step (upward positive; 0 when unset)"),
    (:tiedtke_bechtold_lagged_humidity_flux, "kg/m2/s",
     "Tiedtke-Bechtold D48 lagged area-weighted humidity flux read this step (upward positive; 0 when unset)"),
    (:tiedtke_bechtold_shallow_test_base_half_level, "1", "Tiedtke-Bechtold surface-parcel cloud base (0 = none)"),
    (:tiedtke_bechtold_shallow_test_top_half_level, "1", "Tiedtke-Bechtold surface-parcel top half level (0 = none)"),
    (:tiedtke_bechtold_shallow_test_depth, "Pa", "Tiedtke-Bechtold surface-parcel test-cloud depth"),
    (:tiedtke_bechtold_surface_parcel_wstar, "m/s", "Tiedtke-Bechtold D47 convective velocity scale w*"),
    (:tiedtke_bechtold_surface_parcel_delta_T, "K", "Tiedtke-Bechtold D46 surface-parcel temperature excess"),
    (:tiedtke_bechtold_surface_parcel_delta_q, "kg/kg", "Tiedtke-Bechtold D46 surface-parcel humidity excess"),
    (:tiedtke_bechtold_shallow_source, "W/m2", "Tiedtke-Bechtold shallow closure source S_sub (HOST closure, D71)"),
    (:tiedtke_bechtold_shallow_moist_static_energy_jump, "J/kg", "Tiedtke-Bechtold cloud-base jump Delta h_base"),
    (:tiedtke_bechtold_shallow_closure_mass_flux, "kg/m2/s",
     "Tiedtke-Bechtold shallow M_base = S_sub / Delta h_base before the D41 limiter"),
)
# Stage 1B counters (D53, D57, D73, D81, D83 r5): persistent sums of the per-call increments, in the checkpoint.
const _TIEDTKE_BECHTOLD_STAGE1B_COUNTERS = (
    (:tiedtke_bechtold_lagged_flux_unset_steps, :lagged_flux_unset, "steps run with the unset lagged flux (D81)"),
    (:tiedtke_bechtold_shallow_nonupward_buoyancy_flux_steps, :shallow_nonupward_buoyancy_flux,
     "shallow tests under a zero or downward lagged buoyancy flux (D53)"),
    (:tiedtke_bechtold_shallow_nonpositive_jump_steps, :shallow_nonpositive_jump,
     "shallow closures refused for Delta h_base <= 0 (D73)"),
    (:tiedtke_bechtold_shallow_nonpositive_source_steps, :shallow_nonpositive_source,
     "shallow closures refused for S_sub <= 0 (D73)"),
    (:tiedtke_bechtold_shallow_candidate_replaced_steps, :shallow_candidate_replaced,
     "shallow candidates replaced by a deep cloud (D57)"),
    (:tiedtke_bechtold_shallow_cfl_capped_steps, :shallow_cfl_capped,
     "shallow columns with the D41 CFL factor below 1 (D76, D83 r5)"),
)

SpeedyWeather.variables(convection::TiedtkeBechtoldStage1BHostConvection) = (
    SpeedyWeather.variables(TiedtkeBechtoldHostConvection(convection.convection))...,
    SpeedyWeather.PrognosticVariable(
        :tiedtke_bechtold_lagged_flux_set, SpeedyWeather.Grid2D();
        desc = "Tiedtke-Bechtold D81 flag: 1 once the column has run a full step (lagged flux set)", units = "1",
        namespace = TIEDTKE_BECHTOLD_STAGE1B_FLAG_NAMESPACE),
    (SpeedyWeather.ParameterizationVariable(name, SpeedyWeather.Grid2D(); desc, units)
     for (name, units, desc) in _TIEDTKE_BECHTOLD_STAGE1B_DIAGNOSTICS)...,
    (SpeedyWeather.ParameterizationVariable(name, SpeedyWeather.Grid2D(); desc = "Tiedtke-Bechtold count of $desc",
                                            units = "1")
     for (name, _, desc) in _TIEDTKE_BECHTOLD_STAGE1B_COUNTERS)...,
)

"""
    tiedtke_bechtold_stage1b_host_cloud_top(host_cloud_top, tb_top_full_level, convection_type, nlayers)

D43 and D79 (a) host cloud-top rule for stage 1B (T14, T26): no cloud (type 0, top n+1) and a shallow cloud (type 2,
top 1..n) leave the host value unchanged, and deep convection (type 1, top 1..n) gives min(host cloud_top,
tb_top_full_level). Every other returned index or type is refused, as in the stage 1a rule, which keeps refusing type 2.
"""
@inline function tiedtke_bechtold_stage1b_host_cloud_top(host_cloud_top, tb_top_full_level::Integer,
                                                         convection_type::Integer, nlayers::Integer)
    valid = convection_type == 0 ? tb_top_full_level == nlayers + 1 :
            convection_type in (1, 2) ? 1 <= tb_top_full_level <= nlayers : false
    valid || _tiedtke_bechtold_cloud_top_error()
    convection_type == 1 || return host_cloud_top
    return min(host_cloud_top, convert(typeof(host_cloud_top), tb_top_full_level))
end

@inline function _write_tiedtke_bechtold_stage1b_diagnostics!(p, ij, result, lagged)
    NF = typeof(result.rain_mass_flux_surface)
    p.tiedtke_bechtold_trigger_path[ij] = convert(NF, result.trigger_path)
    p.tiedtke_bechtold_lagged_sensible_heat_flux[ij] = lagged.sensible_heat_flux
    p.tiedtke_bechtold_lagged_humidity_flux[ij] = lagged.humidity_flux
    p.tiedtke_bechtold_shallow_test_base_half_level[ij] = convert(NF, result.shallow_test_base_half_level)
    p.tiedtke_bechtold_shallow_test_top_half_level[ij] = convert(NF, result.shallow_test_top_half_level)
    p.tiedtke_bechtold_shallow_test_depth[ij] = result.shallow_test_depth
    p.tiedtke_bechtold_surface_parcel_wstar[ij] = result.surface_parcel_wstar
    p.tiedtke_bechtold_surface_parcel_delta_T[ij] = result.surface_parcel_delta_T
    p.tiedtke_bechtold_surface_parcel_delta_q[ij] = result.surface_parcel_delta_q
    p.tiedtke_bechtold_shallow_source[ij] = result.shallow_source
    p.tiedtke_bechtold_shallow_moist_static_energy_jump[ij] = result.shallow_moist_static_energy_jump
    p.tiedtke_bechtold_shallow_closure_mass_flux[ij] = result.shallow_closure_mass_flux
    p.tiedtke_bechtold_lagged_flux_unset_steps[ij] += convert(NF, result.lagged_flux_unset)
    p.tiedtke_bechtold_shallow_nonupward_buoyancy_flux_steps[ij] += convert(NF, result.shallow_nonupward_buoyancy_flux)
    p.tiedtke_bechtold_shallow_nonpositive_jump_steps[ij] += convert(NF, result.shallow_nonpositive_jump)
    p.tiedtke_bechtold_shallow_nonpositive_source_steps[ij] += convert(NF, result.shallow_nonpositive_source)
    p.tiedtke_bechtold_shallow_candidate_replaced_steps[ij] += convert(NF, result.shallow_candidate_replaced)
    p.tiedtke_bechtold_shallow_cfl_capped_steps[ij] += convert(NF, result.shallow_cfl_capped)
    return nothing
end

"""
    _tiedtke_bechtold_lagged_flux(ij, vars, model, NF) -> TiedtkeBechtoldLaggedSurfaceFlux{NF}

D48 adapter `coupled_prescribed_inputs_v1` (P1): the land-fraction-weighted prescribed coupled inputs, read exactly as
the K-profile host reads them (`_pbl_area_weighted`), or the declared unset value while the D81 flag of the column is 0.
"""
Base.@propagate_inbounds function _tiedtke_bechtold_lagged_flux(ij, vars, model, ::Type{NF}) where {NF}
    flag = getproperty(vars.prognostic, TIEDTKE_BECHTOLD_STAGE1B_FLAG_NAMESPACE).tiedtke_bechtold_lagged_flux_set
    flag[ij] > 0 || return TiedtkeBechtoldLaggedSurfaceFlux{NF}()
    land_fraction = NF(model.land_sea_mask.mask[ij])
    sensible = NF(_pbl_area_weighted(vars.prognostic.land.sensible_heat_flux,
                                     vars.prognostic.ocean.sensible_heat_flux, ij, land_fraction))
    humidity = NF(_pbl_area_weighted(vars.prognostic.land.surface_humidity_flux,
                                     vars.prognostic.ocean.surface_humidity_flux, ij, land_fraction))
    return TiedtkeBechtoldLaggedSurfaceFlux(sensible, humidity)
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    convection::TiedtkeBechtoldStage1BHostConvection,
    model,
)
    p = vars.parameterizations
    scheme = convection.convection
    NF = typeof(scheme.parameters.gravity)
    lagged = _tiedtke_bechtold_lagged_flux(ij, vars, model, NF)
    dT = view(p.convective_temperature_tendency, ij, :)
    dq = view(p.convective_humidity_tendency, ij, :)
    result = _tiedtke_bechtold_column_call!(dT, dq, ij, vars, scheme, model,
                                            convection.shallow, lagged, TiedtkeBechtoldHostClosure())
    # As the stage 1a wrapper (tiedtke_bechtold_host.jl), with the same expressions in the same order.
    temperature_tendency = vars.tendencies.grid.temperature
    humidity_tendency = vars.tendencies.grid.humidity
    nlayers = length(dT)
    for k in 1:nlayers
        temperature_tendency[ij, k] += dT[k]
        humidity_tendency[ij, k] += dq[k]
    end
    rain_rate_convection = result.rain_mass_flux_surface / convert(NF, _FRESHWATER_DENSITY_KG_M3)
    p.rain_rate_convection[ij] = rain_rate_convection
    p.rain_rate[ij] += rain_rate_convection
    p.rain_convection[ij] += model.time_stepping.Δt_sec * rain_rate_convection
    p.convective_precipitation_surface_pressure[ij] = vars.grid.pressure_prev[ij]
    # D79 (a): shallow leaves cloud_top unchanged, deep as D43.
    p.cloud_top[ij] = tiedtke_bechtold_stage1b_host_cloud_top(
        p.cloud_top[ij], result.tb_top_full_level, result.convection_type, nlayers)
    _write_tiedtke_bechtold_diagnostics!(p, ij, result)
    _write_tiedtke_bechtold_stage1b_diagnostics!(p, ij, result, lagged)
    # D81 / G2: the column has now run a full step; the next reader is TB in step n+1.
    flag = getproperty(vars.prognostic, TIEDTKE_BECHTOLD_STAGE1B_FLAG_NAMESPACE).tiedtke_bechtold_lagged_flux_set
    flag[ij] = one(eltype(flag))
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------
# Restart contract (D82, P5) and provenance
# ----------------------------------------------------------------------------------------------------------------

"""
The stage 1B `convection_restart_contract` (D82): scheme symbol unchanged, stage version stage 1B, the D82 trigger
variant, the D48 adapter version and the frozen_v6 manifest. The per-key string check of `_validate_convection_restart`
therefore refuses stage 1a to 1B, 1B to 1a and either to or from Betts-Miller.
"""
tiedtke_bechtold_restart_contract(::TiedtkeBechtoldStage1BHostConvection) = (;
    schema = CONVECTION_RESTART_CONTRACT_SCHEMA,
    scheme = String(TIEDTKE_BECHTOLD_CONVECTION_SYMBOL),
    stage_version = TIEDTKE_BECHTOLD_STAGE1B_STAGE_VERSION,
    water_contract = TIEDTKE_BECHTOLD_WATER_CONTRACT,
    trigger_variant = TIEDTKE_BECHTOLD_STAGE1B_TRIGGER_VARIANT,
    parcel_adapter = TIEDTKE_BECHTOLD_STAGE1B_PARCEL_ADAPTER,
    contract_manifest_sha256 = TIEDTKE_BECHTOLD_STAGE1B_CONTRACT_MANIFEST_SHA256,
)

_uses_tiedtke_bechtold_convection(::TiedtkeBechtoldStage1BHostConvection) = true

"""Diagnostics provenance of stage 1B (as stage 1a, with the D82 stage, trigger and adapter values)."""
_tiedtke_bechtold_stage1b_metadata() = (
    atmosphere_convection_closure = TIEDTKE_BECHTOLD_CLOSURE_NAME,
    atmosphere_convection_stage = TIEDTKE_BECHTOLD_STAGE1B_STAGE_VERSION,
    atmosphere_convection_water_contract = TIEDTKE_BECHTOLD_WATER_CONTRACT,
    atmosphere_convection_trigger_variant = TIEDTKE_BECHTOLD_STAGE1B_TRIGGER_VARIANT,
    atmosphere_convection_parcel_adapter = TIEDTKE_BECHTOLD_STAGE1B_PARCEL_ADAPTER,
    atmosphere_convection_shallow_closure = "host_lagged_surface_moist_static_energy_flux_v1",
    atmosphere_convection_kinetic_energy_version = TIEDTKE_BECHTOLD_KE_VERSION,
    atmosphere_convection_organized_detrainment_version = TIEDTKE_BECHTOLD_ORGDET_VERSION,
    atmosphere_convection_contract_manifest_sha256 = TIEDTKE_BECHTOLD_STAGE1B_CONTRACT_MANIFEST_SHA256,
    atmosphere_convection_condensate_policy = "detrained_condensate_reevaporated_in_the_step",
    atmosphere_convection_precipitation_phase = "rain_only_no_convective_snow",
)
