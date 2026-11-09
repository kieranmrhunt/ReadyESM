using ReadyESM
using Test

include("precompile_mode.jl")
include("lazy_plotting_first_load.jl")
include("precipitation_dispatch_paths.jl")
include("precipitation_coupled_scheduler.jl")
include("production_config.jl")
include("calendar_centred_ecco.jl")
include("moist_mass_coupling.jl")
include("full_temperature_vertical_advection.jl")
include("ordered_ocean_grid.jl")
include("solar_forcing.jl")
include("standalone_atmosphere_restart.jl")
include("dynamic_restart_identity.jl")
include("dynamic_prognostic_guard.jl")
include("carbon_stage0a.jl")
include("carbon_stage0b.jl")
include("carbon_stage0c.jl")
include("carbon_emissions.jl")
include("carbon_diagnostics.jl")
include("carbon_inventory_geometry.jl")

include("convection_timescale.jl")
include("cloud_relative_humidity_min.jl")
include("solar_constant_option.jl")
include("finite_multiplume_convection.jl")
include("convective_momentum_transport.jl")
include("optional_land.jl")
include("rrtmgp_trace_gases.jl")
include("pbl_k_profile.jl")
include("pbl_k_profile_host.jl")
include("pbl_k_profile_strength.jl")
include("pbl_k_profile_momentum.jl")
include("dissipative_heating.jl")
include("well_mixed_gas_schedule.jl")
include("prescribed_ozone.jl")
include("rrtmgp_band_optics_hook.jl")
include("prescribed_aerosol_cmip.jl")
include("prescribed_aerosol_strat.jl")
include("tropopause_column_loop.jl")

diagnostic_source = read(joinpath(ReadyESM.PROJECT_ROOT, "src", "diagnostics.jl"), String)
dynamic_source = read(joinpath(ReadyESM.PROJECT_ROOT, "src", "dynamic_ocean.jl"), String)
for contract in (
    "RRTMGP.direct_sw_surface_albedo",
    "RRTMGP.diffuse_sw_surface_albedo",
    "RRTMGP direct and diffuse surface albedos differ",
    "RRTMGP surface albedo is not band invariant",
)
    @test occursin(contract, diagnostic_source)
end
for contract in (
    "atmosphere_radiative_surface_albedo",
    "retained_rrtmgp_direct_diffuse_band_invariant_solver_state",
    "radiative surface albedo is outside [0, 1]",
)
    @test occursin(contract, dynamic_source)
end

et = ReadyESM.ERA5PrescribedVegetationEvapotranspiration(Float32)
land_dispatch = which(
    ReadyESM.Terrarium.compute_auxiliary!,
    Tuple{
        Any,
        Any,
        typeof(et),
        typeof(ReadyESM.Terrarium.NoCanopyInterception(Float32)),
        ReadyESM.Terrarium.PhysicalConstants{Float32},
        ReadyESM.Terrarium.AbstractAtmosphere,
        ReadyESM.Terrarium.AbstractSoil,
        Nothing,
        Nothing,
    },
)
@test land_dispatch.module === ReadyESM

include("config_input_files.jl")
include("ecco_centered_velocity.jl")
include("darcy_boundary.jl")
include("sampled_global_diagnostics.jl")
include("sampled_history_output.jl")
include("ocean_history_bounds.jl")
include("native_runoff.jl")
include("catke_tuple_surface_flux.jl")
include("terrarium_surface_flux_average.jl")
include("bs08_runoff_enthalpy.jl")
include("tripolar_runoff_delivery.jl")
include("precipitation_diagnostics.jl")
include("precipitation_pressure.jl")
include("large_scale_melt_cooling.jl")
include("sea_ice_branch_types.jl")
include("sea_ice_lead_closing.jl")
include("sea_ice_lead_closing_hemispheric.jl")
include("sea_ice_branch_observer.jl")
include("conservative_sea_ice_advection.jl")
include("sea_ice_thermal_conservation.jl")
include("land_precision.jl")
include("mixed_phase_saturation.jl")
include("atmosphere_spectral_truncation.jl")
include("scored_outputs.jl")
include("mualem_conductivity.jl")
include("seasonal_albedo.jl")
include("radiative_surface_albedo_diagnostic.jl")
include("seasonal_vegetation.jl")
include("tiedtke_bechtold_column.jl")
include("tiedtke_bechtold_stage1b.jl")
include("climate_writer/runtests.jl")
println("READYESM_CPU_SUITE_PASS")
