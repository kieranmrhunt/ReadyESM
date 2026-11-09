# Included after a test module's R and with(config; kwargs...) definitions.
# Use analytic initial states with no external data. The three new production
# selections are tested separately on their data-backed paths.
function analytic_coupled_test_config(released; output_dir=tempdir())
    # Small analytic tests retain an explicit no-aerosol, 1365 W/m² control.
    # The production PA2 and piControl settings have separate data-backed tests.
    forcing_fields = (; (n => getfield(released.forcing, n)
        for n in fieldnames(typeof(released.forcing)))...)
    forcing = R.ForcingConfig(; merge(forcing_fields, (;
        aerosol_scheme=:none, aerosol_cmip_background_path="",
        aerosol_cmip_plumes_path="", aerosol_cmip_stratospheric_path="",
        solar_constant_w_m2=R.DEFAULT_SOLAR_CONSTANT_W_M2,
        solar_constant_source="", solar_constant_entry=""))...)
    return with(released; device=:cpu, truncation=7, nlayers=8,
        forcing,
        ocean_grid=:idealized, ocean_nlongitude=32, ocean_nlatitude=16,
        ocean_nlayers=4, ocean_initial_conditions=:analytic,
        ecco_initial_conditions_directory="",
        atmosphere_initial_conditions=:analytic,
        terrarium_initial_conditions=:hydrostatic,
        terrarium_evapotranspiration=:bare_ground,
        atmosphere_dissipative_heating=:off,
        terrarium_surface_albedo=:constant, terrarium_albedo_cycle_path="",
        terrarium_albedo_cycle_sha256="", terrarium_vegetation_cycle=:off,
        terrarium_vegetation_cycle_path="", terrarium_vegetation_cycle_sha256="",
        era5_pressure_levels_path="", era5_single_levels_path="", era5_land_state_path="",
        duration_days=4/24, analysis_start_day=0.0, scored_outputs=false, output_dir,
        atmosphere_vertical_advection=:centered_full_temperature)
end
