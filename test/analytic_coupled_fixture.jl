# Included after a test module's R and with(config; kwargs...) definitions.
# Use analytic initial states with no external data. The three new production
# selections are tested separately on their data-backed paths.
function analytic_coupled_test_config(released; output_dir=tempdir())
    return with(released; device=:cpu, truncation=7, nlayers=8,
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
        duration_days=4/24, analysis_start_day=0.0, output_dir,
        atmosphere_vertical_advection=:centered_full_temperature)
end
